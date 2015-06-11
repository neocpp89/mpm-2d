/**
    \file g_nonlocal_mu2.c
    \author Sachith Dunatunga
    \date 10.09.2014

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "interpolate.h"
#include "particle.h"
#include "node.h"
#include "process.h"
#include "material.h"

/* we need the nodal DOFs to use the node number array */
#include "element.h"
#include "exitcodes.h"

#include <assert.h>

#include "cs_cg.h"

#define jp(x) job->particles[i].x

#undef EMOD
#undef NUMOD

#define dense jp(state[0])
#define gf jp(state[3])
#define eta jp(state[4])
#define gf_local jp(state[5])
#define SZZ_STATE 6
#define szz jp(state[SZZ_STATE])
#define xisq jp(state[7])
#define gammap jp(state[9])
#define gammadotp jp(state[10])

#define MAT_VERSION_STRING "1.0 " __DATE__ " " __TIME__

void create_ngf_stiffness_1(cs *triplets, job_t *job, long int *node_map);
void create_ngf_stiffness_2(cs *triplets, job_t *job, long int *node_map, double *g1);

void calculate_stress(job_t *job);
void solve_diffusion_part(job_t *job);
void calculate_stress_threaded(threadtask_t *task);

static double E, nu, G, K, lambda;
static double mu_s, mu_2, I_0, rho_s, rho_c, d, A;
static double *Kp; //vector which is a solution of the matrix-vector product.
static double *ng_loc; //nodal values of g_local
static double *ng; //nodal values of g
static double *nv; //nodal volumes (nodal mass already computed)
static double *ntau_tr; // nodal values of tau_tr

double calculate_g_local(double tau, double p);
double calculate_g_local_from_g(double tau_tr, double p_tr, double g, double delta_t);

double calculate_g_local(double tau, double p)
{
    double g_local = 0;
    const double S0 = mu_s * p;
    if (tau > S0 && p > 0) {
        const double S2 = mu_2 * p;
        const double zeta = I_0 / (d * sqrt(rho_s));

        if (tau >= S2) {
            fprintf(stderr, "\n%g > %g: %g\n", tau, S2, tau / p);
        }
        assert(tau < S2);

        if (tau < S2) {
            g_local = p * sqrt(p) * zeta * (1.0 - S0 / tau) / (S2 - tau);
        } else {
            g_local = 1;
        }

    }
    return g_local;
}

double calculate_xisq(double tau, double p)
{
    const double cap = 100 * d;
    double xisqc = 0;
    if (p > 0) {
        const double S0 = mu_s * p;
        const double S2 = mu_2 * p;

        if (tau > S2) {
            fprintf(stderr, "\n%g > %g: %g\n", tau, S2, tau / p);
            tau = S2;
        }
        assert(tau <= S2);

        if (tau != S0) {
            xisqc  = (fabs(S2 - tau) * A * A * d * d) / (fabs(tau - S0) * (mu_2 - mu_s));
        } else {
            xisqc = cap;
        }

        if (xisqc > cap) {
            xisqc = cap;
        }
    }
    return xisqc;
}

double calculate_g_local_from_g(double tau_tr, double p_tr, double g, double delta_t)
{
    double g_local = 0;
    const double s = g * G * delta_t + p_tr;
    const double tau_s = mu_s * s;
    const double tau_2 = mu_2 * s;
    if (g >= 0 && tau_tr > tau_s && tau_tr < tau_2 && p_tr > 0) {
        const double zeta = I_0 / (d * sqrt(rho_s));
        g_local = (sqrt(p_tr) / tau_tr) * zeta * s * (tau_tr - tau_s) / (tau_2 - tau_tr);

    }
    return g_local;
}

double calculate_deriv_g_local_from_g(double tau_tr, double p_tr, double g, double delta_t);
double calculate_deriv_g_local_from_g(double tau_tr, double p_tr, double g, double delta_t)
{
    double dg_localdg = 0;
    if (g >= 0 && tau_tr > 0 && p_tr > 0) {
        const double p_tilde = (p_tr + g * G * delta_t);
        const double mu_g = tau_tr / p_tilde;
        if (mu_g < mu_2 && mu_g > mu_s) {
            const double zeta = I_0 / (d * sqrt(rho_s));
            dg_localdg = zeta * sqrt(p_tr) * ((mu_s * (mu_2 - 2*mu_g) + mu_g * mu_g) / ((mu_2 - mu_g) * (mu_2 - mu_g))) * (G * delta_t) / tau_tr;
        } else {
            dg_localdg = 0;
        }
    }
    return dg_localdg;
}

double calculate_load_from_g(double tau_tr, double p_tr, double g, double delta_t, double v);
double calculate_load_from_g(double tau_tr, double p_tr, double g, double delta_t, double v)
{
    double load = 0;
    if (g >= 0 && tau_tr > 0 && p_tr > 0) {
        const double zeta = I_0 / (d * sqrt(rho_s));
        const double s = g * G * delta_t + p_tr;
        if ((tau_tr / s) < mu_s) {
            load = 0;
        } else if ((tau_tr / s) < mu_2) {
            const double g_local = (sqrt(p_tr) / tau_tr) * zeta * s * (tau_tr - mu_s * s) / (mu_2 * s - tau_tr);
            load = v * g_local; 
        } else {
            load = 1;
        }
    }
    return load;
}

typedef struct p_trial_s {
    double tau_tr;
    double p_tr;

    double t0xx_tr;
    double t0xy_tr;
    double t0yy_tr;
    double t0zz_tr;

    double tau_tau;
    double p_tau;

    double s; // scaling factor
} trial_t;

void trial_step(const particle_t *p, const double dt, trial_t *trial)
{

    const double trD = p->exx_t + p->eyy_t;
    const double dsjxx = lambda * trD + 2.0 * G * p->exx_t + 2 * p->wxy_t * p->sxy;
    const double dsjxy = 2.0 * G * p->exy_t - p->wxy_t * (p->sxx - p->syy);
    const double dsjyy = lambda * trD + 2.0 * G * p->eyy_t - 2 * p->wxy_t * p->sxy;
    const double dsjzz = lambda * trD;

    const double sxx_tr = p->sxx + dt * dsjxx;
    const double sxy_tr = p->sxy + dt * dsjxy;
    const double syy_tr = p->syy + dt * dsjyy;
    const double szz_tr = p->state[SZZ_STATE] + dt * dsjzz;

    const double p_tr = -(sxx_tr + syy_tr + szz_tr) / 3.0;;
    const double t0xx_tr = sxx_tr + p_tr;
    const double t0xy_tr = sxy_tr;
    const double t0yy_tr = syy_tr + p_tr;
    const double t0zz_tr = szz_tr + p_tr;
    const double tau_tr = sqrt(0.5*(t0xx_tr*t0xx_tr + 2*t0xy_tr*t0xy_tr + t0yy_tr*t0yy_tr + t0zz_tr*t0zz_tr));

    trial->tau_tr = tau_tr;
    trial->p_tr = p_tr;
    trial->t0xx_tr = t0xx_tr;
    trial->t0xy_tr = t0xy_tr;
    trial->t0yy_tr = t0yy_tr;
    trial->t0zz_tr = t0zz_tr;

    return;
}

// assume an active, valid particle is passed in to p
// retuns dot(bar(gamma))^p
//double local_step(particle_t *p, double dt, double * restrict tau_tau, double * restrict p_tau, int * restrict flag);
double local_step(const particle_t *p, double dt, trial_t *trial, int *flag);

// flag is set to 0 if the material is open, 1 if it is dense.
double local_step(const particle_t *p, double dt, trial_t *trial, int *flag)
{
    trial_step(p, dt, trial);

    const double rho = (p->m / p->v);
    const double tau_tr = trial->tau_tr;
    const double p_tr = trial->p_tr;

    double nup_tau;
    double tau_tau;
    trial->s = 0;
    if (rho < rho_c || p_tr <= 0) {
        *flag = 0;
        nup_tau = (tau_tr) / (G * dt);
        tau_tau = 0;
    } else if (p_tr > 0) {
        *flag = 1;
        const double S0 = mu_s * p_tr;
        double scale_factor = 1.0;
        tau_tau = tau_tr;
        if (tau_tr > S0) {
            const double S2 = mu_2 * p_tr;
            const double alpha = G * I_0 * dt * sqrt(p_tr / rho_s) / d;
            const double B = S2 + tau_tr + alpha;
            const double H = S2 * tau_tr + S0 * alpha;
            tau_tau = 2.0 * H / (B + sqrt(B * B - 4 * H));
            scale_factor = (tau_tau / tau_tr);
        }

        assert(scale_factor <= 1.0);
        assert(scale_factor > 0);
        nup_tau = tau_tr * (1.0 - scale_factor) / G / dt;
        trial->s = scale_factor;
    } else {
        fprintf(stderr, "u");
        *flag = 0;
        tau_tau = 0;
        nup_tau = 0;
    }
    trial->tau_tau = tau_tau;
    if (p_tr > 0) {
        trial->p_tau = p_tr;
    } else {
        trial->p_tau = 0;
    }
    return nup_tau;
}

void cs_print_to_file(const cs *A)
{
    int m, n, *Ap, *Ai;
    double *Ax ;
    if (!A) { return ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    
    double *A_dense = calloc(sizeof(double), m*n);
    FILE *fp = fopen("matrix.cs", "w");

    for (int j = 0; j < n; j++) {
        for (int p = Ap [j]; p < Ap [j+1]; p++) {
            A_dense[Ai[p]*n + j] = Ax[p];
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(fp, "%lg,", A_dense[i*n + j]);
        }
        fprintf(fp, "\n");
    }

    if (fp != NULL) {
        fclose(fp);
    }
    free(A_dense);

    return;
}

void create_ngf_stiffness_1(cs *triplets, job_t *job, long int *node_map)
{
    /* create 1st stiffness matrix. */
    for (size_t i = 0; i < job->num_particles; i++) {
        // collapsed other conditionals into the dense flag
        if (dense == 0) {
            continue;
        }
        double s[4];
        double grad_s[4][2];

        s[0] = job->h1[i];
        s[1] = job->h2[i];
        s[2] = job->h3[i];
        s[3] = job->h4[i];

        grad_s[0][0] = job->b11[i];
        grad_s[1][0] = job->b12[i];
        grad_s[2][0] = job->b13[i];
        grad_s[3][0] = job->b14[i];

        grad_s[0][1] = job->b21[i];
        grad_s[1][1] = job->b22[i];
        grad_s[2][1] = job->b23[i];
        grad_s[3][1] = job->b24[i];

        const size_t p = job->in_element[i];
        int *nn = job->elements[p].nodes;

        for (int ei = 0; ei < NODES_PER_ELEMENT; ei++) {
            for (int ej = 0; ej < NODES_PER_ELEMENT; ej++) {
                size_t gi = nn[ei];
                size_t gj = nn[ej];

                gi = (job->node_number_override[NODAL_DOF * gi + 0] - 0) / NODAL_DOF;
                gj = (job->node_number_override[NODAL_DOF * gj + 0] - 0) / NODAL_DOF;

                const size_t sgi = node_map[gi];
                const size_t sgj = node_map[gj];

                assert(isfinite(xisq));
                assert(isfinite(job->particles[i].v));

                const double k_component = -job->particles[i].v * (s[ei] * s[ej] + 
                        xisq * (grad_s[ei][0]*grad_s[ej][0] + grad_s[ei][1]*grad_s[ej][1]));

                assert(isfinite(k_component));
                
                cs_entry(triplets, sgi, sgj, k_component);
            }
        }
    }

    return;
}

void create_ngf_stiffness_2(cs *triplets, job_t *job, long int *node_map, double *g1)
{
    double *ng = g1;

    /* create 1st stiffness matrix. */
    for (size_t i = 0; i < job->num_particles; i++) {
        // collapsed other conditionals into the dense flag
        if (dense == 0) {
            continue;
        }
        trial_t tr;
        trial_step(&(job->particles[i]), job->dt, &tr);

        double s[4];
        double grad_s[4][2];

        s[0] = job->h1[i];
        s[1] = job->h2[i];
        s[2] = job->h3[i];
        s[3] = job->h4[i];

        grad_s[0][0] = job->b11[i];
        grad_s[1][0] = job->b12[i];
        grad_s[2][0] = job->b13[i];
        grad_s[3][0] = job->b14[i];

        grad_s[0][1] = job->b21[i];
        grad_s[1][1] = job->b22[i];
        grad_s[2][1] = job->b23[i];
        grad_s[3][1] = job->b24[i];

        const size_t p = job->in_element[i];
        int *nn = job->elements[p].nodes;

        double reconstructed_g = 0;
        for (int ei = 0; ei < NODES_PER_ELEMENT; ei++) {
            size_t gi = nn[ei];
            gi = (job->node_number_override[NODAL_DOF * gi + 0 ] - 0) / NODAL_DOF;
            const size_t sgi = node_map[gi];

            reconstructed_g += ng[sgi] * s[ei];
        }

        const double gloc_from_g = calculate_g_local_from_g(tr.tau_tr, tr.p_tr, reconstructed_g, job->dt);
        
        for (int ei = 0; ei < NODES_PER_ELEMENT; ei++) {
            size_t gi = nn[ei];
            gi = (job->node_number_override[NODAL_DOF * gi + 0 ] - 0) / NODAL_DOF;
            const size_t sgi = node_map[gi];

            const double ng_loc_component = (job->particles[i].v) * gloc_from_g * s[ei];
            if (ng_loc_component < 0) {
                printf("v: %lg\ng_loc^1: %g\n", job->particles[i].v, ng_loc_component);
                printf("reconstructed_g: %g\n", reconstructed_g);
                printf("g_loc_from_g: %g\n", gloc_from_g);
            }
            assert(ng_loc_component >= 0);

            ng_loc[sgi] += (-ng_loc_component);
        }

        for (int ei = 0; ei < NODES_PER_ELEMENT; ei++) {
            for (int ej = 0; ej < NODES_PER_ELEMENT; ej++) {
                size_t gi = nn[ei];
                size_t gj = nn[ej];
 
                gi = (job->node_number_override[NODAL_DOF * gi + 0] - 0) / NODAL_DOF;
                gj = (job->node_number_override[NODAL_DOF * gj + 0] - 0) / NODAL_DOF;

                const size_t sgi = node_map[gi];
                const size_t sgj = node_map[gj];

                assert(isfinite(xisq));
                assert(isfinite(job->particles[i].v));

                const double dglocdg = calculate_deriv_g_local_from_g(tr.tau_tr, tr.p_tr, reconstructed_g, job->dt);
                const double k_component = -job->particles[i].v * ((1 - dglocdg) * s[ei] * s[ej] + 
                        xisq * (grad_s[ei][0]*grad_s[ej][0] + grad_s[ei][1]*grad_s[ej][1]));

                assert(isfinite(k_component));
                
                cs_entry(triplets, sgi, sgj, k_component);
            }
        }
    }

    return;
}

/*----------------------------------------------------------------------------*/
void material_init(job_t *job)
{
    for (size_t i = 0; i < job->num_particles; i++) {
        for (size_t j = 0; j < DEPVAR; j++) {
            job->particles[i].state[j] = 0;
        }
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        if ((job->particles[i].m / job->particles[i].v) > rho_c) {
            dense = 1;
        } else {
            dense = 0;
        }
        eta = 0;
        gammap = 0;
        gammadotp = 0;
        gf = 0;
        gf_local = 0;
        xisq = 0;
        szz = 0.5 * (job->particles[i].sxx + job->particles[i].syy); // ugly hack
    }

    if (job->material.num_fp64_props < 9) {
        fprintf(stderr,
            "%s:%s: Need at least 9 properties defined (%s).\n",
            "E, nu, mu_s, mu_2, I_0, rho_s, rho_c, d, A",
            __FILE__, __func__);
        exit(EXIT_ERROR_MATERIAL_FILE);
    } else {
        E = job->material.fp64_props[0];
        nu = job->material.fp64_props[1];
        mu_s = job->material.fp64_props[2];
        mu_2 = job->material.fp64_props[3];
        I_0 = job->material.fp64_props[4];
        rho_s = job->material.fp64_props[5];
        rho_c = job->material.fp64_props[6];
        d = job->material.fp64_props[7];
        A = job->material.fp64_props[8];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2*nu));
        lambda = K - 2.0 * G / 3.0;
        printf("%s:%s: properties (E = %g, nu = %g, G = %g, K = %g, mu_s = %g,"
            " mu_2 = %g, I_0 = %g, rho_s = %g, rho_c = %g, d = %g, A = %g).\n",
            __FILE__, __func__, E, nu, G, K, mu_s, mu_2, I_0, rho_s, rho_c, d, A);
    }

    Kp = calloc(job->num_nodes, sizeof(double));
    printf("%s:%s: Done allocating storage for matrix-vector multiply.\n",
        __FILE__, __func__);
    ng = calloc(job->num_nodes, sizeof(double));
    ng_loc = calloc(job->num_nodes, sizeof(double));
    ntau_tr = calloc(job->num_nodes, sizeof(double));
    printf("%s:%s: Done allocating storage for nodal values of g.\n",
        __FILE__, __func__);
    nv = calloc(job->num_nodes, sizeof(double));
    printf("%s:%s: Done allocating storage for nodal values of v.\n",
        __FILE__, __func__);
    printf("%s:%s: (material version %s) done initializing material.\n",
        __FILE__,  __func__, MAT_VERSION_STRING);
    return;
}
/*----------------------------------------------------------------------------*/

/* nonlocal granular fluidity model. */
void calculate_stress_threaded(threadtask_t *task)
{
    job_t *job = task->job;

    /* Only one of the threads can compute the nonlocal solution. */
    if (task->id == 0) {
        calculate_stress(job);
    }

    return;
}

/*----------------------------------------------------------------------------*/
void calculate_stress(job_t *job)
{
    /* g_local is calculated when we create the stiffness matrix. */
    /* build FEM diffusion array/load vector and solve for g_nonlocal */
    solve_diffusion_part(job);

    for (size_t i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        double nup_tau = 0;
        trial_t tr;
        trial_step(&(job->particles[i]), job->dt, &tr);

        const double p_tr = tr.p_tr;
        const double t0xx_tr = tr.t0xx_tr;
        const double t0xy_tr = tr.t0xy_tr;
        const double t0yy_tr = tr.t0yy_tr;
        const double t0zz_tr = tr.t0zz_tr;
        const double tau_tr = tr.tau_tr;

        if (dense) {
            const double scale_factor = p_tr / (G * job->dt * gf + p_tr);
            const double tau_tau = tau_tr * scale_factor;

            nup_tau = ((tau_tr - tau_tau) / G) / job->dt;

            job->particles[i].sxx = scale_factor * t0xx_tr - p_tr;
            job->particles[i].sxy = scale_factor * t0xy_tr;
            job->particles[i].syy = scale_factor * t0yy_tr - p_tr;
            job->particles[i].state[SZZ_STATE] = scale_factor * t0zz_tr - p_tr;
        } else {
            // local step should already take care of this.
            /* 
            nup_tau = (tau_tr / G) / job->dt;
            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
            job->particles[i].state[SZZ_STATE] = 0;
            */
        }

        /* use strain rate to calculate stress increment */
        gammap += nup_tau * job->dt;
        gammadotp = nup_tau;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void solve_diffusion_part(job_t *job)
{
    size_t i;
    size_t i_new;
    size_t ei, ej;
    size_t gi, gj;
    long int sgi, sgj;

    long int *node_map;
    long int *inv_node_map;

    size_t slda;
    size_t nnz;

    cs *triplets;
    cs *smat;

    cs *triplets_hat;
    cs *smat_hat;

    double *gf_nodes;

    double *f;
    double *f_old;
    double *residual;

    double s[4];
    double grad_s[4][2];

    int p;
    int *nn;

    /* initialize node level fluidity arrays */
    gf_nodes = (double *)malloc(job->num_nodes * sizeof(double));
    for (i = 0; i < job->num_nodes; i++) {
        gf_nodes[i] = 0;
        ng[i] = 0;
        ng_loc[i] = 0;
        ntau_tr[i] = 0;
        // clear volumes
        nv[i] = 0;
    }

    /* save xisq from beginning of step so we don't have to recalculate all the time. */
    double *xisq_particles_initial = calloc(job->num_particles, sizeof(double));
    double *dglocdg_particles = calloc(job->num_particles, sizeof(double));


    /* get number of dofs and initialize mapping arrays */
    node_map = malloc(job->num_nodes * sizeof(long int));
    inv_node_map = malloc(job->num_nodes * sizeof(long int));
    slda = 0;
    for (i = 0; i < job->num_nodes; i++) {
        node_map[i] = -1;
        inv_node_map[i] = -1;
    }
    for (i = 0; i < job->num_nodes; i++) {
        i_new = (job->node_number_override[NODAL_DOF * i + 0 ] - 0) / NODAL_DOF;

        if (node_map[i_new] != -1) {
            continue;
        }

        if (job->nodes[i_new].m > 0) {
            node_map[i_new] = slda;
            inv_node_map[slda] = i_new;
            slda++;
        }
    }

    if (slda == 0) {
        fprintf(stderr, "%s:%s: Invalid size, slda = %zu.\n",
            __FILE__, __func__, slda);
        exit(-1);
    }

    f = calloc(sizeof(double), slda);
    f_old = calloc(sizeof(double), slda);
    residual = calloc(sizeof(double), slda);

    /*  calculate number of nonzero elements (before summing duplicates). */
    nnz = 0;
    for (i = 0; i < job->num_particles; i++) {
        if (dense) {
            nnz += (NODES_PER_ELEMENT * NODES_PER_ELEMENT) * 1;
        }
    }
    nnz += slda;

    // initial g_loc from local step
    for (size_t i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            dense = 0;
            continue; 
        }
        p = job->in_element[i];
        if (p == -1) {
            dense = 0;
            continue;
        }
        int flag = 0;
        trial_t tr;
        local_step(&(job->particles[i]), job->dt, &tr, &flag);
        dense = flag;
        if (dense) {
            xisq = calculate_xisq(tr.tau_tau, tr.p_tau);
            const double g_loc = calculate_g_local(tr.tau_tau, tr.p_tau);

            s[0] = job->h1[i];
            s[1] = job->h2[i];
            s[2] = job->h3[i];
            s[3] = job->h4[i];

            p = job->in_element[i];
            nn = job->elements[p].nodes;

            assert(isfinite(gf_local));
            assert(isfinite(xisq));
            
            for (ei = 0; ei < NODES_PER_ELEMENT; ei++) {
                gi = nn[ei];
                gi = (job->node_number_override[NODAL_DOF * gi + 0 ] - 0) / NODAL_DOF;
                sgi = node_map[gi];

                const double ng_loc_component = (job->particles[i].v) * g_loc * s[ei];
                if (ng_loc_component < 0) {
                    printf("v: %lg\ng_loc^1: %g\n", job->particles[i].v, ng_loc_component);
                }
                assert(ng_loc_component >= 0);

                ng_loc[sgi] += (-ng_loc_component);
                nv[sgi] += job->particles[i].v * s[ei];
            }
        } else {
            // set stresses to 0
            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
            job->particles[i].state[SZZ_STATE] = 0;
        }
    }

    double rel_error = 1;
    const double max_rel_error = 1e-6;
    int inner_iterations = 0;
    const int max_inner_iterations = 1000;
    do {
    /* slda contains degrees of freedom of new matrix */
    triplets = cs_spalloc(slda, slda, nnz, 1, 1);
    triplets_hat = cs_spalloc(slda, slda, nnz, 1, 1); 

    for (i = 0; i < slda; i++) {
        f[i] = ng_loc[i]; // copy load for solution
        f_old[i] = ng_loc[i]; // copy load for later
    }

    create_ngf_stiffness_1(triplets, job, node_map);

//#define DIRECT_SOLVE
#ifdef DIRECT_SOLVE
    /* keep matrix from being dengerate when an element is filled with open particles. */
    for (i = 0; i < slda; i++) {
        cs_entry(triplets, i, i, 1e-10);
    }
#endif
    /* create compressed sparse matrix */
    smat = cs_compress(triplets);
    cs_dupl(smat);

#ifndef DIRECT_SOLVE
//    fprintf(stderr, "%d by %d\n", smat->m, smat->n); // print out matrix for debugging
    if (!cs_cg(smat, f, ng_loc, 1e-7)) {
        fprintf(stderr, "cg error!\n");
        //cs_print(smat, 0); // print out matrix for debugging
        cs_print_to_file(smat);
        exit(EXIT_ERROR_CS_SOL);
    }
#else
//    fprintf(stderr, "%d by %d\n", smat->m, smat->n); // print out matrix for debugging
    if (!cs_lusol(1, smat, f, 1e-15)) {
        fprintf(stderr, "lusol error!\n");
        if (cs_qrsol(1, smat, f)) {
            fprintf(stderr, "qrsol error!\n");
            //cs_print(smat, 0); // print out matrix for debugging
            cs_print_to_file(smat);
            exit(EXIT_ERROR_CS_SOL);
        }
    }
#endif

    for (i = 0; i < slda; i++) {
        if (ng[i] < 0) {
            fprintf(stderr, "g1[%d] = %g\n", i, ng[i]);
        }
        assert(ng[i] >= 0);
    }

    for (i = 0; i < slda; i++) {
        ng[i] = f[i]; // solution is g1, save it as ng
        gf_nodes[i] = f[i];
    }

    create_ngf_stiffness_2(triplets_hat, job, node_map, ng);

    // calculate g_loc from g
    for (size_t i = 0; i < job->num_particles; i++) {
        if (dense) {
            trial_t tr;
            trial_step(&(job->particles[i]), job->dt, &tr);

            s[0] = job->h1[i];
            s[1] = job->h2[i];
            s[2] = job->h3[i];
            s[3] = job->h4[i];

            p = job->in_element[i];
            nn = job->elements[p].nodes;

            double reconstructed_g = 0;
            for (ei = 0; ei < NODES_PER_ELEMENT; ei++) {
                gi = nn[ei];
                gi = (job->node_number_override[NODAL_DOF * gi + 0 ] - 0) / NODAL_DOF;
                sgi = node_map[gi];

                reconstructed_g += ng[sgi] * s[ei];
            }

            const double gloc_from_g = calculate_g_local_from_g(tr.tau_tr, tr.p_tr, reconstructed_g, job->dt);
            
            for (ei = 0; ei < NODES_PER_ELEMENT; ei++) {
                gi = nn[ei];
                gi = (job->node_number_override[NODAL_DOF * gi + 0 ] - 0) / NODAL_DOF;
                sgi = node_map[gi];

                const double ng_loc_component = (job->particles[i].v) * gloc_from_g * s[ei];
                if (ng_loc_component < 0) {
                    printf("v: %lg\ng_loc^1: %g\n", job->particles[i].v, ng_loc_component);
                }
                assert(ng_loc_component >= 0);

                ng_loc[sgi] += (-ng_loc_component);
            }
        }
    }

    double *f_hat = malloc(sizeof(double) * slda);
    double *guess = malloc(sizeof(double) * slda);
    for (i = 0; i < slda; i++) {
        f_hat[i] = ng_loc[i] - f_old[i]; // - g_loc_hat(g_1) + g_loc^1
        guess[i] = 0;
    }

#ifdef DIRECT_SOLVE
    /* keep matrix from being dengerate when an element is filled with open particles. */
    for (i = 0; i < slda; i++) {
        cs_entry(triplets_hat, i, i, 1e-10);
    }
#endif
    /* create compressed sparse matrix */
    smat_hat = cs_compress(triplets_hat);
    cs_dupl(smat_hat);

#ifndef DIRECT_SOLVE
//    fprintf(stderr, "%d by %d\n", smat_hat->m, smat_hat->n); // print out matrix for debugging
    if (!cs_cg(smat_hat, f_hat, guess, 1e-7)) {
        fprintf(stderr, "cg error!\n");
        //cs_print(smat_hat, 0); // print out matrix for debugging
        cs_print_to_file(smat_hat);
        exit(EXIT_ERROR_CS_SOL);
    }
#else
//    fprintf(stderr, "%d by %d\n", smat_hat->m, smat_hat->n); // print out matrix for debugging
    if (!cs_lusol(1, smat_hat, f_hat, 1e-15)) {
        fprintf(stderr, "lusol error!\n");
        if (cs_qrsol(1, smat_hat, f_hat)) {
            fprintf(stderr, "qrsol error!\n");
            //cs_print(smat_hat, 0); // print out matrix for debugging
            cs_print_to_file(smat_hat);
            exit(EXIT_ERROR_CS_SOL);
        }
    }
#endif

    for (i = 0; i < slda; i++) {
        ng[i] += f_hat[i]; // ng now contains g1 + delta g2 = g
    }

    cs_spfree(triplets);
    cs_spfree(smat);

    cs_spfree(triplets_hat);
    cs_spfree(smat_hat);

    for (i = 0; i < slda; i++) {
        // f_old and ng_loc are actually -gloc^1
        // residual[i] = ng_loc[i] - f_old[i];
        residual[i] = f_hat[i]; // try to just use size of delta g2
        // fprintf(stderr, "r[%d] = %g\n", i, residual[i]);
    }

    free(f_hat);
    free(guess);

    const double rtr = dot(residual, residual, slda);
    rel_error =  rtr / slda;

    inner_iterations++;
    printf("%d: %d %g %zu %g\n", job->stepcount, inner_iterations, rtr, slda, rel_error);
    } while ((inner_iterations < max_inner_iterations) && (rel_error > max_rel_error));

    /* map nodal g_nonlocal back to particles */
    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        if (dense == 0) {
            continue;
        }

        p = job->in_element[i];
        if (p == -1) {
            continue;
        }

        s[0] = job->h1[i];
        s[1] = job->h2[i];
        s[2] = job->h3[i];
        s[3] = job->h4[i];

        gf = 0;
        nn = job->elements[p].nodes;

        for (ei = 0; ei < NODES_PER_ELEMENT; ei++) {
            gi = nn[ei];
            gi = (job->node_number_override[NODAL_DOF * gi + 0] - 0) / NODAL_DOF;
            sgi = node_map[gi];
/*
            if (gf_nodes[gi] < 0) {
                fprintf(stderr, "%d: %lg\n", gi, gf_nodes[gi]);
                cs_print_to_file(smat);
                FILE *ff = fopen("load.cs", "w");
                for (size_t j = 0; j < slda; j++) {
                    fprintf(ff, "%lg\n", p_tr[j]);
                }
                fclose(ff);
            }
            assert(gf_nodes[gi] >= 0);
*/
            gf += ng[sgi] * s[ei];
            // gf += gf_nodes[gi] * 0.25;
        }

        if (gf < 0) {
            fprintf(stderr, "%zu: %lg\n", i, gf);
            // cs_print_to_file(smat);
            FILE *ff = fopen("load.cs", "w");
            for (size_t j = 0; j < slda; j++) {
                fprintf(ff, "%lg\n", ng[j]);
            }
            fclose(ff);
            gf = 0;
        }
        assert(gf >= 0);
    }

    free(node_map);
    free(inv_node_map);

    free(gf_nodes);

    free(f);
    free(f_old);
    free(residual);

    return;
}
/*----------------------------------------------------------------------------*/


