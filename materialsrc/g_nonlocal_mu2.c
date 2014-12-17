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

#define jp(x) job->particles[i].x

#undef EMOD
#undef NUMOD

#define dense jp(state[0])
#define gf jp(state[3])
#define eta jp(state[4])
#define gf_local jp(state[5])
#define xisq_inv jp(state[6])
#define gammap jp(state[9])
#define gammadotp jp(state[10])

#define MAT_VERSION_STRING "1.0 " __DATE__ " " __TIME__

void calculate_stress(job_t *job);
void solve_diffusion_part(job_t *job);
void calculate_stress_threaded(threadtask_t *task);

static double E, nu, G, K;
static double mu_s, mu_2, I_0, rho_s, rho_c, d, A;
static double *Kp; //vector which is a solution of the matrix-vector product.
static double *g_loc; //nodal values of g_local
static double *g; //nodal values of g
static double *v; //nodal volumes (nodal mass already computed)
static double *tau_tr; // nodal values of tau_tr
static double *p_tr; // nodal values of p_tr

double calculate_g_local(double tau, double p);
double calculate_xisq_inverse(double tau, double p);
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

double calculate_xisq_inverse(double tau, double p)
{
    double xisq_inverse = 0;
    if (p > 0) {
        const double S0 = mu_s * p;
        const double S2 = mu_2 * p;

        if (tau >= S2) {
            fprintf(stderr, "\n%g > %g: %g\n", tau, S2, tau / p);
        }
        assert(tau < S2);

        if (tau < S2) {
            xisq_inverse = fabs(tau - S0) * (mu_2 - mu_s) / (fabs(S2 - tau) * A * A * d * d);
        } else {
            xisq_inverse = 1;
        }
    }
    return xisq_inverse;
}

double calculate_g_local_from_g(double tau_tr, double p_tr, double g, double delta_t)
{
    double g_local = 0;
    if (g >= 0 && tau_tr > 0 && p_tr > 0) {
        const double zeta = I_0 / (d * sqrt(rho_s));
        const double s = g * G * delta_t + p_tr;
        g_local = (sqrt(p_tr) / tau_tr) * zeta * s * (tau_tr - mu_s * s) / (mu_2 * s - tau_tr);

    }
    return g_local;
}

typedef struct p_trial_s {
    double tau_tr;
    double p_tr;

    double t0xx_tr;
    double t0xy_tr;
    double t0yy_tr;

    double tau_tau;
    double p_tau;

    double s; // scaling factor
} trial_t;

void trial_step(const particle_t *p, const double dt, trial_t *trial)
{
    const double lambda = K - 2.0 * G / 3.0;

    const double trD = p->exx_t + p->eyy_t;
    const double dsjxx = lambda * trD + 2.0 * G * p->exx_t + 2 * p->wxy_t * p->sxy;
    const double dsjxy = 2.0 * G * p->exy_t - p->wxy_t * (p->sxx - p->syy);
    const double dsjyy = lambda * trD + 2.0 * G * p->eyy_t - 2 * p->wxy_t * p->sxy;

    const double sxx_tr = p->sxx + dt * dsjxx;
    const double sxy_tr = p->sxy + dt * dsjxy;
    const double syy_tr = p->syy + dt * dsjyy;

    const double p_tr = -0.5 * (sxx_tr + syy_tr);
    const double t0xx_tr = sxx_tr + p_tr;
    const double t0xy_tr = sxy_tr;
    const double t0yy_tr = syy_tr + p_tr;
    const double tau_tr = sqrt(0.5*(t0xx_tr*t0xx_tr + 2*t0xy_tr*t0xy_tr + t0yy_tr*t0yy_tr));

    trial->tau_tr = tau_tr;
    trial->p_tr = p_tr;
    trial->t0xx_tr = t0xx_tr;
    trial->t0xy_tr = t0xy_tr;
    trial->t0yy_tr = t0yy_tr;

    return;
}

// assume an active, valid particle is passed in to p
// retuns dot(bar(gamma))^p
//double local_step(particle_t *p, double dt, double * restrict tau_tau, double * restrict p_tau, int * restrict flag);
double local_step(const particle_t *p, double dt, trial_t *trial, int *flag);

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

    int m, n, nzmax, nz, *Ap, *Ai ;
    double *Ax ;
    if (!A) { return ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    
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

double dot(const double * a, const double * b, size_t n);
double dot(const double * a, const double * b, size_t n)
{
    double s = 0;
    for (size_t i = 0; i < n; i++) {
        s += a[i] * b[i];
    }
    return s;
}

int cs_cg(const cs *K, double *f, const double *u_0, double tol);
int cs_cg(const cs *K, double *f, const double *u_0, double tol)
{
    int converged = 0;
    if (K->m != K->n) {
        return converged; //nonsquare matrix
    }
    const size_t n = K->n;

    double * restrict u = malloc(sizeof(double) * n);
    double * restrict r = malloc(sizeof(double) * n);
    double * restrict p = malloc(sizeof(double) * n);
    const double rsq_tol = tol*tol;

    if (u_0 != NULL) {
        memcpy(u, u_0, sizeof(double) * n);
        for (size_t j = 0; j < n; j++) {
            r[j] = 0;
        }
        cs_gaxpy(K, u, r); // r is now K * u_0
        for (size_t j = 0; j < n; j++) {
            r[j] = f[j] - r[j];
        }
    } else {
        for (size_t j = 0; j < n; j++) {
            u[j] = 0;
        }
        memcpy(r, f, sizeof(double) * n);
    }
    memcpy(p, r, sizeof(double) * n);

    double * restrict Kp = malloc(sizeof(double) * n);
    for (size_t i = 0; i < n; i++) {
        const double rTr = dot(r, r, n);
        // printf("r.r = %lg\n", rTr);
        if (rTr <= rsq_tol) {
            converged = 1;
            // printf("\n");
            break;
        }
        for (size_t j = 0; j < n; j++) {
            Kp[j] = 0;
        }
        cs_gaxpy(K, p, Kp); // Kp is now K * p_{i-1}
        const double pTKp = dot(p, Kp, n);
        assert(pTKp != 0);
        const double alpha = rTr / pTKp;
        for (size_t j = 0; j < n; j++) {
            u[j] = u[j] + alpha * p[j];
            r[j] = r[j] - alpha * Kp[j];
        }
        const double rTr_new = dot(r, r, n);
        const double beta = rTr_new / rTr;
        for (size_t j = 0; j < n; j++) {
            p[j] = r[j] + beta * p[j];
        }
    }

    memcpy(f, u, sizeof(double) * n);
    free(u);
    free(r);
    free(p);
    return converged;
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
        xisq_inv = 0;
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
        printf("%s:%s: properties (E = %g, nu = %g, G = %g, K = %g, mu_s = %g,"
            " mu_2 = %g, I_0 = %g, rho_s = %g, rho_c = %g, d = %g, A = %g).\n",
            __FILE__, __func__, E, nu, G, K, mu_s, mu_2, I_0, rho_s, rho_c, d, A);
    }

    Kp = calloc(job->num_nodes, sizeof(double));
    printf("%s:%s: Done allocating storage for matrix-vector multiply.\n",
        __FILE__, __func__);
    g = calloc(job->num_nodes, sizeof(double));
    g_loc = calloc(job->num_nodes, sizeof(double));
    tau_tr = calloc(job->num_nodes, sizeof(double));
    p_tr = calloc(job->num_nodes, sizeof(double));
    printf("%s:%s: Done allocating storage for nodal values of g.\n",
        __FILE__, __func__);
    v = calloc(job->num_nodes, sizeof(double));
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
        const double tau_tr = tr.tau_tr;

        if (dense) {
            const double scale_factor = p_tr / (G * job->dt * gf + p_tr);
            const double tau_tau = tau_tr * scale_factor;

            nup_tau = ((tau_tr - tau_tau) / G) / job->dt;

            job->particles[i].sxx = scale_factor * t0xx_tr - p_tr;
            job->particles[i].sxy = scale_factor * t0xy_tr;
            job->particles[i].syy = scale_factor * t0yy_tr - p_tr;
        } else {
            nup_tau = (tau_tr / G) / job->dt;
            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
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
    size_t i, j;
    size_t i_new, j_new;
    size_t ei, ej;
    size_t gi, gj;
    size_t sgi, sgj;

    size_t *node_map;
    size_t *inv_node_map;

    size_t slda;
    size_t nnz;

    cs *triplets;
    cs *smat;

    double *gf_nodes;

    double *f;

    double s[4];
    double grad_s[4][2];

    int p;
    int *nn;

    /* initialize node level fluidity arrays */
    gf_nodes = (double *)malloc(job->num_nodes * sizeof(double));
    for (i = 0; i < job->num_nodes; i++) {
        gf_nodes[i] = 0;
        g[i] = 0;
        g_loc[i] = 0;
        tau_tr[i] = 0;
        p_tr[i] = 0;
        // clear volumes
        v[i] = 0;
    }

    /* get number of dofs and initialize mapping arrays */
    node_map = (size_t *)malloc(job->num_nodes * sizeof(size_t));
    inv_node_map = (size_t *)malloc(job->num_nodes * sizeof(size_t));
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

    /*  calculate number of nonzero elements (before summing duplicates). */
    nnz = 0;
    for (i = 0; i < job->num_particles; i++) {
        if (dense) {
            nnz += (NODES_PER_ELEMENT * NODES_PER_ELEMENT) * 1;
        }
    }
    nnz += slda;

    /* slda contains degrees of freedom of new matrix */
    triplets = cs_spalloc(slda, slda, nnz, 1, 1);

    /* project xi and g_local onto background grid (volume-weighted) */
    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        p = job->in_element[i];
        if (p == -1) {
            continue;
        }

        trial_t tr;
        int flag = 0;
        double nu_p = local_step(&(job->particles[i]), job->dt, &tr, &flag);
        dense = flag;
        // printf("%g, %g, %g\n", tr.tau_tau, tr.p_tau, tr.tau_tau / tr.p_tau);
        if (flag == 1) {
            gf_local = tr.tau_tr * ((1.0 / tr.s) - 1.0) / (G * job->dt);
            xisq_inv = calculate_xisq_inverse(tr.tau_tau, tr.p_tau);
        } else {
            gf_local = 0;
            xisq_inv = 0;
            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
            continue; // don't add contribution to stiffness matrix
        }
        assert(gf_local >= 0);
        assert(xisq_inv >= 0);

        if (dense == 0) {
            continue;
        }

        assert(tr.tau_tau >= 0 && tr.tau_tau <= tr.p_tau * mu_2);

        s[0] = job->h1[i];
        s[1] = job->h2[i];
        s[2] = job->h3[i];
        s[3] = job->h4[i];

        nn = job->elements[p].nodes;

        assert(isfinite(gf_local));
        assert(isfinite(xisq_inv));
        
        for (ei = 0; ei < NODES_PER_ELEMENT; ei++) {
            gi = nn[ei];
            gi = (job->node_number_override[NODAL_DOF * gi + 0 ] - 0) / NODAL_DOF;
            sgi = node_map[gi];
            
            assert(sgi != -1);

            const double f_component = (job->particles[i].v) * gf_local * xisq_inv * s[ei];
            const double gf_component = (job->particles[i].v) * gf_local * s[ei];
            if (f_component < 0) {
                printf("v: %lg\nf: %lg\n", job->particles[i].v, f_component);
            }
            assert(f_component >= 0);
            assert(gf_component >= 0);


            f[sgi] += f_component;
            gf_nodes[sgi] += gf_component;
            v[sgi] += s[ei] * job->particles[i].v;
        }
    }
    
    /* create stiffness matrix. */
    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }
        
        if (dense == 0) {
            continue; // don't add stiffness contribution if not dense.
        }

        p = job->in_element[i];
        if (p == -1) {
            continue;
        }

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

        nn = job->elements[p].nodes;

        for (ei = 0; ei < NODES_PER_ELEMENT; ei++) {
            for (ej = 0; ej < NODES_PER_ELEMENT; ej++) {
                gi = nn[ei];
                gj = nn[ej];

                gi = (job->node_number_override[NODAL_DOF * gi + 0] - 0) / NODAL_DOF;
                gj = (job->node_number_override[NODAL_DOF * gj + 0] - 0) / NODAL_DOF;

                sgi = node_map[gi];
                sgj = node_map[gj];

                assert(isfinite(gf_local));
                assert(isfinite(xisq_inv));
                assert(isfinite(job->particles[i].v));

                const double k_component = job->particles[i].v * (xisq_inv * s[ei] * s[ej] + 
                        (grad_s[ei][0]*grad_s[ej][0] + grad_s[ei][1]*grad_s[ej][1]));

                assert(isfinite(k_component));
                
                cs_entry(triplets, sgi, sgj, k_component);
            }
        }
    }

    /* figure out local values of g at nodes */
    for (i = 0; i < slda; i++) {
        if (v[i] > 0) {
            gf_nodes[i] = gf_nodes[i] / v[i];
        } else {
            gf_nodes[i] = 0;
        }
        g_loc[i] = gf_nodes[i];
        assert(gf_nodes[i] >= 0);
    }

    for (i = 0; i < slda; i++) {
        p_tr[i] = f[i]; // copy load;
    }

#define DIRECT_SOLVE
#ifndef DIRECT_SOLVE
    /* create compressed sparse matrix */
    smat = cs_compress(triplets);
    cs_dupl(smat);
//    fprintf(stderr, "%d by %d\n", smat->m, smat->n); // print out matrix for debugging
    if (!cs_cg(smat, f, gf_nodes, 1e-15)) {
        fprintf(stderr, "cg error!\n");
        cs_print(smat, 0); // print out matrix for debugging
        cs_print_to_file(smat);
        exit(EXIT_ERROR_CS_SOL);
    }
#else
    /* keep matrix from being dengerate when an element is filled with open particles. */
    for (i = 0; i < slda; i++) {
        cs_entry(triplets, i, i, 1e-10);
    }
    /* create compressed sparse matrix */
    smat = cs_compress(triplets);
    cs_dupl(smat);
//    fprintf(stderr, "%d by %d\n", smat->m, smat->n); // print out matrix for debugging

    if (!cs_lusol(1, smat, f, 1e-12)) {
        fprintf(stderr, "lusol error!\n");
        if (cs_qrsol(1, smat, f)) {
            fprintf(stderr, "qrsol error!\n");
            cs_print(smat, 0); // print out matrix for debugging
            cs_print_to_file(smat);
            exit(EXIT_ERROR_CS_SOL);
        }
    }
#endif

    for (i = 0; i < job->num_nodes; i++) {
        gf_nodes[i] = 0; // clear guess values
    }

    for (i = 0; i < job->num_nodes; i++) {
        i_new = (job->node_number_override[NODAL_DOF * i + 0] - 0) / NODAL_DOF;

        if (node_map[i_new] != -1) {
            gf_nodes[i] = f[node_map[i_new]];
        }
    }

    /* ensure periodic BCs ok */
    for (i = 0; i < job->num_nodes; i++) {
        i_new = (job->node_number_override[NODAL_DOF * i + 0] - 0) / NODAL_DOF;

        if (i != i_new) {
            assert(gf_nodes[i] == gf_nodes[i_new]);
        }
    }

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
            gf += gf_nodes[gi] * s[ei];
        }

        if (gf < 0) {
            fprintf(stderr, "%d: %lg\n", i, gf);
            cs_print_to_file(smat);
            FILE *ff = fopen("load.cs", "w");
            for (size_t j = 0; j < slda; j++) {
                fprintf(ff, "%lg\n", p_tr[j]);
            }
            fclose(ff);
        }
        assert(gf >= 0);
    }

    cs_spfree(triplets);
    cs_spfree(smat);

    free(node_map);
    free(inv_node_map);

    free(gf_nodes);

    free(f);

    return;
}
/*----------------------------------------------------------------------------*/


