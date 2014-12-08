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

#define TOL 0

#define MAT_VERSION_STRING "1.0 " __DATE__ " " __TIME__

void calculate_stress(job_t *job);
void calculate_bulk_granular_fluidity(job_t *job);
void solve_diffusion_part(job_t *job);
void calculate_stress_threaded(threadtask_t *task);

static double E, nu, G, K;
static double mu_s, mu_2, I_0, rho_s, rho_c, d, A;
static double *Kp; //vector which is a solution of the matrix-vector product.
static double *g; //nodal values of g

void cs_print_to_file(const cs *A)
{

    int p, j, m, n, nzmax, nz, *Ap, *Ai ;
    double *Ax ;
    if (!A) { printf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    
    double *A_dense = calloc(sizeof(double), m*n);
    FILE *fp = fopen("matrix.cs", "w");

    for (j = 0 ; j < n ; j++)
    {
        printf ("    col %g : locations %g to %g\n", (double) j, 
            (double) (Ap [j]), (double) (Ap [j+1]-1)) ;
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            A_dense[Ai[p]*n + j] = Ax[p];
            printf ("      %g : %g\n", (double) (Ai [p]), Ax ? Ax [p] : 1) ;
        }
    }

    for (int i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
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
    printf("%s:%s: Done allocating storage for nodal values of g.\n",
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
    /* value at end of timestep */
    double tau_tau;
    double scale_factor;

    /* increment using jaumann rate */
    double dsjxx, dsjxy, dsjyy;

    /* trial values */
    double sxx_tr, sxy_tr, syy_tr;
    double t0xx_tr, t0xy_tr, t0yy_tr;
    double p_tr, tau_tr;

    double nup_tau;

    double const c = 0;

    /* solve for g_local */
    calculate_bulk_granular_fluidity(job);

    /* build FEM diffusion array/load vector and solve for g_nonlocal */
    solve_diffusion_part(job);

    for (size_t i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }
        const double lambda = K - 2.0 * G / 3.0;

        /* Calculate tau and p trial values. */
        const double trD = job->particles[i].exx_t + job->particles[i].eyy_t;
        dsjxx = lambda * trD + 2.0 * G * job->particles[i].exx_t;
        dsjxy = 2.0 * G * job->particles[i].exy_t;
        dsjyy = lambda * trD + 2.0 * G * job->particles[i].eyy_t;
        dsjxx += 2 * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy -= job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy -= 2 * job->particles[i].wxy_t * job->particles[i].sxy;

        sxx_tr = job->particles[i].sxx + job->dt * dsjxx;
        sxy_tr = job->particles[i].sxy + job->dt * dsjxy;
        syy_tr = job->particles[i].syy + job->dt * dsjyy;

        p_tr = -0.5 * (sxx_tr + syy_tr);
        t0xx_tr = sxx_tr + p_tr;
        t0xy_tr = sxy_tr;
        t0yy_tr = syy_tr + p_tr;
        tau_tr = sqrt(0.5*(t0xx_tr*t0xx_tr + 2*t0xy_tr*t0xy_tr + t0yy_tr*t0yy_tr));

#if 0
        if ((job->particles[i].m / job->particles[i].v) < rho_c) {
            dense = 0;
/*            printf("%4d: density %lf\n", i, (job->particles[i].m / job->particles[i].v));*/
        } else {
            dense = 1;
        }
#endif

        if (dense == 0 || p_tr <= c) {
            nup_tau = (tau_tr / G) / job->dt;
            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
        } else if (p_tr > c) {
            scale_factor = p_tr / (G * job->dt * gf + p_tr);
            tau_tau = tau_tr * scale_factor;

            nup_tau = ((tau_tr - tau_tau) / G) / job->dt;

            job->particles[i].sxx = scale_factor * t0xx_tr - p_tr;
            job->particles[i].sxy = scale_factor * t0xy_tr;
            job->particles[i].syy = scale_factor * t0yy_tr - p_tr;
        } else {
/*            fprintf(stderr, "u %zu %3.3g %3.3g %d ", i, f, p_tr, density_flag);*/
            fprintf(stderr, "u"); 
            nup_tau = 0;
        }

        /* use strain rate to calculate stress increment */
        gammap += nup_tau * job->dt;
        gammadotp = nup_tau;
    }

    return;
}
/*----------------------------------------------------------------------------*/

double negative_root(double a, double b, double c);

double negative_root(double a, double b, double c)
{
    double x;
    if (b > 0) {
        x = (-b - sqrt(b*b - 4*a*c)) / (2*a);
    } else {
        x = (2*c) / (-b + sqrt(b*b - 4*a*c));
    }
    return x;
}

/*----------------------------------------------------------------------------*/
void calculate_bulk_granular_fluidity(job_t *job)
{
    /* value at end of timestep */
    double tau_tau;
    double scale_factor;

    /* increment using jaumann rate */
    double dsjxx, dsjxy, dsjyy;

    /* trial values */
    double sxx_tr, sxy_tr, syy_tr;
    double t0xx_tr, t0xy_tr, t0yy_tr;
    double p_tr, tau_tr;

    double const c = 0;
   
    size_t i = 0;
 
    double trD;
    const double lambda = K - 2.0 * G / 3.0;

    double S0, S2;
    double B, H;
    double alpha;

    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        /* Calculate tau and p trial values. */
        trD = job->particles[i].exx_t + job->particles[i].eyy_t;
        dsjxx = lambda * trD + 2.0 * G * job->particles[i].exx_t;
        dsjxy = 2.0 * G * job->particles[i].exy_t;
        dsjyy = lambda * trD + 2.0 * G * job->particles[i].eyy_t;
        dsjxx += 2 * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy -= job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy -= 2 * job->particles[i].wxy_t * job->particles[i].sxy;

        sxx_tr = job->particles[i].sxx + job->dt * dsjxx;
        sxy_tr = job->particles[i].sxy + job->dt * dsjxy;
        syy_tr = job->particles[i].syy + job->dt * dsjyy;

        p_tr = -0.5 * (sxx_tr + syy_tr);
        t0xx_tr = sxx_tr + p_tr;
        t0xy_tr = sxy_tr;
        t0yy_tr = syy_tr + p_tr;
        tau_tr = sqrt(0.5*(t0xx_tr*t0xx_tr + 2*t0xy_tr*t0xy_tr + t0yy_tr*t0yy_tr));

        if ((job->particles[i].m / job->particles[i].v) < rho_c || p_tr < 0) {
            dense = 0;
/*            printf("%4d: density %lf\n", i, (job->particles[i].m / job->particles[i].v));*/
        } else {
            dense = 1;
        }

        const double p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);
        const double t0xx_t = job->particles[i].sxx + p_tr;
        const double t0xy_t = job->particles[i].sxy;
        const double t0yy_t = job->particles[i].syy + p_tr;
        const double tau_t = sqrt(0.5*(t0xx_t*t0xx_t + 2*t0xy_t*t0xy_t + t0yy_t*t0yy_t));

        if (dense == 0 || p_t <= c) {
            gf_local = 0;
            xisq_inv = 0;
            dense = 0;
        } else if (p_t > c) {
            S0 = mu_s * p_tr;
            if (tau_tr <= S0) {
                tau_tau = tau_tr;
                scale_factor = 1.0;
            } else {
                S2 = mu_2 * p_tr;
                alpha = G * I_0 * job->dt * sqrt(p_tr / rho_s) / d;
                B = -(S2 + tau_tr + alpha);
                H = S2 * tau_tr + S0 * alpha;
                tau_tau = negative_root(1.0, B, H);
                scale_factor = (tau_tau / tau_tr);
            }
            // gf_local = (p_tr / (G * job->dt)) * ((1.0 / scale_factor) - 1.0);
            // xisq_inv = fabs((tau_tau / p_tr) - mu_s) * (mu_2 - mu_s) / (fabs(mu_2 - (tau_tau / p_tr)) * (A * A * d * d));
            double s = 0;
            if (tau_t > S0) {
                s = (tau_t - S0) / (tau_t * (tau_t - S2));
            }
            const double zeta = I_0 / (d * sqrt(rho_s));
            gf_local = p_t * sqrt(p_t) * zeta * s;
            xisq_inv = fabs((tau_t/ p_t) - mu_s) * (mu_2 - mu_s) / (fabs(mu_2 - (tau_t / p_t)) * (A * A * d * d));
            assert(gf_local >= 0);
            assert(xisq_inv >= 0);
            dense = 1;
        } else {
/*            fprintf(stderr, "u %zu %3.3g %3.3g %d ", i, f, p_tr, density_flag);*/
            fprintf(stderr, "u");
            gf_local = 0;
            xisq_inv = 0;
            dense = 0;
        }
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
    double *m_total_nodes;

    double *f;

    double s[4];
    double grad_s[4][2];

/*    double m_frac;*/

    int p;
    int *nn;

    /* initialize node level fluidity arrays */
    gf_nodes = (double *)malloc(job->num_nodes * sizeof(double));
    for (i = 0; i < job->num_nodes; i++) {
        gf_nodes[i] = 0;
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

        if (job->nodes[i_new].m > TOL) {
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

    m_total_nodes = (double *)malloc(slda * sizeof(double));
    f = (double *)malloc(slda * sizeof(double));
    for (i = 0; i < slda; i++) {
        m_total_nodes[i] = 0;
        f[i] = 0;
    }

    /*  calculate number of nonzero elements (before summing duplicates). */
    nnz = 0;
    for (i = 0; i < job->num_particles; i++) {
        if (dense) {
            nnz += (NODES_PER_ELEMENT * NODES_PER_ELEMENT) * 1;
        }
    }
    nnz += slda;
    /*
    for (i = 0; i < job->num_elements; i++) {
        if (job->elements[i].filled != 0) {
            nnz += (NODES_PER_ELEMENT * NODES_PER_ELEMENT) * 1;
        }
    }
    */

    /* slda contains degrees of freedom of new matrix */
    triplets = cs_spalloc(slda, slda, nnz, 1, 1);
/*    fprintf(stderr, "%s:%s: slda = %zu, nnz = %zu.\n",*/
/*        __FILE__, __func__, slda, nnz);*/

    /* project xi and g_local onto background grid (volume-weighted) */
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

        nn = job->elements[p].nodes;

        assert(isfinite(gf_local));
        assert(isfinite(xisq_inv));
        
        for (ei = 0; ei < NODES_PER_ELEMENT; ei++) {
            gi = nn[ei];
            gi = (job->node_number_override[NODAL_DOF * gi + 0 ] - 0) / NODAL_DOF;
            sgi = node_map[gi];
            
            assert(sgi != -1);

            m_total_nodes[sgi] += job->particles[i].m * s[ei];
            f[sgi] += (job->particles[i].v) * gf_local * xisq_inv * s[ei];
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

    /* keep matrix from being dengerate when an element is filled with open particles. */
    for (i = 0; i < slda; i++) {
        cs_entry(triplets, i, i, 1e-10);
    }

    /* create compressed sparse matrix */
    smat = cs_compress(triplets);
    cs_dupl(smat);

    if (!cs_lusol(1, smat, f, 1e-12)) {
        fprintf(stderr, "lusol error!\n");
        if (cs_qrsol(1, smat, f)) {
            fprintf(stderr, "qrsol error!\n");
            cs_print(smat, 0); // print out matrix for debugging
            cs_print_to_file(smat);
            exit(EXIT_ERROR_CS_SOL);
        }
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

            gf += gf_nodes[gi] * s[ei];
        }

        if (gf < 0) {
            fprintf(stderr, "%d: %lg\n", i, gf);
        }
        assert(gf >= 0);
    }

    cs_spfree(triplets);
    cs_spfree(smat);

    free(node_map);
    free(inv_node_map);

    free(gf_nodes);

    free(m_total_nodes);
    free(f);

    return;
}
/*----------------------------------------------------------------------------*/


