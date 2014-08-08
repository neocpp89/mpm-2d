/**
    \file g_local.c
    \author Sachith Dunatunga
    \date 23.10.13

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

#define signum(x) ((int)((0 < x) - (x < 0)))

#define jp(x) job->particles[i].x

#undef EMOD
#undef NUMOD

#define EMOD 1e8
#define NUMOD 0.3

#define G (EMOD / (2.0f * (1.0f + NUMOD)))
#define K (EMOD / (3.0f * (1.0f - 2*NUMOD)))

#undef G
#undef K

#define PI 3.1415926535897932384626433
#define PHI (30.0*PI/180.0)

/* from Dave + Ken's paper (modified B) */
#define MU_S 0.3819
#define GRAINS_RHO 2450
#define B (0.9377 * 20)
#define A 0.48

/* from Jop (modified I_0) */
#define MU_2 0.6435
#define I_0 (0.278 / 1.0)

/* from geometric considerations */
#define GRAINS_D (0.01 * 5)


#define Epxx jp(state[0])
#define Epxy jp(state[1])
#define Epyy jp(state[2])
#define gf jp(state[3])
#define eta jp(state[4])
#define gf_local jp(state[5])
#define xisq_inv jp(state[6])
/*#define Etxy jp(state[7])*/
/*#define Etyy jp(state[8])*/
#define gammap jp(state[9])
#define gammadotp jp(state[10])

#define TOL 1e-10

#define MAT_VERSION_STRING "1.11 " __DATE__ " " __TIME__

void calculate_stress(job_t *job);
void calculate_bulk_granular_fluidity(job_t *job);
void solve_diffusion_part(job_t *job);
void calculate_stress_threaded(threadtask_t *task);

static double E, nu, G, K;

/*----------------------------------------------------------------------------*/
void material_init(job_t *job)
{
    int i, j;

    for (i = 0; i < job->num_particles; i++) {
        for (j = 0; j < DEPVAR; j++) {
            job->particles[i].state[j] = 0;
        }
    }

    for (i = 0; i < job->num_particles; i++) {
        Epxx = 0;
        Epxy = 0;
        Epyy = 0;
/*        Etxx = 0;*/
/*        Etxy = 0;*/
/*        Etyy = 0;*/
        eta = 0;
        gammap = 0;
        gammadotp = 0;
        gf = 0;
        gf_local = 0;
        xisq_inv = 0;
    }

    if (job->material.num_fp64_props < 2) {
        fprintf(stderr,
            "%s:%s: Need at least 2 properties defined (E, nu).\n",
            __FILE__, __func__);
        exit(EXIT_ERROR_MATERIAL_FILE);
    } else {
        E = job->material.fp64_props[0];
        nu = job->material.fp64_props[1];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2*nu));
        printf("%s:%s: properties (E = %g, nu = %g, G = %g, K = %g).\n",
            __FILE__, __func__, E, nu, G , K);
    }

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
    /* values from previous timestep */
    double p_t, mu_t;

    /* value at end of timestep */
    double tau_tau;
    double scale_factor;

    /* increment using jaumann rate */
    double dsjxx, dsjxy, dsjyy;

    /* trial values */
    double sxx_tr, sxy_tr, syy_tr;
    double t0xx_tr, t0xy_tr, t0yy_tr;
    double p_tr, tau_tr;

    double beta, nup_tau;

    double const c = 1e-2;

    double f;
    int density_flag;
    
    size_t i;

    double inertial_num;

/*    fprintf(stderr, "processing particle ids [%zu %zu].\n", p_start, p_stop);*/

    /* solve for g_local */
    calculate_bulk_granular_fluidity(job);

    /* build FEM diffusion array/load vector and solve for g_nonlocal */
    solve_diffusion_part(job);

    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        /* Calculate p at beginning of timestep. */
        p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);

        /* Calculate tau and p trial values. */
        dsjxx = job->dt * (E / (1 - nu*nu)) * ((job->particles[i].exx_t) + nu * (job->particles[i].eyy_t));
        dsjxy = job->dt * (E / (2 *(1 + nu))) * (job->particles[i].exy_t);
        dsjyy = job->dt * (E / (1 - nu*nu)) * ((job->particles[i].eyy_t) + nu * (job->particles[i].exx_t));
        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

        sxx_tr = job->particles[i].sxx + dsjxx;
        sxy_tr = job->particles[i].sxy + dsjxy;
        syy_tr = job->particles[i].syy + dsjyy;

        p_tr = -0.5 * (sxx_tr + syy_tr);
        t0xx_tr = sxx_tr + p_tr;
        t0xy_tr = sxy_tr;
        t0yy_tr = syy_tr + p_tr;
        tau_tr = sqrt(0.5*(t0xx_tr*t0xx_tr + 2*t0xy_tr*t0xy_tr + t0yy_tr*t0yy_tr));

        inertial_num = gammadotp * GRAINS_D * sqrt(GRAINS_RHO);
/*        if (p_t <= 0) {*/
/*            mu_t = MU_2;*/
/*        } else {*/
/*            inertial_num = inertial_num / sqrt(p_t);*/
/*            if (inertial_num <= 0) {*/
/*                mu_t = MU_S;*/
/*            } else {*/
/*                mu_t = MU_S + (MU_2 - MU_S) / ((I_0 / inertial_num) + 1.0);*/
/*            }*/
/*        }*/

        tau_tau = mu_t*(p_tr + c);
        f = tau_tr - tau_tau;

        if ((job->particles[i].m / job->particles[i].v) < 1485.0f) {
            density_flag = 1;
/*            printf("%4d: density %lf\n", i, (job->particles[i].m / job->particles[i].v));*/
        } else {
            density_flag = 0;
        }

        if (density_flag) {
            nup_tau = (tau_tr) / (G * job->dt);
            beta = -p_tr / (K * job->dt);

            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
        } else if (f < 0) {
            nup_tau = 0;
            beta = 0;

            job->particles[i].sxx = t0xx_tr - p_tr;
            job->particles[i].sxy = t0xy_tr;
            job->particles[i].syy = t0yy_tr - p_tr;
        } else if (f >= 0 && p_tr > -c/mu_t) {
            nup_tau = (tau_tr - tau_tau) / (G * job->dt);
            beta = 0;

            scale_factor = (tau_tau / tau_tr);
            job->particles[i].sxx = scale_factor * t0xx_tr - p_tr;
            job->particles[i].sxy = scale_factor * t0xy_tr;
            job->particles[i].syy = scale_factor * t0yy_tr - p_tr;
        } else if (p_tr <= -c/mu_t) {
            nup_tau = tau_tr / (G * job->dt);
            beta = ((c/mu_t) - p_tr) / (K * job->dt);

            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
        } else {
            fprintf(stderr, "u %zu %3.3g %3.3g %d ", i, f, p_tr, density_flag);
/*            fprintf(stderr, "u"); */
            nup_tau = 0;
        }

        /* use strain rate to calculate stress increment */
        gammap += nup_tau * job->dt;
        gammadotp = nup_tau;
    }

    return;
}
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
void calculate_bulk_granular_fluidity(job_t *job)
{
    size_t i;
    double p_t;
    double tau_t;
    double S0;

    double sxx_tr0;
    double sxy_tr0;
    double syy_tr0;

    double const p_min = 1e-2;

    for (i = 0; i < job->num_particles; i++) {

        /* Calculate pressure at beginning of step. */
        p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);

        sxx_tr0 = job->particles[i].sxx + p_t;
        sxy_tr0 = job->particles[i].sxy;
        syy_tr0 = job->particles[i].syy + p_t;

        /* Calculate equivalent shear at beginning of step. */
        tau_t = sqrt(0.5*(sxx_tr0*sxx_tr0 + 2*sxy_tr0*sxy_tr0 + syy_tr0*syy_tr0));

        if (p_t < p_min) {
/*            job->particles[i].sxx = -0.5 * p_cap;*/
/*            job->particles[i].sxy = 0;*/
/*            job->particles[i].syy = -0.5 * p_cap;*/
            p_t = p_min;
/*            printf("%s:%s: Particle %zu pressure less than cohesion, results may be inaccurate.\n", __FILE__, __func__, i);*/
/*            continue;*/
        }

        S0 = p_t * MU_S;
        eta = B * GRAINS_D * sqrt(p_t * GRAINS_RHO);
        xisq_inv = fabs(tau_t / p_t - MU_S)/ (GRAINS_D * GRAINS_D * A * A);

        if (tau_t > S0) { 
            gf_local = (p_t / eta) * (1 - S0 / tau_t);
        } else {
            gf_local = 0;
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

    if (slda <= 0) {
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
    for (i = 0; i < job->num_elements; i++) {
        if (job->elements[i].filled != 0) {
            nnz += (NODES_PER_ELEMENT * NODES_PER_ELEMENT) * 1;
        }
    }

    /* slda contains degrees of freedom of new matrix */
    triplets = cs_spalloc(slda, slda, nnz, 1, 1);
/*    fprintf(stderr, "%s:%s: slda = %zu, nnz = %zu.\n",*/
/*        __FILE__, __func__, slda, nnz);*/

    /* project xi and g_local onto background grid (volume-weighted) */
    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
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

        nn = job->elements[p].nodes;

        for (ei = 0; ei < NODES_PER_ELEMENT; ei++) {
            gi = nn[ei];
            gi = (job->node_number_override[NODAL_DOF * gi + 0 ] - 0) / NODAL_DOF;
            sgi = node_map[gi];

            m_total_nodes[sgi] += job->particles[i].m * s[ei];
/*            f[sgi] += job->particles[i].m * gf_local * xisq_inv * s[ei];*/
            f[sgi] += (job->particles[i].v) * gf_local * xisq_inv * s[ei];
        }
    }

    /* unweight load f */
/*    for (i = 0; i < slda; i++) {*/
/*        f[i] = f[i] / m_total_nodes[i];*/
/*    }*/

    /* create stiffness matrix. */
    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
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

/*                m_frac = (job->particles[i].m * job->particles[i].m) / (m_total_nodes[sgi] * m_total_nodes[sgj]);*/

                cs_entry(triplets, sgi, sgj,
                    job->particles[i].v * (xisq_inv * s[ei] * s[ej] + 
                        (grad_s[ei][0]*grad_s[ej][0] + grad_s[ei][1]*grad_s[ej][1]))
                );
            }
        }
    }

    /* create compressed sparse matrix */
    smat = cs_compress(triplets);
    cs_dupl(smat);

    if (!cs_lusol(1, smat, f, 1e-12)) {
        fprintf(stderr, "lusol error!\n");
        if (cs_qrsol(1, smat, f)) {
            fprintf(stderr, "qrsol error!\n");
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


