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

#define PI 3.1415926535897932384626433
#define PHI (30.0*PI/180.0)

/* from Dave + Ken's paper (modified B) */
#define MU_S 0.3819
#define GRAINS_RHO 2450
#define B (0.9377 * 20)
#define A 0.48

/* from geometric considerations */
#define GRAINS_D (0.01 * 5)


#define Epxx jp(state[0])
#define Epxy jp(state[1])
#define Epyy jp(state[2])
#define gf jp(state[3])
#define eta jp(state[4])
#define gf_bulk jp(state[5])
#define gammap jp(state[9])
#define xisq_inv jp(state[10])
#define Etxx jp(state[6])
#define Etxy jp(state[7])
#define Etyy jp(state[8])

#define TOL 1e-10

#define MAT_VERSION_STRING "1.11 " __DATE__ " " __TIME__

void calculate_bulk_granular_fluidity(job_t *job);
void solve_diffusion_part(job_t *job);

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
        Etxx = 0;
        Etxy = 0;
        Etyy = 0;
        eta = 0;
        gammap = 0;
        gf = 0;
        gf_bulk = 0;
        xisq_inv = 0;
    }

    printf("%s:%s: (material version %s) done initializing material.\n",
        __FILE__,  __func__, MAT_VERSION_STRING);
    return;
}
/*----------------------------------------------------------------------------*/

/* nonlocal granular fluidity model. */
/*----------------------------------------------------------------------------*/
void calculate_stress(job_t *job)
{
    size_t i;
    double p_t;
    double p_tr;
    double tau_t;
    double tau_tr;
    double tau_tau;

    double dsjxx;
    double dsjxy;
    double dsjyy;

    double sxx_tr;
    double sxy_tr;
    double syy_tr;

    double sxx_tr0;
    double sxy_tr0;
    double syy_tr0;

    double sxx_t0;
    double sxy_t0;
    double syy_t0;

    double Npxx_tau;
    double Npxy_tau;
    double Npyy_tau;
    double nup_tau;

    double gammadot;
    double dxx0;
    double dxy0;
    double dyy0;
    double dpxx;
    double dpxy;
    double dpyy;

    double const c = 1e-2;

    /* solve for g_local */
    calculate_bulk_granular_fluidity(job);

    /* build FEM diffusion array/load vector and solve for g_nonlocal */
    solve_diffusion_part(job);

    /* use g_nonlocal to update stress state (in gf variable) */
    for (i = 0; i < job->num_particles; i++) {
        p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);
        sxx_t0 = job->particles[i].sxx + p_t;
        sxy_t0 = job->particles[i].sxy;
        syy_t0 = job->particles[i].syy + p_t;
        tau_t = sqrt(0.5*(sxx_t0*sxx_t0 + 2*sxy_t0*sxy_t0 + syy_t0*syy_t0));

        /* calculate trial elastic strain */
        dsjxx = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].exx_t) + NUMOD * (job->particles[i].eyy_t));
        dsjxy = job->dt * (EMOD / (2 *(1 + NUMOD))) * (job->particles[i].exy_t);
        dsjyy = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].eyy_t) + NUMOD * (job->particles[i].exx_t));

        /* Jaumann spin terms */
        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

        sxx_tr = job->particles[i].sxx + dsjxx;
        sxy_tr = job->particles[i].sxy + dsjxy;
        syy_tr = job->particles[i].syy + dsjyy;

        p_tr = -0.5 * (sxx_tr + syy_tr);

        if (p_tr < c) {
            job->particles[i].sxx = 0.5 * c / MU_S;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0.5 * c / MU_S;
/*            printf("%s:%s: Particle %zu pressure less than cohesion, results may be inaccurate.\n", __FILE__, __func__, i);*/
            continue;
        }

        sxx_tr0 = sxx_tr + p_tr;
        sxy_tr0 = sxy_tr;
        syy_tr0 = syy_tr + p_tr;

        tau_tr = sqrt(0.5*(sxx_tr0*sxx_tr0 + 2*sxy_tr0*sxy_tr0 + syy_tr0*syy_tr0));

        /* direction of plastic flow */
        if (tau_tr > 0.0) {
            Npxx_tau = sxx_tr0 / tau_tr;
            Npxy_tau = sxy_tr0 / tau_tr;
            Npyy_tau = syy_tr0 / tau_tr;
        } else {
/*            printf("%s:%s: Particle %d has no shear.\n", __FILE__, __func__, i);*/
            Npxx_tau = 0;
            Npxy_tau = 0;
            Npyy_tau = 0;
        }

        tau_tau = p_tr * tau_tr / (p_tr + G * gf * job->dt);

        /* restrict tau */
        if (tau_tau < 0.0) {
            tau_tau = 0.0;
        }
        if (tau_tau > tau_tr) {
            tau_tau = tau_tr;
        }

        if (p_tr > c) {
            nup_tau = gf * tau_tau / p_tr;
        } else {
            nup_tau = gf * tau_tau / c;
        }

        /* restrict dp */
        if (nup_tau < 0) {
            nup_tau = 0;
        }

        /* magnitude of plastic flow */
/*        gammap += nup_tau * job->dt;*/

        /* already calculated above */
/*        tau_tau = tau_tr - G * nup_tau * job->dt;*/

        /* adjust stress */
/*        job->particles[i].sxx = sxx_tr - sqrt(2.0) * (tau_tr - tau_tau) * Npxx_tau;*/
/*        job->particles[i].sxy = sxy_tr - sqrt(2.0) * (tau_tr - tau_tau) * Npxy_tau;*/
/*        job->particles[i].syy = syy_tr - sqrt(2.0) * (tau_tr - tau_tau) * Npyy_tau;*/

        /*
        --
            Explicit constitutive update.
        --
        */
#if 1
        p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);
        sxx_t0 = job->particles[i].sxx + p_t;
        sxy_t0 = job->particles[i].sxy;
        syy_t0 = job->particles[i].syy + p_t;
        tau_t = sqrt(0.5*(sxx_t0*sxx_t0 + 2*sxy_t0*sxy_t0 + syy_t0*syy_t0));
        Npxx_tau = sqrt(0.5) * sxx_t0 / tau_t;
        Npxy_tau = sqrt(0.5) * sxy_t0 / tau_t;
        Npyy_tau = sqrt(0.5) * syy_t0 / tau_t;
        dxx0 = 0.5 * (job->particles[i].exx_t - job->particles[i].eyy_t);
        dxy0 = job->particles[i].exy_t;
        dyy0 = 0.5 * (job->particles[i].eyy_t - job->particles[i].exx_t);
        gammadot = sqrt(2 * (dxx0*dxx0 + 2*dxy0*dxy0 + dyy0*dyy0));
        
        if (p_t > c && gf > 0) {
            nup_tau = gf * tau_t / p_t;

/*            fprintf(job->output.log_fd,*/
/*                "%s:%s: [%zu], g_local = %g g_nonlocal %g diff %e\n",*/
/*                __FILE__, __func__, i, gf_bulk, gf, gf - gf_bulk);*/

            if (nup_tau < 0) {
                nup_tau = 0;
            }
            if (nup_tau > gammadot) {
                nup_tau = gammadot;
            }

            dpxx = sqrt(0.5) * (nup_tau) * Npxx_tau;
            dpxy = sqrt(0.5) * (nup_tau) * Npxy_tau;
            dpyy = sqrt(0.5) * (nup_tau) * Npyy_tau;
/*            fprintf(stderr, "dp != 0\n");*/
/*            dpxx = 0;*/
/*            dpxy = 0;*/
/*            dpyy = 0;*/
        } else {
            nup_tau = 0;
            dpxx = 0;
            dpxy = 0;
            dpyy = 0;
        }

        if (job->t == 0 /* || job->step_number == 0 */) {
            nup_tau = 0;
            dpxx = 0;
            dpxy = 0;
            dpyy = 0;
        }

        dsjxx = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].exx_t - dpxx) + NUMOD * (job->particles[i].eyy_t - dpyy));
        dsjxy = job->dt * (EMOD / (2 *(1 + NUMOD))) * (job->particles[i].exy_t - dpxy);
        dsjyy = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].eyy_t - dpyy) + NUMOD * (job->particles[i].exx_t - dpxx));
        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

        gammap += nup_tau * job->dt;
        job->particles[i].sxx += dsjxx;
        job->particles[i].sxy += dsjxy;
        job->particles[i].syy += dsjyy;
#endif

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

    double const p_cap = 1e-2;

    for (i = 0; i < job->num_particles; i++) {

        /* Calculate pressure at beginning of step. */
        p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);

        sxx_tr0 = job->particles[i].sxx + p_t;
        sxy_tr0 = job->particles[i].sxy;
        syy_tr0 = job->particles[i].syy + p_t;

        /* Calculate equivalent shear at beginning of step. */
        tau_t = sqrt(0.5*(sxx_tr0*sxx_tr0 + 2*sxy_tr0*sxy_tr0 + syy_tr0*syy_tr0));

        if (p_t < p_cap) {
/*            job->particles[i].sxx = -0.5 * p_cap;*/
/*            job->particles[i].sxy = 0;*/
/*            job->particles[i].syy = -0.5 * p_cap;*/
            p_t = p_cap;
/*            printf("%s:%s: Particle %zu pressure less than cohesion, results may be inaccurate.\n", __FILE__, __func__, i);*/
/*            continue;*/
        }

        S0 = p_t * MU_S;
        eta = B * GRAINS_D * sqrt(p_t * GRAINS_RHO);
        xisq_inv = fabs(tau_t / p_t - MU_S)/ (GRAINS_D * GRAINS_D * A * A);

        if (tau_t > S0) { 
            gf_bulk = (p_t / eta) * (1 - S0 / tau_t);
        } else {
            gf_bulk = 0;
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
/*            f[sgi] += job->particles[i].m * gf_bulk * xisq_inv * s[ei];*/
            f[sgi] += (job->particles[i].v) * gf_bulk * xisq_inv * s[ei];
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

            // gf += gf_nodes[gi] * s[ei];
            gf += 0.25 * gf_nodes[gi];
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

