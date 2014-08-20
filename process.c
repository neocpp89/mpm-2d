/**
    \file process.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "interpolate.h"
#include "element.h"
#include "particle.h"
#include "node.h"
#include "process.h"
#include "material.h"
#include "exitcodes.h"
#include "map.h"
#include <suitesparse/cs.h>

#if 0
    #undef assert
    #define assert(x)
#else
    #include <assert.h>
#endif

#define TOL 1e-12

#define signum(x) ((int)((0 < x) - (x < 0)))

#define ijton(i,j,N) ((j)*N + (i))
#define N_TO_P(j,tok,i) SMEAR(j,tok,i,h)
#define DX_N_TO_P(j,tok,c,i) (c) * SMEAR(j,tok,i,b1)
#define DY_N_TO_P(j,tok,c,i) (c) * SMEAR(j,tok,i,b2)

#define IS_VALID_ELEMENT_COORD4(r,c,N) \
    (((r) < (N)) && ((r) >= 0) && ((c) < (N)) && ((c) >= 0))

#define FILL_ELEMENT_NEIGHBOR(en,r,c,N) \
    do { \
        if (IS_VALID_ELEMENT_COORD((r),(c),(N))) { \
            en = (r)*(N) + (c); \
        } else { \
            en = -1; \
        } \
    } while(0)

#define SMEAR SMEAR4
#define ACCUMULATE ACCUMULATE4
#define ACCUMULATE_WITH_MUL ACCUMULATE_WITH_MUL4
#define WHICH_ELEMENT WHICH_ELEMENT4
#define IS_VALID_ELEMENT_COORD IS_VALID_ELEMENT_COORD4

#define __E(j,p) j->in_element[p]
#define __N(j,p,n) j->elements[__E(j,p)].nodes[n]
#define __NE(j,e,n) j->elements[e].nodes[n]

/*#define WHICH_ELEMENT4(xp,yp,N,h) \*/
/*    (floor(xp/h) + floor(yp/h)*(N-1))*/

/* XXX: ugly, fix soon */
#define WHICH_ELEMENT4(xp,yp,N,h) \
    (((xp)<1.0 && (xp)>=0.0 && (yp)<1.0 && (yp)>=0.0)?((floor((xp)/(h)) + floor((yp)/(h))*((N)-1))):(-1))

#define ACCUMULATE4(acc_tok,j,tok,i,n,s) \
    j->nodes[__N(j,i,0)].acc_tok += j->s ## 1[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,1)].acc_tok += j->s ## 2[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,2)].acc_tok += j->s ## 3[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,3)].acc_tok += j->s ## 4[i] * j->particles[i].tok;

#define ACCUMULATE_WITH_MUL4(acc_tok,j,tok,i,n,s,c) \
    j->nodes[__N(j,i,0)].acc_tok += j->s ## 1[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,1)].acc_tok += j->s ## 2[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,2)].acc_tok += j->s ## 3[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,3)].acc_tok += j->s ## 4[i] * j->particles[i].tok * (c);

#define SMEAR4(j,tok,i,s) ( \
    j->s ## 1[i] * j->nodes[__N(j,i,0)].tok + \
    j->s ## 2[i] * j->nodes[__N(j,i,1)].tok + \
    j->s ## 3[i] * j->nodes[__N(j,i,2)].tok + \
    j->s ## 4[i] * j->nodes[__N(j,i,3)].tok \
)

#define CHECK_ACTIVE(j,i) if (j->active[i] == 0) { continue; }

/* Lapack function, double precision direct solver for general matrix */
extern void dgesv_(const int *N, const int *nrhs, double *A, const int *lda, int
    *ipiv, double *b, const int *ldb, int *info);

/* Lapack function, double precision direct solver for symmetric matrix */
extern void dsysv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
    int* ipiv, double* b, int* ldb, double* work, int* lwork, int* info);

/* Lapack function, double precision 2-norm of vector */
extern double dnrm2_(int *n, double *x, int *incx);

/* threaded stress calculation */
void pt_calculate_stress(void *_task)
{
    threadtask_t *task = (threadtask_t *)_task;
    return;
}

/*----------------------------------------------------------------------------*/
job_t *mpm_init(int N, double h, particle_t *particles, int num_particles, double t)
{
    int i;
    int n;

    int r,c;

    job_t *job;

    job = (job_t *)malloc(sizeof(job_t));

    job->t = 0;
    job->t_stop = t;
    job->frame = 0;

    job->N = N;
    job->h = h;
    job->num_nodes = N*N;
    job->num_particles = num_particles;

    job->num_elements = (N - 1) * (N - 1);

    /* Copy particles from given ICs. */
    fprintf(stderr, "Each particle is %zu bytes.\n", sizeof(particle_t));
    job->particles = (particle_t *)malloc(num_particles * sizeof(particle_t));
    fprintf(stderr, "%zu bytes (%.2g MB) allocated for %d particles.\n",
        num_particles * sizeof(particle_t),
        num_particles * sizeof(particle_t) / 1048576.0,
        num_particles);
    memcpy(job->particles, particles, num_particles * sizeof(particle_t));

    /* Set stress, strain to zero. */
    for (i = 0; i < job->num_particles; i++) {
/*        job->particles[i].sxx = 0;*/
/*        job->particles[i].sxy = 0;*/
/*        job->particles[i].syy = 0;*/
        job->particles[i].real_sxx = job->particles[i].sxx;
        job->particles[i].real_sxy = job->particles[i].sxy;
        job->particles[i].real_syy = job->particles[i].syy;

        job->particles[i].exx_t = 0;
        job->particles[i].exy_t = 0;
        job->particles[i].eyy_t = 0;
        job->particles[i].wxy_t = 0;

        job->particles[i].Fxx = 1;
        job->particles[i].Fxy = 0;
        job->particles[i].Fyx = 0;
        job->particles[i].Fyy = 1;

        job->particles[i].ux = 0;
        job->particles[i].uy = 0;

        job->particles[i].material_data = NULL;

        job->particles[i].id = i;
    }
    fprintf(stderr, "Done setting initial particle data.\n");

    /* Get node coordinates. */
    fprintf(stderr, "Each node is %zu bytes.\n", sizeof(node_t));
    job->nodes = (node_t *)malloc(job->num_nodes * sizeof(node_t));
    fprintf(stderr, "%zu bytes (%.2g MB) allocated for %d nodes.\n",
        job->num_nodes * sizeof(node_t),
        job->num_nodes * sizeof(node_t) / 1048576.0,
        job->num_nodes);
    for (i = 0; i < job->num_nodes; i++) {
        node_number_to_coords(&(job->nodes[i].x), &(job->nodes[i].y), i, N, h);
        /* printf("preprocessing: xn=%g yn=%g\n", job->nodes[i].x, job->nodes[i].y); */
    }

    /* Get node numbering for elements. */
    fprintf(stderr, "Each element is %zu bytes.\n", sizeof(element_t));
    job->elements = (element_t *)malloc(job->num_elements * sizeof(element_t));
    fprintf(stderr, "%zu bytes (%.2g MB) allocated for %d elements.\n",
        job->num_elements * sizeof(element_t),
        job->num_elements * sizeof(element_t) / 1048576.0,
        job->num_elements);
    for (i = 0; i < job->num_elements; i++) {
        n = ijton(i % (job->N - 1), i / (job->N - 1), job->N);

        job->elements[i].nodes[0] = n + ijton(0,0,job->N);
        job->elements[i].nodes[1] = n + ijton(1,0,job->N);
        job->elements[i].nodes[2] = n + ijton(1,1,job->N);
        job->elements[i].nodes[3] = n + ijton(0,1,job->N);
    }

    /* Find element neighbors and color element -- cartesian grid only! */
    job->num_colors = 4;
    job->color_indices = (size_t *)malloc(job->num_colors * sizeof(size_t));
/*    job->color_list_lengths = (size_t *)malloc(job->num_colors * sizeof(size_t));*/
/*    for (i = 0; i < job->num_colors; i++) {*/
/*        job->color_list_lengths[i] = 0;*/
/*    }*/
    for (i = 0; i < job->num_elements; i++) {
        r = i / (job->N - 1);
        c = i % (job->N - 1);

        job->elements[i].color = 2 * (r % 2)  + ((c % 2) + 1) - 1;
/*        job->color_list_lengths[job->elements[i].color % job->num_colors]++;*/
/*        fprintf(stderr, "job->elements[%d].color = %d\n", i, job->elements[i].color);*/

        /*
            Fill neighbor array starting with element to the right and
           moving in a postive direction (ccw).
        */
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[0], r, c+1, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[1], r+1, c+1, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[2], r+1, c, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[3], r+1, c-1, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[4], r, c-1, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[5], r-1, c-1, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[6], r-1, c, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[7], r-1, c+1, (job->N - 1));

    }
/*    job->element_id_by_color = (size_t **)malloc(job->num_colors * sizeof(size_t *));*/
/*    for (i = 0; i < job->num_colors; i++) {*/
/*        job->element_id_by_color[i] = (size_t *)malloc(job->color_list_lengths[i] * sizeof(size_t));*/
/*    }*/
/*    for (i = 0; i < job->num_colors; i++) {*/
/*        job->color_indices[i] = 0;*/
/*    }*/
/*    for (i = 0; i < job->num_elements; i++) {*/
/*        c = job->elements[i].color;*/
/*        job->element_id_by_color[job->color_indices[c]] = i;*/
/*        job->color_indices[c]++;*/
/*    }*/
    for (i = 0; i < job->num_colors; i++) {
        job->color_indices[i] = 0;
    }

    /* Allocate space for tracking element->particle map. */
    job->in_element =  (int *)malloc(job->num_particles * sizeof(int));

    /* Allocate space for map of active particles. */
    job->active = (int *)malloc(job->num_particles * sizeof(int));

    /* Allocate space for interpolation functions. */
    job->h1 = (double *)malloc(job->num_particles * sizeof(double));
    job->h2 = (double *)malloc(job->num_particles * sizeof(double));
    job->h3 = (double *)malloc(job->num_particles * sizeof(double));
    job->h4 = (double *)malloc(job->num_particles * sizeof(double));

    job->b11 = (double *)malloc(job->num_particles * sizeof(double));
    job->b12 = (double *)malloc(job->num_particles * sizeof(double));
    job->b13 = (double *)malloc(job->num_particles * sizeof(double));
    job->b14 = (double *)malloc(job->num_particles * sizeof(double));

    job->b21 = (double *)malloc(job->num_particles * sizeof(double));
    job->b22 = (double *)malloc(job->num_particles * sizeof(double));
    job->b23 = (double *)malloc(job->num_particles * sizeof(double));
    job->b24 = (double *)malloc(job->num_particles * sizeof(double));

    /* max size of u_grid is NODAL_DOF * number of nodes. */
    job->vec_len = NODAL_DOF * job->num_nodes;
    job->u_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->du_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->v_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->a_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->m_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->f_ext_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->f_int_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->q_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->node_u_map = (int *)malloc(NODAL_DOF * sizeof(int) * job->num_nodes);
    job->inv_node_u_map = (int *)malloc(NODAL_DOF * sizeof(int) * job->num_nodes);

    /* for dirichlet BCs */
    job->u_dirichlet = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->u_dirichlet_mask = (int *)malloc(NODAL_DOF * sizeof(int) * job->num_nodes);

    /* for periodic BCs */
    job->node_number_override = (int *)malloc(job->vec_len * sizeof(int));

    for (i = 0; i < job->vec_len; i++) {
        job->u_grid[i] = 0;
        job->du_grid[i] = 0;
        job->v_grid[i] = 0;
        job->a_grid[i] = 0;
        job->m_grid[i] = 0;
        job->f_ext_grid[i] = 0;
        job->q_grid[i] = 0;
        job->f_int_grid[i] = 0;
    }

    for (i = 0; i < job->num_particles; i++) {
        job->in_element[i] = WHICH_ELEMENT(
            job->particles[i].x, job->particles[i].y, job->N, job->h);
        if (job->in_element[i] < 0 || job->in_element[i] > job->num_elements) {
            job->active[i] = 0;
        } else {
            job->active[i] = 1;
        }
        if ((job->in_element[i] / (job->N - 1)) % 2 == 0) {
            job->particles[i].color = 0;
        } else {
            job->particles[i].color = 1;
        }
    }

    /* set particle domains */
    for (i = 0; i < job->num_particles; i++) {
        /* seems backwards, but only because loader contains
        current particle volume only. */
        job->particles[i].v0 = job->particles[i].v;
        job->particles[i].r1_initial[0] = 0.5 * sqrt(job->particles[i].v0);
        job->particles[i].r1_initial[1] = 0;
        job->particles[i].r2_initial[0] = 0;
        job->particles[i].r2_initial[1] = 0.5 * sqrt(job->particles[i].v0);
    }

    /* intialize state variables */
/*    material_init(job);*/

    /* Set default timestep. */
    job->dt = 0.4 * job->h * sqrt(job->particles[0].m/(job->particles[0].v * EMOD));

    /* make sure corners know where the are on first step */
    update_particle_vectors(job);
    update_corner_positions(job);

    job->use_cpdi = 0;

    /* initialize mapping matricies */
    job->phi = NULL;
    job->grad_phi = NULL;
    job->phi_transpose = NULL;
    job->grad_phi_transpose = NULL;

    /* nodal quantites */
/*    job->m_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/
/*    job->mx_t_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/
/*    job->my_t_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/
/*    job->mx_tt_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/
/*    job->my_tt_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/
/*    job->x_t_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/
/*    job->y_t_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/
/*    job->x_tt_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/
/*    job->y_tt_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/
/*    job->fx_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/
/*    job->fy_nodes = (double *)malloc(sizeof(double) * job->num_nodes);*/

    /* swap this out later with a pointer from dlopen if needed. */
    job->material.num_fp64_props = 0;
    job->material.num_int_props = 0;
    job->material.fp64_props = NULL;
    job->material.int_props = NULL;
    job->material.material_filename = (char *)malloc(16);
    strcpy(job->material.material_filename, "builtin");
    job->material.material_init = &material_init_linear_elastic;
    job->material.calculate_stress = &calculate_stress_linear_elastic;

    /* used to vary loads/bcs */
    job->step_number = 0;
    job->step_start_time = job->t;

    return job;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void zero_momentum_bc(job_t *job)
{
    int i, j, m, n;

    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            n = NODAL_DOF * i + j;
            m = job->node_number_override[n];
            if (job->u_dirichlet_mask[m] != 0) {
                /* only handle 0 displacement right now. */
                if (job->u_dirichlet[m] == 0) {
                    if (j == XDOF_IDX) {
                        job->nodes[i].mx_t = 0;
                    } else if (j == YDOF_IDX) {
                        job->nodes[i].my_t = 0;
                    }
                }
            }
        }
    }

    /* handle periodicity */
    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            if (job->node_number_override[NODAL_DOF * i + j] != NODAL_DOF * i + j) {
                /* tie nodes together. */
                m = job->node_number_override[NODAL_DOF * i + j];
                n = NODAL_DOF * i + j;

                if (n % NODAL_DOF == XDOF_IDX) {
                    if (m % NODAL_DOF == XDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].mx_t += job->nodes[(n / NODAL_DOF)].mx_t;
                        job->nodes[(n / NODAL_DOF)].mx_t = job->nodes[(m / NODAL_DOF)].mx_t;

                        job->nodes[(m / NODAL_DOF)].mx_tt += job->nodes[(n / NODAL_DOF)].mx_tt;
                        job->nodes[(n / NODAL_DOF)].mx_tt = job->nodes[(m / NODAL_DOF)].mx_tt;
                    } else if (m % NODAL_DOF == YDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].my_t += job->nodes[(n / NODAL_DOF)].mx_t;
                        job->nodes[(n / NODAL_DOF)].mx_t = job->nodes[(m / NODAL_DOF)].my_t;

                        job->nodes[(m / NODAL_DOF)].my_tt += job->nodes[(n / NODAL_DOF)].mx_tt;
                        job->nodes[(n / NODAL_DOF)].mx_tt = job->nodes[(m / NODAL_DOF)].my_tt;
                    }
                } else if (n % NODAL_DOF == YDOF_IDX) {
                    if (m % NODAL_DOF == XDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].mx_t += job->nodes[(n / NODAL_DOF)].my_t;
                        job->nodes[(n / NODAL_DOF)].my_t = job->nodes[(m / NODAL_DOF)].mx_t;

                        job->nodes[(m / NODAL_DOF)].mx_tt += job->nodes[(n / NODAL_DOF)].my_tt;
                        job->nodes[(n / NODAL_DOF)].my_tt = job->nodes[(m / NODAL_DOF)].mx_tt;
                    } else if (m % NODAL_DOF == YDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].my_t += job->nodes[(n / NODAL_DOF)].my_t;
                        job->nodes[(n / NODAL_DOF)].my_t = job->nodes[(m / NODAL_DOF)].my_t;

                        job->nodes[(m / NODAL_DOF)].my_tt += job->nodes[(n / NODAL_DOF)].my_tt;
                        job->nodes[(n / NODAL_DOF)].my_tt = job->nodes[(m / NODAL_DOF)].my_tt;
                    }
                }
            }
        }
    }

    /* fix masses as well */
    for (i = 0; i < job->num_nodes; i++) {
        j = 0;
        if (job->node_number_override[NODAL_DOF * i + j] != NODAL_DOF * i + j) {
            /* tie nodes together. */
            m = job->node_number_override[NODAL_DOF * i + j];
            n = NODAL_DOF * i + j;

            job->nodes[(m / NODAL_DOF)].m += job->nodes[(n / NODAL_DOF)].m;
            job->nodes[(n / NODAL_DOF)].m = job->nodes[(m / NODAL_DOF)].m;
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void zero_force_bc(job_t *job)
{
    int i, j, m, n;

    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            n = NODAL_DOF * i + j;
            m = job->node_number_override[n];
            if (job->u_dirichlet_mask[m] != 0) {
                /* only handle 0 displacement right now. */
                if (job->u_dirichlet[m] == 0) {
                    if (j == XDOF_IDX) {
                        job->nodes[i].fx = 0;
                    } else if (j == YDOF_IDX) {
                        job->nodes[i].fy = 0;
                    }
                }
            }
        }
    }

    /* handle periodicity */
    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            if (job->node_number_override[NODAL_DOF * i + j] != NODAL_DOF * i + j) {
                /* tie nodes together. */
                m = job->node_number_override[NODAL_DOF * i + j];
                n = NODAL_DOF * i + j;

                if (n % NODAL_DOF == XDOF_IDX) {
                    if (m % NODAL_DOF == XDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].fx += job->nodes[(n / NODAL_DOF)].fx;
                        job->nodes[(n / NODAL_DOF)].fx = job->nodes[(m / NODAL_DOF)].fx;
                    } else if (m % NODAL_DOF == YDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].fy += job->nodes[(n / NODAL_DOF)].fx;
                        job->nodes[(n / NODAL_DOF)].fx = job->nodes[(m / NODAL_DOF)].fy;
                    }
                } else if (n % NODAL_DOF == YDOF_IDX) {
                    if (m % NODAL_DOF == XDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].fx += job->nodes[(n / NODAL_DOF)].fy;
                        job->nodes[(n / NODAL_DOF)].fy = job->nodes[(m / NODAL_DOF)].fx;
                    } else if (m % NODAL_DOF == YDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].fy += job->nodes[(n / NODAL_DOF)].fy;
                        job->nodes[(n / NODAL_DOF)].fy = job->nodes[(m / NODAL_DOF)].fy;
                    }
                }
            }
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void explicit_mpm_step_usf(job_t *job)
{
    int i;

    /* Clear grid quantites. */
    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].m = 0;
        job->nodes[i].mx_t = 0;
        job->nodes[i].my_t = 0;
        job->nodes[i].mx_tt = 0;
        job->nodes[i].my_tt = 0;
        job->nodes[i].x_t = 0;
        job->nodes[i].y_t = 0;
        job->nodes[i].x_tt = 0;
        job->nodes[i].y_tt = 0;
        job->nodes[i].fx = 0;
        job->nodes[i].fy = 0;
        job->nodes[i].ux = 0;
        job->nodes[i].uy = 0;
    }

    for (i = 0; i < job->num_elements; i++) {
        job->elements[i].filled = 0;
        job->elements[i].n = 0;
        job->elements[i].m = 0;
    }

    /* Update corner coordinates. NOTE: Not needed unless using cpdi. */
    /* update_corner_domains(job); */

    /* Figure out which element each material point is in. */
    create_particle_to_element_map(job);

    /* Calculate shape and gradient of shape functions. */
    calculate_shapefunctions(job);

    /* Generate phi and grad phi maps. */
/*    generate_mappings(job);*/

    /* Create dirichlet and periodic boundary conditions. */
    generate_dirichlet_bcs(job);
    generate_node_number_override(job);

    /* Map particle state to grid quantites. */
    map_to_grid(job);

    /* Zero perpendicular momentum at edge nodes. */
    zero_momentum_bc(job);

    /* Calculate node velocity. */
    calculate_node_velocity(job);

    /* Calculate strain rate. */
    calculate_strainrate(job);

    /* Calculate stress. */
    (*(job->material.calculate_stress))(job);

    /* Calculate stress. */
/*    task.job = job;*/
/*    for (i = 0; i < PT_NUM_THREADS; i++) {*/
/*        task.offset = i;*/
/*        pthread_create(&(job->threads[i]), NULL, &pt_calculate_stress, &task);*/
/*    }*/
/*    */
/*    for (i = 0; i < PT_NUM_THREADS; i++) {*/
/*        task.offset = i;*/
/*        pthread_join(job->threads[i], NULL);*/
/*    }*/

    update_stress(job);

    /* Zero perpendicular forces at edge nodes. */
    zero_force_bc(job);

    /* Update momentum and velocity at nodes. */
    move_grid(job); 

    /* Update particle position and velocity. */
    move_particles_explicit_usf(job);

    /* Update deformation gradient and volume. */
    /* update_deformation_gradient(job); */

    /* update volume */
    update_particle_densities(job);

    /* Update particle domains. NOTE: Not needed unless using cpdi. */
    /* update_particle_domains(job); */

    /* Increment time. */
    job->t += job->dt;

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void explicit_mpm_step_usf_threaded(job_t *job)
{
    int i;

    /* Clear grid quantites. */
    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].m = 0;
        job->nodes[i].mx_t = 0;
        job->nodes[i].my_t = 0;
        job->nodes[i].mx_tt = 0;
        job->nodes[i].my_tt = 0;
        job->nodes[i].x_t = 0;
        job->nodes[i].y_t = 0;
        job->nodes[i].x_tt = 0;
        job->nodes[i].y_tt = 0;
        job->nodes[i].fx = 0;
        job->nodes[i].fy = 0;
        job->nodes[i].ux = 0;
        job->nodes[i].uy = 0;
    }

    for (i = 0; i < job->num_elements; i++) {
        job->elements[i].filled = 0;
        job->elements[i].n = 0;
        job->elements[i].m = 0;
    }

    /* Update corner coordinates. NOTE: Not needed unless using cpdi. */
    /* update_corner_domains(job); */

    /* Figure out which element each material point is in. */
    create_particle_to_element_map(job);

    /* Calculate shape and gradient of shape functions. */
    calculate_shapefunctions(job);

    /* Generate phi and grad phi maps. */
/*    generate_mappings(job);*/

    /* Create dirichlet and periodic boundary conditions. */
    generate_dirichlet_bcs(job);
    generate_node_number_override(job);

    /* Map particle state to grid quantites. */
    map_to_grid(job);

    /* Zero perpendicular momentum at edge nodes. */
    zero_momentum_bc(job);

    /* Calculate node velocity. */
    calculate_node_velocity(job);

    /* Calculate strain rate. */
    calculate_strainrate(job);

    /* Calculate stress. */
    (*(job->material.calculate_stress))(job);

    /* Calculate stress. */
/*    task.job = job;*/
/*    for (i = 0; i < PT_NUM_THREADS; i++) {*/
/*        task.offset = i;*/
/*        pthread_create(&(job->threads[i]), NULL, &pt_calculate_stress, &task);*/
/*    }*/
/*    */
/*    for (i = 0; i < PT_NUM_THREADS; i++) {*/
/*        task.offset = i;*/
/*        pthread_join(job->threads[i], NULL);*/
/*    }*/

    update_stress(job);

    /* Zero perpendicular forces at edge nodes. */
    zero_force_bc(job);

    /* Update momentum and velocity at nodes. */
    move_grid(job); 

    /* Update particle position and velocity. */
    move_particles_explicit_usf(job);

    /* Update deformation gradient and volume. */
    /* update_deformation_gradient(job); */

    /* update volume */
    update_particle_densities(job);

    /* Update particle domains. NOTE: Not needed unless using cpdi. */
    /* update_particle_domains(job); */

    /* Increment time. */
    job->t += job->dt;

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void explicit_mpm_step_usl(job_t *job)
{
    int i;

    /* Clear grid quantites. */
    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].m = 0;
        job->nodes[i].mx_t = 0;
        job->nodes[i].my_t = 0;
        job->nodes[i].mx_tt = 0;
        job->nodes[i].my_tt = 0;
        job->nodes[i].x_t = 0;
        job->nodes[i].y_t = 0;
        job->nodes[i].x_tt = 0;
        job->nodes[i].y_tt = 0;
        job->nodes[i].fx = 0;
        job->nodes[i].fy = 0;
        job->nodes[i].ux = 0;
        job->nodes[i].uy = 0;
    }

    for (i = 0; i < job->num_elements; i++) {
        job->elements[i].filled = 0;
        job->elements[i].n = 0;
        job->elements[i].m = 0;

        job->elements[i].sxx = 0;
        job->elements[i].sxy = 0;
        job->elements[i].syy = 0;
    }

    /* Update corner coordinates. NOTE: Not needed unless using cpdi. */
    /* update_corner_domains(job); */

    /* Figure out which element each material point is in. */
    create_particle_to_element_map(job);

    /* Calculate shape and gradient of shape functions. */
    calculate_shapefunctions(job);

    /* Create dirichlet and periodic boundary conditions. */
    generate_dirichlet_bcs(job);
    generate_node_number_override(job);

    /* Map particle state to grid quantites. */
    map_to_grid_explicit(job);

    /* Zero perpendicular momentum at edge nodes. */
    zero_momentum_bc(job);

    /*
        This is poorly named -- it actually calculates diverence of stress to
        get nodal force.
    */
    update_stress(job);

    /* Zero perpendicular forces at edge nodes. */
    zero_force_bc(job);

    /* Update momentum and velocity at nodes. */
    move_grid(job); 

    /* Update particle position and velocity. */
    move_particles_explicit_usl(job);

    /* update volume */
    update_particle_densities(job);

    /* Calculate node velocity. */
/*    calculate_node_velocity(job);*/
    /* Unnecessary because move_grid already does this! */

    /* Calculate strain rate. */
    calculate_strainrate(job);

    /* Calculate stress. */
    (*(job->material.calculate_stress))(job);

    /* Increment time. */
    job->t += job->dt;

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void explicit_mpm_step_usl_threaded(void *_task)
{
    threadtask_t *task = (threadtask_t *)_task;
    job_t *job = task->job;
    int rc;
    size_t i;

    size_t p_start = task->offset;
    size_t p_stop = task->offset + task->blocksize;

    size_t n_start = task->n_offset;
    size_t n_stop = task->n_offset + task->n_blocksize;

    size_t e_start = task->e_offset;
    size_t e_stop = task->e_offset + task->e_blocksize;

    pthread_barrier_wait(job->serialize_barrier);

    /* Clear grid quantites. */
    for (i = n_start; i < n_stop; i++) {
        job->nodes[i].m = 0;
        job->nodes[i].mx_t = 0;
        job->nodes[i].my_t = 0;
        job->nodes[i].mx_tt = 0;
        job->nodes[i].my_tt = 0;
        job->nodes[i].x_t = 0;
        job->nodes[i].y_t = 0;
        job->nodes[i].x_tt = 0;
        job->nodes[i].y_tt = 0;
        job->nodes[i].fx = 0;
        job->nodes[i].fy = 0;
        job->nodes[i].ux = 0;
        job->nodes[i].uy = 0;
    }

    for (i = e_start; i < e_stop; i++) {
        job->elements[i].filled = 0;
        job->elements[i].n = 0;
        job->elements[i].m = 0;
    }

    /* Figure out which element each material point is in. */
/*    create_particle_to_element_map_split(job, p_start, p_stop);*/
    create_particle_to_element_map_threaded(task);

    /*
        XXX Normally we find which elements are filled here, but we defer
        because we don't need it until later and want to keep the serial
        sections together after one barrier call.
    */

    /* Calculate shape and gradient of shape functions. */
    calculate_shapefunctions_split(job, p_start, p_stop);

    rc = pthread_barrier_wait(job->serialize_barrier);
    if (rc == PTHREAD_BARRIER_SERIAL_THREAD) {
        /* Increment time. (We do this first to avoid more barrier calls.) */
        job->t += job->dt;

/*    }*/
/*    pthread_barrier_wait(job->serialize_barrier);*/

    /* XXX: most of this is serialized in testing!! */
/*    rc = pthread_barrier_wait(job->serialize_barrier);*/
/*    if (rc == PTHREAD_BARRIER_SERIAL_THREAD) {*/

        /* Update corner coordinates. NOTE: Not needed unless using cpdi. */
        /* update_corner_domains(job); */

        /* find which elements are filled */
        job->update_elementlists_flag = 0;
        for (i = 0; i < job->num_threads; i++) {
            job->update_elementlists_flag += job->update_elementlists[i];
        }

        if (job->update_elementlists_flag != 0) {
            find_filled_elements(job);
        }

        /* Create dirichlet and periodic boundary conditions. */
        generate_dirichlet_bcs(job);
        generate_node_number_override(job);
    }

    pthread_barrier_wait(job->serialize_barrier);
    /* Map particle state to grid quantites. */
    map_to_grid_explicit_split(job, task->id);

    rc = pthread_barrier_wait(job->serialize_barrier);
    if (rc == PTHREAD_BARRIER_SERIAL_THREAD) {
        /* Zero perpendicular momentum at edge nodes. */
        zero_momentum_bc(job);

        /*
            This is poorly named -- it actually calculates diverence of stress to
            get nodal force.
        */
        /* XXX: moved to map_to_grid_explicit_split. */
/*        update_stress(job);*/

        /* Zero perpendicular forces at edge nodes. */
        zero_force_bc(job);

    }

    pthread_barrier_wait(job->serialize_barrier);

#if 0
    if (task->id == 0) {
        double dtcalcmin, mmin = 100000;
        dtcalcmin = job->material.fp64_props[0] * job->h;
        for (i = 0; i < job->num_nodes; i++) {
            if (job->nodes[i].m > 0) {
                if (job->nodes[i].m < mmin) {
                    mmin = job->nodes[i].m;
                }
            }
        }
        dtcalcmin = 2.0 * sqrt(mmin / dtcalcmin);
/*        fprintf(stdout, "\n%g, %g, %g, %g", mmin, job->dt, dtcalcmin, job->dt / dtcalcmin);*/
        if (job->dt > 0.4 * dtcalcmin) {
            fprintf(stdout, " XXX ");
/*            job->dt = 0.4 * dtcalcmin;*/
        }
    }
#endif

    /*
        Update momentum and velocity at nodes.
        Wait for all threads to finish updating nodes before calculating
        strainrates.
    */
    move_grid_split(job, n_start, n_stop); 
    pthread_barrier_wait(job->serialize_barrier);

    /* Update particle position and velocity. */
    move_particles_explicit_usl_split(job, p_start, p_stop);

    /* update volume */
    update_particle_densities_split(job, p_start, p_stop);

    /* Calculate strain rate. */
    calculate_strainrate_split(job, p_start, p_stop);

    /* Calculate stress. */
    (*(job->material.calculate_stress_threaded))(task);

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void implicit_mpm_step(job_t *job)
{
    int i, j;

    /* Clear grid quantites. */
    for (i = 0; i < job->num_nodes; i++) {
/*        job->m_nodes[i] = 0;*/
/*        job->mx_t_nodes[i] = 0;*/
/*        job->my_t_nodes[i] = 0;*/
/*        job->mx_tt_nodes[i] = 0;*/
/*        job->my_tt_nodes[i] = 0;*/
/*        job->fx_nodes[i] = 0;*/
/*        job->fy_nodes[i] = 0;*/

        job->nodes[i].m = 0;
        job->nodes[i].mx_t = 0;
        job->nodes[i].my_t = 0;
        job->nodes[i].mx_tt = 0;
        job->nodes[i].my_tt = 0;
        job->nodes[i].x_t = 0;
        job->nodes[i].y_t = 0;
        job->nodes[i].x_tt = 0;
        job->nodes[i].y_tt = 0;
        job->nodes[i].fx = 0;
        job->nodes[i].fy = 0;

        for (j = 0; j < NODAL_DOF; j++) {
            job->u_grid[NODAL_DOF * i + j] = 0;
            job->du_grid[NODAL_DOF * i + j] = 0;
            job->f_int_grid[NODAL_DOF * i + j] = 0;
            job->f_ext_grid[NODAL_DOF * i + j] = 0;
        }
    }

    for (i = 0; i < job->num_elements; i++) {
        job->elements[i].filled = 0;
        job->elements[i].n = 0;
        job->elements[i].m = 0;
    }

    /* Figure out which element each material point is in. */
    create_particle_to_element_map(job);

    /* Calculate shape and gradient of shape functions. */
    calculate_shapefunctions(job);

    /* Generate phi and grad phi mappings. */
    generate_mappings(job);

    /* Map particle state to grid quantites. */
    map_to_grid(job);

    /* Calculate node velocity. */
    calculate_node_velocity(job);

    /* matrix solve */
    implicit_solve(job);

    /* Update particle position and velocity. */
    move_particles(job);

    /* update deformation gradient, vectors and corners. */
    update_deformation_gradient(job);
    update_particle_vectors(job);
    update_corner_positions(job);

    /* update volume */
    update_particle_densities(job);

    /* Increment time. */
    job->t += job->dt;

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void move_grid(job_t *job)
{
    move_grid_split(job, 0, job->num_nodes);
    return;
}

void move_grid_split(job_t *job, size_t n_start, size_t n_stop)
{
    double m;
    size_t i;

    for (i = n_start; i < n_stop; i++) {
        m = job->nodes[i].m;

        if (m > TOL) {
            job->nodes[i].x_tt = job->nodes[i].fx / m;
            job->nodes[i].y_tt = job->nodes[i].fy / m;

            job->nodes[i].mx_t += job->dt * job->nodes[i].fx;
            job->nodes[i].my_t += job->dt * job->nodes[i].fy;

            job->nodes[i].x_t = job->nodes[i].mx_t / m;
            job->nodes[i].y_t = job->nodes[i].my_t / m;

            job->nodes[i].ux = job->dt * job->nodes[i].x_t;
            job->nodes[i].uy = job->dt * job->nodes[i].y_t;
        } else {
            job->nodes[i].x_tt = 0;
            job->nodes[i].y_tt = 0;
            job->nodes[i].x_t = 0;
            job->nodes[i].y_t = 0;
        }

/*        if (i == 0 && fabs(-0.5 *(job->particles[0].sxx + job->particles[0].syy)) > 1) {*/
/*            printf("\n[%zu:%f %f %f %f %f]", i, job->t, job->particles[0].x, job->particles[0].y, job->particles[0].xl, job->particles[0].yl);*/
/*            printf(" sigma=[%f %f %f], p=%f", job->particles[0].sxx, job->particles[0].sxy, job->particles[0].syy, -0.5 *(job->particles[0].sxx + job->particles[0].syy));*/
/*        }*/

/*        if (i == 0 && job->t < 0.2) {*/
/*            printf("\n[%zu:%f] %f %f %f", job->in_element[0], job->t, m, job->nodes[i].x_t, job->nodes[i].y_t);*/
/*            printf(" sigma=[%f %f %f], p=%f, \tLyy[%g], L=[%g %g %g %g]", job->particles[0].sxx, job->particles[0].sxy, job->particles[0].syy, -0.5 *(job->particles[0].sxx + job->particles[0].syy), job->particles[0].eyy_t, job->particles[0].exx_t, job->particles[0].eyy_t, job->particles[0].exy_t, job->particles[0].wxy_t);*/
/*        }*/

/*        if (job->nodes[i].y_t > 1.0 || job->nodes[i].x_t > 1.0) {*/
/*            printf("\n[%zu:%f] %f %f %f", i, job->t, m, job->nodes[i].x_t, job->nodes[i].y_t);*/
/*            printf(" sigma=[%f %f %f], p=%f, L=[%f %f %f %f]", job->particles[0].sxx, job->particles[0].sxy, job->particles[0].syy, -0.5 *(job->particles[0].sxx + job->particles[0].syy), job->particles[0].exx_t, job->particles[0].eyy_t, job->particles[0].exy_t, job->particles[0].wxy_t);*/
/*        }*/
    }


/*    int i_new;*/
/*    for (i = 0; i < job->num_nodes; i++) {*/
/*        i_new = job->node_number_override[NODAL_DOF * i] / NODAL_DOF;*/

/*        if (i == i_new) { continue; }*/

/*        assert(job->nodes[i].fx == job->nodes[i_new].fx);*/
/*        assert(job->nodes[i].fy == job->nodes[i_new].fy);*/
/*        assert(job->nodes[i].x_tt == job->nodes[i_new].x_tt);*/
/*        assert(job->nodes[i].y_tt == job->nodes[i_new].y_tt);*/
/*        assert(job->nodes[i].x_t == job->nodes[i_new].x_t);*/
/*        assert(job->nodes[i].y_t == job->nodes[i_new].y_t);*/
/*        assert(job->nodes[i].ux == job->nodes[i_new].ux);*/
/*        assert(job->nodes[i].uy == job->nodes[i_new].uy);*/
/*    }*/

    return;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
void create_particle_to_element_map(job_t *job)
{
    create_particle_to_element_map_split(job, 0, job->num_particles);
    find_filled_elements(job);

    return;
}

void create_particle_to_element_map_split(job_t *job, size_t p_start, size_t p_stop)
{
    size_t i;
    int p;

    /* All elements must be cleared before entering this routine! */

    for (i = p_start; i < p_stop; i++) {
        CHECK_ACTIVE(job, i);
        p = WHICH_ELEMENT(
            job->particles[i].x, job->particles[i].y, job->N, job->h);

        if (p != job->in_element[i]) {
            fprintf(job->output.log_fd, 
                "[%g] Particle %d @(%g, %g) left element %d, now in element %d.\n",
                job->t, i, job->particles[i].x, job->particles[i].y,
                job->in_element[i], p);
        }

        /* Update particle element. */
        job->in_element[i] = p;

        if (p == -1) {
            fprintf(job->output.log_fd,
                "[%g] Particle %d outside of grid (%g, %g), marking as inactive.\n",
                job->t, i, job->particles[i].x, job->particles[i].y);
            job->active[i] = 0;
            continue;
        }

    }

    return;
}

void create_particle_to_element_map_threaded(threadtask_t *task)
{
    job_t *job = task->job;
    size_t p_start = task->offset;
    size_t p_stop = task->offset + task->blocksize;

    size_t i;
    int p;
    unsigned int changed = 0;

    /* All elements must be cleared before entering this routine! */


    for (i = p_start; i < p_stop; i++) {
        CHECK_ACTIVE(job, i);
        p = WHICH_ELEMENT(
            job->particles[i].x, job->particles[i].y, job->N, job->h);

        if (p != job->in_element[i]) {
            changed = 1;
            fprintf(job->output.log_fd, 
                "[%g] Particle %d @(%g, %g) left element %d, now in element %d.\n",
                job->t, i, job->particles[i].x, job->particles[i].y,
                job->in_element[i], p);
        }

        /* Update particle element. */
        job->in_element[i] = p;

        if (p == -1) {
            fprintf(job->output.log_fd,
                "[%g] Particle %d outside of grid (%g, %g), marking as inactive.\n",
                job->t, i, job->particles[i].x, job->particles[i].y);
            job->active[i] = 0;
            continue;
        }

    }

    /* set elementlist flag */
    job->update_elementlists[task->id] = changed;

    return;
}


void find_filled_elements(job_t *job)
{
    size_t i, p, c_idx, t_idx, tc_idx, curr_len;

    /* This function should be called ONCE per step (serial function). */

    for (i = 0; i < job->num_colors; i++) {
        job->color_indices[i] = 0;
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        p = job->in_element[i];

        /* Mark element as occupied. */
        job->elements[p].filled = 1;
        job->elements[p].n++;
        job->elements[p].m += job->particles[i].m;
    }

    for (i = 0; i < job->num_elements; i++) {
        if (job->elements[i].filled) {
            /* set color-based element id */
            job->elements[i].color_idx = job->color_indices[job->elements[i].color];
            job->color_indices[job->elements[i].color]++;
        }
    }

    /* Sort the particle IDs for easier element coloring. */
    for (i = 0; i < (job->num_threads * job->num_colors); i++) {
        job->particle_by_element_color_lengths[i] = 0;
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        p = job->in_element[i];
        c_idx = job->elements[p].color;
        t_idx = job->elements[p].color_idx % job->num_threads;
        tc_idx = t_idx * job->num_colors + c_idx;
        curr_len = job->particle_by_element_color_lengths[tc_idx];
        job->particle_by_element_color_lists[tc_idx][curr_len] = i;
        job->particle_by_element_color_lengths[tc_idx]++;
    }

    #if 0
    /*derp*/
    printf("filled elem\n");
    for (t_idx = 0; t_idx < job->num_threads; t_idx++) {
        for (c_idx = 0; c_idx < job->num_colors; c_idx++) {
            printf("%zu %zu: ", t_idx, c_idx);
            for (i = 0; i < job->particle_by_element_color_lengths[t_idx * job->num_colors + c_idx]; i++) {
                printf("%zu ", job->particle_by_element_color_lists[t_idx * job->num_colors + c_idx][i]);
            }
            printf("\n");
        }
    }
    #endif

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_shapefunctions(job_t *job)
{
    calculate_shapefunctions_split(job, 0, job->num_particles);
    return;
}

void calculate_shapefunctions_split(job_t *job, size_t p_start, size_t p_stop)
{
    size_t i, p, n;

    double xn;
    double yn;
    double xl;
    double yl;

    for (i = p_start; i < p_stop; i++) {
        CHECK_ACTIVE(job, i);
        p = job->in_element[i];
        n = job->elements[p].nodes[0];
        xn = job->nodes[n].x;
        yn = job->nodes[n].y;
        global_to_local_coords(&xl, &yl,
            job->particles[i].x, job->particles[i].y, 
            xn, yn, job->h);
        tent(&(job->h1[i]), &(job->h2[i]), &(job->h3[i]), &(job->h4[i]),
            xl, yl);
        grad_tent(
            &(job->b11[i]), &(job->b12[i]), &(job->b13[i]), &(job->b14[i]),
            &(job->b21[i]), &(job->b22[i]), &(job->b23[i]), &(job->b24[i]),
            xl, yl, job->h);
        if (xl < 0.0f || xl > 1.0f || yl < 0.0f || yl > 1.0f) {
            fprintf(stderr, "Particle %d outside of element %d (%g, %g).\n", i,
                p, xl, yl);
        }
        job->particles[i].xl = xl;
        job->particles[i].yl = yl;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_node_velocity(job_t *job)
{
    int i;
    double m;

    for (i = 0; i < job->num_nodes; i++) {
        m = job->nodes[i].m;
        if (m > TOL) {
            job->nodes[i].x_t = job->nodes[i].mx_t / m;
            job->nodes[i].y_t = job->nodes[i].my_t / m;
        } else {
            job->nodes[i].x_t = 0;
            job->nodes[i].y_t = 0;
        }
    }

    for (i = 0; i < job->num_nodes; i++) {
        job->v_grid[NODAL_DOF * i + XDOF_IDX] = job->nodes[i].x_t;
        job->v_grid[NODAL_DOF * i + YDOF_IDX] = job->nodes[i].y_t;
    }

    int i_new;
    for (i = 0; i < job->num_nodes; i++) {
        i_new = job->node_number_override[NODAL_DOF * i] / NODAL_DOF;

        if (i == i_new) { continue; }

        assert(job->nodes[i].ux == job->nodes[i_new].ux);
        assert(job->nodes[i].uy == job->nodes[i_new].uy);
        assert(job->nodes[i].mx_t == job->nodes[i_new].mx_t);
        assert(job->nodes[i].my_t == job->nodes[i_new].my_t);
        assert(job->nodes[i].mx_tt == job->nodes[i_new].mx_tt);
        assert(job->nodes[i].my_tt == job->nodes[i_new].my_tt);
        assert(job->nodes[i].x_t == job->nodes[i_new].x_t);
        assert(job->nodes[i].y_t == job->nodes[i_new].y_t);
        assert(job->nodes[i].x_tt == job->nodes[i_new].x_tt);
        assert(job->nodes[i].y_tt == job->nodes[i_new].y_tt);
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_strainrate(job_t *job)
{
    calculate_strainrate_split(job, 0, job->num_particles);
    return;
}

void calculate_strainrate_split(job_t *job, size_t p_start, size_t p_stop)
{
    int i, j, k;
    int ce, nn[4];
    double dx_tdy;
    double dy_tdx;

    for (i = p_start; i < p_stop; i++) {
        CHECK_ACTIVE(job, i);
        job->particles[i].exx_t = 0;
        job->particles[i].exy_t = 0;
        job->particles[i].eyy_t = 0;
        job->particles[i].wxy_t = 0;

        if (job->use_cpdi) {
            /* loop over corners */
            for (j = 0; j < 4; j++) {
                ce = job->particles[i].corner_elements[j];
                /* corner is outside of particle domain. should probably deal with this better... */
                if (ce == -1) {
                    continue;
                }
                for (k = 0; k < 4; k++) {
                    nn[k] = job->elements[ce].nodes[k];

                    /* actual volume of particle here (not averaging volume) */
                    job->particles[i].exx_t += job->nodes[nn[k]].x_t * job->particles[i].grad_sc[j][k][S_XIDX];
                    job->particles[i].exy_t += 0.5 * (job->nodes[nn[k]].x_t * job->particles[i].grad_sc[j][k][S_YIDX] + job->nodes[nn[k]].y_t * job->particles[i].grad_sc[j][k][S_XIDX]);
                    job->particles[i].eyy_t += job->nodes[nn[k]].y_t * job->particles[i].grad_sc[j][k][S_YIDX];
                    job->particles[i].wxy_t += 0.5 * (job->nodes[nn[k]].x_t * job->particles[i].grad_sc[j][k][S_YIDX] - job->nodes[nn[k]].y_t * job->particles[i].grad_sc[j][k][S_XIDX]);
                }
            }
        } else {
            job->particles[i].exx_t = DX_N_TO_P(job, x_t, 1.0, i);
            job->particles[i].eyy_t = DY_N_TO_P(job, y_t, 1.0, i);

            dx_tdy = DY_N_TO_P(job, x_t, 1.0, i);
            dy_tdx = DX_N_TO_P(job, y_t, 1.0, i);

            job->particles[i].exy_t = 0.5 * (dx_tdy + dy_tdx);
            job->particles[i].wxy_t = 0.5 * (dx_tdy - dy_tdx);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_stress(job_t *job)
{
    int i, j, k, method;
    int ce, nn[4];

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        if (job->use_cpdi) {
            /* loop over corners */
            for (j = 0; j < 4; j++) {
                ce = job->particles[i].corner_elements[j];
                for (k = 0; k < 4; k++) {
                    nn[k] = job->elements[ce].nodes[k];
                    
                    /* actual volume of particle here (not averaging volume) */
                    job->nodes[nn[k]].fx += -job->particles[i].v * (
                    job->particles[i].grad_sc[j][k][S_XIDX] * job->particles[i].sxx + 
                    job->particles[i].grad_sc[j][k][S_YIDX] * job->particles[i].sxy);

                    job->nodes[nn[k]].fy += -job->particles[i].v * (
                    job->particles[i].grad_sc[j][k][S_XIDX] * job->particles[i].sxy + 
                    job->particles[i].grad_sc[j][k][S_YIDX] * job->particles[i].syy);
                }
            }
        } else {
            ACCUMULATE_WITH_MUL(fx, job, sxx, i, n, b1, -job->particles[i].v);
            ACCUMULATE_WITH_MUL(fx, job, sxy, i, n, b2, -job->particles[i].v);
            ACCUMULATE_WITH_MUL(fy, job, sxy, i, n, b1, -job->particles[i].v);
            ACCUMULATE_WITH_MUL(fy, job, syy, i, n, b2, -job->particles[i].v);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_stress_threaded(job_t *job, size_t thread_id)
{
    int i, j, k, method;
    int ce, nn[4];

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        if (job->use_cpdi) {
            /* loop over corners */
            for (j = 0; j < 4; j++) {
                ce = job->particles[i].corner_elements[j];
                for (k = 0; k < 4; k++) {
                    nn[k] = job->elements[ce].nodes[k];
                    
                    /* actual volume of particle here (not averaging volume) */
                    job->nodes[nn[k]].fx += -job->particles[i].v * (
                    job->particles[i].grad_sc[j][k][S_XIDX] * job->particles[i].sxx + 
                    job->particles[i].grad_sc[j][k][S_YIDX] * job->particles[i].sxy);

                    job->nodes[nn[k]].fy += -job->particles[i].v * (
                    job->particles[i].grad_sc[j][k][S_XIDX] * job->particles[i].sxy + 
                    job->particles[i].grad_sc[j][k][S_YIDX] * job->particles[i].syy);
                }
            }
        } else {
            ACCUMULATE_WITH_MUL(fx, job, sxx, i, n, b1, -job->particles[i].v);
            ACCUMULATE_WITH_MUL(fx, job, sxy, i, n, b2, -job->particles[i].v);
            ACCUMULATE_WITH_MUL(fy, job, sxy, i, n, b1, -job->particles[i].v);
            ACCUMULATE_WITH_MUL(fy, job, syy, i, n, b2, -job->particles[i].v);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_internal_stress(job_t *job)
{
    int i, j, k, s;
    int p;
    int ce, nn[4];

    double sxx, sxy, syy;
    double pxx, pxy, pyx, pyy;
    double dudx, dudy, dvdx, dvdy, jdet;

    for (i = 0; i < job->vec_len; i++) {
        job->f_int_grid[i] = 0;
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        p = job->in_element[i];
        for (j = 0; j < 4; j++) {
            nn[j]  = job->elements[p].nodes[j];
        }

        if (job->use_cpdi) {
            /* loop over corners */
            for (j = 0; j < 4; j++) {
                ce = job->particles[i].corner_elements[j];
                /* corner is outside of particle domain. should probably deal with this better... */
                if (ce == -1) {
                    continue;
                }
                for (k = 0; k < 4; k++) {
                    nn[k] = job->elements[ce].nodes[k];

                    /* actual volume of particle here (not averaging volume) */
                    job->f_int_grid[NODAL_DOF * nn[k] + XDOF_IDX] += job->particles[i].v * (
                    job->particles[i].grad_sc[j][k][S_XIDX] * job->particles[i].sxx + 
                    job->particles[i].grad_sc[j][k][S_YIDX] * job->particles[i].sxy);

                    job->f_int_grid[NODAL_DOF * nn[k] + YDOF_IDX] += job->particles[i].v * (
                    job->particles[i].grad_sc[j][k][S_XIDX] * job->particles[i].sxy + 
                    job->particles[i].grad_sc[j][k][S_YIDX] * job->particles[i].syy);


/*                    printf("update_stress: %d\n", nn[k]);*/
/*                    printf("vol %g\ngrad_sc [%g %g]\ns [%g %g %g]\n", job->particles[i].v,*/
/*                     job->particles[i].grad_sc[j][k][S_XIDX], job->particles[i].grad_sc[j][k][S_YIDX], job->particles[i].sxx, job->particles[i].sxy, job->particles[i].syy);*/

/*                    for (s = 0; s < job->num_nodes; s++) {*/
/*                        printf("f[%d] = [%g %g]^T\n", s,*/
/*                            job->f_int_grid[NODAL_DOF * s + XDOF_IDX],*/
/*                            job->f_int_grid[NODAL_DOF * s + YDOF_IDX]);*/
/*                    }*/

                }
/*                exit(254);*/
            }
        } else {
#if 0
            dudx = job->u_grid[NODAL_DOF * nn[0] + XDOF_IDX]*job->b11[i];
            dudx += job->u_grid[NODAL_DOF * nn[1] + XDOF_IDX]*job->b12[i];
            dudx += job->u_grid[NODAL_DOF * nn[2] + XDOF_IDX]*job->b13[i];
            dudx += job->u_grid[NODAL_DOF * nn[3] + XDOF_IDX]*job->b14[i];

            dudy = job->u_grid[NODAL_DOF * nn[0] + XDOF_IDX]*job->b21[i];
            dudy += job->u_grid[NODAL_DOF * nn[1] + XDOF_IDX]*job->b22[i];
            dudy += job->u_grid[NODAL_DOF * nn[2] + XDOF_IDX]*job->b23[i];
            dudy += job->u_grid[NODAL_DOF * nn[3] + XDOF_IDX]*job->b24[i];

            dvdx = job->u_grid[NODAL_DOF * nn[0] + YDOF_IDX]*job->b11[i];
            dvdx += job->u_grid[NODAL_DOF * nn[1] + YDOF_IDX]*job->b12[i];
            dvdx += job->u_grid[NODAL_DOF * nn[2] + YDOF_IDX]*job->b13[i];
            dvdx += job->u_grid[NODAL_DOF * nn[3] + YDOF_IDX]*job->b14[i];

            dvdy = job->u_grid[NODAL_DOF * nn[0] + YDOF_IDX]*job->b21[i];
            dvdy += job->u_grid[NODAL_DOF * nn[1] + YDOF_IDX]*job->b22[i];
            dvdy += job->u_grid[NODAL_DOF * nn[2] + YDOF_IDX]*job->b23[i];
            dvdy += job->u_grid[NODAL_DOF * nn[3] + YDOF_IDX]*job->b24[i];

            jdet = 1 + dudx + dvdy;

            sxx = job->particles[i].sxx;
            sxy = job->particles[i].sxy;
            syy = job->particles[i].syy;

            pxx = jdet * ((1 - dudx) * sxx - dudy * sxy);
            pxy = jdet * ((-dvdx) * sxx + (1 - dvdy) * sxy);
            pyx = jdet * ((1 - dudx) * sxy - dudy * syy);
            pyy = jdet * ((-dvdx) * sxy + (1 - dvdy) * syy);
#endif

            sxx = job->particles[i].sxx;
            sxy = job->particles[i].sxy;
            syy = job->particles[i].syy;
            pxx = sxx;
            pxy = sxy;
            pyx = sxy;
            pyy = syy;

            job->f_int_grid[NODAL_DOF * nn[0] + XDOF_IDX] += job->particles[i].v*(job->b11[i]*pxx + job->b21[i]*pxy);
            job->f_int_grid[NODAL_DOF * nn[0] + YDOF_IDX] += job->particles[i].v*(job->b11[i]*pyx + job->b21[i]*pyy);
            job->f_int_grid[NODAL_DOF * nn[1] + XDOF_IDX] += job->particles[i].v*(job->b12[i]*pxx + job->b22[i]*pxy);
            job->f_int_grid[NODAL_DOF * nn[1] + YDOF_IDX] += job->particles[i].v*(job->b12[i]*pyx + job->b22[i]*pyy);
            job->f_int_grid[NODAL_DOF * nn[2] + XDOF_IDX] += job->particles[i].v*(job->b13[i]*pxx + job->b23[i]*pxy);
            job->f_int_grid[NODAL_DOF * nn[2] + YDOF_IDX] += job->particles[i].v*(job->b13[i]*pyx + job->b23[i]*pyy);
            job->f_int_grid[NODAL_DOF * nn[3] + XDOF_IDX] += job->particles[i].v*(job->b14[i]*pxx + job->b24[i]*pxy);
            job->f_int_grid[NODAL_DOF * nn[3] + YDOF_IDX] += job->particles[i].v*(job->b14[i]*pyx + job->b24[i]*pyy);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void implicit_solve(job_t *job)
{
    int i, j;
    int i_new, j_new;
    int r, s;
    int m, n;

    /* global and elemental indicies */
    int gi, gj;
    int ei, ej;
    double ke;

    /* lapack related*/
    int lda, ldb;

    /* dnrm2 */
    int incx = 1;
    double du_norm = 0;
    double du0_norm = 0;
    double q_norm = 0;
    double q0_norm = 0;

    /* For assembly of elements into global matrix. */
    int nn[4];

    /* iteration count */
    int k = 0;
    static int stable_timestep = 0;

    /* old timestep grid velocity and acceleration */
    double *v_grid_t;
    double *a_grid_t;

    /* mass weighted old timestep grid velocity and acceleration */
    double *mv_grid_t;
    double *ma_grid_t;
    double *m_grid_t;

    /* csparse */
    int nnz;
    int slda, sldb;
    cs *triplets;
    cs *smat;
    double *sb;
    int res = 0;

    /* for small timestep */
    double inv_dt = 1.0 / job->dt;
    double inv_dt_sq = inv_dt / job->dt;

    /* Timing code */
    struct timespec requestStart, requestEnd;
    long ns;

    lda = job->vec_len;
    ldb = lda;

    v_grid_t = (double *)malloc(lda * sizeof(double));
    a_grid_t = (double *)malloc(lda * sizeof(double));

    mv_grid_t = (double *)malloc(lda * sizeof(double));
    ma_grid_t = (double *)malloc(lda * sizeof(double));
    m_grid_t = (double *)malloc(lda * sizeof(double));

    /* keep old timestep velocity and acceleration */
    memcpy(v_grid_t, job->v_grid, lda * sizeof(double));
    memcpy(a_grid_t, job->a_grid, lda * sizeof(double));
/*    for (i = 0; i < ldb; i++) {*/
/*        mv_grid_t[i] = v_grid_t[i] * job->m_grid[i];*/
/*        ma_grid_t[i] = a_grid_t[i] * job->m_grid[i];*/
/*    }*/

    /* keep old stresses */
    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        memcpy(job->particles[i].real_state, job->particles[i].state, DEPVAR * sizeof(double));
        job->particles[i].real_sxx = job->particles[i].sxx;
        job->particles[i].real_sxy = job->particles[i].sxy;
        job->particles[i].real_syy = job->particles[i].syy;
    }

start_implicit:
    k = 0;
    /* no initial grid deformation */
    for (i = 0; i < ldb; i++) {
        job->u_grid[i] = 0;
    }

    /* Generate and apply dirichlet boundary conditions given in BC file. */
    generate_dirichlet_bcs(job);

    /* Generate the new nodal numbers according to BCs. */
    generate_node_number_override(job);

    /*
        We have to change the velocities and accelerations if nodes are tied
        together. Weigh each component by nodal mass.
    */
    for (i = 0; i < lda; i++) {
        mv_grid_t[i] = 0;
        ma_grid_t[i] = 0;
        m_grid_t[i] = 0;
    }
    for (i = 0; i < lda; i++) {
        i_new = job->node_number_override[i];
        mv_grid_t[i_new] += v_grid_t[i] * job->m_grid[i];
        ma_grid_t[i_new] += a_grid_t[i] * job->m_grid[i];
        m_grid_t[i_new] += job->m_grid[i];
    }
    for (i = 0; i < lda; i++) {
        i_new = job->node_number_override[i];
        if (fabs(m_grid_t[i_new]) < TOL) {
            v_grid_t[i] = 0;
            a_grid_t[i] = 0;
            continue;
        }
        v_grid_t[i] = mv_grid_t[i_new] / m_grid_t[i_new];
        a_grid_t[i] = ma_grid_t[i_new] / m_grid_t[i_new];
    }

/*    for (i = 0; i < job->num_nodes; i++) {*/
/*        printf("at[%d] = [%g %g]\n",*/
/*            a_grid_t[NODAL_DOF * i + XDOF_IDX], i, a_grid_t[NODAL_DOF * i + YDOF_IDX]);*/
/*        printf("vt[%d] = [%g %g]\n",*/
/*            v_grid_t[NODAL_DOF * i + XDOF_IDX], i, v_grid_t[NODAL_DOF * i + YDOF_IDX]);*/
/*    }*/
/*    getchar();*/

    /* Begin newton iterations */
    do {
        update_internal_stress(job);
        build_elemental_stiffness(job);

        /* starting timer */
        clock_gettime(CLOCK_REALTIME, &requestStart);

        /* build node and inverse maps, taking BCs into account. */
        for (i = 0; i < lda; i++) {
            job->node_u_map[i] = -1;
            job->inv_node_u_map[i] = -1;
        }
        j = 0;
        for (i = 0; i < lda; i++) {
            i_new = job->node_number_override[i];

            /* if the new node dof has already been assigned a place in the
                matrix, don't set it again */
            if (job->node_u_map[i_new] != -1) {
                continue;
            }

            /* otherwise, check if this node has any particles and the value
                has not been set by dirichlet BCs. */
            if (job->m_grid[i] > TOL && (job->u_dirichlet_mask[i] == 0)) {
                job->node_u_map[i_new] = j;
                job->inv_node_u_map[j] = i_new;
                j++;
            }
        }
        slda = j;   /* The number of free DOFs gives the sparse matrix size. */

        /*  calculate number of nonzero elements (before summing duplicates). */
        nnz = 0;
        for (i = 0; i < job->num_elements; i++) {
            if (job->elements[i].filled != 0) {
                nnz += (NODAL_DOF * NODES_PER_ELEMENT) *
                        (NODAL_DOF * NODES_PER_ELEMENT);
                nnz += NODAL_DOF; /* for diagonal mass matrix entries. */
            }
        }

        /* allocate triplet array. */
        fprintf(job->output.log_fd, "%d nonzeros (including duplicates) in K, which has size [%d x %d].\n", nnz, slda, slda);
        triplets = cs_spalloc(slda, slda, nnz, 1, 1);

        /* Calculate right hand side (load):
                Q = f_ext - f_int - M_g * (4 * u / dt^2 - 4 * v / dt - a) */
        for (i = 0; i < lda; i++) {
            job->q_grid[i] = 0;
        }
        for (i = 0; i < lda; i++) {
            i_new = job->node_number_override[i];
            job->q_grid[i_new] = (-4 * m_grid_t[i_new] * job->u_grid[i_new] + 4 * mv_grid_t[i_new] * job->dt + ma_grid_t[i_new] * job->dt * job->dt);
/*            job->q_grid[i_new] = (1 * mv_grid_t[i_new] * job->dt);*/
        }
        for (i = 0; i < lda; i++) {
            i_new = job->node_number_override[i];
            job->q_grid[i_new] += (job->f_ext_grid[i] - job->f_int_grid[i]) * job->dt * job->dt;
        }

        /* assemble elements into global stiffness matrix */
        for (i = 0; i < job->num_elements; i++) {

            /* Ignore empty elements. */
            if (job->elements[i].filled == 0) {
                continue;
            }

            /* scale stiffness matrix */
            for (r = 0; r < (NODAL_DOF * NODES_PER_ELEMENT); r++) {
                for (s = 0; s < (NODAL_DOF * NODES_PER_ELEMENT); s++) {
                    job->elements[i].kku_element[r][s] *= (job->dt * job->dt);
                }
            }
            nn[0] = job->elements[i].nodes[0];
            nn[1] = job->elements[i].nodes[1];
            nn[2] = job->elements[i].nodes[2];
            nn[3] = job->elements[i].nodes[3];

            /*
                This presumes that the element level stiffness matrices are
                written as (x1, y1, c1, .... , x4, y4, c4) where the
                numbers indicate local node and the variable names are
                nodal degrees of freedom.
            */
            for (r = 0; r < NODES_PER_ELEMENT; r++) {
                for (m = 0; m < NODAL_DOF; m++) {
                    for (s = 0; s < NODES_PER_ELEMENT; s++) {
                        for (n = 0; n < NODAL_DOF; n++) {
                            gi = NODAL_DOF * nn[r] + m;
                            gj = NODAL_DOF * nn[s] + n;

                            /* renumber nodes as appropriate */
                            gi = job->node_number_override[gi];
                            gj = job->node_number_override[gj];

                            ei = NODAL_DOF * r + m;
                            ej = NODAL_DOF * s + n;
                            ke = job->elements[i].kku_element[ei][ej];

                            /*
                                If the global index is masked, don't add it
                                to the sparse matrix. We will correct the
                                solution later.
                            */
                            if (job->u_dirichlet_mask[gi] != 0) {
                                continue;
                            }

                            /* Adjust load for dirichlet bcs. */
                            if (job->u_dirichlet_mask[gj] != 0) {
                                job->q_grid[gi] += (-job->u_dirichlet[gj] * ke);
                                continue;
                            }

                            /* Don't add DOFs which are outside the problem. */
                            if (job->node_u_map[gi] == -1 ||
                                job->node_u_map[gj] == -1) {
                                continue;
                            }

                            res = cs_entry(triplets,
                                    job->node_u_map[gi],
                                    job->node_u_map[gj],
                                    ke);

                            if (res == 0) {
                                fprintf(stderr, "error adding stiffness entry\n");
/*                                fprintf(stderr, "gi=%d gj=%d nmap[gi]=%d nmap[gj]=%d\n",*/
/*                                    gi, gj, job->node_u_map[gi], job->node_u_map[gj]);*/
                            }
                        }
                    }
                }
            }
        }

        for (i = 0; i < lda; i++) {
            i_new = job->node_number_override[i];
            if (job->node_u_map[i_new] != -1) {
                /* add diagonal mass matrix entries to
                    displacement dofs. */

                res = cs_entry(triplets,
                        job->node_u_map[i_new],
                        job->node_u_map[i_new],
                        4 * job->m_grid[i]);

/*                res = cs_entry(triplets,*/
/*                        job->node_u_map[i_new],*/
/*                        job->node_u_map[i_new],*/
/*                        1 * job->m_grid[i]);*/

                if (res == 0) {
                    fprintf(stderr, "error adding mass entry\n");
                }
            }
        }

        /* stopping timer */
        clock_gettime(CLOCK_REALTIME, &requestEnd);

        /* Calculate time it took */
    #define NS_PER_S 1E9
        ns = NS_PER_S * (requestEnd.tv_sec - requestStart.tv_sec )
          + ( requestEnd.tv_nsec - requestStart.tv_nsec );
        fprintf(job->output.log_fd, "Global Assembly[%d]: %lf s\n", k, (double)ns/NS_PER_S );
    #undef NS_PER_S

        /* reset stress state at beginning of newton iteration */
        for (i = 0; i < job->num_particles; i++) {
            CHECK_ACTIVE(job, i);
            
            memcpy(job->particles[i].state, job->particles[i].real_state, DEPVAR * sizeof(double));
            job->particles[i].sxx = job->particles[i].real_sxx;
            job->particles[i].sxy = job->particles[i].real_sxy;
            job->particles[i].syy = job->particles[i].real_syy;
        }

        sb = (double *)malloc(slda * sizeof(double));
        for (i = 0; i < slda; i++) {
            sb[i] = job->q_grid[job->inv_node_u_map[i]];
        }

        /* change to compressed format to sum up entries. */
        smat = cs_compress(triplets);

        /* sum entries with same indicies. */
        res = cs_dupl(smat);

        if (res == 0) {
            fprintf(stderr, "error summing entries\n");
            exit(EXIT_ERROR_CS_DUP);
        }

/*        cs_print(smat, 0);*/

        /* starting timer */
        clock_gettime(CLOCK_REALTIME, &requestStart);

        sldb = slda;
        if (k == 0) {
            //q0_norm = dnrm2_(&sldb, sb, &incx);
        }
        //q_norm = dnrm2_(&sldb, sb, &incx);

        if (!cs_lusol(1, smat, sb, 1e-12)) {
            fprintf(stderr, "lusol error!\n");
            if (cs_qrsol(1, smat, sb)) {
                fprintf(stderr, "qrsol error!\n");
                exit(EXIT_ERROR_CS_SOL);
            }
        }

        if (k == 0) {
            //du0_norm = dnrm2_(&sldb, sb, &incx);
        }
        //du_norm = dnrm2_(&sldb, sb, &incx);

        /* stopping timer */
        clock_gettime(CLOCK_REALTIME, &requestEnd);

        /* Calculate time it took */
    #define NS_PER_S 1E9
        ns = NS_PER_S * (requestEnd.tv_sec - requestStart.tv_sec )
          + ( requestEnd.tv_nsec - requestStart.tv_nsec );
        fprintf(job->output.log_fd, "Iteration[%d]: Norm du: %f Norm q: %f Implicit Solve: %lf s\n", k, du_norm, q_norm, (double)ns/NS_PER_S );
    #undef NS_PER_S

        /* map solution back to entire grid */
        for (i = 0; i < ldb; i++) {
            i_new = job->node_number_override[i];
            if ((j = job->node_u_map[i_new]) != -1) {
                job->du_grid[i] = sb[j];
            } else {
                job->du_grid[i] = 0;
            }
        }

        /* apply dirichlet BCs to solution */
        for (i = 0; i < ldb; i++) {
            i_new = job->node_number_override[i];
            if (job->u_dirichlet_mask[i_new] != 0) {
                job->du_grid[i] = job->u_dirichlet[i_new];
            }
        }

        for (i = 0; i < ldb; i++) {
            /* update grid displacement */
            job->u_grid[i] += job->du_grid[i];

            /* update grid velocity */
            job->v_grid[i] = 2 * (job->u_grid[i] / job->dt) - v_grid_t[i];
/*            job->v_grid[i] = (job->u_grid[i] / job->dt);*/
        }

        for (i = 0; i < job->num_nodes; i++) {
            job->nodes[i].x_t = job->v_grid[NODAL_DOF * i + XDOF_IDX];
            job->nodes[i].y_t = job->v_grid[NODAL_DOF * i + YDOF_IDX];
        }

        calculate_strainrate(job);
        (*(job->material.calculate_stress))(job);

        cs_spfree(triplets);
        cs_spfree(smat);
        free(sb);
        k++;

        if (du_norm < job->implicit.du_norm_converged) {
            fprintf(job->output.log_fd, "norm of du  = %e: Accepting as converged.\n", du_norm);
            break;
        }

        if (q_norm*du_norm > 10*(q0_norm*du0_norm)
                || du_norm > 10*du0_norm
                || k > job->implicit.unstable_iteration_count) {
            job->dt = job->dt * 0.8;
            stable_timestep = 0;
            fprintf(job->output.log_fd, "Trouble converging, modifying Timestep to %f.\n", job->dt);
            goto start_implicit;
        }

        if (job->dt < job->timestep.dt_min) {
            fprintf(stderr, "Timestep %g is too small, aborting.\n", job->dt);
            exit(EXIT_ERROR_DT_TOO_SMALL);
        }

    } while ((du_norm / du0_norm) > job->implicit.du_norm_ratio ||
        ((du_norm*q_norm) / (du0_norm*q0_norm)) > job->implicit.q_norm_ratio);

    stable_timestep++;
    if (job->timestep.allow_dt_increase != 0
        && stable_timestep > job->timestep.stable_dt_threshold 
        && job->dt < job->timestep.dt_max) {
        job->dt = job->dt * 1.25;
        if (job->dt > job->timestep.dt_max) {
            job->dt = job->timestep.dt_max;
        }
        fprintf(job->output.log_fd, "Increasing timestep to %f.\n", job->dt);
        stable_timestep = 0;
    }

    /* update grid acceleration */
    for (i = 0; i < ldb; i++) {
        job->a_grid[i] = (4 * job->u_grid[i] * inv_dt_sq - 4 * v_grid_t[i] * inv_dt - a_grid_t[i]);
/*        job->a_grid[i] = (job->u_grid[i] * inv_dt_sq - v_grid_t[i] * inv_dt);*/
    }

    /* repackage for use with macros */
    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].x_tt = job->a_grid[NODAL_DOF * i + XDOF_IDX];
        job->nodes[i].y_tt = job->a_grid[NODAL_DOF * i + YDOF_IDX];

        job->nodes[i].x_t = job->v_grid[NODAL_DOF * i + XDOF_IDX];
        job->nodes[i].y_t = job->v_grid[NODAL_DOF * i + YDOF_IDX];

        job->nodes[i].ux = job->u_grid[NODAL_DOF * i + XDOF_IDX];
        job->nodes[i].uy = job->u_grid[NODAL_DOF * i + YDOF_IDX];
    }

    for (i = 0; i < job->num_nodes; i++) {
        i_new = job->node_number_override[NODAL_DOF * i] / NODAL_DOF;

        if (i == i_new) { continue; }

/*        if (fabs(job->nodes[i].x_t) > TOL) {*/
/*            printf("u[%d] = [%g %g]\n", i, job->nodes[i].ux, job->nodes[i].uy);*/
/*            printf("v[%d] = [%g %g]\n", i, job->nodes[i].x_t, job->nodes[i].y_t);*/
/*            printf("a[%d] = [%g %g]\n", i, job->nodes[i].x_tt, job->nodes[i].y_tt);*/
/*            printf("f_int[%d] = [%g %g]\n", i, job->f_int_grid[NODAL_DOF * i +XDOF_IDX], job->f_int_grid[NODAL_DOF * i + YDOF_IDX]);*/
/*            printf("f_ext[%d] = [%g %g]\n", i, job->f_ext_grid[NODAL_DOF * i +XDOF_IDX], job->f_ext_grid[NODAL_DOF * i + YDOF_IDX]);*/
/*        }*/

        assert(job->nodes[i].ux == job->nodes[i_new].ux);
        assert(job->nodes[i].uy == job->nodes[i_new].uy);
        assert(job->nodes[i].x_t == job->nodes[i_new].x_t);
        assert(job->nodes[i].y_t == job->nodes[i_new].y_t);
        assert(job->nodes[i].x_tt == job->nodes[i_new].x_tt);
        assert(job->nodes[i].y_tt == job->nodes[i_new].y_tt);
    }

    free(v_grid_t);
    free(a_grid_t);
    free(mv_grid_t);
    free(ma_grid_t);
    free(m_grid_t);

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void map_to_grid_explicit(job_t *job)
{
    size_t i;
    double s[NODES_PER_ELEMENT];

    /* kinda a hack */
    const int pdata_len = 7;
    size_t node_field_offsets[pdata_len];
    double pdata[pdata_len];

    int p;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        s[0] = job->h1[i];
        s[1] = job->h2[i];
        s[2] = job->h3[i];
        s[3] = job->h4[i];

        p = job->in_element[i];

        /* Mass. */
        pdata[0] = job->particles[i].m;
        node_field_offsets[0] = offsetof(node_t, m);

        /* Momentum. */
        pdata[1] = job->particles[i].x_t * job->particles[i].m;
        pdata[2] = job->particles[i].y_t * job->particles[i].m;
        node_field_offsets[1] = offsetof(node_t, mx_t);
        node_field_offsets[2] = offsetof(node_t, my_t);

        /* Pseudoforce. */
        pdata[3] = job->particles[i].x_tt * job->particles[i].m;
        pdata[4] = job->particles[i].y_tt * job->particles[i].m;
        node_field_offsets[3] = offsetof(node_t, mx_tt);
        node_field_offsets[4] = offsetof(node_t, my_tt);

        /* Body forces. */
        pdata[5] = job->particles[i].bx * job->particles[i].m;
        pdata[6] = job->particles[i].by * job->particles[i].m;
        node_field_offsets[5] = offsetof(node_t, fx);
        node_field_offsets[6] = offsetof(node_t, fy);

        accumulate_p_to_n_ds_list(job->nodes,
            node_field_offsets, job->elements[p].nodes, s, 
            NODES_PER_ELEMENT, pdata, pdata_len);

    }

    return;
}

void map_to_grid_explicit_split(job_t *job, size_t thread_id)
{
    size_t c, i, p_idx, tc_idx;
    double s[NODES_PER_ELEMENT];
    double ds[NODES_PER_ELEMENT];

    /* Kinda a hack, but used this because C doesn't have reflection. */
    const int pdata_len = 7;
    size_t node_field_offsets[pdata_len];
    double pdata[pdata_len];

    /* also accumulate stress using gradients of shapefunctions. */
    const int stress_len = 3;
    size_t node_d_offsets[stress_len];
    double stressdata[stress_len];

    int p;

    /* Be sure to replicate this order in the particle data array! */
    node_field_offsets[0] = offsetof(node_t, m);
    node_field_offsets[1] = offsetof(node_t, mx_t);
    node_field_offsets[2] = offsetof(node_t, my_t);
    node_field_offsets[3] = offsetof(node_t, mx_tt);
    node_field_offsets[4] = offsetof(node_t, my_tt);
    node_field_offsets[5] = offsetof(node_t, fx);
    node_field_offsets[6] = offsetof(node_t, fy);

    node_d_offsets[0] = offsetof(node_t, fx);
    node_d_offsets[1] = offsetof(node_t, fy);

    for (c = 0; c < job->num_colors; c++) {
        tc_idx = thread_id * job->num_colors + c;
        for (i = 0; i < job->particle_by_element_color_lengths[tc_idx]; i++) {
            p_idx = job->particle_by_element_color_lists[tc_idx][i];

            /*
                Note: We don't need to check if the particle is active since
                the list assembly already checks for this condition.
            */

            p = job->in_element[p_idx];

            s[0] = job->h1[p_idx];
            s[1] = job->h2[p_idx];
            s[2] = job->h3[p_idx];
            s[3] = job->h4[p_idx];

            /* Mass. */
            pdata[0] = job->particles[p_idx].m;

            /* Momentum. */
            pdata[1] = job->particles[p_idx].x_t * job->particles[p_idx].m;
            pdata[2] = job->particles[p_idx].y_t * job->particles[p_idx].m;

            /* Inertia. */
            pdata[3] = job->particles[p_idx].x_tt * job->particles[p_idx].m;
            pdata[4] = job->particles[p_idx].y_tt * job->particles[p_idx].m;

            /* Body forces. */
            pdata[5] = job->particles[p_idx].bx * job->particles[p_idx].m;
            pdata[6] = job->particles[p_idx].by * job->particles[p_idx].m;

            /* Stress. */
            stressdata[0] = -job->particles[p_idx].sxx * job->particles[p_idx].v;
            stressdata[1] = -job->particles[p_idx].sxy * job->particles[p_idx].v;
            stressdata[2] = -job->particles[p_idx].syy * job->particles[p_idx].v;

            accumulate_p_to_n_ds_list(job->nodes,
                node_field_offsets, job->elements[p].nodes, s, 
                NODES_PER_ELEMENT, pdata, pdata_len);

            ds[0] = job->b11[p_idx];
            ds[1] = job->b12[p_idx];
            ds[2] = job->b13[p_idx];
            ds[3] = job->b14[p_idx];

            accumulate_p_to_n_ds_list(job->nodes,
                node_d_offsets, job->elements[p].nodes, ds, 
                NODES_PER_ELEMENT, &(stressdata[0]), 2);

            ds[0] = job->b21[p_idx];
            ds[1] = job->b22[p_idx];
            ds[2] = job->b23[p_idx];
            ds[3] = job->b24[p_idx];

            accumulate_p_to_n_ds_list(job->nodes,
                node_d_offsets, job->elements[p].nodes, ds, 
                NODES_PER_ELEMENT, &(stressdata[1]), 2);
        }

        /*
            Each color can be done simultaneously, but we have to sync between
            colors.
        */
        pthread_barrier_wait(job->serialize_barrier);
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void map_to_grid_implicit(job_t *job)
{
    double mn;
    
    size_t i;
    double s[NODES_PER_ELEMENT];

    /* kinda a hack */
    const int pdata_len = 7;
    size_t node_field_offsets[pdata_len];
    double pdata[pdata_len];

    int p;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        s[0] = job->h1[i];
        s[1] = job->h2[i];
        s[2] = job->h3[i];
        s[3] = job->h4[i];

        p = job->in_element[i];

        /* Mass. */
        pdata[0] = job->particles[i].m;
        node_field_offsets[0] = offsetof(node_t, m);

        /* Momentum. */
        pdata[1] = job->particles[i].x_t * job->particles[i].m;
        pdata[2] = job->particles[i].y_t * job->particles[i].m;
        node_field_offsets[1] = offsetof(node_t, mx_t);
        node_field_offsets[2] = offsetof(node_t, my_t);

        /* Pseudoforce. */
        pdata[3] = job->particles[i].x_tt * job->particles[i].m;
        pdata[4] = job->particles[i].y_tt * job->particles[i].m;
        node_field_offsets[3] = offsetof(node_t, mx_tt);
        node_field_offsets[4] = offsetof(node_t, my_tt);

        /* Body forces. */
        pdata[5] = job->particles[i].bx * job->particles[i].m;
        pdata[6] = job->particles[i].by * job->particles[i].m;
        node_field_offsets[5] = offsetof(node_t, fx);
        node_field_offsets[6] = offsetof(node_t, fy);

        accumulate_p_to_n_ds_list(job->nodes,
            node_field_offsets, job->elements[p].nodes, s, 
            NODES_PER_ELEMENT, pdata, pdata_len);

    }

    /* for implicit method, map to another array */
    for (i = 0; i < job->num_nodes; i++) {
        mn = job->nodes[i].m;
        /* used lumped mass */
        job->m_grid[NODAL_DOF * i + XDOF_IDX] = mn;
        job->m_grid[NODAL_DOF * i + YDOF_IDX] = mn;

        /* need previous timestep's acceleration for implicit method */
        if (mn > TOL) {
            job->a_grid[NODAL_DOF * i + XDOF_IDX] = job->nodes[i].mx_tt / mn;
            job->a_grid[NODAL_DOF * i + YDOF_IDX] = job->nodes[i].my_tt  / mn;
        } else {
            job->a_grid[NODAL_DOF * i + XDOF_IDX] = 0;
            job->a_grid[NODAL_DOF * i + YDOF_IDX] = 0;
        }

        /* external forces */
        job->f_ext_grid[NODAL_DOF * i + XDOF_IDX] = job->nodes[i].fx;
        job->f_ext_grid[NODAL_DOF * i + YDOF_IDX] = job->nodes[i].fy;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void map_to_grid(job_t *job)
{
    double mn;
    size_t i;
    double s[NODES_PER_ELEMENT];

    int p;
/*    double *m_particles;*/

/*    m_particles = malloc(sizeof(double) * job->num_particles);*/
/*    for (i = 0; i < job->num_nodes; i++) {*/
/*        job->m_nodes[i] = 0;*/
/*    }*/
/*    for (i = 0; i < job->num_particles; i++) {*/
/*        m_particles[i] = job->particles[i].m;*/
/*    }*/

/*    cs_gaxpy(job->phi, m_particles, job->m_nodes);*/

/*    for (i = 0; i < job->num_nodes; i++) {*/
/*        job->nodes[i].m = job->m_nodes[i];*/
/*    }*/

    /* accumulate mass */
/*    map_particles_to_nodes_doublescalar(job, offsetof(node_t, m), offsetof(particle_t, m));*/

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        s[0] = job->h1[i];
        s[1] = job->h2[i];
        s[2] = job->h3[i];
        s[3] = job->h4[i];

        p = job->in_element[i];

        /* Mass. */
/*        ACCUMULATE(m,job,m,i,n,h);*/
        accumulate_p_to_n_doublescalar(job->nodes,
            offsetof(node_t, m), job->elements[p].nodes, s, 
            NODES_PER_ELEMENT, job->particles[i].m);

        /* Momentum. */
/*        ACCUMULATE_WITH_MUL(mx_t,job,x_t,i,n,h,job->particles[i].m);*/
/*        ACCUMULATE_WITH_MUL(my_t,job,y_t,i,n,h,job->particles[i].m);*/
        accumulate_p_to_n_doublescalar(job->nodes,
            offsetof(node_t, mx_t), job->elements[p].nodes, s,
            NODES_PER_ELEMENT, job->particles[i].x_t * job->particles[i].m);
        accumulate_p_to_n_doublescalar(job->nodes,
            offsetof(node_t, my_t), job->elements[p].nodes, s, 
            NODES_PER_ELEMENT, job->particles[i].y_t * job->particles[i].m);

        /* Pseudoforce. */
/*        ACCUMULATE_WITH_MUL(mx_tt,job,x_tt,i,n,h,job->particles[i].m);*/
/*        ACCUMULATE_WITH_MUL(my_tt,job,y_tt,i,n,h,job->particles[i].m);*/
        accumulate_p_to_n_doublescalar(job->nodes,
            offsetof(node_t, mx_tt), job->elements[p].nodes, s,
            NODES_PER_ELEMENT, job->particles[i].x_tt * job->particles[i].m);
        accumulate_p_to_n_doublescalar(job->nodes,
            offsetof(node_t, my_tt), job->elements[p].nodes, s, 
            NODES_PER_ELEMENT, job->particles[i].y_tt * job->particles[i].m);


        /* Body forces. */
/*        ACCUMULATE_WITH_MUL(fx,job,bx,i,n,h,job->particles[i].m);*/
/*        ACCUMULATE_WITH_MUL(fy,job,by,i,n,h,job->particles[i].m);*/
        accumulate_p_to_n_doublescalar(job->nodes,
            offsetof(node_t, fx), job->elements[p].nodes, s,
            NODES_PER_ELEMENT, job->particles[i].bx * job->particles[i].m);
        accumulate_p_to_n_doublescalar(job->nodes,
            offsetof(node_t, fy), job->elements[p].nodes, s, 
            NODES_PER_ELEMENT, job->particles[i].by * job->particles[i].m);

    }

    /* for implicit method, map to another array */
    for (i = 0; i < job->num_nodes; i++) {
        mn = job->nodes[i].m;
        /* used lumped mass */
        job->m_grid[NODAL_DOF * i + XDOF_IDX] = mn;
        job->m_grid[NODAL_DOF * i + YDOF_IDX] = mn;

        /* need previous timestep's acceleration for implicit method */
        if (mn > TOL) {
            job->a_grid[NODAL_DOF * i + XDOF_IDX] = job->nodes[i].mx_tt / mn;
            job->a_grid[NODAL_DOF * i + YDOF_IDX] = job->nodes[i].my_tt  / mn;
        } else {
            job->a_grid[NODAL_DOF * i + XDOF_IDX] = 0;
            job->a_grid[NODAL_DOF * i + YDOF_IDX] = 0;
        }

        /* external forces */
        job->f_ext_grid[NODAL_DOF * i + XDOF_IDX] = job->nodes[i].fx;
        job->f_ext_grid[NODAL_DOF * i + YDOF_IDX] = job->nodes[i].fy;
    }

/*    free(m_particles);*/

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void move_particles(job_t *job)
{
    int i;
    double dux, a_x_t;
    double duy, a_y_t;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        dux = (N_TO_P(job, ux, i));
        duy = (N_TO_P(job, uy, i));

        job->particles[i].x += dux;
        job->particles[i].y += duy;
        job->particles[i].ux += dux;
        job->particles[i].uy += duy;

        while(job->particles[i].x < 0) { job->particles[i].x += -floor(job->particles[i].x); }
        while(job->particles[i].x > 1) { job->particles[i].x -= floor(job->particles[i].x); }

        while(job->particles[i].y < 0) { job->particles[i].y += -floor(job->particles[i].y); }
        while(job->particles[i].y > 1) { job->particles[i].y -= floor(job->particles[i].y); }

        a_x_t = job->particles[i].x_tt;
        a_y_t = job->particles[i].y_tt;

        job->particles[i].x_tt = (N_TO_P(job, x_tt, i));
        job->particles[i].y_tt = (N_TO_P(job, y_tt, i));

        job->particles[i].x_t += 0.5 * job->dt * (job->particles[i].x_tt + a_x_t);
        job->particles[i].y_t += 0.5 * job->dt * (job->particles[i].y_tt + a_y_t);
/*        job->particles[i].x_t = (N_TO_P(job, x_t, i));*/
/*        job->particles[i].y_t = (N_TO_P(job, y_t, i));*/
    }
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void move_particles_explicit_usf(job_t *job)
{
    int i;
    double dux, duy;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        job->particles[i].x_tt = (N_TO_P(job, x_tt, i));
        job->particles[i].y_tt = (N_TO_P(job, y_tt, i));

        job->particles[i].x_t = (N_TO_P(job, x_t, i));
        job->particles[i].y_t = (N_TO_P(job, y_t, i));
        dux = job->dt * (N_TO_P(job, x_t, i));
        duy = job->dt * (N_TO_P(job, y_t, i));
        job->particles[i].x += dux;
        job->particles[i].y += duy;
        job->particles[i].ux += dux;
        job->particles[i].uy += duy;

        while(job->particles[i].x < 0) { job->particles[i].x += -floor(job->particles[i].x); }
        while(job->particles[i].x > 1) { job->particles[i].x -= floor(job->particles[i].x); }

        while(job->particles[i].y < 0) { job->particles[i].y += -floor(job->particles[i].y); }
        while(job->particles[i].y > 1) { job->particles[i].y -= floor(job->particles[i].y); }
    }
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void move_particles_explicit_usl(job_t *job)
{
    move_particles_explicit_usl_split(job, 0, job->num_particles);
    return;
}

void move_particles_explicit_usl_split(job_t *job, size_t p_start, size_t p_stop)
{
    size_t i, j, el, n;
    double dux, duy;
    double m, s[4];

    for (i = p_start; i < p_stop; i++) {
        CHECK_ACTIVE(job, i);

#if 1
/*        job->particles[i].x_tt = (N_TO_P(job, x_tt, i));*/
/*        job->particles[i].y_tt = (N_TO_P(job, y_tt, i));*/

        s[0] = job->h1[i];
        s[1] = job->h2[i];
        s[2] = job->h3[i];
        s[3] = job->h4[i];

        el = job->in_element[i];
        job->particles[i].x_tt = 0;
        job->particles[i].y_tt = 0;
        for (j = 0; j < 4; j++) {
            n = job->elements[el].nodes[j];
            m = job->nodes[n].m;
            job->particles[i].x_tt += ((s[j] * job->nodes[n].fx) / m);
            job->particles[i].y_tt += ((s[j] * job->nodes[n].fy) / m);
        }

        job->particles[i].x_t += job->dt * job->particles[i].x_tt;
        job->particles[i].y_t += job->dt * job->particles[i].y_tt;
/*        dux = job->dt * (N_TO_P(job, x_t, i));*/
/*        duy = job->dt * (N_TO_P(job, y_t, i));*/
        dux = 0;
        duy = 0;
        for (j = 0; j < 4; j++) {
            n = job->elements[el].nodes[j];
            m = job->nodes[n].m;
            dux += ((s[j] * job->nodes[n].mx_t) / m);
            duy += ((s[j] * job->nodes[n].my_t) / m);
        }
        dux *= job->dt;
        duy *= job->dt;
        job->particles[i].x += dux;
        job->particles[i].y += duy;
        job->particles[i].ux += dux;
        job->particles[i].uy += duy;
#endif

/*        job->particles[i].x_t = (N_TO_P(job, x_t, i)) + job->dt * job->particles[i].x_tt;*/
/*        job->particles[i].y_t = (N_TO_P(job, y_t, i)) + job->dt * job->particles[i].y_tt;*/

/*        job->particles[i].x_tt = (N_TO_P(job, x_tt, i));*/
/*        job->particles[i].y_tt = (N_TO_P(job, y_tt, i));*/

/*        dux = job->dt * (N_TO_P(job, x_t, i));*/
/*        duy = job->dt * (N_TO_P(job, y_t, i));*/
/*        job->particles[i].x += dux;*/
/*        job->particles[i].y += duy;*/
/*        job->particles[i].ux += dux;*/
/*        job->particles[i].uy += duy;*/

/*        dux = job->dt * (N_TO_P(job, x_t, i));*/
/*        duy = job->dt * (N_TO_P(job, y_t, i));*/
/*        job->particles[i].x_tt = (N_TO_P(job, x_tt, i));*/
/*        job->particles[i].y_tt = (N_TO_P(job, y_tt, i));*/

/*        job->particles[i].x += dux;*/
/*        job->particles[i].y += duy;*/
/*        job->particles[i].x_t += job->dt * job->particles[i].x_tt;*/
/*        job->particles[i].y_t += job->dt * job->particles[i].y_tt;*/

/*        job->particles[i].ux += dux;*/
/*        job->particles[i].uy += duy;*/


        while(job->particles[i].x < 0) { job->particles[i].x += -floor(job->particles[i].x); }
        while(job->particles[i].x > 1) { job->particles[i].x -= floor(job->particles[i].x); }

        while(job->particles[i].y < 0) { job->particles[i].y += -floor(job->particles[i].y); }
        while(job->particles[i].y > 1) { job->particles[i].y -= floor(job->particles[i].y); }

    }
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_deformation_gradient(job_t *job)
{
    int i;

    double delta;
    double delta_inv;
    double m11, m12, m21, m22;
    double a, b, c, d;
    double cd, sd, f;
    double dsq;
    double tmp_Fxx;
    double tmp_Fxy;
    double tmp_Fyx;
    double tmp_Fyy;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        a = job->dt * DX_N_TO_P(job, x_t, 1, i);
        b = job->dt * DY_N_TO_P(job, x_t, 1, i);
        c = job->dt * DX_N_TO_P(job, y_t, 1, i);
        d = job->dt * DY_N_TO_P(job, y_t, 1, i);

        dsq = (a - d) * (a - d) + 4 * b * c;
        delta = 0.5 * sqrt(fabs(dsq));
        delta_inv = (1.0f / delta);
        if (dsq == 0 || delta < 1e-15) {
            cd = 1;
            sd = 1;
            delta_inv = 1;
        } else if (dsq > 0) {
            cd = cosh(delta);
            sd = sinh(delta);
        } else {
            cd = cos(delta);
            sd = sin(delta);
        }

        m11 = cd + 0.5 * (a - d) * sd * delta_inv;
        m12 = b * sd * delta_inv;
        m21 = c * sd * delta_inv;
        m22 = cd - 0.5 * (a - d) * sd * delta_inv;

        f = exp(0.5 * (a + d));
        tmp_Fxx = f * (m11);
        tmp_Fxy = f * (m12);
        tmp_Fyx = f * (m21);
        tmp_Fyy = f * (m22);

        job->particles[i].Fxx = job->particles[i].Fxx * tmp_Fxx + job->particles[i].Fyx * tmp_Fxy;
        job->particles[i].Fxy = job->particles[i].Fxy * tmp_Fxx + job->particles[i].Fyy * tmp_Fxy;
        job->particles[i].Fyx = job->particles[i].Fxx * tmp_Fyx + job->particles[i].Fyx * tmp_Fyy;
        job->particles[i].Fyy = job->particles[i].Fxy * tmp_Fyx + job->particles[i].Fyy * tmp_Fyy;

/*            job->particles[i].v = (job->particles[i].Fxx * job->particles[i].Fyy - job->particles[i].Fxy * job->particles[i].Fyx) * job->particles[i].v0;*/
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_particle_vectors(job_t *job)
{
    int i;
    for (i = 0; i < job->num_particles; i++) {
        job->particles[i].r1[0] =
            job->particles[i].Fxx * job->particles[i].r1_initial[0] +
            job->particles[i].Fxy * job->particles[i].r1_initial[1];
        job->particles[i].r1[1] =
            job->particles[i].Fyx * job->particles[i].r1_initial[0] +
            job->particles[i].Fyy * job->particles[i].r1_initial[1];
        job->particles[i].r2[0] =
            job->particles[i].Fxx * job->particles[i].r2_initial[0] +
            job->particles[i].Fxy * job->particles[i].r2_initial[1];
        job->particles[i].r2[1] =
            job->particles[i].Fyx * job->particles[i].r2_initial[0] +
            job->particles[i].Fyy * job->particles[i].r2_initial[1];
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_corner_positions(job_t *job)
{
    int i, j, k;
    int e, n;

/*    FILE *fd;*/
/*    fd = fopen("pnodes.txt", "a");*/

/*    for (i = 0; i < (job->num_particles * job->num_nodes); i++) {*/
/*        job->phi[i] = 0;*/
/*    }*/

    for (i = 0; i < job->num_particles; i++) {
        /* Find corner positions from (adjusted) vectors. */
        job->particles[i].corners[2][S_XIDX] = job->particles[i].x +
            job->particles[i].r1[S_XIDX] + job->particles[i].r2[S_XIDX];
        job->particles[i].corners[2][S_YIDX] = job->particles[i].y +
            job->particles[i].r1[S_YIDX] + job->particles[i].r2[S_YIDX];

        job->particles[i].corners[3][S_XIDX] = job->particles[i].x -
            job->particles[i].r1[S_XIDX] + job->particles[i].r2[S_XIDX];
        job->particles[i].corners[3][S_YIDX] = job->particles[i].y -
            job->particles[i].r1[S_YIDX] + job->particles[i].r2[S_YIDX];

        job->particles[i].corners[0][S_XIDX] = job->particles[i].x -
            job->particles[i].r1[S_XIDX] - job->particles[i].r2[S_XIDX];
        job->particles[i].corners[0][S_YIDX] = job->particles[i].y -
            job->particles[i].r1[S_YIDX] - job->particles[i].r2[S_YIDX];

        job->particles[i].corners[1][S_XIDX] = job->particles[i].x +
            job->particles[i].r1[S_XIDX] - job->particles[i].r2[S_XIDX];
        job->particles[i].corners[1][S_YIDX] = job->particles[i].y +
            job->particles[i].r1[S_YIDX] - job->particles[i].r2[S_YIDX];

        for (j = 0; j < 4; j++) {
            /* Figure out which element each corner is in. */
            job->particles[i].corner_elements[j] = WHICH_ELEMENT(job->particles[i].corners[j][S_XIDX], job->particles[i].corners[j][S_YIDX], job->N, job->h);

            e = job->particles[i].corner_elements[j];
            n = job->elements[e].nodes[0];

/*            fprintf(fd, "%d:", i);*/
/*            for (k = 0; k < 4; k++) {*/
/*                fprintf(fd, " %d", job->elements[e].nodes[k]);*/
/*            }*/
/*            fprintf(fd, "\n");*/

            /* calculate local coordinates of corners */
            if (job->particles[i].corner_elements[j] != -1) {
                global_to_local_coords(
                    &(job->particles[i].cornersl[j][S_XIDX]),
                    &(job->particles[i].cornersl[j][S_YIDX]),
                    job->particles[i].corners[j][S_XIDX],
                    job->particles[i].corners[j][S_YIDX],
                    job->nodes[n].x,
                    job->nodes[n].y,
                    job->h
                );

                tent(
                    &(job->particles[i].sc[j][0]),
                    &(job->particles[i].sc[j][1]),
                    &(job->particles[i].sc[j][2]),
                    &(job->particles[i].sc[j][3]),
                    job->particles[i].cornersl[j][S_XIDX],
                    job->particles[i].cornersl[j][S_YIDX]
                );
            }
        }

        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {

                switch(k) {
                    case 0:
                        job->particles[i].grad_sc[k][j][S_XIDX] = (job->particles[i].sc[0][j])*(job->particles[i].r1[S_YIDX] - job->particles[i].r2[S_YIDX]);

                        job->particles[i].grad_sc[k][j][S_YIDX] = (job->particles[i].sc[0][j])*(job->particles[i].r2[S_XIDX] - job->particles[i].r1[S_XIDX]);
                    break;
                    case 1:
                        job->particles[i].grad_sc[k][j][S_XIDX] = (job->particles[i].sc[1][j])*(job->particles[i].r1[S_YIDX] + job->particles[i].r2[S_YIDX]);

                        job->particles[i].grad_sc[k][j][S_YIDX] = (job->particles[i].sc[1][j])*(-job->particles[i].r1[S_XIDX] - job->particles[i].r2[S_XIDX]);
                    break;
                    case 2:
                        job->particles[i].grad_sc[k][j][S_XIDX] = (-job->particles[i].sc[2][j])*(job->particles[i].r1[S_YIDX] - job->particles[i].r2[S_YIDX]);

                        job->particles[i].grad_sc[k][j][S_YIDX] = (-job->particles[i].sc[2][j])*(job->particles[i].r2[S_XIDX] - job->particles[i].r1[S_XIDX]);
                    break;
                    case 3:
                        job->particles[i].grad_sc[k][j][S_XIDX] = (-job->particles[i].sc[3][j])*(job->particles[i].r1[S_YIDX] + job->particles[i].r2[S_YIDX]);

                        job->particles[i].grad_sc[k][j][S_YIDX] = (-job->particles[i].sc[3][j])*(-job->particles[i].r1[S_XIDX] - job->particles[i].r2[S_XIDX]);
                    break;
                }
/*                job->particles[i].grad_sc[k][j][S_XIDX] = (job->particles[i].sc[0][j] - job->particles[i].sc[2][j])*(job->particles[i].r1[S_YIDX] - job->particles[i].r2[S_YIDX]) +*/
/*                (job->particles[i].sc[1][j] - job->particles[i].sc[3][j])*(job->particles[i].r1[S_YIDX] + job->particles[i].r2[S_YIDX]);*/

/*                job->particles[i].grad_sc[k][j][S_YIDX] = (job->particles[i].sc[0][j] - job->particles[i].sc[2][j])*(job->particles[i].r2[S_XIDX] - job->particles[i].r1[S_XIDX]) +*/
/*                (job->particles[i].sc[1][j] - job->particles[i].sc[3][j])*(-job->particles[i].r1[S_XIDX] - job->particles[i].r2[S_XIDX]);*/

/*                printf("%d,%d sc [%g %g %g %g]\n", j, k, job->particles[i].sc[k][0], job->particles[i].sc[k][1], job->particles[i].sc[k][2], job->particles[i].sc[k][3]);*/
/*                printf("%d,%d sc [%g %g]\n", j, k, job->particles[i].grad_sc[k][j][S_XIDX],*/
/*                    job->particles[i].grad_sc[k][j][S_YIDX]);*/

                /* should be averaging domain, but use particle volume for now */
                job->particles[i].grad_sc[k][j][S_XIDX] *= (1.0 / job->particles[i].v);
                job->particles[i].grad_sc[k][j][S_YIDX] *= (1.0 / job->particles[i].v);

                if (job->particles[i].corner_elements[k] == -1) {
                    continue;
                }

                /* phi rows are nodal indicies, cols are particle indicies */
/*                job->phi[job->num_particles * job->elements[job->particles[i].corner_elements[k]].nodes[j] + i] += 0.25 * job->particles[i].sc[k][j];*/
            }
        }

    }

/*    for (i = 0; i < job->num_nodes; i++) {*/
/*        for (j = 0; j < job->num_particles; j++) {*/
/*            fprintf(fd, " %g", job->phi[job->num_particles * i + j]);*/
/*        }*/
/*        fprintf(fd, "\n");*/
/*    }*/
/*    fprintf(fd, "\n");*/

/*    fclose(fd);*/

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_particle_densities(job_t *job)
{
    update_particle_densities_split(job, 0, job->num_particles);
    return;
}

void update_particle_densities_split(job_t *job, size_t p_start, size_t p_stop)
{
    size_t i;
#if 0
    int j;

    double left, right, top, bottom;
    double h_spacing, v_spacing;

    for (i = 0; i < job->num_elements; i++) {
        job->elements[i].filled = 0;
        job->elements[i].n = 0;
        job->elements[i].m = 0;
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);
        job->elements[job->in_element[i]].filled = 1;
        job->elements[job->in_element[i]].n++;
        job->elements[job->in_element[i]].m += job->particles[i].m;
        
    }

    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].num_filled_element_neighbors = 0;
        job->nodes[i].mass_filled_element_neighbors = 0;
    }

    for (i = 0; i < job->num_elements; i++) {
        if (job->elements[i].filled) {
            for (j = 0; j < 4; j++) {
                job->nodes[job->elements[i].nodes[j]].num_filled_element_neighbors++;
                job->nodes[job->elements[i].nodes[j]].mass_filled_element_neighbors += job->elements[i].m;
            }
        }
    }

    for (i = 0; i < job->num_elements; i++) {
        h_spacing = 2.0f;
        v_spacing = 2.0f;

        if (job->elements[i].neighbors[0] == -1) {
            right = job->elements[i].n;
            h_spacing = 1.0f;
        } else {
            right = job->elements[job->elements[i].neighbors[0]].n;
        }

        if (job->elements[i].neighbors[4] == -1) {
            left = job->elements[i].n;
            h_spacing = 1.0f;
        } else {
            left = job->elements[job->elements[i].neighbors[4]].n;
        }

        if (job->elements[i].neighbors[2] == -1) {
            top = job->elements[i].n;
            v_spacing = 1.0f;
        } else {
            top = job->elements[job->elements[i].neighbors[2]].n;
        }

        if (job->elements[i].neighbors[6] == -1) {
            bottom = job->elements[i].n;
            v_spacing = 1.0f;
        } else {
            bottom = job->elements[job->elements[i].neighbors[6]].n;
        }

        job->elements[i].grad_x = (right - left) / h_spacing;
        job->elements[i].grad_y = (top - bottom) / v_spacing;

        job->elements[i].grad_mag = hypot(job->elements[i].grad_x, job->elements[i].grad_y);

        /* point normal outward from body (negative of gradient). */
        if (job->elements[i].grad_mag > 0) {
            job->elements[i].n_x = -job->elements[i].grad_x / job->elements[i].grad_mag;
            job->elements[i].n_y = -job->elements[i].grad_y / job->elements[i].grad_mag;
            job->elements[i].n_theta = atan2(job->elements[i].n_y, job->elements[i].n_x);
        }

        /* printing */
/*        if (job->elements[i].grad_mag > 1.0f) {*/
/*            printf("%d: %f %f %f %f\n", i, left, right, top, bottom);*/
/*            printf("%d: %f\n", i, job->elements[i].grad_mag);*/
/*            printf("%d: %f %f\n", i, job->elements[i].n_x, job->elements[i].n_y);*/
/*        }*/

    }

        #define GRAD_THRESHOLD 0.50f

    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].rho = job->nodes[i].mass_filled_element_neighbors / (job->nodes[i].num_filled_element_neighbors * job->h * job->h);
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        /* Fixed version */
/*        if (job->elements[job->in_element[i]].grad_mag > GRAD_THRESHOLD) {*/
/*            job->particles[i].v = job->particles[i].v * exp (job->dt * (job->particles[i].exx_t + job->particles[i].eyy_t));*/
/*        } else {*/
/*            job->particles[i].v = job->h * job->h / (job->elements[job->in_element[i]].n);*/
/*        }*/

        if (job->use_cpdi) {
            job->particles[i].v = job->particles[i].v0 *
                ((job->particles[i].Fxx * job->particles[i].Fyy) - 
                job->particles[i].Fxy * job->particles[i].Fyx);
        } else {
            job->particles[i].v = job->h * job->h / (job->elements[job->in_element[i]].n);
        }
    }
#endif


    if (job->use_cpdi) {
        for (i = p_start; i < p_stop; i++) {
            CHECK_ACTIVE(job, i);
            job->particles[i].v = job->particles[i].v0 *
                ((job->particles[i].Fxx * job->particles[i].Fyy) - 
                job->particles[i].Fxy * job->particles[i].Fyx);
        }
    } else {
        for (i = p_start; i < p_stop; i++) {
            CHECK_ACTIVE(job, i);
/*            job->particles[i].v = job->h * job->h / (job->elements[job->in_element[i]].n);*/
            job->particles[i].v = job->particles[i].v * exp (job->dt * (job->particles[i].exx_t + job->particles[i].eyy_t));
#if 0
            size_t j, el, n, r, c, nvols;
            double rho, s[4];
            el = job->in_element[i];
            rho = 0;
            s[0] = job->h1[i];
            s[1] = job->h2[i];
            s[2] = job->h3[i];
            s[3] = job->h4[i];
            for (j = 0; j < 4; j++) {
                nvols = 4;
                n = job->elements[el].nodes[j];
                r = n / job->N;
                c = n % job->N;
                if (r == 0 || r == (job->N - 1)) {
                    nvols -= 2;
                } 
                if (c == 0 || c == (job->N / 2) || c == (job->N - 1)) {
                    nvols -= 2;
                }
                if (nvols == 0) {
                    nvols = 1;
                }
                rho += s[j] * job->nodes[n].m / (0.25 * nvols * job->h * job->h);
            }
/*            printf("rho = %f, ", rho);*/
            job->particles[i].v = job->particles[i].m / rho;
#endif
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void generate_mappings(job_t *job)
{
    size_t p_idx;
    size_t e_idx;
    size_t i_idx;

    size_t i;

    size_t nr;
    size_t nc;
    size_t nnz;
    cs *triplets;

    double s[4];

    int res;

    nr = job->num_nodes;
    nc = job->num_particles;
    /* Each particle only affects one element at most. */
    nnz = 4 * job->num_particles;

    triplets = cs_spalloc(nr, nc, nnz, 1, 1);

    for (p_idx = 0; p_idx < job->num_particles; p_idx++) {
        if (job->active[p_idx] == 0) {
            continue;
        }

        e_idx = job->in_element[p_idx];
        tent(&(s[0]), &(s[1]), &(s[2]), &(s[3]),
            job->particles[p_idx].xl, job->particles[p_idx].yl);
        for (i = 0; i < 4; i++) {
            i_idx = job->elements[e_idx].nodes[i];

            res = cs_entry(triplets, i_idx, p_idx, s[i]);
            if (res == 0) {
                fprintf(stderr, "Can't add entry to generate mappings.\n");
                exit(EXIT_ERROR_CS_ENTRY);
            }
        }
    }

    if (job->phi != NULL) {
        cs_spfree(job->phi);
    }
    if (job->phi_transpose != NULL) {
        cs_spfree(job->phi_transpose);
    }

    job->phi = cs_compress(triplets);
    cs_dupl(job->phi);
    job->phi_transpose = cs_transpose(job->phi, 1); /* copy values. */
    cs_spfree(triplets);

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void mpm_cleanup(job_t *job)
{
    free(job->nodes);
    free(job->particles);
    free(job->elements);
    
    free(job->in_element);

/*    free(job->m_nodes);*/
/*    free(job->mx_t_nodes);*/
/*    free(job->my_t_nodes);*/
/*    free(job->mx_tt_nodes);*/
/*    free(job->my_tt_nodes);*/
/*    free(job->x_t_nodes);*/
/*    free(job->y_t_nodes);*/
/*    free(job->x_tt_nodes);*/
/*    free(job->y_tt_nodes);*/
/*    free(job->fx_nodes);*/
/*    free(job->fy_nodes);*/

    free(job->h1);
    free(job->h2);
    free(job->h3);
    free(job->h4);

    free(job->b11);
    free(job->b12);
    free(job->b13);
    free(job->b14);

    free(job->b21);
    free(job->b22);
    free(job->b23);
    free(job->b24);

    free(job);

    return;
}
/*----------------------------------------------------------------------------*/

