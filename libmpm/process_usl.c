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
#include "process_usl.h"
#include "material.h"
#include "exitcodes.h"

#include "spmd.h"

#include <suitesparse/cs.h>

#include <assert.h>

#define TOL 5e-11

#define ijton(i,j,N) ((j)*(N) + (i))

#define WHICH_ELEMENT WHICH_ELEMENT4

/* XXX: ugly, fix soon */
#define WHICH_ELEMENT4(xp,yp,N,Lx,Ly,hx,hy) \
    ((int)(((xp)<Lx && (xp)>=0.0 && (yp)<Ly && (yp)>=0.0)?((floor((xp)/(hx)) + floor((yp)/(hy))*((N)-1))):(-1)))

#define CHECK_ACTIVE(j,i) if (j->active[i] == 0) { continue; }

/*----------------------------------------------------------------------------*/
job_t *mpm_init(int N, double hx, double hy, double Lx, double Ly, particle_t *particles, size_t num_particles, double t)
{
    int n;

    job_t *job;

    job = (job_t *)malloc(sizeof(job_t));

    job->t = 0;
    job->t_stop = t;
    job->frame = 0;

    job->N = N;
    job->hx = hx;
    job->hy = hy;
    job->Lx = Lx;
    job->Ly = Ly;
    job->num_nodes = N*N;
    job->num_particles = num_particles;

    job->num_elements = (N - 1) * (N - 1);

    /* Copy particles from given ICs. */
    job->particles = (particle_t *)calloc(num_particles, sizeof(particle_t));
    memcpy(job->particles, particles, num_particles * sizeof(particle_t));

    /* Set stress, strain to zero. */
    for (size_t i = 0; i < job->num_particles; i++) {
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

        job->particles[i].state[9] = 0;
        job->particles[i].state[10] = 0;

        job->particles[i].corners[0][0] = 0;
        job->particles[i].corners[0][1] = 0;
        job->particles[i].corners[1][0] = 0;
        job->particles[i].corners[1][1] = 0;
        job->particles[i].corners[2][0] = 0;
        job->particles[i].corners[2][1] = 0;
        job->particles[i].corners[3][0] = 0;
        job->particles[i].corners[3][1] = 0;
        job->particles[i].color = 0;
    }

    /* hack it in */
    assert(hx == hy);
    const double h = hx;

    /* Get node coordinates. */
    job->nodes = (node_t *)calloc(job->num_nodes, sizeof(node_t));
    for (size_t i = 0; i < job->num_nodes; i++) {
        node_number_to_coords(&(job->nodes[i].x), &(job->nodes[i].y), i, N, h);
        /* printf("preprocessing: xn=%g yn=%g\n", job->nodes[i].x, job->nodes[i].y); */
    }

    /* Get node numbering for elements. */
    job->elements = (element_t *)calloc(job->num_elements, sizeof(element_t));
    for (size_t i = 0; i < job->num_elements; i++) {
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
    const size_t Ne = job->N - 1;
    for (size_t i = 0; i < job->num_elements; i++) {
        size_t r = i / (job->N - 1);
        size_t c = i % (job->N - 1);

        job->elements[i].color = 2 * (r % 2)  + ((c % 2) + 1) - 1;
/*        job->color_list_lengths[job->elements[i].color % job->num_colors]++;*/
/*        fprintf(stderr, "job->elements[%d].color = %d\n", i, job->elements[i].color);*/

        /*
            Fill neighbor array starting with element to the right and
            moving in a postive direction (ccw). If there is no neighbor,
            set that neighbor to (-1).
        */
        for (size_t j = 0; j < 8; j++) {
            job->elements[i].neighbors[j] = -1;
        }

        if (r != (Ne - 1)) {
            job->elements[i].neighbors[2] = ijton(r+1, c, Ne);
            if (c != (Ne - 1)) {
                job->elements[i].neighbors[1] = ijton(r+1, c+1, Ne);
            }
            if (c != 0) {
                job->elements[i].neighbors[3] = ijton(r+1, c-1, Ne);
            }
        }

        if (r != 0) {
            job->elements[i].neighbors[6] = ijton(r-1, c, Ne);
            if (c != (Ne - 1)) {
                job->elements[i].neighbors[7] = ijton(r-1, c+1, Ne);
            }
            if (c != 0) {
                job->elements[i].neighbors[5] = ijton(r-1, c-1, Ne);
            }
        }

        if (c != (Ne - 1)) {
            job->elements[i].neighbors[0] = ijton(r, c+1, Ne);
        }

        if (c != 0) {
            job->elements[i].neighbors[4] = ijton(r, c-1, Ne);
        }

    }
    for (size_t i = 0; i < job->num_colors; i++) {
        job->color_indices[i] = 0;
    }

    /* Allocate space for tracking element->particle map. */
    job->in_element = (int *)malloc(job->num_particles * sizeof(int));

    /* Allocate space for map of active particles. */
    job->active = (int *)malloc(job->num_particles * sizeof(int));

    /* max size of u_grid is NODAL_DOF * number of nodes. */
    job->vec_len = NODAL_DOF * job->num_nodes;

    /* for dirichlet BCs */
    job->u_dirichlet = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->u_dirichlet_mask = (int *)malloc(NODAL_DOF * sizeof(int) * job->num_nodes);

    /* for periodic BCs */
    job->node_number_override = (int *)malloc(job->vec_len * sizeof(int));

    for (size_t i = 0; i < job->num_particles; i++) {
        job->in_element[i] = WHICH_ELEMENT(
            job->particles[i].x, job->particles[i].y, job->N, job->Lx, job->Ly, job->hx, job->hy);
        if (job->in_element[i] < 0 || (size_t)job->in_element[i] >= job->num_elements) {
            job->active[i] = 0;
        } else {
            job->active[i] = 1;
        }
    }

    /* set particle domains */
    for (size_t i = 0; i < job->num_particles; i++) {
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
    job->dt = 0.4 * h * sqrt(job->particles[0].m/(job->particles[0].v * EMOD));

    /* Don't use cpdi here. */
    job->use_cpdi = 0;

    /* swap this out later with a pointer from dlopen if needed. */
    job->material.num_fp64_props = 0;
    job->material.num_int_props = 0;
    job->material.fp64_props = NULL;
    job->material.int_props = NULL;
    job->material.material_filename = NULL;
    job->material.material_init = &material_init_linear_elastic;
    job->material.calculate_stress = &calculate_stress_linear_elastic;

    /* used to vary loads/bcs */
    job->step_number = 0;
    job->step_start_time = job->t;

    if (job->use_cpdi == 1) {
    } else {
        // 4 entries per particle (1 pt in element, 4 nodes per element).
        const size_t nnz = NODES_PER_ELEMENT * job->num_particles;
        const size_t rows = job->num_particles;
        const size_t cols = job->num_nodes;
        job->phi = spmd_create(nnz, rows, cols);

        // derivative matrices have same structure
        job->dphi_x = spmd_create(nnz, rows, cols);
        job->dphi_y = spmd_create(nnz, rows, cols);
    }

    return job;
}
/*----------------------------------------------------------------------------*/

// XXX: Phi operates from NODES to PARTICLES
// this is the TRANSPOSE of the usual definition
// use spmd_gatxpy to go the other way...
void setup_phi(job_t *job);
void setup_dphi_x(job_t *job);
void setup_dphi_y(job_t *job);

void setup_phi(job_t *job)
{
    struct sparsematrix_double *sp = job->phi;
    size_t n = 0;
    for (size_t i = 0; i < job->num_particles; i++) {
        sp->row_pointer[i] = n;

        const int el = job->in_element[i];
        if (el == -1) {
            continue;
        }

        const int *nn = job->elements[el].nodes;
        for (size_t j = 0; j < NODES_PER_ELEMENT; j++) {
            sp->vals[n+j] = job->particles[i].s[j];
            sp->column_index[n+j] = nn[j];
        }
        n += NODES_PER_ELEMENT;
    }
    sp->nnz = n;
    sp->row_pointer[sp->rows] = n;
    return;
}

void setup_dphi_x(job_t *job)
{
    struct sparsematrix_double *sp = job->dphi_x;
    size_t n = 0;
    for (size_t i = 0; i < job->num_particles; i++) {
        sp->row_pointer[i] = n;

        const int el = job->in_element[i];
        if (el == -1) {
            continue;
        }

        const int *nn = job->elements[el].nodes;
        for (size_t j = 0; j < NODES_PER_ELEMENT; j++) {
            sp->vals[n+j] = job->particles[i].grad_s[j][0];
            sp->column_index[n+j] = nn[j];
        }
        n += NODES_PER_ELEMENT;
    }
    sp->nnz = n;
    sp->row_pointer[sp->rows] = n;
    return;
}

void setup_dphi_y(job_t *job)
{
    struct sparsematrix_double *sp = job->dphi_y;
    size_t n = 0;
    for (size_t i = 0; i < job->num_particles; i++) {
        sp->row_pointer[i] = n;

        const int el = job->in_element[i];
        if (el == -1) {
            continue;
        }

        const int *nn = job->elements[el].nodes;
        for (size_t j = 0; j < NODES_PER_ELEMENT; j++) {
            sp->vals[n+j] = job->particles[i].grad_s[j][1];
            sp->column_index[n+j] = nn[j];
        }
        n += NODES_PER_ELEMENT;
    }
    sp->nnz = n;
    sp->row_pointer[sp->rows] = n;
    return;
}
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
        job->nodes[i].x_t = 0;
        job->nodes[i].y_t = 0;
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
        setup_phi(job);
        setup_dphi_x(job);
        setup_dphi_y(job);
        // spmd_print(job->phi, 1);
        // spmd_print(job->dphi_x, 1);
        // spmd_print(job->dphi_y, 1);

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
        (*(job->boundary.bc_time_varying))(job);
    }

    pthread_barrier_wait(job->serialize_barrier);
    /* Map particle state to grid quantites. */
    map_to_grid_explicit_split(job, task->id);

/* start DCA algorithm */
    #if 0
    rc = pthread_barrier_wait(job->serialize_barrier);
    for (size_t i = n_start; i < n_stop; i++) {
        job->nodes[i].sum_sqrt_m_neighbors = sqrt(job->nodes[i].m);
        job->nodes[i].max_m_neighbors = 0;
        const size_t r = i / job->N;
        const size_t c = i % job->N;
        size_t nidx = 0;
        const size_t max_N_idx = job->N - 1;
        if (r != max_N_idx) {
            nidx = ijton(r+1, c, job->N);
            job->nodes[i].sum_sqrt_m_neighbors += sqrt(job->nodes[nidx].m); 
            if (job->nodes[nidx].m > job->nodes[i].max_m_neighbors) {
                job->nodes[i].max_m_neighbors = job->nodes[nidx].m;
            }
            if (c != max_N_idx) {
                nidx = ijton(r+1, c+1, job->N);
                job->nodes[i].sum_sqrt_m_neighbors += sqrt(job->nodes[nidx].m);
                if (job->nodes[nidx].m > job->nodes[i].max_m_neighbors) {
                    job->nodes[i].max_m_neighbors = job->nodes[nidx].m;
                }
            }
            if (c != 0) {
                nidx = ijton(r+1, c-1, job->N);
                job->nodes[i].sum_sqrt_m_neighbors += sqrt(job->nodes[nidx].m);
                if (job->nodes[nidx].m > job->nodes[i].max_m_neighbors) {
                    job->nodes[i].max_m_neighbors = job->nodes[nidx].m;
                }
            }
        }

        if (r != 0) {
            nidx = ijton(r-1, c, job->N);
            job->nodes[i].sum_sqrt_m_neighbors += sqrt(job->nodes[nidx].m); 
            if (job->nodes[nidx].m > job->nodes[i].max_m_neighbors) {
                job->nodes[i].max_m_neighbors = job->nodes[nidx].m;
            }
            if (c != max_N_idx) {
                nidx = ijton(r-1, c+1, job->N);
                job->nodes[i].sum_sqrt_m_neighbors += sqrt(job->nodes[nidx].m);
                if (job->nodes[nidx].m > job->nodes[i].max_m_neighbors) {
                    job->nodes[i].max_m_neighbors = job->nodes[nidx].m;
                }
            }
            if (c != 0) {
                nidx = ijton(r-1, c-1, job->N);
                job->nodes[i].sum_sqrt_m_neighbors += sqrt(job->nodes[nidx].m);
                if (job->nodes[nidx].m > job->nodes[i].max_m_neighbors) {
                    job->nodes[i].max_m_neighbors = job->nodes[nidx].m;
                }
            }
        }

        if (c != max_N_idx) {
            nidx = ijton(r, c+1, job->N);
            job->nodes[i].sum_sqrt_m_neighbors += sqrt(job->nodes[nidx].m); 
            if (job->nodes[nidx].m > job->nodes[i].max_m_neighbors) {
                job->nodes[i].max_m_neighbors = job->nodes[nidx].m;
            }
        }

        if (c != 0) {
            nidx = ijton(r, c-1, job->N);
            job->nodes[i].sum_sqrt_m_neighbors += sqrt(job->nodes[nidx].m); 
            if (job->nodes[nidx].m > job->nodes[i].max_m_neighbors) {
                job->nodes[i].max_m_neighbors = job->nodes[nidx].m;
            }
        }
    }
    rc = pthread_barrier_wait(job->serialize_barrier);
    if (rc == PTHREAD_BARRIER_SERIAL_THREAD) {
        for (size_t i = 0; i < job->num_nodes; i++) {
            if ((job->nodes[i].m == 0) || (job->nodes[i].m > 0.04*job->nodes[i].max_m_neighbors)) {
                continue;
            }
            const double denom = job->nodes[i].sum_sqrt_m_neighbors;
            const size_t r = i / job->N;
            const size_t c = i % job->N;
            size_t nidx = i;
            double f = 0;
            const size_t max_N_idx = job->N - 1;
            if (r != max_N_idx) {
                nidx = ijton(r+1, c, job->N);
                f = sqrt(job->nodes[nidx].m) / denom;
                job->nodes[nidx].fx += f * job->nodes[i].fx;
                job->nodes[nidx].fy += f * job->nodes[i].fy;
                if (c != max_N_idx) {
                    nidx = ijton(r+1, c+1, job->N);
                    f = sqrt(job->nodes[nidx].m) / denom;
                    job->nodes[nidx].fx += f * job->nodes[i].fx;
                    job->nodes[nidx].fy += f * job->nodes[i].fy;
                }
                if (c != 0) {
                    nidx = ijton(r+1, c-1, job->N);
                    f = sqrt(job->nodes[nidx].m) / denom;
                    job->nodes[nidx].fx += f * job->nodes[i].fx;
                    job->nodes[nidx].fy += f * job->nodes[i].fy;
                }
            }

            if (r != 0) {
                nidx = ijton(r-1, c, job->N);
                f = sqrt(job->nodes[nidx].m) / denom;
                job->nodes[nidx].fx += f * job->nodes[i].fx;
                job->nodes[nidx].fy += f * job->nodes[i].fy;
                if (c != max_N_idx) {
                    nidx = ijton(r-1, c+1, job->N);
                    f = sqrt(job->nodes[nidx].m) / denom;
                    job->nodes[nidx].fx += f * job->nodes[i].fx;
                    job->nodes[nidx].fy += f * job->nodes[i].fy;
                }
                if (c != 0) {
                    nidx = ijton(r-1, c-1, job->N);
                    f = sqrt(job->nodes[nidx].m) / denom;
                    job->nodes[nidx].fx += f * job->nodes[i].fx;
                    job->nodes[nidx].fy += f * job->nodes[i].fy;
                }
            }

            if (c != max_N_idx) {
                nidx = ijton(r, c+1, job->N);
                f = sqrt(job->nodes[nidx].m) / denom;
                job->nodes[nidx].fx += f * job->nodes[i].fx;
                job->nodes[nidx].fy += f * job->nodes[i].fy;
            }

            if (c != 0) {
                nidx = ijton(r, c-1, job->N);
                f = sqrt(job->nodes[nidx].m) / denom;
                job->nodes[nidx].fx += f * job->nodes[i].fx;
                job->nodes[nidx].fy += f * job->nodes[i].fy;
            }
            f = sqrt(job->nodes[i].m) / denom;
            job->nodes[i].fx *= f;
            job->nodes[i].fy *= f;
        }
    }
    #endif
/* end DCA algorithm*/

    rc = pthread_barrier_wait(job->serialize_barrier);
    if (rc == PTHREAD_BARRIER_SERIAL_THREAD) {
        /* Zero perpendicular momentum at edge nodes. */
        (*(job->boundary.bc_momentum))(job);

        /*
            This is poorly named -- it actually calculates diverence of stress to
            get nodal force.
        */
        /* XXX: moved to map_to_grid_explicit_split. */
/*        update_stress(job);*/

        /* Zero perpendicular forces at edge nodes. */
        (*(job->boundary.bc_force))(job);

    }

    pthread_barrier_wait(job->serialize_barrier);

    /*
        Update momentum and velocity at nodes.
        Wait for all threads to finish updating nodes before calculating
        strainrates.
    */
    move_grid_split(job, n_start, n_stop); 
    pthread_barrier_wait(job->serialize_barrier);

    /* Update particle position and velocity. */
    move_particles_explicit_usl_split(job, p_start, p_stop);

    /* Calculate strain rate. */
    calculate_strainrate_split(job, p_start, p_stop);

    /* update volume */
    update_particle_densities_split(job, p_start, p_stop);

    /* Calculate stress. */
    (*(job->material.calculate_stress_threaded))(task);

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void move_grid_split(job_t *job, size_t n_start, size_t n_stop)
{
    for (size_t i = n_start; i < n_stop; i++) {
        const double m = job->nodes[i].m;

        if (m > TOL) {

            job->nodes[i].mx_t += job->dt * job->nodes[i].fx;
            job->nodes[i].my_t += job->dt * job->nodes[i].fy;

            job->nodes[i].x_t = job->nodes[i].mx_t / m;
            job->nodes[i].y_t = job->nodes[i].my_t / m;

            job->nodes[i].ux = job->dt * job->nodes[i].x_t;
            job->nodes[i].uy = job->dt * job->nodes[i].y_t;
        } else {
            // job->nodes[i].mx_t = 0;
            // job->nodes[i].my_t = 0;
            job->nodes[i].x_t = 0;
            job->nodes[i].y_t = 0;
            job->nodes[i].ux = 0;
            job->nodes[i].uy = 0;
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
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
            job->particles[i].x, job->particles[i].y, job->N, job->Lx, job->Ly, job->hx, job->hy);

        if (p != job->in_element[i]) {
            changed = 1;
            /*
            fprintf(job->output.log_fd, 
                "[%g] Particle %zu @(%g, %g) left element %d, now in element %d.\n",
                job->t, i, job->particles[i].x, job->particles[i].y,
                job->in_element[i], p);
            */
        }

        /* Update particle element. */
        job->in_element[i] = p;

        if (p == -1) {
            /*
            fprintf(job->output.log_fd,
                "[%g] Particle %zu outside of grid (%g, %g), marking as inactive.\n",
                job->t, i, job->particles[i].x, job->particles[i].y);
            */
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

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_shapefunctions_split(job_t *job, size_t p_start, size_t p_stop)
{
    /* assert(hx == hy) before this. (hack) */
    const double hx = job->hx;
    const double hy = job->hy;

    for (size_t i = p_start; i < p_stop; i++) {
        CHECK_ACTIVE(job, i);
        int p = job->in_element[i];
        const int *nn = job->elements[p].nodes;
        const double x1 = job->nodes[nn[0]].x;
        const double y1 = job->nodes[nn[0]].y;
        const double x2 = job->nodes[nn[1]].x;
        const double y2 = job->nodes[nn[1]].y;
        const double x3 = job->nodes[nn[2]].x;
        const double y3 = job->nodes[nn[2]].y;
        const double x4 = job->nodes[nn[3]].x;
        const double y4 = job->nodes[nn[3]].y;
        const double xp = job->particles[i].x;
        const double yp = job->particles[i].y;

        double h1 = 0;
        double h2 = 0;
        double h3 = 0;
        double h4 = 0;

        double h1x = 0;
        double h2x = 0;
        double h3x = 0;
        double h4x = 0;

        double h1y = 0;
        double h2y = 0;
        double h3y = 0;
        double h4y = 0;
        tent(&h1, &h2, &h3, &h4,
            x1, y1,
            x2, y2,
            x3, y3,
            x4, y4,
            hx, hy,
            xp, yp);
        grad_tent(
            &h1x, &h2x, &h3x, &h4x,
            &h1y, &h2y, &h3y, &h4y,
            x1, y1,
            x2, y2,
            x3, y3,
            x4, y4,
            hx, hy,
            xp, yp);

        job->particles[i].s[0] = h1;
        job->particles[i].s[1] = h2;
        job->particles[i].s[2] = h3;
        job->particles[i].s[3] = h4;

        job->particles[i].grad_s[0][0] = h1x;
        job->particles[i].grad_s[1][0] = h2x;
        job->particles[i].grad_s[2][0] = h3x;
        job->particles[i].grad_s[3][0] = h4x;

        job->particles[i].grad_s[0][1] = h1y;
        job->particles[i].grad_s[1][1] = h2y;
        job->particles[i].grad_s[2][1] = h3y;
        job->particles[i].grad_s[3][1] = h4y;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_strainrate_split(job_t *job, size_t p_start, size_t p_stop)
{
    if (p_start != 0) {
        return;
    }

    double *pvec = calloc(job->num_particles, sizeof(double));
    double *nvec = calloc(job->num_nodes, sizeof(double));

    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->nodes[i].m != 0) {
            nvec[i] = job->nodes[i].x_t;
        } else {
            nvec[i] = 0;
        }
    }
    spmd_gaxpy(job->dphi_x, nvec, pvec);
    for (size_t i = 0; i < job->num_particles; i++) {
        job->particles[i].L[0]= pvec[i];
        pvec[i] = 0;
    }
    spmd_gaxpy(job->dphi_y, nvec, pvec);
    for (size_t i = 0; i < job->num_particles; i++) {
        job->particles[i].L[1]= pvec[i];
        pvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->nodes[i].m != 0) {
            nvec[i] = job->nodes[i].y_t;
        } else {
            nvec[i] = 0;
        }
    }
    spmd_gaxpy(job->dphi_x, nvec, pvec);
    for (size_t i = 0; i < job->num_particles; i++) {
        job->particles[i].L[2] = pvec[i];
        pvec[i] = 0;
    }
    spmd_gaxpy(job->dphi_y, nvec, pvec);
    for (size_t i = 0; i < job->num_particles; i++) {
        job->particles[i].L[3] = pvec[i];
        pvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        job->particles[i].exx_t = job->particles[i].L[0];
        job->particles[i].exy_t = 0.5*(job->particles[i].L[1] + job->particles[i].L[2]);
        job->particles[i].wxy_t = 0.5*(job->particles[i].L[1] - job->particles[i].L[2]);
        job->particles[i].eyy_t = job->particles[i].L[3];
    }

    free(pvec);
    free(nvec);
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void map_to_grid_explicit_split(job_t *job, size_t thread_id)
{
    if (thread_id != 0) {
        return;
    }

    // dumb, but needed until we switch to SoA.
    double *pvec = calloc(job->num_particles, sizeof(double));
    double *nvec = calloc(job->num_nodes, sizeof(double));

    for (size_t i = 0; i < job->num_particles; i++) {
        pvec[i] = job->particles[i].m;
    }
    spmd_gatxpy(job->phi, pvec, nvec);
    for (size_t i = 0; i < job->num_nodes; i++) {
        job->nodes[i].m = nvec[i];
        nvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        pvec[i] = job->particles[i].m * job->particles[i].x_t;
    }
    spmd_gatxpy(job->phi, pvec, nvec);
    for (size_t i = 0; i < job->num_nodes; i++) {
        job->nodes[i].mx_t = nvec[i];
        nvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        pvec[i] = job->particles[i].m * job->particles[i].y_t;
    }
    spmd_gatxpy(job->phi, pvec, nvec);
    for (size_t i = 0; i < job->num_nodes; i++) {
        job->nodes[i].my_t = nvec[i];
        nvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        pvec[i] = job->particles[i].m * job->particles[i].bx;
    }
    spmd_gatxpy(job->phi, pvec, nvec);
    for (size_t i = 0; i < job->num_nodes; i++) {
        job->nodes[i].fx = nvec[i];
        nvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        pvec[i] = job->particles[i].m * job->particles[i].by;
    }
    spmd_gatxpy(job->phi, pvec, nvec);
    for (size_t i = 0; i < job->num_nodes; i++) {
        job->nodes[i].fy = nvec[i];
        nvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        pvec[i] = -job->particles[i].v * job->particles[i].sxx;
    }
    spmd_gatxpy(job->dphi_x, pvec, nvec);
    for (size_t i = 0; i < job->num_nodes; i++) {
        job->nodes[i].fx += nvec[i];
        nvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        pvec[i] = -job->particles[i].v * job->particles[i].sxy;
    }
    spmd_gatxpy(job->dphi_x, pvec, nvec);
    for (size_t i = 0; i < job->num_nodes; i++) {
        job->nodes[i].fy += nvec[i];
        nvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        pvec[i] = -job->particles[i].v * job->particles[i].sxy;
    }
    spmd_gatxpy(job->dphi_y, pvec, nvec);
    for (size_t i = 0; i < job->num_nodes; i++) {
        job->nodes[i].fx += nvec[i];
        nvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        pvec[i] = -job->particles[i].v * job->particles[i].syy;
    }
    spmd_gatxpy(job->dphi_y, pvec, nvec);
    for (size_t i = 0; i < job->num_nodes; i++) {
        job->nodes[i].fy += nvec[i];
        nvec[i] = 0;
    }

    free(pvec);
    free(nvec);
    return;

}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void move_particles_explicit_usl_split(job_t *job, size_t p_start, size_t p_stop)
{
    if (p_start != 0) {
        return;
    }

    double *pvec = calloc(job->num_particles, sizeof(double));
    double *nvec = calloc(job->num_nodes, sizeof(double));

    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->nodes[i].m != 0) {
            nvec[i] = job->dt * job->nodes[i].mx_t / job->nodes[i].m;
        } else {
            nvec[i] = 0;
        }
    }
    spmd_gaxpy(job->phi, nvec, pvec);
    for (size_t i = 0; i < job->num_particles; i++) {
        job->particles[i].x += pvec[i];
        job->particles[i].ux += pvec[i];
        pvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->nodes[i].m != 0) {
            nvec[i] = job->dt * job->nodes[i].my_t / job->nodes[i].m;
        } else {
            nvec[i] = 0;
        }
    }
    spmd_gaxpy(job->phi, nvec, pvec);
    for (size_t i = 0; i < job->num_particles; i++) {
        job->particles[i].y += pvec[i];
        job->particles[i].uy += pvec[i];
        pvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->nodes[i].m != 0) {
            nvec[i] = job->dt * job->nodes[i].fx / job->nodes[i].m;
        } else {
            nvec[i] = 0;
        }
    }
    spmd_gaxpy(job->phi, nvec, pvec);
    for (size_t i = 0; i < job->num_particles; i++) {
        job->particles[i].x_t += pvec[i];
        pvec[i] = 0;
    }

    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->nodes[i].m != 0) {
            nvec[i] = job->dt * job->nodes[i].fy / job->nodes[i].m;
        } else {
            nvec[i] = 0;
        }
    }
    spmd_gaxpy(job->phi, nvec, pvec);
    for (size_t i = 0; i < job->num_particles; i++) {
        job->particles[i].y_t += pvec[i];
        pvec[i] = 0;
    }

    free(pvec);
    free(nvec);
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_particle_densities_split(job_t *job, size_t p_start, size_t p_stop)
{
    for (size_t i = p_start; i < p_stop; i++) {
        CHECK_ACTIVE(job, i);
        const double trD = (job->particles[i].exx_t + job->particles[i].eyy_t);
        job->particles[i].v *= exp(job->dt * trD);
    }

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
    free(job->active);

    free(job->u_dirichlet);
    free(job->u_dirichlet_mask);
    free(job->node_number_override);
    free(job->color_indices);

    return;
}
/*----------------------------------------------------------------------------*/

