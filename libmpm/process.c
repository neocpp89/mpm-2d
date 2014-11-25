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

#include <assert.h>

#define TOL 1e-10

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
    ((int)(((xp)<1.0 && (xp)>=0.0 && (yp)<1.0 && (yp)>=0.0)?((floor((xp)/(h)) + floor((yp)/(h))*((N)-1))):(-1)))

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
    job->particles = (particle_t *)calloc(num_particles, sizeof(particle_t));
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
    fprintf(stderr, "Done setting initial particle data.\n");

    /* Get node coordinates. */
    fprintf(stderr, "Each node is %zu bytes.\n", sizeof(node_t));
    job->nodes = (node_t *)calloc(job->num_nodes, sizeof(node_t));
    fprintf(stderr, "%zu bytes (%.2g MB) allocated for %zu nodes.\n",
        job->num_nodes * sizeof(node_t),
        job->num_nodes * sizeof(node_t) / 1048576.0,
        job->num_nodes);
    for (i = 0; i < job->num_nodes; i++) {
        node_number_to_coords(&(job->nodes[i].x), &(job->nodes[i].y), i, N, h);
        /* printf("preprocessing: xn=%g yn=%g\n", job->nodes[i].x, job->nodes[i].y); */
    }

    /* Get node numbering for elements. */
    fprintf(stderr, "Each element is %zu bytes.\n", sizeof(element_t));
    job->elements = (element_t *)calloc(job->num_elements, sizeof(element_t));
    fprintf(stderr, "%zu bytes (%.2g MB) allocated for %zu elements.\n",
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

    /* for dirichlet BCs */
    job->u_dirichlet = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->u_dirichlet_mask = (int *)malloc(NODAL_DOF * sizeof(int) * job->num_nodes);

    /* for periodic BCs */
    job->node_number_override = (int *)malloc(job->vec_len * sizeof(int));

    for (i = 0; i < job->num_particles; i++) {
        job->in_element[i] = WHICH_ELEMENT(
            job->particles[i].x, job->particles[i].y, job->N, job->h);
        if (job->in_element[i] < 0 || job->in_element[i] > job->num_elements) {
            job->active[i] = 0;
        } else {
            job->active[i] = 1;
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

    /* Don't use cpdi here. */
    job->use_cpdi = 0;

    /* swap this out later with a pointer from dlopen if needed. */
    job->material.num_fp64_props = 0;
    job->material.num_int_props = 0;
    job->material.fp64_props = NULL;
    job->material.int_props = NULL;
    job->material.material_filename = NULL;
    //job->material.material_filename = (char *)malloc(16);
    //strcpy(job->material.material_filename, "builtin");
    job->material.material_init = &material_init_linear_elastic;
    job->material.calculate_stress = &calculate_stress_linear_elastic;

    /* used to vary loads/bcs */
    job->step_number = 0;
    job->step_start_time = job->t;

    return job;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void explicit_mpm_step_musl_threaded(void *_task)
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
        (*(job->boundary.bc_time_varying))(job);
    }

    pthread_barrier_wait(job->serialize_barrier);
    /* Map particle state to grid quantites. */
    map_to_grid_explicit_split(job, task->id);

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
   
    /* recalculate nodal velocity by projecting from particles again. */ 
  
    /* 
    for (i = n_start; i < n_stop; i++) {
        if (job->nodes[i].m > TOL)
            printf("before %zu: %g, %g\n", i, job->nodes[i].x_t, job->nodes[i].y_t);
    }
    */

    pthread_barrier_wait(job->serialize_barrier);
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
    pthread_barrier_wait(job->serialize_barrier);
    map_to_grid_explicit_split(job, task->id);
    pthread_barrier_wait(job->serialize_barrier);
    for (i = n_start; i < n_stop; i++) {
        if(job->nodes[i].m > TOL) {
            job->nodes[i].x_t = job->nodes[i].mx_t / job->nodes[i].m;
            job->nodes[i].y_t = job->nodes[i].my_t / job->nodes[i].m;
        } else {
            job->nodes[i].x_t = 0;
            job->nodes[i].y_t = 0;
        }
    }
    rc = pthread_barrier_wait(job->serialize_barrier);
    if (rc == PTHREAD_BARRIER_SERIAL_THREAD) {
        /* Zero perpendicular momentum at edge nodes. */
        (*(job->boundary.bc_momentum))(job);
    }
    pthread_barrier_wait(job->serialize_barrier);

    /*
    for (i = n_start; i < n_stop; i++) {
        if (job->nodes[i].m > TOL)
            printf("after %zu: %g, %g\n", i, job->nodes[i].x_t, job->nodes[i].y_t);
    }
    */

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
        (*(job->boundary.bc_time_varying))(job);
    }

    pthread_barrier_wait(job->serialize_barrier);
    /* Map particle state to grid quantites. */
    map_to_grid_explicit_split(job, task->id);

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
            job->particles[i].x, job->particles[i].y, job->N, job->h);

        if (p != job->in_element[i]) {
            changed = 1;
            fprintf(job->output.log_fd, 
                "[%g] Particle %zu @(%g, %g) left element %d, now in element %d.\n",
                job->t, i, job->particles[i].x, job->particles[i].y,
                job->in_element[i], p);
        }

        /* Update particle element. */
        job->in_element[i] = p;

        if (p == -1) {
            fprintf(job->output.log_fd,
                "[%g] Particle %zu outside of grid (%g, %g), marking as inactive.\n",
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

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
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
            fprintf(stderr, "Particle %zu outside of element %zu (%g, %g).\n", i,
                p, xl, yl);
        }
        job->particles[i].xl = xl;
        job->particles[i].yl = yl;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
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

            accumulate_p_to_n_ds_list47(job->nodes,
                node_field_offsets, job->elements[p].nodes, s, 
                pdata);

            ds[0] = job->b11[p_idx];
            ds[1] = job->b12[p_idx];
            ds[2] = job->b13[p_idx];
            ds[3] = job->b14[p_idx];

            accumulate_p_to_n_ds_list42(job->nodes,
                node_d_offsets, job->elements[p].nodes, ds, 
                &(stressdata[0]));

            ds[0] = job->b21[p_idx];
            ds[1] = job->b22[p_idx];
            ds[2] = job->b23[p_idx];
            ds[3] = job->b24[p_idx];

            accumulate_p_to_n_ds_list42(job->nodes,
                node_d_offsets, job->elements[p].nodes, ds, 
                &(stressdata[1]));
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
void move_particles_explicit_usl_split(job_t *job, size_t p_start, size_t p_stop)
{
    for (size_t i = p_start; i < p_stop; i++) {
        CHECK_ACTIVE(job, i);

        double s[4];
        s[0] = job->h1[i];
        s[1] = job->h2[i];
        s[2] = job->h3[i];
        s[3] = job->h4[i];

        size_t el = job->in_element[i];
        job->particles[i].x_tt = 0;
        job->particles[i].y_tt = 0;
        for (size_t j = 0; j < 4; j++) {
            size_t n = job->elements[el].nodes[j];
            double m = job->nodes[n].m;
            if (m < TOL) { continue; }
            job->particles[i].x_tt += ((s[j] * job->nodes[n].fx) / m);
            job->particles[i].y_tt += ((s[j] * job->nodes[n].fy) / m);
        }

        job->particles[i].x_t += job->dt * job->particles[i].x_tt;
        job->particles[i].y_t += job->dt * job->particles[i].y_tt;

        double dux = 0;
        double duy = 0;
        for (size_t j = 0; j < 4; j++) {
            size_t n = job->elements[el].nodes[j];
            double m = job->nodes[n].m;
            if (m < TOL) { continue; }
            dux += ((s[j] * job->nodes[n].mx_t) / m);
            duy += ((s[j] * job->nodes[n].my_t) / m);
        }
        dux *= job->dt;
        duy *= job->dt;
        job->particles[i].x += dux;
        job->particles[i].y += duy;
        job->particles[i].ux += dux;
        job->particles[i].uy += duy;

    }
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_particle_densities_split(job_t *job, size_t p_start, size_t p_stop)
{
    for (size_t i = p_start; i < p_stop; i++) {
        CHECK_ACTIVE(job, i);
        job->particles[i].v = job->particles[i].v *
            exp(job->dt * (job->particles[i].exx_t + job->particles[i].eyy_t));
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

    free(job->u_dirichlet);
    free(job->u_dirichlet_mask);
    free(job->node_number_override);
    free(job->color_indices);

    return;
}
/*----------------------------------------------------------------------------*/

