/**
    \file process.h
    \author Sachith Dunatunga
    \date 04.06.12

    Contains the structure for a single job in MPM.
*/
#ifndef __PROCESS_H__
#define __PROCESS_H__
#include "particle.h"
#include "node.h"
#include "element.h"
#include <stdio.h>
#include <pthread.h>

#include <suitesparse/cs.h>

typedef struct material_s {
    double E;
    double nu;
    double rho;
} material_t;

typedef struct ts_control_s {
    double dt_max;
    double dt_min;
    double dt;

    int automatic_dt;

    int allow_dt_increase;
    int stable_dt_threshold;
} timestep_control_t;

typedef struct im_control_s {
    double du_norm_ratio;
    double q_norm_ratio;
    double du_norm_converged;

    int unstable_iteration_count;
    int stable_step_count;
} implicit_control_t;

typedef struct op_control_s {
    char *directory;
    char *user;
    char *particle_filename;
    char *element_filename;
    char *state_filename;

    char *particle_filename_fullpath;
    char *element_filename_fullpath;
    char *state_filename_fullpath;

    FILE *particle_fd;
    FILE *element_fd;
    FILE *state_fd;

    char *job_name;
    char *job_description;
    int job_id;

    char *date;

    int override_directory_with_id;
    int prepend_date;

    int modified_directory;
} output_control_t;

typedef struct job_s {
    double t;
    double dt;
    double t_stop;

    int num_particles;
    int num_nodes;
    int num_elements;
    int N;
    double h;

    int use_cpdi;

    node_t *nodes;
    particle_t *particles;
    element_t *elements;

    int *in_element;

    /*
        Shape functions (for each particle). Get the node(s) to map it to with
        in_element.
    */
    double *h1;
    double *h2;
    double *h3;
    double *h4;
//    double *h5;
//    double *h6;
//    double *h7;
//    double *h8;
//    double *h9;

    double *b11;
    double *b12;
    double *b13;
    double *b14;
//    double *b15;
//    double *b16;
//    double *b17;
//    double *b18;
//    double *b19;

    double *b21;
    double *b22;
    double *b23;
    double *b24;
//    double *b25;
//    double *b26;
//    double *b27;
//    double *b28;
//    double *b29;

    double *u_grid;
    double *du_grid;
    double *v_grid;
    double *a_grid;
    double *m_grid;
    double *q_grid;
    double *f_ext_grid;
    double *f_int_grid;

    /*
        Maps particle to nodal quantites (and gradient of quantites).
        The matrix is a compressed sparse matrix, and the transpose will
        map nodal to particle quantites.
    */
    cs *phi;
    cs *grad_phi;
    cs *phi_transpose;
    cs *grad_phi_transpose;

    double *m_nodes;
    double *mx_t_nodes;
    double *my_t_nodes;
    double *mx_tt_nodes;
    double *my_tt_nodes;
    double *x_t_nodes;
    double *y_t_nodes;
    double *x_tt_nodes;
    double *y_tt_nodes;
    double *fx_nodes;
    double *fy_nodes;

    /* maps node number to position in u_grid vector. */
    int *node_u_map;
    int *inv_node_u_map;

    /* saves which dofs are controlled by BCs and their values */
    double *u_dirichlet;
    int *u_dirichlet_mask;

    /*
        Renumbers the nodes before feeding to the node map. Used to implement
        periodic boundary conditions.
    */
    int *node_number_override;

    int vec_len;
    double *kku_grid;
    /* kku_grid is vec_len*vec_len in size */

//    FILE *ke_data;
//    FILE *stress_data;

    /* simulation timestep options */
    timestep_control_t timestep;

    /* implicit solver options */
    implicit_control_t implicit;

    /* output control options */
    output_control_t output;

#define PT_NUM_THREADS 64
    pthread_t threads[PT_NUM_THREADS];
} job_t;

typedef struct s_threadtask {
    int id;
    int num_threads;

    job_t *job;
} threadtask_t;

job_t *mpm_init(int N, double h, particle_t *particles, int num_particles, double t);
void mpm_step(job_t *job);
void mpm_cleanup(job_t *job);

typedef struct s_strided_task {
    job_t *job;
    int offset;
} stask_t;
void pt_update_stress(void *args);

void create_particle_to_element_map(job_t *job);
void calculate_shapefunctions(job_t *job);
void calculate_node_velocity(job_t *job);
void calculate_strainrate(job_t *job);
void update_grid(job_t *job);
void move_particles(job_t *job);
void update_deformation_gradient(job_t *job);
void update_particle_domains(job_t *job);
void update_particle_densities(job_t *job);
void update_corner_domains(job_t *job);

void generate_mappings(job_t *job);
void map_particles_to_nodes(double *node_var, cs *phi, double *particle_var);
void map_gradient_particles_to_nodes(double *node_var, cs *grad_phi, double *particle_var);
void map_nodes_to_particles(double *particle_var, cs *phi_transpose, double *node_var);

#endif

