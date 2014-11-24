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

enum solver_e {
    IMPLICIT_SOLVER=0,
    EXPLICIT_SOLVER_USF,
    EXPLICIT_SOLVER_USL,
    NUM_SOLVERS
};

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
    char *log_filename;

    char *particle_filename_fullpath;
    char *element_filename_fullpath;
    char *state_filename_fullpath;
    char *log_filename_fullpath;

    FILE *particle_fd;
    FILE *element_fd;
    FILE *state_fd;
    FILE *log_fd;

    FILE *info_fd;

    char *job_name;
    char *job_description;
    int job_id;

    char *date;

    int override_directory_with_id;
    int prepend_date;

    int modified_directory;

    double sample_rate_hz;
} output_control_t;

// forward declare the job structure
struct job_s;

typedef struct mat_control_s {
    const char *material_filename;
    int use_builtin;
    void (*material_init)(struct job_s *);
    void (*calculate_stress)(struct job_s *);
    void (*calculate_stress_threaded)(void *);

    double *fp64_props;
    int *int_props;
    int num_fp64_props;
    int num_int_props;
} material_control_t;

typedef struct bc_control_s {
    const char *bc_filename;
    int use_builtin;
    void (*bc_init)(void *);
    int (*bc_validate)(void *);
    void (*bc_time_varying)(void *);
    void (*bc_momentum)(void *);
    void (*bc_force)(void *);
    
    double *fp64_props;
    int *int_props;
    int num_fp64_props;
    int num_int_props;
} boundary_control_t;

typedef struct job_s {
    double t;
    double dt;
    double t_stop;
    size_t frame;
    int stepcount;

    size_t num_particles;
    size_t num_nodes;
    size_t num_elements;
    size_t num_colors;
    size_t *color_indices;
//    size_t *color_list_lengths;
//    size_t **element_id_by_color;

    /* offsets are arranged [t0c0, t0c1 ... t1c0, t1c1]. */
//    size_t *particle_by_element_color_offsets;
//    size_t *particle_by_element_color_lengths;
//    size_t *particle_by_element_color_list;
    size_t *particle_by_element_color_lengths;
    size_t **particle_by_element_color_lists;

    int N;
    double h;

    int use_cpdi;
    struct timespec tic, toc;

    node_t *nodes;
    particle_t *particles;
    element_t *elements;

    int *in_element;
    int *active;

    /*
        Shape functions (for each particle). Get the node(s) to map it to with
        in_element.
    */
    double *h1;
    double *h2;
    double *h3;
    double *h4;

    double *b11;
    double *b12;
    double *b13;
    double *b14;

    double *b21;
    double *b22;
    double *b23;
    double *b24;

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

    /* solver type */
    enum solver_e solver;

    /* simulation timestep options */
    timestep_control_t timestep;

    /* implicit solver options */
    implicit_control_t implicit;

    /* output control options */
    output_control_t output;

    /* material model options */
    material_control_t material;

    /* boundary condition options */
    boundary_control_t boundary;

    int step_number;
    double step_start_time;

    pthread_barrier_t *step_barrier;
    pthread_barrier_t *serialize_barrier;
    size_t num_threads;

    int *update_elementlists;
    int *update_elementlists_flag;
} job_t;

typedef struct s_threadtask {
    size_t id;
    size_t num_threads;

    /* global array bounds (may not be used). */
    size_t gidx_min;
    size_t gidx_max;

    /*
       This thread is responsible for particles from [offset] to
       [offset+blocksize].
    */
    size_t offset;
    size_t blocksize;


    /*
       This thread is responsible for nodes from [n_offset] to
       [n_offset+n_blocksize].
    */
    size_t n_offset;
    size_t n_blocksize;


    /*
       This thread is responsible for elements from [e_offset] to
       [offset+blocksize].
    */
    size_t e_offset;
    size_t e_blocksize;

    size_t stride;

    job_t *job;
} threadtask_t;

job_t *mpm_init(int N, double h, particle_t *particles, int num_particles, double t);
void implicit_mpm_step(job_t *job);
void explicit_mpm_step_usf(job_t *job);
void explicit_mpm_step_usl(job_t *job);
void explicit_mpm_step_usl_threaded(void *task);
void explicit_mpm_step_musl_threaded(void *task);
void explicit_mpm_step_usf_threaded(void *task);
void mpm_cleanup(job_t *job);

void pt_update_stress(void *args);

void create_particle_to_element_map(job_t *job);
void calculate_shapefunctions(job_t *job);
void calculate_node_velocity(job_t *job);
void calculate_strainrate(job_t *job);
void calculate_strainrate_split(job_t *job, size_t p_start, size_t p_stop);
void update_grid(job_t *job);
void move_grid(job_t *job);
void move_particles(job_t *job);
void update_deformation_gradient(job_t *job);
void update_particle_domains(job_t *job);
void update_particle_densities(job_t *job);

void update_particle_vectors(job_t *job);
void update_corner_positions(job_t *job);

void generate_mappings(job_t *job);
void map_particles_to_nodes(double *node_var, cs *phi, double *particle_var);
void map_gradient_particles_to_nodes(double *node_var, cs *grad_phi, double *particle_var);
void map_nodes_to_particles(double *particle_var, cs *phi_transpose, double *node_var);

void move_grid_split(job_t *job, size_t n_start, size_t n_stop);
void create_particle_to_element_map_split(job_t *job, size_t p_start, size_t p_stop);
void create_particle_to_element_map_threaded(threadtask_t *task);
void find_filled_elements(job_t *job);
void calculate_shapefunctions_split(job_t *job, size_t p_start, size_t p_stop);
void update_stress(job_t *job);
void implicit_solve(job_t *job);
void map_to_grid_explicit(job_t *job);
void map_to_grid_explicit_split(job_t *job, size_t thread_id);
void map_to_grid(job_t *job);
void move_particles_explicit_usf(job_t *job);
void move_particles_explicit_usl(job_t *job);
void move_particles_explicit_usl_split(job_t *job, size_t p_start, size_t p_stop);
void update_particle_densities_split(job_t *job, size_t p_start, size_t p_stop);

#endif

