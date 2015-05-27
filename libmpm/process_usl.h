/**
    \file process_usl.h
    \author Sachith Dunatunga
    \date 29.10.2014

    Contains the structure for a single job in MPM.
*/
#ifndef __PROCESS_USL_H__
#define __PROCESS_USL_H__
#include "process.h"

job_t *mpm_init(int N, double hx, double hy, double Lx, double Ly, particle_t *particles, size_t num_particles, double t);
void explicit_mpm_step_musl_threaded(void *_task);
void explicit_mpm_step_usl_threaded(void *_task);
void move_grid_split(job_t *job, size_t n_start, size_t n_stop);
void create_particle_to_element_map_threaded(threadtask_t *task);
void find_filled_elements(job_t *job);
void calculate_shapefunctions_split(job_t *job, size_t p_start, size_t p_stop);
void calculate_strainrate_split(job_t *job, size_t p_start, size_t p_stop);
void map_to_grid_explicit_split(job_t *job, size_t thread_id);
void move_particles_explicit_usl_split(job_t *job, size_t p_start, size_t p_stop);
void update_particle_densities_split(job_t *job, size_t p_start, size_t p_stop);
void mpm_cleanup(job_t *job);

#endif //__PROCESS_USL_H__

