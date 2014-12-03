/**
    \file map.h
    \author Sachith Dunatunga
    \date 05.11.2013

    Functions to map scalars, gradients of scalars, and diverence of tensors
    to and from material points to nodes.
*/
#ifndef __MAP_H__
#define __MAP_H__
#include <stddef.h>

/* double scalar arguments */
void map_particles_to_nodes_doublescalar(job_t *job,
    const size_t node_field_offset, const size_t particle_field_offset);
void accumulate_p_to_n_doublescalar(node_t *nodes,
    const size_t node_field_offset, const int * restrict nodelist, const double * restrict sfvalues, 
    const size_t nodelist_len, const double pdata);
void accumulate_p_to_n_ds_list(node_t *nodes,
    const size_t * restrict node_field_offset, const int * restrict nodelist, const double * restrict sfvalues, 
    const size_t nodelist_len, const double * restrict pdata, const size_t pdata_len);
void accumulate_p_to_n_ds_list47(node_t *nodes,
    const size_t * restrict node_field_offset, const int * restrict nodelist, const double * restrict sfvalues,
    const double * restrict pdata);
void accumulate_p_to_n_ds_list42(node_t *nodes,
    const size_t * restrict node_field_offset, const int * restrict nodelist, const double * restrict sfvalues, 
    const double * restrict pdata);
#endif

