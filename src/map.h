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
    size_t node_field_offset, size_t particle_field_offset);
void accumulate_p_to_n_doublescalar(node_t *nodes,
    size_t node_field_offset, int *nodelist, double *sfvalues, 
    int nodelist_len, double pdata);
void accumulate_p_to_n_ds_list(node_t *nodes,
    size_t *node_field_offset, int *nodelist, double *sfvalues, 
    size_t nodelist_len, double *pdata, size_t pdata_len);

void accumulate_p_to_n_ds_list47(node_t *nodes,
    size_t *node_field_offset, int *nodelist, double *sfvalues, 
    double *pdata);
void accumulate_p_to_n_ds_list42(node_t *nodes,
    size_t *node_field_offset, int *nodelist, double *sfvalues, 
    double *pdata);
#endif

