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

void map_particles_to_nodes_doublescalar(job_t *job,
    size_t node_field_offset, size_t particle_field_offset);

#endif

