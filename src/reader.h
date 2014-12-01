/**
    \file reader.h
    \author Sachith Dunatunga
    \date 04.06.12

    Writes particle data to disk.
*/
#include <stdio.h>
#include "process.h"

#ifndef __READER_H__
#define __READER_H__
#include "particle.h"

/*
    Read the grid in as NxN nodes spanning the domain from 0 to len in both x
    and y directions.
*/
typedef struct grid_s {
    int N;
    double len;
} grid_t;

int read_grid_params(grid_t *grid, const char *fname);
int read_particles(particle_t **particles, size_t *num_particles, const char *fname);
job_t *read_state(FILE *fd);

#endif

