/**
    \file writer.h
    \author Sachith Dunatunga
    \date 04.06.12

    Writes particle data to disk.
*/
#include <stdio.h>

#ifndef __WRITER_H__
#define __WRITER_H__
#include "particle.h"
#include "process.h"


size_t v2_write_frame(const char *directory, FILE *metafd, job_t *job,
    size_t (*particle_output_fn)(FILE *, particle_t *),
    size_t (*element_output_fn)(FILE *, element_t *));
size_t v2_write_particle(FILE *fd, particle_t *p);

void write_frame(FILE *fd, size_t frame, double time, job_t *job);
void write_element_frame(FILE *fd, size_t frame, double time, job_t *job);
void write_state(FILE *fd, job_t *job);

#endif

