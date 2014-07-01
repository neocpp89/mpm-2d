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

void write_frame(FILE *fd, int frame, double time, job_t *job);
void write_element_frame(FILE *fd, int frame, double time, job_t *job);
void write_state(FILE *fd, job_t *job);

#endif

