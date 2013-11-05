/**
    \file writer.h
    \author Sachith Dunatunga
    \date 04.06.12

    Writes particle data to disk.
*/
#include <stdio.h>
#include <hdf5.h>

#ifndef __WRITER_H__
#define __WRITER_H__
#include "particle.h"
#include "process.h"

typedef struct s_h5writer {
    hid_t particle_type;
    hid_t element_type;
    hid_t node_type;
    hid_t file;
    hid_t pdataspace;
    hid_t pdataset;
    hid_t edataspace;
    hid_t edataset;
} h5writer_t;

void write_frame(FILE *fd, int frame, double time, job_t *job);
void write_element_frame(FILE *fd, int frame, double time, job_t *job);
void write_state(FILE *fd, job_t *job);

hid_t h5_particle_type(void);
hid_t h5_element_type(void);
h5writer_t *h5_init(const char *fname, const job_t *job);
void h5_write_frame(h5writer_t *h5, int frame, job_t *job);
void h5_write_state(const char *hfilename, job_t *job);
void h5_cleanup(h5writer_t *h5state);

#endif

