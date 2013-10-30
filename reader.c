/**
    \file reader.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#include "reader.h"
#include "writer.h"
#include "process.h"

/*---read_grid_params---------------------------------------------------------*/
void read_grid_params(grid_t *grid, char *fname)
{
    int r;
    FILE *fp;

    fp = fopen(fname, "r");
    r = fscanf(fp, "%d", &(grid->N));
    r += fscanf(fp, "%lf", &(grid->len));
    if (r != 2) {
        printf("error reading grid!\n");
        exit(-1);
    }

    fclose(fp);

    return;
}
/*----------------------------------------------------------------------------*/

/*---h5_read_particles--------------------------------------------------------*/
void h5_read_particles(particle_t **particles, int *num_particles, char *fname, int frame)
{
    hid_t particle_type;
    hid_t file;
    hid_t dataset;
    hid_t dataspace;
    herr_t status = 0;

    /* hyperslab size */
    hsize_t offset[2] = {0, frame};
    hsize_t count[2] = {0, 1};
    hsize_t scdims[2] = {0, 1};

    hid_t memspace;

    hsize_t dims[2] = {0, 0};

    int ndims;

    *particles = NULL;
    *num_particles = 0;

    particle_type = h5_particle_type();

    file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen2(file, "/particles", H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    ndims = H5Sget_simple_extent_dims(dataspace, dims, NULL);

    if (frame != 0 && frame >= dims[1]) {
        printf("ndims = %d\n", ndims);
        printf("can't get frame %d, size is (%d %d)\n", frame, (int)dims[0], (int)dims[1]);
        return;
    }

    *num_particles = dims[0];
    *particles = (particle_t *)malloc(*num_particles * sizeof (particle_t));
    count[0] = *num_particles;
    scdims[0] = *num_particles;

    status |= H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, 
        count, NULL);
    memspace = H5Screate_simple(sizeof(scdims)/sizeof(scdims[0]), scdims, NULL);

    status |= H5Dread(dataset, particle_type, memspace, dataspace, H5P_DEFAULT,
                *particles);

    status |= H5Dclose(dataset);
    status |= H5Sclose(dataspace);
    status |= H5Sclose(memspace);
    status |= H5Tclose(particle_type);
    status |= H5Fclose(file);

    return;
}
/*----------------------------------------------------------------------------*/

/*---read_particles-----------------------------------------------------------*/
void read_particles(particle_t **particles, int *num_particles, char *fname)
{
    FILE *fp;
    int i, r;

    char s[16384];
    char *tok;
    double value;
    int field;

    fp = fopen(fname, "r");

    if (NULL == fgets(s, sizeof(s)/sizeof(char), fp)) {
        fprintf(stderr, "can't read initial particle file (num particles).\n");
    }
    r = sscanf(s, "%d", num_particles);
    if (r != 1) {
        printf("Couldn't read number of particles!\n");
        exit(-1);
    }
    *particles = (particle_t *)malloc(sizeof(particle_t) * (*num_particles));


    for (i = 0; i < *num_particles; i++) {
        if (NULL == fgets(s, sizeof(s)/sizeof(char), fp)) {
            fprintf(stderr, "can't read initial particle file.\n");
            continue;
        }
        field = 0;
        tok = strtok(s, " ,");
        while (tok != NULL) {
            sscanf(tok, "%lf", &value);
            switch (field) {
                case 0:
                    (*particles)[i].m = value;
                    break;
                case 1:
                    (*particles)[i].v = value;
                    break;
                case 2:
                    (*particles)[i].x = value;
                    break;
                case 3:
                    (*particles)[i].y = value;
                    break;
                case 4:
                    (*particles)[i].x_t = value;
                    break;
                case 5:
                    (*particles)[i].y_t = value;
                    break;
                case 6:
                    (*particles)[i].sxx = value;
                    break;
                case 7:
                    (*particles)[i].sxy = value;
                    break;
                case 8:
                    (*particles)[i].syy = value;
                    break;
            }
            tok = strtok(NULL, " ,");
            field++;
        }
        if (field < 9) {
/*            fprintf(stderr, "No initial stress state for %d\n", i);*/
            (*particles)[i].sxx = 0;
            (*particles)[i].sxy = 0;
            (*particles)[i].syy = 0;
        }
        if (field < 6) {
            fprintf(stderr, "error reading particle %d\n", i);
        }
    }

    fclose(fp);

    return;
}
/*----------------------------------------------------------------------------*/

/*---read_state---------------------------------------------------------------*/
job_t *read_state(FILE *fd)
{
    int i, j, r;
    int dep;
    job_t *job;

    job = (job_t *)malloc(sizeof(job_t));

    r = fscanf(fd, "%lf %lf %lf", &(job->t), &(job->dt), &(job->t_stop));
    r += fscanf(fd, "%d %d %d", &(job->num_particles), &(job->num_nodes), &(job->num_elements));
    r += fscanf(fd, "%d %lf", &(job->N), &(job->h));

    if (r != 8) {
        printf("Error reading state header!\n");
        exit(-1);
    }

    job->particles = (particle_t *)malloc(sizeof(particle_t) * job->num_particles);
    job->nodes = (node_t *)malloc(sizeof(node_t) * job->num_nodes);
    job->elements = (element_t *)malloc(sizeof(element_t) * job->num_elements);

    /* Allocate space for tracking element->particle map. */
    job->in_element =  (int *)malloc(job->num_particles * sizeof(int));

    /* Allocate space for interpolation functions. */
    job->h1 = (double *)malloc(job->num_particles * sizeof(double));
    job->h2 = (double *)malloc(job->num_particles * sizeof(double));
    job->h3 = (double *)malloc(job->num_particles * sizeof(double));
    job->h4 = (double *)malloc(job->num_particles * sizeof(double));
/*    job->h5 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->h6 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->h7 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->h8 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->h9 = (double *)malloc(job->num_particles * sizeof(double));*/

    job->b11 = (double *)malloc(job->num_particles * sizeof(double));
    job->b12 = (double *)malloc(job->num_particles * sizeof(double));
    job->b13 = (double *)malloc(job->num_particles * sizeof(double));
    job->b14 = (double *)malloc(job->num_particles * sizeof(double));
/*    job->b15 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b16 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b17 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b18 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b19 = (double *)malloc(job->num_particles * sizeof(double));*/

    job->b21 = (double *)malloc(job->num_particles * sizeof(double));
    job->b22 = (double *)malloc(job->num_particles * sizeof(double));
    job->b23 = (double *)malloc(job->num_particles * sizeof(double));
    job->b24 = (double *)malloc(job->num_particles * sizeof(double));
/*    job->b25 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b26 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b27 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b28 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b29 = (double *)malloc(job->num_particles * sizeof(double));*/

    for (i = 0; i < job->num_particles; i++) {
        /* Position */
        r = fscanf(fd, "%lf %lf", &(job->particles[i].x), &(job->particles[i].y));

        /* Local Coordinates */
        r += fscanf(fd, "%lf %lf", &(job->particles[i].xl), &(job->particles[i].yl));

        /* Velocity */
        r += fscanf(fd, "%lf %lf", &(job->particles[i].x_t), &(job->particles[i].y_t));

        /* Mass */
        r += fscanf(fd, "%lf", &(job->particles[i].m));

        /* Volume */
        r += fscanf(fd, "%lf", &(job->particles[i].v));

        /* Initial volume */
        r += fscanf(fd, "%lf", &(job->particles[i].v0));

        /* Stress */
        r += fscanf(fd, "%lf %lf %lf",
            &(job->particles[i].sxx),
            &(job->particles[i].sxy),
            &(job->particles[i].syy));

        /* Strain rate */
        r += fscanf(fd, "%lf %lf %lf %lf",
            &(job->particles[i].exx_t),
            &(job->particles[i].exy_t),
            &(job->particles[i].eyy_t),
            &(job->particles[i].wxy_t));

        /* Body forces */
        r += fscanf(fd, "%lf %lf", &(job->particles[i].bx), &(job->particles[i].by));

        /* Deformation gradient tensor */
        r += fscanf(fd, "%lf %lf %lf %lf",
            &(job->particles[i].Fxx), &(job->particles[i].Fxy),
            &(job->particles[i].Fyx), &(job->particles[i].Fyy));

        /* Displacements */
        r += fscanf(fd, "%lf %lf", &(job->particles[i].ux), &(job->particles[i].uy));

        /* Color used by splot visualization */
        r += fscanf(fd, "%lf", &(job->particles[i].color));

        /* State Variables (for constitutive law) */
        r += fscanf(fd, "%d", &dep);
        for (j = 0; j < dep; j++) {
            r += fscanf(fd, "%lf", &(job->particles[i].state[j]));
        }

        /* Flag particle as active or not (used for discharge problems). */
        r += fscanf(fd, "%d", &(job->particles[i].active));

        if (r != (26 + dep + 1)) {
            printf("error reading state of particle %d\n", i);
        }
    }

    for (i = 0; i < job->num_nodes; i++) {
        r = fscanf(fd, "%lf %lf", &(job->nodes[i].x), &(job->nodes[i].y));
        if (r != 2) {
            printf("error reading state of node %d\n", i);
        }
    }

    for (i = 0; i < job->num_elements; i++) {
        r = fscanf(fd, "%d %d %d %d",
            &(job->elements[i].nodes[0]),
            &(job->elements[i].nodes[1]),
            &(job->elements[i].nodes[2]),
            &(job->elements[i].nodes[3]));
        if (r != 4) {
            printf("error reading state of element %d\n", i);
        }
    }

    printf("Done loading state.\n");

    return job;
}
/*----------------------------------------------------------------------------*/

