/**
    \file reader.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "reader.h"
#include "writer.h"
#include "process.h"

/*---read_grid_params---------------------------------------------------------*/
int read_grid_params(grid_t *grid, const char *fname)
{
    int r;
    FILE *fp;

    fp = fopen(fname, "r");
    if (fp == NULL) {
        return -1;
    }
    r = fscanf(fp, "%d", &(grid->N));
    r += fscanf(fp, "%lg", &(grid->len));
    if (r != 2) {
        printf("error reading grid!\n");
        exit(-1);
    }

    fclose(fp);

    return 0;
}
/*----------------------------------------------------------------------------*/

/*---read_particles-----------------------------------------------------------*/
int read_particles(particle_t **particles, int *num_particles, const char *fname)
{
    FILE *fp;
    int i, r;

    char s[16384];
    char *tok;
    double value;
    int field;

    fp = fopen(fname, "r");
    if (fp == NULL) {
        return -1;
    }

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
            sscanf(tok, "%lg", &value);
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

    return 0;
}
/*----------------------------------------------------------------------------*/

/*---read_state---------------------------------------------------------------*/
job_t *read_state(FILE *fd)
{
    int i, j, r;
    int dep;
    job_t *job;

    job = (job_t *)malloc(sizeof(job_t));

    r = fscanf(fd, "%lg %lg %lg", &(job->t), &(job->dt), &(job->t_stop));
    r += fscanf(fd, "%zu %zu %zu", &(job->num_particles), &(job->num_nodes), &(job->num_elements));
    r += fscanf(fd, "%d %lg", &(job->N), &(job->h));

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
        r = fscanf(fd, "%lg %lg", &(job->particles[i].x), &(job->particles[i].y));

        /* Local Coordinates */
        r += fscanf(fd, "%lg %lg", &(job->particles[i].xl), &(job->particles[i].yl));

        /* Velocity */
        r += fscanf(fd, "%lg %lg", &(job->particles[i].x_t), &(job->particles[i].y_t));

        /* Mass */
        r += fscanf(fd, "%lg", &(job->particles[i].m));

        /* Volume */
        r += fscanf(fd, "%lg", &(job->particles[i].v));

        /* Initial volume */
        r += fscanf(fd, "%lg", &(job->particles[i].v0));

        /* Stress */
        r += fscanf(fd, "%lg %lg %lg",
            &(job->particles[i].sxx),
            &(job->particles[i].sxy),
            &(job->particles[i].syy));

        /* Strain rate */
        r += fscanf(fd, "%lg %lg %lg %lg",
            &(job->particles[i].exx_t),
            &(job->particles[i].exy_t),
            &(job->particles[i].eyy_t),
            &(job->particles[i].wxy_t));

        /* Body forces */
        r += fscanf(fd, "%lg %lg", &(job->particles[i].bx), &(job->particles[i].by));

        /* Deformation gradient tensor */
        r += fscanf(fd, "%lg %lg %lg %lg",
            &(job->particles[i].Fxx), &(job->particles[i].Fxy),
            &(job->particles[i].Fyx), &(job->particles[i].Fyy));

        /* Displacements */
        r += fscanf(fd, "%lg %lg", &(job->particles[i].ux), &(job->particles[i].uy));

        /* Color used by splot visualization */
        r += fscanf(fd, "%lg", &(job->particles[i].color));

        /* State Variables (for constitutive law) */
        r += fscanf(fd, "%d", &dep);
        for (j = 0; j < dep; j++) {
            r += fscanf(fd, "%lg", &(job->particles[i].state[j]));
        }

        /* Flag particle as active or not (used for discharge problems). */
        r += fscanf(fd, "%d", &(job->active[i]));

        if (r != (26 + dep + 1)) {
            printf("error reading state of particle %d\n", i);
        }
    }

    for (i = 0; i < job->num_nodes; i++) {
        r = fscanf(fd, "%lg %lg", &(job->nodes[i].x), &(job->nodes[i].y));
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

