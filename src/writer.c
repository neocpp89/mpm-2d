/**
    \file writer.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "process.h"
#include "writer.h"

/*
    Since we don't have reflection in C, make a structure to hold relevant data
    to output as header information. Also, the macro is used to generate the
    structure to enforce consistency.
*/
typedef struct {
    size_t offset;
    const char *format_specifier;
    const char *fieldname;
} headerinfo_t;
#define HEADER_FP64(field) { offsetof(particle_t, field), "%lg", #field }
#define HEADERINFO(field, format) { offsetof(particle_t, field), format, #field }
#define NAMEDHEADERINFO(field, format, name) { offsetof(particle_t, field), format, name }

headerinfo_t particle_out_fields[] = {
    /* Position */
    HEADERINFO(x, "%lg"),
    HEADERINFO(y, "%lg"),

    /* Local Coordinates */
/*    HEADERINFO(xl, "%lg"),*/
/*    HEADERINFO(yl, "%lg"),*/

    /* Volume */
    HEADERINFO(v, "%lg"),

    /* Initial volume */
/*    HEADERINFO(v0, "%lg"),*/

    /* Skip shapefunctions in output. */

    /* Mass */
    HEADERINFO(m, "%lg"),

    /* Velocity */
    HEADERINFO(x_t, "%lg"),
    HEADERINFO(y_t, "%lg"),

    /* Acceleration */
/*    HEADERINFO(x_tt, "%lg"),*/
/*    HEADERINFO(y_tt, "%lg"),*/

    /* Body forces */
/*    HEADERINFO(bx, "%lg"),*/
/*    HEADERINFO(by, "%lg"),*/

    /* Stress */
    HEADERINFO(sxx, "%lg"),
    HEADERINFO(sxy, "%lg"),
    HEADERINFO(syy, "%lg"),

    /* Skip full 3D stress tensor */

    /* Strain rate */
/*    HEADERINFO(exx_t, "%lg"),*/
/*    HEADERINFO(exy_t, "%lg"),*/
/*    HEADERINFO(eyy_t, "%lg"),*/
/*    HEADERINFO(wxy_t, "%lg"),*/

    /* Skip full 3D velocity gradient tensor and Deformation Gradient */

    /* Displacements */
/*    HEADERINFO(ux, "%lg"),*/
/*    HEADERINFO(uy, "%lg"),*/

    /* Color used by splot visualization */
    HEADERINFO(color, "%lg"),

    /* State Variables (for constitutive law) */
/*    double state[DEPVAR];*/
    NAMEDHEADERINFO(state[9], "%lg", "gammap"), 
    NAMEDHEADERINFO(state[10], "%lg", "gammadotp") /* gammadotp */
};

/*---Version 2 of output format-----------------------------------------------*/
size_t v2_write_frame(const char *directory, FILE *metafd, job_t *job,
    size_t (*particle_output_fn)(FILE *, job_t *, particle_t *),
    size_t (*element_output_fn)(FILE *, job_t *, element_t *))
{
    size_t i = 0;
    size_t bytes_out = 0;
    size_t particles_written = 0;
    char fp_name[1024];
    char fp_name_nobase[1024];
    FILE *fp = NULL;

    snprintf(fp_name_nobase, 1024, "fp_%d.h.csv", job->frame);
    snprintf(fp_name, 1024, "%s%s", directory, fp_name_nobase);
    fp = fopen(fp_name, "w+");
    if (fp != NULL) {
        if (particle_output_fn != NULL) {

            /* Write header lines for csv file. */
            bytes_out += fprintf(fp, "id");
            for (i = 0; i < sizeof(particle_out_fields)/sizeof(particle_out_fields[0]); i++) {
                bytes_out += fprintf(fp, ",%s",
                    particle_out_fields[i].fieldname);
            }
            bytes_out += fprintf(fp, "\n");

            /* Write particle data from simulation. */
            for (i = 0; i < job->num_particles; i++) {
                if (job->active[i]) {
                    bytes_out += (*particle_output_fn)(fp, job, &(job->particles[i]));
                    particles_written++;
                }
            }
        }
        fclose(fp);
    }

    if (element_output_fn != NULL) {
    }

    /*
        Write metadata to file (if it exists).
    */
    if (metafd != NULL) {
        bytes_out += fprintf(metafd, "frame-%d = %s,%zu,%d,%lg\n",
            job->frame, fp_name_nobase, particles_written, job->frame, job->t);
    }

    return bytes_out;
}

size_t v2_write_particle(FILE *fd, job_t *job, particle_t *p)
{
    size_t i = 0;
    size_t bytes_out = 0;

    bytes_out += fprintf(fd, "%zu", p->id); 
    for (i = 0; i < sizeof(particle_out_fields)/sizeof(particle_out_fields[0]); i++) {
        bytes_out += fprintf(fd, ",");
        bytes_out += fprintf(fd,
            particle_out_fields[i].format_specifier,
            *(double *)((char *)p + particle_out_fields[i].offset));
    }
    bytes_out += fprintf(fd, "\n");

    return bytes_out;
}

#undef HEADERINFO
/*----------------------------------------------------------------------------*/

/*---write_particle-----------------------------------------------------------*/
void write_particle(FILE *fd, particle_t p, double active)
{
    int i, j;
    fprintf(fd, "%lg %lg %lg %lg %lg %lg ",
        p.m, p.v, p.x, p.y, p.x_t, p.y_t);
    fprintf(fd, "%lg %lg %lg %lg %lg %lg %lg %lg %lg",
        p.sxx, p.sxy, p.syy, p.ux, p.uy, p.state[9], p.color, p.state[10],
        (double)active);

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 2; j++) {
            fprintf(fd, " %lg", p.corners[i][j]);
        }
    }
    fprintf(fd, "\n");

    return;
}
/*----------------------------------------------------------------------------*/

/*---write_frame--------------------------------------------------------------*/
void write_frame(FILE *fd, int frame, double time, job_t *job)
{
    int i;

    fprintf(fd, "%d %lg %d\n", frame, time, job->num_particles);
    for (i = 0; i < job->num_particles; i++) {
        write_particle(fd, job->particles[i], (double)(job->active[i]));
    }

    /* dump the entire frame to disk (or wherever) */
    fflush(fd);

    return;
}
/*----------------------------------------------------------------------------*/

/*---write_element_frame------------------------------------------------------*/
void write_element_frame(FILE *fd, int frame, double time, job_t *job)
{
    int i;
    int e;
    double *sxx_acc;
    double *sxy_acc;
    double *syy_acc;
    double *v_acc;

    double x;
    double y;

    sxx_acc = (double *)malloc(sizeof(double) * job->num_elements);
    sxy_acc = (double *)malloc(sizeof(double) * job->num_elements);
    syy_acc = (double *)malloc(sizeof(double) * job->num_elements);
    v_acc = (double *)malloc(sizeof(double) * job->num_elements);

    for (i = 0; i < job->num_elements; i++) {
        sxx_acc[i] = 0.0f;
        sxy_acc[i] = 0.0f;
        syy_acc[i] = 0.0f;
        v_acc[i] = 0.0f;
    }

    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }
        e = job->in_element[i];
        sxx_acc[e] += job->particles[i].v * job->particles[i].sxx;
        sxy_acc[e] += job->particles[i].v * job->particles[i].sxy;
        syy_acc[e] += job->particles[i].v * job->particles[i].syy;
        v_acc[e] += job->particles[i].v;
    }

    fprintf(fd, "%d %lg %d\n", frame, time, job->num_elements);

    for (i = 0; i < job->num_elements; i++) {
        node_number_to_coords(&x, &y, job->elements[i].nodes[0], job->N, job->h);
        fprintf(fd, "%lg %lg ", x, y);
        node_number_to_coords(&x, &y, job->elements[i].nodes[1], job->N, job->h);
        fprintf(fd, "%lg %lg ", x, y);
        node_number_to_coords(&x, &y, job->elements[i].nodes[2], job->N, job->h);
        fprintf(fd, "%lg %lg ", x, y);
        node_number_to_coords(&x, &y, job->elements[i].nodes[3], job->N, job->h);
        fprintf(fd, "%lg %lg ", x, y);
        fprintf(fd, "%lg %lg %lg\n",
            (v_acc[i] > 0) ? (sxx_acc[i] / v_acc[i]) : 0.0f,
            (v_acc[i] > 0) ? (sxy_acc[i] / v_acc[i]) : 0.0f,
            (v_acc[i] > 0) ? (syy_acc[i] / v_acc[i]) : 0.0f);
    }

    free(sxx_acc);
    free(sxy_acc);
    free(syy_acc);
    free(v_acc);

    return;
}
/*----------------------------------------------------------------------------*/

/*---write_state--------------------------------------------------------------*/
void write_state(FILE *fd, job_t *job)
{
    int i, j;

    fprintf(fd, "%lg %lg %lg\n", job->t, job->dt, job->t_stop);
    fprintf(fd, "%d %d %d\n", job->num_particles, job->num_nodes, job->num_elements);
    fprintf(fd, "%d %lg\n", job->N, job->h);

    for (i = 0; i < job->num_particles; i++) {
        /* Position */
        fprintf(fd, "%lg %lg\n", job->particles[i].x, job->particles[i].y);

        /* Local Coordinates */
        fprintf(fd, "%lg %lg\n", job->particles[i].xl, job->particles[i].yl);

        /* Velocity */
        fprintf(fd, "%lg %lg\n", job->particles[i].x_t, job->particles[i].y_t);

        /* Mass */
        fprintf(fd, "%lg\n", job->particles[i].m);

        /* Volume */
        fprintf(fd, "%lg\n", job->particles[i].v);

        /* Initial volume */
        fprintf(fd, "%lg\n", job->particles[i].v0);

        /* Stress */
        fprintf(fd, "%lg %lg %lg\n",
            job->particles[i].sxx,
            job->particles[i].sxy,
            job->particles[i].syy);

        /* Stress rate */
/*        fprintf(fd, "%lg %lg %lg\n",*/
/*            job->particles[i].sxx_t,*/
/*            job->particles[i].sxy_t,*/
/*            job->particles[i].syy_t);*/

        /* Strain rate */
        fprintf(fd, "%lg %lg %lg %lg\n",
            job->particles[i].exx_t,
            job->particles[i].exy_t,
            job->particles[i].eyy_t,
            job->particles[i].wxy_t);

        /* Jaumann stress increment */
/*        fprintf(fd, "%lg %lg %lg\n",*/
/*            job->particles[i].dsjxx,*/
/*            job->particles[i].dsjxy,*/
/*            job->particles[i].dsjyy);*/

        /* Elastic and plastic strain increments */
/*        fprintf(fd, "%lg %lg %lg\n",*/
/*            job->particles[i].deexx,*/
/*            job->particles[i].deexy,*/
/*            job->particles[i].deeyy);*/
/*        fprintf(fd, "%lg %lg %lg\n",*/
/*            job->particles[i].depxx,*/
/*            job->particles[i].depxy,*/
/*            job->particles[i].depyy);*/

        /* Body forces */
        fprintf(fd, "%lg %lg\n", job->particles[i].bx, job->particles[i].by);

        /* Deformation gradient tensor */
        fprintf(fd, "%lg %lg %lg %lg\n",
            job->particles[i].Fxx, job->particles[i].Fxy,
            job->particles[i].Fyx, job->particles[i].Fyy);

        /* Initial particle domain vectors */
/*        fprintf(fd, "%lg %lg %lg %lg\n",*/
/*            job->particles[i].r1x0, job->particles[i].r1y0,*/
/*            job->particles[i].r2x0, job->particles[i].r2y0);*/

        /* Updated particle domain vectors */
/*        fprintf(fd, "%lg %lg %lg %lg\n",*/
/*            job->particles[i].r1xn, job->particles[i].r1yn,*/
/*            job->particles[i].r2xn, job->particles[i].r2yn);*/

        /* Corner positions */
/*        fprintf(fd, "%lg %lg %lg %lg %lg %lg %lg %lg\n",*/
/*            job->particles[i].c1x, job->particles[i].c1y,*/
/*            job->particles[i].c2x, job->particles[i].c2y,*/
/*            job->particles[i].c3x, job->particles[i].c3y,*/
/*            job->particles[i].c4x, job->particles[i].c4y);*/

        /* Displacements */
        fprintf(fd, "%lg %lg\n", job->particles[i].ux, job->particles[i].uy);

        /* Color used by splot visualization */
        fprintf(fd, "%lg\n", job->particles[i].color);

        /* State Variables (for constitutive law) */
        fprintf(fd, "%d\n", DEPVAR);
        for (j = 0; j < DEPVAR; j++) {
            fprintf(fd, "%lg\n", job->particles[i].state[j]);
        }

        /* Flag particle as active or not (used for discharge problems). */
        fprintf(fd, "%d\n", job->active[i]);
    }

    for (i = 0; i < job->num_nodes; i++) {
        fprintf(fd, "%lg %lg\n", job->nodes[i].x, job->nodes[i].y);
    }

    for (i = 0; i < job->num_elements; i++) {
         fprintf(fd, "%d %d %d %d\n",
            job->elements[i].nodes[0],
            job->elements[i].nodes[1],
            job->elements[i].nodes[2],
            job->elements[i].nodes[3]);
    }

    return;
}
/*----------------------------------------------------------------------------*/

