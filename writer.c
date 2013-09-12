/**
    \file writer.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include "process.h"
#include "writer.h"

/*----------------------------------------------------------------------------*/
h5writer_t *h5_init(const char *fname, const job_t *job)
{
    h5writer_t *h5state;
    herr_t status = 0;

    hid_t plist;
    hid_t elist;

    hsize_t pdims[2] = {job->num_particles, 1};
    hsize_t edims[2] = {job->num_elements, 1};
    hsize_t max_pdims[2] = {job->num_particles, H5S_UNLIMITED};
    hsize_t max_edims[2] = {job->num_elements, H5S_UNLIMITED};

    /* For compression. */
    hsize_t chunk_pdims[2] = {job->num_particles, 1};
    hsize_t chunk_edims[2] = {job->num_elements, 1};

    h5state = (h5writer_t *)malloc(sizeof(h5writer_t));
    h5state->particle_type = h5_particle_type();
    h5state->element_type = h5_element_type();
    h5state->file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    h5state->pdataspace = H5Screate_simple(sizeof(pdims)/sizeof(pdims[0]), pdims, max_pdims);
    h5state->edataspace = H5Screate_simple(sizeof(edims)/sizeof(edims[0]), edims, max_edims);

    plist  = H5Pcreate(H5P_DATASET_CREATE);
    status |= H5Pset_chunk(plist, sizeof(chunk_pdims)/sizeof(chunk_pdims[0]), chunk_pdims);
    status |= H5Pset_deflate(plist, 1);
    elist  = H5Pcreate(H5P_DATASET_CREATE);
    status |= H5Pset_chunk(elist, sizeof(chunk_edims)/sizeof(chunk_edims[0]), chunk_edims);
    status |= H5Pset_deflate(elist, 1);

    h5state->pdataset = H5Dcreate2(h5state->file, "/particles",
                        h5state->particle_type, h5state->pdataspace,
                        H5P_DEFAULT, plist, H5P_DEFAULT);
    h5state->edataset = H5Dcreate2(h5state->file, "/elements",
                        h5state->element_type, h5state->edataspace,
                        H5P_DEFAULT, elist, H5P_DEFAULT);

    return h5state;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void h5_cleanup(h5writer_t *h5state)
{
    herr_t status = 0;

    status |= H5Dclose(h5state->pdataset);
    status |= H5Sclose(h5state->pdataspace);
    status |= H5Dclose(h5state->edataset);
    status |= H5Sclose(h5state->edataspace);
    status |= H5Fclose(h5state->file);
    free(h5state);

    return;
}
/*----------------------------------------------------------------------------*/

/*---h5_particle_type---------------------------------------------------------*/
hid_t h5_particle_type()
{
    int i;
    char buf[80];
    hid_t hparticle_type;

    hparticle_type = H5Tcreate(H5T_COMPOUND, sizeof(particle_t));
    H5Tinsert(hparticle_type, "x", HOFFSET(particle_t, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "y", HOFFSET(particle_t, y), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "xl", HOFFSET(particle_t, xl), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "yl", HOFFSET(particle_t, yl), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "x_t", HOFFSET(particle_t, x_t), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "y_t", HOFFSET(particle_t, y_t), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "m", HOFFSET(particle_t, m), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "v", HOFFSET(particle_t, v), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "v0", HOFFSET(particle_t, v0), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "sxx", HOFFSET(particle_t, sxx), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "sxy", HOFFSET(particle_t, sxy), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "syy", HOFFSET(particle_t, syy), H5T_NATIVE_DOUBLE);
/*    H5Tinsert(hparticle_type, "sxx_t", HOFFSET(particle_t, sxx_t), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "sxy_t", HOFFSET(particle_t, sxy_t), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "syy_t", HOFFSET(particle_t, syy_t), H5T_NATIVE_DOUBLE);*/
    H5Tinsert(hparticle_type, "exx_t", HOFFSET(particle_t, exx_t), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "exy_t", HOFFSET(particle_t, exy_t), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "eyy_t", HOFFSET(particle_t, eyy_t), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "wxy_t", HOFFSET(particle_t, wxy_t), H5T_NATIVE_DOUBLE);
/*    H5Tinsert(hparticle_type, "dsjxx", HOFFSET(particle_t, dsjxx), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "dsjxy", HOFFSET(particle_t, dsjxy), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "dsjyy", HOFFSET(particle_t, dsjyy), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "deexx", HOFFSET(particle_t, deexx), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "deexy", HOFFSET(particle_t, deexy), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "deeyy", HOFFSET(particle_t, deeyy), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "depxx", HOFFSET(particle_t, depxx), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "depxy", HOFFSET(particle_t, depxy), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "depyy", HOFFSET(particle_t, depyy), H5T_NATIVE_DOUBLE);*/
    H5Tinsert(hparticle_type, "bx", HOFFSET(particle_t, bx), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "by", HOFFSET(particle_t, by), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "Fxx", HOFFSET(particle_t, Fxx), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "Fxy", HOFFSET(particle_t, Fxy), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "Fyx", HOFFSET(particle_t, Fyx), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "Fyy", HOFFSET(particle_t, Fyy), H5T_NATIVE_DOUBLE);
/*    H5Tinsert(hparticle_type, "r1x0", HOFFSET(particle_t, r1x0), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "r1y0", HOFFSET(particle_t, r1y0), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "r2x0", HOFFSET(particle_t, r2x0), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "r2y0", HOFFSET(particle_t, r2y0), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "r1xn", HOFFSET(particle_t, r1xn), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "r1yn", HOFFSET(particle_t, r1yn), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "r2xn", HOFFSET(particle_t, r2xn), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "r2yn", HOFFSET(particle_t, r2yn), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "c1x", HOFFSET(particle_t, c1x), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "c1y", HOFFSET(particle_t, c1y), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "c2x", HOFFSET(particle_t, c2x), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "c2y", HOFFSET(particle_t, c2y), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "c3x", HOFFSET(particle_t, c3x), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "c3y", HOFFSET(particle_t, c3y), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "c4x", HOFFSET(particle_t, c4x), H5T_NATIVE_DOUBLE);*/
/*    H5Tinsert(hparticle_type, "c4y", HOFFSET(particle_t, c4y), H5T_NATIVE_DOUBLE);*/
    H5Tinsert(hparticle_type, "ux", HOFFSET(particle_t, ux), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "uy", HOFFSET(particle_t, uy), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "color", HOFFSET(particle_t, color), H5T_NATIVE_DOUBLE);
    for (i = 0; i < DEPVAR; i++) {
        snprintf(buf, sizeof(buf), "state%d", i);
        H5Tinsert(hparticle_type, buf, HOFFSET(particle_t, state[i]), H5T_NATIVE_DOUBLE);
    }
    H5Tinsert(hparticle_type, "active", HOFFSET(particle_t, active), H5T_NATIVE_INT);
    H5Tinsert(hparticle_type, "material", HOFFSET(particle_t, material), H5T_NATIVE_INT);
    H5Tinsert(hparticle_type, "x_tt", HOFFSET(particle_t, x_tt), H5T_NATIVE_DOUBLE);
    H5Tinsert(hparticle_type, "y_tt", HOFFSET(particle_t, y_tt), H5T_NATIVE_DOUBLE);

    return hparticle_type;
}
/*----------------------------------------------------------------------------*/

/*---h5_element_type----------------------------------------------------------*/
hid_t h5_element_type()
{
    int i;
    char buf[80];
    hid_t helement_type;

    helement_type = H5Tcreate(H5T_COMPOUND, sizeof(element_t));
    for (i = 0; i < 9; i++) {
        snprintf(buf, sizeof(buf), "nodes%d", i);
        H5Tinsert(helement_type, buf, HOFFSET(element_t, nodes[i]), H5T_NATIVE_INT);
    }
    H5Tinsert(helement_type, "n", HOFFSET(element_t, n), H5T_NATIVE_INT);
    H5Tinsert(helement_type, "m", HOFFSET(element_t, m), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "filled", HOFFSET(element_t, filled), H5T_NATIVE_INT);
    for (i = 0; i < 8; i++) {
        snprintf(buf, sizeof(buf), "neighbors%d", i);
        H5Tinsert(helement_type, buf, HOFFSET(element_t, neighbors[i]), H5T_NATIVE_INT);
    }
    H5Tinsert(helement_type, "grad_x", HOFFSET(element_t, grad_x), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "grad_y", HOFFSET(element_t, grad_y), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "grad_mag", HOFFSET(element_t, grad_mag), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "n_x", HOFFSET(element_t, n_x), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "n_y", HOFFSET(element_t, n_y), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "n_theta", HOFFSET(element_t, n_theta), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "v", HOFFSET(element_t, v), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "p_index", HOFFSET(element_t, p_index), H5T_NATIVE_INT);
    H5Tinsert(helement_type, "p_dist", HOFFSET(element_t, p_dist), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "sxx", HOFFSET(element_t, sxx), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "sxy", HOFFSET(element_t, sxy), H5T_NATIVE_DOUBLE);
    H5Tinsert(helement_type, "syy", HOFFSET(element_t, syy), H5T_NATIVE_DOUBLE);

    return helement_type;
}
/*----------------------------------------------------------------------------*/

/*---write_particle-----------------------------------------------------------*/
void inline write_particle(FILE *fd, particle_t p)
{

    fprintf(fd, "%f %f %f %f %f %f ",
        p.m, p.v, p.x, p.y, p.x_t, p.y_t);
    fprintf(fd, "%f %f %f %18.18lf %18.18lf %18.18lf %f %18.18lf %f\n",
        p.sxx, p.sxy, p.syy, p.ux, p.uy, p.state[9], p.color, p.state[10],
        (double)p.active);

    return;
}
/*----------------------------------------------------------------------------*/

/*---write_frame--------------------------------------------------------------*/
void write_frame(FILE *fd, int frame, double time, job_t *job)
{
    int i;

    fprintf(fd, "%d %f %d\n", frame, time, job->num_particles);
    for (i = 0; i < job->num_particles; i++) {
        write_particle(fd, job->particles[i]);
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*---h5_write_frame-----------------------------------------------------------*/
void h5_write_frame(h5writer_t *h5, int frame, job_t *job)
{
    herr_t status = 0;
    hsize_t pdims[2] = {job->num_particles, frame+1};
    hsize_t edims[2] = {job->num_elements, frame+1};

    /* hyperslab size */
    hsize_t poffset[2] = {0, frame};
    hsize_t pcount[2] = {job->num_particles, 1};
    hsize_t scpdims[2] = {job->num_particles, 1};
    hsize_t eoffset[2] = {0, frame};
    hsize_t ecount[2] = {job->num_elements, 1};
    hsize_t scedims[2] = {job->num_elements, 1};

    hid_t pmemspace;
    hid_t ememspace;

    /* increase size of file to add frame... */
    status |= H5Dset_extent(h5->pdataset, pdims);
    status |= H5Sclose(h5->pdataspace);
    h5->pdataspace = H5Dget_space(h5->pdataset);
    status |= H5Dset_extent(h5->edataset, edims);
    status |= H5Sclose(h5->edataspace);
    h5->edataspace = H5Dget_space(h5->edataset);

    status |= H5Sselect_hyperslab(h5->pdataspace, H5S_SELECT_SET, poffset, NULL, 
        pcount, NULL);
    status |= H5Sselect_hyperslab(h5->edataspace, H5S_SELECT_SET, eoffset, NULL, 
        ecount, NULL);

    /* memory space */
    pmemspace = H5Screate_simple(sizeof(scpdims)/sizeof(scpdims[0]), scpdims, NULL);
    ememspace = H5Screate_simple(sizeof(scedims)/sizeof(scedims[0]), scedims, NULL);

    /* first dim is particle id, second dim is time */
    status |= H5Dwrite(h5->pdataset, h5->particle_type, pmemspace,
                h5->pdataspace, H5P_DEFAULT, job->particles);
    status |= H5Sclose(pmemspace);
    status |= H5Dwrite(h5->edataset, h5->element_type, ememspace,
                h5->edataspace, H5P_DEFAULT, job->elements);
    status |= H5Sclose(ememspace);

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
        e = job->in_element[i];
        sxx_acc[e] += job->particles[i].v * job->particles[i].sxx;
        sxy_acc[e] += job->particles[i].v * job->particles[i].sxy;
        syy_acc[e] += job->particles[i].v * job->particles[i].syy;
        v_acc[e] += job->particles[i].v;
    }

    fprintf(fd, "%d %f %d\n", frame, time, job->num_elements);

    for (i = 0; i < job->num_elements; i++) {
        node_number_to_coords(&x, &y, job->elements[i].nodes[0], job->N, job->h);
        fprintf(fd, "%f %f ", x, y);
        node_number_to_coords(&x, &y, job->elements[i].nodes[1], job->N, job->h);
        fprintf(fd, "%f %f ", x, y);
        node_number_to_coords(&x, &y, job->elements[i].nodes[2], job->N, job->h);
        fprintf(fd, "%f %f ", x, y);
        node_number_to_coords(&x, &y, job->elements[i].nodes[3], job->N, job->h);
        fprintf(fd, "%f %f ", x, y);
        fprintf(fd, "%f %f %f\n",
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

    fprintf(fd, "%f %f %f\n", job->t, job->dt, job->t_stop);
    fprintf(fd, "%d %d %d\n", job->num_particles, job->num_nodes, job->num_elements);
    fprintf(fd, "%d %f\n", job->N, job->h);

    for (i = 0; i < job->num_particles; i++) {
        /* Position */
        fprintf(fd, "%f %f\n", job->particles[i].x, job->particles[i].y);

        /* Local Coordinates */
        fprintf(fd, "%f %f\n", job->particles[i].xl, job->particles[i].yl);

        /* Velocity */
        fprintf(fd, "%f %f\n", job->particles[i].x_t, job->particles[i].y_t);

        /* Mass */
        fprintf(fd, "%f\n", job->particles[i].m);

        /* Volume */
        fprintf(fd, "%f\n", job->particles[i].v);

        /* Initial volume */
        fprintf(fd, "%f\n", job->particles[i].v0);

        /* Stress */
        fprintf(fd, "%f %f %f\n",
            job->particles[i].sxx,
            job->particles[i].sxy,
            job->particles[i].syy);

        /* Stress rate */
/*        fprintf(fd, "%f %f %f\n",*/
/*            job->particles[i].sxx_t,*/
/*            job->particles[i].sxy_t,*/
/*            job->particles[i].syy_t);*/

        /* Strain rate */
        fprintf(fd, "%f %f %f %f\n",
            job->particles[i].exx_t,
            job->particles[i].exy_t,
            job->particles[i].eyy_t,
            job->particles[i].wxy_t);

        /* Jaumann stress increment */
/*        fprintf(fd, "%f %f %f\n",*/
/*            job->particles[i].dsjxx,*/
/*            job->particles[i].dsjxy,*/
/*            job->particles[i].dsjyy);*/

        /* Elastic and plastic strain increments */
/*        fprintf(fd, "%f %f %f\n",*/
/*            job->particles[i].deexx,*/
/*            job->particles[i].deexy,*/
/*            job->particles[i].deeyy);*/
/*        fprintf(fd, "%f %f %f\n",*/
/*            job->particles[i].depxx,*/
/*            job->particles[i].depxy,*/
/*            job->particles[i].depyy);*/

        /* Body forces */
        fprintf(fd, "%f %f\n", job->particles[i].bx, job->particles[i].by);

        /* Deformation gradient tensor */
        fprintf(fd, "%f %f %f %f\n",
            job->particles[i].Fxx, job->particles[i].Fxy,
            job->particles[i].Fyx, job->particles[i].Fyy);

        /* Initial particle domain vectors */
/*        fprintf(fd, "%f %f %f %f\n",*/
/*            job->particles[i].r1x0, job->particles[i].r1y0,*/
/*            job->particles[i].r2x0, job->particles[i].r2y0);*/

        /* Updated particle domain vectors */
/*        fprintf(fd, "%f %f %f %f\n",*/
/*            job->particles[i].r1xn, job->particles[i].r1yn,*/
/*            job->particles[i].r2xn, job->particles[i].r2yn);*/

        /* Corner positions */
/*        fprintf(fd, "%f %f %f %f %f %f %f %f\n",*/
/*            job->particles[i].c1x, job->particles[i].c1y,*/
/*            job->particles[i].c2x, job->particles[i].c2y,*/
/*            job->particles[i].c3x, job->particles[i].c3y,*/
/*            job->particles[i].c4x, job->particles[i].c4y);*/

        /* Displacements */
        fprintf(fd, "%f %f\n", job->particles[i].ux, job->particles[i].uy);

        /* Color used by splot visualization */
        fprintf(fd, "%f\n", job->particles[i].color);

        /* State Variables (for constitutive law) */
        fprintf(fd, "%d\n", DEPVAR);
        for (j = 0; j < DEPVAR; j++) {
            fprintf(fd, "%f\n", job->particles[i].state[j]);
        }

        /* Flag particle as active or not (used for discharge problems). */
        fprintf(fd, "%d\n", job->particles[i].active);
    }

    for (i = 0; i < job->num_nodes; i++) {
        fprintf(fd, "%f %f\n", job->nodes[i].x, job->nodes[i].y);
    }

    for (i = 0; i < job->num_elements; i++) {
         fprintf(fd, "%d %d %d %d %d %d %d %d %d\n",
            job->elements[i].nodes[0],
            job->elements[i].nodes[1],
            job->elements[i].nodes[2],
            job->elements[i].nodes[3],
            job->elements[i].nodes[4],
            job->elements[i].nodes[5],
            job->elements[i].nodes[6],
            job->elements[i].nodes[7],
            job->elements[i].nodes[8]);
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void h5_write_state(const char *hfilename, job_t *job)
{
    hid_t hfile;
    herr_t status = 0;
    hsize_t dims[1];    /* data by particle index */
    hsize_t chunk_dims[1];  /* for compression */

    hid_t hdataspace;
    hid_t hdataset;

    hid_t hplist;
    hid_t hparticle_type;

    dims[0] = job->num_particles;
    chunk_dims[0] = job->num_particles; /* one giant chunk */

    hfile = H5Fcreate(hfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hparticle_type = h5_particle_type();

    hdataspace = H5Screate_simple(sizeof(dims)/sizeof(dims[0]), dims, NULL);

    hplist  = H5Pcreate(H5P_DATASET_CREATE);
    status |= H5Pset_chunk(hplist, sizeof(chunk_dims)/sizeof(chunk_dims[0]), chunk_dims);
    status |= H5Pset_deflate(hplist, 9);
    hdataset = H5Dcreate2(hfile, "/particles", hparticle_type, hdataspace, 
                          H5P_DEFAULT, hplist, H5P_DEFAULT);

    status |= H5Dwrite(hdataset, hparticle_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                     job->particles);

    status |= H5Dclose(hdataset);
    status |= H5Sclose(hdataspace);
    status |= H5Fclose(hfile);

    return;
}
/*----------------------------------------------------------------------------*/

