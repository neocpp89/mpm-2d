/**
    \file node.h
    \author Sachith Dunatunga
    \date 28.07.12

    Contains the structure for nodes in MPM.
*/
#ifndef __NODE_H__
#define __NODE_H__

typedef struct node_s {
    /* Filled element neighbors. */
    /* int element_neighbors[4];
    int filled_element_neighbors[4]; */
    int num_filled_element_neighbors;
    double mass_filled_element_neighbors;

    /* Mass */
    double m;

    /* Position */
    double x;
    double y;

    /* Displacemnent */
    double ux;
    double uy;

    /* Velocity */
    double x_t;
    double y_t;

    /* Acceleration */
    double x_tt;
    double y_tt;

    /* Momentum */
    double mx_t;
    double my_t;

    /* "pseudoforce" */
    double mx_tt;
    double my_tt;

    /* Force */
    double fx;
    double fy;

    /* Density */
    double rho;

    /* marker for displacment/velocity update */
    int velocity_update_flag;
} node_t;

#endif

