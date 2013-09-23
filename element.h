/**
	\file element.h
	\author Sachith Dunatunga
	\date 04.06.12
	
	Contains the structure for elements in MPM.
*/
#ifndef __ELEMENT_H__
#define __ELEMENT_H__

#define NODAL_DOF 3
#define NODES_PER_ELEMENT 4

/*
    node numbering
    1 -> bottom left
    2 -> bottom right
    3 -> top right
    4 -> top left
    5 -> bottom edge
    6 -> right edge
    7 -> top edge
    8 -> left edge
    9 -> center
*/
typedef struct element_s {
    int nodes[9];
    int n;
    double m;
    int filled;
    int neighbors[8];
    double grad_x;
    double grad_y;
    double grad_mag;
    double n_x;
    double n_y;
    double n_theta;
    double v;

    int p_index;
    double p_dist;

    double sxx;
    double sxy;
    double syy;

    double kku_element[NODAL_DOF * NODES_PER_ELEMENT][NODAL_DOF * NODES_PER_ELEMENT]; /* dof order: x1, y1, x2, y2 ... x4 , y4 ... c1, c2, ... c4 */
    double f_element[NODAL_DOF * NODES_PER_ELEMENT];
} element_t;

/*
    Convert the node number into coordinates (x,y).
*/
void node_number_to_coords(double *x, double *y, int num, int N, double h);
#endif

