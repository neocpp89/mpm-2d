/**
    \file process.h
    \author Sachith Dunatunga
    \date 04.06.12

    Contains the structure for a constitutive law in MPM.
*/
#ifndef __MATERIAL_H__
#define __MATERIAL_H__
#include "process.h"

#define EMOD 1e6
#define NUMOD 0.3

//#define calculate_stress calculate_stress_dp_shearonly_indep
#define calculate_stress calculate_stress_dp_shearonly
#define material_tangent material_tangent_dp_shearonly_indep

enum e_material_types {M_RIGID=0, M_DRUCKER_PRAGER};

void material_init(job_t *job);
void calculate_stress(job_t *job);
void material_tangent(job_t *job);

#endif

