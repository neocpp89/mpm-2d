/**
	\file implicit.h
	\author Sachith Dunatunga
	\date 20.10.2014
	
    Functions to create the stiffness matrix needed for implicit MPM.
*/
#ifndef __IMPLICIT_H__
#define __IMPLICIT_H__
#include "process.h"

void add_particle_stiffness(int idx, job_t *job);
void build_elemental_stiffness(job_t *job);
#endif
