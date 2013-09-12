/**
    \file visualization.h
    \author Sachith Dunatunga
    \date 02.06.12

    Writes particle data to disk.
*/
#include <stdio.h>

#ifndef __VISUALIZATION_H__
#define __VISUALIZATION_H__
#include "particle.h"
#include "point.h"
#include "process.h"

FILE *makeplot(int persist);
void plot_particles(FILE *pipe, job_t *job);
void closeplot(FILE *pipe);

#endif

