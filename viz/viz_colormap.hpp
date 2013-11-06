/**
    \file viz_colormap.hpp
    \author Sachith Dunatunga
    \date 2013-09-12

    Colormap functions for visualization program.
*/

#ifndef __VIZ_COLORMAP_HPP__
#define __VIZ_COLORMAP_HPP__

typedef struct rgba_s {
    float r;
    float g;
    float b;
    float a;
} rgba_t;

typedef struct cm_s {
    rgba_t *colors;
    float *anchors;
    int num_colors;
} cm_t;

#endif

