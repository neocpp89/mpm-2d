/**
    \file viz_colormap.hpp
    \author Sachith Dunatunga
    \date 2013-09-12

    Colormap functions for visualization program.
*/
#include <map>
#include <sstream>
#include <fstream>
#include <vector>

#include "viz_reader.hpp"
#include "viz_builtin_colormap.hpp"

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

class rgbaColor
{
    private:
        float r, g, b, a;

    public:
        rgbaColor(rgbaColor const &c) { r = c.getr(); g = c.getg(); b = c.getb(); a = c.geta();}
        rgbaColor(float _r = 0.0f, float _g = 0.0f, float _b = 0.0f, float _a = 1.0f) :
            r(_r), g(_g), b(_b), a(_a) { return; }
        void set(float _r = 0.0f, float _g = 0.0f, float _b = 0.0f, float _a = 1.0f)
        {
            r = _r; g = _g; b = _b; a = _a;
            return;
        }
        float getr() const { return r; }
        float getg() const { return g; }
        float getb() const { return b; }
        float geta() const { return a; }
        void invert() { r = 1.0 - r; g = 1.0 - g; b = 1.0 - b; a = 1.0 - a; return; }

        friend std::ostream& operator<<(std::ostream& os, const rgbaColor& c);
};

class Colormap
{
    private:
        std::map<float, rgbaColor> colormap;
        std::ifstream cmapfile;

        void parseColormap(std::istream &str);

    public:
        Colormap() { std::stringstream s; s << viz_default_colormap_str; parseColormap(s); return; }
        Colormap(std::string colormap_file) { cmapfile.open(colormap_file); parseColormap(cmapfile); return; }
        friend std::ostream& operator<<(std::ostream& os, const Colormap& c)
        {
            os << "colormap = { ";
            for (auto const &kv: c.colormap) {
                os << kv.first << ":" << kv.second << " ";
            }
            os << "}";
            return os;
        }

        rgbaColor interpolatedColor(float point);
};

#endif

