/**
    \file viz_colormap.cpp
    \author Sachith Dunatunga
    \date 2013-09-12

    Colormap functions for visualization program.
*/
#include "viz_colormap.hpp"

std::ostream& operator<<(std::ostream& os, const rgbaColor& c)
{
    os << "[r:" << c.r << ", g:" << c.g << ", b:" << c.b << ", a:" << c.a << "]";
    return os;
}

void Colormap::parseColormap(std::istream &str)
{
    std::vector<std::string> tokens;
    char type;
    float control_point;
    float r, g, b, a;
    long hexcolor;
    
    do {
        tokens.clear();
        tokens = Tokenizer::splitNextLine(str, ',');
        if(str.eof() || tokens.size() <= 0) {
            break;
        }
        type = tokens[0][0];
        switch (type) {
            case 'F':
                if (tokens.size() != 6) {
                    continue;
                }
                control_point = stof(tokens[1]);
                r = stof(tokens[2]);
                g = stof(tokens[3]);
                b = stof(tokens[4]);
                a = stof(tokens[5]);
                break;
            case 'I':
                if (tokens.size() != 6) {
                    continue;
                }
                control_point = stof(tokens[1]);
                r = stoi(tokens[2]) / 255.0;
                g = stoi(tokens[3]) / 255.0;
                b = stoi(tokens[4]) / 255.0;
                a = stoi(tokens[5]) / 255.0;
                break;
            case 'H':
                if (tokens.size() != 3) {
                    continue;
                }
                control_point = stof(tokens[1]);
                hexcolor = stoul(tokens[2]);
                r = ((hexcolor >> 24) & 0xFF) / 255.0;
                g = ((hexcolor >> 16) & 0xFF) / 255.0;
                b = ((hexcolor >> 8) & 0xFF) / 255.0;
                a = ((hexcolor >> 0) & 0xFF) / 255.0;
                break;
            case '#':
            default:
                continue;
        }
        colormap[control_point] = rgbaColor(r,g,b,a);
    } while (true);
    return;
}

rgbaColor Colormap::interpolatedColor(float point)
{
    rgbaColor result;

    if (colormap.empty()) {
        return result;
    }

    float minpt = colormap.begin()->first;
    float maxpt = colormap.rbegin()->first;

    if (point <= minpt) {
        result = rgbaColor(colormap.begin()->second);
    } else if (point >= maxpt) {
        result = rgbaColor(colormap.rbegin()->second);
    } else {
        auto lb = colormap.lower_bound(point);
        float leftpt = lb->first;
        rgbaColor left(lb->second);
        lb++;
        float rightpt = lb->first;
        rgbaColor right(lb->second);
        float diff = rightpt - leftpt;
        result = rgbaColor(
            ((point - leftpt) / diff) * left.getr() + ((rightpt - point) / diff) * right.getr(),
            ((point - leftpt) / diff) * left.getg() + ((rightpt - point) / diff) * right.getg(),
            ((point - leftpt) / diff) * left.getb() + ((rightpt - point) / diff) * right.getb(),
            ((point - leftpt) / diff) * left.geta() + ((rightpt - point) / diff) * right.geta()
        );
    }

    return result;
}

