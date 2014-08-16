#include <vector>
#include "node.hpp"

#ifndef __BACKGROUNDMESH_HPP__
#define __BACKGROUNDMESH_HPP__

class BackgroundMesh {
    private:
        std::vector<Node> nodes;
        size_t xdim;
        size_t ydim;
        double width;
        double height;

        size_t linearIndex(size_t x, size_t y)
        {
            return (y * xdim + x);
        }

    public:
        BackgroundMesh(size_t _xdim = 0, size_t _ydim = 0, double _w = 1.0, double _h = 1.0)
            : xdim(_xdim), ydim(_ydim), width(_w), height(_h)
        {
            double dx, dy;
            nodes.reserve(xdim*ydim);
            dx = width / (xdim - 1);
            dy = height / (ydim - 1);

            for (size_t j = 0; j < ydim; j++) {
                for (size_t i = 0; i < xdim; i++) {
                    nodes.push_back(Node(i * dx, j * dy, dx, dy));
                }
            }

            return;
        }

        void clear()
        {
            for (auto &n : nodes) {
                n.clear();
            }
            return;
        }

        void rescaleStress()
        {
            for (auto &n : nodes) {
                n.RescaleStress();
            }
            return;
        }

        size_t numNodes()
        {
            return nodes.size();
        }

        template <typename Integral>
        Node &operator()(Integral x, Integral y)
        {
            return nodes[linearIndex(x, y)];
        }

        template <typename Integral> 
        Node &operator[](Integral lin)
        {
            return nodes[lin];
        }
};

#endif //__BACKGROUNDMESH_HPP__

