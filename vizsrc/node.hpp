#ifndef __NODE_HPP__
#define __NODE_HPP__

class Node {
    private:
        double sxx, sxy, syy;
        double m;
        double x, y;
        double w, h;

    public:
        Node(double X, double Y, double W, double H) :
            sxx(0), sxy(0), syy(0), m(0), x(X), y(Y), w(W), h(H) { return; }
        ~Node() { return; }

        double ShapeFunction(double xp, double yp)
        {
            double xl, yl, s;
            xl = (xp - x) / w;
            yl = (yp - y) / h;

            s = (1 - fabs(xl)) * (1 - fabs(yl));

            if (s < 0) {
                s = 0;
            }

            return s;
        }
        void clear()
        {
            ClearStress();
            return;
        }
        void ClearStress()
        {
            sxx = 0;
            sxy = 0;
            syy = 0;
            m = 0;
            return;
        }
        void CopyStressFromParticle(aux_particle_t *p)
        {
            double s = ShapeFunction(p->x, p->y);
            this->sxx += s * p->sxx * p->m;
            this->sxy += s * p->sxy * p->m;
            this->syy += s * p->syy * p->m;
            this->m += s * p->m;
            return;
        }
        void CopyStressToParticle(aux_particle_t *p)
        {
            double s = ShapeFunction(p->x, p->y);
            p->sxx += s * this->sxx;
            p->sxy += s * this->sxy;
            p->syy += s * this->syy;
            return;
        }
        void Print()
        {
            printf("%p: %g %g %g %g %g %g\n", this, m, x, y, sxx, sxy, syy);
            fflush(stdout);
            return;
        }
        void RescaleStress()
        {
            if (m != 0) {
                this->sxx = this->sxx / this->m;
                this->sxy = this->sxy / this->m;
                this->syy = this->syy / this->m;
            }
        }
        friend std::ostream& operator<<(std::ostream& os, const Node& c)
        {
            os << "[ m: " << c.m << " x: " << c.x << " y: " << c.y <<
                " sxx: " << c.sxx << " sxy: " << c.sxy <<
                " syy: " << c.syy << " ]";
            return os;
        }
};

#endif //__NODE_HPP__

