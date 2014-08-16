#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <iostream>

#ifndef __VIZ_PARTICLE_HPP__
#define __VIZ_PARTICLE_HPP__

#ifndef M_PI_2
    const double pi = 4 * std::atan(1.0);
    #define M_PI_2 (pi / 2)
#endif

//TODO: add exceptions for key errors

class Particle {
    private:
        int id;
        bool active;
        std::unordered_map<std::string,double> data;

        double eigvals[2] = {0};
        double eigangs[2] = {0}; //eigenvector angles with horizontal
        double eigvecs[2][2] = {{0}, {0}};

        float rgbcolor[3]; //for use with opengl

//        double velocity_mag;
//        double velocity_angle;

    public:
        void setId(int _id) { id = _id; return; }

        Particle(int _id, bool _active) : id(_id), active(_active) { return; }

        friend std::ostream& operator<<(std::ostream& os, Particle& p)
        {
            std::vector<std::string> keys;
            keys.reserve(p.data.size());

            os << p.id << (p.active?": { ":"[inactive]: { ");
            for (int i = 0; i < 2; i++) {
                os << "stress-eigenvalue[" << i << "]:" << p.eigvals[i] << " ";
                os << "stress-eigenvector-angle[" << i << "]:" << p.eigangs[i] << " ";
            }
            os << "data = { ";
            for (auto const &kv: p.data) {
                keys.push_back(kv.first);
            }
            std::sort(keys.begin(), keys.end());
            for (auto const &k: keys) {
                os << k << ":" << p.getValue(k) << " ";
            }
            os << "} ";
            os << "}";

            return os;
        }

        bool keyExists(std::string key) const { return (data.count(key) == 1); }
        void setKeyValuePair(std::string key, double val) { data[key] = val; }
        double getValue(std::string key) { if (keyExists(key)) { return data[key]; } return 0; }
        double &operator[](const std::string key) { return data[key]; }
        double operator[](const std::string key) const { return data.at(key); }

        void setActive() { active = true; return; }
        void setInactive() { active = false; return; }
        bool isActive() const { return active; }

        bool calculatePressure()
        {
            bool set = false;
            if (keyExists("sxx") && keyExists("syy")) {
                data["p"] = -0.5 * (data["sxx"] + data["syy"]);
                set = true;
            }
            return set;
        }

        bool calculateDeviator()
        {
            bool set = false;
            if (keyExists("sxx") && keyExists("sxy") && keyExists("syy") && keyExists("p")) {
                data["s0xx"] = data["sxx"] + data["p"];
                data["s0xy"] = data["sxy"];
                data["s0yy"] = data["syy"] + data["p"];
                set = true;
            }
            return set;
        }

        bool calculateEquivalentShear()
        {
            bool set = false;
            if (keyExists("s0xx") && keyExists("s0xy") && keyExists("s0yy")) {
                data["tau"] = std::sqrt(0.5 * (data["s0xx"]*data["s0xx"]
                    + 2*data["s0xy"]*data["s0xy"]
                    + data["s0yy"]*data["s0yy"]));
                set = true;
            }
            return set;
        }

        bool calculateFrictionCoefficient()
        {
            bool set = false;
            if (keyExists("tau") && keyExists("p")) {
                if (data["p"] > 0) {
                    data["mu"] = data["tau"] / data["p"];
                } else {
                    data["mu"] = 1e16;
                }
                set = true;
            }
            return set;
        }

        bool calculateStressEigenvalues()
        {
            bool set = false;
            if (keyExists("sxx") && keyExists("sxy") && keyExists("syy")) {
                double a, b, c, p, n;
                a = 1.0;
                b = data["sxx"] + data["syy"];
                c = data["sxx"]*data["syy"] - data["sxy"]*data["sxy"];
                if (b >= 0) {
                    double s = (-b - std::sqrt(b * b - 4 * a * c));
                    p = (2.0 * c) / s;
                    n = s / (2.0 * a);
                } else {
                    double s = (-b + std::sqrt(b * b - 4 * a * c));
                    p = s / (2.0 * a);
                    n = (2.0 * c) / s;
                }
                eigvals[0] = p;
                eigvals[1] = n;
                set = true;
            }
            return set;
        }

        bool calculateStressEigenvectors()
        {
            bool set = false;
            if (keyExists("sxx") && keyExists("sxy") && keyExists("syy")) {
                double x, y;
                y = eigvals[0] - data["sxx"];
                x = data["sxy"];
                eigangs[0] = std::atan2(y, x);
                eigangs[1] = M_PI_2 + eigangs[0];
                eigvecs[0][0] = cos(eigangs[0]);
                eigvecs[0][1] = sin(eigangs[0]);
                eigvecs[1][0] = cos(eigangs[1]);
                eigvecs[1][1] = sin(eigangs[1]);
                set = true;
            }
            return set;
        }
};

#endif

