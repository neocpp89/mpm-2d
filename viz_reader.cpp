#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>

#include "viz_reader.hpp"

/* from http://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c */
std::vector<std::string> Tokenizer::splitNextLine(std::istream& str, const char delim)
{
    std::vector<std::string>   result;
    std::string line;
    std::getline(str,line);

    std::stringstream lineStream(line);
    std::string cell;

    while(std::getline(lineStream,cell,delim)) {
        result.push_back(cell);
    }

    return result;
}

std::vector<Particle> TXTReader::nextParticles()
{
    size_t num_particles = 0;

    /*
        Important: this has to match the order we wrote out the fields. I made
        a mistake in the original design in that the field documentation is
        not written along with the file. However, the fields written are
        luckily a subset of the fields below (with order preserved). The second
        version should avoid this problem...
    */
    const std::string known_keys[] = {
        "m", "v", "x", "y", "x_t", "y_t", "sxx", "sxy", "syy", "ux", "uy",
        "gammap", "color", "magEf", "active", "corner_0x", "corner_0y",
        "corner_1x", "corner_1y", "corner_2x", "corner_2y",
        "corner_3x", "corner_3y"
    };

    std::vector<std::string> tokens = Tokenizer::splitNextLine(pfstream);
    if (tokens.size() != 3) {
        //TODO: throw exception
    }
    frame = std::stoull(tokens[0]);
    time = std::stod(tokens[1]);
    num_particles = std::stoull(tokens[2]);

    const size_t num_known_keys = sizeof(known_keys) / sizeof(known_keys[0]);
    std::vector<Particle> next;
    for (size_t j = 0; j < num_particles; j++) {
        tokens.clear();
        tokens = Tokenizer::splitNextLine(pfstream);
        size_t num_keys = std::min(num_known_keys, tokens.size());

        /*
            Default to active state, we will correct later if the version of
            mpm we are using writes this flag.
        */
        next.push_back(Particle(j, true));
        Particle &p = next.back();
        for (size_t k = 0; k < num_keys; k++) {
            p[known_keys[k]] = std::stod(tokens[k]);
        }
        
        if (p.keyExists("active") && p["active"] == 0) {
                p.setInactive();
        }

        /* Compute a few special values, if the prereqs exist. */
        p.calculatePressure();
        p.calculateDeviator();
        p.calculateEquivalentShear();
        p.calculateFrictionCoefficient();
        p.calculateStressEigenvalues();
        p.calculateStressEigenvectors();

//        if (p.keyExists("sxx") && p.keyExists("sxy") && p.keyExists("syy")) {
//            p["p"] = -0.5 * (p["sxx"] + p["syy"]);
//            p["s0xx"] = p["sxx"] + p["p"];
//            p["s0xy"] = p["sxy"];
//            p["s0yy"] = p["syy"] + p["p"];
//            p["tau"] = std::sqrt(0.5 * (p["s0xx"]*p["s0xx"] + 2*p["s0xy"]*p["s0xy"] + p["s0yy"]*p["s0yy"]));
//            if (p["p"] > 0) {
//                p["mu"] = p["tau"] / p["p"];
//            } else {
//                p["mu"] = 1e16;
//            }
//        }

//        std::cout << p << std::endl;

//        p.print();
    }
//    std::cout << "Current Frame: " << frame << std::endl;
//    std::cout << "Current Time: " << time << std::endl;
    return next;
}

std::vector<Element> TXTReader::nextElements()
{
    std::vector<Element> next;
    return next;
}

