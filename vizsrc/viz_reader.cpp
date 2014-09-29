#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <ios>

#include "tokenizer.hpp"
#include "viz_particle.hpp"
#include "viz_reader.hpp"

//std::ifstream::sync_with_stdio(false);
//std::ios::sync_with_stdio(false);

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
    tokens.resize(num_known_keys);

    std::vector<Particle> next;
    for (size_t j = 0; j < num_particles; j++) {
//        tokens = Tokenizer::splitNextLine(pfstream);
//        size_t num_keys = std::min(num_known_keys, tokens.size());

        auto end = Tokenizer::splitNextLineIt(pfstream, ' ', tokens.begin(), tokens.end());
        size_t num_keys = end - tokens.begin();

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

    }

    return next;
}

std::vector<Element> TXTReader::nextElements()
{
    std::vector<Element> next;
    return next;
}

std::vector<Particle> CSVReader::nextParticles()
{
    std::vector<Particle> next;

    if (frame < (ifp.getTotalFrames() - 1)) {
        std::string dataFilename = ifp.getFrameCSVFilename(frame);
        std::cout << "datafile: " << dataFilename << std::endl;
        std::ifstream dataStream(dataFilename);
        readSingleFrame(std::back_inserter(next), dataStream, true);
        frame++;
    }

    return next;
}

template <typename Iterator>
void CSVReader::readSingleFrame(Iterator particleIt, std::ifstream &dataStream, bool has_header)
{
    std::vector<std::string> header;
    bool strict = true;
    size_t line = 0;

    // first row is a label
    if (has_header) {
        header = Tokenizer::splitNextLine(dataStream, ',');
    } else {
        const std::string known_keys[] = {
            "m", "v", "x", "y", "x_t", "y_t", "sxx", "sxy", "syy", "ux", "uy",
            "gammap", "color", "magEf", "active", "corner_0x", "corner_0y",
            "corner_1x", "corner_1y", "corner_2x", "corner_2y",
            "corner_3x", "corner_3y"
        };
        const size_t num_known_keys = sizeof(known_keys) / sizeof(known_keys[0]);
        for (size_t i = 0; i < num_known_keys; i++) {
            header.push_back(known_keys[i]);
        }
    }

//    for (auto const &k : header) {
//        std::cout << k << std::endl;
//    }

    do {
        line++;
        std::vector<std::string> tokens = Tokenizer::splitNextLine(dataStream, ',');
        if (dataStream.eof() || !dataStream.good()) {
            break;
        }
        if (tokens.size() == 0 || tokens[0][0] == '#') {
            // comment or blank line; ignore it
            continue;
        }
        if (tokens.size() <= header.size()) {
            Particle p = Particle(line, true);
            for (size_t i = 0; i < tokens.size(); i++) {
                p[header[i]] = std::stod(tokens[i]);
            }

            if (p.keyExists("id")) {
                p.setId(static_cast<size_t>(p["id"]));
            }

            if (p.keyExists("active")) {
                if (p["active"] == 0) {
                    p.setInactive();
                } else {
                    p.setActive();
                }
            }

            *(particleIt++) = p;
        } else {
            if (strict) {
                std::cerr << "FATAL: line " << line;
                std::cerr << " has too many tokens (";
                std::cerr << tokens.size() << ") [";
                for (auto const &t : tokens) {
                    std::cerr << "," << t;
                }
                std::cerr<< " ]" << std::endl;
                exit(1);
            } else {
                std::cerr << "Ignoring line with " << tokens.size() << " tokens." << std::endl;
            }
        }
    } while (true);
}

std::vector<Particle> CSVReader::loadParticles(size_t frame)
{
    std::vector<Particle> next;
    return next;
}

std::vector<Element> CSVReader::nextElements()
{
    std::vector<Element> next;
    return next;
}

