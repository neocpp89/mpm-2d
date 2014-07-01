#include <fstream>
#include <string>
#include <vector>

#include "viz_particle.hpp"
#include "viz_element.hpp"

#ifndef __VIZ_READER_HPP__
#define __VIZ_READER_HPP__

class SimulationReader
{
    public:
        virtual ~SimulationReader() { return; }
        virtual std::vector<Particle> nextParticles() = 0;
        virtual std::vector<Element> nextElements() = 0;
        virtual double currentTime() = 0;
        virtual size_t currentFrame() = 0;
};

class TXTReader : public SimulationReader
{
    public:
        TXTReader(std::string _pf, std::string _ef) :
            particle_filename(_pf), element_filename(_ef),
            frame(0), time(0), state(FRAME)
        {
            pfstream.open(particle_filename);
            efstream.open(element_filename);
            return;
        }
        std::vector<Particle> nextParticles();
        std::vector<Element> nextElements();
        double currentTime() { return time; }
        size_t currentFrame() { return frame; }

    private:
        std::string particle_filename;
        std::string element_filename;

        double frame;
        double time;

        enum ReaderState { FRAME, PARTICLE, ELEMENT };
        enum ReaderState state;

        std::ifstream pfstream;
        std::ifstream efstream;

};

class Tokenizer
{
    public:
        Tokenizer() { return; }
        static std::vector<std::string> splitNextLine(std::istream& str, const char delim = ' ');
};

#endif

