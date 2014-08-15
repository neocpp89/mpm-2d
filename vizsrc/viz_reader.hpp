#include <fstream>
#include <string>
#include <vector>

#include "infoparser.hpp"
#include "tokenizer.hpp"
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

class RandomAccessSimulationReader: public SimulationReader
{
    public:
        virtual ~RandomAccessSimulationReader() { return; }
        virtual std::vector<Particle> loadParticles(size_t frame) = 0;
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

class CSVReader : public RandomAccessSimulationReader
{
    public:
        CSVReader(std::string const & _infoFile) :
            infoFile(_infoFile)
        {
            infoStream.open(infoFile);
            return;
        }
        std::vector<Particle> nextParticles();
        std::vector<Element> nextElements();

    private:
        std::string infoFile;
        std::ifstream infoStream;
};
#endif

