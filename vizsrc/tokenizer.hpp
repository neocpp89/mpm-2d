#include <fstream>
#include <string>
#include <vector>

#ifndef __TOKENIZER_HPP__
#define __TOKENIZER_HPP__

class Tokenizer
{
    public:
        Tokenizer() { return; }
        static std::vector<std::string>
            splitNextLine(std::istream& str, const char delim = ' ');
        static std::vector<std::string>
            splitString(std::string& str, const char delim = ' ');
        template <typename Iterator>
        static Iterator splitNextLineIt(std::istream& strm, const char delim,
            Iterator begin, Iterator end = nullptr)
        {
            std::string line;
            std::getline(strm, line);
            std::stringstream lineStream(line);
            std::string cell;

            while(std::getline(lineStream,cell,delim) && begin != end) {
                *begin = cell;
                begin++;
            }

            return begin; // new end
        }
};

#endif

