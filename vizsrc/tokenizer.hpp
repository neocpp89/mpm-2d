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
};

#endif

