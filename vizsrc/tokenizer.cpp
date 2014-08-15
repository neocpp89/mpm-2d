#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "tokenizer.hpp"

/* from http://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c */
std::vector<std::string> Tokenizer::splitNextLine(std::istream& strm, const char delim)
{
    std::string line;
    std::getline(strm,line);

    return splitString(line, delim);
}

std::vector<std::string> Tokenizer::splitString(std::string& str, const char delim)
{
    std::vector<std::string> result;
    std::stringstream lineStream(str);
    std::string cell;

    while(std::getline(lineStream,cell,delim)) {
        result.push_back(cell);
    }

    return result;
}

