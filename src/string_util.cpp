#include "string_util.hpp"
#include <algorithm>
#include <iterator>
#include <iostream>

namespace propane::util 
{
void Check(bool cond, const char* err) {
    if(!cond) {
        std::cout << err << std::endl;
        exit(EXIT_FAILURE);
    }
}

std::vector<std::string> Split(std::string target, char seperator, StringSplitOptions options) {
    std::vector<std::string> out;
    std::string::iterator substring_begin, substring_end;
    substring_begin = target.begin();
    
    if(*substring_begin == seperator) {
        ++substring_begin;
    }

    do {
        substring_end = std::find(substring_begin, target.end(), seperator);
        std::string substring(std::distance(substring_begin, substring_end), ' ');
        std::copy(substring_begin, substring_end, substring.begin());

        if(substring.size() != 0 || options == kStringSplitOptions_None) {
            out.push_back(substring);
        }

        substring_begin = substring_end;
        ++substring_begin;
    } while(substring_end != target.end());

    return out;
}

std::vector<std::string> Split(std::string target, std::vector<char> seperator, StringSplitOptions options) {
    std::vector<std::string> out;
    std::string::iterator substring_begin, substring_end;
    substring_begin = target.begin();
    
    for(auto c : seperator) {
        if(*substring_begin == c) {
            ++substring_begin;
            break;
        }
    }
    
    do {
        substring_end = std::find_if(substring_begin, target.end(), [&](char c)
            { return std::find(seperator.begin(), seperator.end(), c) != seperator.end(); });
        std::string substring(std::distance(substring_begin, substring_end), ' ');
        std::copy(substring_begin, substring_end, substring.begin());

        if(substring.size() != 0 || options == kStringSplitOptions_None) {
            out.push_back(substring);
        }

        substring_begin = substring_end;
        ++substring_begin;
    } while(substring_end != target.end());

    return out;
}
}