/* Copyright (c) 2016 C. Pattison
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
 
#include "string_util.hpp"
#include <algorithm>
#include <iterator>
#include <iostream>

namespace psqa { namespace util  {

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
}}