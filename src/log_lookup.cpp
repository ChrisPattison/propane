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
 
#include "log_lookup.hpp"
#include <cmath>
#include <limits>

namespace psqa { namespace util {
    
LogLookup::LogLookup() {
    lookup_table_.resize(kLookupTableSize);
    for(int k = 0; k < lookup_table_.size(); ++k) {
        lookup_table_[k] = std::log(static_cast<double>(k) / (lookup_table_.size() - 1));
    }
}

LogLookup::Bound LogLookup::GetBound(double value) {
    int lower_index = std::floor(value * (lookup_table_.size() - 1));
    // Bounds checking
    if(lower_index < 0 || lower_index >= lookup_table_.size()-1) {
        return {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    }
    return {lookup_table_[lower_index+1], lookup_table_[lower_index]};
}

LogLookup::Bound LogLookup::operator()(double value) {
    return GetBound(value);
}
}}