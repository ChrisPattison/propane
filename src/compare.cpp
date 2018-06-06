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
 
#include <cstdint>
#include <cmath>

namespace psqa { namespace util {
    
bool FuzzyUlpCompare(const float& a, const float& b, const int err = 10) {
    if(std::signbit(a)!=std::signbit(b)) {
        return a == b;
    }
    auto ai = reinterpret_cast<const std::uint32_t&>(a);
    auto bi = reinterpret_cast<const std::uint32_t&>(b);

    // std::abs(ai-bi) is ambiguous with clang
    return (ai > bi ? ai-bi : bi-ai) < err;
}

bool FuzzyUlpCompare(const double& a, const double& b, const int err = 10) {
    if(std::signbit(a)!=std::signbit(b)) {
        return a == b;
    }
    auto ai = reinterpret_cast<const std::uint64_t&>(a);
    auto bi = reinterpret_cast<const std::uint64_t&>(b);

    return (ai > bi ? ai-bi : bi-ai) < err;
}

bool FuzzyEpsCompare(const float& a, const float& b, const float err = 1e-6) {
    return std::abs(a - b) < err;
}

bool FuzzyEpsCompare(const double& a, const double& b, const double err = 1e-15) {
    return std::abs(a - b) < err;
}

bool FuzzyCompare(const double& a, const double& b) {
    return FuzzyEpsCompare(a,b) || FuzzyUlpCompare(a,b);
}

bool FuzzyCompare(const float& a, const float& b) {
    return FuzzyEpsCompare(a,b) || FuzzyUlpCompare(a,b);
}
}}