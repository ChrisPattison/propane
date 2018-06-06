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
 
#pragma once
#include <cstdlib>
#include <cstdint>
#include <type_traits>
#include <climits>
#include "dSFMT.h"

namespace psqa {
/** Wrapper around the refreence implementation of dSFMT19937 by Saito and Matsumoto.
 */
class RandomNumberGenerator {
    dsfmt_t* state_;
    std::uint64_t seed_;
public:
/** Seed generator with given value.
 */
    RandomNumberGenerator(std::uint64_t seed);
/** Seed generator with value from RandomSeed.
 */
    RandomNumberGenerator();
    ~RandomNumberGenerator();

    RandomNumberGenerator(const RandomNumberGenerator& other);

    RandomNumberGenerator(RandomNumberGenerator&& other);

    RandomNumberGenerator& operator=(const RandomNumberGenerator& other);

    RandomNumberGenerator& operator=(RandomNumberGenerator&& other);
/** Returns a random value suitable for seeding.
 */
    std::uint64_t RandomSeed();
/** Returns the value used to seed this instance.
 */
    std::uint64_t GetSeed();
/** Returns a double uniformly distributed in [0,1).
 */
    double Probability();
/** Returns an integer uniformly distributed in [0,N).
 */
    int Range(int N);
/** Return random 32 bits
 */
    std::uint32_t operator()();
};
}