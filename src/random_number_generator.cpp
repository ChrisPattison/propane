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
 
#include "random_number_generator.hpp"
#include <random>
#include <thread>
#include <ctime>
#include <chrono>
#include <cassert>

namespace psqa {
    
RandomNumberGenerator::RandomNumberGenerator(std::uint64_t seed) {
    seed_ = seed;
    state_ = new dsfmt_t;
    dsfmt_init_gen_rand(state_, static_cast<std::uint32_t>(seed));
}

RandomNumberGenerator::RandomNumberGenerator() : RandomNumberGenerator(RandomSeed()) {};

RandomNumberGenerator::~RandomNumberGenerator() {
    delete state_;
    state_ = nullptr;
}

RandomNumberGenerator::RandomNumberGenerator(const RandomNumberGenerator& other) {
    seed_ = other.seed_;
    state_ = new dsfmt_t;
    *state_ = *(other.state_);
}

RandomNumberGenerator::RandomNumberGenerator(RandomNumberGenerator&& other) {
    if(this != &other) {
        seed_ = other.seed_;
        state_ = other.state_;
        other.state_ = nullptr;
    }
}

RandomNumberGenerator& RandomNumberGenerator::operator=(const RandomNumberGenerator& other) {
    seed_ = other.seed_;
    *state_ = *(other.state_);
    return *this;
}

RandomNumberGenerator& RandomNumberGenerator::operator=(RandomNumberGenerator&& other) {
    if(this != &other) {
        delete state_;
        seed_ = other.seed_;
        state_ = other.state_;
        other.state_ = nullptr;
    }
    return *this;
}

double RandomNumberGenerator::Probability() {
    return dsfmt_genrand_close_open(state_);
}

std::uint64_t RandomNumberGenerator::RandomSeed() {
    return time(NULL) ^
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count() ^
        std::hash<std::thread::id>()(std::this_thread::get_id()) ^
        rand() ^
        std::random_device()();
}

std::uint64_t RandomNumberGenerator::GetSeed() {
    return seed_;
}

int RandomNumberGenerator::Range(int N) {
    return std::floor(Probability() * N);
}

std::uint32_t RandomNumberGenerator::operator()() {
    return dsfmt_genrand_uint32(state_);
}
}