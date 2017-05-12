#include "random_number_generator.hpp"
#include <random>
#include <thread>
#include <ctime>
#include <chrono>
#include <cassert>

namespace propane
{
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
}