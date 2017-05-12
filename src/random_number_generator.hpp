#pragma once
#include <cstdlib>
#include <cstdint>
#include <type_traits>
#include <climits>
#include "dSFMT.h"

namespace propane
{
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
};
}