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
#include <cstdint>
#include <type_traits>

namespace psqa {

#ifndef PSQA_TROTTERSLICES
    #warning "PSQA_TROTTERSLICES not defined. Assuming 64. Valid values are {16, 32, 64}"
    #define PSQA_TROTTERSLICES 64
#endif

static constexpr std::uint32_t ktrotter_slices = PSQA_TROTTERSLICES;

using VertexType = 
    std::conditional_t<ktrotter_slices == 32, std::uint32_t,
    std::conditional_t<ktrotter_slices == 16, std::uint16_t, 
    std::uint64_t>>;
static_assert(ktrotter_slices == 16 || ktrotter_slices == 32 || ktrotter_slices == 64, "PSQA_TROTTERSLICES must be in {16, 32, 64}");

using EdgeType = double;
using IndexType = std::size_t;
const double kEpsilon = 1e-13;


/** Map bit masked value to +/- 1
 */
template<typename value_type>
auto GetValue(value_type v) { return v ? 1 : -1; }
/** Call functor on each value of the multi spin coded bit vector
 */
template<typename functor, typename bitvector, typename accumulate>
auto EvalFunctor(bitvector value, functor function, accumulate init) {
    auto accumulator = init;
    for(std::uint32_t k = 0; k < ktrotter_slices; ++k) {
        accumulator += function(GetValue(value & 1U));
        value >>= 1;
    }
    return accumulator;
}
/** Call functor on each value of the multi spin coded bit vector
 */
template<typename functor, typename bitvector>
void EvalFunctor(bitvector value, functor function) {
    for(std::uint32_t k = 0; k < ktrotter_slices; ++k) {
        function(GetValue(value & 1U));
        value >>= 1;
    }
}
/** Rotate value
 */
template<typename value_type, typename = std::enable_if_t<std::is_unsigned<value_type>::value>>
value_type RotateL(value_type x, int n) {
    return (x << n) | (x >> ( -n & (ktrotter_slices - 1)));
}

template<typename value_type, typename = std::enable_if_t<std::is_unsigned<value_type>::value>>
value_type RotateR(value_type x, int n) {
    return (x >> -n) | (x << ( n & (ktrotter_slices - 1)));
}
template<typename value_type>
value_type PopCount(value_type x) {
    return __builtin_popcount(x);
}
}