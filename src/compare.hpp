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

namespace psqa { namespace util {
/** Fuzzy Comparison of floats by units in last place.
 * Adjacent floats have adjacent integer representation.
 * The different in this representation is a scalable way to compare for equality.
 * Compares the absolute value of the difference of the integer representation of the inputs to err.
 */
bool FuzzyUlpCompare(const float& a, const float& b, const int err = 10);

/** Double precision version of FuzzyULPCompare.
 */
bool FuzzyUlpCompare(const double& a, const double& b, const int err = 10);

/** Comparison of floats by an epsilon value.
 * ULP comparison breaks down near zero so comparison by a small number is necessary.
 * compares that absolute value of the different to the epsilon value.
 */
bool FuzzyEpsCompare(const float& a, const float& b, const float err = 1e-6);

/** double version of FuzzyEpsCompare.
 */
bool FuzzyEpsCompare(const double& a, const double& b, const double err = 1e-15);

/** Double version of FuzzyCompare.
 */
bool FuzzyCompare(const double& a, const double& b);

/** General purpose fuzzy comparison.
 * Logical OR of the outputs of Epsilon and ULP comparison.
 * More permissive than either one.
 */
bool FuzzyCompare(const float& a, const float& b);
}}