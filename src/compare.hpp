#pragma once

namespace propane::util
{
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
}