#pragma once
#include <vector>

namespace propane::util
{
/** Utility to quickly obtain a bound for the natural logarithm of a value
 */
class LogLookup {
    int const kLookupTableSize = 1024;
    std::vector<double> lookup_table_;
public:
    struct Bound {
        double upper;
        double lower;
    };
/** Allocates and fills lookup table
 */
    LogLookup();
/** Returns the upper and lower bound of the natural log of a number in [0,1]
 * Returns +/- inf if out of bounds
 */
    Bound GetBound(double value);

    Bound operator()(double value);
};
}