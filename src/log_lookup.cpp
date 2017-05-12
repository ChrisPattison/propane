#include "log_lookup.hpp"
#include <cmath>
#include <limits>

namespace propane::util
{
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
}