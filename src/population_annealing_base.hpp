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
#include <vector>
#include <limits>

namespace propane  {
    
class PopulationAnnealingBase {
    public:
/** Observables for a single step.
 */
    struct Result {
        struct Histogram {
            double bin;
            double value;
        };
        std::vector<Histogram> overlap;
        std::vector<Histogram> link_overlap;
        std::vector<Histogram> energy_distribution;
        std::vector<Histogram> ground_distribution;
        double beta = std::numeric_limits<double>::quiet_NaN();
        int population = -1;
        double norm_factor = std::numeric_limits<double>::quiet_NaN();
        double average_energy = std::numeric_limits<double>::quiet_NaN();
        double average_squared_energy = std::numeric_limits<double>::quiet_NaN();
        double ground_energy = std::numeric_limits<double>::quiet_NaN();
        int grounded_replicas = 0;
        double entropy = std::numeric_limits<double>::quiet_NaN();
        std::uint64_t seed = 0;
        int sweeps = -1;
        double mean_square_family_size = std::numeric_limits<double>::quiet_NaN();
        double montecarlo_walltime = std::numeric_limits<double>::quiet_NaN();
        unsigned long long int total_sweeps = 0;
        double total_time = std::numeric_limits<double>::quiet_NaN();
    };
/** Parameters for a single step.
 */
    struct Schedule {
        double beta;
        int sweeps = 10;
        bool heat_bath = false;
        bool overlap_dist = false;
        bool energy_dist = false;
        bool ground_dist = false;
        bool compute_observables = false;
        double population_fraction = 1.0;
    };
/** Parameters for entire run.
 */
    struct Config {
        int population;
        std::uint64_t seed;
        std::vector<PopulationAnnealingBase::Schedule> schedule;
        bool solver_mode = false;
    };
};
}