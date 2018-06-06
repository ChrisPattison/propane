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
#include <utility>
#include <algorithm>
#include "types.hpp"
#include "graph.hpp"
#include "random_number_generator.hpp"
#include "population_annealing_base.hpp"
#include "log_lookup.hpp"

namespace psqa {
/** Implementation of Population Annealing Monte Carlo.
 * Replicas have an associated entry in the family vector indicating lineage.
 */
class PopulationAnnealing : public PopulationAnnealingBase {
protected:

    util::LogLookup log_lookup_;

    Graph structure_;

    RandomNumberGenerator rng_;

    using StateVector = std::vector<VertexType>;
    std::vector<StateVector> replicas_;
    std::vector<int> replica_families_;

    int init_population_;
    int average_population_;
    int trotter_slices_;
    std::vector<Schedule> schedule_;
    // This should have getters and setters
    // Write only through setparams
    double beta_;
    double coeff_P_;
    double coeff_D_;
    bool solver_mode_;

/** Gets the number of replicas in each family as a fraction of the total population
 */
    std::vector<double> FamilyCount();
/**
 * Uses a look up table to compute a bound on the logarithm of a random number 
 * and compares to the exponent of the acceptance probability.
 * If the probability is inside the bound given by the look table, 
 * true exponential is computed and compared.
 */
    bool AcceptedMove(double log_probability);
/** Returns the energy of a replica
 * Implemented as the sum of elementwise multiplication of the replica vector with the 
 * product of matrix multiplication between the upper half of the adjacency matrix
 * and the replica.
 */
    double Hamiltonian(const StateVector& replica);
/** Returns problem energy without weight factor (zeta)
 */
    double ProblemHamiltonian(const StateVector& replica);
/** Returns driver energy without weight factor (gamma)
 */
    double DriverHamiltonian(const StateVector& replica);
/** Carries out moves*quantum sites Wolff cluster moves of replica
 */
    void WolffSweep(StateVector& replica, std::size_t moves);
/** Resamples population according to the Boltzmann distribution.
 * Attempts to maintain approximately the same population as detailed in arXiv:1508.05647
 * Returns the normalization factor Q as a byproduct.
 */
    double Resample(double new_beta, double new_gamma, double new_lambda, double new_population_fraction);
/** Changes coeff_P, coeff_D, and beta for a particular transverse field
 * Do not change beta and gamma directly
 */
    void SetParams(double new_beta, double new_gamma, double new_lambda);
public:

    PopulationAnnealing() = delete;
/** Intializes solver.
 * The inputs will be replaced by a struct in the future.
 * schedule specifies the annealing schedule, sweep counts, and histogram generation at each step.
 * seed may be zero in which case one will be generated.
 */
    PopulationAnnealing(Graph& structure, Config schedule);
/** Run solver and return results.
 */
    std::vector<Result> Run();
/** Return the spatial site energy for each trotter slice
 */
template<typename random_access_iterator>
void SpatialSiteEnergy(const StateVector& replica, std::uint32_t site, random_access_iterator it) {
    std::fill(it, it + ktrotter_slices, 0);
    for(std::size_t edge = 0; edge < structure_.adjacent()[site].size(); ++edge) {
        auto parity = replica[site] ^ replica[structure_.adjacent()[site][edge]];
        auto weight = structure_.weights()[site][edge];
        auto temporal = it;
        EvalFunctor(parity, [&](const auto& v) { *temporal++ += v*weight; });
    }

    auto temporal = it;
    auto field = structure_.fields()[site];
    EvalFunctor(replica[site], [&](const auto& v) { *temporal++ += v*field; });
}
/** Return the problem energy for each trotter slice
 */
template<typename random_access_iterator>
void SpatialProblemHamiltonain(const StateVector& replica, random_access_iterator it) {
    std::fill(it, it + ktrotter_slices, 0);
    for(std::size_t site = 0; site < structure_.size(); ++site) {
        for(std::size_t edge = 0; edge < structure_.adjacent()[site].size(); ++edge) {
            auto parity = replica[site] ^ replica[structure_.adjacent()[site][edge]];
            auto weight = structure_.weights()[site][edge];
            auto temporal = it;
            EvalFunctor(parity, [&](const auto& v) { *temporal++ += v*weight; });
        }
    }

    std::transform(it, it + ktrotter_slices, it, [](const auto& v) { return v / 2; });

    for(std::size_t site = 0; site < structure_.size(); ++site) {
        auto field = structure_.fields()[site];
        auto temporal = it;
        EvalFunctor(replica[site], [&](const auto& v) { *temporal++ += v*field; });
    }
}

};
}