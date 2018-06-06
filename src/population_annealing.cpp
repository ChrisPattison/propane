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
 
#include "population_annealing.hpp"
#include "compare.hpp"
#include "log_sum_exp.hpp"
#include <cmath>
#include <cassert>
#include <iostream>
#include <limits>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <chrono>
#include <array>

namespace psqa {
std::vector<double> PopulationAnnealing::FamilyCount() {
    std::vector<double> count;
    count.reserve(replica_families_.size());
    auto i = replica_families_.begin();
    do {
        auto i_next = std::find_if(i, replica_families_.end(), [&](const int& v){return v != *i;});
        count.push_back(static_cast<double>(std::distance(i, i_next)));
        i = i_next;
    }while(i != replica_families_.end());
    return count;
}

PopulationAnnealing::PopulationAnnealing(Graph& structure, Config config) {
    if(config.seed != 0) {
        rng_ = RandomNumberGenerator(config.seed);
    }

    schedule_ = config.schedule;
    schedule_.back().compute_observables = true;
    beta_ = NAN;
    coeff_D_ = NAN;
    coeff_P_ = NAN;
    structure_ = structure;
    structure_.Compress();
    average_population_ = config.population;
    init_population_ = average_population_;
 }

double PopulationAnnealing::Hamiltonian(const StateVector& replica) {
    return coeff_P_ * ProblemHamiltonian(replica) + coeff_D_ * DriverHamiltonian(replica);
}

double PopulationAnnealing::DriverHamiltonian(const StateVector& replica) {
    double energy = 0;
    for(std::size_t site = 0; site < structure_.size(); ++site) {
        auto parity = replica[site] ^ Rotate(replica[site], 1);
        energy += EvalFunctor(parity, [&](const auto& v) { return v; }, 0);
    }
    return energy;
}

double PopulationAnnealing::ProblemHamiltonian(const StateVector& replica) {
    double energy = 0;
    for(std::size_t site = 0; site < structure_.size(); ++site) {
        for(std::size_t edge = 0; edge < structure_.adjacent()[site].size(); ++edge) {
            auto parity = replica[site] ^ replica[structure_.adjacent()[site][edge]];
            auto weight = structure_.weights()[site][edge];
            energy += EvalFunctor(parity, [&](const auto& v) { return v*weight; }, static_cast<decltype(energy)>(0));
        }
    }
    energy /= 2;

    for(std::size_t site = 0; site < structure_.size(); ++site) {
        auto field = structure_.fields()[site];
        energy += EvalFunctor(replica[site], [&](const auto& v) { return v*field; }, static_cast<decltype(energy)>(0));
    }
    return energy;
}

// double PopulationAnnealing::SliceProblemHamiltonian(StateVector& replica, int slice) {
//     return (structure_.Adjacent().triangularView<Eigen::Upper>() * GetTrotterSlice(replica, slice).cast<EdgeType>()).dot(GetTrotterSlice(replica, slice).cast<EdgeType>());
// }

// double PopulationAnnealing::DeltaProblemEnergy(StateVector& replica, int vertex) {
//     auto slice = vertex / structure_.size();
//     return -2 * replica(vertex) * (structure_.Adjacent().innerVector(vertex % structure_.size()).dot(GetTrotterSlice(replica, slice).cast<EdgeType>()) - structure_.Fields()(vertex % structure_.size()));
// }

// double PopulationAnnealing::DeltaDriverEnergy(StateVector& replica, int vertex) {
//     auto slice = vertex / structure_.size();
//     auto stride = structure_.size();
//     return -2 * replica(vertex) * ( replica((slice + stride)%replica.size()) + replica(((slice + replica.size()) - stride)%replica.size()) );
// }

// void PopulationAnnealing::WolffSweep(StateVector& replica, int moves) {
//     std::vector<char> cluster(ktrotter_slices);
//     double growth_prob = 1.0 - std::exp(2.0 * beta_ * coeff_D_);
//     for(std::size_t k = 0; k < moves * structure_.size(); ++k) {
//         // Setup
//         std::fill(cluster.begin(), cluster.end(), 1);
//         int seed = rng_.Range(replica.size());
//         int site = seed % structure_.size();
//         int seed_slice = seed / structure_.size();
//         double delta_energy = DeltaProblemEnergy(replica, seed);
//         cluster[seed_slice] = -1;

//         // Grow cluster upwards and downwards
//         for(auto direction : {1, -1}) {
//             int prev_slice = seed_slice;
//             for(int next_slice = (prev_slice + direction + ktrotter_slices)%ktrotter_slices; cluster[next_slice] == 1;
//                 next_slice = (prev_slice + direction + ktrotter_slices)%ktrotter_slices) {

//                 int prev_spin = prev_slice * structure_.size() + site;
//                 int next_spin = next_slice * structure_.size() + site;

//                 if(replica(next_spin) == replica(prev_spin) && rng_.Probability() < growth_prob) {
//                     cluster[next_slice] = -1;
//                     delta_energy += DeltaProblemEnergy(replica, next_spin);
//                 }else {
//                     break;
//                 }
//                 prev_slice = next_slice;
//             }
//         }
        
//         // Attempt to flip
//         delta_energy *= coeff_P_;
//         if(MetropolisAcceptedMove(delta_energy)) {
//             for(std::size_t i = 0; i < cluster.size(); ++i) {
//                 replica(i * structure_.size() + site) *= cluster[i];
//             }
//         }
//     }
// }

void PopulationAnnealing::WolffSweep(StateVector& replica, std::size_t moves) {
    double growth_prob = 1.0 - std::exp(2.0 * beta_ * coeff_D_);
    for(std::size_t k = 0; k < moves * structure_.size(); ++k) {

        std::uint32_t seed = rng_.Range(replica.size());
        std::uint32_t site = seed % structure_.size();
        std::uint32_t seed_slice = seed / structure_.size();
        VertexType cluster = 1 << seed_slice;
        VertexType spins = replica[site];
        // Spins of the same parity evaluate to true
        spins = (spins & cluster) ? spins : !spins;

        // Build Cluster
        // Note: it may be faster to evaluate growth prob all at once
        for(auto direction : {1, -1}) {
            VertexType mask = 1 << seed_slice;
            for(std::size_t i = 0; i < ktrotter_slices; ++i) {
                mask = Rotate(mask, direction);
                if(spins & mask && rng_.Probability() < growth_prob) {
                    cluster |= mask;
                } else {
                    break;
                }
                if(cluster == std::numeric_limits<decltype(cluster)>::max()) {
                    break;
                }
            }
        }

        // Compute spatial energy delta
        std::array<double, ktrotter_slices> site_delta_energy;
        SpatialSiteEnergy(replica, site, site_delta_energy.begin());
        double delta_energy = 0;
        VertexType mask = 1;
        for(std::size_t i = 0; i < ktrotter_slices; ++i) {
            delta_energy += mask & cluster ? site_delta_energy[i] : 0;
        }
        delta_energy += structure_.fields()[site] * GetValue(cluster & replica[site]) * PopCount(spins & cluster);
        delta_energy *= -2 * coeff_P_;

        // Flip cluster
        if(AcceptedMove(-delta_energy*beta_)) {
            replica[site] ^= cluster;
        }
    }
}

std::vector<PopulationAnnealing::Result> PopulationAnnealing::Run() {
    std::vector<Result> results;

    replicas_.resize(average_population_);
    replica_families_.resize(average_population_);
    
    for(auto& r : replicas_) {
        r = StateVector();
        r.resize(structure_.size());
        for(std::size_t k = 0; k < r.size(); ++k) {
            static_assert(sizeof(decltype(r[k])) == 8);
            r[k] = rng_() ^ (static_cast<std::uint64_t>(rng_()) << 32);
        }
    }

    std::iota(replica_families_.begin(), replica_families_.end(), 0);
    SetParams(schedule_.front().beta, schedule_.front().gamma, schedule_.front().lambda);
    std::vector<double> energy;

    auto total_time_start = std::chrono::high_resolution_clock::now();
    unsigned long long int total_sweeps = 0;

    for(auto step : schedule_) {
        Result observables;

        auto time_start = std::chrono::high_resolution_clock::now();
        observables.norm_factor = Resample(step.beta, step.gamma, step.lambda, step.population_fraction);
        for(std::size_t k = 0; k < replicas_.size(); ++k) {
            WolffSweep(replicas_[k], step.wolff_sweeps);
        }
        total_sweeps += replicas_.size() * step.wolff_sweeps;
        
        observables.montecarlo_walltime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - time_start).count();

        observables.beta = beta_;
        observables.gamma = step.gamma;
        observables.lambda = step.lambda;
        observables.d = coeff_D_;
        observables.p = coeff_P_;
        observables.population = replicas_.size();

        if(step.compute_observables) {
            energy.resize(replicas_.size() * ktrotter_slices);
            for(std::size_t k = 0; k < replicas_.size(); ++k) {
                SpatialProblemHamiltonain(replicas_[k], energy.begin() + k * ktrotter_slices);
            }
            // Basic observables
            observables.average_energy = std::accumulate(energy.begin(), energy.end(), 0.0)/energy.size();
            observables.average_squared_energy = std::accumulate(energy.begin(), energy.end(), 0.0, [](const auto& acc, const auto& v) { return acc + v*v; })/energy.size();
            observables.ground_energy = *std::min_element(energy.begin(), energy.end());
            // Round-off /probably/ isn't an issue here
            observables.grounded_replicas = std::accumulate(energy.begin(), energy.end(), 0, [&](const auto& acc, const auto& v) {
                return acc + (util::FuzzyCompare(v, observables.ground_energy) ? 1 : 0);
            });
            // Family statistics
            std::vector<double> family_size = FamilyCount();
            std::transform(family_size.begin(),family_size.end(),family_size.begin(),
                [&](double n) -> double {return n /= observables.population;});
            // Entropy
            observables.entropy = -std::accumulate(family_size.begin(), family_size.end(), 0.0, 
                [](double acc, double n) {return acc + n*std::log(n); });
            // Mean Square Family Size
            observables.mean_square_family_size = observables.population * 
                std::accumulate(family_size.begin(), family_size.end(), 0.0, [](double acc, double n) {return acc + n*n; });
        }

        observables.seed = rng_.GetSeed();
        observables.sweeps = step.wolff_sweeps;
        observables.total_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - total_time_start).count();
        observables.total_sweeps = total_sweeps;
        results.push_back(observables);
    }
    return results;
}

// bool PopulationAnnealing::MetropolisAcceptedMove(double delta_energy) {
//     if(delta_energy < 0.0) {
//         return true;
//     }
    
//     double acceptance_prob_exp = -delta_energy*beta_;
//     return AcceptedMove(acceptance_prob_exp);
// }

bool PopulationAnnealing::AcceptedMove(double log_probability) {
    double test = rng_.Probability();
    // Compute bound on log of test number
    auto bound = log_lookup_(test);

    if(bound.upper < log_probability) {
        return true;
    }else if(bound.lower > log_probability) {
        return false;
    }
    // Compute exp if LUT can't resolve it
    return std::exp(log_probability) > test;
}

double PopulationAnnealing::Resample(double new_beta, double new_gamma, double new_lambda, double new_population_fraction) {
    auto old_driver_coeff = coeff_D_;
    auto old_problem_coeff = coeff_P_;
    auto old_beta = beta_;
    SetParams(new_beta, new_gamma, new_lambda);
    auto new_driver_coeff = coeff_D_;
    auto new_problem_coeff = coeff_P_;

    std::vector<StateVector> resampled_replicas;
    std::vector<int> resampled_families;
    resampled_replicas.reserve(replicas_.size());
    resampled_families.reserve(replicas_.size());

    average_population_ = new_population_fraction * init_population_;
    
    if(average_population_ == 1.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> log_weights(replicas_.size());
    std::transform(replicas_.begin(), replicas_.end(), log_weights.begin(), [&](auto& r) {
        auto driver_energy = this->DriverHamiltonian(r);
        auto problem_energy = this->ProblemHamiltonian(r);
        return -(new_beta*new_problem_coeff-old_beta*old_problem_coeff) * problem_energy - (new_beta*new_driver_coeff-old_beta*old_driver_coeff) * driver_energy;
    });

    auto log_norm = util::LogSumExp(log_weights.begin(), log_weights.end());

    for(std::size_t k = 0; k < replicas_.size(); ++k) {
        double weight = average_population_ * std::exp(log_weights[k] - log_norm);
        unsigned int n = (weight - std::floor(weight)) > rng_.Probability() ? static_cast<unsigned int>(std::ceil(weight)) : static_cast<unsigned int>(std::floor(weight));
        for(std::size_t i = 0; i < n; ++i) {
            resampled_replicas.push_back(replicas_[k]);
            resampled_families.push_back(replica_families_[k]);
        }
    }
    replicas_ = resampled_replicas;
    replica_families_ = resampled_families;
    return std::exp(log_norm);
}

void PopulationAnnealing::SetParams(double new_beta, double new_gamma, double new_lambda) {
    beta_ = new_beta;
    coeff_D_ = std::log(std::tanh(beta_ * new_gamma / ktrotter_slices)) / (2.0 * beta_);
    coeff_P_ = new_lambda/ktrotter_slices;
}
}
