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

namespace propane {
// Fix this to do something with the single replicas
std::vector<std::pair<int, int>> PopulationAnnealing::BuildReplicaPairs() {
    std::vector<std::pair<int, int>> pairs;
    int num_pairs = replica_families_.size()/2;
    pairs.reserve(num_pairs);

    for(int k = 0; k < num_pairs; ++k) {
        pairs.push_back({k, k + num_pairs});
    }

    // Consistency check

    for(auto& it : pairs) {
        // this is highly unlikely to be true
        if(replica_families_[it.first] == replica_families_[it.second]) {
            for(auto& it_other : pairs) {
                if(replica_families_[it_other.first] != replica_families_[it.second] && 
                    replica_families_[it.first] != replica_families_[it_other.second]) {
                    it.swap(it_other);
                    break;
                }
            }
        }
    }
    return pairs;
}

std::vector<PopulationAnnealing::Result::Histogram> PopulationAnnealing::BuildHistogram(const std::vector<double>& samples) {
    std::vector<Result::Histogram> hist;
    for(auto v : samples) {
        auto it = std::lower_bound(hist.begin(), hist.end(), v, [&](const Result::Histogram& a, const double& b) { return a.bin < b && !util::FuzzyCompare(a.bin, b); });
        if(it==hist.end() || !util::FuzzyCompare(v, it->bin)) {
            hist.insert(it, {v, 1.0});
        } else {
            it->value++;
        }
    }
    for(auto& bin : hist) {
        bin.value /= samples.size();
    }
    return hist;
}

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
    beta_ = NAN;
    coeff_D_ = NAN;
    coeff_P_ = NAN;
    structure_ = structure;
    structure_.Adjacent().makeCompressed();
    average_population_ = config.population;
    init_population_ = average_population_;
    solver_mode_ = config.solver_mode;
    trotter_slices_ = config.trotter_slices;
 }

double PopulationAnnealing::Hamiltonian(StateVector& replica) {
    return coeff_P_ * ProblemHamiltonian(replica) + coeff_D_ * DriverHamiltonian(replica);
}

double PopulationAnnealing::DriverHamiltonian(StateVector& replica) {
    double energy = 0;
    for(int qvertex = 0; qvertex < structure_.size(); ++qvertex) {
        for(int slice = 0; slice < trotter_slices_; ++slice) {
            energy += replica(qvertex + slice * structure_.size()) * replica(qvertex + (slice+1+trotter_slices_)%trotter_slices_ * structure_.size());
        }
    }
    return energy;
}

double PopulationAnnealing::ProblemHamiltonian(StateVector& replica) {
    // There should be a way to construct a view on a full block diagonal adjacency matrix
    double energy = 0;
    for(auto slice = 0; slice < trotter_slices_; ++slice) {
        if(structure_.has_field()) {
            energy += (structure_.Adjacent().triangularView<Eigen::Upper>() * GetTrotterSlice(replica, slice).cast<EdgeType>() - structure_.Fields()).dot(GetTrotterSlice(replica, slice).cast<EdgeType>());
        }else {
            energy += (structure_.Adjacent().triangularView<Eigen::Upper>() * GetTrotterSlice(replica, slice).cast<EdgeType>()).dot(GetTrotterSlice(replica, slice).cast<EdgeType>());
        }
    }
    return energy;
}

double PopulationAnnealing::DeltaProblemEnergy(StateVector& replica, int vertex) {
    auto slice = vertex / structure_.size();
    return -2 * replica(vertex) * (structure_.Adjacent().innerVector(vertex % structure_.size()).dot(GetTrotterSlice(replica, slice).cast<EdgeType>()) - structure_.Fields()(vertex % structure_.size()));
}

double PopulationAnnealing::DeltaDriverEnergy(StateVector& replica, int vertex) {
    auto slice = vertex / structure_.size();
    auto stride = structure_.size();
    return -2 * replica(vertex) * ( replica((slice + stride)%replica.size()) + replica(((slice + replica.size()) - stride)%replica.size()) );
}

void PopulationAnnealing::MetropolisSweep(StateVector& replica, int sweeps) {
    for(std::size_t k = 0; k < sweeps; ++k) {
        for(std::size_t i = 0; i < replica.size(); ++i) {
            int vertex = i;
            double delta_problem_energy = DeltaProblemEnergy(replica, vertex);
            double delta_driver_energy = DeltaDriverEnergy(replica, vertex);
            double delta_energy = coeff_P_ * delta_problem_energy + coeff_D_ * delta_driver_energy;

            //round-off isn't a concern here
            if(MetropolisAcceptedMove(delta_energy)) {
                replica(vertex) *= -1;
            }
        }
    }
}

void PopulationAnnealing::HeatbathSweep(StateVector& replica, int sweeps) {
    for(std::size_t k = 0; k < sweeps; ++k) {
        for(std::size_t vertex = 0; vertex < replica.size(); ++vertex) {
            double delta_problem_energy = DeltaProblemEnergy(replica, vertex);
            double delta_driver_energy = DeltaDriverEnergy(replica, vertex);
            double delta_energy = coeff_P_ * delta_problem_energy + coeff_D_ * delta_driver_energy;

            if(HeatbathAcceptedMove(delta_energy)) {
                replica(vertex) *= -1;
            }
        }
    }
}

bool PopulationAnnealing::IsLocalMinimum(StateVector& replica) {
    for(int k = 0; k < replica.size(); ++k) {
        if(DeltaProblemEnergy(replica, k) < 0) {
            return false;
        }
    }
    return true;
}

double PopulationAnnealing::Overlap(StateVector& alpha, StateVector& beta) {
    return (alpha.array() * beta.array()).cast<double>().sum() / structure_.size();
}

double PopulationAnnealing::LinkOverlap(StateVector& alpha, StateVector& beta) {
    double ql = 0;
    for(std::size_t k = 0; k < structure_.Adjacent().outerSize(); ++k) {
        for(Eigen::SparseTriangularView<Eigen::SparseMatrix<EdgeType>,Eigen::Upper>::InnerIterator 
            it(structure_.Adjacent().triangularView<Eigen::Upper>(), k); it; ++it) {
            ql += alpha(k) * beta(k) * alpha(it.index()) * beta(it.index());
        }
    }
    return ql / structure_.edges();
}

std::vector<PopulationAnnealing::Result> PopulationAnnealing::Run() {
    std::vector<Result> results;

    replicas_.resize(average_population_);
    replica_families_.resize(average_population_);
    
    for(auto& r : replicas_) {
        r = StateVector();
        r.resize(structure_.size()*trotter_slices_);
        for(std::size_t k = 0; k < r.size(); ++k) {
            r(k) = rng_.Probability() < 0.5 ? 1 : -1;
        }
    }

    std::iota(replica_families_.begin(), replica_families_.end(), 0);
    SetParams(schedule_.front().gamma, schedule_.front().beta);
    std::vector<double> energy;

    auto total_time_start = std::chrono::high_resolution_clock::now();
    unsigned long long int total_sweeps = 0;

    for(auto step : schedule_) {
        Result observables;

        auto time_start = std::chrono::high_resolution_clock::now();
        observables.norm_factor = Resample(step.beta, step.gamma, step.population_fraction);
        for(std::size_t k = 0; k < replicas_.size(); ++k) {
            if(step.heat_bath) {
                HeatbathSweep(replicas_[k], step.sweeps);
            }else {
                MetropolisSweep(replicas_[k], step.sweeps);
            }
        }
        total_sweeps += replicas_.size() * step.sweeps;
        
        observables.montecarlo_walltime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - time_start).count();

        observables.beta = beta_;
        observables.gamma = step.gamma;
        observables.d = coeff_D_;
        observables.p = coeff_P_;
        observables.population = replicas_.size();

        if(!solver_mode_ || beta_ == schedule_.back().beta) {
            if(step.compute_observables) {
                energy.resize(replicas_.size());
                for(std::size_t k = 0; k < replicas_.size(); ++k) {
                    energy[k] = ProblemHamiltonian(replicas_[k]);
                }
                Eigen::Map<Eigen::VectorXd> energy_map(energy.data(), energy.size());
                // Basic observables
                observables.average_energy = energy_map.mean();
                observables.average_squared_energy = energy_map.array().pow(2).mean();
                observables.ground_energy = energy_map.minCoeff();
                // Round-off /probably/ isn't an issue here
                observables.grounded_replicas = energy_map.array().unaryExpr(
                    [&](double E){return util::FuzzyCompare(E, observables.ground_energy) ? 1 : 0;}).sum();
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

                if(step.energy_dist) {
                    // Energy
                    observables.energy_distribution = BuildHistogram(energy);
                }
                
                if(step.ground_dist) {
                    std::vector<StateVector> ground_states;
                    std::vector<double> samples;
                    samples.reserve(observables.grounded_replicas);
                    for(std::size_t k = 0; k < replicas_.size(); ++k) {
                        if(util::FuzzyCompare(observables.ground_energy, energy[k])) {
                            auto state = std::distance(ground_states.begin(), std::find(ground_states.begin(), ground_states.end(), replicas_[k]));
                            if(state == ground_states.size()) {
                                ground_states.push_back(replicas_[k]);
                            }
                            
                            samples.push_back(state);
                        }
                    }
                    observables.ground_distribution = BuildHistogram(samples);
                }
            }

            if(step.overlap_dist) {
                // Overlap
                std::vector<std::pair<int, int>> overlap_pairs = BuildReplicaPairs();
                std::vector<double> overlap_samples(overlap_pairs.size());

                std::transform(overlap_pairs.begin(), overlap_pairs.end(), overlap_samples.begin(),
                    [&](std::pair<int, int> p){return Overlap(replicas_[p.first], replicas_[p.second]);});
                observables.overlap = BuildHistogram(overlap_samples);
                // Link Overlap
                std::transform(overlap_pairs.begin(), overlap_pairs.end(), overlap_samples.begin(),
                    [&](std::pair<int, int> p){return LinkOverlap(replicas_[p.first], replicas_[p.second]);});
                observables.link_overlap = BuildHistogram(overlap_samples);
            }

            observables.seed = rng_.GetSeed();
            observables.sweeps = step.sweeps;
            observables.total_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - total_time_start).count();
            observables.total_sweeps = total_sweeps;
            results.push_back(observables);
        }
    }
    return results;
}

bool PopulationAnnealing::HeatbathAcceptedMove(double delta_energy) {
    double accept_prob = 1.0/(1.0 + std::exp(-delta_energy*beta_));
    return rng_.Probability() < accept_prob;
}

bool PopulationAnnealing::MetropolisAcceptedMove(double delta_energy) {
    if(delta_energy < 0.0) {
        return true;
    }
    
    double acceptance_prob_exp = -delta_energy*beta_;
    return AcceptedMove(acceptance_prob_exp);
}

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

double PopulationAnnealing::Resample(double new_beta, double new_gamma, double new_population_fraction) {
    auto old_driver_coeff = coeff_D_;
    auto old_problem_coeff = coeff_P_;
    auto old_beta = beta_;
    SetParams(new_gamma, new_beta);
    auto new_driver_coeff = coeff_D_;
    auto new_problem_coeff = coeff_P_;

    std::vector<StateVector> resampled_replicas;
    std::vector<int> resampled_families;
    resampled_replicas.reserve(replicas_.size());
    resampled_families.reserve(replicas_.size());

    average_population_ = new_population_fraction * init_population_;
    
    std::vector<double> log_weights(replicas_.size());
    std::transform(replicas_.begin(), replicas_.end(), log_weights.begin(), [&](auto& r) {
        auto driver_energy = DriverHamiltonian(r);
        auto problem_energy = ProblemHamiltonian(r);
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

void PopulationAnnealing::SetParams(double new_gamma, double new_beta) {
    beta_ = new_beta;
    coeff_D_ = std::log(std::tanh(beta_ * new_gamma / trotter_slices_)) / (2.0 * beta_);
    coeff_P_ = 1.0/trotter_slices_;
}
}
