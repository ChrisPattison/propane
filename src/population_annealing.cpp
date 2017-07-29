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
    structure_ = structure;
    structure_.Adjacent().makeCompressed();
    average_population_ = config.population;
    init_population_ = average_population_;
    solver_mode_ = config.solver_mode;
    uniform_init_ = config.uniform_init;
 }

double PopulationAnnealing::Hamiltonian(const StateVector& replica) {
    double energy = 0.0;
    for(std::size_t k = 0; k < structure_.Adjacent().outerSize(); ++k) {
        for(Eigen::SparseTriangularView<Eigen::SparseMatrix<EdgeType>,Eigen::Upper>::InnerIterator 
            it(structure_.Adjacent().triangularView<Eigen::Upper>(), k); it; ++it) {
            energy += replica[k] * it.value() * replica[it.index()];
        }
        energy -= replica[k] * field_[k];
    }
    return energy;
}

PopulationAnnealing::StateVector PopulationAnnealing::Project(const StateVector& replica) {
    StateVector projected;
    projected.resize(replica.size());
    for( std::size_t k = 0; k < replica.size(); ++k ) {
        projected[k] = replica[k] * FieldType(1.0, 0.0) > 0 ? FieldType(1.0, 0.0) : FieldType(-1.0, 0.0);
    }
    return projected;
}

double PopulationAnnealing::ProjectedHamiltonian(const StateVector& projected) {
    double energy = 0.0;
    for(std::size_t k = 0; k < structure_.Adjacent().outerSize(); ++k) {
        for(Eigen::SparseTriangularView<Eigen::SparseMatrix<EdgeType>,Eigen::Upper>::InnerIterator 
            it(structure_.Adjacent().triangularView<Eigen::Upper>(), k); it; ++it) {
            energy += projected[k] * it.value() * projected[it.index()];
        }
        energy -= projected[k] * FieldType(field_[k][0], 0.);
    }
    return energy;
}

FieldType PopulationAnnealing::LocalField(StateVector& replica, int vertex) {
    FieldType h;
    for(Eigen::SparseMatrix<EdgeType>::InnerIterator it(structure_.Adjacent(), vertex); it; ++it) {
        h += it.value() * replica[it.index()];
    }
    h -= field_[vertex];
    return h;
}

double PopulationAnnealing::DeltaEnergy(StateVector& replica, int vertex, FieldType new_value) {
    return (new_value - replica[vertex]) * LocalField(replica, vertex);
}

void PopulationAnnealing::MicroCanonicalSweep(StateVector& replica, int sweeps) {
    for(std::size_t k = 0; k < sweeps; ++k) {
        for(std::size_t i = 0; i < replica.size(); ++i) {
            // pick random site
            auto vertex = rng_.CheapRange(replica.size());
            // get local field
            auto h = LocalField(replica, vertex);
            // normalize
            h /= std::sqrt(h*h);
            // reflect about h
            auto prev = replica[i];
            replica[i] = FieldType(prev * FieldType(h[0]*h[0] - h[1]*h[1], 2*h[0]*h[1]), prev * FieldType(2*h[0]*h[1], h[1]*h[1]-h[0]*h[0]));
        }
    }
}

void PopulationAnnealing::MetropolisSweep(StateVector& replica, int sweeps) {
    for(std::size_t k = 0; k < sweeps; ++k) {
        for(std::size_t i = 0; i < replica.size(); ++i) {
            int vertex = rng_.Range(replica.size());
            auto new_value = VertexType(rng_.Probability());
            double delta_energy = DeltaEnergy(replica, vertex, new_value);
            
            //round-off isn't a concern here
            if(MetropolisAcceptedMove(delta_energy)) {
                replica[vertex] = new_value;
            }
        }
    }
}

void PopulationAnnealing::HeatbathSweep(StateVector& replica, int sweeps) {
    for(std::size_t k = 0; k < sweeps; ++k) {
        for(std::size_t vertex = 0; vertex < replica.size(); ++vertex) {
            auto new_value = VertexType(rng_.Probability());            
            double delta_energy = DeltaEnergy(replica, vertex, new_value);

            if(HeatbathAcceptedMove(delta_energy)) {
                replica[vertex] = new_value;
            }
        }
    }
}

double PopulationAnnealing::Overlap(StateVector& alpha, StateVector& beta) {
    return std::inner_product(alpha.begin(), alpha.end(), beta.begin(), 0.0) / structure_.size();
}

double PopulationAnnealing::LinkOverlap(StateVector& alpha, StateVector& beta) {
    double ql = 0;
    for(std::size_t k = 0; k < structure_.Adjacent().outerSize(); ++k) {
        for(Eigen::SparseTriangularView<Eigen::SparseMatrix<EdgeType>,Eigen::Upper>::InnerIterator 
            it(structure_.Adjacent().triangularView<Eigen::Upper>(), k); it; ++it) {
            ql += alpha[k] * beta[k] * alpha[it.index()] * beta[it.index()];
        }
    }
    return ql / structure_.edges();
}

std::vector<PopulationAnnealing::Result> PopulationAnnealing::Run() {
    // Initialize
    std::vector<Result> results;

    replicas_.resize(average_population_);
    replica_families_.resize(average_population_);
    
    field_.resize(structure_.Fields().size());
    for(int k = 0; k < structure_.Fields().size(); ++k) {
        // No transverse field for now
        field_[k] = FieldType(structure_.Fields()(k), 0.);
    }

    for(auto& r : replicas_) {
        r = StateVector();
        r.resize(structure_.size());
        for(std::size_t k = 0; k < r.size(); ++k) {
            if(uniform_init_) {
                r[k] = FieldType(0.,1.);
            }else {
                // Aligned with sigma_z
                r[k] = FieldType(1.,0.) * (rng_.Probability() < 0.5 ? 1 : -1);
            }
        }
    }

    std::iota(replica_families_.begin(), replica_families_.end(), 0);
    beta_ = schedule_.front().beta;
    TransverseField(schedule_.front().gamma);
    std::vector<double> energy;

    auto total_time_start = std::chrono::high_resolution_clock::now();
    unsigned long long int total_sweeps = 0;

    for(auto step : schedule_) {
        Result observables;
        // Thermalize
        auto time_start = std::chrono::high_resolution_clock::now();
        if(step.beta != beta_ || step.gamma != gamma_) {
            // Set new field
            TransverseField(step.gamma);
            // Resample
            if(step.resample) {
                observables.norm_factor = Resample(step.beta, step.population_fraction);
            }
        }
        for(std::size_t k = 0; k < replicas_.size(); ++k) {
            MicroCanonicalSweep(replicas_[k], step.microcanonical);
            // Monte Carlo Sweeps
            if(step.heat_bath) {
                HeatbathSweep(replicas_[k], step.sweeps);
            }else {
                MetropolisSweep(replicas_[k], step.sweeps);
            }
        }
        total_sweeps += replicas_.size() * step.sweeps;
        
        observables.montecarlo_walltime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - time_start).count();

        observables.beta = beta_;
        observables.gamma = gamma_;
        observables.population = replicas_.size();
        bool report_results = (beta_ == schedule_.back().beta && gamma_ == schedule_.back().gamma);
        if(!solver_mode_ || report_results || step.compute_observables) {
            if(step.compute_observables) {
                std::vector<StateVector> projected_replicas;
                projected_replicas.reserve(replicas_.size());
                std::for_each(replicas_.begin(), replicas_.end(), [&](const StateVector& v) { projected_replicas.push_back(Project(v)); });

                energy.resize(projected_replicas.size());
                for(std::size_t k = 0; k < replicas_.size(); ++k) {
                    energy[k] = ProjectedHamiltonian(projected_replicas[k]);
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
                            auto state = std::distance(ground_states.begin(), std::find(ground_states.begin(), ground_states.end(), projected_replicas[k]));
                            if(state == ground_states.size()) {
                                ground_states.push_back(projected_replicas[k]);
                            }
                            
                            samples.push_back(state);
                        }
                    }
                    observables.ground_distribution = BuildHistogram(samples);
                }

                if(step.overlap_dist) {
                    // Overlap
                    std::vector<std::pair<int, int>> overlap_pairs = BuildReplicaPairs();
                    std::vector<double> overlap_samples(overlap_pairs.size());

                    std::transform(overlap_pairs.begin(), overlap_pairs.end(), overlap_samples.begin(),
                        [&](std::pair<int, int> p){return Overlap(projected_replicas[p.first], projected_replicas[p.second]);});
                    observables.overlap = BuildHistogram(overlap_samples);
                    // Link Overlap
                    std::transform(overlap_pairs.begin(), overlap_pairs.end(), overlap_samples.begin(),
                        [&](std::pair<int, int> p){return LinkOverlap(projected_replicas[p.first], projected_replicas[p.second]);});
                    observables.link_overlap = BuildHistogram(overlap_samples);
                }
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

double PopulationAnnealing::Resample(double new_beta, double new_population_fraction) {
    std::vector<StateVector> resampled_replicas;
    std::vector<int> resampled_families;
    resampled_replicas.reserve(replicas_.size());
    resampled_families.reserve(replicas_.size());

    average_population_ = new_population_fraction * init_population_;

    Eigen::VectorXd weighting(replicas_.size());
    for(std::size_t k = 0; k < replicas_.size(); ++k) {
        weighting(k) = std::exp(-(new_beta-beta_) * Hamiltonian(replicas_[k]));
    }
    
    double summed_weights = weighting.sum();
    double normalize = average_population_ / summed_weights;
    for(std::size_t k = 0; k < replicas_.size(); ++k) {
        double weight = normalize * weighting(k);
        unsigned int n = (weight - std::floor(weight)) > rng_.Probability() ? std::ceil(weight) : std::floor(weight);
        for(std::size_t i = 0; i < n; ++i) {
            resampled_replicas.push_back(replicas_[k]);
            resampled_families.push_back(replica_families_[k]);
        }
    }

    beta_ = new_beta;
    replicas_ = resampled_replicas;
    replica_families_ = resampled_families;
    return summed_weights;
}

void PopulationAnnealing::TransverseField(double magnitude) {
    gamma_ = magnitude;
    for(int k = 0; k < field_.size(); ++k) {
        field_[k][1] = magnitude;
    }
}
}
