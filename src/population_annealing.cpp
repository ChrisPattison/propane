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

namespace propane
{
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
 }

double PopulationAnnealing::Hamiltonian(StateVector& replica) {
    if(structure_.has_field()) {
        return (structure_.Adjacent().triangularView<Eigen::Upper>() * replica.cast<EdgeType>() - structure_.Fields()).dot(replica.cast<EdgeType>());
    }else {
        return (structure_.Adjacent().triangularView<Eigen::Upper>() * replica.cast<EdgeType>()).dot(replica.cast<EdgeType>());
    }
}

double PopulationAnnealing::DeltaEnergy(StateVector& replica, int vertex) {
    return -2 * replica(vertex) * (structure_.Adjacent().innerVector(vertex).dot(replica.cast<EdgeType>()) - structure_.Fields()(vertex));
}

void PopulationAnnealing::MonteCarloSweep(StateVector& replica, int sweeps) {
    for(std::size_t k = 0; k < sweeps; ++k) {
        for(std::size_t i = 0; i < replica.size(); ++i) {
            int vertex = i;
            double delta_energy = DeltaEnergy(replica, vertex);
            
            //round-off isn't a concern here
            if(AcceptedMove(delta_energy)) {
                replica(vertex) *= -1;
            }
        }
    }
}

bool PopulationAnnealing::IsLocalMinimum(StateVector& replica) {
    for(int k = 0; k < replica.size(); ++k) {
        if(DeltaEnergy(replica, k) < 0) {
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
        r.resize(structure_.size());
        for(std::size_t k = 0; k < r.size(); ++k) {
            r(k) = rng_.Probability() < 0.5 ? 1 : -1;
        }
    }

    std::iota(replica_families_.begin(), replica_families_.end(), 0);
    beta_ = schedule_.front().beta;
    std::vector<double> energy;

    for(auto step : schedule_) {
        Result observables;

        auto time_start = std::chrono::high_resolution_clock::now();
        if(step.beta != beta_) {
            observables.norm_factor = Resample(step.beta, step.population_fraction);
        }
        for(std::size_t k = 0; k < replicas_.size(); ++k) {
            MonteCarloSweep(replicas_[k], step.sweeps);
        }
        observables.montecarlo_walltime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - time_start).count();

        energy.resize(replicas_.size());
        for(std::size_t k = 0; k < replicas_.size(); ++k) {
            energy[k] = Hamiltonian(replicas_[k]);
        }
        Eigen::Map<Eigen::VectorXd> energy_map(energy.data(), energy.size());
        // Basic observables
        observables.beta = beta_;
        observables.population = replicas_.size();
        observables.average_energy = energy_map.mean();
        observables.average_squared_energy = energy_map.array().pow(2).mean();
        observables.ground_energy = energy_map.minCoeff();
        // Round-off /probably/ isn't an issue here
        observables.grounded_replicas = energy_map.array().unaryExpr(
            [&](double E){return E == observables.ground_energy ? 1 : 0;}).sum();
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
        results.push_back(observables);
    }
    return results;
}

bool PopulationAnnealing::AcceptedMove(double delta_energy) {
    if(delta_energy < 0.0) {
        return true;
    }
    // Get probability exponent and test probability
    double acceptance_prob_exp = -delta_energy*beta_;
    double test = rng_.Probability();

    // Compute bound on log of test number
    auto bound = log_lookup_(test);

    // return acceptance_prob_exp > log(test);
    if(bound.upper < acceptance_prob_exp) {
        return true;
    }else if(bound.lower > acceptance_prob_exp) {
        return false;
    }
    // Compute exp if LUT can't resolve it
    return std::exp(acceptance_prob_exp) > test;
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
}