#include "fpga_population_annealing.hpp"
#include "graph.hpp"
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace propane
{

std::vector<FpgaPopulationAnnealing::Result> FpgaPopulationAnnealing::Run() {
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
            SweepPopulation(step.sweeps);
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
>
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

FpgaPopulationAnnealing::FpgaPopulationAnnealing(Graph& structure, Config config) {
    // mmap bars
    cfg_disc_ = open("/sys/bus/pci/devices/0000:81:00.0/resource0", O_RDWR | O_SYNC);
    cfg_bar_ = reinterpret_cast<volatile std::std::uint32_t*>(mmap(nullptr, , PROT_READ | PROT_WRITE, MAP_SHARED, cfg_disc_, 0));
    if(!cfg_bar_) {
        throw std::exception("mmap failed.");
    }
    
    mdl_disc_ = open("/sys/bus/pci/devices/0000:81:00.0/resource0", O_RDWR | O_SYNC);
    mdl_bar_ = reinterpret_cast<volatile std::std::uint32_t*>(mmap(nullptr, , PROT_READ | PROT_WRITE, MAP_SHARED, mdl_disc_, 0));
    if(!mdl_bar_) {
        throw std::exception("mmap failed.");
    }

    mem_disc = open("/sys/bus/pci/devices/0000:81:00.0/resource0", O_RDWR | O_SYNC);
    mem_bar_ = reinterpret_cast<volatile std::std::uint32_t*>(mmap(nullptr, , PROT_READ | PROT_WRITE, MAP_SHARED, mem_disc_, 0));
    if(!mem_bar_) {
        throw std::exception("mmap failed.");
    }
    // load configuration

}

FpgaPopulationAnnealing::SweepPopulation(int sweeps) {
    int base_pointer = 0;

    for(int k = 0; k <= replicas_.size() / 64; k++) {
        std::vector<std::uint64_t> replica_pack(replicas_.front().size());
        // Package replicas
        for(int i = 0; i < replica_pack.size(); ++i) {
            replica_pack[i] = 0;
            // Assemble single package
            for(int j = 0; j < 64; ++j) {
                replica_pack[j] <<= 1;
                if(j + k*64 < replicas.size()) {
                    replica_pack[j] |= replicas_[k*64+j](i) == 1 ? 1 : 0;
                }
            }
        }

        // Copy to accelerator memory
        for(int i = 0; i < replica_pack.size(); ++i) {
            mem_bar_[base_pointer + i + replica_pack.size() * k] = replica_pack[i]; 
        }
    }

// -------

    // Configure and run
    cfg_bar_[kReplicaSetCount] = (replicas_.size()+63) / 64;
    cfg_bar_[kBasePointer] = base_pointer;
    cfg_bar_[kBeta] = beta_;
    cfg_bar_[kSweeps] = sweeps;

    // Wait for completion
    while(cfg_bar_[kSweeps] > 0) { ; }

// -------

    // Unpack
    for(int k = 0; k <= replicas_.size() / 64; k++) {
        std::vector<std::uint64_t> replica_pack(replicas_.front().size());
        
        // Copy froms accelerator memory
        for(int i = 0; i < replica_pack.size(); ++i) {
            replica_pack[i] = mem_bar_[base_pointer + i + replica_pack.size() * k];
        }
        
        // Unpack replicas
        for(int i = replica_pack.size()-1; i >= 0; --i) {
            for(int j = 0; j < 64; ++j) {
                replica_pack[j] >>= 1;
                if(j + k*64 < replicas.size()) {
                    replicas_[k*64+j](i) = replica_pack[j] & 1 ? 1 : -1;
                }
            }
        }
    }
}
}