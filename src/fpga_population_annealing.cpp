#include "fpga_population_annealing.hpp"
#include "graph.hpp"
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace propane
{
FpgaPopulationAnnealing::Resample(double new_beta, double population_fraction) {
    // TODO: add in support for padding on the end
    
    std::vector<double> energy(population_.size() * 64);
    std::fill(energy.begin(), energy.end(), 0.0);

    // Compute energy of each replica
    for(int k = 0; k < structure_.Adjacent().outerSize(); ++k) {
        for(Eigen::SparseTriangularView<Eigen::SparseMatrix<EdgeType>,Eigen::Upper>::InnerIterator 
            it(structure_.Adjacent().triangularView<Eigen::Upper>(), k); it; ++it) {
            for(i = 0; i < population_.begin().size(); ++i) {
                for(int j = 0; j < 64; ++j) {
                    index = i*64 + j;
                    energy[index] += 
                        ((population_[k][i] >> j) | 0x01 ? 1 : -1)
                        * ((population_[it.index()][i] >> j) | 0x01 ? 1 : -1)
                        * it.value();
                }
            }
        }
    }

    // Compute resampling weights (suitability)
    std::vector<double> weights(energy.size());
    std::transform(energy.begin(), energy.end(), weights.begin(), [&](double energy) 
        {return std::exp(-(new_beta-beta_) * energy);});

    double summed_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
    double normalize = average_population_ / summed_weights;

    // Determine number of copies of each replica
    std::vector<int> num_replicas(energy.size());
    std::transform(weights.begin(), weights.end(), num_replicas.begin(), [&](double weight) 
        {return weight - std::floor(weight)) > rng_.Probability() ? std::ceil(weight) : std::floor(weight);});

    // Copy spins to new vector, copy families to new vector
    int new_size = std::accumulate(num_replicas.begin(), num_replicas.end(), 0);
    // Number of uint64_t required to store packed spins
    int packed_size = new_size/64 + (new_size % 64 != 0);
    // Make allocations
    StateVectorPack new_population(population_.size());
    std::for_each(new_population.begin(), new_population.end(), [&](auto& v) { v.resize(packed_size); });
    std::vector<int> new_families(new_size);

    // Make copies
    for(i = 0; i < new_population.size(); ++i) {
        for(j = 0; j < packed_size; ++j) {
            std::uint64_t new_pack
            for(k = 0; k < 64; ++k) {

            }
        }
    }
}
}