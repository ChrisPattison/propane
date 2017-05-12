#pragma once
#include <vector>
#include <utility>
#include "types.hpp"
#include "graph.hpp"
#include "types.hpp"
#include "random_number_generator.hpp"
#include "population_annealing.hpp"
#include "parallel.hpp"
#include <Eigen/Dense>

namespace propane
{
/** Parallel version of PopulationAnnealing using MPI.
 */
class ParallelPopulationAnnealing : protected PopulationAnnealing {
protected:
/** Minimum exceeded fraction of average node population before redistribution will take place.
  */
    double kMaxPopulation = 1.10;

    parallel::Mpi parallel_;

    int average_node_population_;
/** Resamples population according to the Boltzmann distribution.
 * Attempts to maintain approximately the same population as detailed in arXiv:1508.05647
 * Returns the normalization factor Q as a byproduct.
 */
    double Resample(double new_beta, double population_fraction);
/** Redistributes population evenly among all processes.
 * Ordering of the replicas must remain preserved at all times.
 * Packs replicas in Packets to send to other processes.
 */
    void Redistribute();
/** Combines source with target additively.
 */
    void CombineHistogram(std::vector<Result::Histogram>& target, const std::vector<Result::Histogram>& source);
/** Vector of the integer population of all surviving families originating from current process.
 * Uses Packets to send the count of families to the originator who compiles full count.
 * Invalid if families are store non-contiguously.
 */
    std::vector<double> FamilyCount();
/** Packets replicas into a single vector
 */
    std::vector<VertexType> Pack(const std::vector<StateVector>& source);
/** Packets replicas into a single vector
 */
    std::vector<VertexType> Pack(const std::vector<StateVector>::const_iterator begin_iterator, const std::vector<StateVector>::const_iterator end_iterator);
/** Unpacks a single vector into a vector of replicas
 */
    std::vector<StateVector> Unpack(std::vector<VertexType>& source);
public:
/** Data relavent only to the parallel implementation
 */
    struct Result : PopulationAnnealing::Result {
        double redist_walltime;
        double observables_walltime;
        int max_family_size;
    };
/** Config parameters relavent to only the parallel implementation
 */
    struct Config : PopulationAnnealing::Config {
        double max_population = 1.10;
    };
    ParallelPopulationAnnealing(const ParallelPopulationAnnealing&) = delete;
/** Intializes solver.
 * schedule specifies the annealing schedule, sweep counts, and histogram generation at each step.
 * seed may be zero in which case one will be generated.
 */
    ParallelPopulationAnnealing(Graph& structure, Config config);
/** Run solver and return results.
 */
    std::vector<Result> Run();
};
}