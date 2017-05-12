#pragma once
#include <vector>
#include <utility>
#include "types.hpp"
#include "graph.hpp"
#include "random_number_generator.hpp"
#include "population_annealing_base.hpp"
#include "log_lookup.hpp"
#include <Eigen/Dense>

namespace propane 
{    
/** Implementation of Population Annealing Monte Carlo.
 * Replicas have an associated entry in the family vector indicating lineage.
 */
class PopulationAnnealing : public PopulationAnnealingBase {
protected:
    util::LogLookup log_lookup_;

protected:
    Graph structure_;

    RandomNumberGenerator rng_;

    using StateVector = Eigen::Matrix<VertexType, Eigen::Dynamic, 1>;
    std::vector<StateVector> replicas_;
    std::vector<int> replica_families_;

    int init_population_;
    int average_population_;
    std::vector<Schedule> schedule_;
    double beta_;

/** Determinstically builds a list of replicas with different Markov Chains.
 */
    std::vector<std::pair<int, int>> BuildReplicaPairs();
/** Gets the number of replicas in each family as a fraction of the total population
 */
    std::vector<double> FamilyCount();
/** Returns true if a move is accepted according to detailed balance.
 * Uses a look up table to compute a bound on the logarithm of a random number 
 * and compares to the exponent of the acceptance probability.
 * If the probability is inside the bound given by the look table, 
 * true exponetnial is computed and compared.
 */
    virtual bool AcceptedMove(double delta_energy);
/** Returns the energy of a replica
 * Implemented as the sum of elementwise multiplication of the replica vector with the 
 * product of matrix multiplication between the upper half of the adjacency matrix
 * and the replica.
 */
    virtual double Hamiltonian(StateVector& replica);
/** Returns the energy change associated with flipping spin vertex.
 * Implemented as the dot product of row vertex of the adjacency matrix 
 * with the replica vector multiplied by the spin at vertex.
 */
    virtual double DeltaEnergy(StateVector& replica, int vertex); 
/** Carries out moves monte carlo sweeps of replica.
 */
    virtual void MonteCarloSweep(StateVector& replica, int moves);
/** Returns true if a move may be made that reduces the total energy.
 */
    bool IsLocalMinimum(StateVector& replica);
/** Returns the overlap between replicas alpha and beta.
 */
    double Overlap(StateVector& alpha, StateVector& beta);
/** Returns the link overlap between replicas alpha and beta.
 */
    double LinkOverlap(StateVector& alpha, StateVector& beta);
/** Given computes a histogram of the samples.
 * Does not attempt to find the bins that are zero.
 * Normalizes the values so that the sum of the values is 1.
 */
    std::vector<Result::Histogram> BuildHistogram(const std::vector<double>& samples);
/** Resamples population according to the Boltzmann distribution.
 * Attempts to maintain approximately the same population as detailed in arXiv:1508.05647
 * Returns the normalization factor Q as a byproduct.
 */
    virtual double Resample(double new_beta, double new_population_fraction);
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
};
}