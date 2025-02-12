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
#include "types.hpp"
#include "graph.hpp"
#include "random_number_generator.hpp"
#include "population_annealing_base.hpp"
#include "log_lookup.hpp"
#include <Eigen/Dense>

namespace propane {
/** Implementation of Population Annealing Monte Carlo.
 * Replicas have an associated entry in the family vector indicating lineage.
 */
class PopulationAnnealing : public PopulationAnnealingBase {
protected:
    util::LogLookup log_lookup_;

    Graph structure_;

    RandomNumberGenerator rng_;

    using StateVector = Eigen::Matrix<VertexType, Eigen::Dynamic, 1>;
    std::vector<StateVector> replicas_;
    std::vector<int> replica_families_;

    int init_population_;
    int average_population_;
    std::vector<Schedule> schedule_;
    double beta_;
    bool solver_mode_;

/** Determinstically builds a list of replicas with different Markov Chains.
 */
    std::vector<std::pair<int, int>> BuildReplicaPairs();
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
/** Returns true if a move is accepted according to the Metropolis algorithm.
 */
    virtual bool MetropolisAcceptedMove(double delta_energy);

/** Uses the Heat Bath algorithm. See MetropolisAcceptedMove.
 */
 
    virtual bool HeatbathAcceptedMove(double delta_energy);
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
/** Carries out moves monte carlo sweeps of replica using the Metropolis algorithm.
 */
    virtual void MetropolisSweep(StateVector& replica, int moves);
/** Carries out moves monte carlo sweeps of replica using the Heatbath algorithm.
 */
    virtual void HeatbathSweep(StateVector& replica, int moves);
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