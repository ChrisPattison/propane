#pragma once
#include <vector>
#include <utility>
#include "types.hpp"
#include "graph.hpp"
#include "types.hpp"
#include "random_number_generator.hpp"
#include "Eigen/Dense"

class PopulationAnnealing {
public:
    struct Result {
        struct Histogram {
            double bin;
            double value;
        };
        std::vector<Histogram> overlap;
        std::vector<Histogram> link_overlap;
        double beta;
        std::size_t population;
        double average_energy;
        double average_energy_squared;
        double ground_energy;
        std::size_t grounded_replicas;
        double entropy;
    };

    
protected:
    Graph structure_;

    RandomNumberGenerator rng_;

    using StateVector = Eigen::Matrix<VertexType, Eigen::Dynamic, 1>;
    std::vector<StateVector> replicas_;
    std::vector<int> replica_families_;

    int average_population_;
    std::vector<double> betalist_;
    double beta_;

    // Builds a list of replicas with different Markov Chains
    std::vector<std::pair<int, int>> BuildReplicaPairs();
    // Gets the number of replicas in each family as a fraction of the total population
    std::vector<double> FamilyCount();

    double AcceptanceProbability(double delta_energy) const;
    // Gives the Hamiltonian of the given state
    double Hamiltonian(StateVector& replica);
    // Change in energy associated with a single change in state
    double DeltaEnergy(StateVector& replica, int vertex); 

    void MonteCarloSweep(StateVector& replica, int moves);

    bool IsLocalMinimum(StateVector& replica);

    StateVector Quench(const StateVector& replica);

    double Overlap(StateVector& alpha, StateVector& beta);

    double LinkOverlap(StateVector& alpha, StateVector& beta);

    std::vector<Result::Histogram> BuildHistogram(const std::vector<double>& samples);

    void Resample(double new_beta);
public:

    PopulationAnnealing() = delete;

    PopulationAnnealing(Graph& structure, std::vector<double> betalist, int average_population);

    void Run(std::vector<Result>& results);
};