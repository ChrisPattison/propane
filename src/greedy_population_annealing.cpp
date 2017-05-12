#include "greedy_population_annealing.hpp"

namespace propane
{
double GreedyPopulationAnnealing::DeltaEnergy(StateVector& replica, int vertex) {
    replica(vertex) *= -1;
    double new_energy = Hamiltonian(replica);

    replica(vertex) *= -1;
    double old_energy = Hamiltonian(replica);
    
    return new_energy - old_energy;
}

double GreedyPopulationAnnealing::Hamiltonian(StateVector& replica) {
    StateVector reduced_replica = replica;
    for(;;) {
        Eigen::Matrix<EdgeType, Eigen::Dynamic, 1>::Index index;
        auto delta_energy = (-(structure_.Adjacent() * reduced_replica.cast<EdgeType>()  - structure_.Fields()).array() * reduced_replica.cast<EdgeType>().array()).minCoeff(&index);
        if(delta_energy >= 0) {
            break;
        } else {
            reduced_replica(index) *= -1;
        }
    }
    return PopulationAnnealing::Hamiltonian(reduced_replica);
}
}