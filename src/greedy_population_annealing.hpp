#pragma once
#include "population_annealing.hpp"

namespace propane
{
/** Population Annealing using a redefined Hamiltonian based on a greedy descent to the minimum energy
 */
class GreedyPopulationAnnealing : public PopulationAnnealing {
    using PopulationAnnealing::PopulationAnnealing;
/** Returns the energy of a replica
 * Defined as the result of a greedy search for the ground state
 */
    virtual double Hamiltonian(StateVector& replica);
/** Returns the energy change associated with flipping spin vertex.
 * Naively calls Hamiltonian twice for now.
 */
    virtual double DeltaEnergy(StateVector& replica, int vertex);
};
}