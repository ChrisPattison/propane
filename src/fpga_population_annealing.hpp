#pragma once
#include "population_annealing.hpp"
#include "types.hpp"
#include "graph.hpp"
#include "monte_carlo_driver.hpp"
#include <Eigen/Dense>

/** Accelerated Implementation of Population Annealing Monte Carlo.
 * This implementation stores spin data very differently than the others due to the need to offload to the FPGA
 * population_ has dimensions num(spins) x num(replicas/64)
 */
class FpgaPopulationAnnealing : public PopulationAnnealingBase {
    MonteCarloDriver driver_;
public:
protected:
    Graph structure_;
    RandomNumberGenerator rng_;

    using SpinPack = std::vector<std::uint64_t>;
    using StateVectorPack = std::vector<SpinPack>;

    StateVectorPack population_;

    std::vector<Schedule> schedule_;

/** Carries out moves monte carlo sweeps of all replicas on the accelerator.
 */
    void Sweep(int moves);
/** Returns true if a move may be made that reduces the total energy.
 */
    double Resample(double new_beta);
/** Returns new population size
 * Uses a logistic curve with parameters given in input file
 * This probably will be removed in the future with preference given 
 * to the current method of specifying beta schedules (one per temperature)
 */
public:

    FpgaPopulationAnnealing() = delete;
/** Intializes solver.
 * The inputs will be replaced by a struct in the future.
 * schedule specifies the annealing schedule, sweep counts, and histogram generation at each step.
 * seed may be zero in which case one will be generated.
 */
    FpgaPopulationAnnealing(Graph& structure, Config schedule);
/** Run solver and return results.
 */
};