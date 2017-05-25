#pragma once
#include "population_annealing.hpp"
#include "types.hpp"
#include "graph.hpp"
#include "monte_carlo_driver.hpp"
#include <Eigen/Dense>

namespace propane
{
/** Accelerated Implementation of Population Annealing Monte Carlo.
 * This implementation stores spin data very differently than the others due to the need to offload to the FPGA
 * population_ has dimensions num(spins) x num(replicas/64)
 */
class FpgaPopulationAnnealing : public PopulationAnnealing {
private:
    const int kVendorID = 0x7075;
    const int kDeviceID = 0x1337;

    const std::size_t kMemSize = 0x200'000;

    const int kPackSize = 256;

    volatile std::uint32_t* cfg_bar_;
    volatile std::uint32_t* mdl_bar_;
    volatile std::uint64_t* mem_bar_;

    const int kReplicaSetCount = 0;
    const int kSpinCount = 1;
    const int kBasePointer = 2;
    const int kSeed = 3;
    const int kSweeps = 4;
    const int kBeta = 5;
protected:
    

    Graph structure_;
    RandomNumberGenerator rng_;

    StateVectorPack population_;

    std::vector<Schedule> schedule_;

/** Carries out moves monte carlo sweeps of all replicas on the accelerator.
 */
    void SweepPopulation(int sweeps);
/** Returns true if a move may be made that reduces the total energy.
 */
    // double Resample(double new_beta, double population_fraction);
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
    std::vector<Result> Run();
};
}