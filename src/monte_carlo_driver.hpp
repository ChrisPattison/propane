// Call with graph 
// Call into object with spin vector vector
//
#include <vector>
#include <functional>
#include "graph.hpp"

namespace propane
{
/** Driver for Hardware Monte Carlo accelerator
 */
class MonteCarloDriver {
    const int kVendorID = 0x7075;
    const int kDeviceID = 0x1337;

    const std::size_t kMemSize = 0x10'000'000;

    const std::size_t kSpinBase = 0x4000;
    const std::size_t kLutAddr = 0x3FC0;

    const std::size_t kDimensionAddr = 0x0000;
    const std::size_t kSpinAddr = 0x0004;
    const std::size_t kSweepAddr = 0x0008;
    const std::size_t kSeedAddr = 0x000c;
    const std::size_t kLutEntriesAddr = 0x0010;

    volatile std::uint8_t* bar0_;

    std::vector<double> local_field_;
    int bar0_disc_;

    int bus_;
    int device_;

    bool MapBar();
public:
    MonteCarloDriver();

    ~MonteCarloDriver();

    MonteCarloDriver& operator=(const MonteCarloDriver&) = delete;

    MonteCarloDriver& operator=(MonteCarloDriver&&);

    MonteCarloDriver(const MonteCarloDriver&) = delete;

    MonteCarloDriver(MonteCarloDriver&&);
/** Seed RNGs
 */
    void SeedRng(std::uint32_t seed);
/** Perform sweeps Monte Carlo sweeps on replica
 */
    void Sweep(std::vector<std::uint32_t>* replica, std::uint32_t sweeps);

/** Wait for completion of sweeps
 */
    void CompletionWait();

/** Set dE tables
 */
    void SetGraph(Graph& structure);
/** Set acceptance probabilities
 * This should not assume a Hamiltonian in the future
 * Assumes that indices match machine for now
 */
    void SetProb(double beta);
};
}