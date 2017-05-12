#pragma once
#include <vector>

namespace propane 
{
class PopulationAnnealingBase {
    public:
/** Observables for a single step.
 */
    struct Result {
        struct Histogram {
            double bin;
            double value;
        };
        std::vector<Histogram> overlap;
        std::vector<Histogram> link_overlap;
        std::vector<Histogram> energy_distribution;
        double beta;
        int population;
        double norm_factor;
        double average_energy;
        double average_squared_energy;
        double ground_energy;
        int grounded_replicas;
        double entropy;
        std::uint64_t seed;
        int sweeps;
        double mean_square_family_size;
        double montecarlo_walltime;
    };
/** Parameters for a single step.
 */
    struct Schedule {
        double beta;
        int sweeps = 10;
        bool overlap_dist = false;
        bool energy_dist = false;
        double population_fraction = 1.0;
    };
/** Parameters for entire run.
 */
    struct Config {
        int population;
        std::uint64_t seed;
        std::vector<PopulationAnnealingBase::Schedule> schedule;
        double population_ratio;
        double population_slope;
        double population_shift;
    };
};
}