#include "gtest/gtest.h"
#include "population_annealing.hpp"
#include <memory>
#include <random>

namespace {
using namespace psqa;

class WolffTest : public ::testing::Test {
    struct UnitUnderTest : PopulationAnnealing {
        using PopulationAnnealing::PopulationAnnealing;
        using PopulationAnnealing::StateVector;

        using PopulationAnnealing::rng_;
        using PopulationAnnealing::structure_;

        using PopulationAnnealing::beta_;
        using PopulationAnnealing::coeff_P_;
        using PopulationAnnealing::coeff_D_;

        using PopulationAnnealing::WolffSweep;
        using PopulationAnnealing::SpatialSiteEnergy;
        using PopulationAnnealing::AcceptedMove;
    };
protected:
    std::mt19937_64 random_;
    std::unique_ptr<UnitUnderTest> uut;
    UnitUnderTest::StateVector test_state;

    /** Reference version of Wolff updates
     * Uses RNG in same way as the production version
     * Avoids fancy bit hacks to find boundaries
     */
    void WolffReference(UnitUnderTest::StateVector& replica) {
        std::vector<double> random(replica.size() * ktrotter_slices);
        std::generate(random.begin(), random.end(), [&]() { return uut->rng_(); });
        double growth_prob = 1.0 - std::exp(2.0 * uut->beta_ * uut->coeff_D_);

        for(std::size_t k = 0; k < uut->structure_.size(); ++k) {
            std::size_t random_offset = ktrotter_slices * k;
            
            std::array<bool, ktrotter_slices> grow;
            for(std::size_t i = 0; i < ktrotter_slices; ++i) {
                grow[random_offset + i] = uut->rng_() < growth_prob;
            }

            std::uint32_t seed = uut->rng_.Range(uut->structure_.size() * ktrotter_slices);
            std::uint32_t site = seed % uut->structure_.size();
            std::uint32_t seed_slice = seed / uut->structure_.size();
            VertexType site_value = replica[site];

            VertexType cluster = 1U << seed_slice;
            VertexType prev_mask =  cluster;
            // Build Cluster
            for(std::size_t i = (seed_slice + 1)%ktrotter_slices; i != seed_slice; i = (i+1)%ktrotter_slices) {
                VertexType next_mask = 1U << i;
                if(grow[i] && (static_cast<bool>(next_mask & site_value) == (prev_mask & site_value))) {
                    cluster |= next_mask;
                }
            }
            for(std::size_t i = (seed_slice + ktrotter_slices - 1)%ktrotter_slices; 
                i != seed_slice; i = (i+ ktrotter_slices - 1)%ktrotter_slices) {
                VertexType next_mask = 1U << i;
                if(grow[i] && (static_cast<bool>(next_mask & site_value) == (prev_mask & site_value))) {
                    cluster |= next_mask;
                }
            }
            
            std::array<EnergyType, ktrotter_slices> site_delta_energy;
            uut->SpatialSiteEnergy(replica, site, site_delta_energy.begin());

            EnergyType delta_energy = 0;
            for(std::size_t i = 0; i < ktrotter_slices; ++i) {
                VertexType mask = 1U << i;
                delta_energy += mask & cluster ? site_delta_energy[i] : 0;
            }
            delta_energy *= -2 * uut->coeff_P_;
            // Attempt flip
            if(uut->AcceptedMove(-delta_energy*uut->beta_)) {
                replica[site] ^= cluster;
            }
        }
    }
public:

    WolffTest() {
        Graph graph;
        graph.Resize(2);
        graph.AddEdge(0, 1, -1.0);
        graph.AddEdge(1, 0, -1.0);
        graph.SetField(1, 2.0);
        PopulationAnnealingBase::Config config;
        config.population = 1;
        config.seed = 0;
        config.schedule.emplace_back();
        config.schedule.back().beta = 10.0;
        config.schedule.back().gamma = 1e-3;
        config.schedule.back().lambda = 1.0;
        uut = std::make_unique<UnitUnderTest>(graph,  config);

        random_ = std::mt19937_64(std::random_device()());
    }

    virtual ~WolffTest() {

    }

    virtual void SetUp() {
        test_state.resize(2, 0);
        uut->beta_ = 1.0;
        uut->coeff_P_ = 1.0;
        uut->coeff_D_ = -1.0;
    }

    virtual void TearDown() {

    }
};

/** At low temperature in the GS, nothing should change
 */
TEST_F(WolffTest, ClusterFreeze) {
    uut->beta_ = 10.0;
    test_state[0] = ~test_state[0];
    uut->WolffSweep(test_state, 10);
    for(auto& v : test_state) {
        ASSERT_EQ(v & 1 ? ~v : v, 0);
    }
}

/** Quench a random config with a strong field
 */
TEST_F(WolffTest, ClusterQuench) {
    uut->beta_ = 10.0;
    test_state[0] = random();
    uut->WolffSweep(test_state, 1000);
    ASSERT_EQ(test_state[0] & 1 ? ~test_state[0] : test_state[0], 0);
}

/** Check against the reference implementation
 */
TEST_F(WolffTest, ClusterReference) {
    auto seed = std::random_device()();

    std::generate(test_state.begin(), test_state.end(), random);
    auto ref_state = test_state;
    
    uut->rng_ = RandomNumberGenerator(seed);
    WolffReference(ref_state);

    uut->rng_ = RandomNumberGenerator(seed);
    uut->WolffSweep(test_state, 1);

    ASSERT_EQ(test_state[0], ref_state[0]);
}
}