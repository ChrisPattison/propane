#include "output.hpp"
#include <iomanip>
#include <iostream>
#include <cfenv>
#include <string>
#include "compare.hpp"
#include "parallel.hpp"

namespace propane::io {   
void Header(Graph& model, std::string config_path, std::string bond_path) {
    std::cout << "# Parallel Optimized Population Annealing V" << version::kMajor << "." << version::kMinor << std::endl;
    std::cout << "# C. Pattison" << std::endl;
    std::cout << "# Branch: " << version::kRefSpec << std::endl;
    std::cout << "# Commit: " << version::kCommitHash << std::endl;
    std::cout << "# Built: " << version::kBuildTime << std::endl;
    std::cout << "# Config: " << config_path << std::endl;
    std::cout << "# Input: " << bond_path << std::endl;
    std::cout << "# Spins: " << model.size() << std::endl;
}

void MpiHeader(parallel::Mpi& parallel) {
    util::Check(parallel.size() % 4 == 0, "Number of cores must be a multiple of 4");
    std::cout << "# Cores: " << parallel.size() << std::endl;
}

void ColumnNames() {
    std::cout << std::right << std::setw(kHeaderWidth)
        << "Beta" << std::setw(kHeaderWidth)
        << "Sweeps" << std::setw(kHeaderWidth)
        << "<E>" << std::setw(kHeaderWidth) 
        << "<E^2>" << std::setw(kHeaderWidth) 
        << "QR" << std::setw(kHeaderWidth) 
        << "R" << std::setw(kHeaderWidth) 
        << "E_MIN" << std::setw(kHeaderWidth) 
        << "R_MIN" << std::setw(kHeaderWidth) 
        << "S_f" << std::setw(kHeaderWidth) 
        << "rho_t" << std::setw(kHeaderWidth) 
        << "MC_Walltime";
}

void MpiColumnNames() {
    std::cout << std::setw(kHeaderWidth) 
        << "Redist_Walltime" << std::setw(kHeaderWidth) 
        << "Obs_Walltime" << std::setw(kHeaderWidth) 
        << "R_f_MAX";
}

void Results(PopulationAnnealing::Result& r) {
    std::cout << std::setprecision(10) << std::scientific << std::setw(kWidth)
        << r.beta << " " << std::setw(kWidth) 
        << r.sweeps << " " << std::setw(kWidth)
        << r.average_energy << " " << std::setw(kWidth) 
        << r.average_squared_energy << " " << std::setw(kWidth) 
        << r.norm_factor << " " << std::setw(kWidth) 
        << r.population << " " << std::setw(kWidth) 
        << r.ground_energy << " " << std::setw(kWidth) 
        << r.grounded_replicas << " " << std::setw(kWidth) 
        << r.entropy << " " << std::setw(kWidth) 
        << r.mean_square_family_size << " " << std::setw(kWidth)
        << r.montecarlo_walltime << " ";
}

void MpiResults(ParallelPopulationAnnealing::Result& r) {
    std::cout << std::setprecision(10) << std::scientific << std::setw(kWidth)
        << r.redist_walltime << " " << std::setw(kWidth) 
        << r.observables_walltime << " " << std::setw(kWidth)
        << r.max_family_size;
}

void Histograms(std::vector<PopulationAnnealing::Result>& results) {
    std::cout << std::endl << kMagicString << std::endl << "# Histograms" << std::endl;
    for(auto r : results) {
        for(auto q : r.overlap) {
            std::cout << kHistChar << "q " 
                << std::setprecision(4) << std::fixed << r.beta << " "
                << std::setprecision(10) << std::scientific << std::setw(kWidth) << q.bin << " "
                << std::setw(kWidth) << q.value << std::endl;
        }
        for(auto ql : r.link_overlap) {
            std::cout << kHistChar << "ql " 
                << std::setprecision(4) << std::fixed << r.beta << " "
                << std::setprecision(10) << std::scientific << std::setw(kWidth) << ql.bin << " "
                << std::setw(kWidth) << ql.value << std::endl;
        }
        for(auto E : r.energy_distribution) {
            std::cout << kHistChar << "E "
                << std::setprecision(4) << std::fixed << r.beta << " "
                << std::setprecision(10) << std::scientific << std::setw(kWidth) << E.bin << " "
                << std::setw(kWidth) << E.value << std::endl;
        }
    }
}


void IjjDump(Graph& model, std::ostream& stream) {
    stream << std::endl << kMagicString << std::endl << "# Input" << std::endl;
    auto round_mode = std::fegetround();
    std::fesetround(FE_TONEAREST);
    stream << model.size() << kOutputSeperator << kOutputCouplerCoeff << std::endl;
    // Coefficients
    for(std::size_t k = 0; k < model.Adjacent().outerSize(); ++k) {
        for(Eigen::SparseTriangularView<Eigen::SparseMatrix<EdgeType>,Eigen::Upper>::InnerIterator 
            it(model.Adjacent().triangularView<Eigen::Upper>(), k); it; ++it) {
            double value = kOutputCouplerCoeff * it.value();
            stream << k << kOutputSeperator << it.index() << kOutputSeperator;
            if(util::FuzzyUlpCompare(value, std::lrint(value), 100)) {
                stream << std::lrint(value);
            }else {
                stream << value;
            }
            stream << std::endl;
        }
    }
    // Fields
    if(model.has_field()) {
        for(std::size_t k = 0; k < model.Fields().size(); ++k) {
            stream << k << kOutputSeperator << k << kOutputSeperator;
            double value = kOutputCouplerCoeff * model.Fields()[k];
            if(util::FuzzyUlpCompare(value, std::lrint(value), 100)) {
                stream << std::lrint(value);
            }else {
                stream << value;
            }
            stream << std::endl;
        }
    }
    std::fesetround(round_mode);
}

}
