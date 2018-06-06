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
 
#include "output.hpp"
#include <iomanip>
#include <iostream>
#include <cfenv>
#include <string>
#include "compare.hpp"

namespace psqa { namespace io {   
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

void ColumnNames() {
    std::cout << std::right << std::setw(kHeaderWidth)
        << "Beta" << std::setw(kHeaderWidth)
        << "Gamma" << std::setw(kHeaderWidth)
        << "Lambda" << std::setw(kHeaderWidth)
        << "D" << std::setw(kHeaderWidth)
        << "P" << std::setw(kHeaderWidth)
        << "Sweeps" << std::setw(kHeaderWidth)
        << "<E>" << std::setw(kHeaderWidth) 
        << "<E^2>" << std::setw(kHeaderWidth) 
        << "QR" << std::setw(kHeaderWidth) 
        << "R" << std::setw(kHeaderWidth) 
        << "E_MIN" << std::setw(kHeaderWidth) 
        << "R_MIN" << std::setw(kHeaderWidth) 
        << "S_f" << std::setw(kHeaderWidth) 
        << "rho_t" << std::setw(kHeaderWidth) 
        << "MC_Walltime" << std::setw(kHeaderWidth) 
        << "Total_Walltime" << std::setw(kHeaderWidth) 
        << "Total_Sweeps";
}

void Results(PopulationAnnealing::Result& r) {
    std::cout << std::setprecision(10) << std::scientific << std::setw(kWidth)
        << r.beta << " " << std::setw(kWidth) 
        << r.gamma << " " << std::setw(kWidth) 
        << r.lambda << " " << std::setw(kWidth) 
        << r.d << " " << std::setw(kWidth) 
        << r.p << " " << std::setw(kWidth) 
        << r.sweeps << " " << std::setw(kWidth)
        << r.average_energy << " " << std::setw(kWidth) 
        << r.average_squared_energy << " " << std::setw(kWidth) 
        << r.norm_factor << " " << std::setw(kWidth) 
        << r.population << " " << std::setw(kWidth) 
        << r.ground_energy << " " << std::setw(kWidth) 
        << r.grounded_replicas << " " << std::setw(kWidth) 
        << r.entropy << " " << std::setw(kWidth) 
        << r.mean_square_family_size << " " << std::setw(kWidth)
        << r.montecarlo_walltime << " " << std::setw(kWidth)
        << r.total_time << " " << std::setw(kWidth)
        << r.total_sweeps << " ";
}
}}
