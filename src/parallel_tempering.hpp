/* Copyright (c) 2017 C. Pattison
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
 
#pragma once
#include "types.hpp"
#include "graph.hpp"
#include "random_number_generator.hpp"
#include "parallel_tempering_base.hpp"
#include "spin_vector_monte_carlo.hpp"
#include <vector>
#include <utility>
#include <chrono>
#include <Eigen/Dense>

namespace propane {
class ParallelTempering : public ParallelTemperingBase, protected SpinVectorMonteCarlo {
protected:
    std::vector<StateVector> replicas_;
    std::vector<int> replica_families_;

    std::vector<Schedule> schedule_;
    std::vector<std::size_t> bin_set_;
    std::vector<double> profile_;
    std::size_t sweeps_;
    bool solver_mode_;
    bool uniform_init_;
    double planted_energy_;
public:
/** Intializes solver.
 * schedule specifies the temperature set and sweep types to do at each temperature
 * seed may be zero in which case one will be generated.
 */
    ParallelTempering(const Graph& structure, Config config);
/** Run solver and return results.
 */
    std::vector<ParallelTempering::Result> Run();
private:
/** Carry out replica exchange on replicas
 */
    void ReplicaExchange(std::vector<StateVector>& replica_set);
/** Record observables
 */
    Bin Observables(const StateVector& replica);
};
}