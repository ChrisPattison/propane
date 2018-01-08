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

#include "spin_vector_monte_carlo.hpp"
#include "compare.hpp"

namespace propane {

double SpinVectorMonteCarlo::Hamiltonian(const StateVector& replica) {
    double energy = 0.0;
    for(std::size_t k = 0; k < structure_.Adjacent().outerSize(); ++k) {
        for(Eigen::SparseTriangularView<Eigen::SparseMatrix<EdgeType>,Eigen::Upper>::InnerIterator 
            it(structure_.Adjacent().triangularView<Eigen::Upper>(), k); it; ++it) {
            energy += replica[k][0] * it.value() * replica[it.index()][0];
        }
        energy -= replica[k] * field_[k];
        energy *= replica.lambda; // Problem Hamiltonian Strength

        energy -= replica[k][0] * replica.gamma;
    }
    return energy;
}

SpinVectorMonteCarlo::StateVector SpinVectorMonteCarlo::Project(const StateVector& replica) {
    StateVector projected;
    projected.resize(replica.size());
    for( std::size_t k = 0; k < replica.size(); ++k ) {
        projected[k] = replica[k] * FieldType(1.0, 0.0) > 0 ? FieldType(1.0, 0.0) : FieldType(-1.0, 0.0);
    }
    return projected;
}

double SpinVectorMonteCarlo::ProjectedHamiltonian(const StateVector& projected) {
    double energy = 0.0;
    for(std::size_t k = 0; k < structure_.Adjacent().outerSize(); ++k) {
        for(Eigen::SparseTriangularView<Eigen::SparseMatrix<EdgeType>,Eigen::Upper>::InnerIterator 
            it(structure_.Adjacent().triangularView<Eigen::Upper>(), k); it; ++it) {
            energy += projected[k][0] * it.value() * projected[it.index()][0];
        }
        energy -= projected[k] * FieldType(field_[k][0], 0.);
    }
    return energy;
}

FieldType SpinVectorMonteCarlo::LocalField(StateVector& replica, int vertex) {
    FieldType h;
    for(Eigen::SparseMatrix<EdgeType>::InnerIterator it(structure_.Adjacent(), vertex); it; ++it) {
        h += FieldType(it.value() * replica[it.index()][0], 0.0);
    }
    h -= field_[vertex];
    h *= replica.lambda;
    h -= FieldType(0., replica.gamma);
    return h;
}

double SpinVectorMonteCarlo::DeltaEnergy(StateVector& replica, int vertex, FieldType new_value) {
    return (new_value - replica[vertex]) * LocalField(replica, vertex);
}

void SpinVectorMonteCarlo::MicroCanonicalSweep(StateVector& replica, int sweeps) {
    for(std::size_t k = 0; k < sweeps; ++k) {
        for(std::size_t i = 0; i < replica.size(); ++i) {
            auto vertex = i;
            
            // get local field
            auto h = LocalField(replica, vertex);
            auto new_spin = ((2*h*(replica[vertex]*h))/(h*h))-replica[vertex];
            replica[vertex] = new_spin;
        }
    }
}

void SpinVectorMonteCarlo::MetropolisSweep(StateVector& replica, int sweeps) {
    for(std::size_t k = 0; k < sweeps; ++k) {
        for(std::size_t i = 0; i < replica.size(); ++i) {
            int vertex = i;
            auto new_value = VertexType(rng_.Probability());
            double delta_energy = DeltaEnergy(replica, vertex, new_value);
            
            //round-off isn't a concern here
            if(MetropolisAcceptedMove(delta_energy, replica.beta)) {
                replica[vertex] = new_value;
            }
        }
    }
}

void SpinVectorMonteCarlo::HeatbathSweep(StateVector& replica, int sweeps) {
    for(std::size_t k = 0; k < sweeps; ++k) {
        for(std::size_t i = 0; i < replica.size(); ++i) {
            int vertex = rng_.Range(replica.size());
            auto h = LocalField(replica, vertex);
            auto h_mag = std::sqrt(h*h);
            auto h_unit = h / h_mag;
            double x = -std::log(1 + rng_.Probability() * (std::exp(-2 * replica.beta * h_mag) - 1)) / (replica.beta * h_mag) - 1.;
            auto h_perp = FieldType(h_unit[1], -h_unit[0]);
            h_perp *= rng_.Probability() < 0.5 ? -1 : 1;
            replica[vertex] = h_unit * x + h_perp * std::sqrt(1.0-x*x);
        }
    }
}

double SpinVectorMonteCarlo::Overlap(StateVector& alpha, StateVector& beta) {
    return std::inner_product(alpha.begin(), alpha.end(), beta.begin(), 0.0) / structure_.size();
}

double SpinVectorMonteCarlo::LinkOverlap(StateVector& alpha, StateVector& beta) {
    double ql = 0;
    for(std::size_t k = 0; k < structure_.Adjacent().outerSize(); ++k) {
        for(Eigen::SparseTriangularView<Eigen::SparseMatrix<EdgeType>,Eigen::Upper>::InnerIterator 
            it(structure_.Adjacent().triangularView<Eigen::Upper>(), k); it; ++it) {
            ql += alpha[k] * beta[k] * alpha[it.index()] * beta[it.index()];
        }
    }
    return ql / structure_.edges();
}

bool SpinVectorMonteCarlo::MetropolisAcceptedMove(double delta_energy, double beta) {
    if(delta_energy < 0.0) {
        return true;
    }
    
    double acceptance_prob_exp = -delta_energy*beta;
    return AcceptedMove(acceptance_prob_exp);
}

bool SpinVectorMonteCarlo::AcceptedMove(double log_probability) {
    double test = rng_.Probability();
    // Compute bound on log of test number
    auto bound = log_lookup_(test);

    if(bound.upper < log_probability) {
        return true;
    }else if(bound.lower > log_probability) {
        return false;
    }
    // Compute exp if LUT can't resolve it
    return std::exp(log_probability) > test;
}

void SpinVectorMonteCarlo::TransverseField(StateVector& replica, double magnitude, double p_magnitude) {
    replica.gamma = magnitude;
}
}