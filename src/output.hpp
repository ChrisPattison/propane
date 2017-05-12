#pragma once
#include "graph.hpp"
#include "parallel_population_annealing.hpp"
#include "population_annealing.hpp"
#include "fpga_population_annealing.hpp"
#include <iostream>
#include <string>
#include <cstdint>
#include "version.hpp"
#include "string_util.hpp"
#include "parallel.hpp"

namespace propane::io {
constexpr int kWidth = 18;
constexpr int kHeaderWidth = kWidth + 1;
const auto kMagicString = "%%%---%%%";
const auto kHistChar = "|";
static constexpr char kOutputSeperator = ' ';
static constexpr int kOutputCouplerCoeff = 1;

void Header(Graph& model, std::string config_path, std::string bond_path);

void MpiHeader(parallel::Mpi& parallel);

void ColumnNames();

void MpiColumnNames();

void Results(PopulationAnnealing::Result& result);

void MpiResults(ParallelPopulationAnnealing::Result& result);

void Histograms(std::vector<PopulationAnnealing::Result>& results);

/** Opposite of IjjParse: Dumps bonds to file.
 */
void IjjDump(Graph& model, std::ostream& stream);
}
