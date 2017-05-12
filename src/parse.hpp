#pragma once
#include "graph.hpp"
#include "parallel_population_annealing.hpp"
#include <iostream>
#include <string>
#include <cstdint>

namespace propane::io 
{
static constexpr std::size_t kBufferSize = 1024; // size of input buffer
static constexpr std::size_t kReservedVertices = 6; // average number of edges per vertex
static constexpr char kCommentToken = {'#'};
static const std::vector<char> kWhitespaceTokens = {' ', '\t'};

/** Reads the bonds from file and constructs the Graph accordingly.
 * Expects a header with the number of spins and scaling factor.
 * Bond format is
 * i j J_ij
 */
Graph IjjParse(std::istream& file);
/** Parses configuration file for single threaded PA and populates a Config.
 * Uses boost::property_tree to read the JSON.
 */ 
void ConfigParse(std::istream& file, PopulationAnnealing::Config* config);
}
