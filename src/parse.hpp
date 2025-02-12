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
 
#pragma once
#include "graph.hpp"
#include "parallel_population_annealing.hpp"
#include <iostream>
#include <string>
#include <cstdint>

namespace propane { namespace io {

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
}}
