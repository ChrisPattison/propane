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
 
#include "parse.hpp"
#include <vector>
#include <algorithm>
#include <iterator>
#include <exception>
#include <tuple>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "string_util.hpp"
#include "compare.hpp"

namespace propane { namespace io {
Graph IjjParse(std::istream& file) {
    Graph model;
    try {
        std::string buffer(kBufferSize, '\0');
        std::vector<std::string> tokens;
        
        do {
            std::getline(file, buffer);
        } while(util::Split(buffer, kWhitespaceTokens, util::kStringSplitOptions_RemoveEmptyEntries)[0][0]=='#');

        auto header = util::Split(buffer, kWhitespaceTokens, util::kStringSplitOptions_RemoveEmptyEntries);
        int coupler_multiplier = std::stoi(header[1]);
        std::vector<std::tuple<double, double, double>> values;

        do {
            std::getline(file, buffer);
            buffer = util::Split(buffer, kCommentToken)[0];
            tokens = util::Split(buffer, kWhitespaceTokens, util::kStringSplitOptions_RemoveEmptyEntries);
            if(tokens.size() >= 3) {
                values.push_back(std::make_tuple(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]) / coupler_multiplier));
            }
        } while(!file.eof());
        
        model.Resize(stoi(header[0]), values.size());
        for(auto& v : values) {
            int i = std::get<0>(v);
            int j = std::get<1>(v);
            EdgeType coeff = std::get<2>(v);
            if(i != j) {
                model.AddEdge(i, j, coeff);
                model.AddEdge(j, i, coeff);
            }else {
                model.SetField(i, coeff);
            }
        }
        util::Check(model.IsConsistent(), "Missing edge somewhere.");
        util::Check(model.Adjacent().size() > 0, "No Elements.");
    } catch(std::exception& e) {
        util::Check(false, "Input parsing failed.");
    }
    return model;
}

void ConfigParse(std::istream& file, PopulationAnnealing::Config* config) {
    try {
        boost::property_tree::ptree tree;
        boost::property_tree::read_json(file, tree);

        config->population = tree.get<int>("population");
        std::stringstream converter(tree.get<std::string>("seed", "0"));
        converter >> std::hex >> config->seed;
        int default_sweeps = tree.get<int>("default_sweeps", 10);
        int default_wolff = tree.get<int>("default_wolff", 0);
        config->solver_mode = tree.get<bool>("solver_mode", false);
        config->trotter_slices = tree.get<int>("trotter_slices");
        for(auto& item : tree.get_child("schedule")) {
            config->schedule.emplace_back();
            config->schedule.back().beta = item.second.get<double>("beta");
            config->schedule.back().gamma = item.second.get<double>("gamma");
            config->schedule.back().lambda = item.second.get<double>("lambda", 1.0);
            config->schedule.back().population_fraction = item.second.get<double>("population_fraction", 1.0);
            config->schedule.back().sweeps = item.second.get("sweeps", default_sweeps);
            config->schedule.back().wolff_sweeps = item.second.get("wolff_sweeps", default_wolff);
            config->schedule.back().heat_bath = item.second.get("heat_bath", false);
            config->schedule.back().compute_observables = item.second.get("compute_observables", true);
        }
    } catch(std::exception& e) {
        util::Check(false, "Config parsing failed.");
    }
    std::stable_sort(config->schedule.begin(), config->schedule.end(), [](const auto& left, const auto& right) {return (left.beta != right.beta) ?  (left.beta < right.beta) : left.gamma > right.gamma;});
}
}}