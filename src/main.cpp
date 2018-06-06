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
#include "graph.hpp"
#include "types.hpp"
#include "string_util.hpp"
#include "version.hpp"
#include "output.hpp"
#include "population_annealing.hpp"
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <map>

/** Read model and config for regular PA
 */
void SinglePaPre(std::string& config_path, std::string& bond_path, psqa::Graph* model, psqa::PopulationAnnealing::Config* config) {
    auto file = std::ifstream(config_path);
    psqa::io::ConfigParse(file, config);
    file.close();

    file = std::ifstream(bond_path);
    *model = psqa::io::IjjParse(file);
    file.close();

    psqa::io::Header(*model, config_path, bond_path);
}

/** Output Data for regular PA
 */
void SinglePaPost(std::vector<psqa::PopulationAnnealing::Result>& results, psqa::Graph& model) {
    psqa::io::ColumnNames();
    std::cout << std::endl;
    for(auto& r : results) {
        psqa::io::Results(r);
        std::cout << std::endl;
    }
}

void SinglePa(std::string config_path, std::string bond_path) {
    psqa::Graph model;
    psqa::PopulationAnnealing::Config config;
    SinglePaPre(config_path, bond_path, &model, &config);

    psqa::PopulationAnnealing population_annealing(model, config);
    auto results = population_annealing.Run();

    SinglePaPost(results, model);
}

enum ModeOption{
    kModeOptionSingle
};

int main(int argc, char** argv) {
    // Parse Arguments
    boost::program_options::options_description description("Options");
    boost::program_options::positional_options_description positional_description;
    positional_description.add("config", 1);
    positional_description.add("bondfile", 1);

    description.add_options()
        ("help,h", "help message")
        ("config", "configuration file")
        ("version,v", "version number")
        ("bondfile", "file containing graph and couplers")
        ("mode,m", boost::program_options::value<std::string>()->default_value("1"), "select run mode <1>");
    boost::program_options::variables_map var_map;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv)
        .options(description).positional(positional_description).run(), var_map);
    boost::program_options::notify(var_map);

    // Print Help
    if(var_map.count("help") || argc == 1) {
        std::cout << "Population Simulated Quantum Annealing";
        std::cout << "C. Pattison" << std::endl << std::endl;
        std::cout << "Usage: " << argv[0] << " [options] <config> <bondfile>" << std::endl;
        std::cout << description << std::endl;
        return EXIT_SUCCESS;
    }

    if(var_map.count("version")) {
        std::cout << "Population Simulated Quantum Annealing";
        std::cout << "Branch: " << psqa::version::kRefSpec << std::endl;
        std::cout << "Commit: " << std::string(psqa::version::kCommitHash).substr(0, 8) << std::endl;
        std::cout << "Build:  " << psqa::version::kBuildType << std::endl;
        std::cout << "Built:  " << psqa::version::kBuildTime << std::endl;
        return EXIT_SUCCESS;
    }

    // Select PA implementation
    std::map<std::string, ModeOption> selector_map;
    ModeOption selection;
    if(selector_map.count(var_map["mode"].as<std::string>()) == 0) {
        selection = kModeOptionSingle;
    }else {
        selection = selector_map.at(var_map["mode"].as<std::string>());
    }

    switch(selection) {
        case kModeOptionSingle : SinglePa(var_map["config"].as<std::string>(), var_map["bondfile"].as<std::string>()); break;
    }

    return EXIT_SUCCESS;
}
