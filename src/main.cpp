#include "parse.hpp"
#include "graph.hpp"
#include "parallel_population_annealing.hpp"
#include "greedy_population_annealing.hpp"
#include "types.hpp"
#include "string_util.hpp"
#include "version.hpp"
#include "output.hpp"
#include <boost/program_options.hpp>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <map>

void MpiPa(std::string config_path, std::string bond_path) {
    auto file = std::ifstream(config_path);
    propane::ParallelPopulationAnnealing::Config config;
    propane::io::ConfigParse(file, &config);
    file.close();

    file = std::ifstream(bond_path);
    propane::Graph model = propane::io::IjjParse(file);
    file.close();

    parallel::Mpi parallel;
    parallel.ExecRoot([&]() {
        propane::io::Header(model, config_path, bond_path);
        propane::io::MpiHeader(parallel);
    });

    propane::ParallelPopulationAnnealing population_annealing(model, config);
    auto results = population_annealing.Run();

    parallel.ExecRoot([&]() {
        propane::io::ColumnNames();
        propane::io::MpiColumnNames();
        std::cout << std::endl;
        for(auto& r : results) {
            propane::io::Results(r);
            propane::io::MpiResults(r);
            std::cout << std::endl;
        }
        std::vector<propane::PopulationAnnealing::Result> basic_result(results.size());
        std::transform(results.begin(), results.end(), basic_result.begin(), [] (propane::ParallelPopulationAnnealing::Result& r) 
            {return static_cast<propane::PopulationAnnealing::Result>(r);});

        propane::io::IjjDump(model, std::cout);
        propane::io::Histograms(basic_result);
    });
}

/** Read model and config for regular PA
 */
void SinglePaPre(std::string& config_path, std::string& bond_path, propane::Graph* model, propane::PopulationAnnealing::Config* config) {
    auto file = std::ifstream(config_path);
    propane::io::ConfigParse(file, config);
    file.close();

    file = std::ifstream(bond_path);
    *model = propane::io::IjjParse(file);
    file.close();

    propane::io::Header(*model, config_path, bond_path);
}

/** Output Data for regular PA
 */
void SinglePaPost(std::vector<propane::PopulationAnnealing::Result>& results, propane::Graph& model) {
    propane::io::ColumnNames();
    std::cout << std::endl;
    for(auto& r : results) {
        propane::io::Results(r);
        std::cout << std::endl;
    }
    propane::io::IjjDump(model, std::cout);
    propane::io::Histograms(results);
}

void SinglePa(std::string config_path, std::string bond_path) {
    propane::Graph model;
    propane::PopulationAnnealing::Config config;
    SinglePaPre(config_path, bond_path, &model, &config);

    propane::PopulationAnnealing population_annealing(model, config);
    auto results = population_annealing.Run();

    SinglePaPost(results, model);
}

void GreedyPa(std::string config_path, std::string bond_path) {
    propane::Graph model;
    propane::GreedyPopulationAnnealing::Config config;
    SinglePaPre(config_path, bond_path, &model, &config);

    propane::GreedyPopulationAnnealing population_annealing(model, config);
    auto results = population_annealing.Run();

    SinglePaPost(results, model);
}

void FpgaPa(std::string config_path, std::string bond_path) {
    propane::Graph model;
    propane::PopulationAnnealing::Config config;
    SinglePaPre(config_path, bond_path, &model, &config);

    propane::FpgaPopulationAnnealing population_annealing(model, config);
    auto results = population_annealing.Run();

    SinglePaPost(results, model);
}

enum ModeOption{
    kModeOptionFpga,
    kModeOptionMpi,
    kModeOptionSingle,
    kModeOptionGreedy
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
        ("mode,m", boost::program_options::value<std::string>()->default_value("1"), "select run mode <1/mpi/fpga/greedy>");
    boost::program_options::variables_map var_map;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv)
        .options(description).positional(positional_description).run(), var_map);
    boost::program_options::notify(var_map);

    // Print Help
    if(var_map.count("help") || argc == 1) {
        std::cout << "Parallel Optimized Population Annealing V" << propane::version::kMajor << "." << propane::version::kMinor << std::endl;
        std::cout << "C. Pattison" << std::endl << std::endl;
        std::cout << "Usage: " << argv[0] << " [options] <config> <bondfile>" << std::endl;
        std::cout << description << std::endl;
        return EXIT_SUCCESS;
    }

    if(var_map.count("version")) {
        std::cout << "Parallel Optimized Population Annealing V" << propane::version::kMajor << "." << propane::version::kMinor << std::endl;
        std::cout << "Branch: " << propane::version::kRefSpec << std::endl;
        std::cout << "Commit: " << std::string(propane::version::kCommitHash).substr(0, 8) << std::endl;
        std::cout << "Build:  " << propane::version::kBuildType << std::endl;
        std::cout << "Built:  " << propane::version::kBuildTime << std::endl;
        return EXIT_SUCCESS;
    }

    // Select PA implementation
    std::map<std::string, ModeOption> selector_map;
    selector_map.insert({"mpi", kModeOptionMpi});
    selector_map.insert({"fpga", kModeOptionFpga});
    selector_map.insert({"greedy", kModeOptionGreedy});
    ModeOption selection;
    if(selector_map.count(var_map["mode"].as<std::string>()) == 0) {
        selection = kModeOptionSingle;
    }else {
        selection = selector_map.at(var_map["mode"].as<std::string>());
    }

    switch(selection) {
        case kModeOptionMpi : MpiPa(var_map["config"].as<std::string>(), var_map["bondfile"].as<std::string>()); break;
        case kModeOptionFpga : FpgaPa(var_map["config"].as<std::string>(), var_map["bondfile"].as<std::string>()); break;
        case kModeOptionSingle : SinglePa(var_map["config"].as<std::string>(), var_map["bondfile"].as<std::string>()); break;
        case kModeOptionGreedy : GreedyPa(var_map["config"].as<std::string>(), var_map["bondfile"].as<std::string>()); break;
    }

    return EXIT_SUCCESS;
}
