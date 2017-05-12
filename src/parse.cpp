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

namespace propane::io 
{
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
        for(auto& item : tree.get_child("schedule")) {
            config->schedule.emplace_back();
            config->schedule.back().beta = item.second.get<double>("beta");
            config->schedule.back().population_fraction = item.second.get<double>("population_fraction", 1.0);
            config->schedule.back().sweeps = item.second.get("sweeps", 10);
            config->schedule.back().overlap_dist = item.second.get("overlap_hist", false);
            config->schedule.back().energy_dist = item.second.get("energy_hist", false);
        }
    } catch(std::exception& e) {
        util::Check(false, "Config parsing failed.");
    }
    std::sort(config->schedule.begin(), config->schedule.end(), [](const auto& left, const auto& right) {return left.beta < right.beta;});
}
}