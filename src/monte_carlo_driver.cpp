#include "monte_carlo_driver.hpp"
#include "string_util.hpp"
#include <sys/mman.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <exception>
#include <string>
#include <cmath>
#include <limits>

namespace propane
{
MonteCarloDriver::MonteCarloDriver() {
    bus_ = -1;
    device_ = -1;
    // std::fstream devices_file("/proc/bus/pci/devices", std::ios::in);
    // Grab the first valid device for now.
    // This will be moved to a seperate class eventually
    // do {
    //     std::string line;
    //     std::getline(devices_file, line);
    //     std::vector<std::string> tokens = utilities::Split(line, ' ', utilities::kStringSplitOptions_RemoveEmptyEntries);
    //     if(tokens.size() >= 2) {
    //         int vendor = std::stoi(tokens[1].substr(0, 4), 0, 16);
    //         int id = std::stoi(tokens[1].substr(4, 4), 0, 16);
    //         if(vendor == kVendorID && id == kDeviceID) {
    //             bus_ = std::stoi(tokens[0].substr(0,2), 0, 16);
    //             device_ = std::stoi(tokens[0].substr(2,2), 0, 16);
    //             break;
    //         }
    //     }
    // }while(!(devices_file.rdstate & (std::ios::eofbit | std::ios::failbit)));

    // if(bus_ == -1 && device_ == -1) {
    //     throw std::exception("Valid Accelerator not found");
    // }

    // Testing hardcode
    bar0_disc_ = open("/sys/bus/pci/devices/0000:81:00.0/resource0", O_RDWR | O_SYNC);
    bar0_ = reinterpret_cast<volatile std::uint8_t*>(mmap(nullptr, kMemSize, PROT_READ | PROT_WRITE, MAP_SHARED, bar0_disc_, 0));
    if(!bar0_) {
        // throw std::exception("mmap failed.");
    }
}

MonteCarloDriver::~MonteCarloDriver() {
    close(bar0_disc_);
    munmap((void*)bar0_, kMemSize);
}

void MonteCarloDriver::SeedRng(std::uint32_t seed) {
    volatile std::uint32_t* seed_register = reinterpret_cast<volatile std::uint32_t*>(&bar0_[kSeedAddr]);
    // Force reseed
    *seed_register = seed;
}

void MonteCarloDriver::Sweep(std::vector<std::uint32_t>* replica, std::uint32_t sweeps) {
    volatile std::uint32_t* spin_base = reinterpret_cast<volatile std::uint32_t*>(&bar0_[kSpinBase]);
    volatile std::uint32_t* sweep_register = reinterpret_cast<volatile std::uint32_t*>(&bar0_[kSweepAddr]);
    
    for(int i = 0; i < replica->size(); ++i) {
        spin_base[i] = (*replica)[i];
    }
    
    *sweep_register = sweeps;
    CompletionWait();

    for(int i = 0; i < replica->size(); ++i) {
         (*replica)[i] = spin_base[i];
    }
}

void MonteCarloDriver::CompletionWait() {
    volatile std::uint32_t* sweep_register = reinterpret_cast<volatile std::uint32_t*>(&bar0_[kSweepAddr]);
    do {
        ;
    }while(*sweep_register > 0);
}

void MonteCarloDriver::SetGraph(Graph& structure) {
    std::uint32_t spins = *reinterpret_cast<volatile std::uint32_t*>(&bar0_[kSpinAddr]);
    std::uint32_t dimension = *reinterpret_cast<volatile std::uint32_t*>(&bar0_[kDimensionAddr]);
    std::uint32_t lutentries = *reinterpret_cast<volatile std::uint32_t*>(&bar0_[kLutEntriesAddr]);
    std::uint32_t length = std::rint(std::pow(spins, 1.0/dimension));

    assert(spins == structure.size());
    assert(dimension == 2);
    int lut_size = 1<<(2*dimension);
    assert(lut_size == lutentries);
    local_field_.resize(spins * lut_size);

    // Site
    for(int vertex = 0; vertex < structure.size(); ++vertex) {
        // Single LUT entry
        for(int entry = 0; entry < lut_size; ++entry) {
            double sum = 0.0;
            // sum over each dimension
            for(int dim = 0; dim < dimension; ++dim) {
                int stride = std::pow(length, dim);
                int neighbor_vertex = (vertex - stride + structure.size()) % structure.size();
//                sum += structure.Adjacent().innerVector(vertex).coeff(neighbor_vertex) * (entry & (1<<(dim*2)) ? 1 : -1);
                sum += -1.0 * (entry & (1<<(dim*2))) == 0 ? -1 : 1;

                neighbor_vertex = (vertex + stride + structure.size()) % structure.size();
//                sum += structure.Adjacent().innerVector(vertex).coeff(neighbor_vertex) * (entry & (1<<(dim*2+1)) ? 1 : -1);
                sum += -1.0 * (entry & (1<<(dim*2+1))) == 0 ? -1 : 1;
            }
            // sum -= structure.Fields()(vertex);
            local_field_[vertex*lut_size + entry] = sum;
            // Table order?????!!!
            // local_field_.at((vertex+1)*lut_size - entry - 1) = sum;
        }
    }
}

void MonteCarloDriver::SetProb(double beta) {
    volatile std::uint32_t* lut_register = reinterpret_cast<volatile std::uint32_t*>(&bar0_[kLutAddr]);
    // 31 bits set 0x7FFFFFFF
    std::uint32_t fixed_mask = 0x7FFFFFFF;// ~0U>>1;

    for(int i = 0; i < local_field_.size(); ++i) {
        // set MSB is 1
        // MSB is s where s -> s' is the non-unity transition

        // local_field_[i] is for the 1 -> -1 transition
        // 31 bit fixed point probability in [0,1)
        auto probability = static_cast<std::uint32_t>(std::exp(-2*std::abs(local_field_[i]*beta)) * static_cast<double>(fixed_mask));


        assert(probability < 0x80000000);
        if(local_field_[i] > 0) {
            probability |= ~fixed_mask;
        }

        *lut_register = probability;
    }
}
}