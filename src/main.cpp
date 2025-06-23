#include <math.h>
#include <random>
#include <sstream>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <array>
#include <thread>
#include <sys/stat.h>
#include <string>
#include <chrono>
#include <thread>

#include "Electron.h"
#include "ProgressBar.h"

int main()
{
    std::cout << "V/cm (multiple values separated by spaces): ";
    std::string line;
    double volts;
    std::vector<double> volts_list;
    std::getline(std::cin, line);
    std::istringstream stream(line);
    while (stream >> volts)
        volts_list.push_back(volts);

    double elec_energy;
    std::cout << "Electron starting energy [eV]:";
    std::cin >> elec_energy;
    while (std::cin.fail()) {
      std::cin.clear();
      std::cin.ignore();
      std::cin >> elec_energy;
    }
    

    double angle;
    std::cout << "Angle of ionization track with Electric field direction [Degrees]: ";
    std::cin >> angle;
    while (std::cin.fail()){
        std::cin.clear();
        std::cin.ignore();
        std::cin >> angle;
    }

    double density;
    std::cout << "Distance between ions [nm]: ";
    std::cin >> density;
    while (std::cin.fail()){
        std::cin.clear();
        std::cin.ignore();
        std::cin >> density;
    }


    double cutoff;
    std::cout << "Stop After (s): ";
    std::cin >> cutoff;
    while (std::cin.fail()) {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> cutoff;
    }

    int write_every;
    std::cout << "Record Every n Measurements: ";
    std::cin >> write_every;
    while (std::cin.fail()) {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> write_every;
    }

    int cores = std::thread::hardware_concurrency();
    if (cores < 1) cores = 1;
    cores = 1;
    int batches = 1;

    std::cout << "Batches (in groups of " + std::to_string(cores) + "): ";
    std::cin >> batches;
    while (std::cin.fail()) {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> batches;
    }


    int debug = 0;
    int status = 0;
    /*
    std::cout << "debug? ";
    std::cin >> debug;
    while (std::cin.fail()) {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> debug;
    }
    */

    // Distributes the requested number of simulations among the processor cores
    for (int j = 0; j < volts_list.size(); j++) {
        std::thread branches[cores];
        ProgressBar bar(cores);

        for (int k = 0; k < cores; k++) {
	  branches[k] = std::thread(generate_plot, int(volts_list[j]), elec_energy, angle, density, cutoff, cores, write_every, k, batches, debug, status, std::ref(bar));   
        }
    
        for (int k = 0; k < cores; k++) {
            branches[k].join();
        }

        std::cout << "\n";
    }

    std::cout << "Complete." << std::endl;

    return 0;
}
