/*!
* \file GrapheneFloquet.cpp
*
*
* \author Author: D. Gagnon <denisg6@hotmail.com>
*/

// Include some headers
#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <complex>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "GrapheneFloquet.hpp"
#include "Utils.hpp"

/// Main function
int main(int argc, char *argv[])
{

    // Parse parameter file
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("GrapheneFloquet.ini", pt);

    // Problem parameters
    double freq  =          std::stof(pt.get<std::string>("Parameters.frequency"));  // Angular frequency in Hz
    double E0 =             std::stof(pt.get<std::string>("Parameters.E0"));         // Peak electric field in V/m
    size_t nblocks =  std::stoi(pt.get<std::string>("Parameters.nblocks"));    // Envelope frequency, in units of freq

    // Parameter sweep
    double xmin =   std::stof(pt.get<std::string>("Sweep.xmin")); // Minimum frequency in sweep
    double xmax =   std::stof(pt.get<std::string>("Sweep.xmax")); // Minimum frequency in sweep
    size_t x_elem =    std::stoi(pt.get<std::string>("Sweep.xnum"));

    // Parameter sweep
    double ymin =   std::stof(pt.get<std::string>("Sweep.ymin")); // Minimum frequency in sweep
    double ymax =   std::stof(pt.get<std::string>("Sweep.ymax")); // Minimum frequency in sweep
    size_t y_elem =    std::stoi(pt.get<std::string>("Sweep.ynum"));

    // Process input variables (angular frequency units)
    double omega = 2.0*M_PI*freq;     // Angular frequency

    // Meshgrid (vectors of parameters)
    auto xvec = arma::linspace(xmin, xmax, x_elem);
    auto yvec = arma::linspace(ymin, ymax, y_elem);

    // Prepare output file
    // std::ofstream outfile;
    // outfile.open("quasienergies.dat");

    // Prepare output file for transition probabilities
    std::ofstream outfile2;
    // outfile2.open("kmap.dat");

    // Variables to store probability values
    auto probvec = arma::mat(y_elem,x_elem);
    auto probvec_sigma = arma::mat(y_elem,x_elem);

    // Variables to store quasienergy values
    auto quasi0 = arma::mat(y_elem,x_elem);
    auto quasi1 = arma::mat(y_elem,x_elem);
    auto quasi2 = arma::mat(y_elem,x_elem);

    // Variables for loops
    size_t id, id2;
    double prob = 0.0;
    double prob_sigma = 0.0; // Trans. prob. between sigma_z eigenstates

    // Loop for each K value and compute probability
    # pragma omp parallel for default(shared) private (id, id2, prob, prob_sigma)
    for (id=0; id < y_elem; id++)
    {
        // Initialize probabilities (will be passed by reference)
        for (id2=0; id2 < x_elem; id2++)
        {
            // Compute energies and transition probability
            auto Energies = QuasiEnergies(xvec[id2],yvec[id],omega,E0,nblocks, prob, prob_sigma);

            // Store in arrays
            probvec(id,id2) = prob;
            probvec_sigma(id,id2) = prob_sigma;

            quasi0(id,id2) = Energies[Energies.n_elem/2 - 1]; // Quasi-energ.
            quasi1(id,id2) = Energies[Energies.n_elem/2];
            quasi2(id,id2) = Energies[Energies.n_elem/2 + 1];
        }

    }

    // Save Floquet Data
    probvec.save("probability.dat", arma::raw_ascii);

    // Save parameters vector
    xvec.save("xvec.dat", arma::raw_ascii);
    yvec.save("yvec.dat", arma::raw_ascii);


    quasi0.save("quasi0.dat", arma::raw_ascii);
    quasi1.save("quasi1.dat", arma::raw_ascii);
    quasi2.save("quasi2.dat", arma::raw_ascii);

    // Save parameters vector
    xvec.save("xvec.dat", arma::raw_ascii);
    yvec.save("yvec.dat", arma::raw_ascii);

    return 0;

    /*
    auto Gamma = GammaValues(2.30,0.0,316227766016837.94,1e10,1);

    for(auto it = Gamma.begin(); it != Gamma.end(); it++)
    {
        std::cout << *it << std::endl;
    }

    double prob = 0.0;
    double prob_sigma = 0.0;
    auto Quasi = QuasiEnergies(2.30,0.0,316227766016837.94,1e10,1, prob, prob_sigma);
    */
}
