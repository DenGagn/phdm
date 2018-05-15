/*!
* \file GrapheneExp.cpp
*
* \brief Tight-binding model for graphene photo-excitation
* Solves the ODE system for a wave with a few-cycle envelope (exp)
*
* The excitation
* is a slowly varying envelope with exponential shape
* See e.g. http://link.aps.org/doi/10.1103/PhysRevB.91.045439
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
#include <boost/numeric/odeint.hpp>

#include "GrapheneExp.hpp"
#include "Utils.hpp"

/// Main function
int main(int argc, char *argv[])
{
    // Parse parameter file
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("GrapheneExp.ini", pt);

    // Problem parameters
    double E0 = std::stof(pt.get<std::string>("Parameters.E0"));            // Peak electric field in V/m
    double tau = std::stof(pt.get<std::string>("Parameters.tau"));        // Pulse duration in fs

    // Parameter sweep
    double xmin = std::stof(pt.get<std::string>("Sweep.xmin")); // Minimum frequency in sweep
    double xmax = std::stof(pt.get<std::string>("Sweep.xmax")); // Minimum frequency in sweep
    int x_elem = std::stoi(pt.get<std::string>("Sweep.xnum"));

    // Parameter sweep
    double ymin = std::stof(pt.get<std::string>("Sweep.ymin")); // Minimum frequency in sweep
    double ymax = std::stof(pt.get<std::string>("Sweep.ymax")); // Minimum frequency in sweep
    int y_elem = std::stoi(pt.get<std::string>("Sweep.ynum"));


    // Meshgrid (vectors of parameters)
    auto xvec = arma::linspace(xmin, xmax, x_elem);
    auto yvec = arma::linspace(ymin, ymax, y_elem);

    // Variables to store probability values
    auto prob_mat = arma::mat(y_elem,x_elem);

    // Initial value of time, integration time and interval

    double tinit = -4.0*tau;
    double inttime = 8.0*tau;
    double dt = inttime/100;

    // Variables for loops
    unsigned id, id2;

    // Loop for each K value and compute probability
    for (id=0; id < y_elem; id++)
    {
        for (id2=0; id2 < x_elem; id2++)
        {

            // Initialize tight-binding model
            tight_binding_exp model(xvec[id2],yvec[id],tau,E0);

            // Prepare initial states (negative and positive eigenstates)
            state_type psi = EigenState(model, xvec[id2], yvec[id], -1.0);
            state_type eigen_p = EigenState(model, xvec[id2], yvec[id], 1.0);

            // Integrate
            size_t steps = boost::numeric::odeint::integrate( model,
                           psi, tinit, inttime, dt );

            // Calculate the projection
            std::complex<double> projection = Projection(eigen_p,psi);

            prob_mat(id,id2) = std::abs(projection)*std::abs(projection);


        }

    }

    // Save ODE integration data
    prob_mat.save("probability.dat", arma::raw_ascii);

    // Save parameters vector
    xvec.save("xvec.dat", arma::raw_ascii);
    yvec.save("yvec.dat", arma::raw_ascii);

    // For post-processing: time evolution of vector potential after ODE integration
    int numtimes = ceil(inttime/dt);
    auto times = arma::linspace(tinit,tinit+inttime,numtimes);

    // Initialize output file
    std::ofstream outfile("potential.dat");
    outfile << std::scientific << std::setprecision(10);

    // Create tight binding object
    tight_binding_exp model(0.0,0.0,tau,E0);

    for (size_t i=0; i < numtimes; i++)
    {
        outfile << times(i) << " "
                << model.Gx(times(i)) << " "
                << model.Gy(times(i)) << std::endl;
    }


    return 0;

}
