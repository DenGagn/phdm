/*!
* \file PhosphoreneSine.cpp
*
* \brief Tight-binding model for phosphorene photo-excitation
* Solves the ODE system for a wave with a slowly varying envelope (sine)
*
* The excitation
* is a few-cycle sine wave with variable duration and amplitude, but zero
* carrier-envelope phase
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

#include "PhosphoreneSine.hpp"
#include "UtilsComplex.hpp"

/// Main function
int main(int argc, char *argv[])
{
    // Parse parameter file
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("PhosphoreneSine.ini", pt);

    // Problem parameters
    double freq  = std::stof(pt.get<std::string>("Parameters.frequency"));  // Angular frequency in Hz
    double E0 = std::stof(pt.get<std::string>("Parameters.E0"));            // Peak electric field in V/m
    double afreq = std::stof(pt.get<std::string>("Parameters.afreq"));         // Envelope frequency, in units of freq

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

    // Process input variables (angular frequency units)
    double omega = 2.0*M_PI*freq;     // Angular frequency
    double a = omega*afreq;         // Envelope frequency

    // Initial value of time, integration time and interval
    double tinit = 0.0;
    double inttime = M_PI/a;
    double dt = inttime/100.0;

    // Variables for loops
    unsigned id, id2;

    // Loop for each K value and compute probability
    for (id=0; id < y_elem; id++)
    {
        for (id2=0; id2 < x_elem; id2++)
        {

            // Initialize tight-binding model
            tight_binding_sine model(xvec[id2],yvec[id],omega,a,E0);

            // Prepare initial states (negative and positive eigenstates)
            state_type psi = EigenState(model, xvec[id2], yvec[id], false);
            state_type eigen_p = EigenState(model, xvec[id2], yvec[id], true);

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
    tight_binding_sine model(0.0,0.0,omega,a,E0);

    for (size_t i=0; i < numtimes; i++)
    {
        outfile << times(i) << " "
                << model.Gx(times(i)) << " "
                << model.Gy(times(i)) << std::endl;
    }


    return 0;

}
