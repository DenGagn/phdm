/*!
* \file GrapheneSine.cpp
*
* \brief Tight-binding model for graphene photo-excitation
* Solves the ODE system for a wave with a slowly varying envelope
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

typedef std::vector< double > state_type; // Type of container used to hold the state vector
using namespace std::complex_literals; // Complex numbers


#ifndef I_INCLUDE
#define I_INCLUDE
/// Imaginary unit
arma::cx_double I(0.0,1.0);
#endif

#include "GrapheneSine.hpp"

/// Main function
int main(int argc, char *argv[])
{
    // Parse parameter file
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("GrapheneSine.ini", pt);

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

    // Loop for each K va1ilue and compute probability
    for (id=0; id < y_elem; id++)
    {
        for (id2=0; id2 < x_elem; id2++)
        {

            // Prepare initial state
            double Re_Gamma = 1.0 + 2.0*cos(0.5*sqrt(3.0)*xvec[id2])*cos(0.5*3.0*yvec[id]);
            double Im_Gamma = 2.0*cos(0.5*sqrt(3.0)*xvec[id2])*sin(0.5*3.0*yvec[id]);
            double angle_Gamma = atan2(Im_Gamma,Re_Gamma);

            state_type init_psi = {0.0,0.0,0.0,0.0};
            state_type fina_psi = {0.0,0.0,0.0,0.0};

            // Assign values of initial state (the ket)
            init_psi[0] = sqrt(0.5);
            init_psi[1] = 0.0;
            init_psi[2] = -sqrt(0.5)*(cos(angle_Gamma));
            init_psi[3] = -sqrt(0.5)*(sin(angle_Gamma));

            // Assign values of final state (the bra)
            fina_psi[0] = sqrt(0.5);
            fina_psi[1] = 0.0;
            fina_psi[2] = sqrt(0.5)*(cos(angle_Gamma));
            fina_psi[3] = sqrt(0.5)*(sin(angle_Gamma));

            // Initialize tight-binding model
            tight_binding_sine tb(xvec[id2],yvec[id],omega,a,E0);

            // Integrate
            size_t steps = boost::numeric::odeint::integrate( tb,
                           init_psi, tinit, inttime, dt );

            // Calculate the projection
            std::complex<double> projection =
                (fina_psi[0] - 1i*fina_psi[1])*(init_psi[0] + 1i*init_psi[1]) +
                (fina_psi[2] - 1i*fina_psi[3])*(init_psi[2] + 1i*init_psi[3]);

            prob_mat(id,id2) = std::abs(projection)*std::abs(projection);


        }

    }

    // Save data
    prob_mat.save("probability.dat", arma::raw_ascii);

    // Save parameters vector
    xvec.save("xvec.dat", arma::raw_ascii);
    yvec.save("yvec.dat", arma::raw_ascii);

    return 0;

}
