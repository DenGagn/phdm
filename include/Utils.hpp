/*!
* \file Utils.hpp
*
* \brief Various utility functions for phdm routines
*
* \author Author: D. Gagnon <denisg6@hotmail.com>
*/
// Include some headers
#include <cmath>
#include <complex>

using namespace std::complex_literals; // Complex numbers
typedef std::vector< double > state_type; // Type of container used to hold the state vector

/// Function computing an eigenstate for a given model
/// @param model T, the Hamiltonian model
/// @param kx float, x-component of the momentum
/// @param ky float, y-component of the momentum
/// @param sign float, sign in front of the last 2 components (should be -1.0 for negative energy)
template <typename T>
state_type EigenState(T model, double kx, double ky, double sign=1.0)
{
        // Prepare initial state (calculate gamma factor and its angles)
        double Re_Gamma = model.Re_Gamma(kx, ky);
        double Im_Gamma = model.Im_Gamma(kx, ky);
        double angle_Gamma = atan2(Im_Gamma,Re_Gamma);

        state_type eigen = {0.0,0.0,0.0,0.0};

        // Assign values of initial state (the ket)
        eigen[0] = sqrt(0.5);
        // eigen[1] = 0.0;
        eigen[2] = sign*sqrt(0.5)*(cos(angle_Gamma));
        eigen[3] = sign*sqrt(0.5)*(sin(angle_Gamma));
        
        return eigen;
}

/// Function computing the projection of two eigenvectors (returns a complex number)
template <typename T>
std::complex<T> Projection(std::vector<T> vec1, std::vector<T> vec2)
{
        return (vec1[0] - 1i*vec1[1])*(vec2[0] + 1i*vec2[1]) +
                (vec1[2] - 1i*vec1[3])*(vec2[2] + 1i*vec2[3]);
}
