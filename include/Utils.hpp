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

/// Function computing an eigenstate for a given model
/// with a complex 2 x 2 Hamiltonian
/// @param model T, the Hamiltonian model
/// @param kx float, x-component of the momentum
/// @param ky float, y-component of the momentum
/// @param bool pos, true if positive eigenstate, false otherwise
template <typename T>
state_type EigenStateGeneral(T model, double kx, double ky, bool pos=true)
{
    // Calculate matrix elements of the Hamiltonian
    std::complex<double> Gamma0 = model.Gamma0(kx,ky);
    std::complex<double> Gamma1 = model.Gamma1(kx,ky);

    // Define Hamiltonian matrix
    arma::cx_mat Ham;
    Ham.resize(2,2);

    Ham(0,0) = Gamma0;
    Ham(0,1) = Gamma1;
    Ham(1,0) = std::conj(Gamma1);
    Ham(1,1) = Gamma0;

    // Initialize eigenvalue vector and eigenvector matrix and
    arma::vec eigval;
    arma::cx_mat eigvecs;

    // Compute eigenvalues
    bool status = arma::eig_sym(eigval, eigvecs, Ham);
    if (status == false)
    {
        std::cout << "Eigenvalue decomposition failed" << std::endl;
    }

    // Control loop for positive/negative energy eigenstate
    arma::cx_vec the_eigenvector;

    if (pos == true)
    {
        the_eigenvector = eigvecs.col(0);
    }
    else
    {
        the_eigenvector = eigvecs.col(1);
    }

    // Normalisation constant
    double norm_const
        = 1.0/sqrt(std::norm(the_eigenvector(0)) + std::norm(the_eigenvector(1)));

    //std::cout << the_eigenvector(0) << std::endl;


    // Assign values of the eigenstate
    state_type eigen = {0.0,0.0,0.0,0.0};

    eigen[0] = norm_const*std::real(the_eigenvector(0));
    eigen[1] = norm_const*std::imag(the_eigenvector(0));
    eigen[2] = norm_const*std::real(the_eigenvector(1));
    eigen[3] = norm_const*std::imag(the_eigenvector(1));

    return eigen;
}

/// Function computing the projection of two eigenvectors (returns a complex number)
template <typename T>
std::complex<T> Projection(std::vector<T> vec1, std::vector<T> vec2)
{
    return (vec1[0] - 1i*vec1[1])*(vec2[0] + 1i*vec2[1]) +
           (vec1[2] - 1i*vec1[3])*(vec2[2] + 1i*vec2[3]);
}
