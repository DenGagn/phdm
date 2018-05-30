#pragma once
/*!
* \file GrapheneFloquet.hpp
*
* \brief Functions calculating gamma factors (Fourier coefficients) for Floquet Hamiltonian
*
* \author Author: D. Gagnon <denisg6@hotmail.com>
*/
// Include some headers
#include <complex>
#include <armadillo>
#include <boost/numeric/odeint.hpp>

using namespace std::complex_literals; // Complex literals

typedef std::complex<double> c_state_type; // Type of container used to hold the state vector
typedef boost::numeric::odeint::runge_kutta_cash_karp54< c_state_type > error_stepper_type; // Error stepper for odeint

/// Real part of gamma factor appearing in the tight-binding Hamiltonian
double Re_Gamma (double kx, double ky)
{
    return 1.0 + 2.0*cos(0.5*sqrt(3.0)*kx)*cos(0.5*3.0*ky);
}

/// Imaginary part of gamma factor appearing in the tight-binding Hamiltonian
double Im_Gamma (double kx, double ky)
{
    return 2.0*cos(0.5*sqrt(3.0)*kx)*sin(0.5*3.0*ky);
}

/*!
* \class gamma_integrand
*
* \brief Base class for the integrand used in gamma factor calculations
*
* \author Author: D. Gagnon <denisg6@hotmail.com>
*/
class gamma_integrand {

    // Parameters of the Hamiltonian
    double m_kx,m_ky,m_omega,m_E0,m_index;

private:

    // Physical constants definitions
    double hbar = 6.5821195140e-16;                             // Planck constant over 2 pi in eV s
    double lat_constant = 2.4e-10/sqrt(3.0);                    // Lattice constant in meters
    double v_Fermi = 0.003646*2.9979245800e+08;                 // Fermi velocity in graphene

public:

    /// Constructor taking as input the parameters of the Hamiltonian
    /// @param kx float, x-component of the momentum
    /// @param ky float, y-component of the momentum
    /// @param omega float, angular frequency of the field in rad/s
    /// @param E0 float, electric field peak value in V/m
    /// @param index int, index of Fourier coefficients
    gamma_integrand(double kx, double ky, double omega, double E0, int index)
        : m_kx(kx)
        , m_ky(ky)
        , m_omega(omega)
        , m_E0(E0)
        , m_index(index) { }

    /// Overload of operator() for ODE integration
    void operator() ( const c_state_type &z, c_state_type &dzdt, const double t)
    {

        // Conversion factor
        double factor = lat_constant*m_E0/(hbar*m_omega);

        // Peierls substitution
        double kx_field = m_kx + factor*std::sin(t);
        double ky_field = m_ky;

        // Compute gamma factors
        double Re_Gamma_t = Re_Gamma(kx_field, ky_field);
        double Im_Gamma_t = Im_Gamma(kx_field, ky_field);

        // Rhs of the ODE system
        dzdt = (Re_Gamma_t + 1i*Im_Gamma_t)*std::exp(-1i*m_index*t);

    }

    /// Function updating the value of the index
    void SetIndex(int m)
    {
        this->m_index = m;
    }

};

/// Function returning the Fourier coefficients to be used in Floquet Hamiltonian
std::vector< c_state_type > GammaValues(double kx, double ky, double omega, double E0, int max_index)
{
    std::vector< c_state_type > Gamma(2*max_index + 1, 0.0);

    for (int vec_index = 0; vec_index < 2*max_index + 1; vec_index ++)
    {

        // State variable to solve result
        c_state_type z(0.0,0.0);

        // Define integrand
        gamma_integrand integrand(kx,ky,omega,E0, vec_index - max_index);

        // Adaptive integration
        boost::numeric::odeint::integrate_adaptive(
            boost::numeric::odeint::make_controlled< error_stepper_type >( 1.0e-10, 1.0e-6 ),
            integrand,
            z, 0.0, 2.0*M_PI, 2.0*M_PI/10000.0 );

        // Assign value to "Gamma" vector
        Gamma[vec_index] = (0.5/M_PI)*z;

    }

    return Gamma;
}

/// Compute Floquet eigen-energies and probabilities
arma::vec QuasiEnergies(double kx, double ky, double omega, double E0,
                        int blocks, double &prob, double &prob_sigma)
{

    // Physical constants definitions
    double hbar = 6.5821195140e-16;                             // Planck constant over 2 pi in eV s
    double lat_constant = 2.4e-10/sqrt(3.0);                    // Lattice constant in meters
    double v_Fermi = 0.003646*2.9979245800e+08;                 // Fermi velocity in graphene

    double freq = hbar*omega;                        // Frequency in eV
    double tb = - 2*hbar/(3.0*v_Fermi*lat_constant); // Tight-binding energy of graphene

    // Total number of blocks including
    // zeroth block
    int totalblocks = 2*blocks + 1;

    // Pre-calculate Fourier coefficients of Hamiltonian
    auto Gamma = GammaValues(kx,ky,omega,E0,2*blocks);

    // Initialize Tmatrix which is to be filled by [totalblocks] 2 x 2 matrices
    arma::Mat<double> zeromat(2*totalblocks,2*totalblocks, arma::fill::zeros);
    arma::cx_mat Ham(zeromat,zeromat); // Hamiltonian

    // FILL MAIN DIAGONAL OF HAMILTONIAN
    // Initialize index m
    int m = -blocks;

    for (int j_= 0; j_ < 2*totalblocks; j_ = j_ + 2)
    {
        Ham.diag(0)[j_] = m*freq;
        Ham.diag(0)[j_ + 1] = m*freq;

        // Increment m
        m++;
    }

    // FILL UPPER AND LOWER DIAGONAL OF HAMILTONIAN

    int nm = 0; // Initialize index (n-m)
    int index_shift = 2*blocks; // Index shift (for gamma vector)

    for (int i_= 1; i_ < 2*totalblocks; i_ = i_ + 2) // Loop on diagonals
    {
        int index0 = nm;
        int index1 = nm + 1;

        //std::cout << i_ << std::endl;

        for (int j_ = 0; j_ < 2*totalblocks - i_; j_ ++) // Loop on elements of diagonals
        {

            //std::cout << j_ << std::endl;

            if (j_ % 2 == 0) // If j_ is even do
            {
                Ham.diag( i_)[j_] = tb*std::conj(Gamma[-index0 + index_shift]);
                Ham.diag(-i_)[j_] = tb*Gamma[index0 + index_shift];

                //std::cout << "Even j" << std::endl;
                //std::cout << Ham.diag( i_)[j_] << std::endl;
                //std::cout << Ham.diag( -i_)[j_] << std::endl;

            }
            else // If j_ is odd do
            {
                Ham.diag( i_)[j_] = tb*Gamma[-index1 + index_shift];
                Ham.diag(-i_)[j_] = tb*std::conj(Gamma[index1 + index_shift]);

                //std::cout << "Odd j" << std::endl;
                //std::cout << Ham.diag( i_)[j_] << std::endl;
                //std::cout << Ham.diag( -i_)[j_] << std::endl;
            }
        }

        nm++; // Increment (n-m) as we shift to the next non-zero diagonal

    }

    // std::cout << Ham << std::endl;

    // arma::mat realpart = arma::real(Ham);
    // arma::mat imagpart = arma::imag(Ham);

    // realpart.save("Hamreal.dat",arma::raw_ascii);
    // imagpart.save("Hamimag.dat",arma::raw_ascii);

    // Initialize eigenvalue vector and eigenvector matrix and
    arma::vec eigval;
    arma::cx_mat eigvec;

    // Compute eigenvalues
    bool status = arma::eig_sym(eigval, eigvec, Ham);
    if (status == false)
    {
        std::cout << "Eigenvalue decomposition failed" << std::endl;
    }

    std::cout << eigvec << std::endl;

    // COMPUTE TRANSITION PROBABILITY

    // Variables to be used in loop
    prob = 0.0;
    prob_sigma = 0.0; // Passed by reference, will change
    std::complex<double> alphazero, betazero; // Ground state amplitude
    std::complex<double> alphan, betan; // Excited state amplitude
    std::complex<double> pos; // Positive energy state amplitude
    std::complex<double> neg; // Negative energy state amplitude
    arma::cx_vec the_eigenvec;

    // Compute "no-field" phase factor
    double angle = std::atan2(Im_Gamma(kx,ky), Re_Gamma(kx,ky) );
    std::complex<double> phase_factor = std::exp(1i*angle);

    // Loop on every eigenvector
    for (size_t i_= 0; i_ < eigvec.n_cols; i_ ++)
    {
        the_eigenvec = eigvec.col(i_); // Slice eigenvectors

        alphazero = the_eigenvec[the_eigenvec.n_elem/2 - 1];
        betazero = the_eigenvec[the_eigenvec.n_elem/2];

        // Negative energy state
        neg = 0.5*std::sqrt(2.0)*(alphazero - phase_factor*betazero);

        // Loop on elements
        for (size_t j_= 0; j_ < the_eigenvec.n_elem; j_ ++)
        {

            if (j_ % 2 != 0) // If j is odd do, else do nothing
            {
                alphan = the_eigenvec[j_ - 1];
                betan  = the_eigenvec[j_];

                // Positive energy state
                pos = 0.5*std::sqrt(2.0)*(alphan + phase_factor*betan);

                // Fermi's rule
                prob += std::norm(std::conj(pos)*neg)*std::norm(std::conj(pos)*neg); // Transition between eigenstates
                prob_sigma += std::norm(std::conj(betan)*alphazero)*std::norm(std::conj(betan)*alphazero); // Transition between sigma_z eigenstates
            }
        }
    }


    return eigval;



}
