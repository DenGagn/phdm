#pragma once

#include <armadillo>

using namespace std::complex_literals; // Complex numbers
typedef std::vector< double > state_type; // Type of container used to hold the state vector

namespace phosphorene {

/*!
* \class base_sine
*
* \brief Base class for a few-cycle sine wave with variable duration and
* amplitude, but zero carrier-envelope phase (version for phosphorene)
*
* \author Author: D. Gagnon <denisg6@hotmail.com>
*/
class base_sine {

    // Parameters of the Hamiltonian
    double m_kx,m_ky,m_omega,m_a,m_E0,m_ellip;

private:

    // Physical constants definitions
    double hbar = 6.5821195140e-16; // Planck constant over 2 pi in eV s

public:

    /// Constructor taking as input the parameters of the Hamiltonian
    /// @param kx float, x-component of the momentum
    /// @param ky float, y-component of the momentum
    /// @param omega float, angular frequency of the field in rad/s
    /// @param a float, envelope frequency (normalized)
    /// @param E0 float, electric field peak value in V/m
    /// @param ellip float, ellipticity of the field in the circularly polarized case
    base_sine(double kx, double ky, double omega, double a, double E0, double ellip=0.0 )
        : m_kx(kx)
        , m_ky(ky)
        , m_omega(omega)
        , m_a(a)
        , m_E0(E0)
        , m_ellip(ellip) { }

    /// Real part of gamma factor appearing in the Hamiltonian
    virtual std::complex<double> Gamma0 (double, double) = 0;

    /// Imaginary part of gamma factor appearing in the Hamiltonian
    virtual std::complex<double> Gamma1 (double, double) = 0;

    /// Overload of operator() for ODE integration
    void operator() ( const state_type &z, state_type &dzdt, const double t)
    {

        // Conversion factor
        double factor = 1.0e-10*m_E0/(hbar*m_omega); // Multiply by a characteristic length of 1 A

        // Peierls substitution
        double kx_field = m_kx + factor*Gx(t);
        double ky_field = m_ky + m_ellip*factor*Gy(t);

        // Compute gamma factors
        auto Gamma0_field = Gamma0(kx_field, ky_field);
        auto Gamma1_field = Gamma1(kx_field, ky_field);

        // The ODE system
        std::complex<double> z0 = z[0] + 1i*z[1];
        std::complex<double> z1 = z[2] + 1i*z[3];

        std::complex<double> dzdt0 = -1i*(1.0/hbar)*(Gamma0_field*z0 + Gamma1_field*z1);
        std::complex<double> dzdt1 = -1i*(1.0/hbar)*(std::conj(Gamma1_field)*z0 + Gamma0_field*z1);

        dzdt[0] = std::real(dzdt0);
        dzdt[1] = std::imag(dzdt0);
        dzdt[2] = std::real(dzdt1);
        dzdt[3] = std::imag(dzdt1);

    }

    /// Normalized vector potential, x-component
    double Gx (double t)
    {

        return 0.25*m_omega*(sin((2.0*m_a-m_omega)*t)/(2.0*m_a-m_omega)
                             + sin((2.0*m_a+m_omega)*t)/(2.0*m_a+m_omega)
                             - 2.0*sin(m_omega*t)/m_omega);

//         // Version used in 10.1364/JOSAB.35.000958
//         return sin(m_omega*t)*sin(m_a*t)*sin(m_a*t);
    }

    /// Normalized vector potential, y-component
    double Gy (double t)
    {
//         return cos(m_omega*t)*sin(m_a*t)*sin(m_a*t);
        return 0.0;
    }

};

/*!
* \class tight_binding_sine
*
* \brief Tight-binding model for phosphorene (RHS of ODE system)
*
* The excitation
* is a few-cycle sine wave with variable duration and amplitude, but zero
* carrier-envelope phase
*
*/
class tight_binding_sine : public base_sine {

private:
    // Tight-binding constant definitions (see http://arxiv.org/abs/1710.05808)
    double eta = 0.58; // eV A^2
    double nu  = 1.01; // eV A^2
    double gx =  3.93; // eV A^2
    double gy =  3.83; // eV A^2
    double d0 = -0.42; // eV
    double d1 =  0.76; // eV
    double chi = 5.25; // eV A

public:
    /// Constructor taking as input the parameters of the  Hamiltonian (tight-binding specialization)
    /// @param kx float, x-component of the momentum
    /// @param ky float, y-component of the momentum
    /// @param omega float, angular frequency of the field in rad/s
    /// @param a float, envelope frequency (normalized)
    /// @param E0 float, electric field peak value in V/m
    /// @param ellip float, ellipticity of the field in the circularly polarized case
    tight_binding_sine(double kx, double ky, double omega, double a, double E0, double ellip=0.0 ) :
        base_sine(kx, ky, omega, a, E0, ellip ) {}

    /// Gamma0 factor appearing in the tight-binding Hamiltonian
    std::complex<double> Gamma0(double kx, double ky)
    {
        return d0 + eta*kx*kx + nu*ky*ky;
    }

    /// Gamma1 factor appearing in the tight-binding Hamiltonian
    std::complex<double> Gamma1(double kx, double ky)
    {
        return d1 + 1i*chi*ky + gx*kx*kx + gy*ky*ky;
    }


};

/*!
* \class kp_sine
*
* \brief k.p. model for phosphorene (RHS of ODE system)
*
* NOT YET IMPLEMENTED, TO BE ADDED LATER
*
* The excitation
* is a few-cycle sine wave with variable duration and amplitude, but zero
* carrier-envelope phase
*
*/
class kp_sine : public base_sine {

private:
    // Tight-binding constant definitions (see http://arxiv.org/abs/1710.05808)
    // NOT YET IMPLEMENTED, TO BE ADDED LATER

public:
    /// Constructor taking as input the parameters of the  Hamiltonian (tight-binding specialization)
    /// @param kx float, x-component of the momentum
    /// @param ky float, y-component of the momentum
    /// @param omega float, angular frequency of the field in rad/s
    /// @param a float, envelope frequency (normalized)
    /// @param E0 float, electric field peak value in V/m
    /// @param ellip float, ellipticity of the field in the circularly polarized case
    kp_sine(double kx, double ky, double omega, double a, double E0, double ellip=0.0 ) :
        base_sine(kx, ky, omega, a, E0, ellip ) {}

    /// Gamma0 factor appearing in the tight-binding Hamiltonian
    std::complex<double> Gamma0(double kx, double ky)
    {
        return 0.0;
    }

    /// Gamma1 factor appearing in the tight-binding Hamiltonian
    std::complex<double> Gamma1(double kx, double ky)
    {
        return 0.0;
    }


};

}
