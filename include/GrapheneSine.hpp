#pragma once

typedef std::vector< double > state_type; // Type of container used to hold the state vector

/*!
* \class base_sine
*
* \brief Base class for a few-cycle sine wave with variable duration and
* amplitude, but zero carrier-envelope phase
*
* \author Author: D. Gagnon <denisg6@hotmail.com>
*/

namespace graphene {

class base_sine {

    // Parameters of the Hamiltonian
    double m_kx,m_ky,m_omega,m_a,m_E0,m_ellip;

private:

    // Physical constants definitions
    double hbar = 6.5821195140e-16;                             // Planck constant over 2 pi in eV s
    double lat_constant = 2.4e-10/sqrt(3.0);                    // Lattice constant in meters
    double v_Fermi = 0.003646*2.9979245800e+08;                 // Fermi velocity in graphene
    double time_constant = -2.0*v_Fermi / (3.0*lat_constant);    // Time constant in the differential equation (in reciprocal seconds)

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
    virtual double Re_Gamma (double, double) = 0;

    /// Imaginary part of gamma factor appearing in the Hamiltonian
    virtual double Im_Gamma (double, double) = 0;

    /// Overload of operator() for ODE integration
    void operator() ( const state_type &z, state_type &dzdt, const double t)
    {

        // Conversion factor
        double factor = lat_constant*m_E0/(hbar*m_omega);

        // Peierls substitution
        double kx_field = m_kx + factor*Gx(t);
        double ky_field = m_ky + m_ellip*factor*Gy(t);

        // Compute gamma factors
        double Re_Gamma_t = Re_Gamma(kx_field, ky_field);
        double Im_Gamma_t = Im_Gamma(kx_field, ky_field);

        // The ODE system
        dzdt[0] = time_constant*( Re_Gamma_t*z[3] - Im_Gamma_t*z[2]); // RePsi_A
        dzdt[1] = time_constant*(-Re_Gamma_t*z[2] - Im_Gamma_t*z[3]); // ImPsi_A
        dzdt[2] = time_constant*( Re_Gamma_t*z[1] + Im_Gamma_t*z[0]); // RePsi_B
        dzdt[3] = time_constant*(-Re_Gamma_t*z[0] + Im_Gamma_t*z[1]); // ImPsi_B

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
* \brief Tight-binding model for graphene (RHS of ODE system)
*
* The excitation
* is a few-cycle sine wave with variable duration and amplitude, but zero
* carrier-envelope phase
*
*/
class tight_binding_sine : public base_sine {

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


};

/*!
* \class dirac_sine
*
* \brief Dirac model for graphene (RHS of ODE system)
*
* The excitation
* is a few-cycle sine wave with variable duration and amplitude, but zero
* carrier-envelope phase
*
*/
class dirac_sine : public base_sine {

public:
    /// Constructor taking as input the parameters of the Hamiltonian (Dirac specialization)
    /// @param kx float, x-component of the momentum
    /// @param ky float, y-component of the momentum
    /// @param omega float, angular frequency of the field in rad/s
    /// @param a float, envelope frequency (normalized)
    /// @param E0 float, electric field peak value in V/m
    /// @param ellip float, ellipticity of the field in the circularly polarized case
    dirac_sine(double kx, double ky, double omega, double a, double E0, double ellip=0.0 ) :
        base_sine(kx, ky, omega, a, E0, ellip ) {}

    /// Real part of gamma factor appearing in the Dirac Hamiltonian
    double Re_Gamma (double kx, double)
    {
        return -1.5*kx;
    }

    /// Imaginary part of gamma factor appearing in the Dirac Hamiltonian
    double Im_Gamma (double, double ky)
    {
        return -1.5*ky;
    }


};

}
