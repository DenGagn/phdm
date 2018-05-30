#pragma once

typedef std::vector< double > state_type; // Type of container used to hold the state vector

/*!
* \class base_exp
*
* \brief Base class for a slowly varying envelope with exponential shape
* See e.g. http://link.aps.org/doi/10.1103/PhysRevB.91.045439
*
* \author Author: D. Gagnon <denisg6@hotmail.com>
*/
class base_exp {

    // Parameters of the Hamiltonian
    double m_kx,m_ky,m_tau,m_E0;

private:

    // Physical constants definitions
    double hbar = 6.5821195140e-16;                             // Planck constant over 2 pi in eV s
    double lat_constant = 2.4e-10/sqrt(3.0);                    // Lattice constant in meters
    double v_Fermi = 0.003646*2.9979245800e+08;                 // Fermi velocity in graphene
    double time_constant = -2.0*v_Fermi / (3.0*lat_constant);   // Time constant in the differential equation (in reciprocal seconds)

public:
    /// Constructor taking as input the parameters of the Hamiltonian
    /// @param kx float, x-component of the momentum
    /// @param ky float, y-component of the momentum
    /// @param tau float, duration of the exponential pulse in fs
    /// @param E0 float, electric field peak value in V/m
    base_exp(double kx, double ky, double tau, double E0 )
        : m_kx(kx)
        , m_ky(ky)
        , m_tau(tau)
        , m_E0(E0) { }

    /// Real part of gamma factor appearing in the Hamiltonian
    virtual double Re_Gamma (double, double) = 0;

    /// Imaginary part of gamma factor appearing in the Hamiltonian
    virtual double Im_Gamma (double, double) = 0;

    /// Overload of operator() for ODE integration
    void operator() ( const state_type &z, state_type &dzdt, const double t)
    {

        // Conversion factor
        double factor = lat_constant*m_E0*m_tau/hbar;

        // Peierls substitution
        double kx_field = m_kx + factor*Gx(t);
        double ky_field = m_ky + factor*Gy(t);

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
        return 0.0;
    }

    /// Normalized vector potential, y-component
    double Gy (double t)
    {
        return -exp(-(t/m_tau)*(t/m_tau))*(t/m_tau);
    }


};

/*!
* \class tight_binding_exp
*
* \brief Tight-binding model for graphene (RHS of ODE system)
*
* The excitation
* is a slowly varying envelope with exponential shape
* See e.g. http://link.aps.org/doi/10.1103/PhysRevB.91.045439
*
* \author Author: D. Gagnon <denisg6@hotmail.com>
*/
class tight_binding_exp : public base_exp {

public:
    /// Constructor taking as input the parameters of the Hamiltonian (tight-binding specialization)
    /// @param kx float, x-component of the momentum
    /// @param ky float, y-component of the momentum
    /// @param tau float, duration of the exponential pulse in fs
    /// @param E0 float, electric field peak value in V/m
    tight_binding_exp(double kx, double ky, double tau, double E0 ) :
        base_exp(kx, ky, tau, E0 ) {}

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
* \class dirac_exp
*
* \brief Dirac model for graphene (RHS of ODE system)
*
* The excitation
* is a slowly varying envelope with exponential shape
* See e.g. http://link.aps.org/doi/10.1103/PhysRevB.91.045439
*
* \author Author: D. Gagnon <denisg6@hotmail.com>
*/
class dirac_exp : public base_exp {

public:
    /// Constructor taking as input the parameters of the Hamiltonian (Dirac specialization)
    /// @param kx float, x-component of the momentum
    /// @param ky float, y-component of the momentum
    /// @param tau float, duration of the exponential pulse in fs
    /// @param E0 float, electric field peak value in V/m
    dirac_exp(double kx, double ky, double tau, double E0 ) :
        base_exp(kx, ky, tau, E0 ) {}

    /// Real part of gamma factor appearing in the Dirac Hamiltonian
    double Re_Gamma (double kx, double ky)
    {
        return -1.5*kx;
    }

    /// Imaginary part of gamma factor appearing in the Dirac Hamiltonian
    double Im_Gamma (double kx, double ky)
    {
        return -1.5*ky;
    }
};
