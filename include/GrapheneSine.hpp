/// Tight-binding model for graphene (RHS of ODE system)
class tight_binding {

    double m_kx,m_ky,m_omega,m_a,m_E0;

private:
    double hbar = 6.5821195140e-16;                             // Planck constant over 2 pi in eV s
    double lat_constant = 2.4e-10/sqrt(3.0);                    // Lattice constant in meters
    double v_Fermi = 0.003646*2.9979245800e+08;                 // Fermi velocity in graphene
    double time_constant = 2.0*v_Fermi / (3.0*lat_constant);    // Time constant in the differential equation (in reciprocal seconds)

public:
    tight_binding(double kx, double ky, double omega, double a, double E0 )
        : m_kx(kx)
        , m_ky(ky)
        , m_omega(omega)
        , m_a(a)
        , m_E0(E0) { }

    void operator() ( const state_type &z, state_type &dzdt, const double t)
    {

        // Compute the normalized vector potential (polarization in x direction)
        double Gx = 0.25*m_omega*(sin((2.0*m_a-m_omega)*t)/(2.0*m_a-m_omega)
                                  + sin((2.0*m_a+m_omega)*t)/(2.0*m_a+m_omega)
                                  - 2.0*sin(m_omega*t)/m_omega);

        // Conversion factor
        double factor = lat_constant*m_E0/(hbar*m_omega);

        // Peierls substitution
        double kx_field = m_kx + factor*Gx;
        double ky_field = m_ky + 0.0;

        // Compute gamma factors
        double Re_Gamma = 1.0 + 2.0*cos(0.5*sqrt(3.0)*kx_field)*cos(0.5*3.0*ky_field);
        double Im_Gamma = 2.0*cos(0.5*sqrt(3.0)*kx_field)*sin(0.5*3.0*ky_field);

        // The ODE system
        dzdt[0] = time_constant*( Re_Gamma*z[3] - Im_Gamma*z[2]); // RePsi_A
        dzdt[1] = time_constant*(-Re_Gamma*z[2] - Im_Gamma*z[3]); // ImPsi_A
        dzdt[2] = time_constant*( Re_Gamma*z[1] + Im_Gamma*z[0]); // RePsi_B
        dzdt[3] = time_constant*(-Re_Gamma*z[0] + Im_Gamma*z[1]); // ImPsi_B

    }
};
