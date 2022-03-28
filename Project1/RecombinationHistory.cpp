#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================

RecombinationHistory::RecombinationHistory(
    BackgroundCosmology* cosmo,
    double Yp) :
    cosmo(cosmo),
    Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve() {

    // Compute and spline Xe, ne
    solve_number_density_electrons();

    // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
    solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons() {
    Utils::StartTiming("Xe");

    //=============================================================================
    // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
    //=============================================================================
    Vector x_array = Utils::linspace(x_start,x_end, npts_rec_arrays);
    Vector Xe_arr(npts_rec_arrays);
    Vector ne_arr(npts_rec_arrays);
    Vector Xe_saha(npts_rec_arrays);
    Vector ne_saha(npts_rec_arrays);

    for (int i = 0; i < npts_rec_arrays; i++) {
        auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

        Xe_saha[i] = Xe_ne_data.first;
        ne_saha[i] = Xe_ne_data.second;
    }
    Xe_Saha_spline.create(x_array, Xe_saha, "saha");
    ne_Saha_spline.create(x_array, ne_saha, "ne_saha");

    // Calculate recombination history
    bool saha_regime = true;
    for (int i = 0; i < npts_rec_arrays; i++) {

        //==============================================================
        // TODO: Get X_e from solving the Saha equation so
        // implement the function electron_fraction_from_saha_equation
        //==============================================================
        auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
        

        // Electron fraction and number density
        const double Xe_current = Xe_ne_data.first;
        const double ne_current = Xe_ne_data.second;
        // Are we still in the Saha regime?
        if (Xe_current < Xe_saha_limit)
            saha_regime = false;

        if (saha_regime) {
            Xe_arr[i] = Xe_current;
            ne_arr[i] = ne_current;

        }
        else {

            auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
            const double Xe_init = Xe_ne_data.first;
            std::cout << "end of saha at x=" << x_array[i] << "\n";
            Vector Xe_ini{ Xe_init};
            Vector Peeb_x_array= Utils::linspace(x_array[i],x_end, npts_rec_arrays-i) ;
            // The Peebles ODE equation
            ODESolver peebles_Xe_ode;
            ODEFunction dXedx = [&](double x, const double* Xe, double* dXedx) {
                return rhs_peebles_ode(x, Xe, dXedx);
            };
            peebles_Xe_ode.solve(dXedx, Peeb_x_array, Xe_ini);
            //=============================================================================
            // TODO: Set up IC, solve the ODE and fetch the result 
            //=============================================================================
            //...
            //...
            auto Peeb_array = peebles_Xe_ode.get_data_by_component(0);
            const double OmegaB = cosmo->get_OmegaB(0);
            const double H0 = cosmo->get_H0();
            const double m_H = Constants.m_H;
            const double G = Constants.G;
            for (int j = i; j < npts_rec_arrays; j++) {
                double a = exp(x_array[j]);
                double nH = OmegaB * (3 * pow(H0, 2) / (m_H * pow(a, 3)* 8 * G * M_PI));
                Xe_arr[j] = Peeb_array[j - i];
                ne_arr[j] = Peeb_array[j - i]*nH;
            }
            break;

        }

    }
    //=============================================================================
    // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
    // functions are working
    //=============================================================================
    //...
    //...

    // Compute time when Xe = 0.5
    for (int k = 0; k < npts_rec_arrays; k++) {
        if (Xe_arr[k] < 0.5) {
            double x = x_array[k - 1];
                std::cout << "Xe is 0.5 at x = " << x << "  and z = " << cosmo->get_z(x) << "\n";
                std::cout << "\n";

                break;
        }
    }

    for (int k = 0; k < npts_rec_arrays; k++) {
        if (Xe_saha[k] < 0.5) {
            double x = x_array[k - 1];
            std::cout << "Xe from Saha equation is 0.5 at x = " << x << "  and z = " << cosmo->get_z(x) << "\n";
            std::cout << "\n";

            break;
        }
    }


    auto log_ne_arr = log(ne_arr);

    Xe_of_x_spline.create(x_array, Xe_arr, "Xe");
    log_ne_of_x_spline.create(x_array, log_ne_arr, "ne");
    std::cout << "Freeze-out abundance of free electrons is Xe(0) = " << Xe_of_x_spline(0) <<"\n";
    std::cout << "Freeze-out abundance of free electrons from the Saha equation is Xe(0) = " << Xe_Saha_spline(0) << "\n";


    Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double, double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const {
    const double a = exp(x);

    // Physical constants
    const double k_b = Constants.k_b;
    const double G = Constants.G;
    const double m_e = Constants.m_e;
    const double hbar = Constants.hbar;
    const double m_H = Constants.m_H;
    const double epsilon_0 = Constants.epsilon_0;
    const double H0_over_h = Constants.H0_over_h;
    // Fetch cosmological parameters
    const double Tb = cosmo->get_TCMB(x);
    const double OmegaB = cosmo->get_OmegaB(0);
    const double H0 = cosmo->get_H0() ;
    double Xe; double ne;

    // Electron fraction and number density
    double nH = OmegaB * (3 * pow(H0, 2) / (m_H * pow(a, 3) * 8 * G * M_PI));;
    double b = pow(m_e*Tb*2.*M_PI*k_b /(pow(hbar,2)),3/2.)*exp(-(epsilon_0)/(k_b*Tb))/nH;
    if (log10(b) > 8) {
        Xe = 1. ;
        ne = Xe*nH;
    }
    else {
        Xe = (-b + sqrt(pow(b,2) + 4 * b)) / 2.;
        ne = Xe * nH;
    }


    return std::pair<double, double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double* Xe, double* dXedx) {

    // Current value of a and X_e
    const double X_e = Xe[0];
    const double a = exp(x);
    // Physical constants in SI units
    const double k_b = Constants.k_b;
    const double G = Constants.G;
    const double c = Constants.c;
    const double m_e = Constants.m_e;
    const double hbar = Constants.hbar;
    const double m_H = Constants.m_H;
    const double sigma_T = Constants.sigma_T;
    const double lambda_2s1s = Constants.lambda_2s1s;
    const double epsilon_0 = Constants.epsilon_0;
    const double alpha = 1 / 137.0359992;
    // Cosmological parameters
    // const double OmegaB      = cosmo->get_OmegaB();
    const double H = cosmo->H_of_x(x);
    const double H0 = cosmo->get_H0();
    const double Tb = cosmo->get_TCMB(x);
    const double OmegaB0 = cosmo->get_OmegaB(0);
    double phi2 = 0.448 * log(epsilon_0 / (k_b*Tb));
    double alpha2 = (64*M_PI*pow(alpha,2))/(sqrt(27*M_PI)*pow(m_e,2))*sqrt(epsilon_0/(k_b*Tb))*phi2*pow(hbar,2)/c;
    double nH = OmegaB0 * (3 * pow(H0, 2) / (m_H * pow(a, 3) * 8 * G * M_PI));
    double n1s = (1 - X_e) * nH;
    double beta2 = alpha2*pow(m_e*Tb*k_b/(2*M_PI) ,3/2.)*exp(-epsilon_0/(4.0*k_b*Tb))/pow(hbar,3.);
    double beta = alpha2 * pow(m_e * Tb *k_b/ (2 * M_PI), 3 / 2.) * exp(-epsilon_0 / (k_b*Tb))/pow(hbar,3.);
    // std::cout << X_e << "\n";
    double lambda_a = H*pow(3*epsilon_0,3)/(pow(8*M_PI,2)*n1s);
    double C_r = (lambda_2s1s + lambda_a) / (lambda_2s1s + lambda_a + beta2);
    // ...

    //=============================================================================
    // TODO: Write the expression for dXedx
    //=============================================================================
    //...
    //...

    dXedx[0] = C_r/H *(beta*(1-X_e) - nH*alpha2*pow(X_e,2 ));

    return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau() {
    Utils::StartTiming("opticaldepth");
    const double c = Constants.c;
    const double sigma_T = Constants.sigma_T;
    // const int npts = 1000;
    Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
    ODEFunction dtaudx = [&](double x, const double* tau, double* dtaudx) {
        dtaudx[0] = dtaudx_of_x(x);

        return GSL_SUCCESS;
    };
    Vector tau_init{ 0.0 };
    ODESolver tau_ODE;
    Vector tau_deriv(npts_rec_arrays);
    tau_ODE.solve(dtaudx,x_array,tau_init);
    auto tau_arr = tau_ODE.get_data_by_component(0);
    //

    for (int j = 0; j < npts_rec_arrays; j++) {
        const double tol = 1E-8;
        tau_arr[j] -= tau_arr[npts_rec_arrays - 1]; 
        if (tau_arr[j] < tol) {
            tau_arr[j] = 0.0;
        }
        tau_deriv[j] = dtaudx_of_x(x_array[j]);
        
    }
    // Find the x and z for decoupling
    for (int k = 0; k < npts_rec_arrays; k++) {
        if (tau_arr[k] < 1.) {
            std::cout << "Tau = 1 at x = " << x_array[k - 1] << "  , z = " << cosmo->get_z(x_array[k-1]) << "\n";
            break;
        }
    }
    Vector saha_init{ 0.0 };
    // Find x and z for decoupling when we only use Saha
    ODEFunction sahadtaudx = [&](double x, const double* tau, double* sahadtaudx) {
        sahadtaudx[0] = saha_dtaudx_of_x(x);

        return GSL_SUCCESS;
    };
    ODESolver saha_tau_ODE;
    saha_tau_ODE.solve(sahadtaudx, x_array, saha_init);
    auto saha_tau_arr = saha_tau_ODE.get_data_by_component(0);

    for (int j = 0; j < npts_rec_arrays; j++) {
        saha_tau_arr[j] -= saha_tau_arr[npts_rec_arrays - 1];
    }
    for (int b = 1; b < npts_rec_arrays; b++) {
        if (saha_tau_arr[b] < 1.) {
            std::cout << "For the Saha equation, Tau = 1 at x = " << x_array[b - 1] << "  , z = " << cosmo->get_z(x_array[b - 1]) << "\n";
            break;
        }
    }

    // Spline results
    tau_of_x_spline.create(x_array, tau_arr, "tau");
    dtaudx_of_x_spline.create(x_array, tau_deriv, "dtau");
    // Calculate g tilde

    Vector g_tilde(npts_rec_arrays) ;
    for (int i = 0; i < npts_rec_arrays; i++) {
        g_tilde[i] = -dtaudx_of_x(x_array[i]) * exp(-tau_of_x(x_array[i]));
    }
    g_tilde_of_x_spline.create(x_array, g_tilde, "g");

    Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const {
    return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const {
    const double c = Constants.c;
    const double sigma_T = Constants.sigma_T;
    double H = cosmo->H_of_x(x);
    return -c * sigma_T * ne_of_x(x) / H;;
}

double RecombinationHistory::saha_dtaudx_of_x(double x) const {
    const double c = Constants.c;
    const double sigma_T = Constants.sigma_T;
    double H = cosmo->H_of_x(x);
    return -c * sigma_T * Saha_ne_of_x(x) / H;;
}


double RecombinationHistory::ddtauddx_of_x(double x) const {

    return dtaudx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const {
    return g_tilde_of_x_spline(x);
}
double RecombinationHistory::Saha_Xe_of_x(double x) const {
    return Xe_Saha_spline(x);
}
double RecombinationHistory::Saha_ne_of_x(double x) const {
    return ne_Saha_spline(x);
}



double RecombinationHistory::dgdx_tilde_of_x(double x) const {
    return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const {

    return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const {
    return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const {
    return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const {
    return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const {
    std::cout << "\n";
    std::cout << "Info about recombination/reionization history class:\n";
    std::cout << "Yp:          " << Yp << "\n";
    std::cout << std::endl;
}

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const {
    std::ofstream fp(filename.c_str());
    const int npts = 5000;
    const double x_min = x_start;
    const double x_max = x_end;

    Vector x_array = Utils::linspace(x_min, x_max, npts);
    auto print_data = [&](const double x) {
        fp << x << " ";
        fp << Xe_of_x(x) << " ";
        fp << ne_of_x(x) << " ";
        fp << tau_of_x(x) << " ";
        fp << dtaudx_of_x(x) << " ";
        fp << ddtauddx_of_x(x) << " ";
        fp << g_tilde_of_x(x) << " ";
        fp << dgdx_tilde_of_x(x) << " ";
        fp << ddgddx_tilde_of_x(x) << " ";
        fp << Saha_Xe_of_x(x) << " ";
        fp << "\n";
    };
    std::for_each(x_array.begin(), x_array.end(), print_data);
}
