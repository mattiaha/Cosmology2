#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology* cosmo,
    RecombinationHistory* rec) :
    cosmo(cosmo),
    rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve() {

    // Integrate all the perturbation equation and spline the result
    integrate_perturbations();

    // Compute source functions and spline the result
    //compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations() {
    Utils::StartTiming("integrateperturbation");

    //===================================================================
    // TODO: Set up the k-array for the k's we are going to integrate over
    // Start at k_min end at k_max with n_k points with either a
    // quadratic or a logarithmic spacing
    //===================================================================
    const int n_points = n_x;
    const int n_2 = 2 * n_points;
    Vector2D vb_array(n_2, Vector(n_k, 0.0));         Vector2D vcdm_array(n_2, Vector(n_k, 0.0));      Vector2D Phi_array(n_2, Vector(n_k, 0.0));
    Vector2D Psi_array(n_2, Vector(n_k, 0.0));        Vector2D delta_b_array(n_2, Vector(n_k, 0.0));   Vector2D delta_cdm_array(n_2, Vector(n_k, 0.0));
    Vector2D Theta0_array(n_2, Vector(n_k, 0.0));   Vector2D Theta1_array(n_2, Vector(n_k, 0.0));    Vector2D Theta2_array(n_2, Vector(n_k, 0.0));
    Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));


    // Loop over all wavenumbers
    for (int ik = 0; ik < n_k; ik++) {

        // Progress bar...
        if ((10 * ik) / n_k != (10 * ik + 10) / n_k) {
            std::cout << (100 * ik + 100) / n_k << "% " << std::flush;
            if (ik == n_k - 1) std::cout << std::endl;
        }
        const double c = Constants.c;

        // Current value of k
        double k = k_array[ik];
        // Find value to integrate to
        double x_end_tight = get_tight_coupling_time(k);
        Vector x_array_tight = Utils::linspace(x_start, x_end_tight, n_points);
        //===================================================================
        // TODO: Tight coupling integration
        // Remember to implement the routines:
        // set_ic : The IC at the start
        // rhs_tight_coupling_ode : The dydx for our coupled ODE system
        //===================================================================

        // Set up initial conditions in the tight coupling regime
        auto y_tight_coupling_ini = set_ic(x_start, k);
        ODESolver ode_tight;
        // The tight coupling ODE system
        ODEFunction dydx_tight_coupling = [&](double x, const double* y, double* dydx) {
            return rhs_tight_coupling_ode(x, k, y, dydx);
        };
        // Integrate from x_start -> x_end_tight
        ode_tight.solve(dydx_tight_coupling, x_array_tight, y_tight_coupling_ini);
        auto y_tight_coupling = ode_tight.get_data();

        Vector y_tc = y_tight_coupling[n_points - 1];


        //====i===============================================================
        // TODO: Full equation integration
        // Remember to implement the routines:
        // set_ic_after_tight_coupling : The IC after tight coupling ends
        // rhs_full_ode : The dydx for our coupled ODE system
        //===================================================================

        // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
        auto y_full_ini = set_ic_after_tight_coupling(y_tc, x_end_tight, k);

        // The full ODE system
        ODEFunction dydx_full = [&](double x, const double* y, double* dydx) {
            return rhs_full_ode(x, k, y, dydx);
        };
        ODESolver ode_full;
        Vector x_array_full = Utils::linspace(x_end_tight, x_end, n_points);
        ode_full.solve(dydx_full, x_array_full, y_full_ini);
        auto y_full = ode_full.get_data();

        for (int ix = 0; ix < n_points; ix++) {
            const double x = x_array_tight[ix];
            const double Hp = cosmo->Hp_of_x(x);
            const double dtau = rec->dtaudx_of_x(x);
            delta_cdm_array[ix][ik] = y_tight_coupling[ix][Constants.ind_deltacdm];
            delta_b_array[ix][ik] = y_tight_coupling[ix][Constants.ind_deltab];
            vcdm_array[ix][ik] = y_tight_coupling[ix][Constants.ind_vcdm];
            vb_array[ix][ik] = y_tight_coupling[ix][Constants.ind_vb];
            Phi_array[ix][ik] = y_tight_coupling[ix][Constants.ind_Phi];
            Theta0_array[ix][ik] = y_tight_coupling[ix][Constants.ind_start_theta];
            Theta1_array[ix][ik] = y_tight_coupling[ix][Constants.ind_start_theta + 1];
            Theta2_array[ix][ik] = -20 * c * k / (45.0 * Hp * dtau) * Theta1_array[ix][ik];
            Psi_array[ix][ik] = -Phi_array[ix][ik];


        }
        const double H0 = cosmo->get_H0();
        const double Omega_R0 = cosmo->get_OmegaR(0);
        for (int ix = n_points; ix < n_2; ix++) {
            const double x = x_array_full[ix - n_points];
            const double dtau = rec->dtaudx_of_x(x);

            delta_cdm_array[ix][ik] = y_full[ix - n_points][Constants.ind_deltacdm];
            delta_b_array[ix][ik] = y_full[ix - n_points][Constants.ind_deltab];
            vcdm_array[ix][ik] = y_full[ix - n_points][Constants.ind_vcdm];
            vb_array[ix][ik] = y_full[ix - n_points][Constants.ind_vb];
            Phi_array[ix][ik] = y_full[ix - n_points][Constants.ind_Phi];
            Theta0_array[ix][ik] = y_full[ix - n_points][Constants.ind_start_theta];
            Theta1_array[ix][ik] = y_full[ix - n_points][Constants.ind_start_theta + 1];
            Theta2_array[ix][ik] = y_full[ix - n_points][Constants.ind_start_theta + 2];
            Psi_array[ix][ik] = -Phi_array[ix][ik] - Omega_R0 * Theta2_array[ix][ik] * 12 * pow(H0/(c*k), 2)* exp(-2 * x);
        }


    }
    Vector x_full = Utils::linspace(x_start, x_end, n_2);
    delta_cdm_spline.create(x_full, k_array, delta_cdm_array, "delta_cdm");
    delta_b_spline.create(x_full, k_array, delta_b_array, "delta_b");
    v_cdm_spline.create(x_full, k_array, vcdm_array, "v_cdm");
    v_b_spline.create(x_full, k_array, vb_array, "v_b");
    Phi_spline.create(x_full, k_array, Phi_array, "Phi");
    Psi_spline.create(x_full, k_array, Psi_array, "Psi");
    Theta0_spline.create(x_full, k_array, Theta0_array, "Theta0");
    Theta1_spline.create(x_full, k_array, Theta1_array, "Theta1");
    Theta2_spline.create(x_full, k_array, Theta2_array, "Theta2");

    Utils::EndTiming("integrateperturbation");

    //=============================================================================
    // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
    //=============================================================================
    // ...
    // ...
    // ...
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const {

    // The vector we are going to fill
    Vector y_tc(Constants.n_ell_tot_tc);

    //=============================================================================
    // Compute where in the y_tc array each component belongs
    // This is just an example of how to do it to make it easier
    // Feel free to organize the component any way you like
    //=============================================================================

    // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
    const int n_ell_theta_tc = Constants.n_ell_theta_tc;
    //const int n_ell_neutrinos_tc = Constants.n_ell_neutrinos_tc;
    const int n_ell_tot_tc = Constants.n_ell_tot_tc;
    const bool polarization = Constants.polarization;
    const bool neutrinos = Constants.neutrinos;

    // References to the tight coupling quantities
    double& delta_cdm = y_tc[Constants.ind_deltacdm_tc];
    double& delta_b = y_tc[Constants.ind_deltab_tc];
    double& v_cdm = y_tc[Constants.ind_vcdm_tc];
    double& v_b = y_tc[Constants.ind_vb_tc];
    double& Phi = y_tc[Constants.ind_Phi_tc];
    double* Theta = &y_tc[Constants.ind_start_theta_tc];
    //double* Nu = &y_tc[Constants.ind_start_nu_tc];

    const double c = Constants.c;
    const double Hp = cosmo->Hp_of_x(x);
    const double dtau = rec->dtaudx_of_x(x);

    //=============================================================================
    // TODO: Set the initial conditions in the tight coupling regime
    //=============================================================================
    Phi = 2.0 / 3.0;
    delta_cdm = 1.0;
    delta_b = 1.0;
    v_cdm = c * k / (3.0 * Hp);
    v_b = c * k / (3.0 * Hp);
    Theta[0] = 1.0 / 3.0;
    Theta[1] = -c * k / (9.0 * Hp);
    return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector& y_tc,
    const double x,
    const double k) const {

    // Make the vector we are going to fill
    Vector y(Constants.n_ell_tot_full);

    //=============================================================================
    // Compute where in the y array each component belongs and where corresponding
    // components are located in the y_tc array
    // This is just an example of how to do it to make it easier
    // Feel free to organize the component any way you like
    //=============================================================================

    // Number of multipoles we have in the full regime
    const int n_ell_theta = Constants.n_ell_theta;
    //const int n_ell_thetap = Constants.n_ell_thetap;
    //const int n_ell_neutrinos = Constants.n_ell_neutrinos;
    //const bool polarization = Constants.polarization;
    //const bool neutrinos = Constants.neutrinos;

    // Number of multipoles we have in the tight coupling regime
    const int n_ell_theta_tc = Constants.n_ell_theta_tc;
    //const int n_ell_neutrinos_tc = Constants.n_ell_neutrinos_tc;

    // References to the tight coupling quantities
    const double& delta_cdm_tc = y_tc[Constants.ind_deltacdm_tc];
    const double& delta_b_tc = y_tc[Constants.ind_deltab_tc];
    const double& v_cdm_tc = y_tc[Constants.ind_vcdm_tc];
    const double& v_b_tc = y_tc[Constants.ind_vb_tc];
    const double& Phi_tc = y_tc[Constants.ind_Phi_tc];
    const double* Theta_tc = &y_tc[Constants.ind_start_theta_tc];
    //const double* Nu_tc = &y_tc[Constants.ind_start_nu_tc];

    // References to the quantities we are going to set
    double& delta_cdm = y[Constants.ind_deltacdm_tc];
    double& delta_b = y[Constants.ind_deltab_tc];
    double& v_cdm = y[Constants.ind_vcdm_tc];
    double& v_b = y[Constants.ind_vb_tc];
    double& Phi = y[Constants.ind_Phi_tc];
    double* Theta = &y[Constants.ind_start_theta_tc];
    //double* Theta_p = &y[Constants.ind_start_thetap_tc];
    //double* Nu = &y[Constants.ind_start_nu_tc];
    const double c = Constants.c;
    const double Hp = cosmo->Hp_of_x(x);
    const double dtau = rec->dtaudx_of_x(x);
    //=============================================================================
    // TODO: fill in the initial conditions for the full equation system below
    // NB: remember that we have different number of multipoles in the two
    // regimes so be careful when assigning from the tc array
    //=============================================================================
    delta_cdm = delta_cdm_tc;
    delta_b = delta_b_tc;
    v_cdm = v_cdm_tc;
    v_b = v_b_tc;
    Phi = Phi_tc;
    Theta[0] = Theta_tc[0];
    Theta[1] = Theta_tc[1];
    Theta[2] = -20 * c * k * Theta[1] / (45 * Hp * dtau);

    for (int i = 3; i < n_ell_theta; i++) {
        Theta[i] = -i * c * k * Theta[i - 1] / ((2.0 * i + 1.0) * Hp * dtau);
    }


    return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const {
    double x_tight_coupling_end = 0.0;

    //=============================================================================
    // TODO: compute and return x for when tight coupling ends
    // Remember all the three conditions in Callin
    //=============================================================================
    // ...
    int n_time = 1000;
    double const rec_start = -8.0;
    double const c = Constants.c;
    Vector x_array = Utils::linspace(x_start, rec_start, n_time);

    for (int i = 0; i < n_time - 1; i++) {
        double x = x_array[i];
        double dtau = rec->dtaudx_of_x(x);
        double Hp = cosmo->Hp_of_x(x);
        if (abs(dtau) < 10.0 && abs(dtau) < abs(c * k / Hp)) {
            x_tight_coupling_end = x;
            break;

        }

        x_tight_coupling_end = rec_start;
    }
    //std::cout << x_tight_coupling_end << "\n";
    //std::cout << abs(rec->dtaudx_of_x(rec_start)) << "\n";
    return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions() {
    Utils::StartTiming("source");

    //=============================================================================
    // TODO: Make the x and k arrays to evaluate over and use to make the splines
    //=============================================================================
    // ...
    // ...
    Vector k_array;
    Vector x_array;

    // Make storage for the source functions (in 1D array to be able to pass it to the spline)
    Vector ST_array(k_array.size() * x_array.size());
    Vector SE_array(k_array.size() * x_array.size());

    // Compute source functions
    for (auto ix = 0; ix < x_array.size(); ix++) {
        const double x = x_array[ix];
        for (auto ik = 0; ik < k_array.size(); ik++) {
            const double k = k_array[ik];

            // NB: This is the format the data needs to be stored 
            // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
            const int index = ix + n_x * ik;

            //=============================================================================
            // TODO: Compute the source functions
            //=============================================================================
            // Fetch all the things we need...
            // const double Hp       = cosmo->Hp_of_x(x);
            // const double tau      = rec->tau_of_x(x);
            // ...
            // ...

            // Temperatur source



            ST_array[index] = 0.0;

            // Polarization source
            if (Constants.polarization) {
                SE_array[index] = 0.0;
            }
        }
    }

    // Spline the source functions
    //ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");
    if (Constants.polarization) {
        SE_spline.create(x_array, k_array, SE_array, "Source_Pol_x_k");
    }

    Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double* y, double* dydx) {

    //=============================================================================
    // Compute where in the y / dydx array each component belongs
    // This is just an example of how to do it to make it easier
    // Feel free to organize the component any way you like
    //=============================================================================
    // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
    const int n_ell_theta_tc = Constants.n_ell_theta_tc;
    // const int n_ell_neutrinos_tc = Constants.n_ell_neutrinos_tc;
    //const bool neutrinos = Constants.neutrinos;

    // The different quantities in the y array
    const double& delta_cdm = y[Constants.ind_deltacdm_tc];
    const double& delta_b = y[Constants.ind_deltab_tc];
    const double& v_cdm = y[Constants.ind_vcdm_tc];
    const double& v_b = y[Constants.ind_vb_tc];
    const double& Phi = y[Constants.ind_Phi_tc];
    const double* Theta = &y[Constants.ind_start_theta_tc];
    //const double* Nu = &y[Constants.ind_start_nu_tc];

    // References to the quantities we are going to set in the dydx array
    double& ddelta_cdmdx = dydx[Constants.ind_deltacdm_tc];
    double& ddelta_bdx = dydx[Constants.ind_deltab_tc];
    double& dv_cdmdx = dydx[Constants.ind_vcdm_tc];
    double& dv_bdx = dydx[Constants.ind_vb_tc];
    double& dPhidx = dydx[Constants.ind_Phi_tc];
    double* dThetadx = &dydx[Constants.ind_start_theta_tc];
    //double* dNudx = &dydx[Constants.ind_start_nu_tc];

    const double H0 = cosmo->get_H0();
    const double Hp = cosmo->Hp_of_x(x);
    const double dHp = cosmo->dHpdx_of_x(x);
    const double c = Constants.c;
    const double dtaudx = rec->dtaudx_of_x(x);
    const double ddtauddx = rec->ddtauddx_of_x(x);
    const double OmegaCDM0 = cosmo->get_OmegaCDM(0);
    const double OmegaB0 = cosmo->get_OmegaB(0);
    const double OmegaR0 = cosmo->get_OmegaR(0);
    const double Theta2 = 0;

    const double Psi = -Phi - 12 * pow(H0, 2) * exp(-2 * x) / pow(c * k, 2) * OmegaR0 * Theta2;
    const double R = 4.0 * OmegaR0 / (3.0 * OmegaB0 * exp(x));

    dPhidx = Psi - (pow(c*k/Hp, 2) * Phi)/3.0 + (pow(H0/ Hp, 2) / (2.0)) * (exp(-x) * (OmegaCDM0 * delta_cdm + OmegaB0 * delta_b) + 4 * OmegaR0 * exp(-2 * x) * Theta[0]);
    dThetadx[0] = -c * k * Theta[1] / Hp - dPhidx;
    const double q = (-(3.0 * Theta[1] + v_b) * ((1 - R) * dtaudx + (1 + R) * ddtauddx) + (c*k/Hp)*(-Psi+(1-dHp/Hp)*(-Theta[0]) -dThetadx[0]))/((1+R)*dtaudx + dHp/Hp -1.0);
    dv_bdx = (-v_b-(c*k/Hp)*Psi +R*(q -(c*k/Hp)*( Theta[0]+Psi)))/ (1.0 + R);

    dThetadx[1] = (1.0 / 3.0) * (q - dv_bdx);
    ddelta_cdmdx = c * k * v_cdm - 3 * dPhidx;
    dv_cdmdx = -v_cdm - c * k * Psi / Hp;
    ddelta_bdx = c * k * v_b / Hp - 3 * dPhidx;

    return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double* y, double* dydx) {

    //=============================================================================
    // Compute where in the y / dydx array each component belongs
    // This is just an example of how to do it to make it easier
    // Feel free to organize the component any way you like
    //=============================================================================

    // Index and number of the different quantities
    const int n_ell_theta = Constants.n_ell_theta;


    // The different quantities in the y array
    const double& delta_cdm = y[Constants.ind_deltacdm];
    const double& delta_b = y[Constants.ind_deltab];
    const double& v_cdm = y[Constants.ind_vcdm];
    const double& v_b = y[Constants.ind_vb];
    const double& Phi = y[Constants.ind_Phi];
    const double* Theta = &y[Constants.ind_start_theta];


    // References to the quantities we are going to set in the dydx array
    double& ddelta_cdmdx = dydx[Constants.ind_deltacdm];
    double& ddelta_bdx = dydx[Constants.ind_deltab];
    double& dv_cdmdx = dydx[Constants.ind_vcdm];
    double& dv_bdx = dydx[Constants.ind_vb];
    double& dPhidx = dydx[Constants.ind_Phi];
    double* dThetadx = &dydx[Constants.ind_start_theta];


    const double H0 = cosmo->get_H0();
    const double Hp = cosmo->Hp_of_x(x);
    const double dHp = cosmo->dHpdx_of_x(x);
    const double c = Constants.c;
    const double dtaudx = rec->dtaudx_of_x(x);
    const double ddtauddx = rec->ddtauddx_of_x(x);
    const double OmegaCDM0 = cosmo->get_OmegaCDM(0);
    const double OmegaB0 = cosmo->get_OmegaB(0);
    const double OmegaR0 = cosmo->get_OmegaR(0);
    const double eta = cosmo->eta_of_x(x);
    // eta_of_x is defined without constant c
    const double cHe = 8.0 / (Hp * eta);
    //std::cout << cHe << "\n";
    const double Psi = -Phi - 12 * pow(H0, 2) * exp(-x) * OmegaR0 * Theta[2] / pow(c * k, 2);
    const double R = 4.0 * OmegaR0 / (3.0 * OmegaB0 * exp(x));


    dPhidx = Psi - (pow(c, 2) * pow(k, 2) * Phi) / (3 * pow(Hp, 2)) + pow(H0, 2) / (2 * pow(Hp, 2)) * (exp(-x) * (OmegaCDM0 * delta_cdm + OmegaB0 * delta_b) + 4 * OmegaR0 * exp(-2 * x) * Theta[0]);
    dThetadx[0] = -c * k * Theta[1] / Hp - dPhidx;
    dThetadx[1] = (c * k / (3 * Hp)) * (Theta[0] - 2 * Theta[2] + Psi) + dtaudx * (Theta[1] + 1 / 3.0 * v_b);
    dThetadx[2] = 2 * c * k * Theta[1] / (5.0 * Hp) - 3 * c * k * Theta[3] / (5.0 * Hp) + dtaudx * 9 / 10.0 * Theta[2];
    ddelta_cdmdx = c * k * v_cdm / Hp - 3 * dPhidx;
    dv_cdmdx = -v_cdm - c * k * Psi / Hp;
    dv_bdx = -v_b - c * k * Psi / Hp + dtaudx * R * (3 * Theta[1] + v_b);
    ddelta_bdx = c * k * v_b / Hp - 3 * dPhidx;


    for (int l = 3; l < n_ell_theta - 1; l++)
    {
        dThetadx[l] = l * c * k * Theta[l - 1] / ((2.0 * l + 1) * Hp) - (l + 1) * c * k * Theta[l + 1] / ((2.0 * l + 1) * Hp) + dtaudx * Theta[l];
    }

    dThetadx[n_ell_theta - 1] = (c * k / Hp) * Theta[n_ell_theta - 2] - cHe * Theta[n_ell_theta - 1] + dtaudx * Theta[n_ell_theta - 1];

    return GSL_SUCCESS;
}


//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const {
    return delta_cdm_spline(x, k);
}
double Perturbations::get_delta_b(const double x, const double k) const {
    return delta_b_spline(x, k);
}
double Perturbations::get_v_cdm(const double x, const double k) const {
    return v_cdm_spline(x, k);
}
double Perturbations::get_v_b(const double x, const double k) const {
    return v_b_spline(x, k);
}
double Perturbations::get_Phi(const double x, const double k) const {
    return Phi_spline(x, k);
}
double Perturbations::get_Psi(const double x, const double k) const {
    return Psi_spline(x, k);
}
double Perturbations::get_Source_T(const double x, const double k) const {
    return ST_spline(x, k);
}
double Perturbations::get_Source_E(const double x, const double k) const {
    return SE_spline(x, k);
}
double Perturbations::get_Theta0(const double x, const double k) const {
    return Theta0_spline(x, k);
}
double Perturbations::get_Theta1(const double x, const double k) const {
    return Theta1_spline(x, k);
}
double Perturbations::get_Theta2(const double x, const double k) const {
    return Theta2_spline(x, k);
}
//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const {
    std::cout << "\n";
    std::cout << "Info about perturbations class:\n";
    std::cout << "x_start:       " << x_start << "\n";
    std::cout << "x_end:         " << x_end << "\n";
    std::cout << "n_x:     " << 2 * n_x << "\n";
    std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc << "\n";
    std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc << "\n";
    std::cout << "n_k:     " << n_k << "\n";
    if (Constants.polarization)
        std::cout << "We include polarization\n";
    else
        std::cout << "We do not include polarization\n";
    if (Constants.neutrinos)
        std::cout << "We include neutrinos\n";
    else
        std::cout << "We do not include neutrinos\n";

    std::cout << "Information about the perturbation system:\n";
    std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm << "\n";
    std::cout << "ind_deltab:         " << Constants.ind_deltab << "\n";
    std::cout << "ind_v_cdm:          " << Constants.ind_vcdm << "\n";
    std::cout << "ind_v_b:            " << Constants.ind_vb << "\n";
    std::cout << "ind_Phi:            " << Constants.ind_Phi << "\n";
    std::cout << "ind_start_theta:    " << Constants.ind_start_theta << "\n";
    std::cout << "n_ell_theta:        " << Constants.n_ell_theta << "\n";
    std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full << "\n";

    std::cout << "Information about the perturbation system in tight coupling:\n";
    std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc << "\n";
    std::cout << "ind_deltab:         " << Constants.ind_deltab_tc << "\n";
    std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc << "\n";
    std::cout << "ind_v_b:            " << Constants.ind_vb_tc << "\n";
    std::cout << "ind_Phi:            " << Constants.ind_Phi_tc << "\n";
    std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc << "\n";
    std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc << "\n";
    std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc << "\n";
    std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const {
    std::ofstream fp(filename.c_str());
    const int npts = 2*n_x;
    auto x_array = Utils::linspace(x_start, x_end, npts);
    auto print_data = [&](const double x) {
        double arg = k * Constants.c * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
        fp << x << " ";
        fp << get_delta_cdm(x, k) << " ";
        fp << get_delta_b(x, k) << " ";
        fp << get_v_cdm(x, k) << " ";
        fp << get_v_b(x, k) << " ";
        fp << get_Theta0(x, k) << " ";
        fp << get_Theta1(x, k) << " ";
        fp << get_Phi(x, k) << " ";
        fp << get_Psi(x, k) << " ";

        //fp << get_Source_T(x, k) << " ";
        //fp << get_Source_T(x, k) * Utils::j_ell(5, arg) << " ";
        //fp << get_Source_T(x, k) * Utils::j_ell(50, arg) << " ";
        //fp << get_Source_T(x, k) * Utils::j_ell(500, arg) << " ";
        fp << "\n";
    };
    std::for_each(x_array.begin(), x_array.end(), print_data);
}
