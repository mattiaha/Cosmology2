#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology* cosmo,
    RecombinationHistory* rec,
    Perturbations* pert,
    double A_s,
    double n_s,
    double kpivot_mpc) :
    cosmo(cosmo),
    rec(rec),
    pert(pert),
    A_s(A_s),
    n_s(n_s),
    kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve() {
    /*const double eta_0 = cosmo->eta_of_x(0);
    const double c = Constants.c;
    const double k_step = 2 * M_PI / (32.0 * eta_0 * c);
    const int n_k = k_max / k_step;*/
    //=========================================================================
    // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
    //=========================================================================
    const double eta_0 = cosmo->eta_of_x(0);
    const double k_step = 2 * M_PI / (32.0 * eta_0*Constants.c);
    const double nk = 1000;//(k_max - k_min) / k_step;

    std::cout << eta_0 * Constants.c*k_max << "\n";
    std::cout << nk << "\n";
    Vector k_array = Utils::linspace(k_min,k_max,nk);
    Vector log_k_array = Utils::linspace(log(k_min), log(k_max), nk);
    

    //=========================================================================
    // TODO: Make splines for j_ell. 
    // Implement generate_bessel_function_splines
    //=========================================================================
    generate_bessel_function_splines();

    //=========================================================================
    // TODO: Line of sight integration to get Theta_ell(k)
    // Implement line_of_sight_integration
    //=========================================================================
    line_of_sight_integration(k_array);

    //=========================================================================
    // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
    // Implement solve_for_cell
    //=========================================================================
    auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
    cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");

    
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines() {
    Utils::StartTiming("besselspline");
    const double c = Constants.c;
    const double eta_0 = cosmo->eta_of_x(0)*c;
    // Make storage for the splines
    j_ell_splines = std::vector<Spline>(ells.size());
    const double z_min = 0.0;
    const double z_max = k_max * eta_0;
    const double z_step = 2.0 * M_PI / 16.0;
    const int n_points = z_max / z_step;

    Vector z_array = Utils::linspace(z_min, z_max, n_points);
    Vector j_ell_array(n_points);
    //=============================================================================
    // TODO: Compute splines for bessel functions j_ell(z)
    // Choose a suitable range for each ell
    // NB: you don't want to go larger than z ~ 40000, then the bessel routines
    // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
    //=============================================================================
    for (size_t i = 0; i < ells.size(); i++) {
        const int ell = ells[i];
        for (int j = 0; j < n_points ; j++) {
            //const double z = z_array[j];
            j_ell_array[j] = Utils::j_ell(ell, z_array[j]);
            //std::cout << Utils::j_ell(ell, z) << "\n";
        }
        
        j_ell_splines[i].create(z_array, j_ell_array);

    }
    /*for (size_t i = 0; i < ells.size(); i++) {
        const int ell = ells[i];
        std::cout << j_ell_splines[5](400) << "\n";
        std::cout << j_ell_splines[5](800) << "\n";

    }*/

    Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector& k_array,
    std::function<double(double, double)>& source_function) {
    Utils::StartTiming("lineofsight");
    // Make storage for the results
    
    Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

    
    const int x_points = 2200;
    const double x_step = abs(x_start) / 2200.0;
    Vector x_array = Utils::linspace(x_start, x_end, x_points);

    for (size_t ik = 0; ik < k_array.size(); ik++) {
        
        const double k = k_array[ik];

        for (int i_ell = 0; i_ell < ells.size(); i_ell++) {
            double trapez = 0.0;
            for (int ix = 1; ix < x_array.size()-1; ix++)
            {
                const double x = x_array[ix];

                double arg = k * Constants.c * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
                //std::cout << "arg  " << arg << "\n";
                //std::cout << "source  " << source_function(x, k);
                //std::cout << j_ell_splines[i_ell] << "\n";
                trapez += 2*source_function(x, k)*j_ell_splines[i_ell](arg)*x_step;

            }
            //std::cout << trapez << "\n";
            const double arg0 = k * Constants.c * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x_array[0]));
            const double arge = k * Constants.c * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x_array[x_array.size()-1]));
            trapez += source_function(x_array[0], k) * j_ell_splines[i_ell](arg0);
            trapez += source_function(x_array[x_array.size()-1], k) * j_ell_splines[i_ell](arge);

            result[i_ell][ik] = trapez*x_step/2.0;
        }
        
    }
    
    Utils::EndTiming("lineofsight");
    return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector& k_array) {
    const int n_k = k_array.size();
    //const int n = 100;
    const int nells = ells.size();

    // Make storage for the splines we are to create
    //thetaT_ell_of_k_spline = std::vector<Spline>(nells);

    //============================================================================
    // TODO: Solve for Theta_ell(k) and spline the result
    //============================================================================

    // Make a function returning the source function
    std::function<double(double, double)> source_function_T = [&](double x, double k) {
        return pert->get_Source_T(x, k);
    };

    // Do the line of sight integration
    Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);
    thetaT_ell_of_k_spline.create(ells, k_array, thetaT_ell_of_k);
    
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector& log_k_array,
    Spline2D& f_ell_spline,
    Spline2D& g_ell_spline) {
    const int nells = ells.size();
    const double k_step = (k_max-k_min)/log_k_array.size();
    Vector result = Vector(nells);
    for (size_t i = 0; i < ells.size(); i++) {
        const int ell = ells[i];
        double trapez_cell = 0;
        for (size_t j = 0; j < log_k_array.size(); j++)
        {
            const double k = exp(log_k_array[j]); 
            const double f_ell = f_ell_spline(ell, k);
            const double g_ell = g_ell_spline(ell, k);
            trapez_cell += 4 * M_PI * primordial_power_spectrum(k)*f_ell*g_ell*k_step/k;
        }
        result[i] = trapez_cell;
        //std::cout << trapez_cell << "\n";
    }
    //============================================================================
    // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
    // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
    //============================================================================

    // ...
    // ...
    // ...
    // ...

    

    return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const {
    return A_s * pow(Constants.Mpc * k / kpivot_mpc, n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const {
    const double Phi = pert->get_Phi(x, k_mpc);
    const double H0 = cosmo->get_H0();
    const double OmegaM = cosmo->get_OmegaB(0) + cosmo->get_OmegaCDM(0);
    const double ckH0 = c * k_mpc / H0;
    const double d_m = pow(ckH0, 2) * Phi * 2.0 * exp(-x) / (3.0 * OmegaM);

    double pofk = pow(d_m, 2) * primordial_power_spectrum(k_mpc);
    

    return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const {
    return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const {
    return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const {
    return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const {
    // Output in standard units of muK^2
    std::ofstream fp(filename.c_str());
    const int ellmax = int(ells[ells.size() - 1]);
    auto ellvalues = Utils::linspace(2, ellmax, ellmax - 1);
    auto print_data = [&](const double ell) {
        double normfactor = (ell * (ell + 1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
        double normfactorN = (ell * (ell + 1)) / (2.0 * M_PI)
            * pow(1e6 * cosmo->get_TCMB() * pow(4.0 / 11.0, 1.0 / 3.0), 2);
        double normfactorL = (ell * (ell + 1)) * (ell * (ell + 1)) / (2.0 * M_PI);
        fp << ell << " ";
        fp << cell_TT_spline(ell) * normfactor << " ";
        if (Constants.polarization) {
            fp << cell_EE_spline(ell) * normfactor << " ";
            fp << cell_TE_spline(ell) * normfactor << " ";
        }
        fp << "\n";
    };
    std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}


void PowerSpectrum::MPSoutput(std::string filename) const {
    std::ofstream fp(filename.c_str());
    const double eta_0 = cosmo->eta_of_x(0);
    const double k_step = 2 * M_PI / (32.0 * eta_0 * Constants.c);
    const int nk = (k_max - k_min) / k_step;
    Vector k_array = Utils::linspace(k_min, k_max, nk);
    auto print_data = [&](const double k) {
        fp << k << " ";
        fp << get_matter_power_spectrum(0, k) << " ";
        fp << thetaT_ell_of_k_spline(6, k) << " ";
        fp << thetaT_ell_of_k_spline(100, k) << " ";
        fp << thetaT_ell_of_k_spline(200, k) << " ";
        fp << thetaT_ell_of_k_spline(500, k) << " ";
        fp << thetaT_ell_of_k_spline(1000, k) << " ";


        fp << "\n";
    };
    std::for_each(k_array.begin(), k_array.end(), print_data);
}

