#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================

BackgroundCosmology::BackgroundCosmology(
    double h,
    double OmegaB,
    double OmegaCDM,
    double OmegaK,
    double Neff,
    double TCMB) :
    h(h),
    OmegaB(OmegaB),
    OmegaCDM(OmegaCDM),
    OmegaK(OmegaK),
    Neff(Neff),
    TCMB(TCMB)
{
    
    //=============================================================================
    // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
    //=============================================================================
  
    H0 = 100*h*Constants.km/Constants.s/Constants.Mpc;                      // The Hubble parameter today H0 = 100h km/s/Mpc
    OmegaR = 2*(std::pow(M_PI,3))*Constants.G*(std::pow(Constants.k_b*TCMB,4))*8/(90.0*pow(Constants.c,5)*pow(Constants.hbar,3)*(std::pow(H0,2)));                  // Photon density today (follows from TCMB)
    OmegaNu = 0 ;                 // Neutrino density today (follows from TCMB and Neff)
    OmegaLambda = 1- OmegaB - OmegaCDM - OmegaK- OmegaR - OmegaNu;                  // // Dark energy density today
    
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve() {
    Utils::StartTiming("Eta");
    int npts = 10001;
    double x_start = -18.0; double x_end =  1.0;
    //=============================================================================
    // TODO: Set the range of x and the number of points for the splines
    // For this Utils::linspace(x_start, x_end, npts) is useful
    //=============================================================================
    Vector x_array = Utils::linspace(x_start, x_end, npts);

    // The ODE for deta/dx
    ODEFunction detadx = [&](double x, const double* eta, double* detadx) {
        detadx[0] = 1/Hp_of_x(x);
        return GSL_SUCCESS;
    };
    double etaini = 1/Hp_of_x(x_start);
    Vector eta{ etaini };
    ODESolver ode;
    ode.solve(detadx, x_array, eta);
    auto eta_array = ode.get_data_by_component(0);
    
    eta_of_x_spline.create(x_array, eta_array, "eta_of_x");
    std::cout << "eta at -7 is " << eta_of_x_spline(-7.0) << "\n";  
    std::cout << "Hp at -7 is " << Hp_of_x(-7.0) << "\n" ;

    //=============================================================================
    // TODO: Set the initial condition, set up the ODE system, solve and make
    // the spline eta_of_x_spline 
    //=============================================================================
    

    Utils::EndTiming("Eta");
    
    ODEFunction dtdx = [&](double x, const double* time, double* dtdx) {
        dtdx[0] = 1 / H_of_x(x);
        return GSL_SUCCESS;
    };
    Vector time{1/(2.*H_of_x(x_start))};
    ode.solve(dtdx, x_array, time);
    auto time_array = ode.get_data_by_component(0);
    t_of_x_spline.create(x_array, time_array, "t_of_x");
    std::cout << "The age of the universe is:" << t_of_x_spline(0)/GYear << " Gigayears" << "\n";
    std::cout << "The conformal time today is: " << eta_of_x_spline(0) / GYear << " Gigayears"<< "\n";

}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const {

    return H0*sqrt((OmegaB+OmegaCDM)*exp(-3*x)+(OmegaR)*exp(-4*x) + OmegaLambda   );
}

double BackgroundCosmology::Hp_of_x(double x) const {

    return H_of_x(x)*exp(x);
}

double BackgroundCosmology::dHpdx_of_x(double x) const {

    return H0*(-exp(-2.0*x)*(OmegaB+OmegaCDM) -2.0*exp(-3.0*x)*OmegaR+2.0*exp(x)*OmegaLambda)/(2.0*sqrt((OmegaB + OmegaCDM) * exp(-3.0*x) + (OmegaR) * exp(-4.0 * x) + OmegaLambda));
}

double BackgroundCosmology::ddHpddx_of_x(double x) const {
    
    return H0*(exp(-3 * x) * (exp(2 * x) * (pow(OmegaB + OmegaCDM, 2) + 14 * (OmegaB + OmegaCDM) * OmegaLambda * exp(3 * x) + 4 * pow(OmegaLambda, 2) * exp(6 * x)) + 6 * OmegaR * exp(x) * ((OmegaB + OmegaCDM) + 4 * OmegaLambda * exp(3 * x)) + 4 * OmegaR * OmegaR)) / (4 * ((OmegaR + (OmegaB + OmegaCDM) * exp(x)) + OmegaLambda * exp(4 * x)) * sqrt(OmegaLambda + exp(-4 * x) * (OmegaR + (OmegaB + OmegaCDM) * exp(x))));
}


double BackgroundCosmology::get_OmegaB(double x) const {
    if (x == 0.0) return OmegaB;
 
    return OmegaB*exp(-3*x)*std::pow(H0,2)/std::pow(H_of_x(x),2);
}


double BackgroundCosmology::get_OmegaR(double x) const {
    if (x == 0.0) return OmegaR;

    return OmegaR*exp(-4*x)*std::pow(H0,2)/std::pow(H_of_x(x),2);
}

double BackgroundCosmology::get_OmegaNu(double x) const {
    if (x == 0.0) return OmegaNu;
 
    return 0.0;
}

double BackgroundCosmology::get_OmegaCDM(double x) const {
    if (x == 0.0) return OmegaCDM;
 
    return OmegaCDM*exp(-3*x)*std::pow(H0,2)/std::pow(H_of_x(x),2);
}

double BackgroundCosmology::get_OmegaLambda(double x) const {
    if (x == 0.0) return OmegaLambda;
 
    return OmegaLambda*std::pow(H0,2)/std::pow(H_of_x(x),2);
}

double BackgroundCosmology::get_OmegaK(double x) const {
    if (x == 0.0) return OmegaK;

    return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const {
    return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const {
    return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const {
    return H0;
}

double BackgroundCosmology::get_h() const {
    return h;
}

double BackgroundCosmology::get_Neff() const {
    return Neff;
}

double BackgroundCosmology::get_TCMB(double x) const {
    if (x == 0.0) return TCMB;
    return TCMB * exp(-x);
}




double BackgroundCosmology::get_z(double x) const {
    return exp(-x) - 1.0;
}

double BackgroundCosmology::get_dL(double x) const {
    return (eta_of_x_spline(0) - eta_of_x_spline(x)) * exp(-x)*Constants.c/(1000*Constants.Mpc);
}
//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const {
    std::cout << "\n";
    std::cout << "Info about cosmology class:\n";
    std::cout << "OmegaB:      " << OmegaB << "\n";
    std::cout << "OmegaCDM:    " << OmegaCDM << "\n";
    std::cout << "OmegaLambda: " << OmegaLambda << "\n";
    std::cout << "OmegaK:      " << OmegaK << "\n";
    std::cout << "OmegaNu:     " << OmegaNu << "\n";
    std::cout << "OmegaR:      " << OmegaR << "\n";
    std::cout << "Neff:        " << Neff << "\n";
    std::cout << "h:           " << h << "\n";
    std::cout << "TCMB:        " << TCMB << "\n";
    std::cout << "H0:          " << H0 << "\n";
    std::cout << "z for RM-equality, acceleration, MDE equality:" << get_z(-8.65) <<"  " << get_z(-0.49) << "  " <<get_z(-0.26) << "\n";
    std::cout << "t for RM-equality, acceleration, MDE equality:" << t_of_x(-8.65)/GYear << "  " << t_of_x(-0.49)/GYear <<"  " << t_of_x(-0.26)/GYear << "\n";
    std::cout << std::endl;
}
//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const {
    const double x_min = -18.0;
    const double x_max = 2.0;
    const int    n_pts = 10001;

    Vector x_array = Utils::linspace(x_min, x_max, n_pts);
    std::cout << "K_equality is " << exp(-8.657) * H_of_x(-8.657) / Constants.c << "\n";

    std::ofstream fp(filename.c_str());
    auto print_data = [&](const double x) {
        test1(x);
        fp << x << " ";
        fp << eta_of_x(x)*Constants.c/Constants.Mpc << " ";
        fp << Hp_of_x(x)/H0 << " ";
        fp << dHpdx_of_x(x)/H0 << " ";
        fp << get_OmegaB(x) << " ";
        fp << get_OmegaCDM(x) << " ";
        fp << get_OmegaLambda(x) << " ";
        fp << get_OmegaR(x) << " ";
        fp << H_of_x(x) / H0<< " ";
        fp << ddHpddx_of_x(x)/H0 << " ";
        fp << t_of_x(x)/GYear << " ";
        fp << get_z(x) << " ";
        fp << get_dL(x) << " ";
        fp << "\n";
    };
    std::for_each(x_array.begin(), x_array.end(), print_data);
}



void BackgroundCosmology::test1(double x) const{
    double err = 1E-3;
    if (abs(1 - (Hp_of_x(x) / (exp(x) * H_of_x(x))) > err)) {
        std::cout << "Hp / (exp(x)*H != 1 for x = " << x << "\n";
        exit(1);
    }
}
