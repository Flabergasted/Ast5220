#include <math.h>
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
  double k_b = Constants.k_b;
  double hbar = Constants.hbar;
  double c = Constants.c;
  double G = Constants.G;
  double H0_over_h = Constants.H0_over_h;
  
  H0 = H0_over_h * h;

  OmegaR = 2. * M_PI*M_PI/30. 
          * std::pow((k_b*TCMB), 4)/(std::pow(hbar,3)*std::pow(c,5))
          * (8.*M_PI*G)/(3.*H0*H0);

  OmegaNu = Neff * 7./8. * std::pow(4./11.,4./3.) * OmegaR;

  OmegaLambda = 1. - (OmegaK + OmegaB + OmegaCDM + OmegaR + OmegaNu);
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");

  double c = Constants.c;
  const double x_start = -20.0;
  const double x_end = 5.0;
  const int    npts = 1000;
  
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = c/Hp_of_x(x);
    return GSL_SUCCESS;
  };

  double etaini = c/Hp_of_x(x_start);
  Vector eta_ic {etaini};

  ODESolver ode;
  ode.solve(detadx, x_array, eta_ic);

  auto eta_array = ode.get_data_by_component(0);

  // create spline of eta
  eta_of_x_spline.create(x_array,eta_array,"eta_of_x");

  Utils::EndTiming("Eta");

  Utils::StartTiming("t");

  // The ODE for dt/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
    dtdx[0] = 1./H_of_x(x);
    return GSL_SUCCESS;
  };

  double tini = 1./(2.*H_of_x(x_start));
  Vector t_ic {tini};

  //ODESolver ode_t;
  ode.solve(dtdx, x_array, t_ic);

  auto t_array = ode.get_data_by_component(0);

  // create spline of t
  t_of_x_spline.create(x_array, t_array,"t_of_x");

  Utils::EndTiming("t");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  double density_sum_root = std::sqrt((OmegaB + OmegaCDM) * std::exp(-3.*x)
                + (OmegaR + OmegaNu) * std::exp(-4.*x)
                + OmegaK * std::exp(-2.*x) + OmegaLambda);
  
  return H0 * density_sum_root;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  double exp_x = std::exp(x);
  return H_of_x(x) * exp_x;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  double exp_x = std::exp(x);
  double density_sum_root = std::sqrt((OmegaB + OmegaCDM) * std::exp(-3.*x)
                + (OmegaR + OmegaNu) * std::exp(-4.*x)
                + OmegaK * std::exp(-2.*x) + OmegaLambda);
  double density_sum_1st_derivative = - 3. *(OmegaB + OmegaCDM) * std::exp(-3.*x)
                                      - 4. * (OmegaR + OmegaNu) * std::exp(-4.*x)
                                      - 2. * OmegaK * std::exp(-2.*x);
  return Hp_of_x(x) + (H0*exp_x * (density_sum_1st_derivative)) / (
                2. * density_sum_root);
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  double exp_x = std::exp(x);
  double density_sum_root = std::sqrt((OmegaB + OmegaCDM) * std::exp(-3.*x)
                + (OmegaR + OmegaNu) * std::exp(-4.*x)
                + OmegaK * std::exp(-2.*x) + OmegaLambda);
  double density_sum_1st_derivative = - 3. *(OmegaB + OmegaCDM) * std::exp(-3.*x)
                                      - 4. * (OmegaR + OmegaNu) * std::exp(-4.*x)
                                      - 2. * OmegaK * std::exp(-2.*x);
  double density_sum_2nd_derivative = 9. *(OmegaB + OmegaCDM) * std::exp(-3.*x)
                                      + 16. * (OmegaR + OmegaNu) * std::exp(-4.*x)
                                      + 4. * OmegaK * std::exp(-2.*x);
  return Hp_of_x(x) + H0 * exp_x * (density_sum_1st_derivative)/(density_sum_root)
                + H0 * exp_x * ((density_sum_2nd_derivative)/(2. * density_sum_root)
                - std::pow(density_sum_1st_derivative,2) / 
                (4. * std::pow(density_sum_root,3)));
}

double BackgroundCosmology::get_OmegaB(double x) const{
  double a = std::exp(x); 
  if(x == 0.0) return OmegaB;
  else;  
    return OmegaB/((std::pow(a, 3)*H_of_x(x)*H_of_x(x))/(H0*H0));
}

double BackgroundCosmology::get_OmegaR(double x) const{
  double a = std::exp(x); 
  if(x == 0.0) return OmegaR;
  else;
    return OmegaR/((std::pow(a, 4)*H_of_x(x)*H_of_x(x))/(H0*H0));
}

double BackgroundCosmology::get_OmegaNu(double x) const{
  double a = std::exp(x); 
  if(x == 0.0) return OmegaNu;
  else;
    return OmegaNu/((std::pow(a, 4)*H_of_x(x)*H_of_x(x))/(H0*H0));
}

double BackgroundCosmology::get_OmegaCDM(double x) const{
  double a = std::exp(x); 
  if(x == 0.0) return OmegaCDM;
  else;
    return OmegaCDM/((std::pow(a, 3)*H_of_x(x)*H_of_x(x))/(H0*H0));
}

double BackgroundCosmology::get_OmegaLambda(double x) const{
  double a = std::exp(x); 
  if(x == 0.0) return OmegaLambda;
  else;
    return OmegaLambda/((H_of_x(x)*H_of_x(x))/(H0*H0));
}

double BackgroundCosmology::get_OmegaK(double x) const{
  double a = std::exp(x); 
  if(x == 0.0) return OmegaK;
  else;
    return OmegaK/((std::pow(a, 2)*H_of_x(x)*H_of_x(x))/(H0*H0));
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  double exp_x = std::exp(x);
  //double Omega_k = get_OmegaK(x);
  double chi = get_comoving_distance_of_x(x); 
  double OmegaK_variable = std::sqrt(std::abs(OmegaK))*H0*chi/Constants.c;
  if(OmegaK < 0)
    return chi * (std::sin(OmegaK_variable))/(OmegaK_variable) / exp_x;
  if(OmegaK > 0)
    return chi * (std::sinh(OmegaK_variable))/(OmegaK_variable) / exp_x;
  if(OmegaK == 0)
    return chi/exp_x;
}

double BackgroundCosmology::get_angular_diameter_distance_of_x(double x) const{
  double exp_x = std::exp(x);
  //double Omega_k = get_OmegaK(x);
  double chi = get_comoving_distance_of_x(x); 
  double OmegaK_variable = std::sqrt(std::abs(OmegaK))*H0*chi/Constants.c;
  if(OmegaK < 0)
    return chi * (std::sin(OmegaK_variable))/(OmegaK_variable) * exp_x;
  if(OmegaK > 0)
    return chi * (std::sinh(OmegaK_variable))/(OmegaK_variable) * exp_x;
  if(OmegaK == 0)
    return chi*exp_x;
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  return eta_of_x(0) - eta_of_x(x);
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  double Gyr = std::pow(10., 9.)*365.*24.*60.*60.;
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << "Ages:        " <<                "\n";
  std::cout << "t_RM_eq      " << t_of_x(-8.132)/Gyr << " Gyr \n";
  std::cout << "t_MDE_eq     " << t_of_x(-0.255)/Gyr << " Gyr \n";
  std::cout << "t_accelerate " << t_of_x(-0.485)/Gyr << " Gyr \n";
  std::cout << "t_aou        " << t_of_x(0)/Gyr   << " Gyr \n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  2.;
  const int    n_pts =  1e6;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << t_of_x(x)          << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << get_luminosity_distance_of_x(x)       << " ";
    fp << get_comoving_distance_of_x(x)         << " ";
    fp << get_angular_diameter_distance_of_x(x) << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

