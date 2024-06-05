#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){  
  // Compute and spline Xe, ne
  solve_number_density_electrons(); 
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  Vector x_array(npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);
  Vector Saha_Xe(npts_rec_arrays);

  double peebles_pts;
  double peeblesXe;
  double xi;

  x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Adding the Saha regime to Xe that will remain here, for later comparison.
    Saha_Xe[i] = Xe_current;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;

    } 
    else {
      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      Vector x_peebles_array = {x_array[i-1], x_array[i]};
      double Xeini = Xe_arr[i-1];
      Vector Xe_ic = {Xeini};
      peebles_Xe_ode.solve(dXedx, x_peebles_array, Xe_ic);
      auto Xe_peebles = peebles_Xe_ode.get_data_by_component(0);
      Xe_arr[i] = Xe_peebles[1];

      const double a           = exp(x_array[i]);
 
      // Physical constants
      const double k_b         = Constants.k_b;
      const double G           = Constants.G;
      const double c           = Constants.c;
      const double m_e         = Constants.m_e;
      const double hbar        = Constants.hbar;
      const double m_H         = Constants.m_H;
      const double epsilon_0   = Constants.epsilon_0;

      // Fetch cosmological parameters
      const double TCMB = cosmo->get_TCMB();
      const double OmegaB = cosmo->get_OmegaB();
      const double H0 = cosmo->get_H0();
      double rho_C0 = 3.*H0*H0/(8.*M_PI*G);
      double n_b = OmegaB*rho_C0/(m_H*std::pow(a, 3));

      ne_arr[i] = Xe_arr[i]/n_b;
      }
    }
  Vector Xe_log_arr(npts_rec_arrays);
  Vector Saha_Xe_log_arr(npts_rec_arrays);
  for (int i = 0; i < npts_rec_arrays; i++){
    Xe_log_arr[i] = std::log(Xe_arr[i]);
    Saha_Xe_log_arr[i] = std::log(Saha_Xe[i]);
    if (std::log(Saha_Xe[i]) < -350.){
      Saha_Xe_log_arr[i] = -350.;
    }
  }
  log_Xe_of_x_spline.create(x_array, Xe_log_arr, "Function Xe");
  log_Saha_Xe_of_x_spline.create(x_array, Saha_Xe_log_arr, "Saha Xe");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double TCMB = cosmo->get_TCMB();
  const double OmegaB = cosmo->get_OmegaB();
  const double H0 = cosmo->get_H0();

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  double rho_C0 = 3.*H0*H0/(8.*M_PI*G);
  double n_b = OmegaB*rho_C0/(m_H*std::pow(a, 3));
  double Tb = TCMB/a;
  double Const1 = 1./n_b*std::pow(m_e*(c*c/hbar)*Tb*(k_b/hbar)/(2*M_PI),3./2.);   //m^3/s^3
  double Const2 = std::exp(-epsilon_0/(Tb*k_b));
  double Const = Const1*Const2*1/std::pow(c,3.);
  double A = 1.;
  double B = Const;
  double C = -B;
  if (Const > 4*1e8){
    Xe = 1.;
  }
  else{
    Xe = (-B + std::sqrt(B*B - 4.*A*C))/(2.*A);
    if (Xe != Xe)
    exit(0);
  }
  ne = Xe*n_b;
  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;                                     // J/K
  const double G           = Constants.G;                                       // m^3/(kg s^2)
  const double c           = Constants.c;                                       // m/s
  const double m_e         = Constants.m_e;                                     // kg
  const double hbar        = Constants.hbar;                                    // J/s
  const double m_H         = Constants.m_H;                                     // kg
  const double sigma_T     = Constants.sigma_T;                                 // m^2
  const double lambda_2s1s = Constants.lambda_2s1s;                             // 1/s
  const double epsilon_0   = Constants.epsilon_0;                               // J

  // Cosmological parameters
  const double OmegaB = cosmo->get_OmegaB();
  const double H0 = cosmo->get_H0();                                            // 1/s
  const double TCMB = cosmo->get_TCMB();                                        // K
  const double H = cosmo->H_of_x(x);                                            // 1/s
  double Tb = TCMB/a;                                                           // K
  double alpha = 1./137.0359992;                                                // kg^2/J^2
  double Yp = 0.;

  double phi2 = 0.448*std::log(epsilon_0/(Tb*k_b));                               // log(J/K) -> dim.less
  double alpha2 = 64.*M_PI/std::sqrt(27.*M_PI)*alpha*alpha/(m_e*m_e)
                  *(hbar*hbar/c)*std::sqrt(epsilon_0/(Tb*k_b))*phi2; // 1/J -> m^3/s
  double beta = alpha2*std::pow(m_e*Tb*k_b/(2.*M_PI)*(c*c),3./2.)
                  *std::pow(1./(hbar*c), 3.)*std::exp(-epsilon_0/(Tb*k_b));        // exp(J/K) m^3/s (kgJ)^3/2 -> 1/s
  double beta2 = 0.;
  if (epsilon_0/(Tb*k_b) < 100.){                                                  
    beta2 = beta*std::exp(3.*epsilon_0/(4.*Tb*k_b));                              // exp(J/K) 1/s -> 1/s
  }
  double nH = (1.-Yp)*(3.*H0*H0*OmegaB)/(8.*M_PI*G*m_H*std::pow(a,3.));            // 1/m^3 (no change needed)
  double n1s = (1. - X_e)*nH;                                                     // 1/m^3 (no change needed)
  double lambda_alpha = H*std::pow(3.*epsilon_0,3.)/(std::pow(8.*M_PI,2.)*n1s)
                  *std::pow(1/(hbar*c), 3.);                                      // 1/s J^3 m^3 -> 1/s
  double Cr = (lambda_2s1s+lambda_alpha)/(lambda_2s1s+lambda_alpha+beta2);        // dim.less (no change needed)
  dXedx[0] = Cr/H*(beta*(1.-X_e)-nH*alpha2*X_e*X_e);
  

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);
  Vector x_reverse_array = Utils::linspace(x_end, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    const double a           = exp(x);
    const double c           = Constants.c;
    const double sigma_T     = Constants.sigma_T;
    const double H = cosmo->H_of_x(x);

    // Set the derivative for photon optical depth
    dtaudx[0] = -(c*ne_of_x(x)*sigma_T/H);

    return GSL_SUCCESS;
  };
  const double a_ini       = exp(x_start);
  const double c           = Constants.c;
  const double sigma_T     = Constants.sigma_T;
  const double H_ic     = cosmo->H_of_x(x_end);
  
  double tauini = 0; //-(c*ne_of_x(x_end)*sigma_T/H_ic);
  Vector tau_ic {tauini};
  
  ODESolver ode;
  ode.solve(dtaudx, x_reverse_array, tau_ic);

  auto tau_array = ode.get_data_by_component(0);

  // Reverse the tau array
  std::reverse(tau_array.begin(),tau_array.end());
  // create a tau spline
  tau_of_x_spline.create(x_array, tau_array, "tau_of_x");
  //Visibility function
  Vector g_tilde_array(npts);
  for (int i = 0; i < npts; i++){
    g_tilde_array[i] = -dtaudx_of_x(x_array[i])*std::exp(-tau_of_x(x_array[i]));
  }

  // create g_tilde spline
  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g_tilde");

  Utils::EndTiming("opticaldepth");

  double Omegab0  = cosmo->get_OmegaB();
  double OmegaR0  = cosmo->get_OmegaR();
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
    double a = exp(x);
    double Hp_of_x  = cosmo->Hp_of_x(x);
    double R = 4.*OmegaR0/(3.*Omegab0*a);
    double cs = c*std::sqrt(R/(3.*(1+R)));
    dsdx[0] = cs/Hp_of_x;
    return GSL_SUCCESS;
  };

double R_ini = 4.*OmegaR0/(3.*Omegab0*a_ini);
double Hp_ini = cosmo->Hp_of_x(x_start);
double cs_ini = c*std::sqrt(R_ini/(3.*(1+R_ini)));

double sini = cs_ini/Hp_ini;
Vector s_ic {sini};

// Solve the ode
ode.solve(dsdx, x_array, s_ic);

auto s_array = ode.get_data_by_component(0);

// create the spline
s_of_x_spline.create(x_array, s_array, "s_of_x");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return std::exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::Saha_Xe_of_x(double x) const{
  return std::exp(log_Saha_Xe_of_x_spline(x));
}

double RecombinationHistory::s_of_x(double x) const{
  return s_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  const double a           = exp(x);
  const double G           = Constants.G;
  const double m_H         = Constants.m_H;
  const double OmegaB = cosmo->get_OmegaB();
  const double H0 = cosmo->get_H0();
  
  double rho_C0 = 3.*H0*H0/(8.*M_PI*G);
  double n_b = OmegaB*rho_C0/(m_H*std::pow(a, 3.));
  return Xe_of_x(x)*n_b;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 60000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << cosmo->t_of_x(x)     << " ";
    fp << s_of_x(x)            << " ";
    fp << Saha_Xe_of_x(x)      << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

