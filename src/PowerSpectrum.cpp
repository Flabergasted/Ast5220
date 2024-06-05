#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
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
void PowerSpectrum::solve(){

  Vector k_array(n_k);
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  for (int ik = 0; ik < n_k; ik++) {
    k_array[ik] = exp(log_k_array[ik]);
  }
  // calls and splines the bessel function
  generate_bessel_function_splines();
  // Line of sight integration to get Theta_ell(k)
  line_of_sight_integration(k_array);
  // Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  //=============================================================================
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================
  int zmax = 30000;
  const int npts = 30000;
  Vector z_array = Utils::linspace(0, zmax, npts);
  Vector j_ell_array(npts);

  // Run through all the ell values
  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    // Construct an array over all z's
    for(int j = 0; j < npts; j++){
      j_ell_array[j] = Utils::j_ell(ell, z_array[j]);
    }
    // Make the j_ell_splines[i] spline
    j_ell_splines[i].create(z_array, j_ell_array, "j_ell");
  }
  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size(), 0.0));
  double dx = 0.05;
  double x_start = Constants.x_start; 
  double x_end = Constants.x_end;
  int npts = abs(x_end-x_start)/dx;
  double eta0 = cosmo->eta_of_x(0);
  double eta_start = cosmo->eta_of_x(x_start);

  for(size_t ik = 0; ik < k_array.size(); ik++){
    //=============================================================================
    // Solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell 
    // values for the given value of k
    //=============================================================================
    double k = k_array[ik];
    for(size_t ell = 0; ell < ells.size(); ell++){
      double s = source_function(x_start, k)*j_ell_splines[ell](k*(eta0-eta_start));
      
      // integrating over values of x
      for(int ix = 1; ix < npts; ix++){
        double x = x_start + ix*dx;
        double eta = cosmo->eta_of_x(x);
        s += 2*source_function(x,k)*j_ell_splines[ell](k*(eta0 - eta));
      }
      // Store the result for Source_ell(k) in results[ell][ik]
      result[ell][ik] = dx/2.*s;
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100; //?
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);
  // Spline the result and store it in thetaT_ell_of_k_spline
  for(int ell = 0; ell < nells; ell++){
    thetaT_ell_of_k_spline[ell].create(k_array,thetaT_ell_of_k[ell],"line of sight result");
  }
  //if(Constants.polarization){
  //}
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();
  // Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k

  // define the k-array
  Vector result(ells.size(), 0.0); 
  Vector k_array = exp(log_k_array);
  double k_start = k_array[0];
  double k_end = k_array[k_array.size() - 1];

  // define step size and P(k) values for the first step of the trapazoidal integration
  double dk = abs(k_end-k_start)/(n_k);
  double P_start = primordial_power_spectrum(k_start);
  double P_end = primordial_power_spectrum(k_end);
  // loop over the ell values
  for(size_t ell = 0; ell < ells.size(); ell++){
    // Set the initial value of the integral
    double temp = 4.*M_PI*P_start*f_ell_spline[ell](k_start)*g_ell_spline[ell](k_start)/k_start
                + 4.*M_PI*P_end*f_ell_spline[ell](k_end)*g_ell_spline[ell](k_end)/k_end;
    // Loop over the remaining values of k                
    for(int i = 1; i < (n_k); i++){
      double k = k_start + i*dk;
      double P = primordial_power_spectrum(k);
      temp += 2.*4.*M_PI*P*f_ell_spline[ell](k)*g_ell_spline[ell](k)/k;
    }
    // Perform the final step of the trapezoidal solution on before defining the result
    result[ell] = temp*dk/2.;
  }
  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{

  // Defining constants and imports from previous milestones
  double c = Constants.c;
  double Mpc = Constants.Mpc;
  double a = exp(x);
  double H0 = cosmo->get_H0();
  double OmegaCDM = cosmo->get_OmegaCDM();
  double OmegaB = cosmo->get_OmegaB();
  double OmegaM = OmegaB + OmegaCDM;
  double Phi = pert->get_Phi(x, k_mpc);
  // Define the matter overdensity
  double Deta_m = (c*c*k_mpc*k_mpc*Phi)/(3./2.*OmegaM*1./a*H0*H0);
  // Implement the formula of the primordial power spectrum
  double pofk = abs(Deta_m*Deta_m)*primordial_power_spectrum(k_mpc)*(2*M_PI*M_PI)/(pow(k_mpc, 3.));
  // normalise units
  pofk *= pow(cosmo->get_h()/Mpc, 3.);

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    /*
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    */
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}
//====================================================
// Output the k-dependent values
//====================================================

void PowerSpectrum::output_Pk(std::string filename) const {
  std::ofstream fp(filename.c_str());
  double Mpc = Constants.Mpc;
  const int npts       = 60000;
  
  Vector k_mpc(npts);
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), npts);
  for (int ik = 0; ik < npts; ik++) {
    k_mpc[ik] = exp(log_k_array[ik]);
  }
  auto print_Pk_data = [&] (const double k) {
    fp << k                               << " ";
    fp << get_matter_power_spectrum(0, k) << " ";
    fp << thetaT_ell_of_k_spline[0](k)    << " ";
    fp << thetaT_ell_of_k_spline[3](k)    << " ";
    fp << thetaT_ell_of_k_spline[14](k)   << " ";
    fp << thetaT_ell_of_k_spline[19](k)   << " ";
    fp << thetaT_ell_of_k_spline[32](k)   << " ";
    fp << "\n";
  };
  std::for_each(k_mpc.begin(), k_mpc.end(), print_Pk_data);
}
