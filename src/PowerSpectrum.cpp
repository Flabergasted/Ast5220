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

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  //Vector k_array;
  //Vector k_array = Utils::linspace(k_min, k_max, n_k);
  Vector k_array(n_k);
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  for (int ik = 0; ik < n_k; ik++) {
    k_array[ik] = exp(log_k_array[ik]);
  }

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
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================
  int zmax = 30000;
  const int npts = 50000;
  Vector z_array = Utils::linspace(0, zmax, npts);
  Vector j_ell_array(zmax);

  // Run through all the ell values
  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    // Construct an array over all z's
    for(int j = 0; j < npts; z++){
      j_ell_array[j] = Utils::j_ell(ell, z_array[j]);
    }

    // Make a spline of the j_ell(z) array
    //Spline j_ell_spline;
    //j_ell_spline.create(z_array, j_ell_array, "j_ell");

    // Make the j_ell_splines[i] spline
    j_ell_splines[i].create(z_array, j_ell_array, "j_ell");//= j_ell_spline;
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
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  Spline jell;
  double eta0 = cosmo->eta_of_x(0);
  int npts = 2000;
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, npts);
  double dx = 1./npts;
  ODESolver ode;

  for(size_t ik = 0; ik < k_array.size(); ik++){
    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    double k = k_array[ik];
    // integrating over values of x
    for(int ix = 0; ix < npts; ix++){
      double x = x_array[ix];
      double eta = cosmo->eta_of_x(x);
      // sums over values of ell
      for(size_t ell = 0; ell < ells.size(); ell++){
        //std::cout << "ik: " << ik << "\n";
        result[ell][ik] += source_function(x, k)*j_ell_splines[ell](k*(eta0-eta))*dx;//(k)*(k*(eta0-eta))*dx;
        
        //auto initial_conditions = source_function(Constants.x_start, k)*j_ell_splines[ell](k*(eta0-cosmo->eta_of_x(Constants.x_start)));
        //auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);
        //ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
        //  return rhs_full_ode(x, k, y, dydx);
        //};
        //ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
        //  dtdx[0] = 1./H_of_x(x);
        //  return GSL_SUCCESS;
        //};
        //Vector y_full_ic {y_full_ini};
        //ode.solve(dydx_full, x_arr_full, y_full_ic);
        //auto y_full = ode.get_data();
      }      
    }
    // ...
    // ...
    // ...
    // Store the result for Source_ell(k) in results[ell][ik]
  }
  std::cout << "source function(0,0)" << source_function(x_array[0], k_array[0]) << "\n";
  std::cout << "j_ell(k=0)" << j_ell_splines[0](k_array[0]) << "\n";
  std::cout << "k[0]" << k_array[0] << "\n";
  //(eta0-eta))*dx

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

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);
  // Spline the result and store it in thetaT_ell_of_k_spline
  for(int ell = 0; ell < nells; ell++){
    Spline result;
    result.create(k_array,thetaT_ell_of_k[ell],"line of sight result");
    thetaT_ell_of_k_spline[ell] = result;
  }

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
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

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================

  Vector result(ells.size()); 

  for(size_t ell = 0; ell < ells.size(); ell++){
    double dk = abs(log_k_array[1]-log_k_array[0]);
    for(int i = 0; i < log_k_array.size(); i++){
      if(i != 0){
        dk = abs(log_k_array[i]-log_k_array[i-1]);
      }
      //std::cout << "ik: " << ik << "\n";
      double k = exp(log_k_array[i]);
      double P = primordial_power_spectrum(k);
      result[ell] += 4.*M_PI*P*f_ell_spline[ell](k)*g_ell_spline[ell](k)*dk;
      //std::cout << "f_ell " << f_ell_spline[ell](k) << "\n";
      //std::cout << "g_ell " << g_ell_spline[ell](k) << "\n";
    }
  }
  /*
  ODEFunction dCelldlogk = [&](double k, const double *Cell, double *dCelldlogk){
    P = A_s * std::pow(k/k_pivot, n_s-1);
    dCelldlogk[0] = 4. * M_PI * P * f_ell_spline(k) * g_ell_spline(k);
    return GSL_SUCCESS;
  }
  double Cellini = 4. * M_PI * P * f_ell_spline(k) * g_ell_spline(k);
  Vector Cell_ic {Cellini};
  ODESolver ode;
  ode.solve(dCelldlogk,,Cell_ic) 
  */
  // The ODE for deta/dx
  /*
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
  */
  // ...
  // ...
  // ...
  // ...

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
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================
  /*
  double c = Constants.c;
  double H0 = cosmo->get_H0();
  double OmegaM = cosmo->get_OmegaM();
  double Phi = pert->get_Phi(x, k_mpc);
  double a = exp(x);
  double Deta_m = (c*c*k_mpc*k_mpc*Phi)/(3./2.*OmegaM*1./a*H0*H0);
  pofk = abs(Deta_m*Deta_m)*primordial_power_spectrum(k_mpc);
  */
  // ...
  // ...
  // ...

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
    //fp << get_matter_power_spectrum(x, k_mpc) << " ";
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

