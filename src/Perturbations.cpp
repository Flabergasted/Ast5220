#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  const double OmegaR  = cosmo->get_OmegaR();
  Vector k_array(n_k);
  Vector x_arr = Utils::linspace(x_start, x_end, n_x);
  // Set up a logaritmic linspace
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  // Loop through and perform the exponential of each wich should give logarithmic spacing.
  for (int ik = 0; ik < n_k; ik++) {
    k_array[ik] = exp(log_k_array[ik]);
  }

  // Allocate the vectors for each of the values we want to spline.
  //Vector2D ;
  Vector2D delta_cdm_arr(n_x,Vector(n_k, 0.0));
  Vector2D delta_b_arr(n_x,Vector(n_k, 0.0));
  Vector2D v_cdm_arr(n_x,Vector(n_k, 0.0));
  Vector2D v_b_arr(n_x,Vector(n_k, 0.0));
  Vector2D Phi_arr(n_x,Vector(n_k, 0.0));
  Vector2D Psi_arr(n_x,Vector(n_k, 0.0));
  Vector2D Theta_0_arr(n_x,Vector(n_k, 0.0));
  Vector2D Theta_1_arr(n_x,Vector(n_k, 0.0));
  Vector2D Theta_2_arr(n_x,Vector(n_k, 0.0));

  Vector2D Theta_3_arr(n_x,Vector(n_k, 0.0));

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);
    double x_low = std::lower_bound(x_arr.begin(), x_arr.end(), x_end_tight)- x_arr.begin();
    Vector x_arr_tight = Utils::linspace(x_start, x_arr[x_low], x_low);
    Vector x_arr_full = Utils::linspace(x_arr[x_low],x_end, n_x-x_low);

    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    // ...
    ODESolver ode;
    Vector y_tight_ic {y_tight_coupling_ini};
    ode.solve(dydx_tight_coupling, x_arr_tight, y_tight_ic);
    auto y_tight = ode.get_data();
    auto y_tight_coupling = ode.get_final_data();

    //====i===============================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================
    
    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);
    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    // ...
    Vector y_full_ic {y_full_ini};
    ode.solve(dydx_full, x_arr_full, y_full_ic);
    auto y_full = ode.get_data();

    for(int ix = 0; ix < n_x; ix++){
      // Outside tight coupling
      if (ix < x_low){
        double x = x_arr_tight[ix];
        double a = exp(x);
        auto y = y_tight[ix];
        delta_cdm_arr[ix][ik] = y[0];

        delta_b_arr[ix][ik]   = y[1];
        v_cdm_arr[ix][ik]     = y[2];
        v_b_arr[ix][ik]       = y[3];
        Phi_arr[ix][ik]       = y[4];
        Theta_0_arr[ix][ik]   = y[5];
        Theta_1_arr[ix][ik]   = y[6];
        Theta_2_arr[ix][ik]   = -20.*Constants.c*k/(45.*cosmo->Hp_of_x(x)*rec->dtaudx_of_x(x))
                                *Theta_1_arr[ix][ik];
        Psi_arr[ix][ik]       = -Phi_arr[ix][ik]-12.*cosmo->get_H0()*cosmo->get_H0()
                                /(std::pow(Constants.c*k*a,2.))
                                *(OmegaR*Theta_2_arr[ix][ik]);                      
      }
      else {
        double a = exp(x_arr_full[ix-x_low]);
        auto y = y_full[ix-x_low];
        delta_cdm_arr[ix][ik] = y[0];
        delta_b_arr[ix][ik]   = y[1];
        v_cdm_arr[ix][ik]     = y[2];
        v_b_arr[ix][ik]       = y[3];
        Phi_arr[ix][ik]       = y[4];
        Theta_0_arr[ix][ik]   = y[5];
        Theta_1_arr[ix][ik]   = y[6];
        Theta_2_arr[ix][ik]   = y[7];
        Psi_arr[ix][ik]       = -Phi_arr[ix][ik]-12.*cosmo->get_H0()*cosmo->get_H0()
                                /(std::pow(Constants.c*k*a,2.))
                                *(OmegaR*Theta_2_arr[ix][ik]);
      }
    }
  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  delta_cdm_spline.create(x_arr,k_array,delta_cdm_arr);
  delta_b_spline.create(x_arr,k_array,delta_b_arr);
  v_cdm_spline.create(x_arr,k_array,v_cdm_arr);
  v_b_spline.create(x_arr,k_array,v_b_arr);
  Phi_spline.create(x_arr,k_array,Phi_arr);
  Psi_spline.create(x_arr,k_array, Psi_arr);
  Spline2D Theta_0_spline;
  Theta_0_spline.create(x_arr,k_array,Theta_0_arr);
  Spline2D Theta_1_spline;
  Theta_1_spline.create(x_arr,k_array,Theta_1_arr);
  Spline2D Theta_2_spline;
  Theta_2_spline.create(x_arr,k_array,Theta_2_arr);
  Theta_spline = {Theta_0_spline, Theta_1_spline, Theta_2_spline};
  Pi_spline.create(x_arr,k_array,Theta_2_arr);
  // ...
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================
  // ...
  // ...
  const double c            = Constants.c;
  const double Hp_of_x      = cosmo->Hp_of_x(x);
  //const double dHp_of_x     = cosmo->dHpdx_of_x(x);
  //const double H0           = cosmo->get_H0();
  //const double OmegaCDM     = cosmo->get_OmegaCDM();
  //const double OmegaB       = cosmo->get_OmegaB();
  //const double OmegaR  = cosmo->get_OmegaR();
  //const double dtaudx       = rec->dtaudx_of_x(x);
  //const double ddtauddx     = rec->ddtauddx_of_x(x);
  //double a = exp(x);

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  // ...
  // ...
  double Psi = -2./3.;
  Phi = -Psi;
  delta_b = -3./2.*Psi;
  delta_cdm = delta_b;
  v_b = -c*k/(2.*Hp_of_x)*Psi;
  v_cdm = v_b;

  // SET: Photon temperature perturbations (Theta_ell)
  // ...
  Theta[0] = -1./2.*Psi;
  Theta[1] = c*k/(6*Hp_of_x)*Psi;
  //Theta[2] = -20.*c*k/(45.*Hp_of_x*dtaudx)*Theta[1];

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // ...
    // ...
  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  // ...
  const double c            = Constants.c;
  const double Hp_of_x      = cosmo->Hp_of_x(x);
  const double dHp_of_x     = cosmo->dHpdx_of_x(x);
  const double H0           = cosmo->get_H0();
  const double OmegaCDM     = cosmo->get_OmegaCDM();
  const double OmegaB       = cosmo->get_OmegaB();
  const double OmegaR  = cosmo->get_OmegaR();
  const double dtaudx       = rec->dtaudx_of_x(x);
  const double ddtauddx     = rec->ddtauddx_of_x(x);
  double a = exp(x);

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  // ...
  //double Psi = -2./3.; //?
  Phi = Phi_tc; //-Psi;
  delta_b = delta_b_tc; //-3./2.*Psi;
  delta_cdm = delta_cdm_tc; //delta_b;
  v_b = v_b_tc; //-c*k/(2.*Hp_of_x)*Psi;
  v_cdm = v_cdm_tc; //v_b;

  // SET: Photon temperature perturbations (Theta_ell)
  // ...
  int l_max = n_ell_theta-1;
  Theta[0] = Theta_tc[0];//-1./2.*Psi;
  Theta[1] = Theta_tc[1];//c*k/(6*Hp_of_x)*Psi;
  Theta[2] = -20.*c*k/(45.*Hp_of_x*dtaudx)*Theta[1];
   
  for (int l = 3; l <= l_max; l++){
    Theta[l] = -l/(2*l + 1.)*c*k/(Hp_of_x*dtaudx)*Theta[l-1];
  }

  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    // ...
    // ...
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // ...
    // ...
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end = 0.0;

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  // ...
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  double Hp_of_x;
  double dtaudx;
  double x;
  double c = Constants.c;

  for (int ix = 0; ix < n_x; ix ++){
    x = x_array[ix];
    dtaudx = rec->dtaudx_of_x(x);
    Hp_of_x = cosmo->Hp_of_x(x);
    x_tight_coupling_end = x_array[ix];
    if (x > -8.3 || abs(dtaudx)<10. || abs(c*k/(Hp_of_x*dtaudx))>1./10.){
      break;
    }
  }
  return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  //Vector k_array;
  //Vector x_array;
  Vector k_array(n_k);
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  const double c = Constants.c;
  // Set up a logaritmic linspace
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  // Loop through and perform the exponential of each wich should give logarithmic spacing.
  for (int ik = 0; ik < n_k; ik++) {
    k_array[ik] = exp(log_k_array[ik]);
  }

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      const double Hp         = cosmo->Hp_of_x(x);
      const double dHp        = cosmo->dHpdx_of_x(x);
      const double ddHp       = cosmo->ddHpddx_of_x(x);
      const double tau        = rec->tau_of_x(x);
      const double dtau       = rec->dtaudx_of_x(x);
      const double ddtau      = rec->ddtauddx_of_x(x);
      const double g_tilde    = rec->g_tilde_of_x(x);
      const double dg_tilde   = rec->dgdx_tilde_of_x(x);
      const double ddg_tilde  = rec->ddgddx_tilde_of_x(x);

      double Theta_0          = Theta_spline[0](x, k);
      double Theta_1          = Theta_spline[1](x, k);
      double dTheta_1         = Theta_spline[1].deriv_x(x, k);
      //double Theta_3          = Theta_spline[3](x, k);
      //double dTheta_3         = Theta_spline[3].deriv_x(x, k);
      double PI               = Theta_spline[2](x, k);
      double dPI              = Theta_spline[2].deriv_x(x, k);
      double ddPI             = Theta_spline[2].deriv_xx(x, k);
      double Psi              = Psi_spline(x, k);
      double Phi              = Phi_spline(x, k);
      double dPsi             = Psi_spline.deriv_x(x, k);
      double dPhi             = Phi_spline.deriv_x(x, k);
      double v_b              = v_b_spline(x, k);
      double dv_b             = v_b_spline.deriv_x(x, k);
      
      //double ddPI = 2.*k/(5.*Hp)*(-dHp/Hp*Theta_1 + dTheta_1) + 3./10.*(ddtau*PI + dtau*dPI) 
      //            - 3.*k/(5.*Hp)*(-dHp/Hp*Theta_3 + dTheta_3);
      double dHpgv = (dHp*g_tilde*v_b+Hp*dg_tilde*v_b+Hp*g_tilde*dv_b);
      double dHdHgP = g_tilde*PI*(dHp*dHp+ddHp*Hp) + 3.*Hp*dHp*(dg_tilde*PI+g_tilde*dPI)
                    + Hp*Hp*(ddg_tilde*PI + 2.*dg_tilde*dPI+g_tilde*ddPI);
      // ...

      // Temperatur source
      //ST_array[index] = 0.0;
      ST_array[index] = g_tilde * (Theta_0+Psi+1./4.*PI)+exp(-tau)*(dPsi-dPhi)-1./(c*k)*dHpgv
                      + 3./(4.*c*c*k*k)*dHdHgP;

      // Polarization source
      //if(Constants.polarization){
      //  SE_array[index] = 0.0;
      //}
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  //if(Constants.polarization){
  //  SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  //}

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  const double c            = Constants.c;
  const double Hp_of_x      = cosmo->Hp_of_x(x);
  const double dHp_of_x     = cosmo->dHpdx_of_x(x);
  const double H0           = cosmo->get_H0();
  const double OmegaCDM     = cosmo->get_OmegaCDM();
  const double OmegaB       = cosmo->get_OmegaB();
  const double OmegaR       = cosmo->get_OmegaR();
  const double dtaudx       = rec->dtaudx_of_x(x);
  const double ddtauddx     = rec->ddtauddx_of_x(x);
  double a = exp(x);

  // SET: Scalar quantities (Phi, delta, v, ...)
  double Theta_2 = -20.*c*k/(45.*Hp_of_x*dtaudx)*Theta[1];
  double Psi = -Phi - 12.*H0*H0/(std::pow(c*k*a,2.))*(OmegaR*Theta_2);
  dPhidx = Psi - c*c*k*k/(3.*Hp_of_x*Hp_of_x)*Phi + H0*H0/(2.*Hp_of_x*Hp_of_x)
                *(OmegaCDM/a*delta_cdm + OmegaB/a*delta_b + 4.*OmegaR/(a*a)*Theta[0]);
  ddelta_cdmdx = c*k/Hp_of_x*v_cdm - 3.*dPhidx;
  dv_cdmdx = -v_cdm - c*k/Hp_of_x*Psi;
  ddelta_bdx = c*k/Hp_of_x*v_b - 3.*dPhidx;
  
  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -k*c/Hp_of_x * Theta[1] - dPhidx;

  double R = 4.*OmegaR/(3*OmegaB*a);
  double q = (-((1.-R)*dtaudx+(1.+R)*ddtauddx)*(3.*Theta[1]+v_b)-c*k/Hp_of_x*Psi
              +(1.-dHp_of_x/Hp_of_x)*c*k/Hp_of_x*(-Theta[0]+2.*Theta_2)-c*k/Hp_of_x*dThetadx[0])
              /((1.+R)*dtaudx+dHp_of_x/Hp_of_x-1.);
  dv_bdx = 1./(1.+R)*(-v_b-c*k/Hp_of_x*Psi+R*(q+c*k/Hp_of_x*(-Theta[0]+2.*Theta_2)
                -c*k/Hp_of_x*Psi));

  dThetadx[1] = 1./3.*(q-dv_bdx);

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // ...
    // ...
    // ...
  }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  const double Hp_of_x      = cosmo->Hp_of_x(x);
  const double H0           = cosmo->get_H0();
  const double OmegaCDM     = cosmo->get_OmegaCDM();
  const double OmegaB       = cosmo->get_OmegaB();
  const double OmegaR       = cosmo->get_OmegaR();
  const double eta          = cosmo->eta_of_x(x);
  // ...

  // Recombination variables
  const double dtaudx       = rec->dtaudx_of_x(x);

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  const double c            = Constants.c;
  double a = exp(x);

  // SET: Scalar quantities (Phi, delta, v, ...)
  // ...
  double Psi = -Phi - 12.*H0*H0/(c*c*k*k*a*a)*(OmegaR*Theta[2]);
  //std::cout << "-Phi " << -Phi << "\n";
  dPhidx = Psi - c*c*k*k/(3.*Hp_of_x*Hp_of_x)*Phi + H0*H0/(2.*Hp_of_x*Hp_of_x)
                *(OmegaCDM/a*delta_cdm + OmegaB/a*delta_b + 4.*OmegaR/(a*a)*Theta[0]);

  // SET: Photon multipoles (Theta_ell)
  // ...
  //std::cout << "Theta 0 " << Theta[0] << "\n";
  //std::cout << "Theta 1 " << Theta[1] << "\n";
  //std::cout << "Theta 2 " << Theta[2] << "\n";
  double PI = Theta[2];
  int l_max = n_ell_theta-1;
  //std::cout << "" << << "\n";
  dThetadx[l_max] = c*k/Hp_of_x*Theta[l_max-1] 
                        - c*(l_max+1.)/(Hp_of_x*eta)*Theta[l_max]
                        + dtaudx*Theta[l_max];
  //std::cout << "" << << "\n";
  //*
  for (int l = 3; l < l_max; l++){
    dThetadx[l] = l*c*k/((2.*l+1.)*Hp_of_x)*Theta[l-1] 
                - (l+1.)*c*k/((2.*l+1.)*Hp_of_x)*Theta[l+1]
                + dtaudx*(Theta[l]);
    //std::cout << "" << << "\n";
  }
  //*/
  //std::cout << "" << << "\n";
  dThetadx[2] = 2.*c*k/((2.*2.+1.)*Hp_of_x)*Theta[2-1] 
                - (2.+1.)*c*k/((2.*2.+1.)*Hp_of_x)*Theta[2+1] 
                + dtaudx*(Theta[2]-1./10.*PI);
  //std::cout << "" << << "\n";
  /*
  for (int l = 2; l < l_max; l++){
    dThetadx[l] = l*c*k/((2.*l+1.)*Hp_of_x)*Theta[l-1] 
                - (l+1.)*c*k/((2.*l+1.)*Hp_of_x)*Theta[l+1]
                + dtaudx*(Theta[l] - 1./10.*PI*(l == 2));
  }
  */
  dThetadx[1] = k*c/(3.*Hp_of_x)*Theta[0] - 2.*c*k/(3.*Hp_of_x)*Theta[2] 
              + c*k/(3.*Hp_of_x)*Psi + dtaudx*(Theta[1]+1./3.*v_b);
  //std::cout << "" << << "\n";
  dThetadx[0] = -k*c/Hp_of_x*Theta[1] - dPhidx;
  // Metric perturbations
  double R = 4.*OmegaR/(3.*OmegaB*a);
  //std::cout << "" << << "\n";

  // CDM and baryons
  
  //These look ok
  //std::cout << "" << << "\n";
  dv_bdx = - v_b - c*k/Hp_of_x*Psi + dtaudx*R*(3.*Theta[1]+v_b);
  ddelta_bdx = c*k/Hp_of_x*v_b - 3.*dPhidx;
  dv_cdmdx = - v_cdm - c*k/Hp_of_x*Psi;
  ddelta_cdmdx = c*k/Hp_of_x*v_cdm - 3.*dPhidx;

  // SET: Photon polarization multipoles (Theta_p_ell)
  if(polarization){
    // ...
    // ...
    // ...
  }

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // ...
    // ...
    // ...
  }

  //return GSL_SUCCESS;
  return 0;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_v_b(x,k)       << " ";

    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

