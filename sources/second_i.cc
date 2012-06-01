static double K0_tau_, ES_, E_;

static double integral_cal_In_integrand(double rS)
    {
double C_radial,I_radial,C_angular, weight_SIMSPON,dzS,zS,I_V_rS,R,thetaS_max,
       dthetaS,thetaS,thetaS_inf,I_angular_rS;
int i_zS,n_thetaS,i_thetaS;    
int      n_zS   = 2 * 10;
double      zS_max   = sqrt(4.0 * K0_tau_ * log(10.0) * DIGIT);
      if (zS_max >= L_dif){zS_max = L_dif;}
      dzS    = zS_max / ((double)n_zS);
      zS     = 0.0;
      I_V_rS = 0.0; 

      for (i_zS=0;i_zS<=n_zS;i_zS++)
	{
	  if (i_zS==0 || i_zS==n_zS) {weight_SIMSPON = 1./3.;}
	  else {weight_SIMSPON = (1. + (double)(i_zS % 2)) * 2. / 3.;}
	  R = sqrt(rS*rS + zS*zS);
//	  rho=Rho(R);

	  I_V_rS += weight_SIMSPON * vertical_positron(zS,ES_,0.,E_)*(rhoQ_(R)/Rhosun/Rhosun)*dzS;
	  zS     += dzS;
	}
      I_V_rS *= 2.0;
      n_thetaS		 = 2 * 10;
      C_angular		 = Rsun * rS / (2.0 * K0_tau_);
      thetaS_max   = sqrt((DIGIT*log(10.0)) / (2.0*C_angular));
      if (thetaS_max >= 1.0){thetaS_max = 1.0;}
      thetaS_max   = 2.0 * asin(thetaS_max);
      dthetaS      = thetaS_max / ((double)n_thetaS);
      thetaS       = 0.0;
      thetaS_inf   = 0.0;
      I_angular_rS = 0.0;
      for (i_thetaS=0;i_thetaS<=n_thetaS;i_thetaS++)
	{
	  if (i_thetaS==0 || i_thetaS==n_thetaS) {weight_SIMSPON = 1./3.;}
	  else {weight_SIMSPON = (1. + (double)(i_thetaS % 2)) * 2. / 3.;}
	  I_angular_rS += weight_SIMSPON * exp(C_angular*(cos(thetaS) - 1.0)) * dthetaS;
	  thetaS       += dthetaS;
	}
      I_angular_rS *= 2.0;

      C_radial = (Rsun - rS)*(Rsun - rS) / (4.0 * K0_tau_);
      C_radial = exp((-1.0)*C_radial) / (4.0 * M_PI * K0_tau_);
      
      return  C_radial * I_V_rS * I_angular_rS * rS ;
        
    }


double integral_cal_In(double ES,double E)
{
  double epsilon_t;
  double rho,weight_SIMSPON,RHO_sun,Ediff;
  long i_zS,n_zS,i_thetaS,n_thetaS,i_rS,n_rS;
  double z,t,zS,tS,tau,K0_tau;
  double dzS,zS_max;
  double thetaS,dthetaS,thetaS_max,thetaS_inf;
  double rS,drS,rS_min,rS_max,R;
  double I_V_rS,C_angular,I_angular_rS,C_radial,I_radial;

  RHO_sun=Rhosun;
  Ediff=L_dif;
  
  tS        = tau_dif * function_v(ES);	     /* [sec] */
  t         = tau_dif * function_v(E);             /* [sec] */
  epsilon_t = (0.1 * (3600. * 24. * 365.25));			     /* [sec] */
  tau       = (t - tS) + epsilon_t;                                  /* [sec] */
  K0_tau    = K_dif * pow(CM_PER_KPC,-2.) * tau;   /* [kpc^{2}] */
  
  z         = 0.0;
  
  n_rS		 = 2 * 40;
  zS_max   = sqrt(4.0 * K0_tau * log(10.0) * DIGIT);
  rS_min   = Rsun - zS_max;
  rS_max   = Rsun + zS_max;
  if (rS_min <= 0.0  ){rS_min = 0.0;  }
  if (rS_max >= R_GAL){rS_max = R_GAL;}
  drS			 = (rS_max - rS_min) / ((double)n_rS);
  rS			 = rS_min;

  ES_=ES;  K0_tau_=K0_tau; E_=E;
  I_radial = simpson(integral_cal_In_integrand,rS_min,rS_max,1.E-3);

  return I_radial;
}
