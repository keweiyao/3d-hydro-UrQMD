#include "trancoeff.h"
#include "eos.h"
#include "inc.h"
#include <iostream>

double Ta = 0.995*0.18, Tb = 1.05*0.18;
double Tc = 0.154;
double C1 = 0.03, C2 = 0.001;
double A0 = -13.45, A1 = 27.55, A2 = -13.77;
double sigma1 = 0.0025, sigma2 = 0.022, sigma3 = 0.025, sigma4 = 0.13;
double lambda1 = 0.9, lambda2 = 0.22, lambda3 = 0.9, lambda4 = 0.25;


TransportCoeff::TransportCoeff(double _etaS_min, double _etaS_slope_QGP, double _etaS_slope_HRG, double _zetaS, EoS *_eos)
{
 etaS_min = _etaS_min ;
 etaS_slope_QGP = _etaS_slope_QGP ;
 etaS_slope_HRG = _etaS_slope_HRG ;
 zetaS = _zetaS ;
 eos = _eos ; 
}

void TransportCoeff::getEta(double e, double T, double &_etaS, double &_zetaS)
{
  //zeta/s
  //  double x = T/0.18;
  // if(T<Ta){ _zetaS = zetaS * (C1 + lambda1*exp((x-1)/sigma1) + lambda2*exp((x-1)/sigma2)) ;}
  //else if(T>Tb){ _zetaS = zetaS * (C2 + lambda3*exp(-(x-1)/sigma3) + lambda4*exp(-(x-1)/sigma4)) ;}
  //else {_zetaS = zetaS * (A0 + A1*x + A2*x*x) ;}
  _zetaS = 0.0;  
  // eta/s
  if (T >= Tc)
    {
      _etaS = etaS_min + (T-Tc)*etaS_slope_QGP;
    }
  if (T < Tc)
    {
      _etaS = etaS_slope_HRG ;
    }

}


void TransportCoeff::getTau(double T, double &_taupi, double &_tauPi)
{
        double etaS; 
	if (T >= 0.154)
	  etaS = etaS_min + (T-0.154)*etaS_slope_QGP;
	else
	  etaS = etaS_slope_HRG;
	//if(T>0.) _taupi=std::max(5./5.068*etaS/T,0.003) ; else _taupi=0.1 ;
	//if(T>0.) _tauPi=std::max(5./5.068*(0.25/C_PI)/T,0.005) ; else _tauPi=0.1 ;
	if(T>0.) 
	  {
	    _taupi=std::max(0.98658*etaS/T,0.003);
	    _tauPi=std::max(0.078509/T,0.005); 
	  }
	else 
	  {
	    _taupi=0.1 ;
	    _tauPi=0.1 ;
	  }
}
