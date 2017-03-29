#include "trancoeff.h"
#include "eos.h"
#include "inc.h"
#include <iostream>

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
  _zetaS= zetaS*(1. - 3.*eos->cs2(e))/(1. - 3.*eos->cs2(0.8))*std::exp(-pow((T-0.175)/0.01, 2.0));  
  if (T >= 0.154)
    {
      _etaS = etaS_min + (T-0.154)*etaS_slope_QGP;
      return;
    }
  if (T < 0.154)
    {
      _etaS = etaS_slope_HRG ;
      return;
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
