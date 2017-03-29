#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "eos.h"
#include "inc.h"
//#include <TGraph.h>

using namespace std ;

const double bagp = pow(247.19/197.32,4)/gevtofm ;
const double bagt = pow(247.19/197.32,4)/gevtofm ;
const double deg = 16.0+3.0*12.0*(7.0/8.0);

double EoS::s(double e, double nb, double nq, double ns)
{
	double T, mub, muq, mus, p ;
	eos(e, nb, nq, ns, T, mub, muq, mus, p) ;
	if(T>0.0) return (e+p-mub*nb-muq*nq-mus*ns)/T ;
	else return 0.;
}
