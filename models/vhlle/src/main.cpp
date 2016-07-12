#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <sstream>
#include "fld.h"
#include "hdo.h"
#include "mcgstream.h"
#include "eos.h"
#include "eos_hotqcd.h"
#include "trancoeff.h"

using namespace std ;

// program parameters, to be read from file
int nx, ny, nz, event_number ;
double xmin, xmax, ymin, ymax, etamin, etamax, tau0, tauMax, dtau ;
char outputDir[255];
char icInputFile [255] ;
double etaS_min, etaS_slope_QGP, etaS_slope_HRG, zetaS, eCrit, epsilon0;
void readParameters(char *parFile)
{
	char parName [255], parValue [255] ;
	ifstream fin(parFile) ;
	if(!fin.is_open()) { cout << "cannot open parameters file " << parFile << endl ; exit(1) ; }
	while(fin.good()){
	 string line ;
	 getline(fin, line) ;
	 istringstream sline (line) ;
	 sline >> parName >> parValue ;
	 if     (strcmp(parName,"outputDir")==0) strcpy(outputDir, parValue) ;
	 else if(strcmp(parName,"icInputFile")==0) strcpy(icInputFile, parValue) ;
	 else if(strcmp(parName,"nx")==0) nx = atoi(parValue) ;
	 else if(strcmp(parName,"ny")==0) ny = atoi(parValue) ;
	 else if(strcmp(parName,"nz")==0) nz = atoi(parValue) ;

	 else if(strcmp(parName,"xmin")==0) xmin = atof(parValue) ;
	 else if(strcmp(parName,"xmax")==0) xmax = atof(parValue) ;
	 else if(strcmp(parName,"ymin")==0) ymin = atof(parValue) ;
	 else if(strcmp(parName,"ymax")==0) ymax = atof(parValue) ;
	 else if(strcmp(parName,"etamin")==0) etamin = atof(parValue) ;
	 else if(strcmp(parName,"etamax")==0) etamax = atof(parValue) ;
	 else if(strcmp(parName,"tau0")==0) tau0 = atof(parValue) ;
	 else if(strcmp(parName,"tauMax")==0) tauMax = atof(parValue) ;
	 else if(strcmp(parName,"dtau")==0) dtau = atof(parValue) ;

         else if(strcmp(parName,"e_crit")==0) eCrit = atof(parValue);
	 else if(strcmp(parName,"etaS_min")==0) etaS_min = atof(parValue) ;
	 else if(strcmp(parName,"etaS_slope_QGP")==0) etaS_slope_QGP = atof(parValue) ;
	 else if(strcmp(parName,"etaS_slope_HRG")==0) etaS_slope_HRG = atof(parValue) ;
	 else if(strcmp(parName,"zetaS")==0) zetaS = atof(parValue) ;

	 else if(strcmp(parName,"epsilon0")==0) epsilon0 = atof(parValue) ;
	 else if(strcmp(parName,"event_number")==0) event_number = atoi(parValue) ;
	 else if(parName[0]=='!') cout << "CCC " << sline.str() << endl ;
	 else cout << "UUU " << sline.str() << endl ;
	}
}


void printParameters()
{
  cout << "====== parameters ======\n" ;
  cout << "outputDir = " << outputDir << endl ;
  cout << "icInputFile = " << icInputFile << endl ;
  cout << "nx = " <<  nx << endl ;
  cout << "ny = " << ny << endl ;
  cout << "nz = " << nz << endl ;
  cout << "xmin = " << xmin << endl ;
  cout << "xmax = " << xmax << endl ;
  cout << "ymin = " << ymin << endl ;
  cout << "ymax = " << ymax << endl ;
  cout << "etamin = " << etamin << endl ;
  cout << "etamax = " << etamax << endl ;
  cout << "tau0 = " << tau0 << endl ;
  cout << "tauMax = " << tauMax << endl ;
  cout << "dtau = " << dtau << endl ;
  cout << "e_crit = "<< eCrit << endl;
  cout << "eta/s = " << etaS_min << "+ (T-Tc)*(T>Tc)*" << etaS_slope_QGP << " + (Tc-T)*(T<Tc)*" << etaS_slope_HRG  << endl ;
  cout << "zeta/s = " << zetaS << endl ;
  cout << "epsilon0 = " << epsilon0 << endl ;
  cout << "event number = " << event_number << endl ;
  cout << "======= end parameters =======\n" ;
}

int main(int argc, char **argv)
{
  // pointers to all the main objects
  EoS *eos ;
  TransportCoeff *trcoeff ;
  Fluid *f ;
  Hydro* h ;
  time_t start=0, end ;
  
  time(&start);
  
  // read parameters from file
  char* parFile ;
  if(argc==1)
  {
     cout << "NO PARAMETERS, exiting\n" ;
     exit(1) ;
  }
  else parFile = argv[1] ;
  readParameters(parFile) ;
  printParameters() ;
 
   // eos
  eos = new Eos_hotqcd();
  
  // transport coefficients
  trcoeff = new TransportCoeff(etaS_min, etaS_slope_QGP, etaS_slope_HRG, zetaS, eos) ;
  
  // fluid
  f = new Fluid(eos, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax, etamin, etamax, dtau, eCrit) ;
  
  // initilal conditions (read from file)
  MCGstreamIC *ic = new MCGstreamIC(epsilon0, icInputFile, eos, event_number);
  ic->setIC(f, tau0);
   
  // hydro init
  h = new Hydro(f, eos, trcoeff, tau0, dtau) ;
  int maxstep = ceil((tauMax-tau0)/dtau) ;
  start = 0;
  time(&start);
  //h->setNSvalues() ; // initialize viscous terms
  h->setQfull() ; // set Qfull in each cell, in order to output IC correctly

  // hllev321v1 = with pre-advection
  f->initOutput(outputDir, maxstep, tau0, 2) ;
  f->outputCorona(tau0) ;

  for(int istep=0; istep<maxstep; istep++)
  {
    // decrease timestep automatically, but use fixed dtau for output
    int nSubSteps = 1 ;
    while(dtau/nSubSteps>0.5*(tau0+dtau*istep))
      {
	nSubSteps *= 2 ; // 0.02 in "old" coordinates
      }
    h->setDtau(dtau/nSubSteps) ;
    for(int j=0; j<nSubSteps; j++)
      { 
	h->performStep() ;
      }
    // cout << tau0 + istep*dtau << endl;
        f->outputSurface(h->getTau()) ;
  }

  end=0 ;
  time(&end); 
  float diff2 = difftime(end, start);
  cout<<"Execution time = "<<diff2<< " [sec]" << endl;
   	
  delete f ;
  delete h ;
  delete eos ;
  return 0;
}
