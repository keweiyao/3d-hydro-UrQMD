#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <sstream>
#include "fld.h"
#include "hdo.h"
#include "IC_reader.h"
#include "eos.h"
#include "eos_hotqcd.h"
#include "trancoeff.h"

using namespace std ;

// program parameters, to be read from file
int nx, ny, nz;
double xmax, ymax, etamax, 
       tau0, tauMax, dtau ;

int ic_nxy, ic_neta;
double ic_dxy, ic_deta;

char outputDir[255];
char icInputFile [255] ;
double etaS_min, etaS_slope_QGP, etaS_HRG, zetaS, eCrit;

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
	 if     (strcmp(parName,"outputdir")==0) strcpy(outputDir, parValue) ;
	 else if(strcmp(parName,"icinputfile")==0) strcpy(icInputFile, parValue) ;
	 else if(strcmp(parName,"nx")==0) nx = atoi(parValue) ;
	 else if(strcmp(parName,"ny")==0) ny = atoi(parValue) ;
	 else if(strcmp(parName,"nz")==0) nz = atoi(parValue) ;
	 else if(strcmp(parName,"xmax")==0) xmax = atof(parValue) ;
	 else if(strcmp(parName,"ymax")==0) ymax = atof(parValue) ;
	 else if(strcmp(parName,"etamax")==0) etamax = atof(parValue) ;
	 else if(strcmp(parName,"tau0")==0) tau0 = atof(parValue) ;
	 else if(strcmp(parName,"taumax")==0) tauMax = atof(parValue) ;
	 else if(strcmp(parName,"dtau")==0) dtau = atof(parValue) ;

	 else if(strcmp(parName,"ic_nxy")==0) ic_nxy = atoi(parValue) ;
     else if(strcmp(parName,"ic_neta")==0) ic_neta = atoi(parValue) ;
     else if(strcmp(parName,"ic_dxy")==0) ic_dxy = atof(parValue) ;
	 else if(strcmp(parName,"ic_deta")==0) ic_deta = atof(parValue) ;

     else if(strcmp(parName,"e_crit")==0) eCrit = atof(parValue);
	 else if(strcmp(parName,"etas_min")==0) etaS_min = atof(parValue) ;
	 else if(strcmp(parName,"etas_slope_qgp")==0) etaS_slope_QGP = atof(parValue) ;
	 else if(strcmp(parName,"etas_hrg")==0) etaS_HRG = atof(parValue) ;
	 else if(strcmp(parName,"zetas")==0) zetaS = atof(parValue) ;

	 else if(parName[0]=='!') cout << "CCC " << sline.str() << endl ;
	 else cout << "UUU " << sline.str() << endl ;
	}
}


void printParameters()
{
  cout << "====== parameters ======\n" ;
  cout << "outputDir = " << outputDir << endl ;
  cout << "icInputFile = " << icInputFile << endl ;
  cout << "hydro grid = " <<  nx << " x " << ny << " x " << nz << endl ;
  cout << "hydro area = " << "[" << -xmax << "," << xmax << "]" << " x " 
                          << "[" << -ymax << "," << ymax << "]" << " x "
                          << "[" << -etamax << "," << etamax << "]" << endl;
  cout << "hydro time = " << tau0 << " : " << tauMax << " : " << dtau << endl ;
  cout << "ic grid = "    << ic_nxy << " x " << ic_nxy << " x " << ic_neta << endl ;
  cout << "ic cell = "    << ic_dxy << " x " << ic_dxy << " x " << ic_deta << endl ;
  cout << "E_switch = "<< eCrit << endl;
  cout << "eta/s = " << etaS_min << "+ (T-Tc)*(T>Tc)*" << etaS_slope_QGP << " + (T<Tc)*" << etaS_HRG  << endl ;
  cout << "zeta/s = " << zetaS << endl ;
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
  trcoeff = new TransportCoeff(etaS_min, etaS_slope_QGP, etaS_HRG, zetaS, eos) ;
  
  // fluid
  f = new Fluid(eos, trcoeff, nx, ny, nz, -xmax, xmax, -ymax, ymax, -etamax, etamax, dtau, eCrit) ;
  
  // initilal conditions (read from file)
  IC_reader *ic = new IC_reader(icInputFile, eos, ic_nxy, ic_neta, ic_dxy, ic_deta);
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
