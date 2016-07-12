#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include "params.h"

using namespace std ;

namespace params{

char sSurface [255], sSpectraDir [255], sMultDir [255] ;
bool bEventGeneration ;
bool weakContribution ;
bool rescatter ;
bool shear ;
//double Temp, mu_b, mu_q, mu_s ;
int NEVENTS ;
double NBINS, QMAX ;
double dx, dy, deta ;
double ecrit ;

// ############ reading and processing the parameters

void readParams(char* filename)
{
	char parName [255], parValue [255] ;
	ifstream fin(filename) ;
	if(!fin.is_open()) { cout << "cannot open parameters file " << filename << endl ; exit(1) ; }
	while(fin.good()){
	 string line ;
	 getline(fin, line) ;
	 istringstream sline (line) ;
	 sline >> parName >> parValue ;
	 if     (strcmp(parName,"surface")==0) strcpy(sSurface, parValue) ;
	 else if(strcmp(parName,"spectra_dir")==0) strcpy(sSpectraDir, parValue) ;
	 //	 else if(strcmp(parName,"fmax_dir")==0) strcpy(sMultDir, parValue) ;
//	 else if(strcmp(parName,"T_ch")==0) Temp = atof(parValue) ;
//	 else if(strcmp(parName,"mu_b")==0) mu_b = atof(parValue) ;
//	 else if(strcmp(parName,"mu_q")==0) mu_q = atof(parValue) ;
//	 else if(strcmp(parName,"mu_s")==0) mu_s = atof(parValue) ;
	 else if(strcmp(parName,"Nbins")==0) NBINS = atoi(parValue) ;
	 else if(strcmp(parName,"q_max")==0) QMAX = atof(parValue) ;
	 else if(strcmp(parName,"number_of_events")==0) NEVENTS = atoi(parValue) ;
	 else if(strcmp(parName,"rescatter")==0) rescatter = atoi(parValue) ;
	 else if(strcmp(parName,"weakContribution")==0) weakContribution = atoi(parValue) ;
   else if(strcmp(parName,"shear")==0) shear = atoi(parValue) ;
   else if(strcmp(parName,"ecrit")==0) ecrit = atof(parValue) ;
	 else if(parName[0]=='!') cout << "CCC " << sline.str() << endl ;
	 else cout << "UUU " << sline.str() << endl ;
	}
 deta=0.05 ; dx=dy=0.0 ; // TODO!
}

void printParameters()
{
  cout << "====== parameters ======\n" ;
  cout << "surface = " << sSurface << endl ;
  cout << "spectraDir = " << sSpectraDir << endl ;
  cout << "multiplicityDir = " << sMultDir << endl ;
//  cout << "T_ch = " << Temp << endl ;
//  cout << "mu_b = " << mu_b << endl ;
//  cout << "mu_q = " << mu_q << endl ;
//  cout << "mu_s = " << mu_s << endl ;
  cout << "eventGeneration = " << bEventGeneration << endl ;
  cout << "numberOfEvents = " << NEVENTS << endl ;
  cout << "isRescatter = " << rescatter << endl ;
  cout << "weakContribution = " << weakContribution << endl ;
  cout << "shear_visc_on = " << shear << endl ;
  cout << "e_critical = " << ecrit << endl ;
  cout << "Nbins = " << NBINS << endl ;
  cout << "q_max = " << QMAX << endl ;
  cout << "======= end parameters =======\n" ;
}

}
