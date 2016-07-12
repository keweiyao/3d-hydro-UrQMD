#include <cstdlib>
#include <iostream>
#include <math.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>

TRandom3 *rndm ;
int ninipart =100 ;
double pmax = 0.5 ;

//--------particle arrays -------
int nfpart ;
int *id ;
float *x, *y, *z, *t, *px, *py, *pz, *E, *mass ;
TTree *tree ;
//-------------------------------

extern"C" {
 void urqmdmain_();

 void cxxinit_(int* index, int* id, float* x, float* y, float* z, float* t, float* px, float* py, float* pz, float* E, float* mass);
 void cxxninit_(int *np) ;
 void storepart_(int *np) ;
 void cxxnfinal_(int *np) ;
 void cxxfinal_(int* index, int* id, float* x, float* y, float* z, float* t, float* px, float* py, float* pz, float* E, float* mass);
}


void cxxinit_(int* index, int* id1, float* x1, float* y1, float* z1, float* t1, float* px1, float* py1, float* pz1, float* E1, float* mass1)
{
 std::cout << "Particle initialized!\n" ;
 *id1 = 2212 ;
 *x1 = 10.*rndm->Rndm() ;
 *y1 = 10.*rndm->Rndm() ;
 *z1 = 10.*rndm->Rndm() ;
 *t1 = 1. ;
 *px1 = -pmax + 2.0*pmax*rndm->Rndm() ;
 *py1 = -pmax + 2.0*pmax*rndm->Rndm() ;
 *pz1 = -pmax + 2.0*pmax*rndm->Rndm() ;
 *mass1 = 0.938 ;
 *E1 = sqrt((*mass1)*(*mass1)+(*px1)*(*px1)+(*py1)*(*py1)+(*pz1)*(*pz1)) ;
 return ;
}


void cxxnfinal_(int *np)
{
 nfpart = *np ;
 std::cout << "cxx : created " << *np << " final particles\n" ;
 id = new int [*np] ;
 x = new float [*np] ;
 y = new float [*np] ;
 z = new float [*np] ;
 t = new float [*np] ;
 px = new float [*np] ;
 py = new float [*np] ;
 pz = new float [*np] ;
 E = new float [*np] ;
 mass = new float [*np] ;
}

void cxxfinal_(int* index, int* id1, float* x1, float* y1, float* z1, float* t1, float* px1, float* py1, float* pz1, float* E1, float* mass1)
{
 id[*index] = *id1 ;
 x[*index] = *x1 ;
 y[*index] = *y1 ;
 z[*index] = *z1 ;
 t[*index] = *t1 ;
 px[*index] = *px1 ;
 py[*index] = *py1 ;
 pz[*index] = *pz1 ;
 E[*index] = *mass1 ;
 mass[*index] = *E1 ;
}


void cxxninit_(int *np)
{
 *np = ninipart ;
}


//----------------- plots
void plots()
{
 TGraph *g1 = new TGraph(nfpart, x, t) ;
 g1->Draw("AP*") ;
}
//-----------------------


int main(int argc, char** argv)
{

std::cout << "Hello, world\n" ;
 rndm = new TRandom3() ;
//--------- writing to file
//storepart_(&nparticles) ;
//exit(1) ;
  TFile outputFile("RunOutput.root", "RECREATE"); 
  outputFile.cd();
  tree = new TTree("ti","Initial");
  tree->Branch("npart",&nfpart,"npart/I");
//  tree->Branch("nev",&nev,"nev/I");
  tree->Branch("x",&x[0],"x[npart]/F");

//----------- UrQMD execution
 urqmdmain_() ;
 tree->Fill() ;

outputFile.Write() ;

//------- graphical part
 TApplication theApp("App", &argc, argv);
 plots() ;
theApp.Run() ;
 return 0 ;
}
