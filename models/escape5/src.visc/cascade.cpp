#include <TRandom3.h>
#include <iostream>

#include "cascade.h"
#include "gen.h"
#include "UKUtility.h"
#include "DatabasePDG2.h"
#include "particle.h"
#include "params.h"

int ievcasc ; // event #

using namespace std ;

void cxxninit_(int *np)
{
 *np = gen::npart[ievcasc] ;
}


void cxxinit_(int* index, int* id1, float* x1, float* y1, float* z1, float* t1, float* px1, float* py1, float* pz1, float* E1, float* mass1)
{
// std::cout << "Particle initialized!\n" ;
 Particle* p = gen::pList[ievcasc][*index-1] ;
 *id1 = p->def->GetPDG() ;
 *x1  = p->x ;
 *y1  = p->y ;
 *z1  = p->z ;
 *t1  = p->t ;
 *px1 = p->px ;
 *py1 = p->py ;
 *pz1 = p->pz ;
 *E1  = p->e ;
 *mass1 = sqrt((*E1)*(*E1)-(*px1)*(*px1)-(*py1)*(*py1)-(*pz1)*(*pz1)) ;
 return ;
}


void cxxnfinal_(int *np)
{
 for(int ip=gen::npart[ievcasc]; ip<*np; ip++)
 gen::pList[ievcasc][ip] = new Particle(0.,0.,0.,0., 0.,0.,0.,0., (ParticlePDG2*)0x0,0) ;
 gen::npart[ievcasc] = *np ;
//! std::cout << "cxx : created " << *np << " final particles\n" ;
}


void cxxfinal_(int* index, int* id1, float* x1, float* y1, float* z1, float* t1, float* px1, float* py1, float* pz1, float* E1, float* mass1)
{
 Particle* p = gen::pList[ievcasc][*index-1] ;
 p->def = gen::database->GetPDGParticle(*id1) ;
 p->mid = 0 ;
 p->x = *x1 ;
 p->y = *y1 ;
 p->z = *z1 ;
 p->t = *t1 ;
 p->px = *px1 ;
 p->py = *py1 ;
 p->pz = *pz1 ;
 p->e = *E1 ;
//cout << "FPART: " << setw(14) << *px1 << setw(14) << *py1 << setw(14) << *pz1 << endl ;
}

namespace gen{

using params::rescatter ;
extern DatabasePDG2 *database ;


// here we involve UrQMD
void urqmd(int iev)
{
 ievcasc = iev ;
 if(rescatter) // involve UrQMD
 urqmdmain_() ;
 // no rescattering - just pass particle arrays ass they are
//===== decay of unstable resonances ========
for(int iiter=0; iiter<3; iiter++){

 for(int ipart=0; ipart<npart[iev]; ipart++){
 Particle* p = pList[iev][ipart] ;
 if(p->def==0) { /*cout << "unknown particle: " << p->def << endl ;*/ continue ; }
 if(p->def->GetWidth()>0. && !isStable(p->def->GetPDG())){
  p->x = p->x  + p->px/p->e*(400. - p->t) ;
  p->y = p->y  + p->py/p->e*(400. - p->t) ;
  p->z = p->z  + p->pz/p->e*(400. - p->t) ;
  p->t = 400. ;
#ifdef DEBUG2
  cout << "------ unstable particle decay " << Id[i] << endl ;
  cout << setw(14) << "px" << setw(14) << "py" << setw(14) << "pz" << setw(14) << "E" << setw(14) << "m" << endl ;
  cout << setw(14) << mom[0] << setw(14) << mom[1] << setw(14) << mom[2] << setw(14) << mom[3] << setw(14) << mom[4] << endl ;
#endif
  int nprod ;
  Particle** daughters ;
  decay(p, nprod, daughters) ;
#ifdef DEBUG2
  cout << "decay into : " ; for(int iprod=0; iprod<nprod; iprod++) cout << "  " << ppid[iprod] ;
  cout << endl ;
  for(int iprod=0; iprod<nprod; iprod++)
   cout << setw(14) << pmom[iprod][0] << setw(14) << pmom[iprod][1] << setw(14) << pmom[iprod][2] << setw(14) << pmom[iprod][3] << setw(14) << pmom[iprod][4] << endl ;
#endif
 //------------------ adding daughters to list (daughter #0 replaces original particle)
  pList[iev][ipart] = daughters[0] ;
  for(int iprod=1; iprod<nprod; iprod++){
    pList[iev][npart[iev]+iprod-1] = daughters[iprod] ;
  }
  npart[iev] += nprod-1 ;
  delete [] daughters ;
//--------------------------------------------
  } // decay procedure
 }
 
 } // decay iteration
}

} // end namespace gen

//##########  resonance decay kinematics from FASTMC ###################

double BreitWigner(double mean, double gamma, TRandom* random)
{
//  Return a number distributed following a BreitWigner function with mean and gamma

   double rval, displ;
   rval = 2*random->Rndm() - 1;
   displ = 0.5*gamma*TMath::Tan(rval*TMath::PiOver2());
   return (mean+displ);
}


void decay(Particle *in, int& nprod, Particle** &out){
  DatabasePDG2 *database = gen::database ;
  ParticlePDG2* pDef = in->def ;
  TRandom3 *random3 = gen::rnd ;
  Int_t nDecayChannel = pDef->GetNDecayChannels();
  
  Double_t fullBranching = pDef->GetFullBranching(); // Only 3 or less body decays
  
  nprod = 0 ;
  if(fullBranching < 0.0000001)
    return;
  
  Double_t initialMass0 = pDef->GetMass();
  //Breit Wigner distributon for masses
M1:  Double_t initialMass = BreitWigner(initialMass0,pDef->GetWidth(), random3); 
	//if(initialMass>2.*initialMass0) initialMass = 2.*initialMass0 ;
	if(initialMass>2.*initialMass0) goto M1 ;

  //ml: BW for Xi*(1530) 
  //if(abs(pDef->GetPDG())==3314||abs(pDef->GetPDG())==3324){ 
 // std::cout<<pDef->GetPDG()<<" mi   "<<initialMass0<<" w "<<pDef->GetWidth()<<" mbw "<<initialMass<<std::endl;
  initialMass0=initialMass;
  //}

  Int_t nAllowedChannels = database->GetNAllowedChannels(pDef, initialMass0);
  if(pDef->GetWidth()>0 && nDecayChannel==0) {
//    cout << "Error in HadronDecayer::Decay() : Particle " << pDef->GetPDG() << " has finite width (" << pDef->GetWidth()
//	 << ") but NO decay channels specified in the database. Check it out!!" << endl;
   initialMass0 = pDef->GetMass();
   goto M1;
//test    exit(0);
  }
  if(nAllowedChannels==0) {
//    cout << "Error in HadronDecayer::Decay() : Particle " << pDef->GetPDG() << " has " << nDecayChannel
//	 << " decay channels but NONE of them is allowed!!" << endl;
   initialMass0 = pDef->GetMass();
   goto M1;
//    exit(0);
  }


  // we need to choose an allowed decay channel
  Double_t randValue = random3->Rndm() * fullBranching;
  Int_t chosenChannel = 1000;
  Bool_t found = kFALSE;
  Int_t iterations = 0;
  while(!found) {
    if(iterations > 100) {
      cout << "Warning in HadronDecayer::Decay() : More than 100 iterations to choose a decay channel. Check it out !!" << endl;
    }
    for(Int_t nChannel = 0; nChannel < nDecayChannel; ++nChannel) {
      randValue -= pDef->GetDecayChannel(nChannel)->GetBranching();
      if(randValue <= 0. && database->IsChannelAllowed(pDef->GetDecayChannel(nChannel), initialMass0)) {
	chosenChannel = nChannel;
	found = kTRUE;
	break;
      }
    }
    iterations++;
  }

  DecayChannel *dc = pDef->GetDecayChannel(chosenChannel);
  
  Int_t nSec = dc->GetNDaughters();

//  Particle parent1(database->GetPDGParticle(parent.Encoding()));
  Double_t E1=TMath::Sqrt(in->px*in->px+in->py*in->py+in->pz*in->pz+initialMass0*initialMass0);
  TLorentzVector parentMom(in->px,in->py,in->pz,E1); 
  
//take into account BW
 // TVector3 velocity(parent.Mom().BoostVector());
  TVector3 velocity(parentMom.BoostVector());
  

  //decay is allowed:
  if(nSec == 1) {
    nprod = 1 ;
    out = new Particle* [1] ;
    out[0] = new Particle(in->x, in->y, in->z, in->t, in->px, in->py, in->pz, E1, pDef, pDef->GetPDG()) ;
    cout<<"Warning: 1-particle decay!, pid= "<<pDef->GetPDG()<<endl ;
    return;
    //store information about mother
  } 
  else if(nSec == 2) {
    //two body decay
    nprod = 2 ;
    ParticlePDG2* daughter1 = database->GetPDGParticle(dc->GetDaughterPDG(0)) ;
    ParticlePDG2* daughter2 = database->GetPDGParticle(dc->GetDaughterPDG(1)) ;
    
    double p1mass = daughter1->GetMass() ;
    double p2mass = daughter2->GetMass() ;
    TLorentzVector p1mom, p2mom ;
    MomAntiMom(p1mom, p1mass, p2mom, p2mass, initialMass0);
    
 //   if(abs(parent.Encoding())==3314||abs(parent.Encoding())==3324){
//    std::cout<<"CMS--- 2 decay --- "<< parent.Encoding()<<" "<<initialMass0<<std::endl;
//    std::cout<<"CMS deltas"<<p1.Mom().X()+p2.Mom().X()<<" "<<p1.Mom().Y()+p2.Mom().Y()<<" "<<
//    p1.Mom().Z()+p2.Mom().Z()<<std::endl;}

    // boost to Lab system
    p1mom.Boost(velocity);
    p2mom.Boost(velocity);
    
    out = new Particle* [2] ;
    out[0] = new Particle(in->x, in->y, in->z, in->t, p1mom.Px(), p1mom.Py(), p1mom.Pz(), p1mom.E(), daughter1, pDef->GetPDG()) ;
    out[1] = new Particle(in->x, in->y, in->z, in->t, p2mom.Px(), p2mom.Py(), p2mom.Pz(), p2mom.E(), daughter2, pDef->GetPDG()) ;

    double delta = TMath::Sqrt(
    (parentMom.X()-p1mom.X()-p2mom.X())*(parentMom.X()-p1mom.X()-p2mom.X())+
    (parentMom.Y()-p1mom.Y()-p2mom.Y())*(parentMom.Y()-p1mom.Y()-p2mom.Y())+
    (parentMom.Z()-p1mom.Z()-p2mom.Z())*(parentMom.Z()-p1mom.Z()-p2mom.Z())+
    (parentMom.E()-p1mom.E()-p2mom.E())*(parentMom.E()-p1mom.E()-p2mom.E()));


    if(delta>0.001){
//    if(abs(parent.Encoding())==3314||abs(parent.Encoding())==3324){
    std::cout<<"after boost--- 2 decay --- "<< delta<<" "<<pDef->GetPDG()<<" "<<initialMass0<<std::endl;
//    std::cout<<" after boost deltas"<<parent.Mom().X()-p1.Mom().X()-p2.Mom().X()<<" "<<parent.Mom().Y()-p1.Mom().Y()-p2.Mom().Y()<<" "<<
//    parent.Mom().Z()-p1.Mom().Z()-p2.Mom().Z()<<std::endl;}

   
/*
    std::cout<<"  "<< parent.Encoding()<<" "<<parent.Mom().X()<<" "<<parent.Mom().Y()<<" "
    <<parent.Mom().Z()<<std::endl;

    std::cout<<"  "<< p1.Encoding()<<" "<<p1.Mom().X()<<" "<<p1.Mom().Y()<<" "
    <<p1.Mom().Z()<<std::endl;

    std::cout<<"  "<< p2.Encoding()<<" "<<p2.Mom().X()<<" "<<p2.Mom().Y()<<" "
    <<p2.Mom().Z()<<std::endl;
*/

   initialMass0 = pDef->GetMass();
   goto M1;

     }



    return;       
  } 
  
  else if(nSec == 3) {
    // three body decay
    nprod = 3 ;
    // calculate daughter Pabs momentum
    Double_t pAbs1 = 0., pAbs2 = 0., pAbs3 = 0., sumPabs = 0., maxPabs = 0.;
    ParticlePDG2* daughter1 = database->GetPDGParticle(dc->GetDaughterPDG(0)) ;
    ParticlePDG2* daughter2 = database->GetPDGParticle(dc->GetDaughterPDG(1)) ;
    ParticlePDG2* daughter3 = database->GetPDGParticle(dc->GetDaughterPDG(2)) ;
    Double_t mass1 = daughter1->GetMass(), 
    mass2 = daughter2->GetMass(), 
    mass3 = daughter3->GetMass();
    TLorentzVector mom1, mom2, mom3; 
    Double_t deltaMass = initialMass0 - mass1 - mass2 - mass3;
    do {
      Double_t rd1 = random3->Rndm();
      Double_t rd2 = random3->Rndm();
      
      if (rd2 > rd1)
	std::swap(rd1, rd2);
      // 1
      Double_t e = rd2*deltaMass;
      pAbs1 = TMath::Sqrt(e*e + 2*e*mass1);
      sumPabs = pAbs1;
      maxPabs = sumPabs;
      // 2
      e = (1-rd1)*deltaMass;
      pAbs2 = TMath::Sqrt(e*e + 2*e*mass2);
      
      if(pAbs2 > maxPabs)
	maxPabs = pAbs2;
      
      sumPabs += pAbs2;
      // 3
      e = (rd1-rd2)*deltaMass;
      pAbs3 = TMath::Sqrt(e*e + 2*e*mass3);
      
      if (pAbs3 > maxPabs)
	maxPabs =  pAbs3;
      sumPabs  +=  pAbs3;
    } while(maxPabs > sumPabs - maxPabs);
    
    // isotropic sample first particle 3-momentum
    Double_t cosTheta = 2 * (random3->Rndm()) - 1;
    Double_t sinTheta = TMath::Sqrt(1 - cosTheta*cosTheta);
    Double_t phi      = TMath::TwoPi()*(random3->Rndm());
    Double_t sinPhi   = TMath::Sin(phi);
    Double_t cosPhi   = TMath::Cos(phi);
    
    mom1.SetPxPyPzE(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta, 0);
    mom1 *= pAbs1;
    // sample rest particle 3-momentum
    Double_t cosThetaN = (pAbs2 * pAbs2 - pAbs3 * pAbs3 - pAbs1 * pAbs1)
      / (2 * pAbs1 * pAbs3);
    Double_t sinThetaN = TMath::Sqrt(1 - cosThetaN * cosThetaN);
    Double_t phiN      = TMath::TwoPi() * (random3->Rndm());
    Double_t sinPhiN   = TMath::Sin(phiN);
    Double_t cosPhiN   = TMath::Cos(phiN);
    
    mom3.SetPxPyPzE(sinThetaN * cosPhiN * cosTheta * cosPhi - sinThetaN * sinPhiN * sinPhi + cosThetaN * sinTheta * cosPhi,
		    sinThetaN * cosPhiN * cosTheta * sinPhi + sinThetaN * sinPhiN * cosPhi + cosThetaN * sinTheta * sinPhi,
		    -sinThetaN * cosPhiN * sinTheta + cosThetaN * cosTheta,
		    0.);
    
    mom3 *= pAbs3*mom3.P();
    mom2 = mom1;
    mom2 += mom3;
    mom2 *= -1.;
    // calculate energy
    mom1.SetE(TMath::Sqrt(mom1.P() * mom1.P() + mass1 * mass1));
    mom2.SetE(TMath::Sqrt(mom2.P() * mom2.P() + mass2 * mass2));
    mom3.SetE(TMath::Sqrt(mom3.P() * mom3.P() + mass3 * mass3));
    // boost to Lab system
    mom1.Boost(velocity);
    mom2.Boost(velocity);
    mom3.Boost(velocity);
    
    
//test energy conservation

    double delta =TMath::Sqrt(
    (parentMom.X()-mom1.X()-mom2.X()-mom3.X())*
    (parentMom.X()-mom1.X()-mom2.X()-mom3.X())
   +
    (parentMom.Y()-mom1.Y()-mom2.Y()-mom3.Y())
    *(parentMom.Y()-mom1.Y()-mom2.Y()-mom3.Y())
   +
    (parentMom.Z()-mom1.Z()-mom2.Z()-mom3.Z())*
    (parentMom.Z()-mom1.Z()-mom2.Z()-mom3.Z())
    +
    (parentMom.E()-mom1.E()-mom2.E()-mom3.E())*
    (parentMom.E()-mom1.E()-mom2.E()-mom3.E())); 

    if(delta>1.E-2){
    std::cout<<"bad decay--- 3 decay --- "<< delta<<" "<<pDef->GetPDG()<<" "<<initialMass0<<std::endl;

/*
    std::cout<<"  "<< parent.Encoding()<<" "<<parent.Mom().X()<<" "<<parent.Mom().Y()<<" "
    <<parent.Mom().Z()<<std::endl;

    std::cout<<"  "<< p1.Encoding()<<" "<<p1.Mom().X()<<" "<<p1.Mom().Y()<<" "
    <<p1.Mom().Z()<<std::endl;

    std::cout<<"  "<< p2.Encoding()<<" "<<p2.Mom().X()<<" "<<p2.Mom().Y()<<" "
    <<p2.Mom().Z()<<std::endl;

    std::cout<<"  "<< p3.Encoding()<<" "<<p3.Mom().X()<<" "<<p3.Mom().Y()<<" "
    <<p3.Mom().Z()<<std::endl;
*/  
   initialMass0 = pDef->GetMass();
   goto M1;

    return;

    }

    out = new Particle* [3] ;
    out[0] = new Particle(in->x, in->y, in->z, in->t, mom1.Px(), mom1.Py(), mom1.Pz(), mom1.E(), daughter1, pDef->GetPDG()) ;
    out[1] = new Particle(in->x, in->y, in->z, in->t, mom2.Px(), mom2.Py(), mom2.Pz(), mom2.E(), daughter2, pDef->GetPDG()) ;
    out[2] = new Particle(in->x, in->y, in->z, in->t, mom3.Px(), mom3.Py(), mom3.Pz(), mom3.E(), daughter3, pDef->GetPDG()) ;

    return;
  }
  	 
}
