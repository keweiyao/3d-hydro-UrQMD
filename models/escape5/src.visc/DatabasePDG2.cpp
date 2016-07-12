
#ifndef DATABASE_PDG2
#include "DatabasePDG2.h"
#endif
#include "params.h"

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>
using std::cout;
using std::endl;
using std::strcmp;
using std::strcpy;
using std::ifstream;
using std::ios;
using std::setw;
using std::ofstream;
using std::strcpy;
using std::strcat;

using params::weakContribution ;
int rnegative ; // negative factorial argument flag - for Clebsch-Gordan coeff. calculation

bool isStable(int pdg) ;

int fact(double n)
{
// cout << "::fact " << n << endl ;
 double result = 1;
 if(n<0) rnegative = 1 ;
 while(n>1){
  result *= n ;
  n-- ;
 }
 return result ;
}


double Clebsch_Gordan(double j1, double m1, double j2, double m2, double J, double M)
{
// return 1. ;
 double result = ((M-m1-m2)==0. ? 1. : 0.)*
  sqrt((2.*J+1.)*fact(J+j1-j2)*fact(J-j1+j2)*fact(j1+j2-J)/fact(j1+j2+J+1.))*
  sqrt(fact(J+M)*fact(J-M)*fact(j1-m1)*fact(j1+m1)*fact(j2-m2)*fact(j2+m2)) ;
 double sum = 0. ;
 for(int k=0; k<10; k++){
  rnegative = 0 ;
  double value = pow(-1., k)/(fact(k)*fact(j1+j2-J-k)*fact(j1-m1-k)*fact(j2+m2-k)*fact(J-j2+m1+k)*fact(J-j1-m2+k)) ;
  if(!rnegative) sum += value ;
 }
 result *= sum ;
 return result ;
}


DatabasePDG2::DatabasePDG2(char *fileparticles, char *filedecay) {
  fNParticles = 0;
  strcpy(fParticleFilename, fileparticles);
  strcpy(fDecayFilename, filedecay);
  for(Int_t i=0; i<kMaxParticles2; i++) {
    fParticles[i] = new ParticlePDG2();
    fStatus[i] = kFALSE;
  }
  fUseCharmParticles = kTRUE;
  fMinimumWidth = 0.;
  fMaximumWidth = 10.;
  fMinimumMass = 0.01;
  fMaximumMass = 10.;
};

DatabasePDG2::~DatabasePDG2() {
  for(Int_t i=0; i<kMaxParticles2; i++)
    delete fParticles[i];
};

void DatabasePDG2::SetParticleFilename(Char_t *filename) {
  strcpy(fParticleFilename, filename);
}

void DatabasePDG2::SetDecayFilename(Char_t *filename) {
  strcpy(fDecayFilename, filename);
}

Bool_t DatabasePDG2::LoadData() {
  return (LoadParticles() && LoadDecays());
//  return LoadParticles() ;
}

Bool_t DatabasePDG2::LoadParticles() {
  ifstream particleFile;
  particleFile.open(fParticleFilename);
  if(!particleFile) {
    cout << "ERROR in DatabasePDG::LoadParticles() : The ASCII file containing the PDG particle list (\""
         << fParticleFilename << "\") was not found !! Aborting..." << endl;
    return kFALSE;
  }
  
  Char_t name [30], recipt [30], gparity[10], cparity[10], pparity[10], antiparticle[10], charge[10], status[10] ;
  Double_t mass, width, spin, isospin;
  Int_t pdg;
  Int_t goodStatusParticles = 0;
  Char_t s_plusplus [] = "++" ;
  Char_t s_plus [] = "+" ;
  Char_t s_minus [] = "-" ;
  Char_t s_minusminus [] = "--" ;
  Char_t s_plusplusbar [] = "++bar" ;
  Char_t s_plusbar [] = "+bar" ;
  Char_t s_minusbar [] = "-bar" ;
  Char_t s_minusminusbar [] = "--bar" ;
  Char_t s_bar [] = "bar" ;

//  cout << "Info in DatabasePDG::LoadParticles() : Start loading particles with the following criteria:" << endl
//       << "       Use particles containing charm quarks (1-yes;0-no) : " << fUseCharmParticles << endl
//       << "       Mass range                                         : (" << fMinimumMass << "; " << fMaximumMass << ")" << endl
//       << "       Width range                                        : (" << fMinimumWidth << "; " << fMaximumWidth << ")" << endl;
  
  particleFile.exceptions(ios::failbit);
  while(!particleFile.eof()) {
    try {
      particleFile >> mass >> width >> isospin >> gparity >> spin >> pparity >> cparity >> antiparticle >> pdg >> charge >> status >> name >> recipt ;
      mass = mass/1000. ;
      width = width/1000. ;
      if(width<0.) width = 0. ;
//      cout << setw(15) << name << setw(12) << mass << setw(12) << width << setw(5) << spin << setw(5) << isospin << setw(12) << pdg ;
    }
    catch (ios::failure const &problem) {
      //cout << "reading particle file:\n" ;
      cout << problem.what() << endl;
      break;
    }
        
    fParticles[fNParticles]->SetName(name);
    fParticles[fNParticles]->SetPDG(pdg);
    fParticles[fNParticles]->SetMass(mass);
    fParticles[fNParticles]->SetWidth(width);
    fParticles[fNParticles]->SetIsospin(isospin) ;
    fParticles[fNParticles]->SetSpin(spin) ;
    goodStatusParticles++;
    fStatus[fNParticles] = kTRUE;

	fParticles[fNParticles]->SetPDGstatus(status[0]) ;
	
	switch(gparity[0]){
	case '+' : fParticles[fNParticles]->SetGparity(1) ; break ;
	case '-' : fParticles[fNParticles]->SetGparity(-1) ; break ;
	default : fParticles[fNParticles]->SetGparity(0) ;
	}
	
	switch(pparity[0]){
	case '+' : fParticles[fNParticles]->SetPparity(1) ; break ;
	case '-' : fParticles[fNParticles]->SetPparity(-1) ; break ;
	default : fParticles[fNParticles]->SetPparity(0) ;
	}
	
	switch(cparity[0]){
	case '+' : fParticles[fNParticles]->SetCparity(1) ; break ;
	case '-' : fParticles[fNParticles]->SetCparity(-1) ; break ;
	default : fParticles[fNParticles]->SetCparity(0) ;
	}
	
	switch(charge[0]){
	case '+' : if(charge[1]=='+') fParticles[fNParticles]->SetElectricCharge(2.) ;
		else  fParticles[fNParticles]->SetElectricCharge(1.) ;
		break ;
	case '-' : fParticles[fNParticles]->SetElectricCharge(-1.) ; break ;
	case '0' : fParticles[fNParticles]->SetElectricCharge(0.) ; break ;
	default : cout << "Unknown charge: " << charge << endl ; exit(1) ;
	}
	
	// calculating flavors
	int baryonnum=0, echarge=0, scharge=0, ccharge=0, bcharge=0, tcharge=0, sqnumber=0, saqnumber=0 ;
	int irec=0 ;
	while(recipt[irec]!='\0'){
		switch(recipt[irec]){
		case 'u' : baryonnum += 1; echarge += 2; break ;
		case 'U' : baryonnum += -1; echarge += -2; break ;
		case 'd' : baryonnum += 1; echarge += -1; break ;
		case 'D' : baryonnum += -1; echarge += 1; break ;
		case 's' : baryonnum += 1; echarge += -1; scharge += -1; sqnumber += 1 ; break ;
		case 'S' : baryonnum += -1; echarge += 1; scharge += 1; saqnumber += 1 ; break ;
		case 'c' : baryonnum += 1; echarge += 2; ccharge += 1 ; break ;
		case 'C' : baryonnum += -1; echarge += -2; ccharge += -1 ; break ;
		case 't' : baryonnum += 1; echarge += 2; tcharge += 1 ; break ;
		case 'T' : baryonnum += -1; echarge += -2; tcharge += -1 ; break ;
		case 'b' : baryonnum += 1; echarge += -1; bcharge += -1 ; break ;
		case 'B' : baryonnum += -1; echarge += 1; bcharge += 1 ; break ;
		}
		irec++ ;

	}
	if(baryonnum%3 != 0) { cout << "Fractonal baryon charge? exiting...\n"; exit(1) ; }
	if(echarge%3 != 0) { cout << "Fractonal charge? exiting...\n"; exit(1) ; }

	baryonnum /= 3. ;
	echarge /= 3. ;
	fParticles[fNParticles]->SetBaryonNumber(baryonnum) ;
	fParticles[fNParticles]->SetStrangeness(scharge) ;
	fParticles[fNParticles]->SetStrangeQNumber(sqnumber) ;
	fParticles[fNParticles]->SetStrangeAQNumber(saqnumber) ;
	fParticles[fNParticles]->SetCharm(ccharge) ;
	fParticles[fNParticles]->SetBeauty(bcharge) ;
	fParticles[fNParticles]->SetTruth(tcharge) ;
	fParticles[fNParticles]->SetIsospinZ(fParticles[fNParticles]->GetElectricCharge()-0.5*(scharge+baryonnum+ccharge+tcharge+bcharge)) ;
	if(fabs(fParticles[fNParticles]->GetIsospinZ())>fParticles[fNParticles]->GetIsospin())
		{ cout << "I_3 > I. Error, exiting... " << fParticles[fNParticles]->GetPDG() << "  " << 
		fParticles[fNParticles]->GetIsospinZ() << "  " << fParticles[fNParticles]->GetIsospin() << "\n"; exit(1) ; }


//     check that the particle mass is inside accepted limits
//     check that the particle width is inside accepted limits
//     check that it has a good PDG status
    if(!(fMinimumMass<=mass && mass<=fMaximumMass) || !(fMinimumWidth<=width && width<=fMaximumWidth) || abs(pdg)<100
	|| status[0]=='S' || status[0]=='F' || bcharge!=0 || ccharge!=0 || tcharge!=0) {
      fStatus[fNParticles] = kFALSE;
      goodStatusParticles--;
    }

    fNParticles++;

    //creating antiparticle
    if(antiparticle[0]=='B' || antiparticle[0]=='F'){
	fParticles[fNParticles]->SetName(name) ;
	fParticles[fNParticles]->SetPDG(-pdg);
	fParticles[fNParticles]->SetMass(mass);
	fParticles[fNParticles]->SetWidth(width);
	fParticles[fNParticles]->SetIsospin(isospin) ;
	fParticles[fNParticles]->SetSpin(spin) ;
	fStatus[fNParticles] = fStatus[fNParticles-1];
	if(fStatus[fNParticles]==kTRUE) goodStatusParticles++;

	if(fParticles[fNParticles-1]->GetBaryonNumber()!=0) 
		fParticles[fNParticles]->SetBaryonNumber(-fParticles[fNParticles-1]->GetBaryonNumber()) ;
		else fParticles[fNParticles]->SetBaryonNumber(0) ;
	if(fParticles[fNParticles-1]->GetElectricCharge()!=0)
		fParticles[fNParticles]->SetElectricCharge(-fParticles[fNParticles-1]->GetElectricCharge()) ;
		else fParticles[fNParticles]->SetElectricCharge(0) ;
	if(fParticles[fNParticles-1]->GetStrangeness()!=0)
		fParticles[fNParticles]->SetStrangeness(-fParticles[fNParticles-1]->GetStrangeness()) ;
		else fParticles[fNParticles]->SetStrangeness(0) ;
	if(fParticles[fNParticles-1]->GetIsospinZ()!=0)
		fParticles[fNParticles]->SetIsospinZ(-fParticles[fNParticles-1]->GetIsospinZ()) ;
		else fParticles[fNParticles]->SetIsospinZ(0) ;
	if(fParticles[fNParticles-1]->GetCharm()!=0)
		fParticles[fNParticles]->SetCharm(-fParticles[fNParticles-1]->GetCharm()) ;
		else fParticles[fNParticles]->SetCharm(0) ;
	if(fParticles[fNParticles-1]->GetTruth()!=0)
		fParticles[fNParticles]->SetTruth(-fParticles[fNParticles-1]->GetTruth()) ;
		else fParticles[fNParticles]->SetTruth(0) ;
	if(fParticles[fNParticles-1]->GetBeauty()!=0)
		fParticles[fNParticles]->SetBeauty(-fParticles[fNParticles-1]->GetBeauty()) ;
		else fParticles[fNParticles]->SetBeauty(0) ;

		fParticles[fNParticles]->SetCparity(fParticles[fNParticles-1]->GetCparity()) ; // the same, undefined value!
		if(int(2.*spin) & 1) fParticles[fNParticles]->SetPparity(-fParticles[fNParticles-1]->GetPparity()) ;
			else fParticles[fNParticles]->SetPparity(fParticles[fNParticles-1]->GetPparity()) ;
		fParticles[fNParticles]->SetGparity(fParticles[fNParticles-1]->GetGparity()) ; // ??
	fNParticles++;

    char tname [30] ;
    strcpy(tname, fParticles[fNParticles-2]->GetName()) ;
    switch(fParticles[fNParticles-2]->GetElectricCharge()){
	case 2 : strcat(tname, s_plusplus) ; break ;
	case 1 : strcat(tname, s_plus) ; break ;
	case -1 : strcat(tname, s_minus) ; break ;
	case -2 : strcat(tname, s_minusminus) ; break ;
    }
    fParticles[fNParticles-2]->SetName(tname) ;

    strcpy(tname, fParticles[fNParticles-1]->GetName()) ;
    if(antiparticle[0]=='B') switch(fParticles[fNParticles-1]->GetElectricCharge()){
	case 2 : strcat(tname, s_plusplus) ; break ;
	case 1 : strcat(tname, s_plus) ; break ;
	case -1 : strcat(tname, s_minus) ; break ;
	case -2 : strcat(tname, s_minusminus) ; break ;
    }
    if(antiparticle[0]=='F') switch(fParticles[fNParticles-1]->GetElectricCharge()){
	case 2 : strcat(tname, s_minusminusbar) ; break ;
	case 1 : strcat(tname, s_minusbar) ; break ;
	case 0 : strcat(tname, s_bar) ; break ;
	case -1 : strcat(tname, s_plusbar) ; break ;
	case -2 : strcat(tname, s_plusplusbar) ; break ;
    }
    fParticles[fNParticles-1]->SetName(tname) ;
    }else{	// no antiparticle
    char tname [30] ;
    strcpy(tname, fParticles[fNParticles-1]->GetName()) ;
    switch(fParticles[fNParticles-1]->GetElectricCharge()){
	case 2 : strcat(tname, s_plusplus) ; break ;
	case 1 : strcat(tname, s_plus) ; break ;
	case -1 : strcat(tname, s_minus) ; break ;
	case -2 : strcat(tname, s_minusminus) ; break ;
    }
    fParticles[fNParticles-1]->SetName(tname) ;
}


  }
  particleFile.close();
  if(fNParticles==0) {
    cout << "Warning in DatabasePDG::LoadParticles(): No particles were found in the file specified!!" << endl;
    return kFALSE;
  }
  SortParticles();
  cout << "Particle definitions found : all " << GetNParticles(true)
  << "  good " << GetNParticles() << endl ;
  //    <<  "Good status particles : " << goodStatusParticles << endl;
  return kTRUE;
};


Bool_t DatabasePDG2::LoadDecays() {
  ifstream decayFile;
  decayFile.open(fDecayFilename);
  if(!decayFile) {
    cout << "ERROR in DatabasePDG::LoadDecays() : The ASCII file containing the decays list (\""
         << fDecayFilename << "\") was not found !! Aborting..." << endl;
    return kFALSE;
  }
  
  Int_t mother_pdg, daughter_pdg[3], daughter_found[3];
  Double_t branching;
  Int_t isCG ;
  
  decayFile.exceptions(ios::failbit);
  while(!decayFile.eof()) {
    mother_pdg = 0;
    for(Int_t i=0; i<3; i++) daughter_pdg[i] = 0;
    for(Int_t i=0; i<3; i++) daughter_found[i] = 1;
    branching = -1.0;
    try {
      decayFile >> mother_pdg;
      for(Int_t i=0; i<3; i++) 
        decayFile >> daughter_pdg[i];
      decayFile >> branching >> isCG ;
    }
    catch (ios::failure const &problem) {
    //cout << "reading decays file:\n" ;
      cout << problem.what() << endl;
      break;
    }
    if((mother_pdg!=0) && (daughter_pdg[0]!=0) && (branching>=0)) {
      Int_t nDaughters = 0;
      for(Int_t i=0; i<3; i++)
        if(daughter_pdg[i]!=0){
          nDaughters++;
	  if(!GetPDGParticle(daughter_pdg[i])) daughter_found[i] = 0 ;
	}
      ParticlePDG2* particle = GetPDGParticle(mother_pdg);
      if(particle && daughter_found[0] && daughter_found[1] && daughter_found[2]){
//--------checking if this decay channel is possible
	double energy = particle->GetMass();
	int Q = particle->GetElectricCharge() ;
	int B = particle->GetBaryonNumber() ;
	int S = particle->GetStrangeness() ;
	for(int i=0; i<nDaughters; i++){
	ParticlePDG2* daughter = GetPDGParticle(daughter_pdg[i]) ;
		energy -= daughter->GetMass();
		Q -= daughter->GetElectricCharge() ;
		B -= daughter->GetBaryonNumber() ;
		S -= daughter->GetStrangeness() ;
	}
//	if(energy<0. || Q!=0 || B!=0){
//		cout << "Impossible decay : " << setw(25) << GetPDGParticle(mother_pdg)->GetName() << " -> " << setw(25) 
//		<< GetPDGParticle(daughter_pdg[0])->GetName() << setw(25) << GetPDGParticle(daughter_pdg[1])->GetName() 
//		<< "  Dm = " << energy << " DQ = " << Q << " DB = " << B << endl ;
//		cout << "Impossible decay : " << mother_pdg << " -> " << daughter_pdg[0]
//		<< " " << daughter_pdg[1]  << " " << daughter_pdg[2] << endl ;
//	}
//--------
      if(isCG){
	ParticlePDG2* daughter1 = GetPDGParticle(daughter_pdg[0]) ;
	ParticlePDG2* daughter2 = GetPDGParticle(daughter_pdg[1]) ;
	branching *= pow(Clebsch_Gordan(daughter1->GetIsospin(), daughter1->GetIsospinZ(),
	daughter2->GetIsospin(), daughter2->GetIsospinZ(),particle->GetIsospin(), particle->GetIsospinZ()) ,2);
	}
      DecayChannel decay(mother_pdg, branching, nDaughters, daughter_pdg);
      if(energy>0. && Q==0 && B==0) particle->AddChannel(&decay);
	}
    }
  }
  decayFile.close();
  Int_t nDecayChannels = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    nDecayChannels += fParticles[i]->GetNDecayChannels();
    if(fParticles[i]->GetNDecayChannels()==0 && fParticles[i]->GetWidth()>0. &&
    !isStable(fParticles[i]->GetPDG())){
    fStatus[i] = kFALSE ;
//    cout << "suspect particle : " << fParticles[i]->GetName() << endl ;
  }
  }
  SortParticles() ;
  cout << "Number of decays found in the database : " << nDecayChannels << endl;
  return kTRUE;
};

ParticlePDG2* DatabasePDG2::GetPDGParticleByIndex(Int_t index) {
  if(index<0 || index>fNParticles) {
    cout << "Warning in DatabasePDG::GetPDGParticleByIndex(Int_t): Particle index is negative or too big !!" << endl
         << " It must be inside this range: [0, " << fNParticles-1 << "]" << endl
         << " Returning null pointer!!" << endl;
    return 0x0;
  }
  return fParticles[index];
}

Bool_t DatabasePDG2::GetPDGParticleStatusByIndex(Int_t index) {
  if(index<0 || index>fNParticles) {
    cout << "Warning in DatabasePDG::GetPDGParticleStatusByIndex(Int_t): Particle index is negative or too big !!" << endl
         << " It must be inside this range: [0, " << fNParticles-1 << "]" << endl
         << " Returning null pointer!!" << endl;
    return kFALSE;
  }
  return fStatus[index];
}

ParticlePDG2* DatabasePDG2::GetPDGParticle(Int_t pdg) {
  Int_t nFindings = 0;
  Int_t firstTimeIndex = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    if(pdg == fParticles[i]->GetPDG()) {
      if(nFindings == 0) firstTimeIndex = i;
      nFindings++;
    }
  }
  if(nFindings == 1) return fParticles[firstTimeIndex];
  if(nFindings == 0) {
    //cout << "Warning in DatabasePDG::GetPDGParticle(Int_t): The particle required with PDG = " << pdg
    //     << " was not found in the database!!" << endl;
    return 0x0;
  }
  if(nFindings >= 2) {
    cout << "Warning in DatabasePDG::GetPDGParticle(Int_t): The particle required with PDG = " << pdg
         << " was found with " << nFindings << " entries in the database. Check it out !!" << endl
	 << "Returning the first instance found" << endl;
    return fParticles[firstTimeIndex];
  }
  return 0x0;
};

ParticlePDG2* DatabasePDG2::GetPDGParticle(Int_t pdg, Int_t& firstTimeIndex) {
  Int_t nFindings = 0;
  firstTimeIndex = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    if(pdg == fParticles[i]->GetPDG()) {
      if(nFindings == 0) firstTimeIndex = i;
      nFindings++;
    }
  }
  if(nFindings == 1) return fParticles[firstTimeIndex];
  if(nFindings == 0) {
    //cout << "Warning in DatabasePDG::GetPDGParticle(Int_t): The particle required with PDG = " << pdg
    //     << " was not found in the database!!" << endl;
    return 0x0;
  }
  if(nFindings >= 2) {
    cout << "Warning in DatabasePDG::GetPDGParticle(Int_t): The particle required with PDG = " << pdg
         << " was found with " << nFindings << " entries in the database. Check it out !!" << endl
	 << "Returning the first instance found" << endl;
    return fParticles[firstTimeIndex];
  }
  return 0x0;
};


int DatabasePDG2::GetIndex(Int_t pdg) {
  Int_t nFindings = 0;
  Int_t firstTimeIndex = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    if(pdg == fParticles[i]->GetPDG()) {
      if(nFindings == 0) firstTimeIndex = i;
      nFindings++;
    }
  }
  if(nFindings == 1) return firstTimeIndex;
  if(nFindings == 0) {
    //cout << "Warning in DatabasePDG::GetPDGParticle(Int_t): The particle required with PDG = " << pdg
    //     << " was not found in the database!!" << endl;
    return -1;
  }
  if(nFindings >= 2) {
    cout << "Warning in DatabasePDG::GetPDGParticle(Int_t): The particle required with PDG = " << pdg
         << " was found with " << nFindings << " entries in the database. Check it out !!" << endl
	 << "Returning the first instance found" << endl;
    return firstTimeIndex;
  }
  return -1;
};


int DatabasePDG2::GetPionIndex()
{ 
  char* piname="pi-" ;
  return GetIndex(GetPDGParticle(piname)->GetPDG()) ; 
}

Bool_t DatabasePDG2::GetPDGParticleStatus(Int_t pdg) {
  Int_t nFindings = 0;
  Int_t firstTimeIndex = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    if(pdg == fParticles[i]->GetPDG()) {
      if(nFindings == 0) firstTimeIndex = i;
      nFindings++;
    }
  }
  if(nFindings == 1) return fStatus[firstTimeIndex];
  if(nFindings == 0) {
    //cout << "Warning in DatabasePDG::GetPDGParticle(Int_t): The particle required with PDG = " << pdg
    //     << " was not found in the database!!" << endl;
    return kFALSE;
  }
  if(nFindings >= 2) {
    cout << "Warning in DatabasePDG::GetPDGParticleStatus(Int_t): The particle status required for PDG = " << pdg
         << " was found with " << nFindings << " entries in the database. Check it out !!" << endl
	 << "Returning the status of first instance found" << endl;
    return fStatus[firstTimeIndex];
  }
  return kFALSE;
};

ParticlePDG2* DatabasePDG2::GetPDGParticle(Char_t* name) {
  Int_t nFindings = 0;
  Int_t firstTimeIndex = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    if(!strcmp(name, fParticles[i]->GetName())) {
      if(nFindings == 0) firstTimeIndex = i;
      nFindings++;
    }
  }
  if(nFindings == 1) return fParticles[firstTimeIndex];
  if(nFindings == 0) {
    //cout << "Warning in DatabasePDG::GetPDGParticle(Char_t*): The particle required with name \"" << name
    //     << "\" was not found in the database!!" << endl;
    return 0x0;
  }
  if(nFindings >= 2) {
    cout << "Warning in DatabasePDG::GetPDGParticle(Char_t*): The particle required with name \"" << name
         << "\" was found with " << nFindings << " entries in the database. Check it out !!" << endl
	 << "Returning the first instance found" << endl;
    return fParticles[firstTimeIndex];
  }
  return 0x0;
};

Bool_t DatabasePDG2::GetPDGParticleStatus(Char_t* name) {
  Int_t nFindings = 0;
  Int_t firstTimeIndex = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    if(!strcmp(name, fParticles[i]->GetName())) {
      if(nFindings == 0) firstTimeIndex = i;
      nFindings++;
    }
  }
  if(nFindings == 1) return fStatus[firstTimeIndex];
  if(nFindings == 0) {
    //cout << "Warning in DatabasePDG::GetPDGParticle(Char_t*): The particle required with name \"" << name
    //     << "\" was not found in the database!!" << endl;
    return kFALSE;
  }
  if(nFindings >= 2) {
    cout << "Warning in DatabasePDG::GetPDGParticleStatus(Char_t*): The particle status required for name \"" << name
         << "\" was found with " << nFindings << " entries in the database. Check it out !!" << endl
	 << "Returning the first instance found" << endl;
    return fStatus[firstTimeIndex];
  }
  return kFALSE;
};

void DatabasePDG2::DumpData(Bool_t dumpAll) {
  ofstream fout("dump.dat") ;
  fout << setw(15) << "name" << setw(12) << "mass" << setw(12) << "width" << setw(5) << "J" << setw(5) << "I" << setw(5) << "I_3" 
	<< setw(12) << "pid"
	<< setw(4) << "B" << setw(4) << "Q" << setw(4) << "S" << setw(4) << "C" << setw(4) << "T" << setw(4) << "B" << endl ;
  for(int i=0; i<GetNParticles(); i++){
	ParticlePDG2* p = GetPDGParticleByIndex(i) ;
	fout << setw(15) << p->GetName() << setw(12) << p->GetMass() << setw(12) << p->GetWidth() << setw(5) << p->GetSpin() 
	<< setw(5) << p->GetIsospin() << setw(5) << p->GetIsospinZ() << setw(12) << p->GetPDG() << setw(4) << p->GetBaryonNumber() 
	<< setw(4) << p->GetElectricCharge() << setw(4) << p->GetStrangeness() << setw(4) << p->GetCharm() 
	<< setw(4) << p->GetTruth() << setw(4) << p->GetBeauty() << endl;
  }
  fout.close() ;
}

void DatabasePDG2::SetUseCharmParticles(Bool_t flag) {
  if(fNParticles>0) {
    fUseCharmParticles = flag;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetCharm()>0.)		  
	fStatus[i] = flag;
    }
    SortParticles();
    return;
  }
  else
    fUseCharmParticles = flag;
  return;
};

void DatabasePDG2::SetMinimumWidth(Double_t value) {
  if(fNParticles>0) {
    fMinimumWidth = value;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetWidth() < fMinimumWidth)		  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else
    fMinimumWidth = value;
  return;
};

void DatabasePDG2::SetMaximumWidth(Double_t value) {
  if(fNParticles>0) {
    fMaximumWidth = value;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetWidth() > fMaximumWidth)		  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else
    fMaximumWidth = value;
  return;
};

void DatabasePDG2::SetWidthRange(Double_t min, Double_t max) {
  if(fNParticles>0) {
    fMinimumWidth = min;
    fMaximumWidth = max;
    for(Int_t i=0; i<fNParticles; i++) {
      if((fParticles[i]->GetWidth()<fMinimumWidth) || (fParticles[i]->GetWidth()>fMaximumWidth))  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else {
    fMinimumWidth = min;
    fMaximumWidth = max;
  }
  return;
};

void DatabasePDG2::SetMinimumMass(Double_t value) {
  if(fNParticles>0) {
    fMinimumMass = value;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetMass() < fMinimumMass)		  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else
    fMinimumMass = value;
  return;
};

void DatabasePDG2::SetMaximumMass(Double_t value) {
  if(fNParticles>0) {
    fMaximumMass = value;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetMass() > fMaximumMass)		  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else
    fMaximumMass = value;
  return;
};

void DatabasePDG2::SetMassRange(Double_t min, Double_t max) {
  if(fNParticles>0) {
    fMinimumMass = min;
    fMaximumMass = max;
    for(Int_t i=0; i<fNParticles; i++) {
      if((fParticles[i]->GetMass()<fMinimumMass) || (fParticles[i]->GetMass()>fMaximumMass))  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else {
    fMinimumMass = min;
    fMaximumMass = max;
  }
  return;
};

void DatabasePDG2::SortParticles() {
  if(fNParticles<2) {
    cout << "Warning in DatabasePDG::SortParticles() : No particles to sort. Load data first!!" << endl;
    return;
  }

  Int_t nGoodStatus = 0;
  for(Int_t i=0; i<fNParticles; i++)
    if(fStatus[i]) nGoodStatus++;
  if(nGoodStatus==fNParticles)    // if all particles have good status then there is nothing to do
    return;
  if(nGoodStatus==0)              // no good status particles, again nothing to do
    return;

  Int_t shifts = 1;
  while(shifts) {
    shifts = 0;
    for(Int_t i=0; i<fNParticles-1; i++) {
      if(!fStatus[i] && fStatus[i+1]) {   // switch if false status is imediately before a true status particle
	ParticlePDG2 *temporaryPointer = fParticles[i];
	fParticles[i] = fParticles[i+1];
	fParticles[i+1] = temporaryPointer;
	Bool_t temporaryStatus = fStatus[i];
	fStatus[i] = fStatus[i+1];
	fStatus[i+1] = temporaryStatus;
	shifts++;
      }
    }
  }
  return;
}


void DatabasePDG2::SortParticlesByMass() {
  if(fNParticles<2) {
    cout << "Warning in DatabasePDG::SortParticles() : No particles to sort. Load data first!!" << endl;
    return;
  }

  Int_t shifts = 1;
  while(shifts) {
    shifts = 0;
    for(Int_t i=0; i<GetNParticles()-1; i++) {
      if(fParticles[i]->GetMass()<fParticles[i+1]->GetMass()) {
	ParticlePDG2 *temporaryPointer = fParticles[i];
	fParticles[i] = fParticles[i+1];
	fParticles[i+1] = temporaryPointer;
	Bool_t temporaryStatus = fStatus[i];
	fStatus[i] = fStatus[i+1];
	fStatus[i+1] = temporaryStatus;
	shifts++;
      }
    }
  }
  return;
}


void DatabasePDG2::CorrectBranching(){
 for(int i=0; i<fNParticles; i++){
   ParticlePDG2 *part = GetPDGParticleByIndex(i) ;
   double fb = part->GetFullBranching() ;
//   if(fb<1.) cout << "EEEfb  "  << setw(24) << part->GetName() << setw(14) << fb << setw(14) << part->GetNDecayChannels() << endl ;
   if(fb<1.) for(int id=0; id<part->GetNDecayChannels(); id++){
    double b = part->GetDecayChannel(id)->GetBranching() ;
    part->GetDecayChannel(id)->SetBranching(b/fb) ;
   }
 }
}


Int_t DatabasePDG2::GetNParticles(Bool_t all) {
  if(all)
    return fNParticles;

  Int_t nGoodStatus = 0;
  for(Int_t i=0; i<fNParticles; i++)
    if(fStatus[i]) nGoodStatus++;
  return nGoodStatus;
}


Bool_t DatabasePDG2::IsChannelAllowed(DecayChannel *channel, Double_t motherMass) {
  Double_t daughtersSumMass = 0.0;
  for(Int_t i=0; i<channel->GetNDaughters(); i++)
    daughtersSumMass += GetPDGParticle(channel->GetDaughterPDG(i))->GetMass();
  if(daughtersSumMass<motherMass)
    return kTRUE;
  return kFALSE;
}

Int_t DatabasePDG2::GetNAllowedChannels(ParticlePDG2 *particle, Double_t motherMass) {
  Int_t nAllowedChannels = 0;
  for(Int_t i=0; i<particle->GetNDecayChannels(); i++)
    nAllowedChannels += (IsChannelAllowed(particle->GetDecayChannel(i), motherMass) ? 1:0);
  return nAllowedChannels;
}


bool isStable(int pdg)
{
  if(weakContribution){
  if(abs(pdg)==211 || abs(pdg)==111 || abs(pdg)==2112 || abs(pdg)==2212 || abs(pdg)==311 || abs(pdg)==321) return true ; // +3212 originally, +3122 for Lambda/K0 analysis
  else return false ;
  }
  else{
  if(abs(pdg)==211 || abs(pdg)==111 || abs(pdg)==2112 || abs(pdg)==2212 || abs(pdg)==311 || abs(pdg)==321 || abs(pdg)==9000221
  || abs(pdg)==3322 || abs(pdg)==3312 || abs(pdg)==3122 || /*abs(pdg)==3222 ||*/ abs(pdg)==3212 || abs(pdg)==3112 || abs(pdg)==3334) return true ;
  else return false ;
  }
}


void DatabasePDG2::SortDecayingResonances() {
    for(Int_t i=0; i<fNParticles; i++) {
      if(!isStable(fParticles[i]->GetPDG()) && fParticles[i]->GetNDecayChannels()==0)		  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
};

