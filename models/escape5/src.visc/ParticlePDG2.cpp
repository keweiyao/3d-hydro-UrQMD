
#ifndef PARTICLE_PDG2
#include "ParticlePDG2.h"
#endif

#include <iostream>
using std::cout;
using std::endl;

ParticlePDG2::ParticlePDG2() {
  fPDG   = kNonsensePDG;
  fMass  = -1.0;
  fWidth = 0.0;
  fNDecayChannels = 0;
  for(Int_t i=0; i<kMaxDecayChannels2; i++)
    fDecayChannels[i] = new DecayChannel();
  ffMax = 0.0 ; fDensity = 0.0 ;
};

ParticlePDG2::ParticlePDG2(Char_t *name, Int_t pdg, Double_t mass, Double_t width) {
  for(Int_t i=0; i<9; i++)
    if(*(name+i) != '\0') fName[i] = *(name+i);
    else break;
  fPDG   = pdg;
  fMass  = mass;
  fWidth = width;
  fNDecayChannels = 0;
  for(Int_t i=0; i<kMaxDecayChannels2; i++)
    fDecayChannels[i] = new DecayChannel();
};

ParticlePDG2::~ParticlePDG2() {
  for(Int_t i=0; i<kMaxDecayChannels2; i++)
    delete fDecayChannels[i];
};

Double_t ParticlePDG2::GetFullBranching() {
  Double_t fullBranching = 0.0;
  for(Int_t i=0; i<fNDecayChannels; i++)
    fullBranching += fDecayChannels[i]->GetBranching();
  return fullBranching;
};

void ParticlePDG2::AddChannel(DecayChannel* channel) {
  if(channel->GetMotherPDG() != fPDG) {
    cout << "ERROR in ParticlePDG::AddChannel() : You try to add a channel which has a different mother PDG" << endl;
    return;
  }
  fDecayChannels[fNDecayChannels]->SetMotherPDG(channel->GetMotherPDG());
  fDecayChannels[fNDecayChannels]->SetBranching(channel->GetBranching());
  fDecayChannels[fNDecayChannels]->SetDaughters(channel->GetDaughters(), channel->GetNDaughters());
  fNDecayChannels++;
};
