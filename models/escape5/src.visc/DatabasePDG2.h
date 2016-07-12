
#ifndef DATABASE_PDG2
#define DATABASE_PDG2

#include "Rtypes.h"
#ifndef PARTICLE_PDG2
#include "ParticlePDG2.h"
#endif

const Int_t kMaxParticles2 = 1000;

bool isStable(int pdg) ;

class DatabasePDG2 {
 private:
  Int_t fNParticles;                        // no. of particles in database
  ParticlePDG2 *fParticles[kMaxParticles2];   // array of particle pointers
  Bool_t fStatus[kMaxParticles2];            // status of each particle
  Char_t fParticleFilename[256];            // particle list filename
  Char_t fDecayFilename[256];               // decay channels filename
  Bool_t fUseCharmParticles;                // flag for using (or not) charm particles
  Double_t fMinimumWidth;                   // minimum allowed width for resonances
  Double_t fMaximumWidth;                   // maximum allowed width for resonances
  Double_t fMinimumMass;                    // minimum allowed mass for resonances
  Double_t fMaximumMass;                    // maximum allowed mass for resonances

  Bool_t LoadParticles();
  Bool_t LoadDecays();
  void SortParticles();                     // put the good status particles at the beggining of the list

 public:
  DatabasePDG2(char *fileparticles, char *filedecay);
  ~DatabasePDG2();
  void SortParticlesByMass();
  void CorrectBranching() ;
  int GetPionIndex() ;
  // Load the particle PDG information from the particle and decay files
  Bool_t LoadData();                        
  
  // Set particle and decay filenames
  void SetParticleFilename(Char_t *filename);
  void SetDecayFilename(Char_t *filename);
  // Set criteria for using particles. Those particle which do not match
  // these criteria will be flagged with FALSE in the fStatus array.
  void SetUseCharmParticles(Bool_t flag);
  void SetMinimumWidth(Double_t value);
  void SetMaximumWidth(Double_t value);
  void SetMinimumMass(Double_t value);
  void SetMaximumMass(Double_t value);
  void SetWidthRange(Double_t min, Double_t max);
  void SetMassRange(Double_t min, Double_t max);
  
  
  Char_t* GetParticleFilename() {return fParticleFilename;}
  Char_t* GetDecayFilename() {return fDecayFilename;}
  Int_t GetNParticles(Bool_t all = kFALSE);      // true - no. of all particles; false - no. of good status particles
  ParticlePDG2* GetPDGParticleByIndex(Int_t index);
  Bool_t GetPDGParticleStatusByIndex(Int_t index);
  ParticlePDG2* GetPDGParticle(Int_t pdg);
  ParticlePDG2* GetPDGParticle(Int_t pdg, Int_t& index);
  int GetIndex(Int_t pdg) ;
  Bool_t GetPDGParticleStatus(Int_t pdg);
  ParticlePDG2* GetPDGParticle(Char_t *name);
  Bool_t GetPDGParticleStatus(Char_t *name);
  Bool_t GetUseCharmParticles() {return fUseCharmParticles;};
  Double_t GetMinimumWidth() {return fMinimumWidth;};
  Double_t GetMaximumWidth() {return fMaximumWidth;};
  Double_t GetMinimumMass() {return fMinimumMass;};
  Double_t GetMaximumMass() {return fMaximumMass;};
  void DumpData(Bool_t dumpAll = kFALSE); // print the PDG information in the console
  Bool_t IsChannelAllowed(DecayChannel *channel, Double_t motherMass);
  Int_t GetNAllowedChannels(ParticlePDG2 *particle, Double_t motherMass);
  void SortDecayingResonances();
};


#endif
