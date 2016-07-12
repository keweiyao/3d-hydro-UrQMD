
#ifndef PARTICLE_PDG2
#define PARTICLE_PDG2

#include "Rtypes.h"
#include <cstring>

#ifndef DECAY_CHANNEL
#include "DecayChannel.h"
#endif

const Int_t kMaxDecayChannels2 = 100;

class ParticlePDG2 {
 private:
  Char_t        fName[30];
  Int_t         fPDG;
  Double_t      fMass;
  Double_t      fWidth;
  Double_t      fSpin;                         // J
  Double_t      fIsospin;                      // I
  Double_t      fIsospinZ;                     // I3
  Double_t      ffMax;                      // f_max, for event generation
  Double_t      fDensity;                   // thermal density, for event generation
  Int_t      fBaryonNumber;                 // 
  Int_t      fElectricCharge;               // 
  Int_t      fStrangeness;                  // 
  Int_t      fCharm;                        // 
  Int_t      fBeauty;                       // 
  Int_t      fTruth;                        //
  Int_t	     fStrangeQuarkNumber;
  Int_t      fAntiStrangeQuarkNumber;
  
  Int_t fGparity;
  Int_t fPparity;
  Int_t fCparity;
  Char_t fPDGstatus;
  
  Int_t         fNDecayChannels;
  DecayChannel* fDecayChannels[kMaxDecayChannels2];
   
 public:
  ParticlePDG2();
  ParticlePDG2(Char_t* name, Int_t pdg, Double_t mass, Double_t width);
  ~ParticlePDG2();
  
  void AddChannel(DecayChannel* channel);
  void SetName(Char_t* name) {
    strcpy(fName, name);
  }
  void SetPDG(Int_t value) {fPDG = value;}
  void SetMass(Double_t value) {fMass = value;}
  void SetWidth(Double_t value) {fWidth = value;}
  void SetSpin(Double_t value) {fSpin = value;}
  void SetIsospin(Double_t value) {fIsospin = value;}
  void SetIsospinZ(Double_t value) {fIsospinZ = value;}
  void SetBaryonNumber(Int_t value) {fBaryonNumber = value;}
  void SetElectricCharge(Int_t value) {fElectricCharge = value;}
  void SetStrangeness(Int_t value) {fStrangeness = value;}
  void SetStrangeQNumber(Int_t value) {fStrangeQuarkNumber = value;}
  void SetStrangeAQNumber(Int_t value) {fAntiStrangeQuarkNumber = value;}
  void SetCharm(Int_t value) {fCharm = value;}
  void SetBeauty(Int_t value) {fBeauty = value;}
  void SetTruth(Int_t value) {fTruth = value;}
  void SetFMax(Double_t value) {ffMax = value;}
  void SetDensity(Double_t value) {fDensity = value;}
  
  void SetGparity(Int_t value) {fGparity = value; }
  void SetCparity(Int_t value) {fCparity = value; }
  void SetPparity(Int_t value) {fPparity = value; }
  void SetPDGstatus(Char_t value) {fPDGstatus = value; }
  
  Char_t* GetName() {return fName;}
  Int_t GetPDG() {return fPDG;}
  Double_t GetMass() {return fMass;}
  Double_t GetWidth() {return fWidth;}
  Int_t GetNDecayChannels() {return fNDecayChannels;}
  Double_t GetSpin() {return fSpin;}
  Double_t GetIsospin() {return fIsospin;}
  Double_t GetIsospinZ() {return fIsospinZ;}
  Int_t GetBaryonNumber() {return fBaryonNumber;}
  Int_t GetElectricCharge() {return fElectricCharge;}
  Int_t GetStrangeness() {return fStrangeness;}
  Int_t GetStrangeQNumber() {return fStrangeQuarkNumber;}
  Int_t GetStrangeAQNumber() {return fAntiStrangeQuarkNumber;}
  Int_t GetCharm() {return fCharm;}
  Int_t GetBeauty() {return fBeauty;}
  Int_t GetTruth() {return fTruth;}
  double GetFMax() {return ffMax;}
  double GetDensity() {return fDensity;}
  
  Int_t GetGparity() {return fGparity;}
  Int_t GetCparity() {return fCparity;}
  Int_t GetPparity() {return fPparity;}
  Char_t GetPDGstatus() {return fPDGstatus;}
  
  Double_t GetFullBranching();
  DecayChannel* GetDecayChannel(Int_t i) {
    if(0<=i && i<fNDecayChannels) 
      return fDecayChannels[i];
    else
      return 0x0;
  }
};

#endif
