// Scorer for HGMTNTuple
//
// ********************************************************************
// *                                                                  *
// *                                                                  *
// * This file was obtained from Topas MC Inc under the license       *
// * agreement set forth at http://www.topasmc.org/registration       *
// * Any use of this file constitutes full acceptance of              *
// * this TOPAS MC license agreement.                                 *
// *                                                                  *
// ********************************************************************
//

#include "HGMTNTuple.hh"

#include "G4Event.hh"
#include "G4PSDirectionFlag.hh"
#include "G4VProcess.hh"
#include "TsTrackInformation.hh"

HGMTNTuple::HGMTNTuple(TsParameterManager *pM, TsMaterialManager *mM,
                       TsGeometryManager *gM, TsScoringManager *scM,
                       TsExtensionManager *eM, G4String scorerName,
                       G4String quantity, G4String outFileName,
                       G4bool isSubScorer)
    : TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName,
                      isSubScorer) {
  // SetSurfaceScorer();

  fNtuple->RegisterColumnI(&fEvent, "Event Number");
  fNtuple->RegisterColumnD(&fEnergy, "Energy", "keV");
  fNtuple->RegisterColumnF(&fPosX, "Position X", "cm");
  fNtuple->RegisterColumnF(&fPosY, "Position Y", "cm");
  fNtuple->RegisterColumnF(&fPosZ, "Position Z", "cm");
  // fNtuple->RegisterColumnF(&fMomentumX, "Momentum X", "");
  // fNtuple->RegisterColumnF(&fMomentumY, "Momentum Y", "");
  // fNtuple->RegisterColumnF(&fMomentumZ, "Momentum Z", "");
  //  fNtuple->RegisterColumnF(&fWeight, "Weight", "");
  fNtuple->RegisterColumnF(&fTimeOfFlight, "Time of Flight", "ns");
  // fNtuple->RegisterColumnI(&fParticleType, "Particle Type (in PDG Format)");
  // fNtuple->RegisterColumnS(&fOriginProcessName, "Origin Process");
  // fNtuple->RegisterColumnI(&fOriginProcessID, "Origin Process (int)");
  fNtuple->RegisterColumnI(&fParentID, "Parent ID");
  fNtuple->RegisterColumnI(&fTrackID, "Particle ID");
}

HGMTNTuple::~HGMTNTuple() { ; }

G4bool HGMTNTuple::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
  fParticleType = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  fParentID = aStep->GetTrack()->GetParentID();
  if (!fIsActive || fParticleType != 11 || (fParentID != 2 && fParentID != 3) ||
      aStep->GetTrack()->GetCurrentStepNumber() != 1) {
    fSkippedWhileInactive++;
    return false;
  }
  fEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
  ResolveSolid(aStep);

  G4StepPoint *theStepPoint = 0;
  // this has weird units or something
  // fDeposit = aStep->GetTotalEnergyDeposit();
  fTimeOfFlight = aStep->GetTrack()->GetGlobalTime();

  // fEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
  // G4ThreeVector momentum = aStep->GetPreStepPoint()->GetMomentumDirection();

  fPosX = pos.x();
  fPosY = pos.y();
  fPosZ = pos.z();
  // fMomentumX = momentum.x();
  // fMomentumY = momentum.y();
  // fMomentumZ = momentum.z();

  fEvent = GetEventID();

  fTrackID = aStep->GetTrack()->GetTrackID();

  const G4VProcess *originProcess = aStep->GetTrack()->GetCreatorProcess();

  // Check if this is a new history
  // fRunID   = GetRunID();
  // fEventID = GetEventID();
  // if (fEventID != fPrevEventID || fRunID != fPrevRunID) {
  // 	fIsNewHistory = true;
  // 	fPrevEventID = fEventID;
  // 	fPrevRunID = fRunID;
  // } else {
  // 	fIsNewHistory = false;
  // }

  fNtuple->Fill();
  return true;
}
