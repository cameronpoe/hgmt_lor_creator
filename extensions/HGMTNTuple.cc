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
#include <cstdlib>
#include <vector>

#include "G4Event.hh"
#include "G4PSDirectionFlag.hh"
#include "G4ThreeVector.hh"
#include "G4TrackStatus.hh"
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
  // fNtuple->RegisterColumnF(&fParentMomentumX, "Momentum X", "");
  // fNtuple->RegisterColumnF(&fParentMomentumY, "Momentum Y", "");
  // fNtuple->RegisterColumnF(&fParentMomentumZ, "Momentum Z", "");
  fNtuple->RegisterColumnF(&fTimeOfFlight, "Time of Flight", "ns");
  // fNtuple->RegisterColumnS(&fOriginProcessName, "Origin Process");
  // fNtuple->RegisterColumnI(&fOriginProcessID, "Origin Process (int)");
  fNtuple->RegisterColumnI(&fParentID, "Parent ID");
  fNtuple->RegisterColumnI(&fTrackID, "Track ID");

  // fNtuple->RegisterColumnI(&fParticleType, "Particle Type");
  // fNtuple->RegisterColumnI(&fStepNumber, "Step Number");
}

HGMTNTuple::~HGMTNTuple() { ; }

G4bool HGMTNTuple::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
  G4Track *track = aStep->GetTrack();
  G4TrackStatus trackStatus = track->GetTrackStatus();
  fParticleType = track->GetDefinition()->GetPDGEncoding();
  fParentID = track->GetParentID();
  fStepNumber = track->GetCurrentStepNumber();
  if ((fStepNumber != 1 || fParticleType != 11 ||
       (fParentID != 2 && fParentID != 3)) &&
      (fParticleType != -11 || fParentID != 0 ||
       (trackStatus != fStopAndKill && trackStatus != fStopButAlive))) {
    fSkippedWhileInactive++;
    return false;
  }
  fEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
  ResolveSolid(aStep);

  fTimeOfFlight = track->GetGlobalTime();

  G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition();
  // G4ThreeVector momentum = aStep->GetPreStepPoint()->GetMomentumDirection();

  fPosX = pos.x();
  fPosY = pos.y();
  fPosZ = pos.z();

  fEvent = GetEventID();

  fTrackID = track->GetTrackID();
  fNtuple->Fill();
  return true;
}
