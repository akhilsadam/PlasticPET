//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm4/include/SteppingAction.hh
/// \brief Definition of the SteppingAction class
//
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "vector"

class B3aEventAction;
class G4VUserDetectorConstruction;
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class B3SteppingAction : public G4UserSteppingAction
{
  public:
    B3SteppingAction(B3aEventAction*,G4VUserDetectorConstruction*);
   ~B3SteppingAction();
    static G4int id;
    virtual void UserSteppingAction(const G4Step*);
    const vector<string> statusEnum = {"Undefined","FresnelRefraction","FresnelReflection","TotalInternalReflection","LambertianReflection","LobeReflection","SpikeReflection","BackScattering","Absorption","Detection","NotAtBoundary","SameMaterial","StepTooSmall","NoRINDEX","PolishedLumirrorAirReflection","PolishedLumirrorGlueReflection","PolishedAirReflection","PolishedTeflonAirReflection","PolishedTiOAirReflection","PolishedTyvekAirReflection","PolishedVM2000AirReflection","PolishedVM2000GlueReflection","EtchedLumirrorAirReflection","EtchedLumirrorGlueReflection","EtchedAirReflection","EtchedTeflonAirReflection","EtchedTiOAirReflection","EtchedTyvekAirReflection","EtchedVM2000AirReflection","EtchedVM2000GlueReflection","GroundLumirrorAirReflection","GroundLumirrorGlueReflection","GroundAirReflection","GroundTeflonAirReflection","GroundTiOAirReflection","GroundTyvekAirReflection","GroundVM2000AirReflection","GroundVM2000GlueReflection","Dichroic"};
    private:
    B3aEventAction* fEventAction;
    G4VUserDetectorConstruction* fpatient;
    
};

#endif
