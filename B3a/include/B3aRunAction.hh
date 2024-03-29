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
//
/// \file B3aRunAction.hh
/// \brief Definition of the B3aRunAction class

#ifndef B3aRunAction_h
#define B3aRunAction_h 1

#include "B3Analysis.hh"

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include "RunTools.hh"

class HistoManager;
class DetectorConstruction;
class B3PrimaryGeneratorAction;
/// Run action class

class B3aRunAction : public G4UserRunAction
{
  public:
    B3aRunAction(DetectorConstruction* patient,B3PrimaryGeneratorAction* kin);
    ~B3aRunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void CountEvent()           { fGoodEvents += 1; RunTool->completedEvent(); };
    void SumDose(G4double dose) { fSumDose += dose; };  

    G4int nevents;
    inline G4int GetNevents() { return nevents;}

private:
    G4Accumulable<G4int>    fGoodEvents;
    G4Accumulable<G4double> fSumDose;  
    HistoManager* fHistoManager;
    DetectorConstruction* fpatient;
    B3PrimaryGeneratorAction* fPrimary;
    static bool CrossSectionTrue;
    static RunTools* RunTool;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

