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
/// \file B3aRunAction.cc
/// \brief Implementation of the B3aRunAction class
#include "Version.hh"
#include "B3aRunAction.hh"
#include "B3PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "B3aRun.hh"//"G4Run.hh"
#include "G4AccumulableManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "B3aHistoManager.hh"
//#include "B3DetectorConstruction.hh"
#include "GDMLDetectorConstruction.hh"

#ifdef CrossSectionTest
	#include "G4ProcessManager.hh"
	#include "G4LossTableManager.hh"
	#include "G4PhysicalConstants.hh"
	#include "G4SystemOfUnits.hh"
	#include "G4EmCalculator.hh"
	#include "CrossSectionTester.hh"
#endif

#include <time.h>

std::mutex foo21;
std::mutex barL21;
bool B3aRunAction::CrossSectionTrue = true;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
B3aRunAction::B3aRunAction(G4VUserDetectorConstruction* patient, B3PrimaryGeneratorAction* kin)
 : G4UserRunAction(),
   fGoodEvents(0),
   fSumDose(0.),
   fpatient(patient),
   fHistoManager(0)
{  
  //add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fGoodEvents);
  accumulableManager->RegisterAccumulable(fSumDose); 
  fHistoManager = new HistoManager(fpatient);
  fPrimary = kin;
  CrossSectionTrue = false;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aRunAction::~B3aRunAction()
{
   delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aRunAction::BeginOfRunAction(const G4Run* run)
{ 
  #ifdef CrossSectionTest
  	std::lock(foo21,barL21);
    if(!CrossSectionTrue)
    {
    	CrossSectionTester* cst = new CrossSectionTester();
    	cst->CSRunAction(run, (GDMLDetectorConstruction*) fpatient, fPrimary);
      G4cout << "### CROSS SECTIONS CREATED " << CrossSectionTrue << G4endl;
      G4cout << "### UPDATING MPT ..." << G4endl;
      ((GDMLDetectorConstruction*) fpatient)->UPDATE_GEO_MPT();
      G4cout << "### MPT UPDATED" << G4endl;
      CrossSectionTrue = true;
    }
    foo21.unlock();
		barL21.unlock();
  #endif


  nevents = (G4int) (run->GetNumberOfEventToBeProcessed());
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  
  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }  
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aRunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();
  // save histograms

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {  
    G4double interacted = (analysisManager->GetH1(10)->bin_entries(0));
    G4cout << "INTERACTED " <<interacted <<G4endl;
    if (IsMaster() && (interacted !=0))
  {
    G4int nx = 3;
    G4int ny = 16;
    G4double nStrips = 1.0*nx*ny;
    analysisManager->GetH2(6)->multiply(1.0/interacted);
    analysisManager->GetH2(12)->multiply(1.0/interacted);
    analysisManager->GetH2(13)->multiply(1.0/interacted);
    //analysisManager->GetH2(23)->multiply(nStrips);
    //analysisManager->GetH2(24)->multiply(nStrips);
  }
    //#ifdef CST2
    //analysisManager->GetH2(19)->multiply(1.0/nofEvents);
    //#endif
    //close
    analysisManager->Write();    
    analysisManager->CloseFile();
  }    
  


  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B3PrimaryGeneratorAction* generatorAction
    = static_cast<const B3PrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String partName;
  if (generatorAction) 
  {
    G4ParticleDefinition* particle 
      = generatorAction->GetParticleGun()->GetParticleDefinition();
    partName = particle->GetParticleName();
  }  
          
  // Print results
  //
  if (IsMaster())
  {
    #ifdef SingleStrip
      #ifndef LEGEND
        G4cout
        << G4endl
        << "--------------------End of Global Run------------------------"
        << G4endl
        << "  The run was " << (nofEvents*100) << " "<< partName;
        #else
        G4cout
        << G4endl
        << "--------------------End of Global Run------------------------"
        << G4endl
        << "  The run was " << (nofEvents) << " "<< partName;
      #endif
    #endif
    #ifdef MultipleStripCell
    G4cout
     << G4endl
     << "--------------------End of Global Run------------------------"
     << G4endl
     << "  The run was " << nofEvents << " "<< partName;
    #endif

    /*G4int nDy = 103;
    G4double scale = 1/nevents;
    for(int y = 0; y<(nDy); y++)
    {
      analysisManager->GetH1(8)->multiply(scale);
      analysisManager->GetH1(9)->multiply(scale);
      analysisManager->GetH1(10)->multiply(scale);
      analysisManager->GetH1(11)->multiply(scale);
    }*/

  }
  else
  {
#ifdef SingleStrip
      #ifndef LEGEND
        G4cout
        << G4endl
        << "--------------------End of Local Run------------------------"
        << G4endl
        << "  The run was " << (nofEvents*100) << " "<< partName;
        #else
        G4cout
        << G4endl
        << "--------------------End of Local Run------------------------"
        << G4endl
        << "  The run was " << (nofEvents) << " "<< partName;
      #endif
    #endif
    #ifdef MultipleStripCell
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------"
     << G4endl
     << "  The run was " << nofEvents << " "<< partName;
    #endif
  }      
  /*G4cout
     << "; Nb of 'good' e+ annihilations: " << fGoodEvents.GetValue()  << G4endl
     << " Total dose in patient : " << G4BestUnit(fSumDose.GetValue(),"Dose") 
     << G4endl 
     << "------------------------------------------------------------" << G4endl 
     << G4endl;*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
