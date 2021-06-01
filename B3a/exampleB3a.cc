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
/// \file exampleB3a.cc
/// \brief Main program of the B3a example
#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "Version.hh"
#include "G4UImanager.hh"
#include "G4UIGAG.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4TScoreNtupleWriter.hh"

#include "Randomize.hh"


#include "B3SteppingAction.hh"
//#include "B3DetectorConstruction.hh"
#include "GDMLDetectorConstruction.hh"
#include "B3PhysicsList.hh"
#include "B3aActionInitialization.hh"
#include "B3Analysis.hh"
#include "B3aHistoManager.hh"

#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
//#include "G4OpWLS.hh"


#include "G4GDMLParser.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  G4GDMLParser parser;

  // Detect interactive mode (if no arguments) and define UI session

  #ifdef VIEWPORT_ONLY
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
  #else
  G4UIsession* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIGAG();
  }
  #endif

  // Optionally: choose a different Random engine...
  //
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(6);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //

  parser.Read("gdml.gdml",false);//"G4EMPHATIC_out.gdml"

  //B3DetectorConstruction* det = new B3DetectorConstruction;
  //B3DetectorConstruction detr = (B3DetectorConstruction) (*(new GDMLDetectorConstruction(parser.GetWorldVolume())));
  //B3DetectorConstruction* det = &detr;
  G4VPhysicalVolume* detVol0 = parser.GetWorldVolume();
  G4VUserDetectorConstruction* det = new GDMLDetectorConstruction(parser, detVol0);
  
  runManager->SetUserInitialization(det);
  //
  //runManager->SetUserInitialization(new B3PhysicsList);

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  #ifdef ScintillationDisable
  opticalPhysics->SetScintillationYieldFactor(0.);
  #else
  opticalPhysics->SetScintillationYieldFactor(1.0);
  #endif
  opticalPhysics->SetScintillationExcitationRatio(0.);

  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
  opticalPhysics->SetScintillationByParticleType(false);
  
  opticalPhysics->SetMaxNumPhotonsPerStep(100);
  opticalPhysics->SetMaxBetaChangePerStep(10.0);

  #ifdef LEGEND
    opticalPhysics->SetWLSTimeProfile("exponential"); // not sure if this should be exponential or delta - need to verify!
  #endif


  physicsList->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physicsList);


  // Set user action initialization
  //
  runManager->SetUserInitialization(new B3aActionInitialization(det));  

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");




  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Activate score ntuple writer
  // The Root output type (Root) is selected in B3Analysis.hh.
  // The verbose level can be also set via UI commands
  // /score/ntuple/writerVerbose level
  G4TScoreNtupleWriter<G4AnalysisManager> scoreNtupleWriter;
  scoreNtupleWriter.SetVerboseLevel(0);

  // Process macro or start UI session
  //
  if ( ! ui ) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    #ifndef SSRefractionTest
      UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 20 20");
    #else
      UImanager->ApplyCommand("/vis/viewer/zoom 1024");
      UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 90 0");
    #endif
    #ifdef SingleStrip
      UImanager->ApplyCommand("/vis/filtering/trajectories/create/particleFilter");
      UImanager->ApplyCommand("/vis/filtering/trajectories/particleFilter-0/add opticalphoton");
      //UImanager->ApplyCommand("/cuts/setLowEdge 0.00001 eV");
      //UImanager->ApplyCommand("/run/initialize"); 
    #endif
    ui->SessionStart();
    delete ui;
  }


  //G4VPhysicalVolume* pWorld = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
  //parser.Write("G4gdml_out.gdml", pWorld);


  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;  //workaround for segfault here!
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
