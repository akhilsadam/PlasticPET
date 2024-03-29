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
#include "G4GeometryTolerance.hh"
#include "Randomize.hh"


#include "B3SteppingAction.hh"
//#include "B3DetectorConstruction.hh"
#include "GDMLDetectorConstruction.hh"
#include "B3PhysicsList.hh"
#include "B3aActionInitialization.hh"
#include "B3Analysis.hh"
//#include "B3aHistoManager.hh"

#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "G4Threading.hh"
//#include "G4OpWLS.hh"
#ifdef G4TACC
  #include "G4OpticalParameters.hh"
#endif
#ifdef HumanoidPhantom
	#ifndef ICRPModel
		#include "G4HumanPhantomConstruction.hh"
	#else
		#include "ICRP110PhantomConstruction.hh"
	#endif
#endif

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
  G4Random::setTheEngine(new CLHEP::RanecuEngine());
  G4long seed = time(NULL);
  G4Random::setTheSeed(seed);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  #ifndef TACC_CORES
	runManager->SetNumberOfThreads(8);
  #else
	runManager->SetNumberOfThreads(68);
  #endif
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  // need object!
  //G4double WorldExtent = 1000*mm;
  //G4GeometryManager::GetInstance()->SetWorldMaximumExtent(WorldExtent);
  parser.Read("gdml.gdml",false);//"G4EMPHATIC_out.gdml"

  //B3DetectorConstruction* det = new B3DetectorConstruction;
  //B3DetectorConstruction detr = (B3DetectorConstruction) (*(new GDMLDetectorConstruction(parser.GetWorldVolume())));
  //B3DetectorConstruction* det = &detr;
  G4VPhysicalVolume* detVol0 = parser.GetWorldVolume();
  cout << detVol0 << endl;
  
  //G4VUserDetectorConstruction* det = new GDMLDetectorConstruction(parser, detVol0);
  //runManager->SetUserInitialization(det);

  #ifdef HumanoidPhantom
	DetectorConstruction* det = new GDMLDetectorConstruction(parser, detVol0);
	#ifndef ICRPModel
		DetectorConstruction* det2 = new G4HumanPhantomConstruction(det);
		runManager->SetUserInitialization(det2);
		cout << "DETECTOR1 NAME : " << det->name << endl;
		cout << "DETECTOR2 NAME : " << det2->name << endl;
	#else
		DetectorConstruction* det2 = new ICRP110PhantomConstruction(det);
		runManager->SetUserInitialization(det2);
	#endif
  #else
	DetectorConstruction* det = new GDMLDetectorConstruction(parser, detVol0);
	runManager->SetUserInitialization(det);
  #endif // HumanoidPhantom

  //
  //runManager->SetUserInitialization(new B3PhysicsList);

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  #ifndef G4TACC
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
  #else
	G4OpticalParameters* optical = G4OpticalParameters::Instance();
	#ifdef ScintillationDisable
	optical->SetScintYieldFactor(0.);
	#else
	optical->SetScintYieldFactor(1.0);
	#endif
	optical->SetScintExcitationRatio(0.);

	optical->SetCerenkovTrackSecondariesFirst(true);
	optical->SetScintTrackSecondariesFirst(true);
	optical->SetScintByParticleType(false);
	
	optical->SetCerenkovMaxPhotonsPerStep(100);
	optical->SetCerenkovMaxBetaChange(10.0);
  #endif

  #ifdef LEGEND
	#ifndef G4TACC
	opticalPhysics->SetWLSTimeProfile("exponential"); // not sure if this should be exponential or delta - need to verify!
	#else
	optical->SetWLSTimeProfile("exponential"); // not sure if this should be exponential or delta - need to verify!
	#endif
  #endif


  physicsList->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physicsList);


  // Set user action initialization
  //
	#ifndef HumanoidPhantom
	  runManager->SetUserInitialization(new B3aActionInitialization(det));  
	#else
	  runManager->SetUserInitialization(new B3aActionInitialization(det2));
	#endif
  //runManager->SetUserInitialization(new B3aActionInitialization(det));

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
  // G4TScoreNtupleWriter<G4AnalysisManager> scoreNtupleWriter;
  // scoreNtupleWriter.SetVerboseLevel(0);

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
	#ifdef HumanoidPhantom
		#ifndef ICRPModel
		// Choose model : ORNL, MIRD, MIX
		// Choose Sex of Phantom : Male or Female
			#ifndef ORNL
				UImanager->ApplyCommand("/phantom/setPhantomModel MIRD");
				UImanager->ApplyCommand("/phantom/setPhantomSex Female");
				UImanager->ApplyCommand("/phantom/buildNewPhantom");
			#else
				cout << "ORNL Phantom not yet properly implemented" << endl; throw;
				UImanager->ApplyCommand("/phantom/setPhantomModel ORNLFemale");
				UImanager->ApplyCommand("/phantom/setPhantomSex Female");
				UImanager->ApplyCommand("/phantom/buildNewPhantom");
			#endif
		#else
			cout << "ICRP Phantom not yet properly implemented" << endl;
			// Choose phantom sex(male or female)
			UImanager->ApplyCommand("/phantom/setPhantomSex female");
			UImanager->ApplyCommand("/phantom/setScoreWriterSex female");
			// Choose phantom section(head, trunk or full)
			UImanager->ApplyCommand("/phantom/setPhantomSection full");
			UImanager->ApplyCommand("/phantom/setScoreWriterSection full");
			//UImanager->ApplyCommand("/vis/drawVolume worlds");
			//UImanager->ApplyCommand("/vis/ogl/set/displayListLimit 4000000");
		#endif
	#endif
	UImanager->ApplyCommand("/control/execute init_vis.mac");
	#ifdef ICRPModel
		UImanager->ApplyCommand("/vis/ogl/set/displayListLimit 4000000");
	#endif
	#ifndef SSRefractionTest
	  #ifndef VIEWPORT_ONLY
	   //UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 20 20");
	  #endif
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

  delete det;
#ifdef HumanoidPhantom
#endif

  delete visManager;
  delete runManager;  //workaround for segfault here!
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
