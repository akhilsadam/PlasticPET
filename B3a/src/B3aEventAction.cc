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
/// \file B3aEventAction.cc
/// \brief Implementation of the B3aEventAction class
#include "Version.hh"
#include "B3aEventAction.hh"
#include "B3aRunAction.hh"
#include "B3Analysis.hh"
#include "B3Hit.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4EventManager.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

/*G4int fCollID_cryst_p;
G4int fCollID_cryst_ep;
G4int fCollID_cryst_en;
G4int fCollID_cryst_y;*/

G4int detL_npho;
G4int detR_npho;
std::mutex foo22;
std::mutex barL22;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::B3aEventAction(B3aRunAction* runAction)
 : G4UserEventAction(), 
   fRunAction(runAction),
   fCollID_cryst(-1),
/* do we need this ??
   fCollID_cryst_p(-1),
   fCollID_cryst_ep(-1),
   fCollID_cryst_en(-1),
*/
   fCollID_patient(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::~B3aEventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::BeginOfEventAction(const G4Event* /*evt*/)
{ 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::EndOfEventAction(const G4Event* evt )
{
  std::lock(foo22,barL22);
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int nx = 3;
  G4int ny = 16;
  G4double Dx = 7.74*cm;
  G4double Dy = 10.32*cm;
  G4int entry;
  G4int left=0;
  G4int right=0;
  G4int nevents = fRunAction->GetNevents();

  G4int scintiPhot=0;
  G4int id = 25;
  G4int id2D = 18;
  G4int leftT = 0;
  G4int rightT = 0;

  for(int x = 0; x<(nx); x++)
  {
    for(int y = 0; y<(ny); y++)
    {
      left = (G4int) analysisManager->GetH2(4)->bin_entries((x),(y));
      right = (G4int) analysisManager->GetH2(5)->bin_entries((x),(y)); //detectors
      scintiPhot = (G4int) analysisManager->GetH2(17)->bin_entries((x),(y)); //scintillator
      entry = left + right; //detectors
      //G4cout << "ENTRY: " << entry << "X " << x << "Y " << y << G4endl;
      if ((entry > 0) && (nevents !=0)) //detectors
      {
        //G4cout << " Filling " << G4endl;
        analysisManager->FillH2(6, (x),(y), (1.0) );
      }

      if(scintiPhot>0) //scintillator
      {
        analysisManager->FillH1(id+5, (scintiPhot) );
        analysisManager->FillH2(id2D, (scintiPhot), ((1.0*entry)/scintiPhot) );
      }
      analysisManager->FillH2(12, (x),(y), (left) );
      analysisManager->FillH2(13, (x),(y), (right) );
      analysisManager->FillH2(11, (right),(left), (1) ); //P@1 VS P@2
      analysisManager->FillH1(id, entry); //(entry)
      leftT+=left;
      rightT+=right;
      id+=6;
      id2D+=1;
    }
  }
  G4int photons = (G4int) analysisManager->GetH1(9)->bin_entries(0);
  if(photons>0)
  {
    analysisManager->FillH1(10, 0.5, (1.0/nevents) );
  }
  analysisManager->FillH1(23, leftT);
  analysisManager->FillH1(24, rightT);

  analysisManager->GetH2(4)->reset();
  analysisManager->GetH2(5)->reset();
  analysisManager->GetH2(17)->reset();
  analysisManager->GetH1(9)->reset();

  foo22.unlock();
  barL22.unlock();

/*
 if ( detL_npho <=0 ) {
   detL_npho = G4SDManager::GetSDMpointer()->GetCollectionID("detL/nPho");
  }
 if ( detR_npho <=0 ) {
   detR_npho = G4SDManager::GetSDMpointer()->GetCollectionID("detR/nPho");
  }



   //Hits collections
  //  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) 
  {
    G4cout << "NO HIT COLLECTIONS!" << G4endl;
    return;
  }
    G4SDManager* SDMan = G4SDManager::GetSDMpointer();     
    SDMan->ListTree();     

*/
    //SDMan->SetVerboseLevel(0);       
   // Get hits collections IDs
  /*if (fCollID_cryst < 0) {
 
    fCollID_cryst   = SDMan->GetCollectionID("crystal/edep");

    fCollID_cryst_p    = SDMan->GetCollectionID("crystal/edep_p");
    fCollID_cryst_ep   = SDMan->GetCollectionID("crystal/edep_ep");
    fCollID_cryst_en   = SDMan->GetCollectionID("crystal/edep_en");
    fCollID_cryst_y   = SDMan->GetCollectionID("crystal/edep_y");

        
  }*/

/*
  
  //Energy in crystals : identify 'good events'
  //
  const G4double eThreshold = 500*keV;
  G4int nbOfFired = 0;
  bool mapL = true;
  bool mapR = true;
  B3Hit* hit;  
 
  B3HitsCollection* HCmapL = (B3HitsCollection*)(HCE->GetHC(detL_npho)); 
  if(!HCmapL)  
  {
    G4cout << "NO HIT COLLECTION L!" << G4endl;
    mapL = false;
  }
  uint lengthL;
  std::vector<B3Hit>*  evtMapL;
  
  B3HitsCollection* HCmapR = (B3HitsCollection*)(HCE->GetHC(detR_npho)); 
  if(!HCmapR)
  {
    G4cout << "NO HIT COLLECTION R!" << G4endl;
    mapR = false;
  }
  uint lengthR;
  std::vector<B3Hit>*  evtMapR;
*/
  // "LIT" deposit in RL detectors
  // 
  /*
  if(mapL)
  {
    try
    { 
     evtMapL =  (std::vector<B3Hit>*) (HCmapL->GetVector());
     lengthL = evtMapL->size(); 
    } 
    catch (...) { G4cerr << "EXCEPTION - DETECTORSL - HELP" << G4endl; }
    std::cout << "LENGTHL" << lengthL <<std::endl;
    for (G4int i = 0; i < lengthL; i++) 
    {
   	  hit = (B3Hit*) &(evtMapL->at(i));
    	G4ThreeVector pos = (hit->GetPos());
	    //G4int prmE = (G4int) analysisManager->GetH2(4)->bin_entries((pos.x()),(pos.y()));
	    //G4cout<< "LIT" << prmE <<G4endl;
	    analysisManager->FillH2(4, (pos.x()), (pos.y()), 1 );
      G4cout << "FilledL" << G4endl;
	    //prmE = (G4int) analysisManager->GetH2(4)->bin_entries((pos.x()),(pos.y()));
	    //G4cout<< "LIT" << prmE <<G4endl;
    }
  }
  
  if(mapR)
  {
    try
    { 
      evtMapR =  (std::vector<B3Hit>*) (HCmapR->GetVector());
      lengthR = evtMapR->size(); 
    } 
    catch (...) { G4cerr << "EXCEPTION - DETECTORSR - HELP" << G4endl; }
    std::cout << "LENGTHR" << lengthR << std::endl;
    for (G4int i = 0; i < 12; i++) 
    {
   	  hit = (B3Hit*) &(evtMapR->at(i));
    	G4ThreeVector pos = (hit->GetPos());
	    //G4int prmE = (G4int) analysisManager->GetH2(4)->bin_entries((pos.x()),(pos.y()));
	    //G4cout<< "LIT" << prmE <<G4endl;
	    //analysisManager->FillH2(4, (pos.x()), (pos.y()), 1 );
      //G4cout << "FilledR" << G4endl;
	    //prmE = (G4int) analysisManager->GetH2(4)->bin_entries((pos.x()),(pos.y()));
	    //G4cout<< "LIT" << prmE <<G4endl;
    }
  
  }







/*
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    //G4int copyNb  = (itr->first);
    G4double edep = *(itr->second);
    if (edep > eThreshold) nbOfFired++;
    //G4cout << "\n  cryst" << copyNb << ": " << edep/keV << " keV ";
    *(itr->second) = edep/MeV; // MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(2, *(itr->second));
  }  
  if (nbOfFired == 2) fRunAction->CountEvent();
 // MODIFIED ----------------------------------------------------------------------------------------
  evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst_p));
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep_p = *(itr->second);
    *(itr->second) = edep_p/MeV; // MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(3, *(itr->second));
  }
  // MODIFIED ----------------------------------------------------------------------------------------
  evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst_ep));
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep_ep = *(itr->second);
    *(itr->second) = edep_ep/MeV; // MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(4, *(itr->second));
  }  
  // MODIFIED ----------------------------------------------------------------------------------------
  evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst_en));
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep_en = *(itr->second);
    *(itr->second) = edep_en/MeV; // MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(5, *(itr->second));
  }
  // MODIFIED ----------------------------------------------------------------------------------------
  evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst_y));
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep_y = *(itr->second);
    *(itr->second) = edep_y/MeV; // MODIFIED TO CONVERT UNITS
	analysisManager->FillH1(6, *(itr->second));
  }
  // ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
  //
  //analysisManager->FillH1(1, fTotalEnergyDeposit/MeV);
*/
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
