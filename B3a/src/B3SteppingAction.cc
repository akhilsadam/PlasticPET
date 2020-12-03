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
/// \file electromagnetic/TestEm4/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "Version.hh"
#include "B3SteppingAction.hh"
#include "B3aEventAction.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
//#include "B3DetectorConstruction.hh"
#include "GDMLDetectorConstruction.hh"
#include "B3aHistoManager.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include  "G4OpBoundaryProcess.hh"

std::mutex foo2;
std::mutex barL2;
G4int B3SteppingAction::id = 0;
G4double lastEnergy=0;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3SteppingAction::B3SteppingAction(B3aEventAction* EvAct, G4VUserDetectorConstruction* patient)
:G4UserSteppingAction(),fEventAction(EvAct),fpatient(patient)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3SteppingAction::~B3SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3SteppingAction::UserSteppingAction(const G4Step* step)
{
 	std::lock(foo2,barL2);
 	G4double edep = (step->GetTotalEnergyDeposit())/keV;
	G4Track* track = step->GetTrack();
 	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();     

 	//longitudinal profile of deposited energy
 	//randomize point of energy deposotion
	 //

 	G4StepPoint* prePoint  = step->GetPreStepPoint();
 	G4StepPoint* postPoint = step->GetPostStepPoint(); 
 	G4ThreeVector P1 = prePoint ->GetPosition();
 	G4ThreeVector P2 = postPoint->GetPosition();
 	G4ThreeVector point = P1 + G4UniformRand()*(P2 - P1);
	if (track->GetDefinition()->GetPDGCharge() == 0.) point = P2;
	G4double x = P2.x();
 	G4double y = P2.y();
 	G4double z = P2.z();
 	G4double xr = point.x();
 	G4double yr = point.y();
 	G4double zr = point.z();
 	G4double xshifted = xr + 0.5*(((DetectorConstruction*) fpatient)->GetWorldSizeXY());
 	G4double yshifted = yr + 0.5*(((DetectorConstruction*) fpatient)->GetWorldSizeXY());
 	G4double zshifted = zr + 0.5*(((DetectorConstruction*) fpatient)->GetWorldSizeZ());

 	analysisManager->FillH1(0, yshifted, edep);//y,edep
 	analysisManager->FillH2(0, xshifted, yshifted, edep);
 	analysisManager->FillH3(0, xshifted, yshifted, zshifted, edep);

	G4double Ox = -2.58*cm;
	G4double Oy = 4.83*cm;
	G4double Dx = 7.74*cm;
	G4double Dy = 10.32*cm;
	G4double Dz = 100*cm;
	G4double Oz = Dz/2;

	if(track->GetParticleDefinition()->GetParticleName() == "opticalphoton")
	{
		if(track->GetTrackStatus() == fStopAndKill)
		{
			G4double x  = track->GetPosition().x();
       		G4double y  = track->GetPosition().y();
       		G4double z  = track->GetPosition().z();
			analysisManager->FillH3(1,x-Ox,y-Oy,z-Oz,1);
			if(prePoint->GetStepStatus() == fGeomBoundary)
			{
				analysisManager->FillH1(2,1);
				const G4VProcess* pds = postPoint->GetProcessDefinedStep();
				if(pds)
				{
					//cout << "PDS: " << pds->GetProcessName() << G4endl;
					const std::string pdsN = pds->GetProcessName();
					if(pdsN.compare("Transportation")==0)
					{
						analysisManager->FillH1(3,0.5);
					}
					else if(pdsN.compare("OpAbsorption")==0)
					{
						analysisManager->FillH1(3,0.25);
					}
					else
					{
						analysisManager->FillH1(3,0.75);
					}
				}
			}
			else
			{
				analysisManager->FillH1(2,0);
			}
		}
	}
	
 	const std::vector< const G4Track* >* secondaries = step->GetSecondaryInCurrentStep();  
	 const int size = secondaries->size();
 	//std::vector< G4double> times = std::vector<G4double>();
 	//std::vector< const G4DynamicParticle*>* secondparticles = new std::vector< const G4DynamicParticle* >(size);

 	//int nOfe = 0;
 	const std::string title[] = { "e+","e-","opticalphoton","gamma","proton","alpha","Li6","Be7","C11","C12","N15","O15","O16"};
 	const int titleSize = sizeof(title)/sizeof(title[0]);

	G4double Esec = (edep);
	const G4VProcess* ps;

	for(int n = 0; n < size; n++)
	{
		G4Track* ptrG4Track = (G4Track*)(*secondaries)[n];
		const G4DynamicParticle* pp = (const G4DynamicParticle*) ptrG4Track->GetDynamicParticle();
		const G4ParticleDefinition* pd = (const G4ParticleDefinition*) pp->GetParticleDefinition();
		const std::string pdVar = pd->GetParticleName();
	
		//G4cout << "-- Secondary: " << pdVar << " || Es: " << (pp->GetKineticEnergy()/MeV) << G4endl;
		//-----------Histograms----------
		if(pdVar.compare("opticalphoton")==0)
		{
			//G4cout << "Filled Photon Deposition" << G4endl;
			analysisManager->FillH1(8, (x-Ox), 1);
			analysisManager->FillH2(1, (x-Ox),(y-Oy), 1);
			analysisManager->FillH2(17, (x-Ox),(y-Oy), 1);
			analysisManager->FillH1(9, 0.5, 1);
			ps = postPoint->GetProcessDefinedStep();
			if(ps)
			{		
				const std::string psN = ps->GetProcessName();
				if(psN.compare("eIoni")==0)
				{
					//Electron Ionization (Photoelectric)
					analysisManager->FillH2(7, (x-Ox),(y-Oy), 1);
					
				}
				else if(psN.compare("msc")==0)
				{
					//COMPTON
					analysisManager->FillH2(8, (x-Ox),(y-Oy), 1);
					//analysisManager->FillH1(10, ((Dy/2) - (y-Oy)), 1);
				}
				else if(psN.compare("Cerenkov")==0)
				{
					//Cerenkov
					analysisManager->FillH2(9, (x-Ox),(y-Oy), 1);
					//analysisManager->FillH1(11, ((Dy/2) - (y-Oy)), 1);
				}
				else if(psN.compare("eBrem")==0)
				{
					//Electron Braking radiation
					analysisManager->FillH2(10, (x-Ox),(y-Oy), 1);
					//analysisManager->FillH1(12, ((Dy/2) - (y-Oy)), 1);
				}
				else if(psN.compare("Transportation")==0)
				{
					//motion - do nothing - un-double-count
					analysisManager->FillH1(8, (x-Ox), -1);
					analysisManager->FillH2(1, (x-Ox),(y-Oy), -1);
					analysisManager->FillH2(17, (x-Ox),(y-Oy), -1);
				}
				else
				{
					G4cout << psN << " not accounted for" << G4endl;
				}

			}
		}
		//------------------------------


		Esec += (pp->GetKineticEnergy()/keV);
	
		bool filled = false;
		for(int i = 0; i < titleSize; i++)
		{
			//G4cout << pdVar << " " << title[i] << " " << titleSize << G4endl;
			if(pdVar.compare(title[i])==0)
			{
				filled = true;
				//add to histogram!!
				analysisManager->FillH1(1, (0.5+i), 1);
				break;
	
			}
		}
			
		if(!filled)
		{
		G4cout
			<< G4endl
			<< "------//\\---------- broke if statement "<< pdVar <<" -///////////\\\\\\\\\\"
			<< G4endl;
		}

	}


	std::string prim = step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
	G4double Eprim = (step->GetTrack()->GetTotalEnergy()/keV); 
	if(prim.compare("gamma")==0)
	{
		/*if(Eprim <=0)
		{
			Eprim = 0;
			G4cout<<"//\\prim <=0"<<G4endl;
		}*/
		//G4cout << "/|\\--- Total (should equal the previous primary): " << (Etot) <<G4endl;
		analysisManager->FillH1(7, (id), (Eprim+Esec-lastEnergy));
		lastEnergy = Eprim;
		id += 1;
		if(size>0)
		{
			ps = postPoint->GetProcessDefinedStep();
			if(ps)
			{
				G4String psNm = ps->GetProcessName();
				//G4cout<<(ps->GetProcessName())<<G4endl;
				analysisManager->FillH2(14, (x-Ox),(y-Oy), 1);
				if(psNm.compare("phot")==0)
				{
					analysisManager->FillH2(15, (x-Ox),(y-Oy), 1);
					
				}
				else if(psNm.compare("compt")==0)
				{
					//COMPTON
					analysisManager->FillH2(16, (x-Ox),(y-Oy), 1);
				}
				else
				{
					G4cout << psNm << " not accounted for" << G4endl;
				}

			}
		}
	}

	 //example of saving random number seed of this event, under condition
	 //// if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();  

	//-----------Histograms----------
	G4String vol = prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
	if (prePoint->GetStepStatus() == fGeomBoundary)
	{
		if(vol.compare("detVOLL")==0)
		{	
			G4double lambdaP = (h*c)/(Eprim*1000*nanop);
			//G4cout << "Filled Photon Left End Counts" << G4endl;
			analysisManager->FillH2(2, (x-Ox), (y-Oy), 1);
			analysisManager->FillH2(4, (x-Ox), (y-Oy), 1);
			analysisManager->FillH1(17,  lambdaP, 1);
		}
		else if (vol.compare("detVOLR")==0)
		{
			G4double lambdaP = (h*c)/(Eprim*1000*nanop);
			//G4cout << "Filled Photon Right End Counts" << G4endl;
			analysisManager->FillH2(3, (x-Ox), (y-Oy), 1);
			analysisManager->FillH2(5, (x-Ox), (y-Oy), 1);
			analysisManager->FillH1(18,  lambdaP, 1);
		}
	}

	/*
	// Retrieve the status of the photon (From Kyle and Chris's code)
  	G4OpBoundaryProcessStatus theStatus = Undefined;
  	G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
	if (OpManager) 
	{
    	G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
    	G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);
		for ( G4int i=0; i<MAXofPostStepLoops; i++) 
		{
			G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
			G4OpBoundaryProcess* fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
			if (fOpProcess)
			{
				theStatus = fOpProcess->GetStatus();
				//if(theStatus == Absorption)
				//G4cout << "------------------------// status: " << statusEnum[theStatus+1] << G4endl;
				break;
			}
    	}
		//G4cout << "NEW STEP" << "- ID : " << step->GetTrack()->GetTrackID() << G4endl;
  	}*/


	
	foo2.unlock();
	barL2.unlock();
}
	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
		
