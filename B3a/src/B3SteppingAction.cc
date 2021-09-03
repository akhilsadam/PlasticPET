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

#include "../cnpy-master/cnpy.h"
//#include "../cnpy-master/cnpy.cpp"
#include <cstdlib>
#include <map>

// std::mutex foo2;
// std::mutex barL2;
// G4int B3SteppingAction::id = 0;
// G4double lastEnergy=0;

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
 	// std::lock(foo2,barL2);

	G4double edep = (step->GetTotalEnergyDeposit())/keV;
	G4Track* track = step->GetTrack();
	std::string prim = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();

	if(prim.compare("gamma")==0)
	{
		bool newGamma = true;
		vector<G4int> gammaIDs = fEventAction->gammaID;
		G4int currentTrackID = track->GetTrackID();
		for(int gi = 0; gi<gammaIDs.size();gi++)
		{
			if (gammaIDs[gi] == currentTrackID)
			{
				newGamma = false;
				break;
			}
		}
		if (newGamma)
		{
			currentTrackID = fEventAction->gammaIDCounter;
			track->SetTrackID(fEventAction->gammaIDCounter);
			 fEventAction->gammaIDCounter =  fEventAction->gammaIDCounter + 1;
			fEventAction->gammaID.push_back(currentTrackID);
		}
	}

	if(fEventAction->parentTrack.find(track->GetTrackID())==fEventAction->parentTrack.end())
	{
		fEventAction->parentTrack[track->GetTrackID()] = track->GetParentID();
		fEventAction->vertexPosition[track->GetTrackID()] = track->GetVertexPosition();
	}
 	// G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();     

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


	// Get Array # (G4 copy)
	G4VPhysicalVolume* currentVolume = prePoint->GetTouchableHandle()->GetVolume();
	G4int copyNumber = currentVolume->GetCopyNo();
	G4int detA = 0;
	if(copyNumber!=0)
	{
		detA = (copyNumber-copyNumberOffset)/(NAssemblyElements+1);
	}
	//cout << currentVolume->GetName() << endl;
	//cout << detA <<endl;
	G4ThreeVector xyzS = fEventAction->GlobalToArrayXY(x,y,z,detA);

 	// analysisManager->FillH1(0, yshifted, edep);//y,edep
 	// analysisManager->FillH2(0, xshifted, yshifted, edep);
 	// analysisManager->FillH3(0, xshifted, yshifted, zshifted, edep);
	 
	G4double Eprim = (track->GetTotalEnergy()/keV); 

	// if(track->GetParticleDefinition()->GetParticleName() == "opticalphoton")
	// {
	// 	if(track->GetTrackStatus() == fStopAndKill)
	// 	{
	// 		G4double x  = track->GetPosition().x();
    //    		G4double y  = track->GetPosition().y();
    //    		G4double z  = track->GetPosition().z();
	// 		// analysisManager->FillH3(1,xyzS.x(),xyzS.y(),xyzS.z(),1);
	// 		if(prePoint->GetStepStatus() == fGeomBoundary)
	// 		{
	// 			// analysisManager->FillH1(2,1);
	// 			const G4VProcess* pds = postPoint->GetProcessDefinedStep();
	// 			if(pds)
	// 			{
	// 				//cout << "PDS: " << pds->GetProcessName() << G4endl;
	// 				const std::string pdsN = pds->GetProcessName();
	// 				if(pdsN.compare("Transportation")==0)
	// 				{
	// 					analysisManager->FillH1(3,0.5);
	// 				}
	// 				else if(pdsN.compare("OpAbsorption")==0)
	// 				{
	// 					analysisManager->FillH1(3,0.25);
	// 				}
	// 				else
	// 				{
	// 					analysisManager->FillH1(3,0.75);
	// 				}
	// 			}
	// 		}
	// 		else
	// 		{
	// 			analysisManager->FillH1(2,0);
	// 		}
	// 	}

	// }
	
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

		
		ptrG4Track->SetTrackID(fEventAction->particleIDnum);
		fEventAction->particleIDnum = fEventAction->particleIDnum + 1;

		#ifdef ElectronPathLength
			if((pdVar.compare("e-")==0) && (prim.compare("gamma")==0))
			{
				fEventAction->electronStartPosition[fEventAction->particleIDnum] = P2;//(G4ThreeVector(xyzS.x(),xyzS.y(),xyzS.z()));
			}
		#endif

		if(pdVar.compare("e-")==0)
		{
			fEventAction->electronProcessName.push_back(ptrG4Track->GetCreatorProcess()->GetProcessName()); // electron process name;
    		fEventAction->electronProcess.push_back(ptrG4Track->GetTrackID()); // eventID, trackID
		}

		//G4cout << "-- Secondary: " << pdVar << " || Es: " << (pp->GetKineticEnergy()/MeV) << G4endl;
		//-----------Histograms----------
		if(pdVar.compare("opticalphoton")==0)
		{
			//G4String creatorName = ptrG4Track->GetCreatorProcess()->GetProcessName();
			//if(creatorName != "Scintillation")
			//G4cout << creatorName << G4endl;
			G4int trackid = ptrG4Track->GetTrackID();
			if(count(fEventAction->photonIDList.begin(),fEventAction->photonIDList.end(),trackid))
			{
				//G4cout << "photonIDList ERROR" << G4endl;
				fEventAction->photonIDList.push_back(ptrG4Track->GetTrackID());
			}
			else
			{
				fEventAction->photonIDList.push_back(ptrG4Track->GetTrackID());
			}


			//G4cout << "Filled Photon Deposition" << G4endl;
			// analysisManager->FillH1(8, (xyzS.x()), 1);
			// analysisManager->FillH2(1, (xyzS.x()),(xyzS.y()), 1);

			fEventAction->fillS(xyzS.x(),xyzS.y(),detA,1); // OBSOLETE(right): // analysisManager->FillH2(17, (xyzS.x()),(xyzS.y()), 1);

			// analysisManager->FillH1(9, 0.5, 1);
			ps = postPoint->GetProcessDefinedStep();
			if(ps)
			{		
				const std::string psN = ps->GetProcessName();
				// if(psN.compare("eIoni")==0)
				// {
				// 	//Electron Ionization (Photoelectric)
				// 	analysisManager->FillH2(7, (xyzS.x()),(xyzS.y()), 1);
					
				// }
				// else if(psN.compare("msc")==0)
				// {
				// 	//COMPTON
				// 	analysisManager->FillH2(8, (xyzS.x()),(xyzS.y()), 1);
				// 	//analysisManager->FillH1(10, ((Dy/2) - (y-Oy)), 1);
				// }
				// else if(psN.compare("Cerenkov")==0)
				// {
				// 	//Cerenkov
				// 	analysisManager->FillH2(9, (xyzS.x()),(xyzS.y()), 1);
				// 	//analysisManager->FillH1(11, ((Dy/2) - (y-Oy)), 1);
				// }
				// else if(psN.compare("eBrem")==0)
				// {
				// 	//Electron Braking radiation
				// 	analysisManager->FillH2(10, (xyzS.x()),(xyzS.y()), 1);
				// 	//analysisManager->FillH1(12, ((Dy/2) - (y-Oy)), 1);
				// }
				if(psN.compare("Transportation")==0) //note removed else
				{
					//motion - do nothing - un-double-count
					// analysisManager->FillH1(8, (xyzS.x()), -1);
					// analysisManager->FillH2(1, (xyzS.x()),(xyzS.y()), -1);
					fEventAction->fillS(xyzS.x(),xyzS.y(),detA,-1); // OBSOLETE(right): // analysisManager->FillH2(17, xyzS.x(),xyzS.y(), -1);
				}
				else
				{
					#ifndef TACC
						G4cout << psN << " not accounted for" << G4endl;
					#endif
				}

			}
		}
		//------------------------------


		Esec += (pp->GetKineticEnergy()/keV);
	
		// bool filled = false;
		// for(int i = 0; i < titleSize; i++)
		// {
		// 	//G4cout << pdVar << " " << title[i] << " " << titleSize << G4endl;
		// 	if(pdVar.compare(title[i])==0)
		// 	{
		// 		filled = true;
		// 		//add to histogram!!
		// 		analysisManager->FillH1(1, (0.5+i), 1);
		// 		break;
	
		// 	}
		// }
			
		// if(!filled)
		// {
		// G4cout
		// 	<< G4endl
		// 	<< "------//\\---------- broke if statement "<< pdVar <<" -///////////\\\\\\\\\\"
		// 	<< G4endl;
		// }

	}

	if(prim.compare("gamma")==0)
	{
		G4int currentTrackID = track->GetTrackID();


		/*if(Eprim <=0)

		{
			Eprim = 0;
			G4cout<<"//\\prim <=0"<<G4endl;
		}*/
		//G4cout << "/|\\--- Total (should equal the previous primary): " << (Etot) <<G4endl;
		// analysisManager->FillH1(7, (id), (Eprim+Esec-lastEnergy));
		// lastEnergy = Eprim;
		// id += 1;
		if(size>0)
		{
			ps = postPoint->GetProcessDefinedStep();
			if((ps) && (step->GetStepLength()>0))
			{
				G4String psNm = ps->GetProcessName();
				
				//G4cout<<(ps->GetProcessName())<<G4endl;
				// analysisManager->FillH2(14, (xyzS.x()),(xyzS.y()), 1);
				fEventAction->interactionPos.push_back({x,y,z,0,postPoint->GetGlobalTime()/ns,track->GetTrackID()});//postPoint->GetLocalTime()
				G4int gammaProcessID = 2;
				if(psNm.compare("phot")==0)
				{
					// analysisManager->FillH2(15, xyzS.x(),xyzS.y(), 1);
					fEventAction->interactionPosPhot.push_back(G4ThreeVector(x,y,z));
					gammaProcessID = 1;
				}
				else if(psNm.compare("compt")==0)
				{
					//COMPTON
					// analysisManager->FillH2(16, xyzS.x(),xyzS.y(), 1);
					fEventAction->interactionPosCompt.push_back(G4ThreeVector(x,y,z));
					gammaProcessID = 0;
				}
				else
				{
					#ifndef TACC
						G4cout << psNm << " not accounted for" << G4endl;
					#endif
				}
				fEventAction->gammaProcessIDList.push_back(gammaProcessID);

			}
		}
	}

	 //example of saving random number seed of this event, under condition
	 //// if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();  

	//-----------Histograms----------
	else if((prim.compare("opticalphoton")==0))
	{
		G4int pid = track->GetTrackID();
		G4String vol = currentVolume->GetLogicalVolume()->GetName();

		if ((vol.find(detectorLeftName) != std::string::npos) || (vol.find(detectorRightName) != std::string::npos))
		{
			track->SetTrackStatus(fStopAndKill);
		}

		const G4VProcess* pds = postPoint->GetProcessDefinedStep();
		if (pds->GetProcessName() == "OpAbsorption") {
			if (postPoint->GetStepStatus() != fGeomBoundary) 
			{
				fEventAction->VolAbsorption = fEventAction->VolAbsorption + 1;
			}
			else
			{
				fEventAction->BoundAbsorption = fEventAction->BoundAbsorption + 1;
			}
		} 
		else if (pds->GetProcessName() == "OpRayleigh") 
		{
			fEventAction->Rayleigh = fEventAction->Rayleigh + 1;
		}


		if ((postPoint->GetStepStatus() == fGeomBoundary) && (step->GetStepLength()>0))
		{
			int rid = 2;
			int reflectProcess = 4000;
			G4OpBoundaryProcessStatus theStatus = Undefined;
			G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
			if (OpManager) 
			{
				//G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
				//G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

				G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
				G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

				//G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[MAXofPostStepLoops-1]; //; last step
				//G4OpBoundaryProcess* fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
				//reflectProcess = (int) fOpProcess->GetStatus();
				for ( G4int i=0; i<MAXofPostStepLoops; i++) 
				{
					G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
					G4OpBoundaryProcess* fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
					if (fOpProcess)
					{
						int notInList = 2;
						
						/*if(fEventAction->detPhotonIDList.size()>2)
						{
							cout << pid <<" "<<fEventAction->detPhotonIDList[fEventAction->detPhotonIDList.size()-1]<< endl;
						}*/
						for(int i =fEventAction->detPhotonIDList.size()-1;i>=0;i--)
						{
							if(fEventAction->detPhotonIDList[i] == pid)
							{
								cout << "!!" << endl;
								notInList =  -1;
								break;
							}
						}
						/*if(find(fEventAction->detPhotonIDList.begin(),fEventAction->detPhotonIDList.end(),pid)== fEventAction->detPhotonIDList.end())
						{
							notInList =  2;
						}
						else
						{

						}*/
						#ifdef  ReflectionTracking
							reflectProcess = (int) fOpProcess->GetStatus();
							//theStatus = fOpProcess->GetStatus();
							//if(theStatus == Absorption)
							//G4cout << "------------------------// status: " << statusEnum[theStatus+1] << G4endl;
							//break;

							G4double angleI = 0.0;
							G4double angleR = 0.0;
							// Incident Angle:
							G4ThreeVector photonDirection = prePoint->GetMomentum() / prePoint->GetMomentum().mag();
							const G4VTouchable *touchable = prePoint->GetTouchable();
							const G4RotationMatrix *rotation = touchable->GetRotation();
							G4RotationMatrix rotation_inv = rotation->inverse();
							G4ThreeVector translation = touchable->GetTranslation();
							G4VSolid *sector = touchable->GetSolid();
							G4ThreeVector posLocal = *rotation * (P1 - translation);

							G4ThreeVector normal =  - sector->SurfaceNormal(posLocal);
							if(sector->Inside(P1) == kOutside)
							{
								normal = - normal;
							}
							G4ThreeVector photonDirectionLocal = *rotation * photonDirection;
							G4double val0 = normal.dot(photonDirectionLocal);
							angleI = (val0/abs(val0))*acos(val0);

							if(track->GetTrackStatus() == fAlive)
							{
							// if the photon has been reflected at a boundary:
								rid = 1;
								// Reflected/Refracted Angle
								G4ThreeVector photonPostDirection = postPoint->GetMomentum() / postPoint->GetMomentum().mag();
								const G4VTouchable *touchable2 = postPoint->GetTouchable();
								const G4RotationMatrix *rotation2 = touchable2->GetRotation();
								G4RotationMatrix rotation_inv2 = rotation2->inverse();
								G4ThreeVector translation2 = touchable2->GetTranslation();
								G4VSolid *sector2 = touchable2->GetSolid();
								G4ThreeVector posLocal2 = *rotation2 * (P2 - translation2);
								G4ThreeVector normal2 = normal;// - sector2->SurfaceNormal(posLocal2);
								//if(sector->Inside(P2) != sector->Inside(P1))
								//{
								//	normal2 = - normal2;
								//}
								
								G4ThreeVector photonDirectionLocal2 = *rotation * photonPostDirection;
								G4double val = normal2.dot(photonDirectionLocal2);
								angleR = (val/abs(val))*acos(val);
							}
							else if ((track->GetTrackStatus() == fStopAndKill) || (track->GetTrackStatus() == fKillTrackAndSecondaries))
							{
								/* code */
								rid = 0;
								if ((notInList>1) && ((vol.find(detectorLeftName) != std::string::npos) || (vol.find(detectorRightName) != std::string::npos)))
								{
									rid = 4;
									G4double lambdaP = (h*c)/(Eprim*1000*nanop);
									fEventAction->photonSiPMData.push_back({xyzS.x(),xyzS.y(),xyzS,z(),prePoint->GetGlobalTime()/ns,detA,lambdaP});
									fEventAction->detPhotonIDList.push_back(pid);
									fEventAction->detectedPosition[pid] = {x,y,z};
									if(vol.find(detectorRightName) == std::string::npos)
									{
										//analysisManager->FillH2(2, xyzS.x(),xyzS.y(), 1);
										//analysisManager->FillH2(4, xyzS.x(),xyzS.y(), 1);
										//analysisManager->FillH1(17,  lambdaP, 1);
										fEventAction->fillL(xyzS.x(),xyzS.y(),detA);
									}
									else
									{
										//analysisManager->FillH2(3, xyzS.x(),xyzS.y(), 1);
										//analysisManager->FillH2(5, xyzS.x(),xyzS.y(), 1);
										//analysisManager->FillH1(18,  lambdaP, 1);	
										fEventAction->fillR(xyzS.x(),xyzS.y(),detA);
									}		
								}
							}
							else if (track->GetTrackStatus() == fSuspend)
							{
								/* code */
								rid = 2;
							}
							else if (track->GetTrackStatus() == fStopButAlive)
							{
								/* code */
								rid = 2;
							}
							else if (track->GetTrackStatus() == fPostponeToNextEvent)
							{
								/* code */
								rid = 2;
							}
							else
							{
								/* code */
								rid = 2;
							}
							fEventAction->photonReflectData.push_back({xyzS.x(),xyzS.y(),xyzS.z(),prePoint->GetGlobalTime()/ns,rid,pid,angleI,angleR,reflectProcess});
							//fEventAction->photonReflectProcess.push_back(reflectProcess);
						#else
						
							if ((track->GetTrackStatus() == fStopAndKill) || (track->GetTrackStatus() == fKillTrackAndSecondaries))
							{
								if ((notInList>1) && ((vol.find(detectorLeftName) != std::string::npos) || (vol.find(detectorRightName) != std::string::npos)))
								{
									G4double lambdaP = (h*c)/(Eprim*1000*nanop);
									fEventAction->photonSiPMData.push_back({xyzS.x(),xyzS.y(),xyzS.z(),prePoint->GetGlobalTime()/ns,detA,lambdaP});
									fEventAction->detPhotonIDList.push_back(pid);
									fEventAction->detectedPosition[pid] = {x,y,z};
									if(vol.find(detectorRightName) == std::string::npos)
									{
										//analysisManager->FillH2(2, xyzS.x(),xyzS.y(), 1);
										//analysisManager->FillH2(4, xyzS.x(),xyzS.y(), 1);
										//analysisManager->FillH1(17,  lambdaP, 1);
										fEventAction->fillL(xyzS.x(),xyzS.y(),detA);
										//cout << vol << endl;
									}
									else
									{
										//analysisManager->FillH2(3, xyzS.x(),xyzS.y(), 1);
										//analysisManager->FillH2(5, xyzS.x(),xyzS.y(), 1);
										//analysisManager->FillH1(18,  lambdaP, 1);	
										fEventAction->fillR(xyzS.x(),xyzS.y(),detA);
									}		
								}	
							}
						#endif	
					}
				}
			//G4cout << "NEW STEP" << "- ID : " << step->GetTrack()->GetTrackID() << G4endl;
			}
		}
	}

	#ifdef ElectronPathLength
		else if((prim.compare("e-")==0) && ((track->GetTrackStatus() == fStopAndKill) || (track->GetTrackStatus() == fKillTrackAndSecondaries)))
		{
			G4ThreeVector startpos = fEventAction->electronStartPosition[track->GetTrackID()];
			fEventAction->electronDisplace.push_back(P1 - startpos);
			fEventAction->electronPath.push_back(track->GetTrackLength());
		}
	#endif

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


	
	// foo2.unlock();
	// barL2.unlock();
}
	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
		
