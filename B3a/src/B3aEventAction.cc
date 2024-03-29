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
// * work  make  aNy representation or  warranty, express or implied, *
// * regarding  this  software system or assume aNy liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * aNy work based  on the software)  you  agree  to acknowledge its *
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
#include "B3PrimaryGeneratorAction.hh"
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
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"

#include <iostream>
#include <fstream>
#include <istream>
//#include "../cnpy-master/cnpy.h"
#include "../cnpy-master/cnpy.cpp"
#include <cstdlib>
#include <map>

/*G4int fCollID_cryst_p;
G4int fCollID_cryst_ep;
G4int fCollID_cryst_en;
G4int fCollID_cryst_y;*/
using namespace std;

// G4int detL_npho;
// G4int detR_npho;
std::mutex foo22;
std::mutex barL22;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::B3aEventAction(B3aRunAction* runAction,B3PrimaryGeneratorAction* pgAction)
 : G4UserEventAction(), 
   fRunAction(runAction),
   fpga(pgAction),
   fCollID_cryst(-1),
/* do we need this ??
   fCollID_cryst_p(-1),
   fCollID_cryst_ep(-1),
   fCollID_cryst_en(-1),
*/
   fCollID_patient(-1)
{

    ///vector<G4ThreeVector> 
    interactionPosPhot.reserve(maxVecSize); // global
    //vector<G4ThreeVector> 
    interactionPosCompt.reserve(maxVecSize); // global
    //vector<vector<double>> 
    //interactionPos; // x,y,z,photonCount,t,gammaID // global
    //vector<vector<double>> 
    //interactionPosTrack;
    //vector<int> 
    gammaProcessIDList.reserve(maxVecSize); // compton,photo, or other (0,1,2)
    //map<G4int, G4ThreeVector> 
    //detectedPosition; // pid->global ThreeVector Position.
    //map<G4int, G4int> 
    //detPhotonType; // photon pid->gamma interaction type.
    //map<G4int, G4int> 
    //parentTrack; // key: currentID, value: parentID
    //map<G4int, G4ThreeVector> 
    //vertexPosition; // key: currentID, value: global vertexPosition
    //vector<G4int> 
    photonIDList.reserve(maxVecSize); // list of all currentIDs
    //vector<G4int> 
    detPhotonIDList.reserve(maxVecSize); // list of all currentIDs detected
    //vector<G4int> 
    gammaID.reserve(maxVecSize);
    //vector<vector<double>> 
    //photonSiPMData; //(X,Y,Z,T,A,L)
    //vector<vector<double>> 
    //photonReflectData; //(X,Y,Z,T,alive/dead, id, incident angle, reflected angle, processName)
    //map<G4int, G4ThreeVector> 
    //electronStartPosition; //(electron id, starting position)
    //vector<double> 
    electronPath.reserve(maxVecSize); // (double pathlength)
    //vector<G4ThreeVector> 
    electronDisplace.reserve(maxVecSize); // (x,y,z displacement)
    //vector<string> 
    electronProcessName.reserve(maxVecSize); // electron process name;
    //vector<int> 
    electronProcess.reserve(maxVecSize); // trackID

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::~B3aEventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::BeginOfEventAction(const G4Event* /*evt*/)
{ 
  VolAbsorption = 0;
  BoundAbsorption = 0;
  Rayleigh = 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::EndOfEventAction(const G4Event* evt)
{
    // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    G4int entry;
    G4int left = 0;
    G4int right = 0;
    G4int nevents = fRunAction->GetNevents();

    G4int totalPhot = 0;
    G4int scintiPhot = 0;
    G4int id = 25;
    G4int id2D = 18;
    G4int leftT = 0;
    G4int rightT = 0;
    //G4double firedZ = evt->GetPrimaryVertex(0)->GetZ0();
    //G4double firedX = evt->GetPrimaryVertex(0)->GetX0()+2.58*cm;
    //G4double firedY = evt->GetPrimaryVertex(0)->GetY0()-4.83*cm;
    G4double firedZ = 0;
    G4double firedX = 0;
    G4double firedY = 0;

    const int interType = 3;
    const int xvl = 3;
    const int yvl = 16;
    const int zvl = 3;
    vector<double> eventData(zvl * yvl * xvl); //shape = {z,y,x}
    vector<double> eventDataType(interType * 2 * yvl * xvl); //shape = {type,z,y,x}
    G4PrimaryVertex* pv = evt->GetPrimaryVertex(0);
    G4ThreeVector pdir = fpga->GetParticleGun()->GetParticleMomentumDirection();
    vector<double> beamData = { pv->GetX0(),pv->GetY0(),pv->GetZ0(),pdir.x(),pdir.y(),pdir.z() }; //shape = {z,y,x}

    int evtId = evt->GetEventID();

    if (interactionPosPhot.size() > 0)
    {
        firedZ = interactionPosPhot[0].z();
        firedX = interactionPosPhot[0].x();
        firedY = interactionPosPhot[0].y();
    }

    G4double delX = (length_X / Nx);
    G4double delY = (length_Y / Ny);
    G4double p0X = -(length_X / 2);
    G4double p0Y = -(length_Y / 2);
    G4double predX = 0;
    G4double predY = 0;


#ifdef ZPredictorTest
    G4cout << "FIREDZ - ZPOS: " << firedZ << G4endl;
#endif
#ifdef YPredictorTest
    G4cout << "FIREDY - YPOS: " << firedY << G4endl;
    G4cout << "FIREDX - XPOS: " << firedX << G4endl;
#endif

#ifndef NoFileWrite
    for (int a = 0; a < (Na); a++)
    {
        for (int x = 0; x < (Nx); x++)
        {
            for (int y = 0; y < (Ny); y++)
            {
                //removing ROOT dependency here
                left = B3aEventAction::leftCount[x][y][a];
                right = B3aEventAction::rightCount[x][y][a]; //detectors
                // cout << "left: " << left << endl;
                // cout << "right: "<< right  << endl;

                scintiPhot = B3aEventAction::stripCount[x][y][a]; //scintillator

                eventData[(0 * Nx * Ny) + ((Ny - 1 - y) * Nx) + x] = (left);
                eventData[(1 * Nx * Ny) + ((Ny - 1 - y) * Nx) + x] = (scintiPhot);
                eventData[(2 * Nx * Ny) + ((Ny - 1 - y) * Nx) + x] = (right);
                entry = left + right; //detectors
                //G4cout << "ENTRY: " << entry << "X " << x << "Y " << y << G4endl;
                // if ((entry > 0) && (nevents !=0)) //detectors
                // {
                //   analysisManager->FillH2(6, (x),(y), (1.0) );
                // }

                if (scintiPhot > 0) //scintillator
                {
                    // analysisManager->FillH1(id+5, (scintiPhot) );
                    // analysisManager->FillH2(id2D, (scintiPhot), ((1.0*entry)/scintiPhot) );
                    predY = predY + (scintiPhot * (p0Y + (delY * y)));
                    predX = predX + (scintiPhot * (p0X + (delX * x)));
                    totalPhot = totalPhot + scintiPhot;
                }
                // analysisManager->FillH2(12, (x),(y), (left) );
                // analysisManager->FillH2(13, (x),(y), (right) );
                // analysisManager->FillH2(11, (right),(left), (1) ); //P@1 VS P@2
                // analysisManager->FillH1(id, entry); //(entry)
                leftT += left;
                rightT += right;
                id += 6;
                id2D += 1;
            }
        }
        G4String outname = "photonCounts" + to_string(a) + ".npy";
        std::lock(foo22, barL22);
        cnpy::npy_save(outpath + (outname), &eventData[0], { zvl,yvl,xvl }, "a"); //event photon counts
        foo22.unlock();
        barL22.unlock();
        //clear eventData
        std::fill(std::begin(eventData), std::end(eventData), 0);
    }
#endif
    // G4int photons = (G4int) analysisManager->GetH1(9)->bin_entries(0);
    // if(photons>0)
    // {
    //   analysisManager->FillH1(10, 0.5, (1.0/nevents) );
    // }
    // analysisManager->FillH1(23, leftT);//left and right aggregate histograms
    // analysisManager->FillH1(24, rightT);
    std::cout << "EVT LEFT TOTALS:" << leftT << std::endl;
    std::cout << "EVT RIGHT TOTALS:" << rightT << std::endl;
    std::cout << "EVT TOTAL uniques Detections:" << B3aEventAction::detPhotonIDList.size() << std::endl;
    // //Predict POSITION (position resolution)

    //   G4double zpred = ((att_len/(16))*log(double(leftT)/double(rightT))  + (length_D/(2)));
    //   predY  = predY/double(totalPhot);
    //   predX  = predX/double(totalPhot);
    //   if(interactionPosPhot.size() > 0)
    // {


    //   analysisManager->FillH2(analysisManager->GetH2Id("(L-R) Photon Ratio as Z-Predictor (ln(LR Ratio) by m)"),firedZ/m, zpred/m);//Z-predictor histogram
    //   analysisManager->FillH1(analysisManager->GetH1Id("(L-R) Photon Ratio as Z-Predictor"), (zpred - firedZ)/m);//Z-predictor histogram
    //   analysisManager->FillH2(analysisManager->GetH2Id("(L-R) Photon Ratio as Y-Predictor (ln(LR Ratio) by m)"),firedY/cm, predY/cm);//Y-predictor histogram
    //   analysisManager->FillH1(analysisManager->GetH1Id("(L-R) Photon Ratio as Y-Predictor"), (predY - firedY)/cm);//Y-predictor histogram
    //   analysisManager->FillH2(analysisManager->GetH2Id("(L-R) Photon Ratio as X-Predictor (ln(LR Ratio) by m)"),firedX/cm, predX/cm);//Z-predictor histogram
    //   analysisManager->FillH1(analysisManager->GetH1Id("(L-R) Photon Ratio as X-Predictor"), (predX - firedX)/cm);//Z-predictor histogram

    // }

    //data dump
    //cnpy::npz_save("posRes.npz","data",&eventData[0],{zvl,yvl,xvl},"a"); //event photon counts
#ifndef NoFileWrite
    std::lock(foo22, barL22);

    cnpy::npy_save(outpath + ("beamData.npy"), &beamData[0], { 1,2,3 }, "a"); //event location and momentum direction

    std::ofstream beamInteract(outpath + ("beamInteract.txt"), std::ios_base::app);
    for (int i = 0; i < 3; i++)
    {
        for (G4ThreeVector pos : interactionPosPhot)
        {
            switch (i)
            {
            case 0:
                beamInteract << pos.x() << " ";
                break;
            case 1:
                beamInteract << pos.y() << " ";
                break;
            case 2:
                beamInteract << pos.z() << " ";
                break;
            }
        }
        beamInteract << endl;
    }
    for (int i = 0; i < 3; i++)
    {
        for (G4ThreeVector pos : interactionPosCompt)
        {
            switch (i)
            {
            case 0:
                beamInteract << pos.x() << " ";
                break;
            case 1:
                beamInteract << pos.y() << " ";
                break;
            case 2:
                beamInteract << pos.z() << " ";
                break;
            }
        }
        beamInteract << endl;
    }


    //INTERACTION pos - photon counting
    G4cout << " GAMMA IDs = " << gammaID[0] << "|" << gammaID[1] << G4endl;
    vector<G4double> diffs;
    for (int i = 0; i < photonIDList.size(); i++)
    {
        G4int photon = photonIDList[i];
        G4int pid_original = photonIDList[i];
        G4int parent = photon;
        while ((parent != gammaID[0]) && (parent != gammaID[1]))
        {
            if (parentTrack.find(parent) != parentTrack.end())
            {
                parent = parentTrack.at(parent);
                photon = parent;
            }
            else
            {
                //G4cout << "Non-gamma parented photon" << G4endl;
                break;
            }
        }
        //now parent = gammaID, so photon is first secondary-track 
        G4ThreeVector vertexpos = vertexPosition[photon];
        G4double x = vertexpos.x(); // removing the OX OY offsets
        G4double y = vertexpos.y();
        G4double z = vertexpos.z();
        //G4cout << "(" << x <<" "<<y<<" "<<z<<")";
        for (int i = 0; i < interactionPos.size(); i++)
        {
            auto pos = interactionPos[i];
            G4double diff = pow((pos[0] - x), 2) + pow((pos[1] - y), 2) + pow((pos[2] - z), 2);
            diffs.push_back(diff);
            //G4cout << diff << " ";
        }
        int min_i = min_element(diffs.begin(), diffs.end()) - diffs.begin();
        interactionPos[min_i][3] = interactionPos[min_i][3] + 1;
        diffs.clear();
        //
        // If photon in list of detected photons, detPhotonType[pid_original] = gammaProcessIDList[min_i];
        //
        for (int k = 0; k < detPhotonIDList.size(); k++)
        {
            if (detPhotonIDList[k] == pid_original)
            {
                //cout << "This photon was detected!" << endl;
                detPhotonType[pid_original] = gammaProcessIDList[min_i];
                break;
            }
        }
    }

#ifndef CompleteScanner

    for (G4int pid_k : detPhotonIDList)
    {
        //for(int tp = 0; tp<interType; tp++)
        //for(int i = 0; i<Nx;i++)
        //for(int j = 0;j<Ny; j++)
        G4ThreeVector vpos = detectedPosition[pid_k];
        G4double x = vpos.x();
        G4double y = vpos.y();
        G4double z = vpos.z();

        int xi = int(((x / length_X) + 0.5) * Nx);
        int yi = int(((y / length_Y) + 0.5) * Ny);
        int zi = 0;
        if (xi == Nx) { xi = Nx - 1; } //right edge !
        if (yi == Ny) { yi = Ny - 1; } //right edge !
        if (z > Oz) { zi = 1; }
        int tp = detPhotonType[pid_k];
        eventDataType[(tp * 2 * Nx * Ny) + (zi * Nx * Ny) + ((Ny - 1 - yi) * Nx) + xi] += 1;
    }
    std::lock(foo22, barL22);
    cnpy::npy_save(outpath + ("photonCountTypes.npy"), &eventDataType[0], { interType,2,yvl,xvl }, "a"); //event photon counts
    foo22.unlock();
    barL22.unlock();
#endif

    string psmFile = outpath + ("photonSiPMData.txt");
    std::ofstream photonSiPM(psmFile, std::ios_base::app);
    for (auto pos : photonSiPMData)
    {
        for (int u = 0; u < pos.size(); u++)
        {
            photonSiPM << pos[u] << " ";
        }
        photonSiPM << "|";
    }
    photonSiPM << endl;
    photonSiPM.close();


#ifdef  ReflectionTracking
    string reflectfile = outpath + ("photonReflectData.txt");
    string reflectfile2 = outpath + ("photonReflectCount.txt");
#ifdef DISABLEVK
    reflectfile = outpath + ("photonReflectData_DISABLE_VK.txt");
    reflectfile2 = outpath + ("photonReflectCount_DISABLE_VK.txt");
#endif
    std::ofstream photonReflectC(reflectfile2, std::ios_base::app);
    photonReflectC << photonIDList.size() << endl;
    photonReflectC.close();
    std::ofstream photonReflect(reflectfile, std::ios_base::app);
    for (auto pos : photonReflectData)
    {
        for (int u = 0; u < pos.size(); u++)
        {
            photonReflect << pos[u] << " ";
        }
        photonReflect << "|";
    }
    photonReflect << endl;
    photonReflect.close();
#endif

    for (int i = 0; i < 6; i++)
    {
        for (auto pos : interactionPos)
        {
            beamInteract << pos[i] << " "; // note that this is in global coordinates!
        }
        beamInteract << endl;
    }
    beamInteract << endl;
    beamInteract.close();

    string volFile = outpath + ("volProcess.txt");
    std::ofstream volStream(volFile, std::ios_base::app);
    volStream << VolAbsorption << " ";
    volStream << BoundAbsorption << " ";
    volStream << Rayleigh << " ";
    volStream << endl;
    volStream.close();
#ifdef ElectronPathLength
    string electFile = outpath + ("electronPath.txt");
    std::ofstream electPathStream(electFile, std::ios_base::app);
    for (int i = 0; i < electronPath.size(); i++)
    {
        electPathStream << (electronPath[i] / mm) << " ";
        G4ThreeVector electronDisplaceTV = electronDisplace[i];
        electPathStream << electronDisplaceTV.x() << " ";
        electPathStream << electronDisplaceTV.y() << " ";
        electPathStream << electronDisplaceTV.z() << " ";
        electPathStream << "|";
    }
    electPathStream << endl;
    electPathStream.close();
#endif

    string elect = outpath + ("electronProcess.txt");
    std::ofstream electStream(elect, std::ios_base::app);
    for (int i = 0; i < electronProcess.size(); i++)
    {
        electStream << electronProcessName[i] << " ";
        electStream << electronProcess[i] << " ";
        electStream << "|";
    }
    electStream << endl;
    electStream.close();



    // scatterfraction
    int gcheck{ 0 };
    for (auto gid : gammaID)
    {
        G4ThreeVector gdir = final_gamma_position[gid] - G4ThreeVector(beamData[0], beamData[1], beamData[2]);
        double angle = gdir.angle(G4ThreeVector(beamData[3], beamData[4], beamData[5]));
        //G4cout << "ANGLE : " << angle << G4endl;
        if (fabs(angle) < eps_sf)
        {
            gcheck += -1;
        }
        else if (fabs(angle - pi) < eps_sf)
        {
            gcheck += 1;
        }
        else
        {
            gcheck += 500;
            break;
        }
    }
    string sfract = outpath + ("scatterfraction.txt");
    std::ofstream sfStream(sfract, std::ios_base::app);
    if ((gcheck >= -10) && (gcheck <= 10))
    {
        sfStream << 0 << endl;
    }
    else
    {
        sfStream << 1 << endl;
    }
    sfStream.close();
      
  foo22.unlock();
  barL22.unlock();
  

  // std::lock(foo22,barL22);  

  // foo22.unlock();
  // barL22.unlock();
  #endif

  B3aEventAction::initializeCount();
  fRunAction->CountEvent();

  // analysisManager->GetH2(4)->reset();
  // analysisManager->GetH2(5)->reset();
  // analysisManager->GetH2(17)->reset();
  // analysisManager->GetH1(9)->reset();

  interactionPos.clear();
  interactionPosPhot.clear();
  interactionPosCompt.clear();
  gammaProcessIDList.clear();
  gammaIDCounter = 0;
  detectedPosition.clear();
  detPhotonType.clear();
  parentTrack.clear();
  vertexPosition.clear();
  photonIDList.clear();
  detPhotonIDList.clear();
  particleIDnum = 0;
  photonSiPMData.clear();
  photonReflectData.clear();
  #ifdef ElectronPathLength
    electronStartPosition.clear();
    electronPath.clear();
    electronDisplace.clear();
  #endif

  electronProcessName.clear();
  electronProcess.clear();

  //photonReflectProcess.clear();

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
