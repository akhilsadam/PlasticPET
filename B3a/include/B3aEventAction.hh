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
/// \file B3aEventAction.hh
/// \brief Definition of the B3aEventAction class

#ifndef B3aEventAction_h
#define B3aEventAction_h 1

#ifndef VERSION
#include "Version.hh"
#endif

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>
#include <array>
#include <unordered_map>
#include "G4ThreeVector.hh"
using namespace std;

class B3aRunAction;
class B3PrimaryGeneratorAction;
/// Event action class
///
/// In EndOfEventAction() there is collected information event per event 
/// from Hits Collections, and accumulated statistic for 
/// B3RunAction::EndOfRunAction().

class B3aEventAction : public G4UserEventAction
{
  public:
    B3aEventAction(B3aRunAction* runAction, B3PrimaryGeneratorAction* pgAction);
    virtual ~B3aEventAction();

    virtual void  BeginOfEventAction(const G4Event*);
    virtual void    EndOfEventAction(const G4Event*);

    
    void AddEdep(G4double Edep)     {fTotalEnergyDeposit += Edep;};      
    G4double GetEnergyDeposit()     {return fTotalEnergyDeposit;};   

    int leftCount[Nx][Ny][Na]={}; //avoid using ROOT as an intermediary
    int rightCount[Nx][Ny][Na]={};
    int stripCount[Nx][Ny][Na]={};

    G4ThreeVector GlobalToArrayXY(G4double x, G4double y, G4double z, G4int i)
    {
      double angle = i*(2*M_PI/NArray);
      G4ThreeVector origin = G4ThreeVector(R*cos(angle),R*sin(angle),0);
      G4ThreeVector translatedPosition = G4ThreeVector(x,y,z) - origin;
      // cout << "XYZ:" << x << "|" << y << "|" << z << endl;
      // cout << "newOrigin:" << newOrigin->x() << "|" << newOrigin->y() << "|" << newOrigin->z() << endl;
      // cout << "TP:" << translatedPosition.x() << "|" << translatedPosition.y() << "|" << translatedPosition.z() << endl;
      double xy = sqrt(pow(translatedPosition.x(),2)+pow(translatedPosition.y(),2));
		  double angle2 = -angle + atan2(translatedPosition.y(),translatedPosition.x());
      return G4ThreeVector(xy*cos(angle2),xy*sin(angle2),z);
    }
    G4ThreeVector ArrayToGlobalXY(G4double x,G4double y,G4double z,G4int i)
    {
      double angle = i*(2*M_PI/NArray);
      double cosine = cos(angle);
      double sine = sin(angle);
      double x_global = (R+x)*cosine - y*sine;
      double y_global = (R+x)*sine + y*cosine;
      return G4ThreeVector(x_global,y_global,z);
    }
    void initializeCount()          
    {
      //reset to 0
      memset(leftCount, 0, sizeof(leftCount)); 
      memset(rightCount, 0, sizeof(rightCount)); 
      memset(stripCount, 0, sizeof(stripCount)); 
    };
    void fillL(G4double x,G4double y,G4int A)
    {
      int xi = int (((x/length_X) + 0.5)*Nx);
      int yi = int (((y/length_Y) + 0.5)*Ny);
      if(xi==Nx){xi=Nx-1;} //right edge !
      if(yi==Ny){yi=Ny-1;} //right edge !
      //cout << "xi,x/lx:"<< xi <<"|" << x/length_X << endl;
      if ((0<=xi) && (xi<Nx) && (0<=yi) && (yi<Ny) && (0<=A) && (A<Na))
      {
        leftCount[xi][yi][A]=leftCount[xi][yi][A] + 1;
        //cout << "L : xi: "<< xi <<"|yi: " << yi << "|A: " << A << "|in: " << x << "|" << y << "|" << A << endl;
      }
      else
      {
        //cout << "---------ERROR L : xi: "<< xi <<"|yi: " << yi << "|A: " << A << "|in: " << x << "|" << y << "|" << A << endl;
      }
    };
    void fillR(G4double x,G4double y,G4int A)
    {
      int xi = int (((x/length_X) + 0.5)*Nx);
      int yi = int (((y/length_Y) + 0.5)*Ny);
      if(xi==Nx){xi=Nx-1;} //right edge !
      if(yi==Ny){yi=Ny-1;} //right edge !

      if ((0<=xi) && (xi<Nx) && (0<=yi) && (yi<Ny) && (0<=A) && (A<Na))
      {
        rightCount[xi][yi][A]=rightCount[xi][yi][A] + 1;
      }
      else
      {
        //cout << "---------ERROR R : xi: "<< xi <<"|yi: " << yi << "|A: " << A << endl;
      }
    };
    void fillLRT(G4double x,G4double y,G4double z,G4double lrt,G4int A)
    {
      int xi = int (((x/length_X) + 0.5)*Nx);
      int yi = int (((y/length_Y) + 0.5)*Ny);
      if(xi==Nx){xi=Nx-1;} //right edge !
      if(yi==Ny){yi=Ny-1;} //right edge !
      if(z<Oz)
      {
        //leftCount[xi][yi]=leftCount[xi][yi] + 1;
      }
      else
      {
        //rightCount[xi][yi]=rightCount[xi][yi] + 1;
      }
    };
    void fillS(G4double x,G4double y,G4int A, G4int val)
    {
      int xi = int (((x/length_X) + 0.5)*Nx);
      int yi = int (((y/length_Y) + 0.5)*Ny);
      if(xi==Nx){xi=Nx-1;} //right edge !
      if(yi==Ny){yi=Ny-1;} //right edge !
      
      if ((0<=xi) && (xi<Nx) && (0<=yi) && (yi<Ny) && (0<=A) && (A<Na))
      {
        stripCount[xi][yi][A]=stripCount[xi][yi][A] + val;
        //cout << "S : xi: "<< xi <<"|yi: " << yi << "|A: " << A << "|val: " << val << endl;
      }
      else
      {
        //cout << "---------ERROR S : xi: "<< xi <<"|yi: " << yi << "|A: " << A << "|val: " << val << endl;
      }
    };


    vector<G4ThreeVector> interactionPosPhot; // global
    vector<G4ThreeVector> interactionPosCompt; // global
    vector<array<double,6>> interactionPos; // x,y,z,photonCount,t,gammaID // global
    //vector<vector<double>> interactionPosTrack; 
    vector<int> gammaProcessIDList; // compton,photo, or other (0,1,2)
    unordered_map<G4int,G4ThreeVector> detectedPosition; // pid->global ThreeVector Position.
    unordered_map<G4int,G4int> detPhotonType; // photon pid->gamma interaction type.
    unordered_map<G4int,G4int> parentTrack; // key: currentID, value: parentID
    unordered_map<G4int,G4ThreeVector> vertexPosition; // key: currentID, value: global vertexPosition
    vector<G4int> photonIDList; // list of all currentIDs
    vector<G4int> detPhotonIDList; // list of all currentIDs detected
    G4int particleIDnum;
    vector<G4int> gammaID;
    G4int gammaIDCounter = 0;
    G4double threshold = 0.00001;
    vector<array<double,6>> photonSiPMData; //(X,Y,Z,T,A,L)
    vector<array<double,9>> photonReflectData; //(X,Y,Z,T,alive/dead, id, incident angle, reflected angle, processName)
    unordered_map<G4int,G4ThreeVector> electronStartPosition; //(electron id, starting position)
    vector<double> electronPath; // (double pathlength)
    vector<G4ThreeVector> electronDisplace; // (x,y,z displacement)
    //vector<int> photonReflectProcess; // as above, but process name.
    vector<string> electronProcessName; // electron process name;
    vector<int> electronProcess; // trackID

    //scatter fraction
    unordered_map<G4int, G4ThreeVector> final_gamma_position;
    // have a 6-vector of pos,dir for gamma_shot.

    G4int VolAbsorption;
    G4int BoundAbsorption;
    G4int Rayleigh;

  private:
    B3aRunAction*  fRunAction;
    B3PrimaryGeneratorAction* fpga;
    G4int fCollID_cryst;
    G4int fCollID_patient;   
    G4double fTotalEnergyDeposit;   // Energy deposited 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
