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

    int leftCount[Nx][Ny]={}; //avoid using ROOT as an intermediary
    int rightCount[Nx][Ny]={};
    int stripCount[Nx][Ny]={};

    void initializeCount()          
    {
      //reset to 0
      memset(leftCount, 0, sizeof(leftCount)); 
      memset(rightCount, 0, sizeof(rightCount)); 
      memset(stripCount, 0, sizeof(stripCount)); 
    };
    void fillL(G4double x,G4double y)
    {
      int xi = int (((x/length_X) + 0.5)*Nx);
      int yi = int (((y/length_Y) + 0.5)*Ny);
      if(xi==Nx){xi=Nx-1;} //right edge !
      if(yi==Ny){yi=Ny-1;} //right edge !
      //cout << "xi,x/lx:"<< xi <<"|" << x/length_X << endl;
      leftCount[xi][yi]=leftCount[xi][yi] + 1;
    };
    void fillR(G4double x,G4double y)
    {
      int xi = int (((x/length_X) + 0.5)*Nx);
      int yi = int (((y/length_Y) + 0.5)*Ny);
      if(xi==Nx){xi=Nx-1;} //right edge !
      if(yi==Ny){yi=Ny-1;} //right edge !
      rightCount[xi][yi]=rightCount[xi][yi] + 1;
    };
    void fillS(G4double x,G4double y,G4int val)
    {
      int xi = int (((x/length_X) + 0.5)*Nx);
      int yi = int (((y/length_Y) + 0.5)*Ny);
      if(xi==Nx){xi=Nx-1;} //right edge !
      if(yi==Ny){yi=Ny-1;} //right edge !
      stripCount[xi][yi]=stripCount[xi][yi] + val;
    };


    vector<G4ThreeVector> interactionPosPhot;
    vector<G4ThreeVector> interactionPosCompt;
    vector<vector<double>> interactionPos;
    vector<vector<double>> interactionPosTrack;
    map<G4int,G4int> parentTrack; // key: currentID, value: parentID
    map<G4int,G4ThreeVector> vertexPosition; // key: currentID, value: vertexPosition
    vector<G4int> photonIDList; // list of all currentIDs
    vector<G4int> detPhotonIDList; // list of all currentIDs detected
    G4int particleIDnum;
    G4int gammaID;
    G4double threshold = 0.00001;
    vector<vector<double>> photonSiPMData; //(X,Y,Z,T,L)
    vector<vector<double>> photonReflectData; //(X,Y,Z,T,alive/dead, id, incident angle, reflected angle)
  private:
    B3aRunAction*  fRunAction;
    B3PrimaryGeneratorAction* fpga;
    G4int fCollID_cryst;
    G4int fCollID_patient;   
    G4double fTotalEnergyDeposit;   // Energy deposited 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
