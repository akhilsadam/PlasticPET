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
/// \file B3PrimaryGeneratorAction.cc
/// \brief Implementation of the B3PrimaryGeneratorAction class
#include "B3PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#ifdef BloodRadiation
    #include "ICRP110PhantomConstruction.hh"
#endif
#include <ctime> 
#include <cstdlib>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3PrimaryGeneratorAction::B3PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
    // default particle kinematic

  #ifdef MultipleStripCell
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);
    G4ParticleDefinition* particle = particleTable->FindParticle("gamma");//FindParticle("chargedgeantino");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,length_Z/2));
    fParticleGun->SetParticleEnergy(y_Energy);  
    G4ThreeVector randDir = G4RandomDirection().unit();
    cout << "RANDOM INIT " << randDir.x() << "|" << randDir.y() << "|" << randDir.z() << endl;
    fParticleGun->SetParticleMomentumDirection(randDir);
    #ifdef CompleteScanner
      fParticleGun2  = new G4ParticleGun(n_particle);
      fParticleGun2->SetParticleDefinition(particle);
      fParticleGun2->SetParticlePosition(G4ThreeVector(0.,0.,length_Z/2));
      fParticleGun2->SetParticleEnergy(y_Energy); 
      G4ThreeVector randDir2 = -randDir;
      cout << "OTHER INIT" << randDir2.x() << "|" << randDir2.y() << "|" << randDir2.z() << endl;
      fParticleGun2->SetParticleMomentumDirection(randDir2);
    #endif
  #endif

  #ifdef SensitivityScan
  		cout << endl;
			cout << "****USING Interaction POSITION AS DEFINED BY COMMANDLINE****: SensitivityScan" << endl;
			cout << endl;
  #endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3PrimaryGeneratorAction::~B3PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B3PrimaryGeneratorAction::CSTtest(double E)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
    // fParticleGun = new G4ParticleGun(1);
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticlePolarization(G4RandomDirection().unit());
    fParticleGun->SetParticleEnergy(E);  //visible spectrum between 400-700 nm  
    //fParticleGun->GeneratePrimaryVertex(anEvent);
}
#ifdef BloodRadiation
G4ThreeVector B3PrimaryGeneratorAction::newPosition()
{
    for (int i = 0; i < loopMax; i++)
    {
        G4double x = (length_X * (G4UniformRand() - 0.5));
        G4double y = (length_Y * (G4UniformRand() - 0.5));
        G4double z = (G4UniformRand()) * length_Z;
        if (isConfined(x, y, z))
        {
            return G4ThreeVector(x, y, z);
        }
    }
    G4cout << "[ERROR] : cannot find a position to fire gammas from!" << G4endl;
    throw("");
}

bool B3PrimaryGeneratorAction::isConfined(G4int x, G4int y, G4int z, G4int matID)
{
    return ICRP110PhantomConstruction::paramPV.InMaterial(x, y, z, matID);
}
#endif
void B3PrimaryGeneratorAction::RESETtest()
{
  #ifdef MultipleStripCell
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // G4int n_particle = 1;
  // fParticleGun  = new G4ParticleGun(n_particle);
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");//FindParticle("chargedgeantino");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(y_Energy);  //SetParticleEnergy(1*eV); between 70-250 MeV   
  #ifndef RadialSource
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
  #else
    G4ThreeVector randDir = G4RandomDirection().unit();
    cout << "RANDOM " << randDir.x() << "|" << randDir.y() << "|" << randDir.z() << endl;
    fParticleGun->SetParticleMomentumDirection(randDir);
    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1,0.1,0));
  #endif
  #endif
}


void B3PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // since seeding seems not to work
  G4long seed = anEvent->GetEventID()+time(0); //anEvent->GetEventID()+1092;
  cout << "SEED: " << seed << endl;
  G4Random::setTheSeed(seed);

 /*G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
if (particle == G4ChargedGeantino::ChargedGeantino()) {
    //fluorine 
    G4int Z = 9, A = 18;
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;
    
    G4ParticleDefinition* ion
       = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
  }
*/


  // randomized position
  //
  ///G4double x0  = 0*cm, y0  = 0*cm, z0  = 0*cm;
  ///G4double dx0 = 0*cm, dy0 = 0*cm, dz0 = 0*cm;   
  //G4double x0  = 4*cm, y0  = 4*cm, z0  = 4*cm; MODIFIED TO PLACE GUN OUTSIDE CHAMBER
  //G4double dx0 = 0*cm, dy0 = 0*cm, dz0 = 0*cm; 
  //x0 += dx0*0.001*(G4UniformRand()-0.5);
  //y0 += dy0*0.001*(G4UniformRand()-0.5);
  //z0 += dz0*0.001*(G4UniformRand()-0.5);
            
  //create vertex
  //
  #ifdef MultipleStripCell
    #ifndef BloodRadiation
        G4double x0  = 0*cm;
        #ifndef ZPredictorTest
        G4double z0  = length_D/2;
        #else
        G4double z0  = (G4UniformRand())*length_Z;
        #endif
        #ifndef YPredictorTest
        G4double y0  = 0*cm;
        #else
        G4double y0  = (length_Y*(G4UniformRand()-0.5));
        #endif
        #ifdef SensitivityScan
          z0  = (G4UniformRand())*length_Z;
          x0  = (G4UniformRand())*R_0;
        #endif
    fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
    #else
    fParticleGun->SetParticlePosition(newPosition());
    #endif
  

  #ifdef RadialSource
    G4ThreeVector randDir = G4RandomDirection().unit();
    cout << "RANDOM " << randDir.x() << "|" << randDir.y() << "|" << randDir.z() << endl;
    fParticleGun->SetParticleMomentumDirection(randDir);
    #ifndef Singles
      G4ThreeVector randDir2 = -randDir;
      cout << "OTHER " << randDir2.x() << "|" << randDir2.y() << "|" << randDir2.z() << endl;
      fParticleGun2->SetParticleMomentumDirection(randDir2);
    #endif
  #else
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  #endif

  #ifndef Singles
    fParticleGun2->GeneratePrimaryVertex(anEvent);
  #endif  

  fParticleGun->GeneratePrimaryVertex(anEvent);

  #endif

  #ifdef SingleStrip
    //G4double x0  = -2.58*2*cm, y0  = 9.5*cm, z0  = 50*cm;
    G4double x0 = -5.14*cm, y0  = 9.66*cm,z0  = 0.5*cm;
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

    #ifdef LEGEND
      // fParticleGun = new G4ParticleGun(1);
      fParticleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
      fParticleGun->SetParticlePosition(G4ThreeVector(-5.14*cm,12*cm,50*cm));
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,-1,0));
      G4double L0 = 128*nm;
      fParticleGun->SetParticleEnergy(EoL/L0);
      fParticleGun->GeneratePrimaryVertex(anEvent);
      G4cout << "DEBUG LOG ::::// Event Fired " <<G4endl;
    #else

      #ifndef CST2
        G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");
      #else
        G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
      #endif
      double deltaL = Lmax - Lmin;
      //G4cout << "Emax " << Emax/eV << " Emin " << Emin/eV << G4endl;
      for(int i = 0; i<1000; i++)
      { 
        // fParticleGun = new G4ParticleGun(1);
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticlePolarization(G4RandomDirection().unit());
        #ifdef SSRefractionTest
          G4double y = tan(G4UniformRand()*pi/2);
          fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0.001*m));
          fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,y,-1).unit());
        #elif  defined(SSReflectionTest)
          #ifndef SSSpecularReflectionTest
          fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
          fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,-1,160.4838).unit()); //0,-1,160.4838 for lower wall testing on 100cm
          #else
          fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
          fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
          #endif
        #endif
        fParticleGun->SetParticleEnergy(EoL/(G4UniformRand()*deltaL + Lmin));  //visible spectrum between 400-700 nm  
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    #endif
  #endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

