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
// Code developed by:
// S.Guatelli, M. Large and A. Malaroda, University of Wollongong
//
#ifndef ICRP110PhantomConstruction_H
#define ICRP110PhantomConstruction_H 1

#include "DetectorConstruction.hh"
#include "ICRP110PhantomMessenger.hh"
#include "ICRP110PhantomNestedParameterisation.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSDoseDeposit3D.hh"

#include "G4SDManager.hh"

#include "globals.hh"
#include <map>

class G4VPhysicalVolume;
class ICRP110PhantomMaterial_Female;
class ICRP110PhantomMaterial_Male;
class G4Material;

class ICRP110PhantomConstruction : public DetectorConstruction
{
  public:
     ICRP110PhantomConstruction(DetectorConstruction* detector);
    ~ICRP110PhantomConstruction();
     G4VPhysicalVolume* Construct();
     
     inline G4int GetNumberVoxelX(){return fNVoxelX;};
     inline G4int GetNumberVoxelY() {return fNVoxelY;};
     inline G4int GetNumberVoxelZ() {return fNVoxelZ;};
     
     void SetPhantomSex(G4String);
     void SetPhantomSection(G4String);

     static ICRP110PhantomNestedParameterisation paramPV;

 private:
  void ReadPhantomData(const G4String& sex, const G4String& section);
  void ReadPhantomDataFile(const G4String& sex, const G4String& fname, G4int);
  ICRP110PhantomMaterial_Female* fMaterial_Female;
  ICRP110PhantomMaterial_Male* fMaterial_Male;
  ICRP110PhantomMessenger* fMessenger;
  DetectorConstruction* detector;
 // std::vector<G4Material*> fMaterials;
  G4int fNVoxelX;
  G4int fNVoxelY; 
  G4int fNVoxelZ;
  G4double fVoxelHalfDimX;
  G4double fVoxelHalfDimY;
  G4double fVoxelHalfDimZ;
  G4double fMinX;
  G4double fMaxX;
  G4double fMinY;
  G4double fMaxY;
  G4double fMinZ;
  G4double fMaxZ;
  G4int fNoFiles;
  G4int fNVoxels;
  size_t* fMateIDs; // index of material of each voxel
  G4String fSex;
  G4String fSection;
};
#endif

