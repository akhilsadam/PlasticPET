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
/// \file B3DetectorConstruction.cc
/// \brief Implementation of the B3DetectorConstruction class

#include "B3DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSDFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4Scintillation.hh"
#include "G4PhysicalConstants.hh"

#include "CADMesh.hh" // Important!

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::B3DetectorConstruction()
: DetectorConstruction(),
  fCheckOverlaps(true)
{
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::~B3DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();
  
  G4bool isotopes = false;
  
  G4Element*  H = man->FindOrBuildElement("H" , isotopes);
  G4Element*  O = man->FindOrBuildElement("O" , isotopes); 
  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);
  G4Element* Y = man->FindOrBuildElement("Y", isotopes);  
  G4Element* Ar = man->FindOrBuildElement("Ar", isotopes);
  G4Element* Bi = man->FindOrBuildElement("Bi");
  G4Element* Ge = man->FindOrBuildElement("Ge");  

  G4Material* LYSO = new G4Material("LYSO", 7.1*g/cm3, 4);
  LYSO->AddElement(Lu, 18);
  LYSO->AddElement(Y, 2);
  LYSO->AddElement(Si, 10);
  LYSO->AddElement(O , 50);  

  G4Material* LHYD = new G4Material("lH",0.07*g/cm3, 1); //should really be 0.07g/cm^3
  LHYD->AddElement(H, 1);

  G4Material* lAr = new G4Material("liquidArgon",1.390*g/cm3, 1); 
  lAr->AddElement(Ar, 1);

  G4Material* BGO = new G4Material("BGO", 7.13 * g/cm3, 3);
  BGO->AddElement(Bi, 4);
  BGO->AddElement(Ge, 3);
  BGO->AddElement(O, 12);

  G4MaterialPropertiesTable* BGO_MPT = new G4MaterialPropertiesTable();  
  BGO->SetMaterialPropertiesTable(BGO_MPT);

 // BGO_MPT->AddProperty("FASTCOMPONENT", Penergy, fast, 2.15); WHAT IS THIS EXACTLY?
  BGO_MPT->AddConstProperty("RINDEX", 2.15);
 // BGO_MPT->AddProperty("ABSLENGTH", Penergy, 1.0*cm, 2.15);
  BGO_MPT->AddConstProperty("SCINTILLATIONYIELD", 8500 / MeV);
  BGO_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  BGO_MPT->AddConstProperty("FASTTIMECONSTANT", 300 * ns);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B3DetectorConstruction::Construct()
{  
  // Detector Parameters
  //
  G4double cryst_dX = 1*cm, cryst_dY = 1*cm, cryst_dZ = 1*cm;
  G4double p = 0.3;
  G4double ringAngle = twopi*p;//twopi;
  G4double offset = twopi*(0.5-p)/2;
  G4int nb_cryst = 32*3*(ringAngle/twopi);
  G4int nb_rings = 9*3;
  //
  G4double dPhi = ringAngle/nb_cryst, half_dPhi = 0.5*dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);
  // 
  G4double ring_R1 = 0.3*cryst_dY/tandPhi;
  G4double ring_R2 = (ring_R1+cryst_dZ)/cosdPhi;
  //
  G4double detector_dZ = nb_rings*cryst_dX;
  //
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* default_mat = nist->FindOrBuildMaterial("lH"); //MODIFIED - swapped out G4_AIR with lAr
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("lH");
        
  //     
  // World
  //
   world_sizeXY = 2.4*ring_R2;
   world_sizeZ  = 1.2*detector_dZ;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size
     //G4cout << "CONSTRUCTING WORLD : " << world_sizeXY << G4endl;
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                         default_mat,        //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps 

 /*                
  //
  // ring
  //
  G4Tubs* solidRing =
    new G4Tubs("Ring", ring_R1, ring_R2, 0.5*cryst_dX, offset, ringAngle);
      
  G4LogicalVolume* logicRing =                         
    new G4LogicalVolume(solidRing,           //its solid
                        default_mat,         //its material
                        "Ring");             //its name
  G4Tubs* solidRing2 =
    new G4Tubs("Ring", ring_R1, ring_R2, 0.5*cryst_dX, 3*offset+ringAngle, ringAngle);
      
  G4LogicalVolume* logicRing2 =                         
    new G4LogicalVolume(solidRing2,           //its solid
                        default_mat,         //its material
                        "Ring");             //its name
                    
  //     
  // define crystal
  //
  G4double gap = 0.5*mm;        //a gap for wrapping
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap;
  G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);
                     
  G4LogicalVolume* logicCryst = 
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name
               
  // place crystals within rings 
  //
  for (G4int icrys = 0; icrys < nb_cryst ; icrys++) {
    G4double phi = icrys*dPhi + offset;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg); 
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
    G4ThreeVector position = (ring_R1+0.5*cryst_dZ)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);
                                    
    new G4PVPlacement(transform,             //rotation,position
                      logicCryst,            //its logical volume
                      "crystal",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }

 for (G4int icrys = 0; icrys < nb_cryst ; icrys++) {
    G4double phi = icrys*dPhi + 3*offset + ringAngle;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg); 
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
    G4ThreeVector position = (ring_R1+0.5*cryst_dZ)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);
                                    
    new G4PVPlacement(transform,             //rotation,position
                      logicCryst,            //its logical volume
                      "crystal",             //its name
                      logicRing2,            //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }
                                                      
  //
  // full detector
  //
  G4Tubs* solidDetector =
    new G4Tubs("Detector", ring_R1, ring_R2, 0.5*detector_dZ, offset, ringAngle);
  G4Tubs* solidDetector2 =
    new G4Tubs("Detector", ring_R1, ring_R2, 0.5*detector_dZ, 3*offset+ringAngle, ringAngle);

  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(solidDetector,       //its solid
                        default_mat,         //its material
                        "Detector");         //its name
  G4LogicalVolume* logicDetector2 =                         
    new G4LogicalVolume(solidDetector2,       //its solid
                        default_mat,         //its material
                        "Detector");         //its name
                                 
  // 
  // place rings within each detector 
  //
  G4double OG = -0.5*(detector_dZ + cryst_dX);
  for (G4int iring = 0; iring < nb_rings ; iring++) {
    OG += cryst_dX;
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,OG), //position
                      logicRing,             //its logical volume
                      "ring",                //its name
                      logicDetector,         //its mother  volume
                      false,                 //no boolean operation
                      iring,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }

  OG = -0.5*(detector_dZ + cryst_dX);
  for (G4int iring = 0; iring < nb_rings ; iring++) {
    OG += cryst_dX;
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,OG), //position
                      logicRing2,            //its logical volume
                      "ring",                //its name
                      logicDetector2,        //its mother  volume
                      false,                 //no boolean operation
                      iring,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }
                       
  //
  // place detector in world
  //                    
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicDetector,           //its logical volume
                    "Detector",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps 
                 
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicDetector2,           //its logical volume
                    "Detector",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps 
   */              
  //
  // patient
  //
/*
  auto mesh = CADMesh::TessellatedMesh::FromOBJ("recompiled.obj");//"EMPHATIC_Setup_1.stl"
  mesh->SetScale(1);
  //mesh->SetOffset(G4ThreeVector(0,40,0));

  G4Material* mat = nist->FindOrBuildMaterial("G4_B-100_BONE"); // does not need to be set here!

  G4LogicalVolume* lv_i;

for (G4int i = 0; i < 109; i++) {

  char a[100];
  sprintf(a,"%d", i);
  strcat(a,"_EmphaticLV");
  G4cout << "Volumes " << a <<G4endl;
    
    auto solid = mesh->GetSolid(i);
    lv_i = new G4LogicalVolume(solid,        //its solid
                                mat,         //its material // you could change material here
                                a);        //its name (needs to be unique!)
                                    
    new G4PVPlacement(0,                       //no rotation
                      G4ThreeVector(0,40,0), 
                      lv_i,            //its logical volume
                      a,             //its name
                      logicWorld,            //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps 

 // G4VisAttributes* lbVis = new G4VisAttributes(G4Colour(1.0*(G4UniformRand()-0.5),1.0*(G4UniformRand()-0.5),1.0*(G4UniformRand()-0.5))); //random color!
  //lv_i->SetVisAttributes (lbVis);
  }*/


/*auto solid = mesh->GetSolid();
    lv_i = new G4LogicalVolume(solid,        //its solid
                                mat,         //its material // you could change material here
                                "_EmphaticLV");        //its name (needs to be unique!)
                                    
 // rotation

    G4RotationMatrix rotm2  = G4RotationMatrix();
    rotm2.rotateX(-90*deg); 
    G4ThreeVector* position2 = new G4ThreeVector(0,40,0);
    G4Transform3D transform2 = G4Transform3D(rotm2,(*position2));


    new G4PVPlacement(transform2,                     //no rotation 
                      lv_i,            //its logical volume
                      "EMPHATIC",             //its name
                      logicWorld,            //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps 


*/

  //
  // place patient in world
  //                    
  /*new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicPatient,            //its logical volume
                    "Patient",               //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps */
                                          
  // Visualization attributes
  //

  //logicRing->SetVisAttributes (G4VisAttributes::GetInvisible());
  //logicDetector->SetVisAttributes (G4VisAttributes::GetInvisible());    

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl; 

  //always return the physical World
  //
  return physWorld;
}

void B3DetectorConstruction::ConstructSDandField()
{
 /* G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  
  // declare crystal as a MultiFunctionalDetector scorer
  //  
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4VPrimitiveScorer* p1 = new G4PSEnergyDeposit("edep");

  G4VPrimitiveScorer* pPr = new G4PSEnergyDeposit("edep_p");
   G4VSDFilter* prFilter = new G4SDParticleFilter("prFilter","proton");
   pPr->SetFilter(prFilter);

  G4VPrimitiveScorer* pPo = new G4PSEnergyDeposit("edep_ep");
   G4VSDFilter* poFilter = new G4SDParticleFilter("poFilter","e+");
   pPo->SetFilter(poFilter); 

  G4VPrimitiveScorer* pEl = new G4PSEnergyDeposit("edep_en");
   G4VSDFilter* elFilter = new G4SDParticleFilter("elFilter","e-");
   pEl->SetFilter(elFilter); 

  G4VPrimitiveScorer* pY = new G4PSEnergyDeposit("edep_y");
   G4VSDFilter* yFilter = new G4SDParticleFilter("yFilter","gamma");
   pY->SetFilter(yFilter); 
  
  cryst->RegisterPrimitive(p1);
  cryst->RegisterPrimitive(pPr);
  cryst->RegisterPrimitive(pPo);
  cryst->RegisterPrimitive(pEl);
  cryst->RegisterPrimitive(pY);
  SetSensitiveDetector("CrystalLV",cryst);*/
  
  // declare patient as a MultiFunctionalDetector scorer
  //  
  G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("patient");
  G4SDManager::GetSDMpointer()->AddNewDetector(patient);
  G4VPrimitiveScorer* sec1 = new G4PSDoseDeposit("dose");
  patient->RegisterPrimitive(sec1);

 /*for (G4int i = 0; i < 109; i++) {
  char a[100];
  sprintf(a,"%d", i);
  strcat(a,"_EmphaticLV");
  SetSensitiveDetector(a,patient);
}*/

 SetSensitiveDetector("world_volume",patient);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
