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
/// \file B3DetectorConstruction.hh
/// \brief Definition of the B3DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#include <Version.hh>
#include "G4VUserDetectorConstruction.hh"
#include <globals.hh>


#include <G4NistManager.hh>
#include <G4RunManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4RotationMatrix.hh>
#include <G4Transform3D.hh>
#include <G4SDManager.hh>
#include <G4MultiFunctionalDetector.hh>
#include <G4VPrimitiveScorer.hh>
#include <G4PSEnergyDeposit.hh>
#include <G4PSDoseDeposit.hh>
#include <G4PSNofCollision.hh>
#include <G4VisAttributes.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4VSDFilter.hh>
#include <G4SDParticleFilter.hh>
#include <G4Scintillation.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4PhysicalConstants.hh>

#include <G4GeometryManager.hh>
#include <G4VVisManager.hh>
#ifdef LYSOTest
#include <LYSO.hh>
#endif

#include <G4AssemblyVolume.hh>
#include <G4PhysicalVolumeStore.hh>

#include <G4GDMLParser.hh>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
using namespace std;

#include "Data.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.
///
/// Crystals are positioned in Ring, with an appropriate rotation matrix. 
/// Several copies of Ring are placed in the full detector.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public: 
	inline G4double GetWorldSizeXY() { return world_sizeXY; };
	inline G4double GetWorldSizeZ() { return world_sizeZ; };
	void UpdateWorldDataWithPhantom(G4VPhysicalVolume* in_world) { delete World; World = in_world; };
	G4VPhysicalVolume* GetWorld() { return World; };
	G4double world_sizeXY;
	G4double world_sizeZ;

	double rho_pvt = 1.023 * g / cm3;

	G4ThreeVector R0 = G4ThreeVector(R, 0, 0);
	G4ThreeVector OFV = G4ThreeVector(Ox, Oy, 0);

	// should move to protected
	G4Material* EJ208;
	G4VPhysicalVolume* World;
	int MPT_UPDATE = 0;
	string name = "Undefined";
protected:
	

	const G4int n3 = 21;

public:
	DetectorConstruction() : G4VUserDetectorConstruction()
	{
		world_sizeXY = 237.958; //HARDCODED -NEED TO FIX
		world_sizeZ = 237.958;
	};    //virtual ~DetectorConstruction();

	void copyConstruct(DetectorConstruction* det)
	{
		this->EJ208 = (det->EJ208);
		this->World = (det->World);
		this->MPT_UPDATE = (det->MPT_UPDATE);
	};
	
	G4Material* GetMaterial()
	{
		return EJ208;
	}

	static G4double Interpolation_Calculate(G4double x, G4int bin, const vector<double>& points, const vector<double>& data)
	{
		//G4cout << "G4LinInterpolation is performed (2 arguments)" << G4endl;
		G4int nBins = data.size() - 1;
		G4double value = 0.;
		if (x < points[0])
		{
			value = 0.;
		}
		else if (bin < nBins)
		{
			G4double e1 = points[bin];
			G4double e2 = points[bin + 1];
			G4double d1 = data[bin];
			G4double d2 = data[bin + 1];
			value = d1 + (d2 - d1) * (x - e1) / (e2 - e1);
		}
		else
		{
			value = data[nBins];
		}
		return value;
	}
	static G4double Interpolate(G4double x, const vector<double>& points, const vector<double>& data)
	{
		int bin = -1;
		for (auto p : points)
		{
			if (x < p)
			{
				break;
			}
			bin++;
		}
		return Interpolation_Calculate(x, bin, points, data);
	}
	static bool sipmQE_Hit(G4double wavelength)
	{
		G4double p = Interpolate(wavelength, sipm_QE_lambda, sipm_QE) * 0.01; //probability
		if (G4UniformRand() > p)
		{
			return false;
		}
		return true;
	}

	G4ThreeVector* GetPosition(int nVertex, string name, G4GDMLParser* parser)
	{
		double x = 0;
		double y = 0;
		double z = 0;
		for (int j = 0; j < nVertex; j++)
		{
			name = name + "_v" + to_string(j);
			G4ThreeVector parsedPosition = parser->GetPosition(name);
			x = x + (parsedPosition.x() / nVertex);
			y = y + (parsedPosition.y() / nVertex);
			z = z + (parsedPosition.z() / nVertex);
		}
		return new G4ThreeVector(x, y, z);
	}
	G4ThreeVector GetOrigin(vector<G4ThreeVector>& positions)
	{
		double x = 0;
		double y = 0;
		double z = 0;
		double nVertex = positions.size();
		for (int i = 0; i < nVertex; i++)
		{
			G4ThreeVector position = positions[i];
			x = x + (double)(position.x() / nVertex);
			y = y + (double)(position.y() / nVertex);
			z = z + (double)(position.z() / nVertex);
			//cout << x << y << z << endl;
		}
		return G4ThreeVector(-x, -y, -z);
	}
	G4ThreeVector* Rotate(int i)
	{
		double xy = sqrt(pow(R0.x(), 2) + pow(R0.y(), 2));
		double angle = i * (2 * M_PI / NArray) + atan2(R0.y(), R0.x());
		return new G4ThreeVector(xy * cos(angle), xy * sin(angle), R0.z());
	}
	G4int copyNumber(G4int i)
	{
		return (NAssemblyElements + 1) * i + copyNumberOffset;
	}
	G4double length(G4double v[], int size) {
		G4double len = 0;
		for (int i = 0; i < size; i++)
		{
			len += (v[i]);
		}
		return len;
	}
	void normalize(G4double v[], int size) {
		G4double len = length(v, size);
		for (int i = 0; i < size; i++)
		{
			v[i] = v[i] / len;
		}
	}
	vector<string> split(const string& s, char delim) {
		vector<string> result;
		stringstream ss(s);
		string item;
		while (getline(ss, item, delim)) {
			result.push_back(item);
		}
		return result;
	}

	void UPDATE_GEO_MPT()
	{
		MPT_UPDATE = 1;
		G4RunManager* runManager = G4RunManager::GetRunManager();
		// Open geometry for the physical volume to be modified ...
		G4GeometryManager* geo = G4GeometryManager::GetInstance();
		geo->OpenGeometry(World);
		// Modify dimension of the solid ...
		//
		//World->GetLogicalVolume()->ClearDaughters();
		//World->GetLogicalVolume()->Clean();
		//Parser->Clear();
		//Parser->Read("gdml.gdml",false);
		//World = Parser->GetWorldVolume();
		//WorldBuild((*Parser), World);
		WorldMod();


		//#ifdef HumanoidPhantom
		//	#ifndef ICRPModel
  // 				phantomWorld = new G4HumanPhantomConstruction();
		//	#else
		//		phantomWorld = new ICRP110PhantomConstruction(det);
		//	#endif
		//	UpdateWorldDataWithPhantom(phantomWorld->Construct());
		//#endif // HumanoidPhantom



		// Close geometry for the portion modified ...
		//
		geo->CloseGeometry(World);
		G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();
		visManager->GeometryHasChanged();
		runManager->PhysicsHasBeenModified();
		runManager->GeometryHasBeenModified();
		//runManager->ReinitializeGeometry(true);		
		//EJ208->GetMaterialPropertiesTable()->DumpTable();
		//Vikuiti->GetMaterialPropertiesTable()->DumpTable();
		//Air->GetMaterialPropertiesTable()->DumpTable();
	}

protected:

	void WorldMod()
	{
		G4cout << "--//-------------------------------------------WORLDMODDING-------------------------------------------//--" << G4endl;

		//if visible, (3.105*eV to 1.774*eV), Att length = 400NM! - override.
	//----------------------------------------------------------------------ATT READIN------------------------------------------------////
		vector<G4double> att_length_ejRV{ vector<G4double>() };
		vector<G4double> PhotonEnergyRV{ vector<G4double>() };
		vector<G4double> NPhotonEnergyRV{ vector<G4double>() };
		vector<G4double> rayScr_length_pvtRV{ vector<G4double>() };
		string line;
		int sz0 = 0;
		char delm = '|';
		string fn = "";
		if (MPT_UPDATE < 0)
		{
			fn = "PVT.txt";
		}
		else
		{
			fn = "lambdas5.txt";
			delm = ' ';
		}
		ifstream myfile(fn);
		G4cout << "OPENING PVT FILE -----------";
		if (myfile.is_open())
		{
			if (MPT_UPDATE < 0)
			{
				getline(myfile, line);
				getline(myfile, line);
				getline(myfile, line);
			}
			G4cout << G4endl;
			while (getline(myfile, line))
			{
				G4double en = 0;
				G4double rayScr_CS = 0;
				G4double att_CS = 0;
				G4double rS_l = 0;
				G4double att_l = 0;
				vector<string> val = split(line, delm);
				if (MPT_UPDATE < 0)
				{
					en = G4double(stod(val[0]));
					PhotonEnergyRV.push_back(en);
					NPhotonEnergyRV.push_back(en);
					rayScr_CS = G4double(stod(val[1])) * cm2 / g;
					att_CS = G4double(stod(val[4])) * cm2 / g;
					rS_l = 1 / (rho_pvt * rayScr_CS);
					att_l = 1 / (rho_pvt * att_CS);
					G4cout << "1ST SET" << G4endl;
				}
				else
				{
					en = G4double(stod(val[0])) * MeV;
					//G4cout << val.size() << " " << en << G4endl; //size
					PhotonEnergyRV.push_back(en);
					NPhotonEnergyRV.push_back(en);
					rayScr_CS = G4double(stod(val[4])) * g / cm2; //
					att_CS = G4double(stod(val[5])) * g / cm2;
					rS_l = rayScr_CS / (rho_pvt);
					att_l = att_CS / (rho_pvt);
					//G4cout << "2nd SET" << G4endl;
					if ((en < (3.105 * eV)) && (en > (1.774 * eV)))
					{
						att_l = 400 * cm;
						//G4cout << "MODIFIED OPTICAL ATTENUATION." << G4endl;
					}


				}

				rayScr_length_pvtRV.push_back(rS_l);
				att_length_ejRV.push_back(att_l);
				//G4cout << "values  : " << val[0] << " "<<rS_l/cm<<" "<<att_l/cm<<G4endl;
				sz0 += 1;
			}
			myfile.close();
			G4cout << "opened file" << G4endl;
		}
		else { G4cout << "Unable to open file" << G4endl; exit(3); }

		G4double att_length_ejR[sz0];
		G4double PhotonEnergyR[sz0];
		G4double rayScr_length_pvtR[sz0];
		G4double NPhotonEnergyR[sz0];
		memcpy(att_length_ejR, &att_length_ejRV[0], sizeof(double) * att_length_ejRV.size());
		memcpy(PhotonEnergyR, &PhotonEnergyRV[0], sizeof(double) * PhotonEnergyRV.size());
		memcpy(rayScr_length_pvtR, &rayScr_length_pvtRV[0], sizeof(double) * rayScr_length_pvtRV.size());
		memcpy(NPhotonEnergyR, &NPhotonEnergyRV[0], sizeof(double) * NPhotonEnergyRV.size());
		//----------------------------------------------------------------------RED MAT----------------------------------------------------////
		EJ208->GetMaterialPropertiesTable()->RemoveProperty("ABSLENGTH");
		EJ208->GetMaterialPropertiesTable()->RemoveProperty("RAYLEIGH");
		EJ208->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH", PhotonEnergyR, att_length_ejR, sz0);// MODIFIED WITH DATA
		EJ208->GetMaterialPropertiesTable()->AddProperty("RAYLEIGH", PhotonEnergyR, rayScr_length_pvtR, sz0);// MODIFIED WITH DATA
	//G4NistManager* man = G4NistManager::Instance();
		for (int i = 0; i < World->GetLogicalVolume()->GetNoDaughters(); i++)
		{
			G4VPhysicalVolume* V = World->GetLogicalVolume()->GetDaughter(i);
			if (V != NULL)
			{
				string name = V->GetName();
#ifdef SingleStrip
				string name0 = "EJ208(1)";
#else
				string name0 = "EJ208";
#endif
				if (name.find(name0) != std::string::npos)
				{
					//V->GetLogicalVolume()->SetMaterial(man->FindMaterial("G4_AIR"));
					V->GetLogicalVolume()->SetMaterial(EJ208);

#ifdef SingleStrip
					G4cout << "updated EJ208!" << G4endl;
#endif
				}
			}
		}
		//World->GetLogicalVolume()->GetDaughter(0)->GetLogicalVolume()->GetMaterial()->GetMaterialPropertiesTable()->DumpTable();
	}

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

