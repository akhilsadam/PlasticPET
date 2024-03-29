#ifndef _GDMLDETECTORCONSTRUCTION_H_
#define _GDMLDETECTORCONSTRUCTION_H_
#include <Version.hh>
#include <DetectorConstruction.hh>
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

//#ifdef HumanoidPhantom
//	#ifndef ICRPModel
//		#include "G4HumanPhantomConstruction.hh"
//	#else
//		#include "ICRP110PhantomConstruction.hh"
//	#endif
//#endif

class G4VPhysicalVolume;
class G4LogicalVolume;


//#include "B3DetectorConstruction.hh"
class GDMLDetectorConstruction : public DetectorConstruction
{
	
  public:
	G4Material* Vikuiti;
	G4Material* Air;
	G4Material* TPB;
	G4Material* detMAT;
	
	void WorldBuild(G4GDMLParser& parser, G4VPhysicalVolume *setWorld )
   {
	   G4cout << "--//-------------------------------------------WORLDBUILDING-------------------------------------------//--" << G4endl;
// ---- get positions
		// vector<G4ThreeVector> positions;
		// int i = 0;
		// while(true)
		// {
		// 	G4VPhysicalVolume* V = World->GetLogicalVolume()->GetDaughter(i);
		// 	if(V!=NULL)
		// 	{			
		// 		cout << V->GetTranslation() << endl;
		// 		positions.push_back(V->GetTranslation());
		// 	}
		// 	else
		// 	{
		// 		break;
		// 	}
		// 	i+=1;
		// }
	//cout << positions.size() << endl;


//----------------------------------------------------------------------DEF MAT----------------------------------------------------////
	G4NistManager* man = G4NistManager::Instance();

	Air = man->FindOrBuildMaterial("G4_AIR");
	detMAT = man->FindOrBuildMaterial("G4_AIR");
  	G4bool isotopes = false;
  
  	G4Element*  H = man->FindOrBuildElement("H" , isotopes);
	G4Element*  C = man->FindOrBuildElement("C" , isotopes);
	G4Element*  Pb = man->FindOrBuildElement("Pb" , isotopes);

	G4Material* CarbonM = new G4Material("CarbonM", 2.26*g/cm3, 1);
	CarbonM->AddElement(C,1);

	G4Material* LeadM = new G4Material("LeadM", 207.2*g/cm3, 1);
	LeadM->AddElement(Pb,1);

	Vikuiti = new G4Material("Vikuiti", 1.07*g/cm3, 2);
	Vikuiti->AddElement(C,8);
	Vikuiti->AddElement(H,8);
	Vikuiti->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

	G4Material* PVT = new G4Material("PVT", rho_pvt, 2);
	PVT->AddElement(C,9);
	PVT->AddElement(H,10);
	PVT->GetIonisation()->SetBirksConstant(0.1955*mm/MeV);

	TPB = new G4Material("TPB", 1.079*g/cm3, 2);
	TPB->AddElement(C,28);
	TPB->AddElement(H,22);
	//TPB->GetIonisation()->SetBirksConstant(0.1955*mm/MeV); need to set Birks' constant?

//----------------------------------------------------------------------MIX MAT----------------------------------------------------////
	#ifndef PVTTest
		EJ208 = new G4Material("EJ208",11.33*g/cm3,2);
		EJ208->AddMaterial(PVT,95*perCent);
		EJ208->AddMaterial(LeadM,5*perCent);
	#else
		delete EJ208;
		EJ208 = PVT;
	#endif

//----------------------------------------------------------------------MAT TABLES------------------------------------------------////
const G4int n = 6;
const G4int n2 = 10;
G4double PhotonEnergy[n] = {3.105*eV,2.95714*eV,2.855*eV,2.7*eV,2.5875*eV,2.388*eV}; //visible spectrum (400,420,435,460,480,520)nm
G4double GPhotonEnergy[n+1] = {511*keV, 3.105 * eV,2.95714 * eV,2.855 * eV,2.7 * eV,2.5875 * eV,2.388 * eV };
G4double TPBPhotonEnergy[n2] = {12.3985*eV,6.19926*eV,4.1328*eV,3.5424*eV,3.105*eV,2.95714*eV,2.855*eV,2.7*eV,2.5875*eV,2.388*eV};//visible spectrum (100,200,300,350,400,420,435,460,480,520)nm

G4double refractive_index_vk[n] = {1.6,1.6,1.6,1.6,1.6,1.6};
G4double att_length_vk[n] = {400*cm,400*cm,400*cm,400*cm,400*cm,400*cm};

G4double refractive_index_ej[n+1] = {1,1.58,1.58,1.58,1.58,1.58,1.58};
G4double att_length_ej[n] = {400*cm,400*cm,400*cm,400*cm,400*cm,400*cm};

//Spectroscopic and travelling-wave lasing characterisation of tetraphenylbenzidine and di-naphtalenyl-diphenylbenzidine
//Measurements of the intrinsic quantum efficiency and absorption length of tetraphenyl butadiene thin films in the vacuum ultraviolet regime
G4double refractive_index_tpb[n2] = {1.67,1.67,2.2,1.79,1.9,1.8,1.8,1.79,1.75,1.72};// the following values are preliminary from articles above - needs updating!
G4double att_length_tpb[n2] = {100*nm,100*nm,100*nm,31.6*nm,2000*nm,10000*nm,20000*nm,1*mm,1*mm,1*mm};
G4double ray_length_tpb[n2] = {10*m,10*m,10*m,2.75*um,2.75*um,2.75*um,2.75*um,2.75*um,2.75*um,2.75*um};
G4double em_tpb[n2] = {0,0,0,0,0.00665,0.013,0.012,0.008,0.004,0.002};
normalize(em_tpb,n2);
G4double tc_tpb = 3*ns;
G4double reflectivity_tpb[n2] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
G4double efficiency2[n2] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

//----------------------------------------------------------------------MAT TABLES------------------------------------------------////
//----------------------------------------------------------------------ATT READIN------------------------------------------------////
	vector<G4double> att_length_ejRV(begin(att_length_ej), end(att_length_ej));
	vector<G4double> PhotonEnergyRV(begin(PhotonEnergy), end(PhotonEnergy));
	vector<G4double> NPhotonEnergyRV = *new vector<G4double>();
	vector<G4double> rayScr_length_pvtRV = *new vector<G4double>();
	string line;
	int sze = n;
	int sz0 = 0;
	char delm = '|';
	string fn ="";
	if(MPT_UPDATE<0)
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
		if(MPT_UPDATE<0)
		{
			getline (myfile,line);
			getline (myfile,line);
			getline (myfile,line);
		}
		G4cout << G4endl;
    	while ( getline (myfile,line) )
    	{
			G4double en = 0;
			G4double rayScr_CS = 0;
			G4double att_CS = 0;
			G4double rS_l = 0;
			G4double att_l = 0;
			vector<string> val = split(line,delm);
			if(MPT_UPDATE<0)
			{
				en = G4double(stod(val[0]));
				PhotonEnergyRV.push_back(en);
				NPhotonEnergyRV.push_back(en);
				rayScr_CS = G4double(stod(val[1]))*cm2/g;
				att_CS = G4double(stod(val[4]))*cm2/g;
				rS_l = 1/(rho_pvt*rayScr_CS);
				att_l = 1/(rho_pvt*att_CS);
				G4cout << "1ST SET" << G4endl;
			}
			else
			{
				en = G4double(stod(val[0]));
				//G4cout << val.size() << " " << en << G4endl; //size
				PhotonEnergyRV.push_back(en);
				NPhotonEnergyRV.push_back(en);
				rayScr_CS = G4double(stod(val[4]))*g/cm2; //
				att_CS = G4double(stod(val[5]))*g/cm2;
				rS_l = rayScr_CS/(rho_pvt);
				att_l = att_CS/(rho_pvt);
				//G4cout << "2nd SET" << G4endl;
			}
						
			rayScr_length_pvtRV.push_back(rS_l);
			att_length_ejRV.push_back(att_l);
			//G4cout << "values  : " << val[0] << " "<<rS_l/cm<<" "<<att_l/cm<<G4endl;
			sze+=1;
			sz0+=1;
    	}
    	myfile.close();
		G4cout << "opened file" << G4endl;
	}
  	else {G4cout << "Unable to open file" << G4endl; exit(3);}

	G4double att_length_ejR[sze];
	G4double PhotonEnergyR[sze];
	G4double rayScr_length_pvtR[sz0];
	G4double NPhotonEnergyR[sz0];
	memcpy( att_length_ejR, &att_length_ejRV[0], sizeof(double) * att_length_ejRV.size() );
	memcpy( PhotonEnergyR, &PhotonEnergyRV[0], sizeof(double) * PhotonEnergyRV.size() );
	memcpy( rayScr_length_pvtR, &rayScr_length_pvtRV[0], sizeof(double) * rayScr_length_pvtRV.size() );
	memcpy( NPhotonEnergyR, &NPhotonEnergyRV[0], sizeof(double) * NPhotonEnergyRV.size() );	
	
//----------------------------------------------------------------------MAT TABLES------------------------------------------------////
	G4double refractive_index_air[n] = {1.00029,1.00029,1.00029,1.00029,1.00029,1.00029};

	G4double fast[n] = {0.032258,0.258064,0.322581,0.225806,0.129032,0.032258};//{0.1,0.8,1,0.7,0.4,0.1}
	G4double reflectivity_det[n] = {0,0,0,0,0,0};

	#ifndef SSRefractionTest
		#ifdef ULTRAreflective
			G4double reflectivity_vk[n] = {ref_VK,ref_VK,ref_VK,ref_VK,ref_VK,ref_VK};//{0.5,0.5,0.5,0.5,0.5,0.5};//
			cout << endl;
			cout << "****USING REFLECTIVITY AS DEFINED BY COMMANDLINE****: ref_VK=" << ref_VK << endl;
			cout << endl;
		#else
			G4double reflectivity_vk[n] = {.9150,.9334,.9452,.9566,.9652,.9698}; //using ESR_Clear. DRP:{0.9643,0.9680,0.9698,0.9743,0.9761,0.9798};
		#endif
	#else
		G4double reflectivity_vk[n] = {1,1,1,1,1,1};
		//G4double reflectivity_vk[n] = {0,0,0,0,0,0}; //using ESR_Clear. DRP:{0.9643,0.9680,0.9698,0.9743,0.9761,0.9798};
	#endif


	#ifndef ULTRAreflective
	//G4double reflectivity_ej[n] =  {.9150,.9334,.9452,.9566,.9652,.9698}; //USING VK for now!!
	G4double reflectivity_ej[n] =  {1,1,1,1,1,1};
	#else
	G4double reflectivity_ej[n] =  {1,1,1,1,1,1};//{0.5,0.5,0.5,0.5,0.5,0.5};//
	#endif
	
	G4double efficiency[n] = {1.0,1.0,1.0,1.0,1.0,1.0};

	G4MaterialPropertiesTable* airMPT = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* vkMPT = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* scintMPT = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* detMPT = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* surfVKMPT = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* surfEJMPT = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* surfEJMPTR = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* tpbMPT = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* surfTPBMPT = new G4MaterialPropertiesTable();

	airMPT->AddProperty("RINDEX", PhotonEnergy,refractive_index_air,n);
	vkMPT->AddProperty("RINDEX", PhotonEnergy,refractive_index_vk,n);
	detMPT->AddProperty("RINDEX", PhotonEnergy,refractive_index_air,n);
	vkMPT->AddProperty("ABSLENGTH", PhotonEnergy,att_length_vk,n);
	scintMPT->AddProperty("RINDEX", GPhotonEnergy,refractive_index_ej,n+1);

	tpbMPT->AddProperty("RINDEX", TPBPhotonEnergy,refractive_index_tpb,n2);
	tpbMPT->AddProperty("WLSABSLENGTH", TPBPhotonEnergy,att_length_tpb,n2);
	tpbMPT->AddProperty("WLSCOMPONENT", TPBPhotonEnergy,em_tpb,n2);
	tpbMPT->AddConstProperty("WLSTIMECONSTANT",tc_tpb);
	//tpbMPT->AddProperty("RAYLEIGH", TPBPhotonEnergy,ray_length_tpb,n2); should we use Rayleigh?

	if(MPT_UPDATE<0)
	{
		scintMPT->AddProperty("ABSLENGTH",PhotonEnergyR,att_length_ejR,sze);// MODIFIED WITH DATA
		//scintMPT->AddProperty("RAYLEIGH",NPhotonEnergyR,rayScr_length_pvtR,sz0); turn this off for now!
	}
	else
	{
		scintMPT->RemoveProperty("ABSLENGTH");
		scintMPT->RemoveProperty("RAYLEIGH");
		scintMPT->AddProperty("ABSLENGTH",PhotonEnergyR,att_length_ejR,sze);// MODIFIED WITH DATA
	}

	#ifndef ScintillationDisable
	scintMPT->AddProperty("FASTCOMPONENT", PhotonEnergy, fast, n);
  	scintMPT->AddConstProperty("SCINTILLATIONYIELD", 9200 / MeV);
  	scintMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
	scintMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 1.0 * ns);
 	scintMPT->AddConstProperty("FASTTIMECONSTANT", 3.3 * ns);
	scintMPT->AddConstProperty("YIELDRATIO", 1.0);
	vkMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity_vk, n);
	scintMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity_ej, n); //WHAT is the EJ reflectivity?
	#endif

	Air->SetMaterialPropertiesTable(airMPT);
	EJ208->SetMaterialPropertiesTable(scintMPT);
	Vikuiti->SetMaterialPropertiesTable(vkMPT);


	detMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity_det, n);
	detMAT->SetMaterialPropertiesTable(detMPT);

	/*if(MPT_UPDATE>0)
	{
		WorldMod();
		return;
	}*/
//----------------------------------------------------------------------SETSURF----------------------------------------------------////
	surfVKMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity_vk, n);
	surfEJMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity_ej, n); // WHAT is the EJ reflectivity?
	surfEJMPTR->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity_ej, n); 
	surfTPBMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity_tpb, n2); 
	surfVKMPT->AddProperty("EFFICIENCY", PhotonEnergy, efficiency, n);
	surfEJMPT->AddProperty("EFFICIENCY", PhotonEnergy, efficiency, n);
	surfEJMPTR->AddProperty("EFFICIENCY", PhotonEnergy, efficiency, n);
	surfTPBMPT->AddProperty("EFFICIENCY", PhotonEnergy, efficiency2, n2); //NO Quantum efficiency implementation yet!
	surfVKMPT->AddProperty("RINDEX", PhotonEnergy,refractive_index_vk,n);
	surfEJMPT->AddProperty("RINDEX", PhotonEnergy,refractive_index_ej,n);
	surfEJMPTR->AddProperty("RINDEX", PhotonEnergy,refractive_index_ej,n);
	surfTPBMPT->AddProperty("RINDEX", PhotonEnergy,refractive_index_tpb,n2);
	
	G4OpticalSurface* surfVK = new G4OpticalSurface("surfVK");
	G4OpticalSurface* surfEJ = new G4OpticalSurface("surfEJ");
	G4OpticalSurface* surfEJR = new G4OpticalSurface("surfEJR");
	G4OpticalSurface* surfTPB = new G4OpticalSurface("surfTPB");

	G4double sigma_alpha_tpb = 0.01;
	#ifndef ReflectionDisable
	G4double sigma_alpha_VK = 0.0226893; //0.209439; //12deg., ground //..  vk
	G4double sigma_alpha_endcap_EJ = 0.0226893; //12deg., ground
	#ifndef SSSpecularReflectionTest
		G4double sigma_alpha_EJ = 0.0226893; //1.3 deg., polished 
	#else
		G4double sigma_alpha_EJ = 0.0;
	#endif

	#ifdef VK_SIGMAALPHA
		sigma_alpha_VK = sga_VK;
		cout << endl;
		cout << "****USING VK SIGMAALPHA AS DEFINED BY COMMANDLINE****: sga_VK=" << sigma_alpha_VK << endl;
		cout << endl;
	#endif

	#ifdef EJ_SIGMAALPHA
		sigma_alpha_EJ = sga_EJ;
		cout << endl;
		cout << "****USING EJ SIGMAALPHA AS DEFINED BY COMMANDLINE****: sga_EJ=" << sigma_alpha_EJ << endl;
		cout << endl;
	#endif

	G4double specularlobe[n] = {1.0,1.0,1.0,1.0,1.0,1.0}; //all specular lobe (microfacet reflections) // test with specular spike (perfectly smooth surface!)
	G4double specularspike[n] = {0.0,0.0,0.0,0.0,0.0,0.0};
	G4double backscatter[n] = {0.0,0.0,0.0,0.0,0.0,0.0};
	G4double specularlobe2[n2] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}; //all specular lobe (microfacet reflections) // test with specular spike (perfectly smooth surface!)
	G4double specularspike2[n2] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	G4double backscatter2[n2] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	surfVKMPT->AddProperty("SPECULARLOBECONSTANT",PhotonEnergy,specularlobe,n);
	surfVKMPT->AddProperty("SPECULARSPIKECONSTANT",PhotonEnergy,specularspike,n);
	surfVKMPT->AddProperty("BACKSCATTERCONSTANT",PhotonEnergy,backscatter,n);
	surfEJMPT->AddProperty("SPECULARLOBECONSTANT",PhotonEnergy,specularlobe,n);
	surfEJMPT->AddProperty("SPECULARSPIKECONSTANT",PhotonEnergy,specularspike,n);
	surfEJMPT->AddProperty("BACKSCATTERCONSTANT",PhotonEnergy,backscatter,n);
	surfEJMPTR->AddProperty("SPECULARLOBECONSTANT",PhotonEnergy,specularlobe,n);
	surfEJMPTR->AddProperty("SPECULARSPIKECONSTANT",PhotonEnergy,specularspike,n);
	surfEJMPTR->AddProperty("BACKSCATTERCONSTANT",PhotonEnergy,backscatter,n);
	surfTPBMPT->AddProperty("SPECULARLOBECONSTANT",PhotonEnergy,specularlobe2,n2);
	surfTPBMPT->AddProperty("SPECULARSPIKECONSTANT",PhotonEnergy,specularspike2,n2);
	surfTPBMPT->AddProperty("BACKSCATTERCONSTANT",PhotonEnergy,backscatter2,n2);
	surfVK->SetType(dielectric_metal);
	surfVK->SetModel(unified);
	surfVK->SetFinish(ground);
	surfVK->SetSigmaAlpha(sigma_alpha_VK);
	surfVK->SetMaterialPropertiesTable(surfVKMPT);
	surfEJ->SetType(dielectric_dielectric);
	surfEJ->SetModel(unified);
	surfEJ->SetFinish(ground);
	surfEJ->SetSigmaAlpha(sigma_alpha_EJ);
	surfEJ->SetMaterialPropertiesTable(surfEJMPT);
	surfEJ->SetType(dielectric_dielectric);
	surfEJR->SetModel(unified);
	surfEJR->SetFinish(ground);
	surfEJR->SetSigmaAlpha(sigma_alpha_endcap_EJ);
	surfEJR->SetMaterialPropertiesTable(surfEJMPTR);
	surfTPB->SetType(dielectric_dielectric);
	surfTPB->SetModel(unified);
	surfTPB->SetFinish(ground);
	surfTPB->SetSigmaAlpha(sigma_alpha_tpb);
	surfTPB->SetMaterialPropertiesTable(surfTPBMPT);
	#endif

	#ifdef LYSOTest
		LYSO* lyse = new LYSO();
		EJ208 = lyse->LYSOmat();
		//printf("\n\n LYSO NAME: %s",EJ208->GetName()); HAS AN ERROR
		surfEJ = lyse->LYSOsurf();
	#endif	
//----------------------------------------------------------------------Geometry---------------------------------------------------////

	G4VisAttributes* red_col = new G4VisAttributes(G4Color(0.6,0.4,0.4,1));
	G4VisAttributes* end_col = new G4VisAttributes(G4Color(0.9,0.4,0.4,1));
	G4VisAttributes* blue_col = new G4VisAttributes(G4Color(0.4,0.4,0.6,1));
	G4VisAttributes* det_col = new G4VisAttributes(G4Color(0.8,0.8,1.0,1));
	G4VisAttributes* tpb_col = new G4VisAttributes(G4Color(0.9,0.0,0.0,1));
	G4VisAttributes* invis_col = new G4VisAttributes(G4Color(0.0,0.0,0.0,0));


	#ifdef SSRefractionTest
		G4Box* world = new G4Box("World",(20*length_D),(20*length_D),(20*length_D));
		G4LogicalVolume* logicW = new G4LogicalVolume(world,EJ208,"World");
		logicW->SetMaterial(Air);
		logicW->SetVisAttributes (invis_col);
		World = new G4PVPlacement(0,G4ThreeVector(), logicW,"World", 0, false,0);
	#else
		G4LogicalVolume* VOL = World->GetLogicalVolume();
		VOL->SetMaterial(Air);
	#endif

	#ifdef MultipleStripCell

		#ifdef RoughEnds
			G4double edX = sx; G4double edY = sy; G4double edZ = 0.0000001*cm;
			G4Box* endcap = new G4Box("E_EJ208", edX/2, edY/2, edZ/2);	

			vector<G4ThreeVector> stripPos;

			for(int y = 0; y<(Ny); y++)
			{
				int numSiPM = (int) (y/4);
				double VKOffset = (2*numSiPM)*VKT;
				double yi = (9.66*cm) - ((y)*py + VKOffset);
				for(int x = 0; x<(Nx); x++)
				{
					VKOffset = (2*x+2)*VKT;
					double xi = (x+1)*px + VKOffset - (Dx);
					stripPos.push_back(G4ThreeVector(xi,yi,-(edZ/2)));
					stripPos.push_back(G4ThreeVector(xi,yi,(length_Z + (edZ/2))));
				}
			}

			G4LogicalSkinSurface* EJsurfR;
			for(G4ThreeVector posR : stripPos)
			{
				G4LogicalVolume* endLV = new G4LogicalVolume(endcap,EJ208,"E_EJ208");
				endLV->SetVisAttributes (end_col);
				new G4PVPlacement(0,posR,"end",endLV,World, false,0,fCheckOverlaps);       // checking overlaps 
				EJsurfR = new G4LogicalSkinSurface("E_EJ208",endLV, surfEJR); 
			}
		#endif

	// using https://geant4.web.cern.ch/sites/geant4.web.cern.ch/files/geant4/collaboration/workshops/users2002/talks/lectures/geoadvanced.pdf
	G4AssemblyVolume* Array =  new G4AssemblyVolume();
	G4ThreeVector position = *(new G4ThreeVector());
	// G4ThreeVector rotation;
	G4ThreeVector origin = G4ThreeVector(0,0,0);//GetOrigin(positions);
	position = origin - OFV;
	// cout << origin.x() << origin.y() << origin.z() << endl;
	// cout << position.x() << position.y() << position.z() << endl;
	G4RotationMatrix* rotateMatrix = new G4RotationMatrix();

	//G4PhysicalVolumeStore* STOCK = G4PhysicalVolumeStore::GetInstance();
	// vector<G4RotationMatrix> rotationMatrices;

	// vector<G4LogicalVolume*> scintillatorVOL;
	// vector<G4LogicalVolume*> separatorVOL;
	// vector<G4LogicalVolume*> sipmVOL;

	//note assuming pre-translated shapes (8 vertex rectangular prism and 16 vertex rectangular pipe)

	G4LogicalVolume* EJvol;
	G4VPhysicalVolume* EJvol_P;
	G4LogicalSkinSurface* EJsurf;
	G4String a;
	for (G4int i = 1; i <= NUnit; i++) 
	{
        //char a[100];
		a = "S_EJ208(";
	 	a = a + to_string(i);
	  	a = a + ")";
        //sprintf(a,"%d", i);
        //strcat(a,")");
        EJvol = parser.GetVolume(a);
		// adding roughened end caps to help with detection?
		// positions.push_back(*(GetPosition(8,"",&parser))); -----
		// position = parser.GetPosition(a);
		// rotation = parser.GetRotation(a);
		// rotateMatrix = new G4RotationMatrix();
		// rotateMatrix->rotateX(rotation.x());
		// rotateMatrix->rotateY(rotation.y());
		// rotateMatrix->rotateZ(rotation.z());
		
		/// 
	  	EJvol->SetMaterial(EJ208);
	  	EJvol->SetVisAttributes (red_col);
	  	EJsurf = new G4LogicalSkinSurface("surfEJ_L",EJvol, surfEJ);

		Array->AddPlacedVolume(EJvol,position,rotateMatrix);
	}

	G4LogicalVolume* VKvol;
	G4VPhysicalVolume* VKvol_P;
	G4LogicalSkinSurface* VKsurf;
	for (G4int i = 1; i <= NUnit; i++) 
	{
        //char a[100];
		a = "S_Vikuiti(";
	  	a = a + to_string(i);
	  	a = a + ")";
        //sprintf(a,"%d", i);
        //strcat(a,")");
        VKvol = parser.GetVolume(a);
		#ifndef DISABLEVK
			// position = parser.GetPosition(a);
			// rotation = parser.GetRotation(a);
			// rotateMatrix = new G4RotationMatrix();
			// rotateMatrix->rotateX(rotation.x());
			// rotateMatrix->rotateY(rotation.y());
			// rotateMatrix->rotateZ(rotation.z());

			VKvol->SetMaterial(Vikuiti);
			VKvol->SetVisAttributes (blue_col);
			VKsurf = new G4LogicalSkinSurface("surfVK_L",VKvol, surfVK);
		#else
			VKvol->SetMaterial(Air);
			VKvol->SetVisAttributes (invis_col);
		#endif

		Array->AddPlacedVolume(VKvol,position,rotateMatrix);
    }

	//--------- Detecting (SiPM like) Endcap Geometry -------------//\\

	G4double dX = 7.74*cm; G4double dY = 10.32*cm; G4double dZ = 1*cm;
	G4Box* det = new G4Box("det", dX/2, dY/2, dZ/2);
	//G4Box* ydet = new G4Box("ydet", 0.5*cm, dY/2, length_D/2);
	G4LogicalVolume* logicDetL = new G4LogicalVolume(det,detMAT,detectorLeftName);
	G4LogicalVolume* logicDetR = new G4LogicalVolume(det,detMAT,detectorRightName);
	//G4LogicalVolume* logicDetY = new G4LogicalVolume(ydet,EJ208,"detVOLY");
	logicDetL->SetVisAttributes (det_col);
	logicDetR->SetVisAttributes (det_col);
	//logicDetY->SetVisAttributes (det_col);
    G4ThreeVector pos = G4ThreeVector(Ox,Oy,-0.5001*cm)+position;
	G4ThreeVector pos2 = G4ThreeVector(Ox,Oy,(length_D+0.5001*cm))+position;
	//G4ThreeVector pos3 = G4ThreeVector((-2.58*cm-(dX/2)-0.5001*cm),4.83*cm,(length_D/2));
	//new G4PVPlacement(0,pos,"det",logicDetR,World, false,0,fCheckOverlaps);       // checking overlaps 
	//new G4PVPlacement(0,pos2,"det",logicDetL,World, false,0,fCheckOverlaps); 
	//new G4PVPlacement(0,pos3,"ydet",logicDetY,World, false,0,fCheckOverlaps); 
	
	//rotateMatrix = new G4RotationMatrix();
	Array->AddPlacedVolume(logicDetR,pos, rotateMatrix);
	Array->AddPlacedVolume(logicDetL,pos2, rotateMatrix);

		
	World->GetLogicalVolume()->ClearDaughters();
	
	for( unsigned int i = 0; i < Na; i++ ) 
	{
		// Translation of the assembly inside the world
		G4ThreeVector arrayPosition = ((*Rotate(i)));
		G4ThreeVector& positionRef = arrayPosition;
		rotateMatrix = new G4RotationMatrix();
		rotateMatrix->rotateZ(i*(double)(2*M_PI/NArray));
		Array->MakeImprint(World->GetLogicalVolume(),positionRef,rotateMatrix,copyNumber(i),false);
	}

	for(int i = 0; i< World->GetLogicalVolume()->GetNoDaughters(); i++)
	{
		G4VPhysicalVolume* V = World->GetLogicalVolume()->GetDaughter(i);
		if(V!=NULL)
		{			
			//cout << V->GetName() << endl;
			string name = V->GetName();
			if(name.find(scintillatorName) != std::string::npos)
			{	
				G4LogicalVolume* logicalVol = V->GetLogicalVolume();
				EJsurf = new G4LogicalSkinSurface("surfEJ_L",logicalVol, surfEJ);
			}
			else if(name.find(vikuitiName) != std::string::npos)
			{
				G4LogicalVolume* logicalVol = V->GetLogicalVolume();
				VKsurf = new G4LogicalSkinSurface("surfVK_L",logicalVol, surfVK);
			}
		}
	}
	// G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();
	// G4RunManager* runManager = G4RunManager::GetRunManager();
	// visManager->GeometryHasChanged();
	// runManager->GeometryHasBeenModified();

	#endif
	#ifdef SingleStrip

		#ifdef LEGEND
			G4Box* tpb= new G4Box("TPB",2.54*cm/2.0,5*um/2.0,50*cm);
			G4LogicalVolume* lc = new G4LogicalVolume(tpb,TPB,"TPB");
			G4ThreeVector postpb = G4ThreeVector(-5.14*cm,(9.97*cm+2.5*um),50*cm);
			lc->SetVisAttributes (tpb_col);
			G4LogicalSkinSurface* TPBsurf = new G4LogicalSkinSurface("surfTPB_L",lc, surfTPB); 	
			G4VPhysicalVolume* pC = new G4PVPlacement(0,postpb,"TPB", lc, World, false,0); //TPB film
		#endif

	#ifndef SSRefractionTest
	

	G4LogicalVolume* EJvol;
	G4VPhysicalVolume* EJvol_P;
	G4LogicalSkinSurface* EJsurf;
	G4String a;
	
	for (G4int i = 1; i <= 1; i++) {
        //char a[100];
		a = "S_EJ208(";
	  	a = a + to_string(i);
	  	a = a + ")";
        //sprintf(a,"%d", i);
        //strcat(a,")");
        EJvol = parser.GetVolume(a);
	  	EJvol->SetMaterial(EJ208);
	  	EJvol->SetVisAttributes (red_col);
	  	EJsurf = new G4LogicalSkinSurface("surfEJ_L",EJvol, surfEJ); 	
    }


	for (G4int i = 2; i <= NUnit; i++) {
        //char a[100];
		a = "S_EJ208(";
	  	a = a + to_string(i);
	  	a = a + ")";
        //sprintf(a,"%d", i);
        //strcat(a,")");
        EJvol = parser.GetVolume(a);
	  	EJvol->SetMaterial(Air);
	  	EJvol->SetVisAttributes (invis_col);
	  	//EJsurf = new G4LogicalSkinSurface("surfEJ_L",EJvol, surfEJ); 	
    }
	G4LogicalVolume* VKvol;
	G4VPhysicalVolume* VKvol_P;
	G4LogicalSkinSurface* VKsurf;
	for (G4int i = 1; i <= NUnit; i++) 
	{
        //char a[100];
	  	a = "S_Vikuiti(";
	  	a = a + to_string(i);
	  	a = a + ")";
        //sprintf(a,"%d", i);
        //strcat(a,")");
        VKvol = parser.GetVolume(a);
	  	VKvol->SetMaterial(Air);
	  	VKvol->SetVisAttributes (invis_col);
	  	//VKsurf = new G4LogicalSkinSurface("surfVK_L",VKvol, surfVK);	
	}
	#else
	G4double len = length_D*2;
	G4Box* testP = new G4Box("S_EJ208",len,len,len);
	G4LogicalVolume* EJvol = new G4LogicalVolume(testP,EJ208,"S_EJ208");
	EJvol->SetMaterial(EJ208);
	EJvol->SetVisAttributes (red_col);
	G4ThreeVector posV = G4ThreeVector(0*cm,0*cm,len);
	G4VPhysicalVolume* pC = new G4PVPlacement(0,posV,"S_EJ208", EJvol, World, false,0);
	#endif

	//--------- Detecting (SiPM like) Endcap Geometry -------------//\\

	#ifndef SSRefractionTest
	G4double dX = 2.54*cm; G4double dY = .62*cm; G4double dZ = 1*cm;
	G4Box* det = new G4Box("det", dX/2, dY/2, dZ/2);	
	G4ThreeVector pos = G4ThreeVector(-5.14*cm,9.66*cm,-0.5001*cm);
	G4ThreeVector pos2 = G4ThreeVector(-5.14*cm,9.66*cm,100.5001*cm);
	#else
	G4double dX = 4*length_D; G4double dY = 4*length_D; G4double dZ = 1*cm;
	G4Box* det = new G4Box("det", dX/2, dY/2, dZ/2);	
	G4ThreeVector pos = G4ThreeVector(0*cm,0*cm,-0.5001*cm);
	G4ThreeVector pos2 = G4ThreeVector(0*cm,0*cm,(dY + 0.5001*cm));
	#endif
	G4LogicalVolume* logicDetR = new G4LogicalVolume(det,detMAT,"detVOLR");
	G4LogicalVolume* logicDetL = new G4LogicalVolume(det,detMAT,"detVOLL");
	logicDetR->SetVisAttributes (det_col);
	logicDetL->SetVisAttributes (det_col);
	new G4PVPlacement(0,pos,"det",logicDetR,World, false,0,fCheckOverlaps);       // checking overlaps 
	new G4PVPlacement(0,pos2,"det",logicDetL,World, false,0,fCheckOverlaps); 
	#endif




	/*G4Tubs* testP = new G4Tubs("Crystal",0*cm,10*cm, 5*cm,0,2*pi);
	G4LogicalVolume* lc = new G4LogicalVolume(testP,EJ208,"Crystal");
	G4ThreeVector pos = G4ThreeVector(0*cm,0*cm,0*cm);
	G4VPhysicalVolume* pC = new G4PVPlacement(0,pos,"Crystal", lc, World, false,0); //test crystal

	G4Tubs* testP2 = new G4Tubs("Crystal2",4*cm,10*cm, 5*cm,0,2*pi);
	G4LogicalVolume* lc2 = new G4LogicalVolume(testP2,Vikuiti,"Crystal2");
	G4ThreeVector pos2 = G4ThreeVector(0*cm,8*cm,8*cm);
	G4VPhysicalVolume* pC2 = new G4PVPlacement(0,pos2,"Crystal2", lc2, World, false,0);*/
   }

    GDMLDetectorConstruction(G4GDMLParser& parser, G4VPhysicalVolume *setWorld=0) : DetectorConstruction()
   {

	// G4cout <<""<< G4endl;
	// G4cout << "--//\\\\Scintillator Volume ---- GDML - : " << parser.GetVolume("S_EJ208(1)")->GetSolid()->GetCubicVolume()/(cm3) << G4endl;
	// G4cout <<""<< G4endl;
	Parser = &parser;
     	#ifndef SSRefractionTest
		World = setWorld;
		#endif
		name = "newWorld";
	//--------------
		MPT_UPDATE = -1;
		//WorldBuild(parser, setWorld);
   }

   G4VPhysicalVolume *Construct()
   {
	 WorldBuild((*Parser), World);
	 name = "Constructed";
     return World;
   }

   virtual ~GDMLDetectorConstruction() {};

  public:
    void ConstructSDandField()
    {
		G4SDParticleFilter* photonFilter = new G4SDParticleFilter("photonFilter");
  		photonFilter->add("opticalphoton");
        G4MultiFunctionalDetector* detL = new G4MultiFunctionalDetector("detL");
		G4MultiFunctionalDetector* detR = new G4MultiFunctionalDetector("detR");
		////G4MultiFunctionalDetector* detY = new G4MultiFunctionalDetector("detY");
        G4SDManager::GetSDMpointer()->AddNewDetector(detL);
		G4SDManager::GetSDMpointer()->AddNewDetector(detR);
		////G4SDManager::GetSDMpointer()->AddNewDetector(detY);
        G4VPrimitiveScorer* npho = new G4PSNofCollision("nPho",0);
		npho->SetFilter(photonFilter);
        detL->RegisterPrimitive(npho);
		detR->RegisterPrimitive(npho);

        SetSensitiveDetector("detVOLL",detL);
		SetSensitiveDetector("detVOLR",detR);
		////SetSensitiveDetector("detVOLY",detY);

    }

  private:
	G4GDMLParser* Parser;
	G4GDMLReadStructure* GDMLReader;
    void DefineMaterials();
    G4bool  fCheckOverlaps;
	//G4VUserDetectorConstruction* phantomWorld;
};

#endif
