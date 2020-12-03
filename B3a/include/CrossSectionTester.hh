#include "B3aRunAction.hh"
#include "GDMLDetectorConstruction.hh"
#include "B3PrimaryGeneratorAction.hh"
#include "B3aHistoManager.hh"

#include "G4Run.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include <vector>
#include <math.h> 
#include <stdio.h>
using namespace std;
class CrossSectionTester
{

 G4double  fRangeCut[3];
 G4double fEnergyCut[3];
 const G4Run* aRun;
 GDMLDetectorConstruction* aDetector;
 B3PrimaryGeneratorAction* aPrimary;
 double energyP;
 G4AnalysisManager* analysisManager;
string fileOutputName = "lambdas5.txt";
char* fileOutputNameC;
public: void CSRunAction(const G4Run* runI,GDMLDetectorConstruction* fDetectorI,B3PrimaryGeneratorAction* fPrimaryI)
{
  //Clearing crossSection file
    fileOutputNameC = &fileOutputName[0];
    FILE* pfile = fopen(fileOutputNameC, "w");
    fclose(pfile);
  aRun = runI;
  aDetector = fDetectorI;
  aPrimary = fPrimaryI;
  analysisManager = G4AnalysisManager::Instance();
  //CSRun();
  	//set up photon tests
 double en = 1*eV;
 double eM = 1*MeV;
 int maxN = 500;
 double expStep = (log(eM/en)/(1.0*maxN));
 for(int i = 0; i<=maxN; i++)
  { 
	//double deltaL = Lmax - Lmin;
	//double L = (((1.0*i)/maxN)*deltaL) + Lmin;
  	//aPrimary->CSTtest(EoL/(L));
	energyP = (en*exp(expStep*i));//(i*(eM-en)/maxN)+en;
	aPrimary->CSTtest(energyP);
	//G4cout << "Wavelength = " << L << " nm" << G4endl;
	CSRun();
  }
	aPrimary->RESETtest();

}

void CSRun()//void RunAction::BeginOfRunAction(const G4Run*)
{

  //set precision for printing
  G4int prec = G4cout.precision(6);
   //instantate EmCalculator
  G4EmCalculator emCal;
  //  emCal.SetVerbose(2);
     
  // get particle 
  G4ParticleDefinition* particle = aPrimary->GetParticleGun()->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double charge   = particle->GetPDGCharge();    
  G4double energy   = aPrimary->GetParticleGun()->GetParticleEnergy();
 
  // get material
  const G4Material* material = aDetector->GetMaterial();
  G4String matName     = material->GetName();
  G4double density     = material->GetDensity();
  G4double radl        = material->GetRadlen();  
  G4double atmsVolumic  = material->GetTotNbOfAtomsPerVolume();

  /*G4cout << "\n " << partName << " ("
         << G4BestUnit(energy,"Energy") << ") in " 
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ";   radiation length: "
         << G4BestUnit(radl,   "Length")       << ")" << G4endl;*/

  // get cuts         
  GetCuts();
  if (charge != 0.) {
   /*G4cout << "\n  Range cuts : \t gamma "  
                      << std::setw(8) << G4BestUnit(fRangeCut[0],"Length")
          << "\t e- " << std::setw(8) << G4BestUnit(fRangeCut[1],"Length");
   G4cout << "\n Energy cuts : \t gamma " 
                      << std::setw(8) << G4BestUnit(fEnergyCut[0],"Energy")
          << "\t e- " << std::setw(8) << G4BestUnit(fEnergyCut[1],"Energy")
          << G4endl;*/
   }
   
  // max energy transfert
  if (charge != 0.) {
  G4double Mass_c2 = particle->GetPDGMass();
  G4double moverM = electron_mass_c2/Mass_c2;
  G4double gamM1 = energy/Mass_c2, gam = gamM1 + 1., gamP1 = gam + 1.;
  G4double Tmax = energy; 
  if(particle == G4Electron::Electron()) { 
    Tmax *= 0.5; 
  } else if(particle != G4Positron::Positron()) { 
    Tmax = (2*electron_mass_c2*gamM1*gamP1)/(1.+2*gam*moverM+moverM*moverM);
  }
  G4double range = emCal.GetCSDARange(Tmax,G4Electron::Electron(),material);
  
  /*G4cout << "\n  Max_energy _transferable  : " << G4BestUnit(Tmax,"Energy")
         << " (" << G4BestUnit(range,"Length") << ")" << G4endl;*/               
  }
               
  // get processList and extract EM processes (but not MultipleScattering)
  G4ProcessVector* plist = particle->GetProcessManager()->GetProcessList();
  G4String procName;
  G4double cut;
  std::vector<G4String> emName;
  std::vector<G4double> enerCut;
  size_t length = plist->size();
  for (size_t j=0; j<length; j++) {
     procName = (*plist)[j]->GetProcessName();
     cut = fEnergyCut[1];
     if ((procName == "eBrem")||(procName == "muBrems")) cut = fEnergyCut[0];
     if (((*plist)[j]->GetProcessType() == fElectromagnetic) &&
         (procName != "msc")) {
       emName.push_back(procName);
       enerCut.push_back(cut);
     }  
  }
  
  // write html documentation, if requested
  char* htmlDocName = std::getenv("G4PhysListName");   // file name 
  char* htmlDocDir  = std::getenv("G4PhysListDocDir"); // directory
  if (htmlDocName && htmlDocDir) {
    G4LossTableManager::Instance()->DumpHtml();
  }
  
  // print list of processes
  /*G4cout << "\n  processes :                ";
  for (size_t j=0; j<emName.size();j++)
    G4cout << "\t" << std::setw(13) << emName[j] << "\t";
  G4cout << "\t" << std::setw(13) <<"total";*/
  
  //compute cross section per atom (only for single material)
  if (material->GetNumberOfElements() == 1) {
    G4double Z = material->GetZ();
    G4double A = material->GetA();
     
    std::vector<G4double> sigma0;
    G4double sig, sigtot = 0.;

    for (size_t j=0; j<emName.size();j++) {
      sig = emCal.ComputeCrossSectionPerAtom
                      (energy,particle,emName[j],Z,A,enerCut[j]);
      sigtot += sig;                              
      sigma0.push_back(sig);                
    }
    sigma0.push_back(sigtot);

    /*G4cout << "\n \n  cross section per atom    : ";
    for (size_t j=0; j<sigma0.size();j++) {             
      G4cout << "\t" << std::setw(13) << G4BestUnit(sigma0[j], "Surface");
    }
    G4cout << G4endl;*/
  }
    
  //get cross section per volume 
  std::vector<G4double> sigma0;
  std::vector<G4double> sigma1;
  std::vector<G4double> sigma2;
  G4double Sig, SigtotComp = 0., Sigtot = 0.;

  for (size_t j=0; j<emName.size();j++) {
    Sig = emCal.ComputeCrossSectionPerVolume
      (energy,particle,emName[j],material,enerCut[j]);  
    SigtotComp += Sig;    
    sigma0.push_back(Sig);
    Sig = emCal.GetCrossSectionPerVolume(energy,particle,emName[j],material);
    Sigtot += Sig;    
    sigma1.push_back(Sig);
    sigma2.push_back(Sig/density);                        
  }
  sigma0.push_back(SigtotComp);
  sigma1.push_back(Sigtot);
  sigma2.push_back(Sigtot/density);          
    
   //print cross sections

  /*G4cout << "\n \n  compCrossSectionPerVolume : ";
  for (size_t j=0; j<sigma0.size();j++) {             
    G4cout << "\t" << std::setw(13) << sigma0[j]*cm << " cm^-1";
  }*/

  /*G4cout << "\n  cross section per volume : ";
  for (size_t j=0; j<sigma1.size();j++) {             
    G4cout << "\t" << std::setw(13) << sigma1[j]*cm << " cm^-1";
  }*/
  
  /*G4cout << "\n  cross section per mass   : ";
  for (size_t j=0; j<sigma2.size();j++) {
    G4cout << "\t" << std::setw(13) 
           << G4BestUnit(sigma2[j], "Surface/Mass");
  }*/

  //G4cout << "\n  cross section            : ";
  for (size_t j=0; j<sigma2.size();j++) {
	G4double barns = (sigma1[j]/(atmsVolumic))/barn;
    /*G4cout << "\t" << std::setw(13) 
           << barns << " barn ";*/
	if(j==0)
	{
		//analysisManager->FillH1(20,energyP, barns);
	}
	else if(j==1)
	{
		//analysisManager->FillH1(21,energyP, barns);
	}
	else if(j==3)
	{
		//analysisManager->FillH1(22,energyP, barns);
	}
	else if(j==4)
	{
		//analysisManager->FillH1(19,energyP, barns);
	}
  }
   
   
  //print mean free path
  
  G4double lambda;
  
  //G4cout << "\n  mean free path           : ";
  for (size_t j=0; j<sigma1.size();j++) {
    lambda = DBL_MAX; 
    if (sigma1[j] > 0.) lambda = 1/sigma1[j];
    //G4cout << "\t" << std::setw(13) << G4BestUnit( lambda, "Length");
  }
  
  //mean free path (g/cm2)
  //G4cout << "\n        (g/cm2)            : ";  
  std::ofstream myfile5(fileOutputName, std::ios_base::app);
  //std::ofstream myfileb("barns20.txt", std::ios_base::app);
  myfile5 << energyP << " ";  
  for (size_t j=0; j<sigma2.size();j++) {
    lambda =  DBL_MAX;
    if (sigma2[j] > 0.) lambda = 1/sigma2[j];   
    double out = lambda / g * cm *cm;
    myfile5 << out << " ";                    
    //G4cout << "\t" << std::setw(13) << G4BestUnit( lambda, "Mass/Surface"); 
        G4double barns = ((1.0/out)*1.023*pow(10,3)/(6.022*19));
    if(j==0)
	{
		analysisManager->FillH1(20,energyP, barns); //photoE
    //myfileb << energyP << " " << barns << std::endl;    
	}
	else if(j==1)
	{
		analysisManager->FillH1(21,energyP, barns);
	}
	else if(j==3)
	{
		analysisManager->FillH1(22,energyP, barns);
	}
	else if(j==4)
	{
		analysisManager->FillH1(19,energyP, barns);
	}   
  }
  myfile5 << std::endl;
  myfile5.close(); 
  //myfileb.close(); 
  //G4cout << G4endl;
  
  if (charge == 0.) {
    /*G4cout.precision(prec);
    G4cout << "\n-----------------------------------------------------------\n"
           << G4endl;*/
    return;
  }
  
  //get stopping power 
  std::vector<G4double> dedx1;
  std::vector<G4double> dedx2;  
  G4double dedx, dedxtot = 0.;
  size_t nproc = emName.size();

  for (size_t j=0; j<nproc; j++) {
    dedx = emCal.ComputeDEDX(energy,particle,emName[j],material,enerCut[j]);
    dedxtot += dedx;
    dedx1.push_back(dedx);
    dedx2.push_back(dedx/density);                        
  }
  dedx1.push_back(dedxtot);
  dedx2.push_back(dedxtot/density);          
    
  //print stopping power
  /*G4cout << "\n \n  restricted dE/dx         : ";
  for (size_t j=0; j<=nproc; j++) {             
    G4cout << "\t" << std::setw(13) 
           << G4BestUnit(dedx1[j],"Energy/Length");
  }
  
  G4cout << "\n      (MeV/g/cm2)          : ";
  for (size_t j=0; j<=nproc; j++) {
    G4cout << "\t" << std::setw(13) 
           << G4BestUnit(dedx2[j],"Energy*Surface/Mass");
  }*/
  dedxtot = 0.;

  for (size_t j=0; j<nproc; j++) {
    dedx = emCal.ComputeDEDX(energy,particle,emName[j],material,energy);
    dedxtot += dedx;
    dedx1[j] = dedx;
    dedx2[j] = dedx/density;                        
  }
  dedx1[nproc] = dedxtot;
  dedx2[nproc] = dedxtot/density;          
    
  //print stopping power
  /*G4cout << "\n \n  unrestricted dE/dx       : ";
  for (size_t j=0; j<=nproc; j++) {             
    G4cout << "\t" << std::setw(13) << G4BestUnit(dedx1[j],"Energy/Length");
  }
  
  G4cout << "\n      (MeV/g/cm2)          : ";
  for (size_t j=0; j<=nproc; j++) {
    G4cout << "\t" << std::setw(13) 
           << G4BestUnit(dedx2[j],"Energy*Surface/Mass");
  }*/
  
  //get range from restricted dedx
  G4double range1 = emCal.GetRangeFromRestricteDEDX(energy,particle,material);
  G4double range2 = range1*density;

  //print range
  /*G4cout << "\n \n  range from restrict dE/dx: " 
         << "\t" << std::setw(8) << G4BestUnit(range1,"Length")
         << " (" << std::setw(8) << G4BestUnit(range2,"Mass/Surface") << ")";*/
  
  //get range from full dedx
  G4double EmaxTable = G4EmParameters::Instance()->MaxEnergyForCSDARange();
  if(energy < EmaxTable) {
    G4double Range1 = emCal.GetCSDARange(energy,particle,material);
    G4double Range2 = Range1*density;
     
    /*G4cout << "\n  range from full dE/dx    : " 
           << "\t" << std::setw(8) << G4BestUnit(Range1,"Length")
           << " (" << std::setw(8) << G4BestUnit(Range2,"Mass/Surface") << ")";*/
  }

  //get transport mean free path (for multiple scattering)
  G4double MSmfp1 = emCal.GetMeanFreePath(energy,particle,"msc",material);
  G4double MSmfp2 = MSmfp1*density;
  
  //print transport mean free path
  /*G4cout << "\n \n  transport mean free path : " 
         << "\t" << std::setw(8) << G4BestUnit(MSmfp1,"Length")
         << " (" << std::setw(8) << G4BestUnit(MSmfp2,"Mass/Surface") << ")";

  if (particle == G4Electron::Electron()) CriticalEnergy();
           
  G4cout << "\n-------------------------------------------------------------\n";
  G4cout << G4endl;*/
       
 // reset default precision
 G4cout.precision(prec);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ProductionCutsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 void GetCuts()
{  
  G4ProductionCutsTable* theCoupleTable =
        G4ProductionCutsTable::GetProductionCutsTable();
        
  size_t numOfCouples = theCoupleTable->GetTableSize();
  const G4MaterialCutsCouple* couple = 0;
  G4int index = 0;
  for (size_t i=0; i<numOfCouples; i++) {
     couple = theCoupleTable->GetMaterialCutsCouple(i);
     if (couple->GetMaterial() == aDetector->GetMaterial()) {index = i; break;}
  }
  
  fRangeCut[0] =
         (*(theCoupleTable->GetRangeCutsVector(idxG4GammaCut)))[index];
  fRangeCut[1] =      
         (*(theCoupleTable->GetRangeCutsVector(idxG4ElectronCut)))[index];
  fRangeCut[2] =      
         (*(theCoupleTable->GetRangeCutsVector(idxG4PositronCut)))[index]; 

  fEnergyCut[0] =
         (*(theCoupleTable->GetEnergyCutsVector(idxG4GammaCut)))[index];
  fEnergyCut[1] =      
         (*(theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut)))[index];
  fEnergyCut[2] =      
         (*(theCoupleTable->GetEnergyCutsVector(idxG4PositronCut)))[index];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 void CriticalEnergy()
{
  // compute e- critical energy (Rossi definition) and Moliere radius.
  // Review of Particle Physics - Eur. Phys. J. C3 (1998) page 147
  //
  G4EmCalculator emCal;
    
  const G4Material* material = aDetector->GetMaterial();
  const G4double radl = material->GetRadlen();
  G4double ekin = 5*MeV;
  G4double deioni;
  G4double err  = 1., errmax = 0.001;
  G4int    iter = 0 , itermax = 10;  
  while (err > errmax && iter < itermax) {
    iter++;          
    deioni  = radl*
              emCal.ComputeDEDX(ekin,G4Electron::Electron(),"eIoni",material);
    err = std::abs(deioni - ekin)/ekin;
    ekin = deioni;
  }
 /* G4cout << "\n \n  critical energy (Rossi)  : " 
         << "\t" << std::setw(8) << G4BestUnit(ekin,"Energy");*/
         
  //Pdg formula (only for single material)
  G4double pdga[2] = { 610*MeV, 710*MeV };
  G4double pdgb[2] = { 1.24, 0.92 };
  G4double EcPdg = 0.;
  
  if (material->GetNumberOfElements() == 1) {
    G4int istat = 0;
    if (material->GetState() == kStateGas) istat = 1;  
    G4double Zeff = material->GetZ() + pdgb[istat];
    EcPdg = pdga[istat]/Zeff;
   /* G4cout << "\t\t\t (from Pdg formula : " 
           << std::setw(8) << G4BestUnit(EcPdg,"Energy") << ")";    */
  }
     
 const G4double Es = 21.2052*MeV;
 G4double rMolier1 = Es/ekin, rMolier2 = rMolier1*radl;
 /*G4cout << "\n  Moliere radius           : "
        << "\t" << std::setw(8) << rMolier1 << " X0 "   
        << "= " << std::setw(8) << G4BestUnit(rMolier2,"Length");*/
        
 if (material->GetNumberOfElements() == 1) {
    G4double rMPdg = radl*Es/EcPdg;
    /*G4cout << "\t (from Pdg formula : " 
           << std::setw(8) << G4BestUnit(rMPdg,"Length") << ")";*/
  }         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
};
