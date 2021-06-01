// by Rylee Ross, Kyle Klein, Ben Fox - adapted for PlasticPET by Akhil Sadam

#include "G4SDManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4GeometryManager.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UserLimits.hh"

#include <iostream>
#include <fstream>
#include <vector>

  class LYSO
  {
    G4Material *fAir;
    G4Material *fWater;
    G4Material *fNaI;
    G4Material *fLYSO;
    G4Material *fTeflon;
    G4Material *fBGO;
    G4Material *fLutetiumOxide;
    G4Material *fYttriumOxide;
    G4Material *fSiliconDioxide;
    G4Material *fGel;
    G4Material *fPolystyrene;
    G4Material *fBariumSulfate;
    G4Material *fTitaniumDioxide;
    G4Material *fResin;
    G4Material *fGlue;
    G4MaterialPropertiesTable *fLYSO_MPT;
    G4MaterialPropertiesTable *fAir_MPT;
    G4MaterialPropertiesTable *fGel_MPT;
    G4MaterialPropertiesTable *fBGO_MPT;
    G4MaterialPropertiesTable *fTeflon_MPT;
    G4MaterialPropertiesTable *fPoly_MPT;
    G4MaterialPropertiesTable *fBariumSulfate_MPT;
    G4MaterialPropertiesTable *fTitaniumDioxide_MPT;
    G4MaterialPropertiesTable *fPlasticSurface_MPT;
    G4MaterialPropertiesTable *fCrystalSurface_MPT;
    G4MaterialPropertiesTable *fBaSO4Surface_MPT;
    G4MaterialPropertiesTable *fTiO2Surface_MPT;
    
    G4OpticalSurface* opCrystalSurface;

    public:
    G4Material* LYSOmat()
    {
        return fLYSO;
    }
    G4OpticalSurface* LYSOsurf()
    {
        return opCrystalSurface;
    }

    LYSO()
    {

        G4NistManager * nist = G4NistManager::Instance();
        G4Element *Hydrogen = nist->FindOrBuildElement("H");
        G4Element *Carbon = nist->FindOrBuildElement("C");
        G4Element *Chlorine = nist->FindOrBuildElement("Cl");
        G4Element *Oxygen = nist->FindOrBuildElement("O");
        G4Element *Lutetium = nist->FindOrBuildElement("Lu");
        G4Element *Yttrium = nist->FindOrBuildElement("Y");
        G4Element *fTeflon = nist->FindOrBuildElement("Teflon");
        G4Element *Na = nist->FindOrBuildElement("Na");
        G4Element *I = nist->FindOrBuildElement("I");
        G4Element *Bismuth = nist->FindOrBuildElement("Bi");
        G4Element *Germanium = nist->FindOrBuildElement("Ge");

        // Air
        fAir = nist->FindOrBuildMaterial("G4_AIR");
        //NaI
        fNaI = new G4Material("SodiumIodine", 3.67*g/cm3, 2);
        fNaI->AddElement(Na,1);
        fNaI->AddElement(I,1);
        // Lu2O3
        fLutetiumOxide = new G4Material("LutetiumOxide", 9.41 * g / cm3, 2);
        fLutetiumOxide->AddElement(Lutetium, 2);
        fLutetiumOxide->AddElement(Oxygen, 3);
        // SiO2
        fSiliconDioxide = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
        fPolystyrene = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
        //BaSO4
        fBariumSulfate = nist->FindOrBuildMaterial("G4_BARIUM_SULFATE");
        //TiO2
        fTitaniumDioxide = nist->FindOrBuildMaterial("G4_TITANIUM_DIOXIDE");
        // Y2O3
        fYttriumOxide = new G4Material("YttriumOxide", 5.01 * g / cm3, 2);
        fYttriumOxide->AddElement(Yttrium, 2);
        fYttriumOxide->AddElement(Oxygen, 3);
        // LYSO
        fLYSO = new G4Material("LYSO", 7.1 * g / cm3, 3);
        fLYSO->AddMaterial(fLutetiumOxide, 81 * perCent);
        fLYSO->AddMaterial(fSiliconDioxide, 14 * perCent);
        fLYSO->AddMaterial(fYttriumOxide, 5 * perCent);
        fBGO = new G4Material("BGO", 7.13 * g/cm3, 3);
            fBGO->AddElement(Bismuth, 4);
            fBGO->AddElement(Germanium, 3);
            fBGO->AddElement(Oxygen, 12);
        fGel = new G4Material("Gel", 7.1 * g / cm3, 3);
        fGel->AddMaterial(fLutetiumOxide, 81 * perCent);
        fGel->AddMaterial(fSiliconDioxide, 14 * perCent);
        fGel->AddMaterial(fYttriumOxide, 5 * perCent);
        fResin = new G4Material("Resin", 1.25 * g/cm3, 4);
        fResin->AddElement(Carbon, 21);
        fResin->AddElement(Hydrogen, 25);
        fResin->AddElement(Chlorine, 1);
        fResin->AddElement(Oxygen, 5);
        fGlue = new G4Material("Glue", 1.25 * g/cm3, 4);
        fGlue->AddElement(Carbon, 21);
        fGlue->AddElement(Hydrogen, 25);
        fGlue->AddElement(Chlorine, 1);
        fGlue->AddElement(Oxygen, 5);

    // Optical Properties of Materials ( nm/1242.38 )
        const G4int num = 120;

        G4double ene[num] = {3.51370304301*eV,3.46284728249*eV,3.4156577087*eV,3.36977603198*eV,3.23635193804*eV,3.22094763375*eV,3.20764136152*eV,3.19637713232*eV,3.18291759837*eV,3.17218298275*eV,3.16123184911*eV,3.1503731311*eV,3.1396057038*eV,3.12928927829*eV,3.12192642432*eV,3.11495009612*eV,3.10797996538*eV,3.10108780463*eV,3.09434823581*eV,3.0876570207*eV,3.0810055618*eV,3.07438540703*eV,3.06780442816*eV,3.06124350731*eV,3.05469721997*eV,3.04817354563*eV,3.04165911826*eV,3.03518039566*eV,3.02874761411*eV,3.02234465979*eV,3.01597393385*eV,3.00960145787*eV,3.00323000843*eV,2.99723759298*eV,2.99057154424*eV,2.984211427*eV,2.97752386539*eV,2.9682041556*eV,2.95859140503*eV,2.95063652972*eV,2.93947102052*eV,2.92513342795*eV,2.89963655932*eV,2.87253927733*eV,2.85352309132*eV,2.83628914618*eV,2.8222957219*eV,2.80993218523*eV,2.79766640437*eV,2.78549390188*eV,2.77344008928*eV,2.76149998282*eV,2.74965574426*eV,2.73791911725*eV,2.72625989445*eV,2.71468687979*eV,2.70321484598*eV,2.62319693516*eV,2.61249533729*eV,2.60188652089*eV,2.59136351654*eV,2.58094246644*eV,2.5706077365*eV,2.56035262437*eV,2.55018739616*eV,2.53884632056*eV,2.46805717364*eV,2.45631868318*eV,2.44469337969*eV,2.43204791693*eV,2.41841286976*eV,2.40494311619*eV,2.38943563283*eV,2.37197339274*eV,2.35478359973*eV,2.33681395534*eV,2.31707351876*eV,2.30377129545*eV,2.25609178885*eV,2.23488272047*eV,2.21410318259*eV,2.19371928099*eV,2.17372057346*eV,2.15408900424*eV,2.13481455788*eV,2.11588408116*eV,2.09729119825*eV,2.07902560898*eV,2.0610754277*eV,2.04343516532*eV,2.02609686477*eV,2.00905157727*eV,1.99229317822*eV,1.97581265517*eV,1.95960675492*eV,1.94367103348*eV,1.92800809024*eV,1.91260128415*eV,1.89744776279*eV,1.88252859884*eV,1.86783894415*eV,1.85337461802*eV,1.83912783245*eV,1.82509840454*eV,1.81128037215*eV,1.7976639421*eV,1.78425469025*eV,1.77643183324*eV,2.64723791997*eV,2.65431419592*eV,2.66836164719*eV,2.6790662172*eV,2.6863139081*eV,2.69716329233*eV,2.53364876322*eV,2.5209805194*eV,2.51147200842*eV,2.49905148118*eV,2.48976188644*eV,2.48661723357*eV};
        G4double BGOatt[num] = {1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm, 1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm, 1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm,1.0*cm};
        G4double fast[num] = {0.013820781627794876, 0.013938358896580905, 0.014050591744058316, 0.01770083149773294, 0.1361291993771121, 0.17951253934849465, 0.22274311794414414, 0.2652312673463606, 0.3114997040874239, 0.35139224620750076, 0.40066157551326187, 0.4516999082721216, 0.5045072444840795, 0.5579031796294753, 0.6034774651162591, 0.6472816782656827, 0.688432386235458, 0.7345952706556798, 0.7938782640197157, 0.8552250947457001, 0.9177512611070833, 0.9805722613773165, 1.0445725972829487, 1.1076884314620314, 1.1693300960968656, 1.2303820929140001, 1.2893702523691863, 1.3492429135509216, 1.4111794120946055, 1.4734107445471392, 1.5362317448173723, 1.5958095720902579, 1.6524390602746462, 1.7103057819919416, 1.766140287506697, 1.8148092601521415, 1.8614449477107844, 1.910241473878124, 1.9568572980042371, 1.9966865978054973, 2.040506844218846, 2.084160879129502, 2.1189188575511695, 2.0948333321877324, 2.053620718004467, 2.012992427217593, 1.9736802001779261, 1.9334781269904011, 1.8919493012130524, 1.848651471982605, 1.807122646205256, 1.7669205730177313, 1.7258339981036572, 1.6856319249161322, 1.642334095685685, 1.597267263002139, 1.5526426811818679, 1.2800838749773036, 1.2372282966101307, 1.1952572199695073, 1.153286143328884, 1.1139685718679084, 1.0750932512702074, 1.035775679809232, 0.9977848609380802, 0.9548461768195611, 0.6689313945441024, 0.6314819655514974, 0.594386337249512, 0.556352584860516, 0.5169089741303492, 0.4798240346709801, 0.4399488619201559, 0.3975782897867255, 0.3587457245594927, 0.32193320213070775, 0.2849053292144217, 0.2621546621464601, 0.2066071763487929, 0.17922479084668286, 0.1590792376527035, 0.1416675988862397, 0.12715069304302862, 0.11392033516570743, 0.1019765252542757, 0.09051517083005267, 0.08017954587598375, 0.07064801340059512, 0.061116480925207384, 0.05222822243276415, 0.043983237923266305, 0.03605989040524049, 0.028779816870159625, 0.02166056183081544, 0.015667036261624467, 0.011442514145532101, 0.01156009141431813, 0.013285853640466527, 0.017584711798394714, 0.020757840486169243, 0.02296605819952635, 0.024531001929938068, 0.024648579198724097, 0.024766156467510125, 0.02456209674482368, 0.022428215073302393, 0.021580881367671, 0.023419362297779855, 1.3330963870244557, 1.3692671524065998, 1.4099353409976352, 1.4506171995680286, 1.486787964950173, 1.5274698235205666, 0.9082334285675229, 0.8675652399764877, 0.8223586182384865, 0.7862151928150589, 0.750058097412273, 0.727447951553593};
        G4double rLyso[num] = {1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82};
        // G4double rGel[num] = {1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82};
        G4double rBaSO4[num] = {1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636, 1.636};
        G4double rPoly[num] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
        G4double rAir[num] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        // G4double rBGO[num] = {2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15};
        G4double rGel[num] = {1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55};
        G4double val = 1.67;
        G4double rBGO[num] = {val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val, val};
        G4double eneWithGamma[num+1] = { 511000 * eV, 3.51370304301*eV,3.46284728249*eV,3.4156577087*eV,3.36977603198*eV,3.23635193804*eV,3.22094763375*eV,3.20764136152*eV,3.19637713232*eV,3.18291759837*eV,3.17218298275*eV,3.16123184911*eV,3.1503731311*eV,3.1396057038*eV,3.12928927829*eV,3.12192642432*eV,3.11495009612*eV,3.10797996538*eV,3.10108780463*eV,3.09434823581*eV,3.0876570207*eV,3.0810055618*eV,3.07438540703*eV,3.06780442816*eV,3.06124350731*eV,3.05469721997*eV,3.04817354563*eV,3.04165911826*eV,3.03518039566*eV,3.02874761411*eV,3.02234465979*eV,3.01597393385*eV,3.00960145787*eV,3.00323000843*eV,2.99723759298*eV,2.99057154424*eV,2.984211427*eV,2.97752386539*eV,2.9682041556*eV,2.95859140503*eV,2.95063652972*eV,2.93947102052*eV,2.92513342795*eV,2.89963655932*eV,2.87253927733*eV,2.85352309132*eV,2.83628914618*eV,2.8222957219*eV,2.80993218523*eV,2.79766640437*eV,2.78549390188*eV,2.77344008928*eV,2.76149998282*eV,2.74965574426*eV,2.73791911725*eV,2.72625989445*eV,2.71468687979*eV,2.70321484598*eV,2.62319693516*eV,2.61249533729*eV,2.60188652089*eV,2.59136351654*eV,2.58094246644*eV,2.5706077365*eV,2.56035262437*eV,2.55018739616*eV,2.53884632056*eV,2.46805717364*eV,2.45631868318*eV,2.44469337969*eV,2.43204791693*eV,2.41841286976*eV,2.40494311619*eV,2.38943563283*eV,2.37197339274*eV,2.35478359973*eV,2.33681395534*eV,2.31707351876*eV,2.30377129545*eV,2.25609178885*eV,2.23488272047*eV,2.21410318259*eV,2.19371928099*eV,2.17372057346*eV,2.15408900424*eV,2.13481455788*eV,2.11588408116*eV,2.09729119825*eV,2.07902560898*eV,2.0610754277*eV,2.04343516532*eV,2.02609686477*eV,2.00905157727*eV,1.99229317822*eV,1.97581265517*eV,1.95960675492*eV,1.94367103348*eV,1.92800809024*eV,1.91260128415*eV,1.89744776279*eV,1.88252859884*eV,1.86783894415*eV,1.85337461802*eV,1.83912783245*eV,1.82509840454*eV,1.81128037215*eV,1.7976639421*eV,1.78425469025*eV,1.77643183324*eV,2.64723791997*eV,2.65431419592*eV,2.66836164719*eV,2.6790662172*eV,2.6863139081*eV,2.69716329233*eV,2.53364876322*eV,2.5209805194*eV,2.51147200842*eV,2.49905148118*eV,2.48976188644*eV,2.48661723357*eV};
        G4double abs[num+1] = {1.2 * cm, 20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm};
        G4double reflecVikuiti[num] = {0.9167, 0.9167, 0.9167, 0.9378, 0.9412, 0.9412, 0.9467, 0.9467, 0.9467, 0.9467, 0.9467, 0.9467, 0.9467, 0.9507, 0.9507, 0.9507, 0.9507, 0.9507, 0.9507, 0.9507, 0.9507, 0.9507, 0.9507, 0.9524, 0.9524, 0.9524, 0.9524, 0.9524, 0.9524, 0.9524, 0.9524, 0.9524, 0.9524, 0.9524, 0.9524, 0.9548, 0.9548, 0.9548, 0.9548, 0.9548, 0.9548, 0.9548, 0.9581, 0.9581, 0.9581, 0.9605, 0.9605, 0.9605, 0.9605, 0.9605, 0.9605, 0.9605, 0.9605, 0.9605, 0.9605, 0.9628, 0.9628, 0.964, 0.964, 0.9647, 0.9647, 0.9647, 0.9647, 0.9647, 0.9659, 0.9659, 0.9661, 0.9661, 0.9665, 0.9665, 0.9665, 0.9676, 0.9676, 0.9676, 0.9676, 0.9676, 0.9674, 0.9674, 0.9677, 0.9677, 0.9678, 0.9683, 0.9683, 0.9684, 0.9684, 0.9683, 0.9683, 0.9686, 0.9686, 0.9683, 0.9683, 0.9681, 0.9681, 0.969, 0.969, 0.9695, 0.9695, 0.9693, 0.9693, 0.969, 0.969, 0.9692, 0.9692, 0.9685, 0.9685, 0.9683, 0.9683, 0.9682, 0.964, 0.964, 0.9628, 0.9628, 0.9628, 0.9628, 0.9659, 0.9659, 0.9659, 0.9661, 0.9661, 0.9661};
        //G4double reflecVikuiti[num] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        //G4double reflecTiO2[num] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        G4double reflecBaSO4[num] = {0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955};
        G4double rTiO2[num] = {1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61};
        G4double reflecTiO2[num] = {0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955};
        G4double transmittanceTiO2[num] = {0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045};

        fLYSO_MPT = new G4MaterialPropertiesTable();
        fAir_MPT = new G4MaterialPropertiesTable();
        fGel_MPT = new G4MaterialPropertiesTable();
        fBGO_MPT = new G4MaterialPropertiesTable();
        fPoly_MPT = new G4MaterialPropertiesTable();
        fBariumSulfate_MPT = new G4MaterialPropertiesTable();
        fTitaniumDioxide_MPT = new G4MaterialPropertiesTable();


        fLYSO_MPT->AddProperty("FASTCOMPONENT", ene, fast, num)->SetSpline(true);
        fLYSO_MPT->AddProperty("RINDEX", ene, rLyso, num);
        fLYSO_MPT->AddProperty("ABSLENGTH", eneWithGamma, abs, num+1);
        fLYSO_MPT->AddConstProperty("SCINTILLATIONYIELD", 33200 / MeV);
        fLYSO_MPT->AddConstProperty("RESOLUTIONSCALE", 1);
        fLYSO_MPT->AddConstProperty("FASTTIMECONSTANT", 43 * ns);

        fAir_MPT->AddProperty("RINDEX", ene, rAir, num);

        fGel_MPT->AddProperty("RINDEX", ene, rGel, num);

        // fBGO_MPT->AddProperty("FASTCOMPONENT", ene, fast, num)->SetSpline(true);
        fBGO_MPT->AddProperty("RINDEX", ene,rBGO,num);
        // fBGO_MPT->AddProperty("ABSLENGTH", ene, BGOatt, num);
        // fBGO_MPT->AddConstProperty("SCINTILLATIONYIELD", 9000 / MeV);
        // fBGO_MPT->AddConstProperty("RESOLUTIONSCALE", 1);
        // fBGO_MPT->AddConstProperty("FASTTIMECONSTANT", 300 * ns);


        fLYSO->SetMaterialPropertiesTable(fLYSO_MPT);
        fAir->SetMaterialPropertiesTable(fAir_MPT);
        fGel->SetMaterialPropertiesTable(fGel_MPT);
        fBGO->SetMaterialPropertiesTable(fBGO_MPT);
        fPolystyrene->SetMaterialPropertiesTable(fPoly_MPT);
        fBariumSulfate->SetMaterialPropertiesTable(fBariumSulfate_MPT);
        fTitaniumDioxide->SetMaterialPropertiesTable(fTitaniumDioxide_MPT);
        //fResin->SetMaterialPropertiesTable(fResin_MPT);
        //fGlue->SetMaterialPropertiesTable(fGlue_MPT);

        //SURFACES
        const G4int n = 6;
        G4double PhotonEnergy[n] = {3.105*eV,2.95714*eV,2.855*eV,2.7*eV,2.5875*eV,2.388*eV};
        G4double specularlobe[n] = {1.0,1.0,1.0,1.0,1.0,1.0};
        G4double specularspike[n] = {0.0,0.0,0.0,0.0,0.0,0.0};
        G4double backscatter[n] = {0.0,0.0,0.0,0.0,0.0,0.0};
        //Defining the surface properties of crystals. Crystals may be ground or polished, so choose the correct sigma alpha.
        G4double sigma_alphaGround = 0.209439; //12 deg., for ground crystal
        G4double sigma_alphaPolished = 0.0226893; //1.3 deg., for polished crystal
        //G4double sigma_alphaTest = 0.209439 * radian;
        
        //G4double crystalTransmission[n] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
        opCrystalSurface = new G4OpticalSurface("CrystalSurface");
        fCrystalSurface_MPT = new G4MaterialPropertiesTable();
        fCrystalSurface_MPT->AddProperty("SPECULARLOBECONSTANT",PhotonEnergy,specularlobe,n);
        fCrystalSurface_MPT->AddProperty("SPECULARSPIKECONSTANT",PhotonEnergy,specularspike,n);
        fCrystalSurface_MPT->AddProperty("BACKSCATTERCONSTANT",PhotonEnergy,backscatter,n);
        //fCrystalSurface_MPT->AddProperty("TRANSMITTANCE", PhotonEnergy, crystalTransmission, n);
        opCrystalSurface->SetFinish(ground);
        opCrystalSurface->SetModel(unified);
        opCrystalSurface->SetSigmaAlpha(sigma_alphaGround);
        opCrystalSurface->SetMaterialPropertiesTable(fCrystalSurface_MPT);
        //ONLY set crystal skin surface if using Vikuiti/air coupled separator, to set sigmaalpha at crystal-air interface.
        //G4LogicalSkinSurface* crystalSurface = new G4LogicalSkinSurface("Crystal Surface", logicCryst, opCrystalSurface);

        
        //Defining surface properties of detectors
        G4OpticalSurface* opDetectorSurface = new G4OpticalSurface("DetectorSurface");
        opDetectorSurface->SetFinish(ground);  
        const G4int effNum = 64;
        G4double ephoton_pmt[effNum] = {4.39937628481263* eV, 4.26793273171382* eV, 4.25415407177944* eV, 4.2526101654458* eV, 4.25002775896025* eV, 4.23740768267528* eV, 4.23740768267528* eV, 4.22619410336501* eV, 4.21111365866894* eV, 4.19468373277856* eV, 4.17562764465* eV, 4.16485707439122* eV, 4.13287618164404* eV, 4.11182706958868* eV, 4.06068564113166* eV, 3.98990143108327* eV, 3.91255690517576* eV, 3.82905214632875* eV, 3.73170831479369* eV, 3.63100753721902* eV, 3.54335760187188* eV, 3.39028291586928* eV, 3.313575740457* eV, 3.26400828279202* eV, 3.2073708649341* eV, 3.13717401566535* eV, 3.06569708729327* eV, 2.99026243581681* eV, 2.92677376304993* eV, 2.81453688978407* eV, 2.71059029663409* eV, 2.61404815600077* eV, 2.51541948296409* eV, 2.47766764433903* eV, 2.43449980055764* eV, 2.38261019441473* eV, 2.31642767655566* eV, 2.28160772652005* eV, 2.24045139601668* eV, 2.20375721015054* eV, 2.16824560915714* eV, 2.13386033457783* eV, 2.09415257211902* eV, 2.06881296944741* eV, 2.03552435424357* eV, 1.99418369175086* eV, 1.9602891257364* eV, 1.92752749589036* eV, 1.89584292993575* eV, 1.8621095224178* eV, 1.82796546757295* eV, 1.78873181823744* eV, 1.75185829972651* eV, 1.71101762636774* eV, 1.67203778810094* eV, 1.63479444340984* eV, 1.60473676359121* eV, 1.56380972550422* eV, 1.53439796153312* eV, 1.50385115860456* eV, 1.47365578936701* eV, 1.44464911982904* eV, 1.416762310119* eV, 1.39353045806523};
        G4double reflectivity_pmt[effNum] = {0.02696598, 0.05893489, 0.07987711, 0.11722792, 0.1014316, 0.12730853, 0.14009744, 0.15859531, 0.18073614, 0.2021246, 0.22482906, 0.24095285, 0.26211995, 0.28273989, 0.30246012, 0.32418864, 0.34063229, 0.35529757, 0.36923624, 0.38282257, 0.39408207, 0.40300594, 0.41655242, 0.4362877, 0.45028576, 0.4658789, 0.48091868, 0.49178681, 0.4997318, 0.50860914, 0.50921274, 0.49967056, 0.47483176, 0.45819545, 0.44008796, 0.42196203, 0.39788722, 0.38184955, 0.36710202, 0.35260655, 0.33865954, 0.32471254, 0.3101095, 0.29717251, 0.28064034, 0.26568978, 0.25062753, 0.23575333, 0.22153729, 0.20680864, 0.19418899, 0.17932864, 0.16402218, 0.14939387, 0.13560323, 0.1227101, 0.11222387, 0.10103327, 0.08902995, 0.07927751, 0.06943588, 0.06025241, 0.05136811, 0.0440509};
        G4MaterialPropertiesTable *detectorMPT = new G4MaterialPropertiesTable();
        //Currently using reflectivity to account for detector efficiency, because as of yet I can't make efficiency work. Photons that are "reflected" off the detector are really just counted as detector hits, and they're killed when they reflect. All other photons hitting the 
        //detector are "absorbed" in Geant4, which for us means that they are lost due to detector efficiency.
        detectorMPT->AddProperty("REFLECTIVITY", ephoton_pmt, reflectivity_pmt, effNum);
        //detectorMPT->AddProperty("EFFICIENCY",   ephoton_pmt, efficiency_pmt,   effNum);
        opDetectorSurface->SetMaterialPropertiesTable(detectorMPT);
        //G4LogicalSkinSurface* DetectorSideSurface = new G4LogicalSkinSurface("Detector Surface", logicDet_outer, opDetectorSurface);


        //Defining surface properties for Vikuiti, BaSO4, and TiO2.
        //Ascribing Vikuiti properties
        fPlasticSurface_MPT = new G4MaterialPropertiesTable();
        //Defining probabilities that the Vikuity surface undergoes certain reflection processes. Since using unified model, should be 100% specular lobe, with sigma_alpha specified.
        fPlasticSurface_MPT->AddProperty("SPECULARLOBECONSTANT",PhotonEnergy,specularlobe,n);
        fPlasticSurface_MPT->AddProperty("SPECULARSPIKECONSTANT",PhotonEnergy,specularspike,n);
        fPlasticSurface_MPT->AddProperty("BACKSCATTERCONSTANT",PhotonEnergy,backscatter,n);
        fPlasticSurface_MPT->AddProperty("REFLECTIVITY", ene, reflecVikuiti, num);
        G4OpticalSurface* opPlasticSurface = new G4OpticalSurface("PlasticSurface");
        opPlasticSurface->SetModel(unified);
        opPlasticSurface->SetType(dielectric_metal);
        opPlasticSurface->SetFinish(ground);
        //opPlasticSurface->SetSigmaAlpha(sigma_alphaPolished);
        opPlasticSurface->SetMaterialPropertiesTable(fPlasticSurface_MPT);
        //G4LogicalSkinSurface* PlasticSurface_T = new G4LogicalSkinSurface("Plastic Surface", logicPolyTop, opPlasticSurface);
        //G4LogicalSkinSurface* PlasticSurface_S = new G4LogicalSkinSurface("Plastic Surface", logicPolySide, opPlasticSurface);
        
        //Ascribing BaSO4 properties
        fBaSO4Surface_MPT = new G4MaterialPropertiesTable();
        //Defining probabilities that the surface undergoes certain reflection processes. Since using unified model, should be 100% specular lobe, with sigma_alpha specified.
        fBaSO4Surface_MPT->AddProperty("SPECULARLOBECONSTANT",PhotonEnergy,specularlobe,n);
        fBaSO4Surface_MPT->AddProperty("SPECULARSPIKECONSTANT",PhotonEnergy,specularspike,n);
        fBaSO4Surface_MPT->AddProperty("BACKSCATTERCONSTANT",PhotonEnergy,backscatter,n);
        fBaSO4Surface_MPT->AddProperty("RINDEX", ene, rBaSO4, num);
        fBaSO4Surface_MPT->AddProperty("REFLECTIVITY", ene, reflecBaSO4, num);
        G4OpticalSurface* opBaSO4Surface = new G4OpticalSurface("BaSO4Surface");
        opBaSO4Surface->SetModel(unified);
        opBaSO4Surface->SetType(dielectric_dielectric);
        opBaSO4Surface->SetFinish(groundbackpainted);
        opBaSO4Surface->SetSigmaAlpha(sigma_alphaGround);
        opBaSO4Surface->SetMaterialPropertiesTable(fBaSO4Surface_MPT);
        //G4LogicalSkinSurface* BaSO4Surface_T = new G4LogicalSkinSurface("BaSO4 Surface", logicBaSO4Top, opBaSO4Surface);
        //G4LogicalSkinSurface* BaSO4Surface_S = new G4LogicalSkinSurface("BaSO4 Surface", logicBaSO4Side, opBaSO4Surface);
        //G4LogicalBorderSurface* BaSO4Surface_T  = new G4LogicalBorderSurface("BaSO4 Top Surface", logicCryst, logicBaSO4Top, opBaSO4Surface);
        //G4LogicalBorderSurface* BaSO4Surface_S  = new G4LogicalBorderSurface("BaSO4 Side Surface", logicCryst, logicBaSO4Side, pBaSO4Surface);
        
        //Ascribing TiO2 properties
        fTiO2Surface_MPT = new G4MaterialPropertiesTable();
        //Defining probabilities that the surface undergoes certain reflection processes. Since using unified model, should be 100% specular lobe, with sigma_alpha specified.
        fTiO2Surface_MPT->AddProperty("SPECULARLOBECONSTANT",PhotonEnergy,specularlobe,n);
        fTiO2Surface_MPT->AddProperty("SPECULARSPIKECONSTANT",PhotonEnergy,specularspike,n);
        fTiO2Surface_MPT->AddProperty("BACKSCATTERCONSTANT",PhotonEnergy,backscatter,n);
        fTiO2Surface_MPT->AddProperty("RINDEX", ene, rTiO2, num);
        fTiO2Surface_MPT->AddProperty("REFLECTIVITY", ene, reflecTiO2, num);
        //fTiO2Surface_MPT->AddProperty("TRANSMITTANCE", ene, transmittanceTiO2, num);
        G4OpticalSurface* opTiO2Surface = new G4OpticalSurface("TiO2Surface");
        opTiO2Surface->SetModel(unified);
        opTiO2Surface->SetType(dielectric_dielectric);
        opTiO2Surface->SetFinish(groundbackpainted);
        opTiO2Surface->SetSigmaAlpha(sigma_alphaGround);
        opTiO2Surface->SetMaterialPropertiesTable(fTiO2Surface_MPT);
        //G4LogicalSkinSurface* TiO2Surface_T = new G4LogicalSkinSurface("TiO2 Surface", logicTiO2Top, opTiO2Surface);
        //G4LogicalSkinSurface* TiO2Surface_S = new G4LogicalSkinSurface("TiO2 Surface", logicTiO2Side, opTiO2Surface);
    }

  };
