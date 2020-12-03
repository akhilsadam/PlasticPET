#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GlobalConfig.hh"

#ifndef VERSION
#define VERSION
const int     c = 299792458;
const double  h = 4.135667696*0.000000000000001;
const double  nanop = 0.000000001;
const double  Lmax  = 700;
const double  Lmin  = 350;
const double  EoL   = h*c*eV/(nanop);
const double  length_D  = 100*cm;

// ....oooOO0OOooo........oooOO0OOooo...|SIM-TYPE\...oooOO0OOooo........oooOO0OOooo......

    //default - define CrossSectionTest
    //#define VIEWPORT_ONLY //enable to see a 3d interactive viewport, disable for console run.

    //#define SSSpecularReflectionTest // reflection test - needs SSReflectionTest
    //#define SSReflectionTest // reflection test - needs SingleStrip define and ScintillationDisable
    //#define SSRefractionTest // refraction test - needs SingleStrip define and ScintillationDisable (NOT ReflectionDisable)
    #define CrossSectionTest // prints default cross sections! MAKE SURE THIS IS ENABLED! unless using PVT MPT
    #define LEGEND
// ....oooOO0OOooo........oooOO0OOooo...|MAT-TYPE\...oooOO0OOooo........oooOO0OOooo......
    //#define LYSOTest // swaps scintillator material to LYSO
    //#define PVTTest // swaps scintillator material to PVT
// ....oooOO0OOooo........oooOO0OOooo...|GEOMETRY\...oooOO0OOooo........oooOO0OOooo......
    //#define MultipleStripCell // all 3x16 strips - default

    #define SingleStrip // one EJ208 strip for debugging

// ....oooOO0OOooo........oooOO0OOooo...|PROCESSES\..oooOO0OOooo........oooOO0OOooo......
    //default - define nothing

    //#define ScintillationDisable // define to disable (at present does not really work!!)
    //#define ReflectionDisable

// ....oooOO0OOooo........oooOO0OOooo...|END\........oooOO0OOooo........oooOO0OOooo......
#endif
