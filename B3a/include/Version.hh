#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4GlobalConfig.hh>
#include <string>
using namespace std;

#ifndef VERSION
#define VERSION
const int     c = 299792458;
const double  h = 4.135667696*0.000000000000001;
const double  nanop = 0.000000001;
const double  Lmax  = 700;
const double  Lmin  = 350;
const double  EoL   = h*c*eV/(nanop);
const double  length_D  = 100*cm; //NOT FULLY IMPLEMENTED YET Z_length
const double  length_X  = 7.74*cm;
const double  length_Y  = 10.32*cm; 
const double  length_Z  = length_D;
const double Dx = length_X;
const double Dy = length_Y;
const double Dz = length_D;
const double Ox = -2.58*cm;
const double Oy = 4.83*cm;
const double Oz = Dz/2;
const double VKT= 0.1*mm;
const double sx = 25.4*mm;
const double sy = 6.2*mm;
const double px = sx + 2*VKT;
const double py = sy + 2*VKT;
const int Nx = 3;
const int Ny = 16;
const int Nz = 1;
const double  att_len = 400*cm;
#ifndef OUTPATH
    #define OUTPATH
    const std::string outpath = "../../data/current/";
#endif

#ifdef TACC
    #define G4TACC
    cout << "--------------------------------------------------------***********************************____________________" <<endl;
    cout << "TACC VERSION IS NOW BEING USED!" <<endl;
    cout <<endl;
    #define G4MULTITHREADED
    #define TACC_CORES
#endif


// ....oooOO0OOooo........oooOO0OOooo...|SIM-TYPE\...oooOO0OOooo........oooOO0OOooo......

    //default - define CrossSectionTest
    //#define VIEWPORT_ONLY //enable to see a 3d interactive viewport, disable for console run.

    //#define DISABLEVK //Disable Vikuiti.
    //#define RoughEnds
    //#define SSSpecularReflectionTest // reflection test - needs SSReflectionTest
    //#define SSReflectionTest // reflection test - needs SingleStrip define and ScintillationDisable
    //#define SSRefractionTest // refraction test - needs SingleStrip define and ScintillationDisable (NOT ReflectionDisable)
    #define CrossSectionTest // prints default cross sections! MAKE SURE THIS IS ENABLED! unless using PVT MPT
    #define ZPredictorTest // plots graphs to find the Z-location of the gamma! Will vary Z-position of Gamma hit!
    #define YPredictorTest // plots graphs to find the (and X) Y-location of the gamma! Will vary Y-position of Gamma hit!
    //#define ReflectionTracking // saves data on photon reflections/absorptions
    #define ElectronPathLength // saves data on electron mean displacement and path length
    //#define ULTRAreflective // makes VK have reflectivity 1!
    //#define LEGEND
// ....oooOO0OOooo........oooOO0OOooo...|MAT-TYPE\...oooOO0OOooo........oooOO0OOooo......
    //#define LYSOTest // swaps scintillator material to LYSO
    //#define PVTTest // swaps scintillator material to PVT
// ....oooOO0OOooo........oooOO0OOooo...|GEOMETRY\...oooOO0OOooo........oooOO0OOooo......
    #define MultipleStripCell // all 3x16 strips - default

    //#define SingleStrip // one EJ208 strip for debugging

// ....oooOO0OOooo........oooOO0OOooo...|PROCESSES\..oooOO0OOooo........oooOO0OOooo......
    //default - define nothing

    //#define ScintillationDisable // define to disable (at present does not really work!!)
    //#define ReflectionDisable

// ....oooOO0OOooo........oooOO0OOooo...|END\........oooOO0OOooo........oooOO0OOooo......
#endif
