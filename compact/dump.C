#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Factories.h"
#include "DDG4/Geant4Particle.h"
#include "DDG4/Geant4Data.h"
#include "../src/DualCrystalHcalHit.h"

#include <vector>

using namespace CalVision;

void dumpDump() {
  gSystem->Load("libDualCrystalHcal.so");
  gSystem->Load("libDDG4Plugins.so");
  int i = Dump::DualCrystalHcaldumpData(10,"testSCEPCAL.root");

}


