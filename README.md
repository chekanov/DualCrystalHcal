# DualCrystalHcal 

This is an example of a single barrel ECAL+HCAL tower.
The HCAl is "sandwitch" calorimeter has abserber made of steel.


The tower consists of 40 layers. Each layer has absorber (steel) with the width of 
1.8 cm. The two active media, positioned after the steel plate, 
are  polystyrene (for scintillation light) and quartz (for Cherenkov light). 
These two layers  are separated by the 0.1 cm of steel. 
Each active media has a width of 0.5 cm. 
The final part of this design is the steel plate of 0.1 cm width that protects the  
quartz layer and forms the wall of the single layer. 
This results in a 2 cm absorber and 1 cm active media for each of the 40 layers.
The interaction length of the 40 x (1.8+0.1+0.1 cm)  steel absorbers  is 
4.77 lambda_I. When adding the active media, the total  interaction length of the module is 5.77 lambda_I.
The transverse size of the module is 20 x 20 cm.

Check the ECAL geometry using this top-level file SimpleECAL.xml that includes "ECalBarrel_DualCrystal.xml":

```bash
geoDisplay DRSingleCrystal.xml 
```

Run this example using the script "A_RUN". It  makes the ROOT file (using Sarah's example). 
To create ROOT files in "output" directory, run "A_RUN_ALL" script (takes very long!)

The main script that runs over data stired in "output" is: 

```bash
A_RUN_RESOLUTION
```

It runs "Resolution.C" which is the main analysis program that fills histograms and put them to "histo" directory.
 
This program is based on the SingleDualCrystal by Sarah Eno (eno@umd.edu).
 

##  Installation 

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
git clone https://github.com/AIDASoft/DD4hep.git
cd DD4hep/examples
git clone git@github.com:chekanov/DualCrystalHcal.git 
# edit CMakeLists.txt and add DualCrystalHcal to
# SET(DD4HEP _EXAMPLES "AlignDet CLICSiD ClientTests Conditions DDCMS DDCodex DDDigi DDG4 DDG4_MySensDet LHeD Optica\
lSurfaces Persistency DDCAD SimpleDetector DualCrystalHcal"
CACHE STRING "List of DD4hep Examples to build")

cd ..
mkdir build
mkdir install
cd build/
cmake -DDD4HEP_USE_GEANT4=ON -DBoost_NO_BOOST_CMAKE=ON -DDD4HEP_USE_LCIO=ON -DBUILD_TESTING=ON -DROOT_DIR=$ROOTSYS -D CMAKE_BUILD_TYPE=Release -DDD4HEP_BUILD_EXAMPLES=ON ..
make -j4
make install
cd ..
source bin/thisdd4hep.sh
```

S.Chekanov (ANL)
