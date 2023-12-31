##==========================================================================
#  AIDA Detector description implementation 
#--------------------------------------------------------------------------
# Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
# All rights reserved.
#
# For the licensing terms see $DD4hepINSTALL/LICENSE.
# For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
#
#==========================================================================
#
# DualCrystalHcal is a detector description of a homogeneous dualreadout calorimeter
# S.Chekanov (ANL)
#==========================================================================
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

IF(NOT TARGET DD4hep::DDCore)
  find_package ( DD4hep REQUIRED )
  include ( ${DD4hep_DIR}/cmake/DD4hep.cmake )
  include ( ${DD4hep_DIR}/cmake/DD4hepBuild.cmake )
  dd4hep_configure_output()
ENDIF()

dd4hep_set_compiler_flags()
#
#==========================================================================
# CLHEP dependent stuff
if(NOT TARGET CLHEP::CLHEP)
  # ensure we only look for CLHEP once
  find_package(CLHEP QUIET)
endif()
if (TARGET CLHEP::CLHEP)
  dd4hep_print("|++> CLHEP PRESENT. Building DualCrystalHcal examples.")
else()
  dd4hep_print("|++> CLHEP is not present. NOT building DualCrystalHcal examples.")
  return()
endif()
#==========================================================================
if(NOT TARGET Geant4::Interface)
  dd4hep_print("Not creating DualCrystalHcal tests [No Geant4 found]")
  return()
endif()

if (DD4HEP_USE_GEANT4)
  #---------------------------  Plugin library for the simulation framework  ---------
  dd4hep_add_dictionary(G__DualCrystalHcal
    SOURCES ${DD4hep_DIR}/include/ROOT/Warnings.h src/DualCrystalHcalHit.h 
    LINKDEF ${DD4hep_DIR}/include/ROOT/LinkDef.h
    OUTPUT  ${LIBRARY_OUTPUT_PATH}
    USES    DD4hep::DDCore DD4hep::DDG4 Geant4::Interface
    )


#  dd4hep_add_dictionary(G__DRcaloSiPM
#    SOURCES ${DD4hep_DIR}/include/ROOT/Warnings.h src/DRcaloSiPMHit.h 
#    LINKDEF ${DD4hep_DIR}/include/ROOT/LinkDef.h
#    OUTPUT  ${LIBRARY_OUTPUT_PATH}
#    USES    DD4hep::DDCore DD4hep::DDG4 Geant4::Interface
#    )



endif()


#==========================================================================
#
dd4hep_configure_output ()
#
#-------------------------------------------------------------------------------
# we only create only library for DualCrystalHcal. The whole package is a single component
# library. A priory there is no need to seperate the implementation from the
# plugins....
#
# If need arises, this can be changed easily.
#
#---DualCrystalHcal plugin library -------------------------------------------------------
dd4hep_add_plugin ( DualCrystalHcal
    GENERATED G__DualCrystalHcal.cxx 
  SOURCES src/*.cpp src/plugins/*.cpp
  USES    DD4hep::DDCore DD4hep::DDAlign DD4hep::DDCond DD4hep::DDG4 DDDetectors
          ROOT::Core ROOT::Geom ROOT::GenVector CLHEP::CLHEP Geant4::Interface
)





target_include_directories(DualCrystalHcal PUBLIC include)
target_link_options(DualCrystalHcal PRIVATE -L${Geant4_DIR}/..)
#message (STATUS "CLHEP: ${CLHEP}")
#
#
install(TARGETS DualCrystalHcal LIBRARY DESTINATION lib)
install(DIRECTORY include/DualCrystalHcal DESTINATION include)
#
#---Package installation procedure(s) -----------------------------------------
install(DIRECTORY data eve scripts DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/DualCrystalHcal )
#
#
#
