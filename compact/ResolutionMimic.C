/*
* S.Chekanov (ANL)
**/ 
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Factories.h"
#include "DDG4/Geant4Particle.h"
#include "DDG4/Geant4Data.h"
#include "../src/DualCrystalHcalHit.h"

#include <vector>
#include <algorithm>


// way too slow if track all photons for now
// so randomly delete photons after creation according to this fraction
// keep 1 photon out of 1000
// const double supressRandom=0.001;


// information about hit channel IT numbers
const int nchan = 4;
const int ichan[nchan] = {64,73,74,75};  // channel 74 is the crystal, 73 and 75 the two kill media
std::string namechan[nchan] = {"air","PD1","crystal","PD2"};

// MeV to GeV
const float MeV2GeV = 0.001;

// (17*0.5) / (17*4+17*0.5)
// calorimeter sampling fraction
const float sampling_fraction=0.11;

void ResolutionMimic(int num_evtsmax, const char* inputfilename, float beamE, const char* outputfilename ) {


  typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
  typedef std::vector<CalVision::DualCrystalHcalHit*> CalHits;

  // read in libraries that define the classes
  Long_t result;
  char text[1024];
  const char* dd4hep = gSystem->Getenv("DD4hepINSTALL");
  snprintf(text,sizeof(text)," -I%s/include -D__DD4HEP_DDEVE_EXCLUSIVE__ -Wno-shadow -g -O0",dd4hep);
  gSystem->AddIncludePath(text);
  TString fname = "libDDG4IO";
  const char* io_lib = gSystem->FindDynamicLibrary(fname,kTRUE);
  result = gSystem->Load("libDDG4IO");
  result = gSystem->Load("libDDEvePlugins");
  result = gSystem->Load("libDDEvePlugins");
  result = gSystem->Load("libDualCrystalHcal");
  result = gSystem->Load("libDDG4Plugins");

  // ************************ build functions ******************
  // use 10 GeV pions 
  float ptmin=0.0; 
  float ptmax=50.0;
  TF1* f_scint= new TF1("f_scint","[0]+x*[1]+[2]*x*x",ptmin,ptmax);

  float a1=50174.562;
  float a2=1211.215;
  float a3=-4.176;
  if (beamE ==1) { a1=6669.876; a2=-1393.655; a3=672.962; }
  if (beamE ==5) { a1= 24868.266; a2= 1262.575; a3=-4.136; } 
  if (beamE ==10)  { a1=50174.562; a2=1211.215; a3=-4.176; } 
  if (beamE ==20)  { a1=8627.738; a2=5271.548; a3=-45.110; } 

  // only 10 GeV for neutrons
  if (strcmp(outputfilename, "histos_mimic/hist_neutron_10gev.root") == 0) {
    if (beamE ==10)  { a1=26239.941; a2=3177.677; a3=-50.972; }
  } 


  f_scint->SetParameter(0, a1); 
  f_scint->SetParameter(1, a2); 
  f_scint->SetParameter(2, a3); 

  TF1* f_cheren= new TF1("f_cheren","[0]+x*[1]+[2]*x*x",ptmin,ptmax);
  if (beamE ==1) { a1=50174.562; a2=1211.215; a3=-4.176; }
  if (beamE ==5) { a1= 33227.497; a2=13222.813; a3=-150.427; }
  if (beamE ==10)  { a1=71083.780; a2=12466.743; a3=-56.730; }
  if (beamE ==20)  { a1=61621.966; a2=16003.502; a3=-75.346; }

  // only 10 GeV for neutrons
  if (strcmp(outputfilename, "histos_mimic/hist_neutron_10gev.root") == 0) {
    if (beamE ==10)  { a1=30704.768; a2=14588.43; a3=-145.048; }
  }



  f_cheren->SetParameter(0, a1); 
  f_cheren->SetParameter(1, a2); 
  f_cheren->SetParameter(2, a3); 
  //***********************************************************

  // define histograms
  //gen particles
  TH1F *hgenPsize = new TH1F("hgenPsize","number of generator particles",600,0.,40000);
  TH1F *hgenPdgID = new TH1F("hgenpdgID","pdgID of generator particles",600,-200,200);


  // calorimeter infor
  TH1F *hchan = new TH1F("hchan","channel ID number",1028,0.,1028);
  TH1F *hcEcalE = new TH1F("hcEcalE","sum crystal ecal energy",100,0.,100.);
  TH1F *hcEcalncer = new TH1F("hcEcalncer","total number of cerenkov",100,0.,10000);
  TH1F *hcEcalncer0 = new TH1F("hcEcalncer0","total number of cerenkov chan 0",100,0.,10000);
  TH1F *hcEcalncer1 = new TH1F("hcEcalncer1","total number of cerenkov chan 1",100,0.,10000);
  TH1F *hcEcalncer2 = new TH1F("hcEcalncer2","total number of cerenkov chan 2",100,0.,10000);
  TH1F *hcEcalncer3 = new TH1F("hcEcalncer3","total number of cerenkov chan 3",100,0.,10000);

  TH1F *hcEcalE0 = new TH1F("hcEcalE0","energy chan 0",100,0.,10000);
  TH1F *hcEcalE1 = new TH1F("hcEcalE1","energy chan 1",100,0.,10000);
  TH1F *hcEcalE2 = new TH1F("hcEcalE2","energy chan 2",100,0.,10000);
  TH1F *hcEcalE3 = new TH1F("hcEcalE3","energy chan 3",100,0.,10000);

  TH1F *wave_cherenk  = new TH1F("wave_cherenk ","wave of cherenkov",100,0.,1000);
  TH1F *wave_scintil  = new TH1F("wave_scintil ","wave of scinittilation",100,0.,1000);

  TH1F *heest = new TH1F("heest","estimated energy",100,0.0,3.0);
  TH1F *heest_meass = new TH1F("heest_meass","Energy from Cherenkov+Scintillation",100,0.0,3.0);
  TH1F *heest_scint = new TH1F("heest_scint","Energy using Scintillation (e- for calibration)",100,0.0,3.0);
  TH1F *heest_cherenk = new TH1F("heest_cherenk","Energy using Cherenkov (e- for calibration)",100,0.0,3.0);

  TH1F *energyECAL = new TH1F("energyECAL","Energy in ECAL [GeV]",1000,0.0,50.0);
  TH1F *energyHCAL = new TH1F("energyHCAL","Energy in HCAL [GeV]",1000,0.0,50.0);
 
  // max value for multiplicity
  float maxCherenkov=5000000.;
  float maxScintil=5000000.;
  if (beamE ==10)  maxScintil=maxScintil*10;
  if (beamE ==20)  maxScintil=maxScintil*100;

  TH2F *depos_scintil = new TH2F("depos_scintil","NScintilation vs Depos", 1000,0.,100., 1000, 0.0, maxScintil);
  TH2F *depos_cherenk = new TH2F("depos_cherenk","Cherenkov vs Depos",     1000,0.,100., 1000, 0.0,  maxCherenkov);
  TH2F *schint_cherenk = new TH2F("schint_cherenk","Scintil vs Cherenkov", 1000,0., maxScintil, 1000, 0.0, maxCherenkov);

  // after calibration
  TH2F *schint_cherenk_calib = new TH2F("schint_cherenk_calib","Scintil vs Cherenkov (e- calibrated)", 3000,0.,30000., 3000, 0.0, 30000.0);
  TH2F *schint_cherenk_calib_rel = new TH2F("schint_cherenk_calib_rel","Scintil vs Cherenkov (e- calibrated, normalized)", 3000,0.,1.2, 3000, 0.0, 1.2);

  TH1F *scintilVsdepos = new TH1F("scintilPerMeV","L: N of scintillation per MeV",1000,0., 1000);
  TH1F *cherenkVsdepos = new TH1F("cherenkPerMeV","L: N of cherenkov per MeV",1000,0., 1000);
  TH1F *scintil_cherenkVsdepos = new TH1F("scintil_cherenkPerMeV","N of cherenkov per MeV",1000,0., 1000);

  TH1F *scintil= new TH1F("scintil","N of scintillation",50000,0., maxScintil);
  TH1F *cherenk= new TH1F("cherenk","N of cherenkov",50000,0., maxCherenkov);
  TH1F *scintil_cherenk= new TH1F("scintil_cherenk","Combined  scintillation and cherenkov",50000,0., maxCherenkov+maxScintil);
  TH1F *cherenkDIVscintil= new TH1F("cherenkDIVscintil","Nr of cherenkov / scintillation",100, 0., 2.0);

  // ECAL
  TH1F *scintil_ECAL= new TH1F("scintil_ECAL","N of scintillation in ECAL",50000,0., maxScintil);
  TH1F *cherenk_ECAL= new TH1F("cherenk_ECAL","N of cherenkov in ECAL",50000,0., maxCherenkov);
  TH1F *scintil_cherenk_ECAL= new TH1F("scintil_cherenk_ECAL","Combined  scintillation and cherenkov in ECAL",50000,0., maxCherenkov+maxScintil);
  TH1F *cherenkDIVscintil_ECAL= new TH1F("cherenkDIVscintil_ECAL","Nr of cherenkov / scintillation in ECAL",100, 0., 2.0);

  // HCAL
  TH1F *scintil_HCAL= new TH1F("scintil_HCAL","N of scintillation in HCAL",50000,0., maxScintil);
  TH1F *cherenk_HCAL= new TH1F("cherenk_HCAL","N of cherenkov in ECAL",50000,0., maxCherenkov);
  TH1F *scintil_cherenk_HCAL= new TH1F("scintil_cherenk_HCAL","Combined  scintillation and cherenkov in HCAL",50000,0., maxCherenkov+maxScintil);
  TH1F *cherenkDIVscintil_HCAL= new TH1F("cherenkDIVscintil_HCAL","Nr of cherenkov / scintillation in HCAL",100, 0., 2.0);

  // errors are RMS
  //TProfile* profS  = new TProfile("scint_vs_energy","Profile of Scint versus E [GeV]",500,0, 5.0, 0.0, 500.0, "S");
  //TProfile* profC  = new TProfile("cher_vs_energy","Profile of Cheren versus E [GeV]",500,0, 5.0, 0.0, 500.0, "S");

 TProfile* profS  = new TProfile("scint_vs_energy","Profile of Scint versus E [GeV]",500,0, 5.0, 0.0, 500.0);
 TProfile* profC  = new TProfile("cher_vs_energy","Profile of Cheren versus E [GeV]",500,0, 5.0, 0.0, 500.0);
 TH1F *hitsHCAL= new TH1F("hits_e_HCAL","Energies of hits in HCAL",500, 0., 5.0);
 TProfile* profSdis  = new TProfile("scint_vs_distance","Profile of Scint versus D [mm]",100,0, 1500.0, 0.0, 500.0);
 TProfile* profCdis  = new TProfile("cher_vs_distance","Profile of Cheren versus D [mm]",100,0, 1500.0, 0.0, 500.0);

 // before calibration
 TProfile* profSS  = new TProfile("scint_vs_energy_before","Profile of Scint versus E [GeV] no calib",250,0, 5.0, 0.0, 800000.0);
 TProfile* profCC  = new TProfile("cher_vs_energy_before","Profile of Cheren versus E [GeV] no calib",250,0, 5.0, 0.0, 800000.0);


 TProfile2D* profSS_CC  = new TProfile2D("scint_cher_vs_energy_before","Profile of Scint versus E [GeV] no calib",250,0, 5.0, 250,0, 5.0, 0.0, 800000.0);
 TH2F *profSS_CC_2D= new TH2F("scint_cher_vs_energy_before","Profile of Scint versus E [GeV] no calib", 250, 0.0, 800000.0, 250, 0.0, 800000.0);

 // with RMS
 TProfile* profSS_SD  = new TProfile("scint_vs_energy_before_SD","Profile of Scint versus E [GeV] no calib +SD",250,0, 5.0, 0.0, 800000.0, "s");
 TProfile* profCC_SD  = new TProfile("cher_vs_energy_before_SD","Profile of Cheren versus E [GeV] no calib +SD",250,0, 5.0, 0.0, 800000.0, "s");

 // total light vs EM fraction 
 TProfile* profSS_EM  = new TProfile("scint_vs_EM_before_SD","Profile of Scint versus EM counts",200,0, 10000, 0.0, 2000000.0, "s");
 TProfile* profCC_EM  = new TProfile("cher_vs_EM_before_SD","Profile of Cheren versus EM counts",200,0, 10000, 0.0, 2000000.0, "s");
 TProfile* profSS_HAD  = new TProfile("scint_vs_HAD_before_SD","Profile of Scint versus HAD counts",200,0, 10000, 0.0, 2000000.0, "s");
 TProfile* profCC_HAD  = new TProfile("cher_vs_HAD_before_SD","Profile of Cheren versus HAD counts",200,0, 10000, 0.0, 2000000.0, "s");
 TProfile* profSS_EMfrac  = new TProfile("scint_vs_EMfrac_before_SD","Profile of Scint versus EM",200,0, 1.0, 0.0, 20000000.0, "s");
 TProfile* profCC_EMfrac  = new TProfile("cher_vs_EMfrac_before_SD","Profile of Cheren versus EM",200,0, 1.0, 0.0, 20000000.0, "s");

 float fmax=100;
 if (beamE ==1)  fmax=5;
 if (beamE ==5)  fmax=20;
 if (beamE ==10)  fmax=50;
 if (beamE ==20)  fmax=120;

 // vs EM energy
 TProfile* profSS_EM_E  = new TProfile("scint_vs_EM_E_before_SD","Profile of Scint versus EM energy",40,0, fmax, 0.0, 2000000.0, "s");
 TProfile* profCC_EM_E  = new TProfile("cher_vs_EM_E_before_SD","Profile of Cheren versus EM energy",40,0, fmax, 0.0, 2000000.0, "s");


 float ymax=10;
 if (beamE ==1)  ymax=1;
 if (beamE ==5)  ymax=10;
 if (beamE ==10)  ymax=20;
 if (beamE ==20)  ymax=40;

 TH1F *scintil_EM= new TH1F("scintil_shape_EM","N of scintillation in ECAL",50000,0., maxScintil);
 TH1F *cherenk_EM= new TH1F("cheren_shape_EM","N of cherenkov in ECAL",50000,0., maxScintil);


 // vs HD energy
 TProfile* profSS_HAD_E  = new TProfile("scint_vs_HAD_E_before_SD","Profile of Scint versus HAD energy",100,0, 5000, 0.0, 2000000.0, "s");
 TProfile* profCC_HAD_E  = new TProfile("cher_vs_HAD_E_before_SD","Profile of Cheren versus HAD energy",100,0, 5000, 0.0, 2000000.0, "s");

 // vs beta EM 
 TProfile* profSS_BetaEM  = new TProfile("scint_vs_BetaEM_before_SD","Profile of Scint versus EM beta",100,0, 500000, 0.0, 2000000.0, "s");
 TProfile* profCC_BetaEM  = new TProfile("cher_vs_BetaEM_before_SD","Profile of Cheren versus EM beta",100,0, 500000, 0.0, 2000000.0, "s");


  //  const char* inputfilename="/data/users/eno/dd4hep/DD4hep/DDDetectors/compact/testSid.root";
  // const char* outputfilename="hist.root";

  // get Tree
  //  TFile *f = new TFile(inputfilename);
  //f->Print();
  GenParts* pgenparts = new GenParts();
  CalHits* pcalhits = new CalHits();
  int num_evt,nbyte;

  TFile* f = TFile::Open(inputfilename);
  TTree* t = (TTree*)f->Get("EVENT;1");
  t->Print();



  
  // loop over events
  TBranch* b_mc = t->GetBranch("MCParticles");
  TBranch* b_ecal = t->GetBranch("DRCNoSegment");
  int ihaha = b_mc->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<" doing "<<b_mc->GetName()<<std::endl;
  std::cout<<"num_evt gen loop is "<<num_evt<<std::endl;
  
  
  if(num_evt>0) {


    // find branches
    GenParts* gens = new GenParts();
    b_mc->SetAddress(&gens);
    CalHits* ecalhits = new CalHits();
    b_ecal->SetAddress(&ecalhits);


    int SCEPRINT2=10;
    for(int ievt=0;ievt<num_evt; ++ievt) {
      if (ievt%10 == 0) std::cout<<std::endl<<std::endl<<"event number is "<<ievt<<std::endl;

      // gen particles
      nbyte = b_mc->GetEntry(ievt);
      if( nbyte>0) {
	if(ievt<SCEPRINT2) std::cout<<gens->size()<<" Gen particles "<<std::endl;
      }
      hgenPsize->Fill(gens->size());
      for(size_t i=0;i<gens->size(); ++i) {
        dd4hep::sim::Geant4Particle* agen =gens->at(i);
        hgenPdgID->Fill(agen->pdgID);
      }


      int em_n=0;
      int had_n=0;
      float em_E=0;
      float had_E=0;
      float beta_em=0;
      float beta_had=0;

// https://github.com/AIDASoft/DD4hep/blob/master/DDG4/plugins/Geant4EventReaderHepMC.cpp

      for(size_t i=0;i<gens->size(); ++i) {
        dd4hep::sim::Geant4Particle* agen =gens->at(i);
        int pdgID=agen->pdgID;

        float psx=agen->psx;
        float psy=agen->psy;
        float psz=agen->psz;
        float mass=agen->mass;
        float e=sqrt(psx*psx+psy*psy+psz*psz+mass*mass);

        float vsx=agen->vsx;
        float vsy=agen->vsy;
        float vsz=agen->vsz;
        float time=agen->time;

        float beta=0;
        if (time>0)  beta=sqrt(vsx*vsx+vsy*vsy+vsz*vsz) / (time); // distance CM / time (ns) 


        int getStatus=agen->genStatus;
        hgenPdgID->Fill( pdgID );
        //cout << getStatus << endl;
        // leptons or photons
        if ( (abs( pdgID ) >10 && abs( pdgID ) < 20) || abs( pdgID )==22)  {
             em_n++;
             em_E=em_E+e;
             beta_em=beta_em+ beta;
        } else {
           had_n++;
           had_E=had_E+e;
           beta_had=beta_had+ beta;

        };
      }

    em_E= em_E*MeV2GeV;
    had_E=had_E*MeV2GeV;

  // average beta looks the same 
   // beta_em = beta_em / em_n;
   // beta_had = beta_had / had_n;

   cout << "EM="<< em_n << " HAD=" << had_n << " E(EM)=" << em_E << " E(HAD)=" << had_E << " beta(EM)=" << beta_em << " beta(HAD)=" <<  beta_had << endl;



    // ECAL hits  
    // there are hits in the crystal and also the photodetectors "kill media"
    // in the crystal, photons created in the crystal are counted and their wavelengths stored
    // in the photodetector, photons that enter are counted, wavelength stored, and then they are killed


      int nbyte = b_ecal->GetEntry(ievt);
      if( nbyte>0) {
        if(ievt<SCEPRINT2) std::cout<<ecalhits->size()<<" Ecal Hits "<<std::endl;
      }
      float esum=0.;
      float esum_ECAL=0.;
      float esum_HCAL=0.;
      int ncertot_ECAL=0;
      int nscinttot_ECAL=0;
      int ncertot_HCAL=0;
      int nscinttot_HCAL=0;

      float esumchan[nchan]={0.,0.,0.,0.};
      int ncerchan[nchan]={0,0,0,0};
      int nscintchan[nchan]={0,0,0,0};
      int ncertot=0;
      int nscinttot=0;
      int SCEPRINT=10;
      for(size_t i=0;i<ecalhits->size(); ++i) {
	CalVision::DualCrystalHcalHit* aecalhit =ecalhits->at(i);
	//	std::cout<<"       "<<i<<" energy "<<aecalhit->energyDeposit<<std::endl;
	esum+=aecalhit->energyDeposit;
        float hit_gev=(aecalhit->energyDeposit)*MeV2GeV / sampling_fraction ;

        ncertot+=aecalhit->ncerenkov;
        nscinttot+=aecalhit->nscintillator;


        /*
        // use gaussian with RMS 0.1
        double x1 =  f_scint->Eval(  hit_gev );
        double x2 =  f_cheren->Eval(  hit_gev );
        int  nscin = (int)gRandom->Gaus(x1,0.1*x1);
        int  ncher = (int)gRandom->Gaus(x2,0.1*x2);

	ncertot+=ncher; // aecalhit->ncerenkov;
	nscinttot+=nscin; // aecalhit->nscintillator;
	if(i<SCEPRINT&&ievt<SCEPRINT2) std::cout<<" hit channel is "<<aecalhit->cellID<<" in hex is "<< std::hex<< aecalhit->cellID<<std::dec<<" "<<aecalhit->energyDeposit<<" "<<aecalhit->ncerenkov<<" "<<aecalhit->nscintillator<<std::endl;
        */

	// see ../src/DRCrystal_geo.cpp to see the assignments

        /*
        int ihitchan=aecalhit->cellID;
	int idet = (ihitchan & 0xC0)>>6;  // this assignment is made in SCEPCALConstants.xml
	int ilayer = (ihitchan & 0x38)>>3; // this is 1 for crystal and detectors, 0 for air around it
	int islice = (ihitchan & 0x07);  //   this is 1 or 4 for photodetectors, 2 for crystal
        */
        int ihitchan=aecalhit->cellID;
        int idet = (ihitchan) & 0x07;
        int ix = (ihitchan >>3) & 0x1F ;  // is this right?
        if(ix>16) ix=ix-32;
        int iy =(ihitchan >>8) & 0x1F ; // is this right?
        if(iy>16) iy=iy-32;
        int  islice = (ihitchan >>13) & 0x07;
        int  ilayer = (ihitchan>> 16) & 0x07;

        //std::cout<<" ihitchan =" << ihitchan << std::endl; 
        //std::cout<<" idet,ilayer,islice is ("<<idet<<","<<ilayer<<","<<islice<<")"<<std::endl;

        dd4hep::Position  pos = aecalhit->position;
        //std::cout<<"DRcalo deposit "<< " position ("<<pos.X()<<","<<pos.Y()<<","<<pos.Z()<<") "<<std::endl;

          float nscin = (float)(aecalhit->ncerenkov);
          float ncher = (float)(aecalhit->nscintillator);
          //cout <<  nscin << " " << ncher << endl;


          profSS -> Fill(hit_gev,   nscin );
          profCC -> Fill(hit_gev,   ncher );
          // 2D
          profSS_CC -> Fill(hit_gev,  nscin, ncher );
          profSS_CC_2D -> Fill(nscin, ncher );

          profSS_SD -> Fill(hit_gev,   nscin );
          profCC_SD -> Fill(hit_gev,   ncher );

          //cout << "nhits=" << hit_gev << " scin=" << nscin << " cher=" << ncher << endl; 

        // the total size of this detector is 114.2 mm 
        // It is calculated as 40*(2+0.1+0.5+0.1+0.5) + 20+0.2. See: SCEPCAL_DRCrystal.xml 
        // the center is in: 148.6/2 = 743 mm
        // The most left ellement is in -743  
        // The the position of ECAL is -743  to (-743 + 202) = -541, i.e. (-743 -  -541) 

        // ECAL:
        if (pos.Z()>=-743 && pos.Z()<= -541)  {
                esum_ECAL+=aecalhit->energyDeposit;
                ncertot_ECAL+=ncher; 
                nscinttot_ECAL+=nscin; 
        } else { // HCAL 
                esum_HCAL+=aecalhit->energyDeposit; 
                ncertot_HCAL+=ncher; 
                nscinttot_HCAL+=nscin; 
         }

	// channels are 64 air
	//             73 75 detectors
	//            74 crystal
	if(i<SCEPRINT&&ievt<SCEPRINT2) std::cout<<" idet,ilayer,islice is ("<<idet<<","<<ilayer<<","<<islice<<")"<<std::endl;

        //std::cout<<" idet,ilayer,islice is ("<<idet<<","<<ilayer<<","<<islice<<")"<<std::endl;

         
	// print out wavelength spectra
	int ijchan=aecalhit->nbin;
        int bmin = aecalhit->wavelenmin;
        int bmax = aecalhit->wavelenmax;
        float binsize=(bmax-bmin)/float(ijchan); 
 
	for (int j=0;j<ijchan;j++) {
	 //std::cout<<"  ncerwave["<<j<<"]="<<(aecalhit->ncerwave)[j]<<std::endl;
	 //std::cout<<"  nscintwave["<<j<<"]="<<(aecalhit->nscintwave)[j]<<std::endl;
         // let's unpack
         // ibin = (wavelength-hit->wavelenmin)/binsize;
         float xw= (j*binsize)+bmin; 
         wave_scintil->Fill(xw,(aecalhit->nscintwave)[j]);
         wave_cherenk->Fill(xw,(aecalhit->ncerwave)[j]); 
	}
	hchan->Fill(aecalhit->cellID);

/*
        std::cout<<"Total number of channels=" << nchan << endl;
 
      // there is a better way to do this
	int jchan=aecalhit->cellID;
	int kchan=-1;
	for( int i=0;i<nchan;i++ ) {
	  if(ichan[i]==jchan) kchan=i;
	}
	if(kchan==-1) {
	  std::cout<<"unknown hit channel is "<< aecalhit->cellID<<std::endl;
	} else {
	  esumchan[kchan]+=aecalhit->energyDeposit;
	  ncerchan[kchan]+=aecalhit->ncerenkov;
	  nscintchan[kchan]+=aecalhit->nscintillator;
	}
*/

      }  // end loop over hits


// ---------------------------------------
        // draw numbers from functions
        // use gaussian with RMS 0.1
        double x1 =  f_scint->Eval( em_E );
        double x2 =  f_cheren->Eval( em_E  );
        float sigma_scint=0.15;
        float sigma_cher=0.15;
        nscinttot = (int)gRandom->Gaus(x1,sigma_scint*x1);
        ncertot = (int)gRandom->Gaus(x2,sigma_cher*x2);
        //cout << "Created. Scin=" << nscinttot << " cher=" <<  ncertot << " EM_E=" <<  em_E << endl;

//-----------------------------------------------

      float EMfrac=((float)em_n)/(em_n+had_n);

      profSS_EM -> Fill(em_n,  nscinttot );
      profCC_EM -> Fill(em_n,  ncertot );

      profSS_HAD -> Fill(had_n,  nscinttot );
      profCC_HAD -> Fill(had_n,  ncertot );

      profSS_EMfrac -> Fill(EMfrac,  nscinttot );
      profCC_EMfrac -> Fill(EMfrac,  ncertot );

      if (em_E>(ymax-ymax*0.2) && em_E<(ymax+ymax*0.2)) {
            scintil_EM->Fill(nscinttot);
            cherenk_EM->Fill(ncertot);
       }


      profSS_EM_E -> Fill(em_E,  nscinttot );
      profCC_EM_E -> Fill(em_E,  ncertot );

      profSS_HAD_E -> Fill(had_E,  nscinttot );
      profCC_HAD_E -> Fill(had_E,  ncertot );

      // vs beta
      profSS_BetaEM -> Fill( beta_em,  nscinttot );
      profCC_BetaEM -> Fill( beta_em,  ncertot );


      // Does not work as expected giving a smaller fraction of photons
      // convert back to the expected rate since we fill 1 photon from every 1000
      //nscinttot = nscinttot /  supressRandom;
      //ncertot = ncertot / supressRandom;

   
      energyECAL->Fill(esum_ECAL/1000.);
      energyHCAL->Fill(esum_HCAL/1000.);

      // redefine with sampling fraction
      esum = esum_ECAL + (esum_HCAL / sampling_fraction); 


      hcEcalE->Fill(esum/1000.);
      hcEcalncer->Fill(ncertot);
      hcEcalncer0->Fill(ncerchan[0]);
      hcEcalncer1->Fill(ncerchan[1]);
      hcEcalncer2->Fill(ncerchan[2]);
      hcEcalncer3->Fill(ncerchan[3]);
      hcEcalE0->Fill(esumchan[0]);
      hcEcalE1->Fill(esumchan[1]);
      hcEcalE2->Fill(esumchan[2]);
      hcEcalE3->Fill(esumchan[3]);

      // kludge for now
      float mainee=beamE*1000;
      heest->Fill((esum)/mainee);
      //heest->Fill(esum*MeV2GeV);

      // reconstructed from Scintillation
      // this calibration found using e- guns 0.5 GeV
      float calibration_scint= (0.5*1000)/328400;
      if (beamE ==1)  calibration_scint= (1.0*1000)/662700;
      if (beamE ==5)  calibration_scint= (5.0*1000)/3337000;
      if (beamE ==10)  calibration_scint= (10.0*1000)/6677000;
      if (beamE ==20)  calibration_scint= (20.0*1000)/13370000;

      float energy_scint = calibration_scint * nscinttot;
      heest_scint->Fill(energy_scint / mainee);
      //heest_scint->Fill(energy_scint*MeV2GeV);


      // corrected energy for Scinitallation + Cherenkov 
      // https://arxiv.org/pdf/0707.4021.pdf
      // this calibration found using e- guns 0.5 GeV 
      float calibration_cherenkov = (0.5*1000)/24460;
      if (beamE ==1)   calibration_cherenkov= (1.0*1000)/49420;
      if (beamE ==5)   calibration_cherenkov= (5.0*1000)/250300;
      if (beamE ==10)  calibration_cherenkov= (10.0*1000)/500000;
      if (beamE ==20)  calibration_cherenkov= (20.0*1000)/1002000;


      float energy_cherenkov = calibration_cherenkov * ncertot;
      heest_cherenk->Fill(energy_cherenkov / mainee);
      //heest_cherenk->Fill(energy_cherenkov*MeV2GeV);


       // associate calibrated light with raw hits
       // MeV vs Calibrated photons
       for(size_t i=0;i<ecalhits->size(); ++i) {
        CalVision::DualCrystalHcalHit* aecalhit =ecalhits->at(i);
        float hit = aecalhit->energyDeposit;
        // remove 0.1 kEV hits
        if (hit<0.0001) continue; 

        float  ncher = (float)aecalhit->ncerenkov;
        float  nscin = (float)aecalhit->nscintillator;

        // use gaussian with RMS 0.5
        float hit_gev=( hit )*MeV2GeV / sampling_fraction ;
 
        /*
        double x1 =  f_scint->Eval(  hit_gev  );
        double x2 =  f_cheren->Eval(  hit_gev  );
        int  nscin = (int)gRandom->Gaus(x1,0.5*x1);
        int  ncher = (int)gRandom->Gaus(x2,0.5*x2);
        */

        //cout << nscin << " " << ncher << endl;

        dd4hep::Position  pos = aecalhit->position;

        // after SF
        //float hit_gev=(aecalhit->energyDeposit)*MeV2GeV / sampling_fraction ;
         //if (pos.Z()>=-743 && pos.Z()<= -541)  {
         //    continue;
         //} else { // HCAL
          profS -> Fill(hit_gev,   nscin * calibration_scint );
          profC -> Fill(hit_gev,    ncher * calibration_cherenkov );

          // distance in mm
          // It is calculated as 40*(2+0.1+0.5+0.1+0.5) + 20+0.2. See: SCEPCAL_DRCrystal.xml 
          // the center is in: 148.6/2 = 743 mm

          // beginning of HCAL:
          //Xbeg= 148.6/2 -  20.2 = -541;
          
          float  Zdistance = pos.Z() + 541;
          float  distance=TMath::Sqrt(pos.X()*pos.X() + pos.Y()*pos.Y() + Zdistance*Zdistance); 

          //cout << distance << " pos=" <<  pos.Z() << endl;
          profSdis -> Fill(distance,   nscin * calibration_scint );
          profCdis -> Fill(distance,    ncher * calibration_cherenkov );

          hitsHCAL-> Fill(hit_gev);
          //if (hit_gev>0.01) cout << "Hit=" << hit_gev << " Ncher=" << ncher * calibration_cherenkov  << " Nscin=" << nscin * calibration_scint  << endl;
          //cout <<  hit << " " <<  nscin * calibration_scint << endl; 
         //};

        }

       //cinrofS -> Fill(esum/1000.0, calibration_scint); 
       //profC -> Fill(esum/1000.0, calibration_cherenkov); 

       //cout << esum/1000.0 << " " << calibration_scint << endl;

      // scatter plot for calibrated signals
      schint_cherenk_calib->Fill(energy_scint, energy_cherenkov);
      schint_cherenk_calib_rel->Fill(energy_scint/mainee, energy_cherenkov/mainee);
 
      float h_over_e_S=0.7111; // h/e for scinittilation 
      float h_over_e_C=0.5186; // h/e for cherenkov 
      //float kappa=(1-h_over_e_S) / (1-h_over_e_C);
 
      float kappa=0.600; 
      float a1=( energy_scint - kappa * energy_cherenkov);
      float a2= 1 - kappa;

      // overall calibration
      float calibration = 1.0;
      float energy_rec = (calibration)*a1/a2;
      heest_meass->Fill(energy_rec/mainee);
      //heest_meass->Fill(energy_rec*MeV2GeV);

       /* debug 
       std::cout<<" True energy =  "<< mainee  << std::endl;
       std::cout<<" Hit energy deposit =  "<<esum << std::endl;
       std::cout<<" kappa = "<< kappa  << std::endl;
       std::cout<<" Scint energy = "<< energy_scint  << std::endl;
       std::cout<<" Cherenkov energy = "<< energy_cherenkov   << std::endl;
       std::cout<<" a1 = "<< a1   << std::endl;
       std::cout<<" a2 = "<< a2   << std::endl;
       std::cout<<" Cherenkov+Scint energy = "<< energy_rec  << std::endl;
       */

      if(ievt<SCEPRINT2) std::cout<<" total energy deposit "<<esum<< " reco energy =" <<  energy_rec << std::endl;
      float check=0.;
      for( int i=0;i<nchan;i++) {
	if(ievt<SCEPRINT2) std::cout<<"esum ["<<namechan[i]<<"]="<<esumchan[i]<<std::endl;
	check+=esumchan[i];
      }
      if(ievt<SCEPRINT2) std::cout<<" check total energy desposit "<<check<<std::endl;

      if(ievt<SCEPRINT2) std::cout<<" total number of cherenkov is "<<ncertot<<std::endl;
      check=0;
      for( int i=0;i<nchan;i++) {
	if(ievt<SCEPRINT2) std::cout<<"ncerenkov ["<<namechan[i]<<"]="<<ncerchan[i]<<std::endl;
	check+=ncerchan[i];
      }
      if(ievt<SCEPRINT2) std::cout<<" check ncerenkov "<<check<<std::endl;


      if(ievt<SCEPRINT2) std::cout<<" total number of scintillator is "<<nscinttot<<std::endl;
      check=0;
      for( int i=0;i<nchan;i++) {
	if(ievt<SCEPRINT2) std::cout<<"nscintillator ["<<namechan[i]<<"]="<<nscintchan[i]<<std::endl;
	check+=nscintchan[i];
      }
      if(ievt<SCEPRINT2) std::cout<<" check nscintillator "<<check<<std::endl;

       float  esumcrystal = esum;

       // look at correlations
       depos_scintil->Fill(esumcrystal, (float)nscinttot);
       depos_cherenk->Fill(esumcrystal, (float)ncertot);
       schint_cherenk->Fill((float)nscinttot, (float)ncertot);

       // light per MeV. Luminosyty Photons/MeV =300 See: http://scintillator.lbl.gov/
       scintilVsdepos->Fill(nscinttot/ (esumcrystal) );
       cherenkVsdepos->Fill(ncertot/ (esumcrystal) );
       scintil_cherenkVsdepos->Fill( (nscinttot+ncertot)/ (esumcrystal) );

       // just count  
       scintil->Fill(nscinttot);
       cherenk->Fill(ncertot);
       scintil_cherenk->Fill(nscinttot+ncertot);
       if (nscinttot>0) cherenkDIVscintil->Fill(ncertot/(float)nscinttot);


       // just count in ECAL 
       scintil_ECAL->Fill(nscinttot_ECAL);
       cherenk_ECAL->Fill(ncertot_ECAL);
       scintil_cherenk_ECAL->Fill(nscinttot_ECAL+ncertot_ECAL);
       if (nscinttot_ECAL>0) cherenkDIVscintil_ECAL->Fill(ncertot_ECAL/(float)nscinttot_ECAL);
       
       // just count in HCAL 
       scintil_HCAL->Fill(nscinttot_HCAL);
       cherenk_HCAL->Fill(ncertot_HCAL);
       scintil_cherenk_HCAL->Fill(nscinttot_HCAL+ncertot_HCAL);
       if (nscinttot_HCAL>0) cherenkDIVscintil_HCAL->Fill(ncertot_HCAL/(float)nscinttot_HCAL);


    }  //end loop over events
  }  // end if no events
    
  


  
 
 
  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hgenPsize->Write();
  hgenPdgID->Write();
  hcEcalE->Write();
  hcEcalncer->Write();
  hcEcalncer0->Write();
  hcEcalncer1->Write();
  hcEcalncer2->Write();
  hcEcalncer3->Write();

  wave_cherenk->Write();
  wave_scintil->Write();

  // corrected ratios
  heest->Write();
  heest_meass->Write();
  heest_scint->Write();
  heest_cherenk->Write();

  // sergei
  energyECAL->Write();
  energyHCAL->Write();
  depos_scintil->Write();
  depos_cherenk->Write();
  schint_cherenk->Write();
  schint_cherenk_calib->Write();
  schint_cherenk_calib_rel->Write();
  scintilVsdepos->Write();
  cherenkVsdepos->Write();

  scintil->Write();
  cherenk->Write();
  cherenkDIVscintil->Write();
  scintil_cherenk->Write();
  scintil_cherenkVsdepos->Write();

  // ECAL
  scintil_ECAL->Write();
  cherenk_ECAL->Write();
  scintil_cherenk_ECAL->Write();
  cherenkDIVscintil_ECAL->Write();
  // HCAL
  scintil_HCAL->Write();
  cherenk_HCAL->Write();
  scintil_cherenk_HCAL->Write();
  cherenkDIVscintil_HCAL->Write();

  // before calib
  profSS->Write();
  profCC->Write();
  profSS_CC->Write();
  profSS_CC_2D->Write();

  profSS_SD->Write();
  profCC_SD->Write();

  profSS_EM->Write();
  profCC_EM->Write();
  profSS_HAD->Write();
  profCC_HAD->Write();

  profSS_EMfrac->Write();
  profCC_EMfrac->Write();

  profSS_EM_E->Write();
  profCC_EM_E->Write();

  profSS_HAD_E->Write();
  profCC_HAD_E->Write();

  scintil_EM->Write();
  cherenk_EM->Write();

  // vs average beta
  profSS_BetaEM->Write();
  profCC_BetaEM->Write();

  // after calib
  profS->Write();
  profC->Write();
  profSdis->Write();
  profCdis->Write();


  hitsHCAL->Write();
  out->Close();

}



