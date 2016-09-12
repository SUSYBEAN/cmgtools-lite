#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TVector2.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <ctime>


int main(int argc, char** argv) {

  if( argc<3 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: skimming... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " inputfile.root outputfile.root" << std::endl ;
     exit(1) ;
  }

  time_t now = time(0);
  char* dt = ctime(&now);
  std::cout << "Start time is: " << dt << std::endl;

  // input file name
  std::string inputfile((const char*)argv[1]); 
  // output file name
  std::string outputfile((const char*)argv[2]);
  // initialize
  // root files
  TFile* finput = TFile::Open(inputfile.c_str());
  TFile* foutput = TFile::Open(outputfile.c_str(), "recreate");

  // tree
  TTree* tree = (TTree*)finput->Get("tree");

  Int_t isData;
  tree->SetBranchAddress("isData",&isData);
  // check if tree has events
  if (tree->GetEntries()<=0) return 0;
  // get isData info
  tree->GetEntry(0);

  // select branches
  tree->SetBranchStatus("genjet_*",0);
  tree->SetBranchStatus("lep_*",0);
  tree->SetBranchStatus("jet_*",0);
  tree->SetBranchStatus("vtx_*",0);
  //tree->SetBranchStatus("triggers*",0);
  tree->SetBranchStatus("nLL",0); 
  tree->SetBranchStatus("nLLNuNu",0);
  tree->SetBranchStatus("llnunu_deltaPhi",0);
  //tree->SetBranchStatus("llnunu_TuneP_deltaPhi",0);
  tree->SetBranchStatus("llnunu_deltaR",0);
  //tree->SetBranchStatus("llnunu_TuneP_deltaR",0);
  //tree->SetBranchStatus("llnunu_TuneP_LV_*",0);
  //tree->SetBranchStatus("llnunu_LV_*",0);
  //tree->SetBranchStatus("llnunu_CosdphiZMet",0);
  //tree->SetBranchStatus("llnunu_CosZdeltaPhi",0);
  //tree->SetBranchStatus("llnunu_dPTPara",0);
  //tree->SetBranchStatus("llnunu_dPTParaRel",0);
  //tree->SetBranchStatus("llnunu_dPTPerp",0);
  //tree->SetBranchStatus("llnunu_dPTPerpRel",0);
  //tree->SetBranchStatus("llnunu_metOvSqSET",0);
  //tree->SetBranchStatus("llnunu_l1_px",0);
  //tree->SetBranchStatus("llnunu_l1_py",0);
  //tree->SetBranchStatus("llnunu_l1_pz",0);
  //tree->SetBranchStatus("llnunu_l1_TuneP_usage",0);
  //tree->SetBranchStatus("llnunu_l1_l*_TuneP_type",0);
  tree->SetBranchStatus("llnunu_l1_l*_dxy",0);
  tree->SetBranchStatus("llnunu_l1_l*_dz",0);
  tree->SetBranchStatus("llnunu_l1_l*_edxy",0);
  tree->SetBranchStatus("llnunu_l1_l*_edz",0);
  tree->SetBranchStatus("llnunu_l1_l*_ip3d",0);
  //tree->SetBranchStatus("llnunu_l1_l*_miniRelIso",0);
  tree->SetBranchStatus("llnunu_l1_l*_muonincore03",0);
  //tree->SetBranchStatus("llnunu_l1_l*_heepV60_noISO",0);
  tree->SetBranchStatus("llnunu_l1_l*_looseelectron",0);
  tree->SetBranchStatus("llnunu_l1_l*_looseelectronnoiso",0);
  tree->SetBranchStatus("llnunu_l1_l*_looseisoelectron",0);
  //tree->SetBranchStatus("llnunu_l1_l*_lepsf",0);
  //tree->SetBranchStatus("llnunu_l1_l*_lepsfUp",0);
  //tree->SetBranchStatus("llnunu_l1_l*_lepsfLo",0);
  //tree->SetBranchStatus("llnunu_l1_l*_px",0);
  //tree->SetBranchStatus("llnunu_l1_l*_py",0);
  //tree->SetBranchStatus("llnunu_l1_l*_pz",0);
  tree->SetBranchStatus("llnunu_l1_l*_chargedHadRelIso03",0);
  tree->SetBranchStatus("llnunu_l1_l*_chargedHadRelIso04",0);
  tree->SetBranchStatus("llnunu_l1_l*_softMuonId",0);
  tree->SetBranchStatus("llnunu_l1_l*_pfMuonId",0);
  tree->SetBranchStatus("llnunu_l1_l*_eleCutId2012_full5x5",0);
  tree->SetBranchStatus("llnunu_l1_l*_trackerLayers",0);
  tree->SetBranchStatus("llnunu_l1_l*_pixelLayers",0);
  tree->SetBranchStatus("llnunu_l1_l*_trackerHits",0);
  tree->SetBranchStatus("llnunu_l1_l*_lostOuterHits",0);
  tree->SetBranchStatus("llnunu_l1_l*_innerTrackValidHitFraction",0);
  tree->SetBranchStatus("llnunu_l1_l*_innerTrackChi2",0);
  tree->SetBranchStatus("llnunu_l1_l*_nStations",0);
  tree->SetBranchStatus("llnunu_l1_l*_caloCompatibility",0);
  tree->SetBranchStatus("llnunu_l1_l*_globalTrackChi2",0);
  tree->SetBranchStatus("llnunu_l1_l*_trkKink",0);
  tree->SetBranchStatus("llnunu_l1_l*_segmentCompatibility",0);
  tree->SetBranchStatus("llnunu_l1_l*_chi2LocalPosition",0);
  tree->SetBranchStatus("llnunu_l1_l*_chi2LocalMomentum",0);
  tree->SetBranchStatus("llnunu_l1_l*_glbTrackProbability",0);
  tree->SetBranchStatus("llnunu_l1_l*_sigmaIEtaIEta",0);
  tree->SetBranchStatus("llnunu_l1_l*_dEtaScTrkIn",0);
  tree->SetBranchStatus("llnunu_l1_l*_dPhiScTrkIn",0);
  tree->SetBranchStatus("llnunu_l1_l*_hadronicOverEm",0);
  tree->SetBranchStatus("llnunu_l1_l*_eInvMinusPInv",0);
  tree->SetBranchStatus("llnunu_l1_l*_eInvMinusPInv_tkMom",0);
  tree->SetBranchStatus("llnunu_l1_l*_etaSc",0);
  //tree->SetBranchStatus("llnunu_l2_px",0);
  //tree->SetBranchStatus("llnunu_l2_py",0);
  //tree->SetBranchStatus("llnunu_l2_pz",0);
  tree->SetBranchStatus("llnunu_l2_eta",0);
  tree->SetBranchStatus("llnunu_l2_rapidity",0);
  tree->SetBranchStatus("llnunu_l2_mass",0);
  //tree->SetBranchStatus("llnunu_l2_metSig",0); 
  if (!isData) {
    tree->SetBranchStatus("HLT_*",0); 
    tree->SetBranchStatus("genX_*",0); 
    tree->SetBranchStatus("genLep_*",0);
    //tree->SetBranchStatus("genZ_px",0); 
    //tree->SetBranchStatus("genZ_py",0); 
    //tree->SetBranchStatus("genZ_pz",0); 
    tree->SetBranchStatus("llnunu_l2_genEta",0);
  }

  TTree* tree_out = tree->CloneTree(-1);


  foutput->cd();
  tree_out->Write();
  foutput->Close();
  finput->Close();


  now = time(0);
  dt = ctime(&now);
  std::cout << "End time is: " << dt << std::endl;

  return 0;

}



