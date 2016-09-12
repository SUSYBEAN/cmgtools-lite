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
  tree->SetBranchStatus("photon_*",0);
  tree->SetBranchStatus("jet_*",0);
  tree->SetBranchStatus("gjet_l*_px",0);
  tree->SetBranchStatus("gjet_l*_py",0);
  tree->SetBranchStatus("gjet_l*_pz",0);
  tree->SetBranchStatus("gjet_l2_eta",0);
  tree->SetBranchStatus("gjet_l2_rapidity",0);
  tree->SetBranchStatus("gjet_l2_mass",0);
  tree->SetBranchStatus("gjet_dPT*",0); 
  tree->SetBranchStatus("gjet_deltaPhi",0); 
  tree->SetBranchStatus("gjet_CosDeltaPhi",0);
  tree->SetBranchStatus("gjet_metOvSqSET",0);
  tree->SetBranchStatus("gjet_l1_mass",0);
  tree->SetBranchStatus("gjet_l1_chHadIso",0);
  tree->SetBranchStatus("gjet_l1_phIso",0);
  tree->SetBranchStatus("gjet_l1_neuHadIso",0);
  tree->SetBranchStatus("gjet_l2_metSig",0);
  if (!isData) {
    tree->SetBranchStatus("HLT_*",0); 
    tree->SetBranchStatus("gjet_l1_mcMatchId",0);
    tree->SetBranchStatus("gjet_l1_mcPt",0);
    tree->SetBranchStatus("gjet_l1_mcEta",0);
    tree->SetBranchStatus("gjet_l1_mcPhi",0);
    tree->SetBranchStatus("gjet_l2_genPt",0);
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



