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

  if( argc<4 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: skimming... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " inputfile.root inputfile_gen.root outputfile.root" << std::endl ;
     exit(1) ;
  }

  time_t now = time(0);
  char* dt = ctime(&now);
  std::cout << "Start time is: " << dt << std::endl;

  // input file name
  std::string inputfile((const char*)argv[1]); 
  // input file gen name
  std::string inputfilegen((const char*)argv[2]); 
  // output file name
  std::string outputfile((const char*)argv[3]);
  // initialize
  // root files
  TFile* finput = TFile::Open(inputfile.c_str());
  TFile* finputgen = TFile::Open(inputfilegen.c_str());
  TFile* foutput = TFile::Open(outputfile.c_str(), "recreate");

  // tree
  TTree* tree = (TTree*)finput->Get("tree");
  TTree* treegen = (TTree*)finputgen->Get("tree");
  treegen->SetName("treegen");

  UInt_t          run;
  UInt_t          lumi;
  ULong64_t       evt;
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("lumi",&lumi);
  tree->SetBranchAddress("evt",&evt);
  treegen->SetBranchAddress("run",&run);
  treegen->SetBranchAddress("lumi",&lumi);
  treegen->SetBranchAddress("evt",&evt);
  
  foutput->cd();
  TTree* tree_out = tree->CloneTree(0);

  // lheNb
  Int_t lheNb, lheNj;
  tree_out->Branch("lheNb", &lheNb, "lheNb/I"); 
  tree_out->Branch("lheNj", &lheNj, "lheNj/I"); 

  Int_t ngenZ; 
  Float_t genZ_pt[100], genZ_px[100], genZ_py[100], genZ_pz[100], genZ_eta[100], genZ_rapidity[100], genZ_phi[100], genZ_energy[100], genZ_mass[100], genZ_mt[100];
  Int_t ngenLep, genLep_pdgId[100], genLep_status[100]; 
  Float_t genLep_pt[100], genLep_px[100], genLep_py[100], genLep_pz[100], genLep_eta[100], genLep_rapidity[100], genLep_phi[100], genLep_mass[100], genLep_charge[100];

  treegen->SetBranchAddress("ngenZ", &ngenZ);
  treegen->SetBranchAddress("genZ_pt", genZ_pt);
  treegen->SetBranchAddress("genZ_px", genZ_px);
  treegen->SetBranchAddress("genZ_py", genZ_py);
  treegen->SetBranchAddress("genZ_pz", genZ_pz);
  treegen->SetBranchAddress("genZ_eta", genZ_eta);
  treegen->SetBranchAddress("genZ_rapidity", genZ_rapidity);
  treegen->SetBranchAddress("genZ_phi", genZ_phi);
  treegen->SetBranchAddress("genZ_energy", genZ_energy);
  treegen->SetBranchAddress("genZ_mass", genZ_mass);
  treegen->SetBranchAddress("genZ_mt", genZ_mt);

  treegen->SetBranchAddress("ngenLep", &ngenLep);
  treegen->SetBranchAddress("genLep_pdgId", genLep_pdgId);
  treegen->SetBranchAddress("genLep_pt", genLep_pt);
  treegen->SetBranchAddress("genLep_px", genLep_px);
  treegen->SetBranchAddress("genLep_py", genLep_py);
  treegen->SetBranchAddress("genLep_pz", genLep_pz);
  treegen->SetBranchAddress("genLep_eta", genLep_eta);
  treegen->SetBranchAddress("genLep_rapidity", genLep_rapidity);
  treegen->SetBranchAddress("genLep_phi", genLep_phi);
  treegen->SetBranchAddress("genLep_mass", genLep_mass);
  treegen->SetBranchAddress("genLep_charge", genLep_charge);
  treegen->SetBranchAddress("genLep_status", genLep_status);

  tree_out->SetBranchAddress("ngenZ", &ngenZ);
  tree_out->SetBranchAddress("genZ_pt", genZ_pt);
  //tree_out->SetBranchAddress("genZ_px", genZ_px);
  //tree_out->SetBranchAddress("genZ_py", genZ_py);
  //tree_out->SetBranchAddress("genZ_pz", genZ_pz);
  tree_out->SetBranchAddress("genZ_eta", genZ_eta);
  tree_out->SetBranchAddress("genZ_rapidity", genZ_rapidity);
  tree_out->SetBranchAddress("genZ_phi", genZ_phi);
  //tree_out->SetBranchAddress("genZ_energy", genZ_energy);
  tree_out->SetBranchAddress("genZ_mass", genZ_mass);
  //tree_out->SetBranchAddress("genZ_mt", genZ_mt);

  //tree_out->SetBranchAddress("ngenLep", &ngenLep);
  //tree_out->SetBranchAddress("genLep_pdgId", genLep_pdgId);
  //tree_out->SetBranchAddress("genLep_pt", genLep_pt);
  //tree_out->SetBranchAddress("genLep_px", genLep_px);
  //tree_out->SetBranchAddress("genLep_py", genLep_py);
  //tree_out->SetBranchAddress("genLep_pz", genLep_pz);
  //tree_out->SetBranchAddress("genLep_eta", genLep_eta);
  //tree_out->SetBranchAddress("genLep_rapidity", genLep_rapidity);
  //tree_out->SetBranchAddress("genLep_phi", genLep_phi);
  //tree_out->SetBranchAddress("genLep_mass", genLep_mass);
  //tree_out->SetBranchAddress("genLep_charge", genLep_charge);
  //tree_out->SetBranchAddress("genLep_status", genLep_status);
 
  Int_t isData;
  tree->SetBranchAddress("isData",&isData);
  // check if tree has events
  if (tree->GetEntries()<=0) return 0;
  // get isData info
  tree->GetEntry(0);

  treegen->BuildIndex("run", "evt");

  for (int i=0; i<(int)tree->GetEntries(); i++){
    tree->GetEntry(i);
    treegen->GetEntryWithIndex(run, evt);
    tree_out->Fill();
  }

  foutput->cd();
  tree_out->Write();
  foutput->Close();
  finput->Close();
  finputgen->Close();


  now = time(0);
  dt = ctime(&now);
  std::cout << "End time is: " << dt << std::endl;

  return 0;

}



