#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TVector2.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

bool debug = false;

int main(int argc, char** argv) {

  if( argc<3 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: apply corrections ... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " inputfile.root outputfile.root " << std::endl ;
     exit(1) ;
  }

  // input file name
  std::string inputfile((const char*)argv[1]); 
  // output file name
  std::string outputfile((const char*)argv[2]);


  // initialize
  // root files
  TFile* finput = new TFile(inputfile.c_str());
  TFile* foutput = new TFile(outputfile.c_str(), "recreate");


  // tree
  TTree* tree = (TTree*)finput->Get("tree");

  // out_tree
  TTree* tree_out = tree->CloneTree(0);

  // trigger eff
  //TFile* ftrig = TFile::Open("data/dataefficiency.root");
  //TF1* ftrg_z_e = (TF1*)ftrig->Get("dataz_e_func");
  //TF1* ftrg_z_m = (TF1*)ftrig->Get("dataz_m_func");
  //TFile* ftrig = new TFile("data/trigEff.root");
  //TH1D* htrg_z_e = (TH1D*)ftrig->Get("heff_el_zpt");
  //TH1D* htrg_z_m = (TH1D*)ftrig->Get("heff_mu_zpt");
  //TFile* ftrig = new TFile("data/trig_l1l2_80x.root");
  TFile* ftrig = new TFile("data/trig_l1l2_80x_2016b_2p6fbm1.root");
  TH2D* htrg_e_1 = (TH2D*)ftrig->Get("ell1pteta");
  TH2D* htrg_e_2 = (TH2D*)ftrig->Get("ell2pteta");
  TH2D* htrg_m_1 = (TH2D*)ftrig->Get("mul1pteta");
  TH2D* htrg_m_2 = (TH2D*)ftrig->Get("mul2pteta");




  Float_t triggersf;
  Float_t llnunu_l1_pt, llnunu_l1_eta, llnunu_l1_l1_pt, llnunu_l1_l2_pt, llnunu_l1_l1_eta,llnunu_l1_l2_eta;
  Int_t llnunu_l1_l1_pdgId, llnunu_l1_l2_pdgId;
  
  tree->SetBranchAddress("triggersf",&triggersf);
  tree->SetBranchAddress("llnunu_l1_pt",&llnunu_l1_pt);
  tree->SetBranchAddress("llnunu_l1_eta",&llnunu_l1_eta);
  tree->SetBranchAddress("llnunu_l1_l1_pt",&llnunu_l1_l1_pt);
  tree->SetBranchAddress("llnunu_l1_l1_eta",&llnunu_l1_l1_eta);
  tree->SetBranchAddress("llnunu_l1_l1_pdgId",&llnunu_l1_l1_pdgId);
  tree->SetBranchAddress("llnunu_l1_l2_pt",&llnunu_l1_l2_pt);
  tree->SetBranchAddress("llnunu_l1_l2_eta",&llnunu_l1_l2_eta);
  tree->SetBranchAddress("llnunu_l1_l2_pdgId",&llnunu_l1_l2_pdgId);



  int nentries = (int)tree->GetEntries();
  if (debug) nentries = 1000;
  for (int i=0; i<nentries; i++){
    tree->GetEntry(i);

    int l1_pt_bin(-1), l1_eta_bin(-1), l2_pt_bin(-1), l2_eta_bin(-1);

    triggersf = 1.0;
    double eff1=0;
    double eff2=0;

    if (fabs(llnunu_l1_l1_pdgId)==11&&fabs(llnunu_l1_l2_pdgId)==11) {
      if (llnunu_l1_l1_pt<htrg_e_1->GetXaxis()->GetXmin()) {
        l1_pt_bin = 1;
      }
      else if (llnunu_l1_l1_pt>htrg_e_1->GetXaxis()->GetXmin() && llnunu_l1_l1_pt<htrg_e_1->GetXaxis()->GetXmax()){
        l1_pt_bin = htrg_e_1->GetXaxis()->FindBin(llnunu_l1_l1_pt);
      }
      else {
        l1_pt_bin = htrg_e_1->GetXaxis()->GetNbins();
      }

      if (llnunu_l1_l1_eta<htrg_e_1->GetYaxis()->GetXmin()) {
        l1_eta_bin = 1;
      }
      else if (llnunu_l1_l1_eta>htrg_e_1->GetYaxis()->GetXmin() && llnunu_l1_l1_eta<htrg_e_1->GetYaxis()->GetXmax()){
        l1_eta_bin = htrg_e_1->GetYaxis()->FindBin(llnunu_l1_l1_eta);
      }
      else {
        l1_eta_bin = htrg_e_1->GetYaxis()->GetNbins();
      }      

      if (llnunu_l1_l2_pt<htrg_e_2->GetXaxis()->GetXmin()) {
        l2_pt_bin = 1;
      }
      else if (llnunu_l1_l2_pt>htrg_e_2->GetXaxis()->GetXmin() && llnunu_l1_l2_pt<htrg_e_2->GetXaxis()->GetXmax()){
        l2_pt_bin = htrg_e_2->GetXaxis()->FindBin(llnunu_l1_l2_pt);
      }
      else {
        l2_pt_bin = htrg_e_2->GetXaxis()->GetNbins();
      }

      if (llnunu_l1_l2_eta<htrg_e_2->GetYaxis()->GetXmin()) {
        l2_eta_bin = 1;
      }
      else if (llnunu_l1_l2_eta>htrg_e_2->GetYaxis()->GetXmin() && llnunu_l1_l2_eta<htrg_e_2->GetYaxis()->GetXmax()){
        l2_eta_bin = htrg_e_2->GetYaxis()->FindBin(llnunu_l1_l2_eta);
      }
      else {
        l2_eta_bin = htrg_e_2->GetYaxis()->GetNbins();
      }

      if (l1_pt_bin!=-1 && l1_eta_bin!=-1 && l2_pt_bin!=-1 && l2_eta_bin!=-1) {
        eff1 = htrg_e_1->GetBinContent(l1_pt_bin, l1_eta_bin);
        eff2 = htrg_e_2->GetBinContent(l2_pt_bin, l2_eta_bin);
        eff1 /=100.0; eff2 /=100.0;
        //triggersf = eff1 + (1.0-eff1)*eff2;
        triggersf = eff1 ; // only use eff1, becaue eff2 (2nd pass incase 1st does not pass) is shot in statistics, only 3 events out of 222521 in 2016B 2.6fb-1 has this situation, 
      }

    }

    if (fabs(llnunu_l1_l1_pdgId)==13&&fabs(llnunu_l1_l2_pdgId)==13) {
      if (llnunu_l1_l1_pt<htrg_m_1->GetXaxis()->GetXmin()) {
        l1_pt_bin = 1;
      }
      else if (llnunu_l1_l1_pt>htrg_m_1->GetXaxis()->GetXmin() && llnunu_l1_l1_pt<htrg_m_1->GetXaxis()->GetXmax()){
        l1_pt_bin = htrg_m_1->GetXaxis()->FindBin(llnunu_l1_l1_pt);
      }
      else {
        l1_pt_bin = htrg_m_1->GetXaxis()->GetNbins();
      }

      if (llnunu_l1_l1_eta<htrg_m_1->GetYaxis()->GetXmin()) {
        l1_eta_bin = 1;
      }
      else if (llnunu_l1_l1_eta>htrg_m_1->GetYaxis()->GetXmin() && llnunu_l1_l1_eta<htrg_m_1->GetYaxis()->GetXmax()){
        l1_eta_bin = htrg_m_1->GetYaxis()->FindBin(llnunu_l1_l1_eta);
      }
      else {
        l1_eta_bin = htrg_m_1->GetYaxis()->GetNbins();
      }

      if (llnunu_l1_l2_pt<htrg_m_2->GetXaxis()->GetXmin()) {
        l2_pt_bin = 1;
      }
      else if (llnunu_l1_l2_pt>htrg_m_2->GetXaxis()->GetXmin() && llnunu_l1_l2_pt<htrg_m_2->GetXaxis()->GetXmax()){
        l2_pt_bin = htrg_m_2->GetXaxis()->FindBin(llnunu_l1_l2_pt);
      }
      else {
        l2_pt_bin = htrg_m_2->GetXaxis()->GetNbins();
      }

      if (llnunu_l1_l2_eta<htrg_m_2->GetYaxis()->GetXmin()) {
        l2_eta_bin = 1;
      }
      else if (llnunu_l1_l2_eta>htrg_m_2->GetYaxis()->GetXmin() && llnunu_l1_l2_eta<htrg_m_2->GetYaxis()->GetXmax()){
        l2_eta_bin = htrg_m_2->GetYaxis()->FindBin(llnunu_l1_l2_eta);
      }
      else {
        l2_eta_bin = htrg_m_2->GetYaxis()->GetNbins();
      }

      if (l1_pt_bin!=-1 && l1_eta_bin!=-1 && l2_pt_bin!=-1 && l2_eta_bin!=-1) {
        eff1 = htrg_m_1->GetBinContent(l1_pt_bin, l1_eta_bin);
        eff2 = htrg_m_2->GetBinContent(l2_pt_bin, l2_eta_bin);
        eff1 /=100.0; eff2 /=100.0;
        triggersf = eff1 + (1.0-eff1)*eff2;
      }

    }

    if (debug) std::cout << " trig eff = " << triggersf << ", eff1 = " << eff1 << ", eff2 = " << eff2 << std::endl; 

/*
    if (fabs(llnunu_l1_l1_pdgId)==11&&fabs(llnunu_l1_l2_pdgId)==11) {
      //triggersf = ftrg_z_e->Eval(llnunu_l1_pt);
      if (llnunu_l1_pt<htrg_z_e->GetXaxis()->GetXmin()) {
        triggersf=0.0;
      }
      else if (llnunu_l1_pt>htrg_z_e->GetXaxis()->GetXmin() && llnunu_l1_pt<htrg_z_e->GetXaxis()->GetXmax()){
        triggersf = htrg_z_e->GetBinContent(htrg_z_e->GetXaxis()->FindBin(llnunu_l1_pt));
      }
      else {
        triggersf = htrg_z_e->GetBinContent(htrg_z_e->GetXaxis()->GetLast());
      }
    }
    else if (fabs(llnunu_l1_l1_pdgId)==13&&fabs(llnunu_l1_l2_pdgId)==13) {
      //triggersf = ftrg_z_m->Eval(llnunu_l1_pt);
      if (llnunu_l1_pt<htrg_z_m->GetXaxis()->GetXmin()) {
        triggersf=0.0;
      }
      else if (llnunu_l1_pt>htrg_z_m->GetXaxis()->GetXmin() && llnunu_l1_pt<htrg_z_m->GetXaxis()->GetXmax()){
        triggersf = htrg_z_m->GetBinContent(htrg_z_m->GetXaxis()->FindBin(llnunu_l1_pt));
      }
      else {
        triggersf = htrg_z_m->GetBinContent(htrg_z_m->GetXaxis()->GetLast());
      }
    }
*/
    tree_out->Fill();

  }

  foutput->cd();
  tree_out->Write();
  foutput->Close();
  finput->Close();

  return 0;

}



