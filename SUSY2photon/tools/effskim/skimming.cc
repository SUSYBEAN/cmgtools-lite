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

bool lightWeight=true;

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

  // light weight skim
  if (lightWeight){
    /*
    tree->SetBranchStatus("lep_*",0);
    tree->SetBranchStatus("jet_*",0);
    tree->SetBranchStatus("llnunu_LV_*",0);
    tree->SetBranchStatus("llnunu_TuneP_LV_*",0);
    tree->SetBranchStatus("triggers*",0);
    tree->SetBranchStatus("llnunu_l1_l1_leps*",0);
    tree->SetBranchStatus("llnunu_l1_l2_leps*",0);
    */
    if(tree->FindBranch("lep_*")) tree->SetBranchStatus("lep_*",0);
    if(tree->FindBranch("jet_*")) tree->SetBranchStatus("jet_*",0);
    if(tree->FindBranch("llnunu_LV_*")) tree->SetBranchStatus("llnunu_LV_*",0);
    if(tree->FindBranch("llnunu_TuneP_LV_*")) tree->SetBranchStatus("llnunu_TuneP_LV_*",0);
    if(tree->FindBranch("triggers*")) tree->SetBranchStatus("triggers*",0);
    if(tree->FindBranch("llnunu_l1_l2_leps*")) tree->SetBranchStatus("llnunu_l1_l1_leps*",0);
    if(tree->FindBranch("llnunu_l1_l2_leps*")) tree->SetBranchStatus("llnunu_l1_l2_leps*",0);
  }

  if(tree->FindBranch("trgsf"))  tree->SetBranchStatus("trgsf",0);
  if(tree->FindBranch("idsf"))  tree->SetBranchStatus("idsf",0);
  if(tree->FindBranch("isosf"))  tree->SetBranchStatus("isosf",0);
  if(tree->FindBranch("trksf"))  tree->SetBranchStatus("trksf",0);
  if(tree->FindBranch("trgsf_err"))  tree->SetBranchStatus("trgsf_err",0);
  if(tree->FindBranch("idsf_err"))  tree->SetBranchStatus("idsf_err",0);
  if(tree->FindBranch("isosf_err"))  tree->SetBranchStatus("isosf_err",0);
  if(tree->FindBranch("trksf_err"))  tree->SetBranchStatus("trksf_err",0);
  // out_tree
  TTree* tree_out = tree->CloneTree(0);

  //TFile* ftkhp = TFile::Open("all80x12p9.root");
  TFile* ftkhp = TFile::Open("muon80x12p9.root");
  TH2F* tkhpdt=(TH2F*)ftkhp->Get("eff_trackHighPt_80Xdata_pteta");
  TH2F* tkhpmc=(TH2F*)ftkhp->Get("eff_trackHighPt_80Xmc_pteta");
  TH2F* hpdt=(TH2F*)ftkhp->Get("eff_HighPt_80Xdata_pteta");
  TH2F* hpmc=(TH2F*)ftkhp->Get("eff_HighPt_80Xmc_pteta");
  TH2F* isopteta=(TH2F*)ftkhp->Get("sf_trackerIso_80X_pteta");


  TFile* ftrg = TFile::Open("trigereff12p9.root");
  TH2D* mul1pteta=(TH2D*)ftrg->Get("mul1pteta");
  TH2D* mul2pteta=(TH2D*)ftrg->Get("mul2pteta");
  TH2D* ell1pteta=(TH2D*)ftrg->Get("ell1pteta");

  TFile* ftrg_mu = TFile::Open("trigeff_mu.root");
  TH2D* htrg_l1_tot = (TH2D*)ftrg_mu->Get("htrg_l1_tot");
  TH2D* htrg_l2_tot = (TH2D*)ftrg_mu->Get("htrg_l2_tot");
  TH2D* htrg_l1_l1p = (TH2D*)ftrg_mu->Get("htrg_l1_l1p");
  TH2D* htrg_l2_l1p = (TH2D*)ftrg_mu->Get("htrg_l2_l1p");
  TH2D* htrg_l1_l1f = (TH2D*)ftrg_mu->Get("htrg_l1_l1f");
  TH2D* htrg_l2_l1f = (TH2D*)ftrg_mu->Get("htrg_l2_l1f");
  TH2D* htrg_l1_l1pl2f = (TH2D*)ftrg_mu->Get("htrg_l1_l1pl2f");
  TH2D* htrg_l1_l1pl2p = (TH2D*)ftrg_mu->Get("htrg_l1_l1pl2p");
  TH2D* htrg_l1_l1fl2p = (TH2D*)ftrg_mu->Get("htrg_l1_l1fl2p");
  TH2D* htrg_l2_l1pl2f = (TH2D*)ftrg_mu->Get("htrg_l2_l1pl2f");
  TH2D* htrg_l2_l1pl2p = (TH2D*)ftrg_mu->Get("htrg_l2_l1pl2p");
  TH2D* htrg_l2_l1fl2p = (TH2D*)ftrg_mu->Get("htrg_l2_l1fl2p");

  Double_t Ntrg_tot;
  Double_t Ntrg_tot_err;
  Double_t Ntrg_l1pl2f;
  Double_t Ntrg_l1pl2f_err;
  Double_t Ntrg_l1pl2p;
  Double_t Ntrg_l1pl2p_err;
  Double_t Ntrg_l1fl2p;
  Double_t Ntrg_l1fl2p_err;
  Double_t Ntrg_l1p;
  Double_t Ntrg_l1p_err;
  Double_t Ntrg_l1f;
  Double_t Ntrg_l1f_err;

  Int_t trg_NPtBins = htrg_l2_tot->GetNbinsX();
  Int_t trg_NEtaBins = htrg_l2_tot->GetNbinsY();
  Ntrg_tot = htrg_l2_tot->IntegralAndError(Int_t(1), trg_NPtBins, Int_t(1), trg_NEtaBins, Ntrg_tot_err);
  Ntrg_l1p = htrg_l2_l1p->IntegralAndError(Int_t(1), trg_NPtBins, Int_t(1), trg_NEtaBins, Ntrg_l1p_err);
  Ntrg_l1f = htrg_l2_l1f->IntegralAndError(Int_t(1), trg_NPtBins, Int_t(1), trg_NEtaBins, Ntrg_l1f_err);
  Ntrg_l1pl2f = htrg_l2_l1pl2f->IntegralAndError(Int_t(1), trg_NPtBins, Int_t(1), trg_NEtaBins, Ntrg_l1pl2f_err);
  Ntrg_l1pl2p = htrg_l2_l1pl2p->IntegralAndError(Int_t(1), trg_NPtBins, Int_t(1), trg_NEtaBins, Ntrg_l1pl2p_err);
  Ntrg_l1fl2p = htrg_l2_l1fl2p->IntegralAndError(Int_t(1), trg_NPtBins, Int_t(1), trg_NEtaBins, Ntrg_l1fl2p_err);

  TH2D* htrg_l1_tot_norm = (TH2D*)htrg_l1_tot->Clone("htrg_l1_tot_norm");
  TH2D* htrg_l2_tot_norm = (TH2D*)htrg_l2_tot->Clone("htrg_l2_tot_norm");
  TH2D* htrg_l1_l1p_norm = (TH2D*)htrg_l1_l1p->Clone("htrg_l1_l1p_norm");
  TH2D* htrg_l1_l1f_norm = (TH2D*)htrg_l1_l1f->Clone("htrg_l1_l1f_norm");
  TH2D* htrg_l2_l1p_norm = (TH2D*)htrg_l2_l1p->Clone("htrg_l2_l1p_norm");
  TH2D* htrg_l2_l1f_norm = (TH2D*)htrg_l2_l1f->Clone("htrg_l2_l1f_norm");
  TH2D* htrg_l1_l1pl2f_norm = (TH2D*)htrg_l1_l1pl2f->Clone("htrg_l1_l1pl2f_norm");
  TH2D* htrg_l1_l1pl2p_norm = (TH2D*)htrg_l1_l1pl2p->Clone("htrg_l1_l1pl2p_norm");
  TH2D* htrg_l1_l1fl2p_norm = (TH2D*)htrg_l1_l1fl2p->Clone("htrg_l1_l1fl2p_norm");
  TH2D* htrg_l2_l1pl2f_norm = (TH2D*)htrg_l2_l1pl2f->Clone("htrg_l2_l1pl2f_norm");
  TH2D* htrg_l2_l1pl2p_norm = (TH2D*)htrg_l2_l1pl2p->Clone("htrg_l2_l1pl2p_norm");
  TH2D* htrg_l2_l1fl2p_norm = (TH2D*)htrg_l2_l1fl2p->Clone("htrg_l2_l1fl2p_norm");

  htrg_l1_tot_norm->Scale(1./Ntrg_tot);
  htrg_l2_tot_norm->Scale(1./Ntrg_tot);
  htrg_l1_l1p_norm->Scale(1./Ntrg_l1p);
  htrg_l1_l1f_norm->Scale(1./Ntrg_l1f);
  htrg_l2_l1p_norm->Scale(1./Ntrg_l1p);
  htrg_l2_l1f_norm->Scale(1./Ntrg_l1f);
  htrg_l1_l1pl2f_norm->Scale(1./Ntrg_l1pl2f);
  htrg_l1_l1pl2p_norm->Scale(1./Ntrg_l1pl2p);
  htrg_l1_l1fl2p_norm->Scale(1./Ntrg_l1fl2p);
  htrg_l2_l1pl2f_norm->Scale(1./Ntrg_l1pl2f);
  htrg_l2_l1pl2p_norm->Scale(1./Ntrg_l1pl2p);
  htrg_l2_l1fl2p_norm->Scale(1./Ntrg_l1fl2p);

  TH2D* htrg_l1_l1p_norm_vs_tot    = (TH2D*)htrg_l1_l1p_norm->Clone("htrg_l1_l1p_norm_vs_tot");
  TH2D* htrg_l1_l1f_norm_vs_tot    = (TH2D*)htrg_l1_l1f_norm->Clone("htrg_l1_l1f_norm_vs_tot");
  TH2D* htrg_l2_l1p_norm_vs_tot    = (TH2D*)htrg_l2_l1p_norm->Clone("htrg_l2_l1p_norm_vs_tot");
  TH2D* htrg_l2_l1f_norm_vs_tot    = (TH2D*)htrg_l2_l1f_norm->Clone("htrg_l2_l1f_norm_vs_tot");
  TH2D* htrg_l1_l1pl2f_norm_vs_tot = (TH2D*)htrg_l1_l1pl2f_norm->Clone("htrg_l1_l1pl2f_norm_vs_tot");
  TH2D* htrg_l1_l1pl2p_norm_vs_tot = (TH2D*)htrg_l1_l1pl2p_norm->Clone("htrg_l1_l1pl2p_norm_vs_tot");
  TH2D* htrg_l1_l1fl2p_norm_vs_tot = (TH2D*)htrg_l1_l1fl2p_norm->Clone("htrg_l1_l1fl2p_norm_vs_tot");
  TH2D* htrg_l2_l1pl2f_norm_vs_tot = (TH2D*)htrg_l2_l1pl2f_norm->Clone("htrg_l2_l1pl2f_norm_vs_tot");
  TH2D* htrg_l2_l1pl2p_norm_vs_tot = (TH2D*)htrg_l2_l1pl2p_norm->Clone("htrg_l2_l1pl2p_norm_vs_tot");
  TH2D* htrg_l2_l1fl2p_norm_vs_tot = (TH2D*)htrg_l2_l1fl2p_norm->Clone("htrg_l2_l1fl2p_norm_vs_tot");
  TH2D* htrg_l1_l1pl2f_norm_vs_l1p = (TH2D*)htrg_l1_l1pl2f_norm->Clone("htrg_l1_l1pl2f_norm_vs_l1p");
  TH2D* htrg_l1_l1pl2p_norm_vs_l1p = (TH2D*)htrg_l1_l1pl2p_norm->Clone("htrg_l1_l1pl2p_norm_vs_l1p");
  TH2D* htrg_l1_l1fl2p_norm_vs_l1f = (TH2D*)htrg_l1_l1fl2p_norm->Clone("htrg_l1_l1fl2p_norm_vs_l1f");
  TH2D* htrg_l2_l1pl2f_norm_vs_l1p = (TH2D*)htrg_l2_l1pl2f_norm->Clone("htrg_l2_l1pl2f_norm_vs_l1p");
  TH2D* htrg_l2_l1pl2p_norm_vs_l1p = (TH2D*)htrg_l2_l1pl2p_norm->Clone("htrg_l2_l1pl2p_norm_vs_l1p");
  TH2D* htrg_l2_l1fl2p_norm_vs_l1f = (TH2D*)htrg_l2_l1fl2p_norm->Clone("htrg_l2_l1fl2p_norm_vs_l1f");

  htrg_l1_l1p_norm_vs_tot->Divide(htrg_l1_tot_norm);
  htrg_l1_l1f_norm_vs_tot->Divide(htrg_l1_tot_norm);
  htrg_l2_l1p_norm_vs_tot->Divide(htrg_l2_tot_norm);
  htrg_l2_l1f_norm_vs_tot->Divide(htrg_l2_tot_norm);
  htrg_l1_l1pl2f_norm_vs_tot->Divide(htrg_l1_tot_norm);
  htrg_l1_l1pl2p_norm_vs_tot->Divide(htrg_l1_tot_norm);
  htrg_l1_l1fl2p_norm_vs_tot->Divide(htrg_l1_tot_norm);
  htrg_l2_l1pl2f_norm_vs_tot->Divide(htrg_l2_tot_norm);
  htrg_l2_l1pl2p_norm_vs_tot->Divide(htrg_l2_tot_norm);
  htrg_l2_l1fl2p_norm_vs_tot->Divide(htrg_l2_tot_norm);
  htrg_l1_l1pl2f_norm_vs_l1p->Divide(htrg_l1_l1p_norm);
  htrg_l1_l1pl2p_norm_vs_l1p->Divide(htrg_l1_l1p_norm);
  htrg_l1_l1fl2p_norm_vs_l1f->Divide(htrg_l1_l1f_norm);
  htrg_l2_l1pl2f_norm_vs_l1p->Divide(htrg_l2_l1p_norm);
  htrg_l2_l1pl2p_norm_vs_l1p->Divide(htrg_l2_l1p_norm);
  htrg_l2_l1fl2p_norm_vs_l1f->Divide(htrg_l2_l1f_norm);



  //TFile* ftrg = new TFile("trgeff_L12p9.root");
  //TH2D* mul1pteta = (TH2D*)ftrg->Get("hr_mu_l1_eta_pt");
  //TH2D* ell1pteta = (TH2D*)ftrg->Get("hr_el_l1_eta_pt");
 
  //TFile* loosE = new TFile("Loose_80X_2ndPeriod.txt_SF2D.root");
  TFile* loosE = TFile::Open("egammaEffi.txt_SF2D.root");
  TH2F* esfh2=(TH2F*)loosE->Get("EGamma_SF2D");

  TFile* fEtk = TFile::Open("egammatracking.root");
  TH2F* etksf=(TH2F*)fEtk->Get("EGamma_SF2D");

  TFile* fMtk = TFile::Open("muontrackingsf.root");
  TH1F* mtksf=(TH1F*)fMtk->Get("hist_ratio_eta");

  Double_t effdt1a,effdt2a,effmc1a,effmc2a,errdt1a,errdt2a,errmc1a,errmc2a,effdt1,effdt2,effmc1,effmc2,errdt1,errdt2,errmc1,errmc2,trgsfall,idsfall,isosfall,trgsfallerr,idsfallerr,isosfallerr,trksfall,trksfallerr,temp1,temp2,tkhighptdt1,tkhighptmc1,highptdt1,highptmc1,tkhighptdt2,tkhighptmc2,highptdt2,highptmc2;
  Double_t trgsfall_up,trgsfall_dn;
  TBranch *b_trgsfall=tree_out->Branch("trgsf",&trgsfall,"trgsf/D");
  TBranch *b_isosfall=tree_out->Branch("isosf",&isosfall,"isosf/D");
  TBranch *b_idsfall=tree_out->Branch("idsf",&idsfall,"idsf/D");
  TBranch *b_trksfall=tree_out->Branch("trksf",&trksfall,"trksf/D");
  TBranch *b_trgsfallerr=tree_out->Branch("trgsf_err",&trgsfallerr,"trgsf_err/D");
  TBranch *b_trgsfall_up=tree_out->Branch("trgsf_up",&trgsfall_up,"trgsf_up/D");
  TBranch *b_trgsfall_dn=tree_out->Branch("trgsf_dn",&trgsfall_dn,"trgsf_dn/D");
  TBranch *b_isosfallerr=tree_out->Branch("isosf_err",&isosfallerr,"isosf_err/D");
  TBranch *b_idsfallerr=tree_out->Branch("idsf_err",&idsfallerr,"idsf_err/D");
  TBranch *b_tkhighpt_dt_eff_m1=tree_out->Branch("tkhighpt_dt_eff_m1",&tkhighptdt1,"tkhighpt_dt_eff_m1/D");
  TBranch *b_tkhighpt_mc_eff_m1=tree_out->Branch("tkhighpt_mc_eff_m1",&tkhighptmc1,"tkhighpt_mc_eff_m1/D");
  TBranch *b_highpt_dt_eff_m1=tree_out->Branch("highpt_dt_eff_m1",&highptdt1,"highpt_dt_eff_m1/D");
  TBranch *b_highpt_mc_eff_m1=tree_out->Branch("highpt_mc_eff_m1",&highptmc1,"highpt_mc_eff_m1/D");
  TBranch *b_tkhighpt_dt_eff_m2=tree_out->Branch("tkhighpt_dt_eff_m2",&tkhighptdt2,"tkhighpt_dt_eff_m2/D");
  TBranch *b_tkhighpt_mc_eff_m2=tree_out->Branch("tkhighpt_mc_eff_m2",&tkhighptmc2,"tkhighpt_mc_eff_m2/D");
  TBranch *b_highpt_dt_eff_m2=tree_out->Branch("highpt_dt_eff_m2",&highptdt2,"highpt_dt_eff_m2/D");
  TBranch *b_highpt_mc_eff_m2=tree_out->Branch("highpt_mc_eff_m2",&highptmc2,"highpt_mc_eff_m2/D");
  TBranch *b_trksfallerr=tree_out->Branch("trksf_err",&trksfallerr,"trksf_err/D");

  Float_t llnunu_l1_l1_phi, llnunu_l1_l1_eta, llnunu_l1_l2_phi, llnunu_l1_l2_eta,llnunu_l1_l1_pt,llnunu_l1_l2_pt,llnunu_l1_l1_eSCeta,llnunu_l1_l2_eSCeta;
  Int_t llnunu_l1_l1_pdgId, llnunu_l1_l2_pdgId;
  tree->SetBranchAddress("llnunu_l1_l1_pdgId",&llnunu_l1_l1_pdgId);
  tree->SetBranchAddress("llnunu_l1_l2_pdgId",&llnunu_l1_l2_pdgId);
  tree->SetBranchAddress("llnunu_l1_l1_phi",&llnunu_l1_l1_phi);
  tree->SetBranchAddress("llnunu_l1_l2_phi",&llnunu_l1_l2_phi);
  tree->SetBranchAddress("llnunu_l1_l1_pt",&llnunu_l1_l1_pt);
  tree->SetBranchAddress("llnunu_l1_l2_pt",&llnunu_l1_l2_pt);
  tree->SetBranchAddress("llnunu_l1_l1_eta",&llnunu_l1_l1_eta);
  tree->SetBranchAddress("llnunu_l1_l2_eta",&llnunu_l1_l2_eta);
  tree->SetBranchAddress("llnunu_l1_l1_eSCeta",&llnunu_l1_l1_eSCeta);
  tree->SetBranchAddress("llnunu_l1_l2_eSCeta",&llnunu_l1_l2_eSCeta);

  for (int i=0; i<(int)tree->GetEntries(); i++){
    if (i%10000==0) std::cout << " Events " << i << std::endl;
    tree->GetEntry(i);
    if(abs(llnunu_l1_l1_pdgId)==13&&abs(llnunu_l1_l2_pdgId)==13){
      tkhighptdt1=tkhpdt->GetBinContent(tkhpdt->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      tkhighptmc1=tkhpmc->GetBinContent(tkhpmc->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      highptdt1=hpdt->GetBinContent(hpdt->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      highptmc1=hpmc->GetBinContent(hpmc->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      tkhighptdt2=tkhpdt->GetBinContent(tkhpdt->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      tkhighptmc2=tkhpmc->GetBinContent(tkhpmc->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      highptdt2=hpdt->GetBinContent(hpdt->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      highptmc2=hpmc->GetBinContent(hpmc->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      effdt1=tkhpdt->GetBinContent(tkhpdt->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      effmc1=tkhpmc->GetBinContent(tkhpmc->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      errdt1=tkhpdt->GetBinError(tkhpdt->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      errmc1=tkhpmc->GetBinError(tkhpmc->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      effdt1a=hpdt->GetBinContent(hpdt->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      effmc1a=hpmc->GetBinContent(hpmc->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      errdt1a=hpdt->GetBinError(hpdt->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      errmc1a=hpmc->GetBinError(hpmc->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      effdt2=tkhpdt->GetBinContent(tkhpdt->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      effmc2=tkhpmc->GetBinContent(tkhpmc->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      errdt2=tkhpdt->GetBinError(tkhpdt->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      errmc2=tkhpmc->GetBinError(tkhpmc->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      effdt2a=hpdt->GetBinContent(hpdt->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      effmc2a=hpmc->GetBinContent(hpmc->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      errdt2a=hpdt->GetBinError(hpdt->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      errmc2a=hpmc->GetBinError(hpmc->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      temp1=effdt1*effdt2a+effdt1a*effdt2-effdt1a*effdt2a;
      temp2=effmc1*effmc2a+effmc1a*effmc2-effmc1a*effmc2a;
      if(temp1&&temp2){
      idsfall=temp1/temp2;
      idsfallerr=(TMath::Power((effdt2-effdt2a)*errdt1a,2)+TMath::Power((effdt1-effdt1a)*errdt2a,2)+TMath::Power(effdt1a*errdt2,2)+TMath::Power(effdt2a*errdt1,2))/TMath::Power(temp1,2)+(TMath::Power((effmc2-effmc2a)*errmc1a,2)+TMath::Power((effmc1-effmc1a)*errmc2a,2)+TMath::Power(effmc1a*errmc2,2)+TMath::Power(effmc2a*errmc1,2))/TMath::Power(temp2,2);
      idsfallerr=TMath::Power(idsfallerr,.5)*idsfall;
      }
      else{
	idsfall=1;
	idsfallerr=1;}

      effdt1=isopteta->GetBinContent(isopteta->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      effdt2=isopteta->GetBinContent(isopteta->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      errdt1=isopteta->GetBinError(isopteta->FindBin(llnunu_l1_l1_eta,llnunu_l1_l1_pt));
      errdt2=isopteta->GetBinError(isopteta->FindBin(llnunu_l1_l2_eta,llnunu_l1_l2_pt));
      isosfall=effdt1*effdt2;
      isosfallerr=TMath::Power((TMath::Power(effdt1*errdt2,2)+TMath::Power(errdt1*effdt2,2)),.5);

      effdt1=mtksf->GetBinContent(mtksf->FindBin(llnunu_l1_l1_eta));
      effdt2=mtksf->GetBinContent(mtksf->FindBin(llnunu_l1_l2_eta));
      errdt1=mtksf->GetBinError(mtksf->FindBin(llnunu_l1_l1_eta));
      errdt2=mtksf->GetBinError(mtksf->FindBin(llnunu_l1_l2_eta));
      trksfall=effdt1*effdt2;
      trksfallerr=TMath::Power((TMath::Power(effdt1*errdt2,2)+TMath::Power(errdt1*effdt2,2)),.5);

      //effdt1=mul1pteta->GetBinContent(mul1pteta->FindBin(llnunu_l1_l1_pt,abs(llnunu_l1_l1_eta)))/100;
      //effdt2=mul2pteta->GetBinContent(mul2pteta->FindBin(llnunu_l1_l2_pt,abs(llnunu_l1_l2_eta)))/100;
      //errdt1=mul1pteta->GetBinError(mul1pteta->FindBin(llnunu_l1_l1_pt,abs(llnunu_l1_l1_eta)))/100;
      //errdt2=mul2pteta->GetBinError(mul2pteta->FindBin(llnunu_l1_l2_pt,abs(llnunu_l1_l2_eta)))/100;
      //trgsfall=effdt1+effdt2-effdt1*effdt2;
      //trgsfallerr=TMath::Power((TMath::Power((1-effdt1)*errdt2,2)+TMath::Power((1-effdt2)*errdt1,2)),.5);

      //      effdt1=mul1pteta->GetBinContent(mul1pteta->FindBin(llnunu_l1_l1_pt,abs(llnunu_l1_l1_eta)));
      //errdt1=mul1pteta->GetBinError(mul1pteta->FindBin(llnunu_l1_l1_pt,abs(llnunu_l1_l1_eta)));
      //trgsfall=effdt1;
      //trgsfallerr=errdt1;

      int trg_bin_l1 = htrg_l1_l1p_norm_vs_tot->FindBin(llnunu_l1_l1_pt,fabs(llnunu_l1_l1_eta));
      int trg_bin_l2 = htrg_l2_l1pl2f_norm_vs_l1p->FindBin(llnunu_l1_l2_pt,fabs(llnunu_l1_l2_eta));
      double trg_sc_l1_l1p_vs_tot = htrg_l1_l1p_norm_vs_tot->GetBinContent(trg_bin_l1);
      double trg_sc_l2_l1pl2f_vs_l1p = htrg_l2_l1pl2f_norm_vs_l1p->GetBinContent(trg_bin_l2);
      double trg_sc_l2_l1pl2p_vs_l1p = htrg_l2_l1pl2p_norm_vs_l1p->GetBinContent(trg_bin_l2);
      double trg_sc_l2_l1fl2p_vs_tot = htrg_l2_l1fl2p_norm_vs_tot->GetBinContent(trg_bin_l2);
      double trg_sc_l1_l1p_vs_tot_err = htrg_l1_l1p_norm_vs_tot->GetBinError(trg_bin_l1);
      double trg_sc_l2_l1pl2f_vs_l1p_err = htrg_l2_l1pl2f_norm_vs_l1p->GetBinError(trg_bin_l2);
      double trg_sc_l2_l1pl2p_vs_l1p_err = htrg_l2_l1pl2p_norm_vs_l1p->GetBinError(trg_bin_l2);
      double trg_sc_l2_l1fl2p_vs_tot_err = htrg_l2_l1fl2p_norm_vs_tot->GetBinError(trg_bin_l2);

      double trg_npass = Ntrg_l1pl2f*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2f_vs_l1p
                       + Ntrg_l1pl2p*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2p_vs_l1p
                       + Ntrg_l1fl2p*trg_sc_l2_l1fl2p_vs_tot
                       ;
      double trg_npass_err = pow(Ntrg_l1pl2f_err*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2f_vs_l1p,2)
                           + pow(Ntrg_l1pl2f*trg_sc_l1_l1p_vs_tot_err*trg_sc_l2_l1pl2f_vs_l1p,2)
                           + pow(Ntrg_l1pl2f*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2f_vs_l1p_err,2)
                           + pow(Ntrg_l1pl2p_err*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2p_vs_l1p,2)
                           + pow(Ntrg_l1pl2p*trg_sc_l1_l1p_vs_tot_err*trg_sc_l2_l1pl2p_vs_l1p,2)
                           + pow(Ntrg_l1pl2p*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2p_vs_l1p_err,2)
                           + pow(Ntrg_l1fl2p_err*trg_sc_l2_l1fl2p_vs_tot,2)
                           + pow(Ntrg_l1fl2p*trg_sc_l2_l1fl2p_vs_tot_err,2)
                           ;
      trg_npass_err = sqrt(trg_npass_err);

      double trg_nfail = Ntrg_tot-trg_npass;
      double trg_nfail_err = sqrt(Ntrg_tot_err*Ntrg_tot_err-Ntrg_l1pl2f_err*Ntrg_l1pl2f_err-Ntrg_l1pl2p_err*Ntrg_l1pl2p_err-Ntrg_l1fl2p_err*Ntrg_l1fl2p_err);

      double trg_eff = trg_npass/(trg_npass+trg_nfail);
      double trg_eff_err = (pow(trg_nfail*trg_npass_err,2)+pow(trg_npass*trg_nfail_err,2))/pow(trg_npass+trg_nfail,4);
      trg_eff_err = sqrt(trg_eff_err);

      double trg_eff_up = trg_eff+0.5*trg_eff_err;
      double trg_eff_dn = trg_eff-0.5*trg_eff_err;

      if (trg_eff>=1) trg_eff=1; 
      if (trg_eff<=0) trg_eff=0; 
      if (trg_eff_up>=1) trg_eff_up=1;
      if (trg_eff_dn>=1) trg_eff_dn=1;
      if (trg_eff_up<=0) trg_eff_up=0;
      if (trg_eff_dn<=0) trg_eff_dn=0;
      trg_eff_err = fabs(trg_eff_up-trg_eff_dn);

      trgsfall = trg_eff;
      trgsfallerr = trg_eff_err;
      trgsfall_up = trg_eff_up;
      trgsfall_dn = trg_eff_dn;

    }
    else if(abs(llnunu_l1_l1_pdgId)==11&&abs(llnunu_l1_l2_pdgId)==11){
      tkhighptdt1=1;
      tkhighptmc1=1;
      highptdt1=1;
      highptmc1=1;
      tkhighptdt2=1;
      tkhighptmc2=1;
      highptdt2=1;
      highptmc2=1;
      if(llnunu_l1_l1_pt>200) effdt1=1;
      else effdt1=esfh2->GetBinContent(esfh2->FindBin(llnunu_l1_l1_eSCeta,llnunu_l1_l1_pt));
      if(llnunu_l1_l2_pt>200) effdt2=1;
      else effdt2=esfh2->GetBinContent(esfh2->FindBin(llnunu_l1_l2_eSCeta,llnunu_l1_l2_pt));
      errdt1=esfh2->GetBinError(esfh2->FindBin(llnunu_l1_l1_eSCeta,llnunu_l1_l1_pt));
      errdt2=esfh2->GetBinError(esfh2->FindBin(llnunu_l1_l2_eSCeta,llnunu_l1_l2_pt));
      idsfall=effdt1*effdt2;
      idsfallerr=TMath::Power((TMath::Power(effdt1*errdt2,2)+TMath::Power(errdt1*effdt2,2)),.5);

      effdt1=etksf->GetBinContent(etksf->FindBin(llnunu_l1_l1_eSCeta,100));
      effdt2=etksf->GetBinContent(etksf->FindBin(llnunu_l1_l2_eSCeta,100));
      errdt1=etksf->GetBinError(etksf->FindBin(llnunu_l1_l1_eSCeta,100));
      errdt2=etksf->GetBinError(etksf->FindBin(llnunu_l1_l2_eSCeta,100));
      trksfall=effdt1*effdt2;
      trksfallerr=TMath::Power((TMath::Power(effdt1*errdt2,2)+TMath::Power(errdt1*effdt2,2)),.5);

      isosfall=1;
      isosfallerr=0;
      trgsfall=ell1pteta->GetBinContent(ell1pteta->FindBin(llnunu_l1_l1_pt,abs(llnunu_l1_l1_eta)))/100;
      trgsfallerr=ell1pteta->GetBinError(ell1pteta->FindBin(llnunu_l1_l1_pt,abs(llnunu_l1_l1_eta)))/100;
      //trgsfall=ell1pteta->GetBinContent(ell1pteta->FindBin(llnunu_l1_l1_pt,abs(llnunu_l1_l1_eta)));
      //trgsfallerr=ell1pteta->GetBinError(ell1pteta->FindBin(llnunu_l1_l1_pt,abs(llnunu_l1_l1_eta)));
      trgsfall_up = trgsfall+0.5*trgsfallerr;
      trgsfall_dn = trgsfall-0.5*trgsfallerr;
      if (trgsfall>=1) trgsfall=1;
      if (trgsfall<=0) trgsfall=0;
      if (trgsfall_up>=1) trgsfall_up=1;
      if (trgsfall_up<=0) trgsfall_up=0;
      if (trgsfall_dn>=1) trgsfall_dn=1;
      if (trgsfall_dn<=0) trgsfall_dn=0;
 
      trgsfallerr = fabs(trgsfall_up-trgsfall_dn);     
      
    }

    tree_out->Fill();

  }

  foutput->cd();
  tree_out->Write();
  foutput->Close();
  finput->Close();


  now = time(0);
  dt = ctime(&now);
  std::cout << "End time is: " << dt << std::endl;

  return 0;

}



