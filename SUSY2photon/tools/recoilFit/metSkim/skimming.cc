#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TMath.h"
#include "TVector2.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

int main(int argc, char** argv) {

  if( argc<5 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: skimming... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " inputfile.root outputfile.root Nevts SumWeights " << std::endl ;
     exit(1) ;
  }

  // input file name
  std::string inputfile((const char*)argv[1]); 
  // output file name
  std::string outputfile((const char*)argv[2]);


  // nevts 
  double SumEvents = atof(argv[3]);  
  // sum weights 
  double SumWeights = atof(argv[4]);  

  // initialize
  // root files
  TFile* finput = new TFile(inputfile.c_str());
  TFile* foutput = new TFile(outputfile.c_str(), "recreate");


  // tree
  TTree* tree = (TTree*)finput->Get("tree");

  // out_tree
  TTree* tree_out = tree->CloneTree(0);

//  TBranch *b_SumEvents=tree_out->Branch("SumEvents",&SumEvents,"SumEvents/D");
//  TBranch *b_SumWeights=tree_out->Branch("SumWeights",&SumWeights,"SumWeights/D");


  // ZPT PhiStar reweight model
/*
root [4] TH1D* hdt = (TH1D*)gDirectory->Get("hUnfold");
root [5] TH1D* hmc = (TH1D*)gDirectory->Get("hTruth");
root [6] TH1D* hratio = (TH1D*)hmc->Clone("hmcdt");
root [7] hratio->Divid(hdt);
*/
  TFile* fdyzpt = new TFile("UnfoldingOutputZPt.root");
  TH1D* hdyzptdt = (TH1D*)fdyzpt->Get("hUnfold");  
  TH1D* hdyzptmc = (TH1D*)fdyzpt->Get("hTruth");  
  TH1D* hdyzpt_dtmc_ratio = (TH1D*)hdyzptdt->Clone("hdyzpt_dtmc_ratio");
  hdyzpt_dtmc_ratio->Divide(hdyzptmc);

  TGraphErrors* gdyzpt_dtmc_ratio = new TGraphErrors(hdyzpt_dtmc_ratio);

  TFile* fdyphistar = new TFile("UnfoldingOutputPhiStar.root");
  TH1D* hdyphistardt = (TH1D*)fdyphistar->Get("hUnfold");
  TH1D* hdyphistarmc = (TH1D*)fdyphistar->Get("hTruth");
  TH1D* hdyphistar_dtmc_ratio = (TH1D*)hdyphistardt->Clone("hdyphistar_dtmc_ratio");
  hdyphistar_dtmc_ratio->Divide(hdyphistarmc);

  TGraphErrors* gdyphistar_dtmc_ratio = new TGraphErrors(hdyphistar_dtmc_ratio);

  //
  Double_t ZPtWeight;
  TBranch *b_ZPtWeight=tree_out->Branch("ZPtWeight",&ZPtWeight,"ZPtWeight/D");
  Double_t PhiStarWeight;
  TBranch *b_PhiStarWeight=tree_out->Branch("PhiStarWeight",&PhiStarWeight,"PhiStarWeight/D");

  // metskim
/*  
  TFile* fpar = new TFile("fit_recoil_par.root");
  TH1D* h_peak_upara = (TH1D*)fpar->Get("h_peak_upara");
  TH1D* h_shift_upara = (TH1D*)fpar->Get("h_shift_upara");
  TH1D* h_sigma_upara = (TH1D*)fpar->Get("h_sigma_upara");
  TH1D* h_peak_uperp = (TH1D*)fpar->Get("h_peak_uperp");
  TH1D* h_shift_uperp = (TH1D*)fpar->Get("h_shift_uperp");
  TH1D* h_sigma_uperp = (TH1D*)fpar->Get("h_sigma_uperp");

  h_peak_upara->Smooth(1);
  h_shift_upara->Smooth(1);
  h_sigma_upara->Smooth(1);
  h_peak_uperp->Smooth(1);
  h_shift_uperp->Smooth(1);
  h_sigma_uperp->Smooth(1);

  TGraphErrors* g_peak_upara = new TGraphErrors(h_peak_upara);
  TGraphErrors* g_shift_upara = new TGraphErrors(h_shift_upara);
  TGraphErrors* g_sigma_upara = new TGraphErrors(h_sigma_upara);
  TGraphErrors* g_peak_uperp = new TGraphErrors(h_peak_uperp);
  TGraphErrors* g_shift_uperp = new TGraphErrors(h_shift_uperp);
  TGraphErrors* g_sigma_uperp = new TGraphErrors(h_sigma_uperp);

  g_peak_upara->SetName("g_peak_upara");
  g_shift_upara->SetName("g_shift_upara");
  g_sigma_upara->SetName("g_sigma_upara");
  g_peak_uperp->SetName("g_peak_uperp");
  g_shift_uperp->SetName("g_shift_uperp");
  g_sigma_uperp->SetName("g_sigma_uperp");
*/
/*
 llnunu_l2_pt    = 31.7863
 llnunu_l2_phi   = 1.49733

 llnunu_l1_pt    = 31.9389
 llnunu_l1_phi   = 0.681347

 llnunu_deltaPhi = 0.81598
 llnunu_mta      = 180.101
 llnunu_CosdphiZMet = 0.685155
 llnunu_dPTPara  = 53.7175
 llnunu_dPTParaRel = 1.68188
 llnunu_dPTPerp  = 23.1531
 llnunu_dPTPerpRel = 0.724918

*/
  Float_t llnunu_l1_pt, llnunu_l2_pt, llnunu_l1_phi, llnunu_l2_phi, llnunu_l1_mass;
  Float_t llnunu_mta, llnunu_deltaPhi, llnunu_CosdphiZMet, llnunu_dPTPara, llnunu_dPTParaRel, llnunu_dPTPerp, llnunu_dPTPerpRel;
  Float_t llnunu_l1_l1_phi, llnunu_l1_l1_eta, llnunu_l1_l2_phi, llnunu_l1_l2_eta;
  Float_t llnunu_l1_deltaPhi;
  tree->SetBranchAddress("llnunu_l1_mass",&llnunu_l1_mass);
  tree->SetBranchAddress("llnunu_l1_pt",&llnunu_l1_pt);
  tree->SetBranchAddress("llnunu_l2_pt",&llnunu_l2_pt);
  tree->SetBranchAddress("llnunu_l1_phi",&llnunu_l1_phi);
  tree->SetBranchAddress("llnunu_l2_phi",&llnunu_l2_phi);
  tree->SetBranchAddress("llnunu_l1_l1_phi",&llnunu_l1_l1_phi);
  tree->SetBranchAddress("llnunu_l1_l2_phi",&llnunu_l1_l2_phi);
  tree->SetBranchAddress("llnunu_l1_l1_eta",&llnunu_l1_l1_eta);
  tree->SetBranchAddress("llnunu_l1_l2_eta",&llnunu_l1_l2_eta);
  tree->SetBranchAddress("llnunu_l1_deltaPhi",&llnunu_l1_deltaPhi);
  tree->SetBranchAddress("llnunu_mta",&llnunu_mta);
  tree->SetBranchAddress("llnunu_deltaPhi", &llnunu_deltaPhi);
  //tree->SetBranchAddress("llnunu_CosdphiZMet",&llnunu_CosdphiZMet);
  //tree->SetBranchAddress("llnunu_dPTPara",&llnunu_dPTPara);
  //tree->SetBranchAddress("llnunu_dPTParaRel",&llnunu_dPTParaRel);
  //tree->SetBranchAddress("llnunu_dPTPerp",&llnunu_dPTPerp);
  //tree->SetBranchAddress("llnunu_dPTPerpRel",&llnunu_dPTPerpRel);

  Float_t genZ_pt[10], genLep_eta[10], genLep_phi[10]; 
  tree->SetBranchAddress("genZ_pt", genZ_pt);
  tree->SetBranchAddress("genLep_eta", genLep_eta);
  tree->SetBranchAddress("genLep_phi", genLep_phi);


  for (int i=0; i<(int)tree->GetEntries(); i++){
  //for (int i=0; i<10000; i++){
    tree->GetEntry(i);
    //b_SumEvents->Fill();
    //b_SumWeights->Fill();

    // zpt weight
    //ZPtWeight = gdyzpt_dtmc_ratio->Eval(llnunu_l1_pt);
    ZPtWeight = gdyzpt_dtmc_ratio->Eval(genZ_pt[0]);

    // phistar weight
    //double phistar = TMath::Tan((TMath::Pi()-TMath::Abs(llnunu_l1_deltaPhi))/2.0)*TMath::Sin(TMath::ACos(TMath::TanH((llnunu_l1_l1_eta-llnunu_l1_l2_eta)/2.0)));
    double phistar = TMath::Tan((TMath::Pi()-TMath::Abs(TVector2::Phi_mpi_pi(genLep_phi[0]-genLep_phi[1])))/2.0)*TMath::Sin(TMath::ACos(TMath::TanH((genLep_eta[0]-genLep_eta[1])/2.0)));
    PhiStarWeight = gdyphistar_dtmc_ratio->Eval(phistar);
/*
    // 
    double zmass = llnunu_l1_mass;
    double zpt0 = llnunu_l1_pt;
    double met0 = llnunu_l2_pt;
    double zptphi0 = llnunu_l1_phi;
    double metphi0 = llnunu_l2_phi;
    double zpx0 = zpt0*TMath::Cos(zptphi0);
    double zpy0 = zpt0*TMath::Sin(zptphi0);

    double cosdphi0 = TMath::Cos(metphi0-zptphi0);
    double sindphi0 = TMath::Sin(metphi0-zptphi0);

    double metPara0 = met0*cosdphi0;
    double metPerp0 = met0*sindphi0;

    double peak_upara = g_peak_upara->Eval(zpt0);
    double shift_upara = g_shift_upara->Eval(zpt0);
    double sigma_upara = g_sigma_upara->Eval(zpt0);
    double peak_uperp = g_peak_uperp->Eval(zpt0);
    double shift_uperp = g_shift_uperp->Eval(zpt0);
    double sigma_uperp = g_sigma_uperp->Eval(zpt0);

    // apply para corr
    double metPara1 = metPara0 + shift_upara*zpt0;
    // apply para smear
    metPara1 = metPara1 + sigma_upara * (metPara1 - peak_upara*zpt0);
 
    // apply perp corr
    double metPerp1 = metPerp0 + shift_uperp*zpt0;
    //double metPerp1 = metPerp0;// + shift_uperp*zpt0;

    // met1 in z frame 
    TVector2 vect_met1(metPara1, metPerp1);

    // rotate to original frame
    vect_met1 = vect_met1.Rotate(zptphi0);

    //
    double met1 = vect_met1.Mod();
    double metphi1 = TVector2::Phi_mpi_pi(vect_met1.Phi());
    double metx1 = vect_met1.Px();
    double mety1 = vect_met1.Py();
    double dphi1 = TMath::Abs(TVector2::Phi_mpi_pi(metphi1-zptphi0));  
    double cosdphi1 = TMath::Cos(dphi1);
    double sindphi1 = TMath::Sin(dphi1);
   
    double et1 = TMath::Sqrt(zmass*zmass+zpt0*zpt0);
    double et2 = TMath::Sqrt(91.188*91.188+met1*met1);
    double mta = TMath::Sqrt(zmass*zmass + 91.188*91.188 + 2.0* (et1*et2 - zpx0*metx1 - zpy0*mety1));

    llnunu_l2_pt    = met1;
    llnunu_l2_phi   = metphi1;

    llnunu_deltaPhi = dphi1;
    llnunu_mta      = mta;
    llnunu_CosdphiZMet = cosdphi1;
    llnunu_dPTPara  = zpt0 + met1*cosdphi1;
    llnunu_dPTParaRel = 1.0 + met1*cosdphi1/zpt0;
    llnunu_dPTPerp  = met1*sindphi1;
    llnunu_dPTPerpRel = met1*sindphi1/zpt0;
*/ 
    tree_out->Fill();

  }

  foutput->cd();
  tree_out->Write();
  foutput->Close();
  finput->Close();

  return 0;

}



