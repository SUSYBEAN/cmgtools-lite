
#include<string>
#include<iostream>
#include"TFile.h"
#include "TROOT.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "tdrstyle.h"
#include "TLegend.h"
#include "TMath.h"

int main(){

  tdrstyle();
  //TFile* fin= new TFile("/data2/XZZ2/singlemuonnewtrg12p9.root");
  //TFile* fin= new TFile("/home/heli/XZZ/singlemuonnewtrg12p9.root");
  //TFile* fin= new TFile("/home/yanchu/singlemuonnewtrg12p9.root");
  TFile* fin= new TFile("/home/heli/XZZ/80X_20160825_light/SingleEMU_Run2016BCD_PromptReco/vvTreeProducer/tree.root");

  TFile* fout = new TFile("trigeff_mu_20160830.root", "recreate");
  char outputplot[] = "trigeff_mu_20160830_plots.pdf";
  TCanvas* plots = new TCanvas("plots", "plots");
  char name[10000];
  sprintf(name, "%s[", outputplot);
  plots->Print(name);
   

  std::string sel_l1_base = "((llnunu_l1_mass>70&&llnunu_l1_mass<110)&&abs(llnunu_l1_l1_pdgId)==13&&abs(llnunu_l1_l2_pdgId)==13&&llnunu_l1_l1_pt>20&&fabs(llnunu_l1_l1_eta)<2.4&&llnunu_l1_l2_pt>20&&fabs(llnunu_l1_l2_eta)<2.4&&(llnunu_l1_l1_highPtID>=0.99||llnunu_l1_l2_highPtID>=0.99))";
  std::string sel_l2_base = sel_l1_base+"&&(llnunu_l1_l1_pt>50)";
  //std::string sel_l2_base = sel_l1_base;
  std::string sel_l1_l1pl2f = sel_l1_base+"&&((llnunu_l1_l1_trigerob_HLTbit&12)!=0)&&((llnunu_l1_l2_trigerob_HLTbit&12)==0)";
  std::string sel_l1_l1pl2p = sel_l1_base+"&&((llnunu_l1_l1_trigerob_HLTbit&12)!=0)&&((llnunu_l1_l2_trigerob_HLTbit&12)!=0)";
  std::string sel_l1_l1fl2p = sel_l1_base+"&&((llnunu_l1_l1_trigerob_HLTbit&12)==0)&&((llnunu_l1_l2_trigerob_HLTbit&12)!=0)";
  std::string sel_l2_l1pl2f = sel_l2_base+"&&((llnunu_l1_l1_trigerob_HLTbit&12)!=0)&&((llnunu_l1_l2_trigerob_HLTbit&12)==0)";
  std::string sel_l2_l1pl2p = sel_l2_base+"&&((llnunu_l1_l1_trigerob_HLTbit&12)!=0)&&((llnunu_l1_l2_trigerob_HLTbit&12)!=0)";
  std::string sel_l2_l1fl2p = sel_l2_base+"&&((llnunu_l1_l1_trigerob_HLTbit&12)==0)&&((llnunu_l1_l2_trigerob_HLTbit&12)!=0)";

  std::string sel_l1_l1p = sel_l2_base+"&&((llnunu_l1_l1_trigerob_HLTbit&12)!=0)";
  std::string sel_l1_l1f = sel_l2_base+"&&((llnunu_l1_l1_trigerob_HLTbit&12)==0)";

  std::string sel_base = sel_l1_base+"&&(llnunu_l1_l1_pt>50)";
  //std::string sel_pass = sel_base+"&&(HLT_MUv2||HLT_ELEv2)";
  std::string sel_pass = sel_base+"&&(HLT_MUv2)";

  TTree* tree= (TTree*)fin->Get("tree");




  Double_t PtBins[] = {20, 25, 30, 32, 34, 36, 38, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 58, 60, 66, 72, 78, 84, 90, 96, 112, 118, 124, 130, 140, 160, 180, 200, 3000};
  Int_t NPtBins = sizeof(PtBins)/sizeof(PtBins[0]) - 1;
  
  Double_t EtaBins[] = {-0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.02,2.04,2.06,2.08, 2.1, 2.12,2.14, 2.16, 2.18, 2.2, 2.25, 2.5};
  Int_t NEtaBins = sizeof(EtaBins)/sizeof(EtaBins[0]) - 1;

  TH2D* htrg_l1_tot    = new TH2D("htrg_l1_tot",    "htrg_l1_tot",    NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l1_l1pl2f = new TH2D("htrg_l1_l1pl2f", "htrg_l1_l1pl2f", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l1_l1pl2p = new TH2D("htrg_l1_l1pl2p", "htrg_l1_l1pl2p", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l1_l1fl2p = new TH2D("htrg_l1_l1fl2p", "htrg_l1_l1fl2p", NPtBins, PtBins, NEtaBins, EtaBins);

  TH2D* htrg_l2_tot    = new TH2D("htrg_l2_tot",    "htrg_l2_tot",    NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l2_l1pl2f = new TH2D("htrg_l2_l1pl2f", "htrg_l2_l1pl2f", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l2_l1pl2p = new TH2D("htrg_l2_l1pl2p", "htrg_l2_l1pl2p", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l2_l1fl2p = new TH2D("htrg_l2_l1fl2p", "htrg_l2_l1fl2p", NPtBins, PtBins, NEtaBins, EtaBins);

  TH2D* htrg_l1_l1p = new TH2D("htrg_l1_l1p", "htrg_l1_l1p", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l1_l1f = new TH2D("htrg_l1_l1f", "htrg_l1_l1f", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l2_l1p = new TH2D("htrg_l2_l1p", "htrg_l2_l1p", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l2_l1f = new TH2D("htrg_l2_l1f", "htrg_l2_l1f", NPtBins, PtBins, NEtaBins, EtaBins);

  htrg_l1_tot->Sumw2();
  htrg_l1_l1pl2f->Sumw2();
  htrg_l1_l1pl2p->Sumw2();
  htrg_l1_l1fl2p->Sumw2();
  htrg_l2_tot->Sumw2();
  htrg_l2_l1pl2f->Sumw2();
  htrg_l2_l1pl2p->Sumw2();
  htrg_l2_l1fl2p->Sumw2();
  htrg_l1_l1p->Sumw2();
  htrg_l1_l1f->Sumw2();
  htrg_l2_l1p->Sumw2();
  htrg_l2_l1f->Sumw2();


  tree->Draw("fabs(llnunu_l1_l1_eta):llnunu_l1_l1_pt>>htrg_l1_tot",    sel_l1_base.c_str());
  tree->Draw("fabs(llnunu_l1_l1_eta):llnunu_l1_l1_pt>>htrg_l1_l1pl2f", sel_l1_l1pl2f.c_str());
  tree->Draw("fabs(llnunu_l1_l1_eta):llnunu_l1_l1_pt>>htrg_l1_l1pl2p", sel_l1_l1pl2p.c_str());
  tree->Draw("fabs(llnunu_l1_l1_eta):llnunu_l1_l1_pt>>htrg_l1_l1fl2p", sel_l1_l1fl2p.c_str());

  tree->Draw("fabs(llnunu_l1_l2_eta):llnunu_l1_l2_pt>>htrg_l2_tot",    sel_l2_base.c_str());
  tree->Draw("fabs(llnunu_l1_l2_eta):llnunu_l1_l2_pt>>htrg_l2_l1pl2f", sel_l2_l1pl2f.c_str());
  tree->Draw("fabs(llnunu_l1_l2_eta):llnunu_l1_l2_pt>>htrg_l2_l1pl2p", sel_l2_l1pl2p.c_str());
  tree->Draw("fabs(llnunu_l1_l2_eta):llnunu_l1_l2_pt>>htrg_l2_l1fl2p", sel_l2_l1fl2p.c_str());

  tree->Draw("fabs(llnunu_l1_l1_eta):llnunu_l1_l1_pt>>htrg_l1_l1p", sel_l1_l1p.c_str());
  tree->Draw("fabs(llnunu_l1_l1_eta):llnunu_l1_l1_pt>>htrg_l1_l1f", sel_l1_l1f.c_str());
  tree->Draw("fabs(llnunu_l1_l2_eta):llnunu_l1_l2_pt>>htrg_l2_l1p", sel_l1_l1p.c_str());
  tree->Draw("fabs(llnunu_l1_l2_eta):llnunu_l1_l2_pt>>htrg_l2_l1f", sel_l1_l1f.c_str());

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

  Ntrg_tot = htrg_l2_tot->IntegralAndError(Int_t(1), NPtBins, Int_t(1), NEtaBins, Ntrg_tot_err);
  Ntrg_l1p = htrg_l2_l1p->IntegralAndError(Int_t(1), NPtBins, Int_t(1), NEtaBins, Ntrg_l1p_err);
  Ntrg_l1f = htrg_l2_l1f->IntegralAndError(Int_t(1), NPtBins, Int_t(1), NEtaBins, Ntrg_l1f_err);
  Ntrg_l1pl2f = htrg_l2_l1pl2f->IntegralAndError(Int_t(1), NPtBins, Int_t(1), NEtaBins, Ntrg_l1pl2f_err);
  Ntrg_l1pl2p = htrg_l2_l1pl2p->IntegralAndError(Int_t(1), NPtBins, Int_t(1), NEtaBins, Ntrg_l1pl2p_err);
  Ntrg_l1fl2p = htrg_l2_l1fl2p->IntegralAndError(Int_t(1), NPtBins, Int_t(1), NEtaBins, Ntrg_l1fl2p_err);

  //Double_t Ntrg_tot    = htrg_l2_tot->Integral();
  //Double_t Ntrg_l1pl2f = htrg_l2_l1pl2f->Integral();
  //Double_t Ntrg_l1pl2p = htrg_l2_l1pl2p->Integral();
  //Double_t Ntrg_l1fl2p = htrg_l2_l1fl2p->Integral();
  //Double_t Ntrg_l1p = htrg_l2_l1p->Integral();
  //Double_t Ntrg_l1f = htrg_l2_l1f->Integral();

  Double_t eff_l1pl2f = Ntrg_l1pl2f/Ntrg_tot;
  Double_t eff_l1pl2p = Ntrg_l1pl2p/Ntrg_tot;
  Double_t eff_l1fl2p = Ntrg_l1fl2p/Ntrg_tot;
  Double_t eff_tot = (Ntrg_l1pl2f+Ntrg_l1pl2p+Ntrg_l1fl2p)/Ntrg_tot;

  std::cout << " Ntrg_l1p = " << Ntrg_l1p << "; Ntrg_l1f = " << Ntrg_l1f << std::endl;
  std::cout << " Ntrg_l1pl2f = " << Ntrg_l1pl2f << "; Ntrg_l1pl2p = " << Ntrg_l1pl2p << "; Ntrg_l1fl2p = " << Ntrg_l1fl2p << "; Ntrg_tot = " << Ntrg_tot << std::endl;
  std::cout << " Ntrg_l1p/Ntrg_tot = " << Ntrg_l1p/Ntrg_tot << "; Ntrg_l1f/Ntrg_tot = " << Ntrg_l1f/Ntrg_tot << std::endl;
  std::cout << " Ntrg_l1pl2f/Ntrg_tot  = " << Ntrg_l1pl2f/Ntrg_tot  << "; Ntrg_l1pl2p/Ntrg_tot  = " << Ntrg_l1pl2p/Ntrg_tot  << "; Ntrg_l1fl2p/Ntrg_tot  = " << Ntrg_l1fl2p/Ntrg_tot << std::endl;
  std::cout << " (Ntrg_l1pl2f+Ntrg_l1pl2p+Ntrg_l1fl2p)/Ntrg_tot = " << (Ntrg_l1pl2f+Ntrg_l1pl2p+Ntrg_l1fl2p)/Ntrg_tot << std::endl;


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

  // validation
  TH2D* htrg_l1_base = new TH2D("htrg_l1_base", "htrg_l1_base", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l2_base = new TH2D("htrg_l2_base", "htrg_l1_base", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l1_pass = new TH2D("htrg_l1_pass", "htrg_l1_pass", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l2_pass = new TH2D("htrg_l2_pass", "htrg_l1_pass", NPtBins, PtBins, NEtaBins, EtaBins);
  htrg_l1_base->Sumw2();
  htrg_l2_base->Sumw2();
  htrg_l1_pass->Sumw2();
  htrg_l2_pass->Sumw2();

  TH1D* htrg_l1_pt_base  = new TH1D("htrg_l1_pt_base", "htrg_l1_pt_base", 100, 0, 1000);
  TH1D* htrg_l2_pt_base  = new TH1D("htrg_l2_pt_base", "htrg_l2_pt_base", 100, 0, 500);
  TH1D* htrg_l1_pt_pass  = new TH1D("htrg_l1_pt_pass", "htrg_l1_pt_pass", 100, 0, 1000);
  TH1D* htrg_l2_pt_pass  = new TH1D("htrg_l2_pt_pass", "htrg_l2_pt_pass", 100, 0, 500);

  htrg_l1_pt_base->Sumw2();
  htrg_l2_pt_base->Sumw2();
  htrg_l1_pt_pass->Sumw2();
  htrg_l2_pt_pass->Sumw2();

  TH1D* htrg_l1_eta_base  = new TH1D("htrg_l1_eta_base", "htrg_l1_eta_base", 30, -3, 3);
  TH1D* htrg_l2_eta_base  = new TH1D("htrg_l2_eta_base", "htrg_l2_eta_base", 30, -3, 3);
  TH1D* htrg_l1_eta_pass  = new TH1D("htrg_l1_eta_pass", "htrg_l1_eta_pass", 30, -3, 3);
  TH1D* htrg_l2_eta_pass  = new TH1D("htrg_l2_eta_pass", "htrg_l2_eta_pass", 30, -3, 3);

  htrg_l1_eta_base->Sumw2();
  htrg_l2_eta_base->Sumw2();
  htrg_l1_eta_pass->Sumw2();
  htrg_l2_eta_pass->Sumw2();


  // get trig eff from data sample
  tree->Draw("fabs(llnunu_l1_l1_eta):llnunu_l1_l1_pt>>htrg_l1_base",  sel_base.c_str());
  tree->Draw("fabs(llnunu_l1_l2_eta):llnunu_l1_l2_pt>>htrg_l2_base",  sel_base.c_str());
  tree->Draw("fabs(llnunu_l1_l1_eta):llnunu_l1_l1_pt>>htrg_l1_pass",  sel_pass.c_str());
  tree->Draw("fabs(llnunu_l1_l2_eta):llnunu_l1_l2_pt>>htrg_l2_pass",  sel_pass.c_str());

  tree->Draw("llnunu_l1_l1_pt>>htrg_l1_pt_base",  sel_base.c_str());
  tree->Draw("llnunu_l1_l2_pt>>htrg_l2_pt_base",  sel_base.c_str());
  tree->Draw("llnunu_l1_l1_pt>>htrg_l1_pt_pass",  sel_pass.c_str());
  tree->Draw("llnunu_l1_l2_pt>>htrg_l2_pt_pass",  sel_pass.c_str());

  tree->Draw("llnunu_l1_l1_eta>>htrg_l1_eta_base",  sel_base.c_str());
  tree->Draw("llnunu_l1_l2_eta>>htrg_l2_eta_base",  sel_base.c_str());
  tree->Draw("llnunu_l1_l1_eta>>htrg_l1_eta_pass",  sel_pass.c_str());
  tree->Draw("llnunu_l1_l2_eta>>htrg_l2_eta_pass",  sel_pass.c_str());

  // trig eff from data
  TH2D* htrg_eff_l1_trig = (TH2D*)htrg_l1_pass->Clone("htrg_eff_l1_trig");
  TH2D* htrg_eff_l2_trig = (TH2D*)htrg_l2_pass->Clone("htrg_eff_l2_trig");

  htrg_eff_l1_trig->Divide(htrg_l1_base);
  htrg_eff_l2_trig->Divide(htrg_l2_base);

  TH1D* htrg_eff_l1_pt_trig = (TH1D*)htrg_l1_pt_pass->Clone("htrg_eff_l1_pt_trig");
  TH1D* htrg_eff_l2_pt_trig = (TH1D*)htrg_l2_pt_pass->Clone("htrg_eff_l2_pt_trig");
  TH1D* htrg_eff_l1_eta_trig = (TH1D*)htrg_l1_eta_pass->Clone("htrg_eff_l1_eta_trig");
  TH1D* htrg_eff_l2_eta_trig = (TH1D*)htrg_l2_eta_pass->Clone("htrg_eff_l2_eta_trig");

  htrg_eff_l1_pt_trig->Divide(htrg_l1_pt_base);
  htrg_eff_l2_pt_trig->Divide(htrg_l2_pt_base);
  htrg_eff_l1_eta_trig->Divide(htrg_l1_eta_base);
  htrg_eff_l2_eta_trig->Divide(htrg_l2_eta_base);

  // test trig eff by calculation
  TH2D* htrg_eff_l1_test = new TH2D("htrg_eff_l1_test", "htrg_eff_l1_test", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_eff_l2_test = new TH2D("htrg_eff_l2_test", "htrg_eff_l1_test", NPtBins, PtBins, NEtaBins, EtaBins);
  htrg_eff_l1_test->Sumw2();
  htrg_eff_l2_test->Sumw2();

  for (int ix=1; ix<=htrg_eff_l1_test->GetNbinsX(); ix++){
    for (int iy=1; iy<=htrg_eff_l1_test->GetNbinsY(); iy++){
      double sf1 = Ntrg_l1pl2f/Ntrg_tot
                  *htrg_l1_l1pl2f_norm_vs_tot->GetBinContent(ix,iy)
                 + Ntrg_l1pl2p/Ntrg_tot
                  *htrg_l1_l1pl2p_norm_vs_tot->GetBinContent(ix,iy)
                 + Ntrg_l1fl2p/Ntrg_tot
                  *htrg_l1_l1fl2p_norm_vs_tot->GetBinContent(ix,iy)
                 ;
      double sf2 = Ntrg_l1pl2f/Ntrg_tot
                  *htrg_l2_l1pl2f_norm_vs_tot->GetBinContent(ix,iy)
                 + Ntrg_l1pl2p/Ntrg_tot
                  *htrg_l2_l1pl2p_norm_vs_tot->GetBinContent(ix,iy)
                 + Ntrg_l1fl2p/Ntrg_tot
                  *htrg_l2_l1fl2p_norm_vs_tot->GetBinContent(ix,iy)
                 ;
      double er1 = Ntrg_l1pl2f/Ntrg_tot
                  *htrg_l1_l1pl2f_norm_vs_tot->GetBinError(ix,iy)
                 + Ntrg_l1pl2p/Ntrg_tot
                  *htrg_l1_l1pl2p_norm_vs_tot->GetBinError(ix,iy)
                 + Ntrg_l1fl2p/Ntrg_tot
                  *htrg_l1_l1fl2p_norm_vs_tot->GetBinError(ix,iy)
                 ;
      double er2 = Ntrg_l1pl2f/Ntrg_tot
                  *htrg_l2_l1pl2f_norm_vs_tot->GetBinError(ix,iy)
                 + Ntrg_l1pl2p/Ntrg_tot
                  *htrg_l2_l1pl2p_norm_vs_tot->GetBinError(ix,iy)
                 + Ntrg_l1fl2p/Ntrg_tot
                  *htrg_l2_l1fl2p_norm_vs_tot->GetBinError(ix,iy)
                 ;
      htrg_eff_l1_test->SetBinContent(ix, iy, sf1);
      htrg_eff_l2_test->SetBinContent(ix, iy, sf2);
      htrg_eff_l1_test->SetBinError(ix, iy, er1);
      htrg_eff_l2_test->SetBinError(ix, iy, er2);

    }
  }

 
  // check real sf
  Float_t llnunu_l1_mass;
  Float_t llnunu_l1_l1_phi, llnunu_l1_l1_eta, llnunu_l1_l2_phi, llnunu_l1_l2_eta,llnunu_l1_l1_pt,llnunu_l1_l2_pt;
  Int_t llnunu_l1_l1_pdgId, llnunu_l1_l2_pdgId;
  Float_t llnunu_l1_l1_highPtID, llnunu_l1_l2_highPtID;
  tree->SetBranchAddress("llnunu_l1_mass",&llnunu_l1_mass);
  tree->SetBranchAddress("llnunu_l1_l1_pdgId",&llnunu_l1_l1_pdgId);
  tree->SetBranchAddress("llnunu_l1_l2_pdgId",&llnunu_l1_l2_pdgId);
  tree->SetBranchAddress("llnunu_l1_l1_highPtID",&llnunu_l1_l1_highPtID);
  tree->SetBranchAddress("llnunu_l1_l2_highPtID",&llnunu_l1_l2_highPtID);
  tree->SetBranchAddress("llnunu_l1_l1_phi",&llnunu_l1_l1_phi);
  tree->SetBranchAddress("llnunu_l1_l2_phi",&llnunu_l1_l2_phi);
  tree->SetBranchAddress("llnunu_l1_l1_pt",&llnunu_l1_l1_pt);
  tree->SetBranchAddress("llnunu_l1_l2_pt",&llnunu_l1_l2_pt);
  tree->SetBranchAddress("llnunu_l1_l1_eta",&llnunu_l1_l1_eta);
  tree->SetBranchAddress("llnunu_l1_l2_eta",&llnunu_l1_l2_eta);

  TH2D* htrg_l1_tgsf = new TH2D("htrg_l1_tgsf", "htrg_l1_tgsf", NPtBins, PtBins, NEtaBins, EtaBins);
  TH2D* htrg_l2_tgsf = new TH2D("htrg_l2_tgsf", "htrg_l1_tgsf", NPtBins, PtBins, NEtaBins, EtaBins);
  htrg_l1_tgsf->Sumw2();
  htrg_l2_tgsf->Sumw2();

  TH1D* htrg_l1_pt_tgsf  = new TH1D("htrg_l1_pt_tgsf", "htrg_l1_pt_tgsf", 100, 0, 1000);
  TH1D* htrg_l2_pt_tgsf  = new TH1D("htrg_l2_pt_tgsf", "htrg_l2_pt_tgsf", 100, 0, 500);
  htrg_l1_pt_tgsf->Sumw2();
  htrg_l2_pt_tgsf->Sumw2();

  TH1D* htrg_l1_eta_tgsf  = new TH1D("htrg_l1_eta_tgsf", "htrg_l1_eta_tgsf", 30, -3, 3);
  TH1D* htrg_l2_eta_tgsf  = new TH1D("htrg_l2_eta_tgsf", "htrg_l2_eta_tgsf", 30, -3, 3);
  htrg_l1_eta_tgsf->Sumw2();
  htrg_l2_eta_tgsf->Sumw2();

  for (int i=0; i<(int)tree->GetEntries(); i++){
    tree->GetEntry(i);
    if ( !((llnunu_l1_mass>70&&llnunu_l1_mass<110)&&abs(llnunu_l1_l1_pdgId)==13&&abs(llnunu_l1_l2_pdgId)==13&&llnunu_l1_l1_pt>50&&fabs(llnunu_l1_l1_eta)<2.4&&llnunu_l1_l2_pt>20&&fabs(llnunu_l1_l2_eta)<2.4&&(llnunu_l1_l1_highPtID>=0.991||llnunu_l1_l2_highPtID>=0.99)) ) {
      continue;
    }

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
 
    htrg_l1_tgsf->Fill(llnunu_l1_l1_pt, fabs(llnunu_l1_l1_eta), trg_eff);
    htrg_l2_tgsf->Fill(llnunu_l1_l2_pt, fabs(llnunu_l1_l2_eta), trg_eff);
    htrg_l1_pt_tgsf->Fill(llnunu_l1_l1_pt, trg_eff);
    htrg_l2_pt_tgsf->Fill(llnunu_l1_l2_pt, trg_eff);
    htrg_l1_eta_tgsf->Fill(llnunu_l1_l1_eta, trg_eff);
    htrg_l2_eta_tgsf->Fill(llnunu_l1_l2_eta, trg_eff);

  }

  // eff by sf
  TH2D* htrg_eff_l1_tgsf = (TH2D*)htrg_l1_tgsf->Clone("htrg_eff_l1_tgsf");
  TH2D* htrg_eff_l2_tgsf = (TH2D*)htrg_l2_tgsf->Clone("htrg_eff_l2_tgsf"); 
  htrg_eff_l1_tgsf->Divide(htrg_l1_base);
  htrg_eff_l2_tgsf->Divide(htrg_l2_base);

  TH1D* htrg_eff_l1_pt_tgsf = (TH1D*)htrg_l1_pt_tgsf->Clone("htrg_eff_l1_pt_tgsf");
  TH1D* htrg_eff_l2_pt_tgsf = (TH1D*)htrg_l2_pt_tgsf->Clone("htrg_eff_l2_pt_tgsf");
  htrg_eff_l1_pt_tgsf->Divide(htrg_l1_pt_base);
  htrg_eff_l2_pt_tgsf->Divide(htrg_l2_pt_base);

  TH1D* htrg_eff_l1_eta_tgsf = (TH1D*)htrg_l1_eta_tgsf->Clone("htrg_eff_l1_eta_tgsf");
  TH1D* htrg_eff_l2_eta_tgsf = (TH1D*)htrg_l2_eta_tgsf->Clone("htrg_eff_l2_eta_tgsf");
  htrg_eff_l1_eta_tgsf->Divide(htrg_l1_eta_base);
  htrg_eff_l2_eta_tgsf->Divide(htrg_l2_eta_base);

  Double_t N_base = htrg_l1_base->Integral();
  Double_t N_pass = htrg_l1_pass->Integral();
  Double_t N_tgsf = htrg_l1_tgsf->Integral();
  std::cout << "N_pass = " << N_pass << "; N_tgsf = " << N_tgsf << "; N_base = " << N_base << std::endl;
  std::cout << "N_pass/N_tgsf = " << N_pass/N_tgsf << "; N_pass/N_base = " << N_pass/N_base << "; N_tgsf/N_base = " << N_tgsf/N_base << std::endl;


  TH1D* htrg_l1_pt_pass_vs_tgsf = (TH1D*)htrg_l1_pt_pass->Clone("htrg_l1_pt_pass_vs_tgsf");
  htrg_l1_pt_pass_vs_tgsf->Divide(htrg_l1_pt_tgsf);
  TH1D* htrg_l2_pt_pass_vs_tgsf = (TH1D*)htrg_l2_pt_pass->Clone("htrg_l2_pt_pass_vs_tgsf");
  htrg_l2_pt_pass_vs_tgsf->Divide(htrg_l2_pt_tgsf);
  TH1D* htrg_l1_eta_pass_vs_tgsf = (TH1D*)htrg_l1_eta_pass->Clone("htrg_l1_eta_pass_vs_tgsf");
  htrg_l1_eta_pass_vs_tgsf->Divide(htrg_l1_eta_tgsf);
  TH1D* htrg_l2_eta_pass_vs_tgsf = (TH1D*)htrg_l2_eta_pass->Clone("htrg_l2_eta_pass_vs_tgsf");
  htrg_l2_eta_pass_vs_tgsf->Divide(htrg_l2_eta_tgsf);
  TH1D* htrg_eff_l1_pt_trig_vs_tgsf = (TH1D*)htrg_eff_l1_pt_trig->Clone("htrg_eff_l1_pt_trig_vs_tgsf");
  htrg_eff_l1_pt_trig_vs_tgsf->Divide(htrg_eff_l1_pt_tgsf);
  TH1D* htrg_eff_l2_pt_trig_vs_tgsf = (TH1D*)htrg_eff_l2_pt_trig->Clone("htrg_eff_l2_pt_trig_vs_tgsf");
  htrg_eff_l2_pt_trig_vs_tgsf->Divide(htrg_eff_l2_pt_tgsf);
  TH1D* htrg_eff_l1_eta_trig_vs_tgsf = (TH1D*)htrg_eff_l1_eta_trig->Clone("htrg_eff_l1_eta_trig_vs_tgsf");
  htrg_eff_l1_eta_trig_vs_tgsf->Divide(htrg_eff_l1_eta_tgsf);
  TH1D* htrg_eff_l2_eta_trig_vs_tgsf = (TH1D*)htrg_eff_l2_eta_trig->Clone("htrg_eff_l2_eta_trig_vs_tgsf");
  htrg_eff_l2_eta_trig_vs_tgsf->Divide(htrg_eff_l2_eta_tgsf);


  htrg_l1_pt_base->SetMarkerStyle(20);
  htrg_l2_pt_base->SetMarkerStyle(20);
  htrg_l1_pt_pass->SetMarkerStyle(20);
  htrg_l2_pt_pass->SetMarkerStyle(20);
  htrg_l1_eta_base->SetMarkerStyle(20);
  htrg_l2_eta_base->SetMarkerStyle(20);
  htrg_l1_eta_pass->SetMarkerStyle(20);
  htrg_l2_eta_pass->SetMarkerStyle(20);
  htrg_eff_l1_pt_trig->SetMarkerStyle(20);
  htrg_eff_l2_pt_trig->SetMarkerStyle(20);
  htrg_eff_l1_eta_trig->SetMarkerStyle(20);
  htrg_eff_l2_eta_trig->SetMarkerStyle(20);
  htrg_l1_pt_tgsf->SetMarkerStyle(20);
  htrg_l2_pt_tgsf->SetMarkerStyle(20);
  htrg_l1_eta_tgsf->SetMarkerStyle(20);
  htrg_l2_eta_tgsf->SetMarkerStyle(20);
  htrg_eff_l1_pt_tgsf->SetMarkerStyle(20);
  htrg_eff_l2_pt_tgsf->SetMarkerStyle(20);
  htrg_eff_l1_eta_tgsf->SetMarkerStyle(20);
  htrg_eff_l2_eta_tgsf->SetMarkerStyle(20);

  htrg_l1_pt_base->SetMarkerSize(0.5);
  htrg_l2_pt_base->SetMarkerSize(0.5);
  htrg_l1_pt_pass->SetMarkerSize(0.5);
  htrg_l2_pt_pass->SetMarkerSize(0.5);
  htrg_l1_eta_base->SetMarkerSize(0.5);
  htrg_l2_eta_base->SetMarkerSize(0.5);
  htrg_l1_eta_pass->SetMarkerSize(0.5);
  htrg_l2_eta_pass->SetMarkerSize(0.5);
  htrg_eff_l1_pt_trig->SetMarkerSize(0.5);
  htrg_eff_l2_pt_trig->SetMarkerSize(0.5);
  htrg_eff_l1_eta_trig->SetMarkerSize(0.5);
  htrg_eff_l2_eta_trig->SetMarkerSize(0.5);
  htrg_l1_pt_tgsf->SetMarkerSize(0.5);
  htrg_l2_pt_tgsf->SetMarkerSize(0.5);
  htrg_l1_eta_tgsf->SetMarkerSize(0.5);
  htrg_l2_eta_tgsf->SetMarkerSize(0.5);
  htrg_eff_l1_pt_tgsf->SetMarkerSize(0.5);
  htrg_eff_l2_pt_tgsf->SetMarkerSize(0.5);
  htrg_eff_l1_eta_tgsf->SetMarkerSize(0.5);
  htrg_eff_l2_eta_tgsf->SetMarkerSize(0.5);

  htrg_l1_pt_base->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_pt_base->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_pt_pass->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_pt_pass->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l1_pt_trig->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l2_pt_trig->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_pt_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_pt_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l1_pt_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l2_pt_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_eta_base->GetXaxis()->SetTitle("#eta");
  htrg_l2_eta_base->GetXaxis()->SetTitle("#eta");
  htrg_l1_eta_pass->GetXaxis()->SetTitle("#eta");
  htrg_l2_eta_pass->GetXaxis()->SetTitle("#eta");
  htrg_eff_l1_eta_trig->GetXaxis()->SetTitle("#eta");
  htrg_eff_l2_eta_trig->GetXaxis()->SetTitle("#eta");
  htrg_l1_eta_tgsf->GetXaxis()->SetTitle("#eta");
  htrg_l2_eta_tgsf->GetXaxis()->SetTitle("#eta");
  htrg_eff_l1_eta_tgsf->GetXaxis()->SetTitle("#eta");
  htrg_eff_l2_eta_tgsf->GetXaxis()->SetTitle("#eta");

  htrg_l1_pt_base->GetYaxis()->SetTitle("Entries");
  htrg_l2_pt_base->GetYaxis()->SetTitle("Entries");
  htrg_l1_pt_pass->GetYaxis()->SetTitle("Entries");
  htrg_l2_pt_pass->GetYaxis()->SetTitle("Entries");
  htrg_l1_eta_base->GetYaxis()->SetTitle("Entries");
  htrg_l2_eta_base->GetYaxis()->SetTitle("Entries");
  htrg_l1_eta_pass->GetYaxis()->SetTitle("Entries");
  htrg_l2_eta_pass->GetYaxis()->SetTitle("Entries");
  htrg_l1_pt_tgsf->GetYaxis()->SetTitle("Entries");
  htrg_l2_pt_tgsf->GetYaxis()->SetTitle("Entries");
  htrg_l1_eta_tgsf->GetYaxis()->SetTitle("Entries");
  htrg_l2_eta_tgsf->GetYaxis()->SetTitle("Entries");
  htrg_eff_l1_pt_trig->GetYaxis()->SetTitle("eff.");
  htrg_eff_l2_pt_trig->GetYaxis()->SetTitle("eff.");
  htrg_eff_l1_eta_trig->GetYaxis()->SetTitle("eff.");
  htrg_eff_l2_eta_trig->GetYaxis()->SetTitle("eff.");
  htrg_eff_l1_pt_tgsf->GetYaxis()->SetTitle("eff.");
  htrg_eff_l2_pt_tgsf->GetYaxis()->SetTitle("eff.");
  htrg_eff_l1_eta_tgsf->GetYaxis()->SetTitle("eff.");
  htrg_eff_l2_eta_tgsf->GetYaxis()->SetTitle("eff.");

  htrg_l1_pt_pass_vs_tgsf->SetMarkerStyle(20);
  htrg_l2_pt_pass_vs_tgsf->SetMarkerStyle(20);
  htrg_l1_eta_pass_vs_tgsf->SetMarkerStyle(20);
  htrg_l2_eta_pass_vs_tgsf->SetMarkerStyle(20);
  htrg_eff_l1_pt_trig_vs_tgsf->SetMarkerStyle(20);
  htrg_eff_l2_pt_trig_vs_tgsf->SetMarkerStyle(20);
  htrg_eff_l1_eta_trig_vs_tgsf->SetMarkerStyle(20);
  htrg_eff_l2_eta_trig_vs_tgsf->SetMarkerStyle(20);
  htrg_l1_pt_pass_vs_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_pt_pass_vs_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l1_pt_trig_vs_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l2_pt_trig_vs_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_eta_pass_vs_tgsf->GetXaxis()->SetTitle("#eta");
  htrg_l2_eta_pass_vs_tgsf->GetXaxis()->SetTitle("#eta");
  htrg_eff_l1_eta_trig_vs_tgsf->GetXaxis()->SetTitle("#eta");
  htrg_eff_l2_eta_trig_vs_tgsf->GetXaxis()->SetTitle("#eta");

  htrg_l1_tot_norm->SetTitle("htrg_l1_tot_norm");
  htrg_l2_tot_norm->SetTitle("htrg_l2_tot_norm");
  htrg_l1_l1p_norm->SetTitle("htrg_l1_l1p_norm");
  htrg_l1_l1f_norm->SetTitle("htrg_l1_l1f_norm");
  htrg_l2_l1p_norm->SetTitle("htrg_l2_l1p_norm");
  htrg_l2_l1f_norm->SetTitle("htrg_l2_l1f_norm");
  htrg_l1_l1pl2f_norm->SetTitle("htrg_l1_l1pl2f_norm");
  htrg_l1_l1pl2p_norm->SetTitle("htrg_l1_l1pl2p_norm");
  htrg_l1_l1fl2p_norm->SetTitle("htrg_l1_l1fl2p_norm");
  htrg_l2_l1pl2f_norm->SetTitle("htrg_l2_l1pl2f_norm");
  htrg_l2_l1pl2p_norm->SetTitle("htrg_l2_l1pl2p_norm");
  htrg_l2_l1fl2p_norm->SetTitle("htrg_l2_l1fl2p_norm");
  htrg_l1_l1p_norm_vs_tot->SetTitle("htrg_l1_l1p_norm_vs_tot");
  htrg_l1_l1f_norm_vs_tot->SetTitle("htrg_l1_l1f_norm_vs_tot");
  htrg_l2_l1p_norm_vs_tot->SetTitle("htrg_l2_l1p_norm_vs_tot");
  htrg_l2_l1f_norm_vs_tot->SetTitle("htrg_l2_l1f_norm_vs_tot");
  htrg_l1_l1pl2f_norm_vs_tot->SetTitle("htrg_l1_l1pl2f_norm_vs_tot");
  htrg_l1_l1pl2p_norm_vs_tot->SetTitle("htrg_l1_l1pl2p_norm_vs_tot");
  htrg_l1_l1fl2p_norm_vs_tot->SetTitle("htrg_l1_l1fl2p_norm_vs_tot");
  htrg_l2_l1pl2f_norm_vs_tot->SetTitle("htrg_l2_l1pl2f_norm_vs_tot");
  htrg_l2_l1pl2p_norm_vs_tot->SetTitle("htrg_l2_l1pl2p_norm_vs_tot");
  htrg_l2_l1fl2p_norm_vs_tot->SetTitle("htrg_l2_l1fl2p_norm_vs_tot");
  htrg_l1_l1pl2f_norm_vs_l1p->SetTitle("htrg_l1_l1pl2f_norm_vs_l1p");
  htrg_l1_l1pl2p_norm_vs_l1p->SetTitle("htrg_l1_l1pl2p_norm_vs_l1p");
  htrg_l1_l1fl2p_norm_vs_l1f->SetTitle("htrg_l1_l1fl2p_norm_vs_l1f");
  htrg_l2_l1pl2f_norm_vs_l1p->SetTitle("htrg_l2_l1pl2f_norm_vs_l1p");
  htrg_l2_l1pl2p_norm_vs_l1p->SetTitle("htrg_l2_l1pl2p_norm_vs_l1p");
  htrg_l2_l1fl2p_norm_vs_l1f->SetTitle("htrg_l2_l1fl2p_norm_vs_l1f");
  htrg_eff_l1_trig->SetTitle("htrg_eff_l1_trig");
  htrg_eff_l2_trig->SetTitle("htrg_eff_l2_trig");
  htrg_eff_l1_pt_trig->SetTitle("htrg_eff_l1_pt_trig");
  htrg_eff_l2_pt_trig->SetTitle("htrg_eff_l2_pt_trig");
  htrg_eff_l1_eta_trig->SetTitle("htrg_eff_l1_eta_trig");
  htrg_eff_l2_eta_trig->SetTitle("htrg_eff_l2_eta_trig");
  htrg_eff_l1_tgsf->SetTitle("htrg_eff_l1_tgsf");
  htrg_eff_l2_tgsf->SetTitle("htrg_eff_l2_tgsf");
  htrg_eff_l1_pt_tgsf->SetTitle("htrg_eff_l1_pt_tgsf");
  htrg_eff_l2_pt_tgsf->SetTitle("htrg_eff_l2_pt_tgsf");
  htrg_eff_l1_eta_tgsf->SetTitle("htrg_eff_l1_eta_tgsf");
  htrg_eff_l2_eta_tgsf->SetTitle("htrg_eff_l2_eta_tgsf");
  htrg_l1_pt_pass_vs_tgsf->SetTitle("htrg_l1_pt_pass_vs_tgsf");
  htrg_l2_pt_pass_vs_tgsf->SetTitle("htrg_l2_pt_pass_vs_tgsf");
  htrg_l1_eta_pass_vs_tgsf->SetTitle("htrg_l1_eta_pass_vs_tgsf");
  htrg_l2_eta_pass_vs_tgsf->SetTitle("htrg_l2_eta_pass_vs_tgsf");
  htrg_eff_l1_pt_trig_vs_tgsf->SetTitle("htrg_eff_l1_pt_trig_vs_tgsf");
  htrg_eff_l2_pt_trig_vs_tgsf->SetTitle("htrg_eff_l2_pt_trig_vs_tgsf");
  htrg_eff_l1_eta_trig_vs_tgsf->SetTitle("htrg_eff_l1_eta_trig_vs_tgsf");
  htrg_eff_l2_eta_trig_vs_tgsf->SetTitle("htrg_eff_l2_eta_trig_vs_tgsf");

  TLegend* lg[100];

  lg[0] = new TLegend(0.3,0.3,0.7,0.5);  
  lg[0]->SetName("lg_htrg_pt_base_l1_l2");
  htrg_l1_pt_base->SetLineColor(2);
  htrg_l2_pt_base->SetLineColor(4);
  htrg_l1_pt_base->SetMarkerColor(2);
  htrg_l2_pt_base->SetMarkerColor(4);
  lg[0]->AddEntry(htrg_l1_pt_base, "l_{1} no trig.", "pl");
  lg[0]->AddEntry(htrg_l2_pt_base, "l_{2} no trig.", "pl");

  plots->cd();
  plots->Clear();
  htrg_l1_pt_base->Draw();
  htrg_l2_pt_base->Draw("same");
  lg[0]->Draw("same");
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();
 

  lg[1] = new TLegend(0.3,0.3,0.7,0.5);
  lg[1]->SetName("lg_htrg_l1_pt_pass_tgsf");
  htrg_l1_pt_pass->SetLineColor(2);
  htrg_l1_pt_tgsf->SetLineColor(4);
  htrg_l1_pt_pass->SetMarkerColor(2);
  htrg_l1_pt_tgsf->SetMarkerColor(4);
  lg[1]->AddEntry(htrg_l1_pt_pass, "l_{1} pass trig.", "pl");
  lg[1]->AddEntry(htrg_l1_pt_tgsf, "l_{1} trig. s.f.", "pl");

  plots->cd();
  plots->Clear();
  htrg_l1_pt_pass->Draw();
  htrg_l1_pt_tgsf->Draw("same");
  lg[1]->Draw("same");
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();

  lg[2] = new TLegend(0.3,0.3,0.7,0.5);
  lg[2]->SetName("lg_htrg_l2_pt_pass_tgsf");
  htrg_l2_pt_pass->SetLineColor(2);
  htrg_l2_pt_tgsf->SetLineColor(4);
  htrg_l2_pt_pass->SetMarkerColor(2);
  htrg_l2_pt_tgsf->SetMarkerColor(4);
  lg[2]->AddEntry(htrg_l2_pt_pass, "l_{2} pass trig.", "pl");
  lg[2]->AddEntry(htrg_l2_pt_tgsf, "l_{2} trig. s.f.", "pl");

  plots->cd();
  plots->Clear();
  htrg_l2_pt_pass->Draw();
  htrg_l2_pt_tgsf->Draw("same");
  lg[2]->Draw("same");
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();

  lg[3] = new TLegend(0.3,0.3,0.7,0.5);
  lg[3]->SetName("lg_htrg_eff_l1_pt_trig_tgsf");
  htrg_eff_l1_pt_trig->SetLineColor(2);
  htrg_eff_l1_pt_tgsf->SetLineColor(4);
  htrg_eff_l1_pt_trig->SetMarkerColor(2);
  htrg_eff_l1_pt_tgsf->SetMarkerColor(4);
  lg[3]->AddEntry(htrg_eff_l1_pt_trig, "l_{1} pass trig.", "pl");
  lg[3]->AddEntry(htrg_eff_l1_pt_tgsf, "l_{1} trig. s.f.", "pl");

  plots->cd();
  plots->Clear();
  htrg_eff_l1_pt_trig->Draw();
  htrg_eff_l1_pt_tgsf->Draw("same");
  lg[3]->Draw("same");
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();

  lg[4] = new TLegend(0.3,0.3,0.7,0.5);
  lg[4]->SetName("lg_htrg_eff_l2_pt_trig_tgsf");
  htrg_eff_l2_pt_trig->SetLineColor(2);
  htrg_eff_l2_pt_tgsf->SetLineColor(4);
  htrg_eff_l2_pt_trig->SetMarkerColor(2);
  htrg_eff_l2_pt_tgsf->SetMarkerColor(4);
  lg[4]->AddEntry(htrg_eff_l2_pt_trig, "l_{2} pass trig.", "pl");
  lg[4]->AddEntry(htrg_eff_l2_pt_tgsf, "l_{2} trig. s.f.", "pl");

  plots->cd();
  plots->Clear();
  htrg_eff_l2_pt_trig->Draw();
  htrg_eff_l2_pt_tgsf->Draw("same");
  lg[4]->Draw("same");
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();


  // eta
  lg[5] = new TLegend(0.3,0.3,0.7,0.5);
  lg[5]->SetName("lg_htrg_eta_base_l1_l2");
  htrg_l1_eta_base->SetLineColor(2);
  htrg_l2_eta_base->SetLineColor(4);
  htrg_l1_eta_base->SetMarkerColor(2);
  htrg_l2_eta_base->SetMarkerColor(4);
  lg[5]->AddEntry(htrg_l1_eta_base, "l_{1} no trig.", "pl");
  lg[5]->AddEntry(htrg_l2_eta_base, "l_{2} no trig.", "pl");

  plots->cd();
  plots->Clear();
  htrg_l1_eta_base->Draw();
  htrg_l2_eta_base->Draw("same");
  lg[5]->Draw("same");
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();


  lg[6] = new TLegend(0.3,0.3,0.7,0.5);
  lg[6]->SetName("lg_htrg_l1_eta_pass_tgsf");
  htrg_l1_eta_pass->SetLineColor(2);
  htrg_l1_eta_tgsf->SetLineColor(4);
  htrg_l1_eta_pass->SetMarkerColor(2);
  htrg_l1_eta_tgsf->SetMarkerColor(4);
  lg[6]->AddEntry(htrg_l1_eta_pass, "l_{1} pass trig.", "pl");
  lg[6]->AddEntry(htrg_l1_eta_tgsf, "l_{1} trig. s.f.", "pl");

  plots->cd();
  plots->Clear();
  htrg_l1_eta_pass->Draw();
  htrg_l1_eta_tgsf->Draw("same");
  lg[6]->Draw("same");
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();

  lg[7] = new TLegend(0.3,0.3,0.7,0.5);
  lg[7]->SetName("lg_htrg_l2_eta_pass_tgsf");
  htrg_l2_eta_pass->SetLineColor(2);
  htrg_l2_eta_tgsf->SetLineColor(4);
  htrg_l2_eta_pass->SetMarkerColor(2);
  htrg_l2_eta_tgsf->SetMarkerColor(4);
  lg[7]->AddEntry(htrg_l2_eta_pass, "l_{2} pass trig.", "pl");
  lg[7]->AddEntry(htrg_l2_eta_tgsf, "l_{2} trig. s.f.", "pl");

  plots->cd();
  plots->Clear();
  htrg_l2_eta_pass->Draw();
  htrg_l2_eta_tgsf->Draw("same");
  lg[7]->Draw("same");
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();
 
  lg[8] = new TLegend(0.3,0.3,0.7,0.5);
  lg[8]->SetName("lg_htrg_eff_l1_eta_trig_tgsf");
  htrg_eff_l1_eta_trig->SetLineColor(2);
  htrg_eff_l1_eta_tgsf->SetLineColor(4);
  htrg_eff_l1_eta_trig->SetMarkerColor(2);
  htrg_eff_l1_eta_tgsf->SetMarkerColor(4);
  lg[8]->AddEntry(htrg_eff_l1_eta_trig, "l_{1} pass trig.", "pl");
  lg[8]->AddEntry(htrg_eff_l1_eta_tgsf, "l_{1} trig. s.f.", "pl");

  plots->cd();
  plots->Clear();
  htrg_eff_l1_eta_trig->Draw();
  htrg_eff_l1_eta_tgsf->Draw("same");
  lg[8]->Draw("same");
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();

  lg[9] = new TLegend(0.3,0.3,0.7,0.5);
  lg[9]->SetName("lg_htrg_eff_l2_eta_trig_tgsf");
  htrg_eff_l2_eta_trig->SetLineColor(2);
  htrg_eff_l2_eta_tgsf->SetLineColor(4);
  htrg_eff_l2_eta_trig->SetMarkerColor(2);
  htrg_eff_l2_eta_tgsf->SetMarkerColor(4);
  lg[9]->AddEntry(htrg_eff_l2_eta_trig, "l_{2} pass trig.", "pl");
  lg[9]->AddEntry(htrg_eff_l2_eta_tgsf, "l_{2} trig. s.f.", "pl");

  plots->cd();
  plots->Clear();
  htrg_eff_l2_eta_trig->Draw();
  htrg_eff_l2_eta_tgsf->Draw("same");
  lg[9]->Draw("same");
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();

  plots->cd();
  plots->Clear();
  htrg_l1_pt_pass_vs_tgsf->Draw();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_pt_pass_vs_tgsf->Draw();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_eta_pass_vs_tgsf->Draw();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_eta_pass_vs_tgsf->Draw();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_eff_l1_pt_trig_vs_tgsf->Draw();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_eff_l2_pt_trig_vs_tgsf->Draw();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_eff_l1_eta_trig_vs_tgsf->Draw();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_eff_l2_eta_trig_vs_tgsf->Draw();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->Clear();

 
  // 2D
  htrg_l1_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1pl2f->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1pl2p->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1fl2p->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1pl2f->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1pl2p->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1fl2p->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1p->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1f->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1p->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1f->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_tot_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_tot_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1p_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1f_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1p_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1f_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1pl2f_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1pl2p_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1fl2p_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1pl2f_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1pl2p_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1fl2p_norm->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1p_norm_vs_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1f_norm_vs_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1p_norm_vs_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1f_norm_vs_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1pl2f_norm_vs_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1pl2p_norm_vs_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1fl2p_norm_vs_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1pl2f_norm_vs_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1pl2p_norm_vs_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1fl2p_norm_vs_tot->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1pl2f_norm_vs_l1p->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1pl2p_norm_vs_l1p->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_l1fl2p_norm_vs_l1f->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1pl2f_norm_vs_l1p->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1pl2p_norm_vs_l1p->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_l1fl2p_norm_vs_l1f->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_base->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_base->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_pass->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_pass->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l1_trig->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l2_trig->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l1_test->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l2_test->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l2_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l1_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_eff_l2_tgsf->GetXaxis()->SetTitle("p_{T} (GeV)");
  htrg_l1_tot->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1pl2f->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1pl2p->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1fl2p->GetYaxis()->SetTitle("#eta");
  htrg_l2_tot->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1pl2f->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1pl2p->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1fl2p->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1p->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1f->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1p->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1f->GetYaxis()->SetTitle("#eta");
  htrg_l1_tot_norm->GetYaxis()->SetTitle("#eta");
  htrg_l2_tot_norm->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1p_norm->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1f_norm->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1p_norm->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1f_norm->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1pl2f_norm->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1pl2p_norm->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1fl2p_norm->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1pl2f_norm->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1pl2p_norm->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1fl2p_norm->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1p_norm_vs_tot->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1f_norm_vs_tot->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1p_norm_vs_tot->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1f_norm_vs_tot->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1pl2f_norm_vs_tot->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1pl2p_norm_vs_tot->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1fl2p_norm_vs_tot->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1pl2f_norm_vs_tot->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1pl2p_norm_vs_tot->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1fl2p_norm_vs_tot->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1pl2f_norm_vs_l1p->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1pl2p_norm_vs_l1p->GetYaxis()->SetTitle("#eta");
  htrg_l1_l1fl2p_norm_vs_l1f->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1pl2f_norm_vs_l1p->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1pl2p_norm_vs_l1p->GetYaxis()->SetTitle("#eta");
  htrg_l2_l1fl2p_norm_vs_l1f->GetYaxis()->SetTitle("#eta");
  htrg_l1_base->GetYaxis()->SetTitle("#eta");
  htrg_l2_base->GetYaxis()->SetTitle("#eta");
  htrg_l1_pass->GetYaxis()->SetTitle("#eta");
  htrg_l2_pass->GetYaxis()->SetTitle("#eta");
  htrg_eff_l1_trig->GetYaxis()->SetTitle("#eta");
  htrg_eff_l2_trig->GetYaxis()->SetTitle("#eta");
  htrg_eff_l1_test->GetYaxis()->SetTitle("#eta");
  htrg_eff_l2_test->GetYaxis()->SetTitle("#eta");
  htrg_l1_tgsf->GetYaxis()->SetTitle("#eta");
  htrg_l2_tgsf->GetYaxis()->SetTitle("#eta");
  htrg_eff_l1_tgsf->GetYaxis()->SetTitle("#eta");
  htrg_eff_l2_tgsf->GetYaxis()->SetTitle("#eta");

  plots->cd();
  plots->Clear();
  htrg_l1_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1pl2f->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1pl2p->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1fl2p->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1pl2f->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1pl2p->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1fl2p->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1p->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1f->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1p->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1f->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_tot_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_tot_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1p_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1f_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1p_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1f_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1pl2f_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1pl2p_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1fl2p_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1pl2f_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1pl2p_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1fl2p_norm->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1p_norm_vs_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1f_norm_vs_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1p_norm_vs_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1f_norm_vs_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1pl2f_norm_vs_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1pl2p_norm_vs_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1fl2p_norm_vs_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1pl2f_norm_vs_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1pl2p_norm_vs_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1fl2p_norm_vs_tot->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1pl2f_norm_vs_l1p->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1pl2p_norm_vs_l1p->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_l1fl2p_norm_vs_l1f->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1pl2f_norm_vs_l1p->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1pl2p_norm_vs_l1p->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_l1fl2p_norm_vs_l1f->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_base->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_base->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_pass->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_pass->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_eff_l1_trig->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_eff_l2_trig->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_eff_l1_test->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_eff_l2_test->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l1_tgsf->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_l2_tgsf->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_eff_l1_tgsf->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
 
  plots->cd();
  plots->Clear();
  htrg_eff_l2_tgsf->Draw("colz text");
  plots->SetLogx();
  sprintf(name, outputplot);
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();

  // end
  sprintf(name, "%s]", outputplot);
  plots->Print(name);

  // save 

  fout->cd();
  htrg_l1_pt_base->Write();
  htrg_l2_pt_base->Write();
  htrg_l1_pt_pass->Write();
  htrg_l2_pt_pass->Write();
  htrg_l1_eta_base->Write();
  htrg_l2_eta_base->Write();
  htrg_l1_eta_pass->Write();
  htrg_l2_eta_pass->Write();
  htrg_eff_l1_pt_trig->Write();
  htrg_eff_l2_pt_trig->Write();
  htrg_eff_l1_eta_trig->Write();
  htrg_eff_l2_eta_trig->Write();
  htrg_l1_pt_tgsf->Write();
  htrg_l2_pt_tgsf->Write();
  htrg_l1_eta_tgsf->Write();
  htrg_l2_eta_tgsf->Write();
  htrg_eff_l1_pt_tgsf->Write();
  htrg_eff_l2_pt_tgsf->Write();
  htrg_eff_l1_eta_tgsf->Write();
  htrg_eff_l2_eta_tgsf->Write();
  htrg_l1_tot->Write();
  htrg_l1_l1pl2f->Write();
  htrg_l1_l1pl2p->Write();
  htrg_l1_l1fl2p->Write();
  htrg_l2_tot->Write();
  htrg_l2_l1pl2f->Write();
  htrg_l2_l1pl2p->Write();
  htrg_l2_l1fl2p->Write();
  htrg_l1_l1p->Write();
  htrg_l1_l1f->Write();
  htrg_l2_l1p->Write();
  htrg_l2_l1f->Write();
  htrg_l1_tot_norm->Write();
  htrg_l2_tot_norm->Write();
  htrg_l1_l1p_norm->Write();
  htrg_l1_l1f_norm->Write();
  htrg_l2_l1p_norm->Write();
  htrg_l2_l1f_norm->Write();
  htrg_l1_l1pl2f_norm->Write();
  htrg_l1_l1pl2p_norm->Write();
  htrg_l1_l1fl2p_norm->Write();
  htrg_l2_l1pl2f_norm->Write();
  htrg_l2_l1pl2p_norm->Write();
  htrg_l2_l1fl2p_norm->Write();
  htrg_l1_l1p_norm_vs_tot->Write();
  htrg_l1_l1f_norm_vs_tot->Write();
  htrg_l2_l1p_norm_vs_tot->Write();
  htrg_l2_l1f_norm_vs_tot->Write();
  htrg_l1_l1pl2f_norm_vs_tot->Write();
  htrg_l1_l1pl2p_norm_vs_tot->Write();
  htrg_l1_l1fl2p_norm_vs_tot->Write();
  htrg_l2_l1pl2f_norm_vs_tot->Write();
  htrg_l2_l1pl2p_norm_vs_tot->Write();
  htrg_l2_l1fl2p_norm_vs_tot->Write();
  htrg_l1_l1pl2f_norm_vs_l1p->Write();
  htrg_l1_l1pl2p_norm_vs_l1p->Write();
  htrg_l1_l1fl2p_norm_vs_l1f->Write();
  htrg_l2_l1pl2f_norm_vs_l1p->Write();
  htrg_l2_l1pl2p_norm_vs_l1p->Write();
  htrg_l2_l1fl2p_norm_vs_l1f->Write();
  htrg_l1_base->Write();
  htrg_l2_base->Write();
  htrg_l1_pass->Write();
  htrg_l2_pass->Write();
  htrg_eff_l1_trig->Write();
  htrg_eff_l2_trig->Write();
  htrg_eff_l1_test->Write();
  htrg_eff_l2_test->Write();
  htrg_l1_tgsf->Write();
  htrg_l2_tgsf->Write();
  htrg_eff_l1_tgsf->Write();
  htrg_eff_l2_tgsf->Write();
  fout->Close();

// if use pt1>50 for Npass/Nall, gives:
// N_pass = 2.34669e+06; N_tgsf = 2.40924e+06; N_pass/N_tgsf = 0.974037
// 
}
