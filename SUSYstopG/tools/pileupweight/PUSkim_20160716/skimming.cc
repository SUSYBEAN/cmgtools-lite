#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TVector2.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>


bool doDyJets = true;
 
int main(int argc, char** argv) {

  if( argc<6 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: skimming... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " inputfile.root outputfile.root Nevts SumWeights cut_string " << std::endl ;
     exit(1) ;
  }

  // input file name
  std::string inputfile((const char*)argv[1]); 
  // output file name
  std::string outputfile((const char*)argv[2]);

  bool isDyJets = (outputfile.find("DYJets")!=std::string::npos);


  // nevts 
  double SumEvents = atof(argv[3]);  
  // sum weights 
  double SumWeights = atof(argv[4]);  
  // cuts
  std::string Cuts((const char*)argv[5]);

  // initialize
  // root files
  TFile* finput = new TFile(inputfile.c_str());
  TFile* foutput = new TFile(outputfile.c_str(), "recreate");


  // tree
  TTree* tree = (TTree*)finput->Get("tree");

  // out_tree
  //TTree* tree_out = tree->CopyTree(Cuts.c_str());
  TTree* tree_out = tree->CloneTree(0);;

  tree_out->Branch("SumEvents",&SumEvents,"SumEvents/D");
  tree_out->Branch("SumWeights",&SumWeights,"SumWeights/D");


  char name[3000];
  char name1[3000];


  Int_t isData;
  tree->SetBranchAddress("isData",&isData);
  
  // check if tree has events
  if (tree->GetEntries()<=0) return 0;

  tree->GetEntry(0);
  

  UInt_t          run;
  UInt_t          lumi;
  ULong64_t       evt;
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("lumi",&lumi);
  tree->SetBranchAddress("evt",&evt);
  Float_t genWeight;
  tree->SetBranchAddress("genWeight",&genWeight);

  // pileup file tags
  std::vector<std::string> pileup_tags = {
   //"68715", "68883", 
   "62154",
   "61665",
 //   "78000", "78500", "79000", "79500", "80000", "80500", "81000", "81500", "82000", "82500", "83000", "83500", "84000", "84500", "85000", "85500", "86000", "86500", "87000", "87500", "88000", "88500", "89000", "89500", "90000", "90500", "91000", "91500", "92000", "92500", "93000", "93500", "94000", "94500", "95000",
    "50000", "50500", "51000", "51500", "52000", "52500", "53000", "53500", "54000", "54500", "55000", "55500", "56000", "56500", "57000", 
    "61100", "61200", "61300", "61400",  "61600", "61700", "61800", "61900",  "62100", "62200", "62300", "62400",  "62600", "62700", "62800", "62900", 
     "57500", "58000", "58500", "59000", "59500", "60000", "60500", "61000", "61500", "62000", "62500", "63000", "63500", "64000", "64500", "65000", "65500", "66000", "66500", "67000", "67500", "68000", "68500", "69000", "69500", "70000", "70500", "71000", "71500", "72000", "72500", "73000", "73500", "74000", "74500", "75000", "75500", "76000", "76500", "77000", "77500",
 //"96000", "97000", "98000", "99000", "100000", "101000", "102000", "103000", "104000", "105000", "106000", "107000", "108000", "109000", "110000", "111000", "112000", "113000", "114000", "115000", "116000", "117000", "118000", "119000", "120000", "121000", "122000", "123000", "124000", "125000", "126000", "127000", "128000", "129000", "130000", "131000", "132000", "133000", "134000", "135000", "136000", "137000", "138000", "139000", "140000", "141000", "142000", "143000", "144000", "145000", "146000", "147000", "148000", "149000", "150000"
  //    "69452"
  //  "69100", "69200", "69300", "69400", "69500", "69600", "69700", "69800", "69900", "70000", "70100", "70200", "70300", "70400", "70500"
    };
  // pileup inputs
  std::vector<TFile*> pileup_files;
  std::vector<TH1D*> pileup_hists;
  std::vector<Double_t*> pileup_weights;
  std::vector<TBranch*> pileup_branches;

  if (!isData) {
    for (int i=0; i<(int)pileup_tags.size(); i++) {
      sprintf(name, "results/pileup_MC_80x_271036-276384_%s.root",pileup_tags.at(i).c_str());
      pileup_files.push_back(new TFile(name));
      pileup_hists.push_back((TH1D*)pileup_files.at(i)->Get("puweight_dtmc"));
      pileup_weights.push_back(new Double_t(1.0));
      sprintf(name, "puWeight%s", pileup_tags.at(i).c_str());
      sprintf(name1, "puWeight%s/D", pileup_tags.at(i).c_str());
      pileup_branches.push_back(tree_out->Branch(name,pileup_weights.at(i),name1));
    }
  }

  // doDyJetsSigma, correct dyjets mc to match data
  TFile* file_dt_sigma;
  TFile* file_mc_sigma;
  TH1D* h_dt_met_para_shift;
  TH1D* h_mc_met_para_shift;
  TH1D* h_dt_met_para_sigma;
  TH1D* h_dt_met_perp_sigma;
  TH1D* h_mc_met_para_sigma;
  TH1D* h_mc_met_perp_sigma;
  TH1D* h_ratio_met_para_sigma_dtmc;
  TH1D* h_ratio_met_perp_sigma_dtmc;

  Float_t llnunu_mt;
  Float_t llnunu_l1_mass;
  Float_t llnunu_l1_pt;
  Float_t llnunu_l1_phi;
  Float_t llnunu_l1_eta;
  Float_t llnunu_l2_pt;
  Float_t llnunu_l2_phi;
  Float_t llnunu_l1_l1_pt;
  Float_t llnunu_l1_l1_phi;
  Float_t llnunu_l1_l1_eta;
  Float_t llnunu_l1_l2_pt;
  Float_t llnunu_l1_l2_phi;
  Float_t llnunu_l1_l2_eta;

  tree->SetBranchAddress("llnunu_mt",&llnunu_mt);
  tree->SetBranchAddress("llnunu_l1_mass",&llnunu_l1_mass);
  tree->SetBranchAddress("llnunu_l1_pt",&llnunu_l1_pt);
  tree->SetBranchAddress("llnunu_l1_phi",&llnunu_l1_phi);
  tree->SetBranchAddress("llnunu_l1_eta",&llnunu_l1_eta);
  tree->SetBranchAddress("llnunu_l2_pt",&llnunu_l2_pt);
  tree->SetBranchAddress("llnunu_l2_phi",&llnunu_l2_phi);
  tree->SetBranchAddress("llnunu_l1_l1_pt",&llnunu_l1_l1_pt);
  tree->SetBranchAddress("llnunu_l1_l1_phi",&llnunu_l1_l1_phi);
  tree->SetBranchAddress("llnunu_l1_l1_eta",&llnunu_l1_l1_eta);
  tree->SetBranchAddress("llnunu_l1_l2_pt",&llnunu_l1_l2_pt);
  tree->SetBranchAddress("llnunu_l1_l2_phi",&llnunu_l1_l2_phi);
  tree->SetBranchAddress("llnunu_l1_l2_eta",&llnunu_l1_l2_eta);

  //
  Double_t ZPtWeight;
  Double_t PhiStarWeight;
  tree_out->Branch("ZPtWeight",&ZPtWeight,"ZPtWeight/D");
  tree_out->Branch("PhiStarWeight",&PhiStarWeight,"PhiStarWeight/D");

  TFile* fdyzpt;
  TH1D* hdyzptdt;
  TH1D* hdyzptmc;
  TH1D* hdyzpt_dtmc_ratio;
  TGraphErrors* gdyzpt_dtmc_ratio;

  TFile* fdyphistar;
  TH1D* hdyphistardt;
  TH1D* hdyphistarmc;
  TH1D* hdyphistar_dtmc_ratio;
  TGraphErrors* gdyphistar_dtmc_ratio;

  Int_t ngenZ;
  Int_t ngenLep;
  Float_t genZ_pt[10], genLep_eta[10], genLep_phi[10];
  if (!isData) {
    tree->SetBranchAddress("ngenZ", &ngenZ);
    tree->SetBranchAddress("genZ_pt", genZ_pt);
    tree->SetBranchAddress("ngenLep", &ngenLep);
    tree->SetBranchAddress("genLep_eta", genLep_eta);
    tree->SetBranchAddress("genLep_phi", genLep_phi);
  }


  // met shift  sigma
  file_dt_sigma = new TFile("SingleEMU_Run2016BCD_PromptReco_met_para_study.root");
  file_mc_sigma = new TFile("DYJetsToLL_M50_met_para_study.root");

  h_dt_met_para_shift = (TH1D*)file_dt_sigma->Get("h_met_para_vs_zpt_mean");
  h_mc_met_para_shift = (TH1D*)file_mc_sigma->Get("h_met_para_vs_zpt_mean");

  h_dt_met_para_sigma = (TH1D*)file_dt_sigma->Get("h_met_para_vs_zpt_sigma");
  h_dt_met_perp_sigma = (TH1D*)file_dt_sigma->Get("h_met_perp_vs_zpt_sigma");
  h_mc_met_para_sigma = (TH1D*)file_mc_sigma->Get("h_met_para_vs_zpt_sigma");
  h_mc_met_perp_sigma = (TH1D*)file_mc_sigma->Get("h_met_perp_vs_zpt_sigma");

  h_ratio_met_para_sigma_dtmc = (TH1D*)h_dt_met_para_sigma->Clone("h_ratio_met_para_sigma_dtmc");
  h_ratio_met_perp_sigma_dtmc = (TH1D*)h_dt_met_perp_sigma->Clone("h_ratio_met_perp_sigma_dtmc");
  h_ratio_met_para_sigma_dtmc->Divide(h_mc_met_para_sigma);
  h_ratio_met_perp_sigma_dtmc->Divide(h_mc_met_perp_sigma);

  fdyzpt = new TFile("UnfoldingOutputZPt.root");
  hdyzptdt = (TH1D*)fdyzpt->Get("hUnfold");
  hdyzptmc = (TH1D*)fdyzpt->Get("hTruth");
  hdyzpt_dtmc_ratio = (TH1D*)hdyzptdt->Clone("hdyzpt_dtmc_ratio");
  hdyzpt_dtmc_ratio->Divide(hdyzptmc);
  gdyzpt_dtmc_ratio = new TGraphErrors(hdyzpt_dtmc_ratio);

  fdyphistar = new TFile("UnfoldingOutputPhiStar.root");
  hdyphistardt = (TH1D*)fdyphistar->Get("hUnfold");
  hdyphistarmc = (TH1D*)fdyphistar->Get("hTruth");
  hdyphistar_dtmc_ratio = (TH1D*)hdyphistardt->Clone("hdyphistar_dtmc_ratio");
  hdyphistar_dtmc_ratio->Divide(hdyphistarmc);
  gdyphistar_dtmc_ratio = new TGraphErrors(hdyphistar_dtmc_ratio);





  //
  int n_interval = 5000;
  Int_t nTrueInt; 
  tree->SetBranchAddress( "nTrueInt", &nTrueInt );


  TTree* tree_out_duplic;
  if (isData){
    tree_out_duplic = tree->CloneTree(0);
    tree_out_duplic->SetName("tree_duplic");
  }

  for (int i=0; i<(int)tree->GetEntries(); i++){
  //for (int i=0; i<(int)1000; i++){
    tree->GetEntry(i);

    if ( i%n_interval == 0 ) {
      std::cout << "Event " << i << "   " << std::endl;
    }

    if (!isData) {
      for (int j=0; j<(int)pileup_tags.size(); j++){
        *(pileup_weights.at(j)) = pileup_hists.at(j)->GetBinContent(pileup_hists.at(j)->FindBin(nTrueInt));
      }
    }

    double met_para = llnunu_l2_pt*cos(llnunu_l2_phi-llnunu_l1_phi);
    double met_perp = llnunu_l2_pt*sin(llnunu_l2_phi-llnunu_l1_phi);
    if (isData) {
      met_para -= h_dt_met_para_shift->GetBinContent(h_dt_met_para_shift->FindBin(llnunu_l1_pt));
    }
    else {
      met_para -= h_mc_met_para_shift->GetBinContent(h_mc_met_para_shift->FindBin(llnunu_l1_pt));
    }

    if (!isData && isDyJets){
      met_para *= h_ratio_met_para_sigma_dtmc->GetBinContent(h_ratio_met_para_sigma_dtmc->FindBin(llnunu_l1_pt));
      met_perp *= h_ratio_met_perp_sigma_dtmc->GetBinContent(h_ratio_met_perp_sigma_dtmc->FindBin(llnunu_l1_pt));
    }
    double met_x = met_para*cos(llnunu_l1_phi)-met_perp*sin(llnunu_l1_phi);
    double met_y = met_para*sin(llnunu_l1_phi)+met_perp*cos(llnunu_l1_phi);
    TVector2 vec_met;
    vec_met.Set(met_x, met_y);
    llnunu_l2_pt = vec_met.Mod();
    llnunu_l2_phi = TVector2::Phi_mpi_pi(vec_met.Phi());

    double et1 = TMath::Sqrt(llnunu_l1_mass*llnunu_l1_mass + llnunu_l1_pt*llnunu_l1_pt);
    double et2 = TMath::Sqrt(llnunu_l1_mass*llnunu_l1_mass);
    llnunu_mt = TMath::Sqrt(2.0*llnunu_l1_mass*llnunu_l1_mass + 2.0* (et1*et2 - llnunu_l1_pt*cos(llnunu_l1_phi)*met_x - llnunu_l1_pt*sin(llnunu_l1_phi)*met_y));

    //
    if (!isData&&isDyJets) {
      // zpt weight
      if (ngenZ>0) ZPtWeight = gdyzpt_dtmc_ratio->Eval(genZ_pt[0]);
      else ZPtWeight = gdyzpt_dtmc_ratio->Eval(llnunu_l1_pt);
      // put zpt weight in genWeight
      genWeight *= ZPtWeight;

      // phistar weight
      double phistar(0.0);
      if (ngenLep>1) phistar = TMath::Tan((TMath::Pi()-TMath::Abs(TVector2::Phi_mpi_pi(genLep_phi[0]-genLep_phi[1])))/2.0)*TMath::Sin(TMath::ACos(TMath::TanH((genLep_eta[0]-genLep_eta[1])/2.0)));
      else phistar = TMath::Tan((TMath::Pi()-TMath::Abs(TVector2::Phi_mpi_pi(llnunu_l1_l1_phi-llnunu_l1_l2_phi)))/2.0)*TMath::Sin(TMath::ACos(TMath::TanH((llnunu_l1_l1_eta-llnunu_l1_l2_eta)/2.0)));
      PhiStarWeight = gdyphistar_dtmc_ratio->Eval(phistar);
    }

    sprintf(name, "run==%u&&lumi==%u&&evt==%llu", run,lumi,evt);
    if ( isData && (tree_out->GetEntries(name)) >0 ) {
      tree_out_duplic->Fill();
    }
    else {
      tree_out->Fill();
    }
  }

  foutput->cd();
  tree_out->Write();
  tree_out_duplic->Write();
  foutput->Close();
  finput->Close();

  return 0;

}



