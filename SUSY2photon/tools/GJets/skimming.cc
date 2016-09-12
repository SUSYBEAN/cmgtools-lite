#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TVector2.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <utility>
#include "TLorentzVector.h"
#include "KalmanMuonCalibrator.h"

bool doDyJets = true;
bool doDyJetsUseSmooth = false;
bool doDyJetsUseFunction = true;
bool doDyJetsLOUseSmooth = false;
bool doDyJetsLOUseFunction = true;
bool addZrapidity = false;
bool doRecoil = false;
bool correctData = false;
bool doRecoilUseSmooth = true;
bool doRecoilUseGraph = true;

bool correctMuonPt = false;

bool lightWeight = false;
bool addZjetsNewGenWeight = true;

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
  bool isDyJetsLO = (outputfile.find("DYJets")!=std::string::npos && outputfile.find("MGMLM")!=std::string::npos);


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

  // light weight skim
  if (lightWeight){
    tree->SetBranchStatus("lep_*",0);
    tree->SetBranchStatus("jet_*",0);
    tree->SetBranchStatus("llnunu_LV_*",0);
    tree->SetBranchStatus("llnunu_TuneP_LV_*",0);
    tree->SetBranchStatus("triggers*",0);
    tree->SetBranchStatus("llnunu_l1_l1_leps*",0);
    tree->SetBranchStatus("llnunu_l1_l2_leps*",0);
  }

  // out_tree
  //TTree* tree_out = tree->CopyTree(Cuts.c_str());
  TTree* tree_out = tree->CloneTree(0);

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

  Float_t ZJetsLOSumWeights(49877138);
  Float_t ZJetsNLOSumWeights(450670522117);
  Float_t ZJetsLOSumEvents(49877138);
  Float_t ZJetsNLOSumEvents(28696958);
  Float_t genWeight;
  if (!isData) tree->SetBranchAddress("genWeight",&genWeight);


  Float_t ZJetsGenWeight;
  if (addZjetsNewGenWeight and !isData and isDyJets) {
    tree_out->Branch("ZJetsGenWeight", &ZJetsGenWeight, "ZJetsGenWeight/F");  
  }
  // pileup file tags
  std::vector<std::string> pileup_tags = {
    "69200", "68075", "67921",
   // "62194", "61674", "62118", "61651", "61658",
   //"68715", "68883", 
  // "62154",
  // "61665",
 //   "78000", "78500", "79000", "79500", "80000", "80500", "81000", "81500", "82000", "82500", "83000", "83500", "84000", "84500", "85000", "85500", "86000", "86500", "87000", "87500", "88000", "88500", "89000", "89500", "90000", "90500", "91000", "91500", "92000", "92500", "93000", "93500", "94000", "94500", "95000",
  //  "50000", "50500", "51000", "51500", "52000", "52500", "53000", "53500", "54000", "54500", "55000", "55500", "56000", "56500", "57000", 
//    "61100", "61200", "61300", "61400",  "61600", "61700", "61800", "61900",  "62100", "62200", "62300", "62400",  "62600", "62700", "62800", "62900", 
//     "57500", "58000", "58500", "59000", "59500", "60000", "60500", "61000", "61500", "62000", "62500", "63000", "63500", "64000", "64500", "65000", "65500", "66000", "66500", "67000", "67500", "68000", "68500", "69000", "69500", "70000", "70500", "71000", "71500", "72000", "72500", "73000", "73500", "74000", "74500", "75000", "75500", "76000", "76500", "77000", "77500",

// "68100", "68200", "68300", "68400", "68600", "68700", "68800", "68900", "69100", "69200", "69300", "69400", "69600", "69700", "69800", "69900", "70100", "70200", "70300", "70400", "70600", "70700", "70800", "70900",
// "63100", "63200", "63300", "63400", "63600", "63700", "63800", "63900", "64100", "64200", "64300", "64400", "64600", "64700", "64800", "64900", "65100", "65200", "65300", "65400", "65600", "65700", "65800", "65900", "66100", "66200", "66300", "66400", "66600", "66700", "66800", "66900", "67100", "67200", "67300", "67400", "67600", "67700", "67800", "67900"
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
      sprintf(name, "results/pileup_MC_80x_271036-276811_%s.root",pileup_tags.at(i).c_str());
      pileup_files.push_back(new TFile(name));
      pileup_hists.push_back((TH1D*)pileup_files.at(i)->Get("puweight_dtmc"));
      pileup_weights.push_back(new Double_t(1.0));
      sprintf(name, "puWeight%s", pileup_tags.at(i).c_str());
      sprintf(name1, "puWeight%s/D", pileup_tags.at(i).c_str());
      pileup_branches.push_back(tree_out->Branch(name,pileup_weights.at(i),name1));
    }
  }

  // doDyJetsSigma, correct dyjets mc to match data
  TFile* file_dt_sigma[10];
  TFile* file_mc_sigma[10];
  TH1D* h_dt_met_para_shift[10];
  TH1D* h_mc_met_para_shift[10];
  TH1D* h_met_para_shift_dtmc[10];
  TH1D* h_dt_met_para_sigma[10];
  TH1D* h_dt_met_perp_sigma[10];
  TH1D* h_mc_met_para_sigma[10];
  TH1D* h_mc_met_perp_sigma[10];
  TH1D* h_ratio_met_para_sigma_dtmc[10];
  TH1D* h_ratio_met_perp_sigma_dtmc[10];
  TGraphErrors* gr_dt_met_para_shift[10];
  TGraphErrors* gr_mc_met_para_shift[10];
  TGraphErrors* gr_met_para_shift_dtmc[10];
  TGraphErrors* gr_ratio_met_para_sigma_dtmc[10];
  TGraphErrors* gr_ratio_met_perp_sigma_dtmc[10];

  Float_t llnunu_mt,llnunu_l1_mass, llnunu_l1_pt, llnunu_l1_phi, llnunu_l1_eta;
  Float_t llnunu_l2_px, llnunu_l2_py, llnunu_l2_pt, llnunu_l2_phi;
  Float_t llnunu_l1_l1_pt, llnunu_l1_l1_phi, llnunu_l1_l1_eta, llnunu_l1_l2_pt, llnunu_l1_l2_phi, llnunu_l1_l2_eta;
  Int_t llnunu_l1_l1_pdgId, llnunu_l1_l2_pdgId;

  Float_t llnunu_deltaPhi, llnunu_TuneP_deltaPhi, llnunu_TuneP_mt, llnunu_CosdphiZMet;
  Float_t llnunu_dPTPara, llnunu_dPTParaRel, llnunu_dPTPerp, llnunu_dPTPerpRel, llnunu_metOvSqSET, llnunu_l2_sumEt;
  Float_t llnunu_l1_TuneP_pt, llnunu_l1_TuneP_eta, llnunu_l1_TuneP_phi, llnunu_l1_TuneP_mass; 

  Float_t llnunu_l1_px, llnunu_l1_py, llnunu_l1_pz, llnunu_l1_mt, llnunu_l1_deltaPhi, llnunu_l1_deltaR, llnunu_l1_TuneP_mt, llnunu_l1_TuneP_deltaPhi, llnunu_l1_TuneP_deltaR;
  Float_t llnunu_l1_l1_ptErr, llnunu_l1_l1_TuneP_pt, llnunu_l1_l1_TuneP_ptErr, llnunu_l1_l1_TuneP_eta, llnunu_l1_l1_TuneP_phi, llnunu_l1_l1_px, llnunu_l1_l1_py, llnunu_l1_l1_pz, llnunu_l1_l1_mass;
  Float_t llnunu_l1_l2_ptErr, llnunu_l1_l2_TuneP_pt, llnunu_l1_l2_TuneP_ptErr, llnunu_l1_l2_TuneP_eta, llnunu_l1_l2_TuneP_phi, llnunu_l1_l2_px, llnunu_l1_l2_py, llnunu_l1_l2_pz, llnunu_l1_l2_mass;
  Int_t llnunu_l1_l1_charge, llnunu_l1_l2_charge;

  Float_t llnunu_l1_l1_rapidity, llnunu_l1_l2_rapidity, llnunu_l1_rapidity;

  tree->SetBranchAddress("llnunu_mt",&llnunu_mt);
  tree->SetBranchAddress("llnunu_l1_mass",&llnunu_l1_mass);
  tree->SetBranchAddress("llnunu_l1_pt",&llnunu_l1_pt);
  tree->SetBranchAddress("llnunu_l1_phi",&llnunu_l1_phi);
  tree->SetBranchAddress("llnunu_l1_eta",&llnunu_l1_eta);
  //tree->SetBranchAddress("llnunu_l2_px",&llnunu_l2_px);
  //tree->SetBranchAddress("llnunu_l2_py",&llnunu_l2_py);
  tree->SetBranchAddress("llnunu_l2_pt",&llnunu_l2_pt);
  tree->SetBranchAddress("llnunu_l2_phi",&llnunu_l2_phi);
  tree->SetBranchAddress("llnunu_l1_l1_pt",&llnunu_l1_l1_pt);
  tree->SetBranchAddress("llnunu_l1_l1_phi",&llnunu_l1_l1_phi);
  tree->SetBranchAddress("llnunu_l1_l1_eta",&llnunu_l1_l1_eta);
  tree->SetBranchAddress("llnunu_l1_l2_pt",&llnunu_l1_l2_pt);
  tree->SetBranchAddress("llnunu_l1_l2_phi",&llnunu_l1_l2_phi);
  tree->SetBranchAddress("llnunu_l1_l2_eta",&llnunu_l1_l2_eta);
  tree->SetBranchAddress("llnunu_l1_l1_pdgId",&llnunu_l1_l1_pdgId);
  tree->SetBranchAddress("llnunu_l1_l2_pdgId",&llnunu_l1_l2_pdgId);

  //tree->SetBranchAddress("llnunu_deltaPhi", &llnunu_deltaPhi);
  //tree->SetBranchAddress("llnunu_TuneP_deltaPhi", &llnunu_TuneP_deltaPhi);
  //tree->SetBranchAddress("llnunu_TuneP_mt", &llnunu_TuneP_mt);
  //tree->SetBranchAddress("llnunu_CosdphiZMet", &llnunu_CosdphiZMet);
  //tree->SetBranchAddress("llnunu_dPTPara", &llnunu_dPTPara);
  //tree->SetBranchAddress("llnunu_dPTParaRel", &llnunu_dPTParaRel);
  //tree->SetBranchAddress("llnunu_dPTPerp", &llnunu_dPTPerp);
  //tree->SetBranchAddress("llnunu_dPTPerpRel", &llnunu_dPTPerpRel);
  //tree->SetBranchAddress("llnunu_metOvSqSET", &llnunu_metOvSqSET);
  //tree->SetBranchAddress("llnunu_l2_sumEt", &llnunu_l2_sumEt);
  //tree->SetBranchAddress("llnunu_l1_TuneP_pt", &llnunu_l1_TuneP_pt);
  //tree->SetBranchAddress("llnunu_l1_TuneP_eta", &llnunu_l1_TuneP_eta);
  //tree->SetBranchAddress("llnunu_l1_TuneP_phi", &llnunu_l1_TuneP_phi);

  //tree->SetBranchAddress("llnunu_l1_px", &llnunu_l1_px);
  //tree->SetBranchAddress("llnunu_l1_py", &llnunu_l1_py);
  //tree->SetBranchAddress("llnunu_l1_pz", &llnunu_l1_pz);
  tree->SetBranchAddress("llnunu_l1_mt", &llnunu_l1_mt);
  tree->SetBranchAddress("llnunu_l1_deltaPhi", &llnunu_l1_deltaPhi);
  tree->SetBranchAddress("llnunu_l1_deltaR", &llnunu_l1_deltaR);
  //tree->SetBranchAddress("llnunu_l1_TuneP_mt", &llnunu_l1_TuneP_mt);
  //tree->SetBranchAddress("llnunu_l1_TuneP_deltaPhi", &llnunu_l1_TuneP_deltaPhi);
  //tree->SetBranchAddress("llnunu_l1_TuneP_deltaR", &llnunu_l1_TuneP_deltaR);
  tree->SetBranchAddress("llnunu_l1_l1_charge", &llnunu_l1_l1_charge);
  tree->SetBranchAddress("llnunu_l1_l1_ptErr", &llnunu_l1_l1_ptErr);
  //tree->SetBranchAddress("llnunu_l1_l1_TuneP_pt", &llnunu_l1_l1_TuneP_pt);
  //tree->SetBranchAddress("llnunu_l1_l1_TuneP_ptErr", &llnunu_l1_l1_TuneP_ptErr);
  //tree->SetBranchAddress("llnunu_l1_l1_TuneP_eta", &llnunu_l1_l1_TuneP_eta);
  //tree->SetBranchAddress("llnunu_l1_l1_TuneP_phi", &llnunu_l1_l1_TuneP_phi);
  //tree->SetBranchAddress("llnunu_l1_l1_px", &llnunu_l1_l1_px);
  //tree->SetBranchAddress("llnunu_l1_l1_py", &llnunu_l1_l1_py);
  //tree->SetBranchAddress("llnunu_l1_l1_pz", &llnunu_l1_l1_pz);
  tree->SetBranchAddress("llnunu_l1_l1_mass", &llnunu_l1_l1_mass);
  tree->SetBranchAddress("llnunu_l1_l2_charge", &llnunu_l1_l2_charge);
  tree->SetBranchAddress("llnunu_l1_l2_ptErr", &llnunu_l1_l2_ptErr);
  //tree->SetBranchAddress("llnunu_l1_l2_TuneP_pt", &llnunu_l1_l2_TuneP_pt);
  //tree->SetBranchAddress("llnunu_l1_l2_TuneP_ptErr", &llnunu_l1_l2_TuneP_ptErr);
  //tree->SetBranchAddress("llnunu_l1_l2_TuneP_eta", &llnunu_l1_l2_TuneP_eta);
  //tree->SetBranchAddress("llnunu_l1_l2_TuneP_phi", &llnunu_l1_l2_TuneP_phi);
  //tree->SetBranchAddress("llnunu_l1_l2_px", &llnunu_l1_l2_px);
  //tree->SetBranchAddress("llnunu_l1_l2_py", &llnunu_l1_l2_py);
  //tree->SetBranchAddress("llnunu_l1_l2_pz", &llnunu_l1_l2_pz);
  tree->SetBranchAddress("llnunu_l1_l2_mass", &llnunu_l1_l2_mass);

  //tree->SetBranchAddress("llnunu_l1_l1_rapidity", &llnunu_l1_l1_rapidity);
  //tree->SetBranchAddress("llnunu_l1_l2_rapidity", &llnunu_l1_l2_rapidity);
  tree->SetBranchAddress("llnunu_l1_rapidity", &llnunu_l1_rapidity);

  // old met 
  Float_t llnunu_l2_pt_old, llnunu_l2_phi_old, llnunu_mt_old;
  if (isDyJets|| (isData && correctData) ) {
    std::cout << " -- isDyJets " << std::endl;
    tree_out->Branch("llnunu_mt_old", &llnunu_mt_old, "llnunu_mt_old/F");
    tree_out->Branch("llnunu_l2_pt_old", &llnunu_l2_pt_old, "llnunu_l2_pt_old/F");
    tree_out->Branch("llnunu_l2_phi_old", &llnunu_l2_phi_old, "llnunu_l2_phi_old/F");
  }

  // add rapidity
  Float_t llnunu_l1_l1_TuneP_rapidity, llnunu_l1_l2_TuneP_rapidity, llnunu_l1_TuneP_rapidity;
  if (addZrapidity) {
    tree_out->Branch("llnunu_l1_l1_TuneP_rapidity", &llnunu_l1_l1_TuneP_rapidity, "llnunu_l1_l1_TuneP_rapidity/F");
    tree_out->Branch("llnunu_l1_l2_TuneP_rapidity", &llnunu_l1_l2_TuneP_rapidity, "llnunu_l1_l2_TuneP_rapidity/F");
    tree_out->Branch("llnunu_l1_TuneP_rapidity", &llnunu_l1_TuneP_rapidity, "llnunu_l1_TuneP_rapidity/F");
  }

  //
  Float_t ZPtWeight, ZPtWeight_up, ZPtWeight_dn;
  //Float_t PhiStarWeight, PhiStarWeight_up, PhiStarWeight_dn;
  if (doDyJets && !isData && isDyJets) {
    tree_out->Branch("ZPtWeight",&ZPtWeight,"ZPtWeight/F");
    tree_out->Branch("ZPtWeight_up",&ZPtWeight_up,"ZPtWeight_up/F");
    tree_out->Branch("ZPtWeight_dn",&ZPtWeight_dn,"ZPtWeight_dn/F");
    //tree_out->Branch("PhiStarWeight",&PhiStarWeight,"PhiStarWeight/F");
    //tree_out->Branch("PhiStarWeight_up",&PhiStarWeight_up,"PhiStarWeight_up/F");
    //tree_out->Branch("PhiStarWeight_dn",&PhiStarWeight_dn,"PhiStarWeight_dn/F");
  }

  TFile* fdyzpt;
  TH1D* hdyzpt_dtmc_ratio;
  TH1D* hdyzpt_dtmc_lo_ratio;
  TH1D* hdyzpt_dtmc_ratio_smooth;
  TH1D* hdyzpt_dtmc_lo_ratio_smooth;
  TF1* fcdyzpt_dtmc_ratio;
  TF1* fcdyzpt_dtmc_lo_ratio;
  TH1D* hdyzpt_mc_nlo_lo_ratio;
  TF1* fcdyzpt_mc_nlo_lo_ratio;

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
    //tree->SetBranchAddress("ngenLep", &ngenLep);
    //tree->SetBranchAddress("genLep_eta", genLep_eta);
    //tree->SetBranchAddress("genLep_phi", genLep_phi);
  }


  if (doRecoil && ((isData && correctData) || isDyJets)) {
    // met shift  sigma
    file_dt_sigma[0] = new TFile("SingleEMU_Run2016BCD_PromptReco_met_para_study.root");
    file_dt_sigma[1] = new TFile("SingleEMU_Run2016BCD_PromptReco_met_para_study_mu.root");
    file_dt_sigma[2] = new TFile("SingleEMU_Run2016BCD_PromptReco_met_para_study_el.root");
    if (isDyJetsLO) {
      file_mc_sigma[0] = new TFile("DYJetsToLL_M50_MGMLM_Ext1_NoRecoil_met_para_study.root");
      file_mc_sigma[1] = new TFile("DYJetsToLL_M50_MGMLM_Ext1_NoRecoil_met_para_study_mu.root");
      file_mc_sigma[2] = new TFile("DYJetsToLL_M50_MGMLM_Ext1_NoRecoil_met_para_study_el.root");
    }
    else {
      file_mc_sigma[0] = new TFile("DYJetsToLL_M50_NoRecoil_met_para_study.root");
      file_mc_sigma[1] = new TFile("DYJetsToLL_M50_NoRecoil_met_para_study_mu.root");
      file_mc_sigma[2] = new TFile("DYJetsToLL_M50_NoRecoil_met_para_study_el.root");
    }
    h_dt_met_para_shift[0] = (TH1D*)file_dt_sigma[0]->Get("h_met_para_vs_zpt_mean");
    h_mc_met_para_shift[0] = (TH1D*)file_mc_sigma[0]->Get("h_met_para_vs_zpt_mean");
    h_met_para_shift_dtmc[0] = (TH1D*)h_dt_met_para_shift[0]->Clone("h_met_para_shift_dtmc_all");
    h_met_para_shift_dtmc[0]->Add(h_mc_met_para_shift[0], -1);

    h_dt_met_para_sigma[0] = (TH1D*)file_dt_sigma[0]->Get("h_met_para_vs_zpt_sigma");
    h_dt_met_perp_sigma[0] = (TH1D*)file_dt_sigma[0]->Get("h_met_perp_vs_zpt_sigma");
    h_mc_met_para_sigma[0] = (TH1D*)file_mc_sigma[0]->Get("h_met_para_vs_zpt_sigma");
    h_mc_met_perp_sigma[0] = (TH1D*)file_mc_sigma[0]->Get("h_met_perp_vs_zpt_sigma");

    h_ratio_met_para_sigma_dtmc[0] = (TH1D*)h_dt_met_para_sigma[0]->Clone("h_ratio_met_para_sigma_dtmc_all");
    h_ratio_met_perp_sigma_dtmc[0] = (TH1D*)h_dt_met_perp_sigma[0]->Clone("h_ratio_met_perp_sigma_dtmc_all");
    h_ratio_met_para_sigma_dtmc[0]->Divide(h_mc_met_para_sigma[0]);
    h_ratio_met_perp_sigma_dtmc[0]->Divide(h_mc_met_perp_sigma[0]);

    h_dt_met_para_shift[1] = (TH1D*)file_dt_sigma[1]->Get("h_met_para_vs_zpt_mean");
    h_mc_met_para_shift[1] = (TH1D*)file_mc_sigma[1]->Get("h_met_para_vs_zpt_mean");
    h_met_para_shift_dtmc[1] = (TH1D*)h_dt_met_para_shift[1]->Clone("h_met_para_shift_dtmc_mu");
    h_met_para_shift_dtmc[1]->Add(h_mc_met_para_shift[1], -1);

    h_dt_met_para_sigma[1] = (TH1D*)file_dt_sigma[1]->Get("h_met_para_vs_zpt_sigma");
    h_dt_met_perp_sigma[1] = (TH1D*)file_dt_sigma[1]->Get("h_met_perp_vs_zpt_sigma");
    h_mc_met_para_sigma[1] = (TH1D*)file_mc_sigma[1]->Get("h_met_para_vs_zpt_sigma");
    h_mc_met_perp_sigma[1] = (TH1D*)file_mc_sigma[1]->Get("h_met_perp_vs_zpt_sigma");

    h_ratio_met_para_sigma_dtmc[1] = (TH1D*)h_dt_met_para_sigma[1]->Clone("h_ratio_met_para_sigma_dtmc_mu");
    h_ratio_met_perp_sigma_dtmc[1] = (TH1D*)h_dt_met_perp_sigma[1]->Clone("h_ratio_met_perp_sigma_dtmc_mu");
    h_ratio_met_para_sigma_dtmc[1]->Divide(h_mc_met_para_sigma[1]);
    h_ratio_met_perp_sigma_dtmc[1]->Divide(h_mc_met_perp_sigma[1]);

    h_dt_met_para_shift[2] = (TH1D*)file_dt_sigma[2]->Get("h_met_para_vs_zpt_mean");
    h_mc_met_para_shift[2] = (TH1D*)file_mc_sigma[2]->Get("h_met_para_vs_zpt_mean");
    h_met_para_shift_dtmc[2] = (TH1D*)h_dt_met_para_shift[2]->Clone("h_met_para_shift_dtmc_el");
    h_met_para_shift_dtmc[2]->Add(h_mc_met_para_shift[2], -1);

    h_dt_met_para_sigma[2] = (TH1D*)file_dt_sigma[2]->Get("h_met_para_vs_zpt_sigma");
    h_dt_met_perp_sigma[2] = (TH1D*)file_dt_sigma[2]->Get("h_met_perp_vs_zpt_sigma");
    h_mc_met_para_sigma[2] = (TH1D*)file_mc_sigma[2]->Get("h_met_para_vs_zpt_sigma");
    h_mc_met_perp_sigma[2] = (TH1D*)file_mc_sigma[2]->Get("h_met_perp_vs_zpt_sigma");

    h_ratio_met_para_sigma_dtmc[2] = (TH1D*)h_dt_met_para_sigma[2]->Clone("h_ratio_met_para_sigma_dtmc_el");
    h_ratio_met_perp_sigma_dtmc[2] = (TH1D*)h_dt_met_perp_sigma[2]->Clone("h_ratio_met_perp_sigma_dtmc_el");
    h_ratio_met_para_sigma_dtmc[2]->Divide(h_mc_met_para_sigma[2]);
    h_ratio_met_perp_sigma_dtmc[2]->Divide(h_mc_met_perp_sigma[2]);

    // smooth functions
    h_dt_met_para_shift[0]->SetName("h_dt_met_para_shift_all");
    h_mc_met_para_shift[0]->SetName("h_mc_met_para_shift_all");
    h_dt_met_para_shift[1]->SetName("h_dt_met_para_shift_mu");
    h_mc_met_para_shift[1]->SetName("h_mc_met_para_shift_mu");
    h_dt_met_para_shift[2]->SetName("h_dt_met_para_shift_el");
    h_mc_met_para_shift[2]->SetName("h_mc_met_para_shift_el");
    h_dt_met_para_shift[3] = (TH1D*)h_dt_met_para_shift[0]->Clone("h_dt_met_para_shift_all_smooth");
    h_mc_met_para_shift[3] = (TH1D*)h_mc_met_para_shift[0]->Clone("h_mc_met_para_shift_all_smooth");
    h_dt_met_para_shift[4] = (TH1D*)h_dt_met_para_shift[1]->Clone("h_dt_met_para_shift_mu_smooth");
    h_mc_met_para_shift[4] = (TH1D*)h_mc_met_para_shift[1]->Clone("h_mc_met_para_shift_mu_smooth");
    h_dt_met_para_shift[5] = (TH1D*)h_dt_met_para_shift[2]->Clone("h_dt_met_para_shift_el_smooth");
    h_mc_met_para_shift[5] = (TH1D*)h_mc_met_para_shift[2]->Clone("h_mc_met_para_shift_el_smooth");
    h_dt_met_para_shift[3]->Smooth();
    h_mc_met_para_shift[3]->Smooth();
    h_dt_met_para_shift[4]->Smooth();
    h_mc_met_para_shift[4]->Smooth();
    h_dt_met_para_shift[5]->Smooth();
    h_mc_met_para_shift[5]->Smooth();
 
    gr_dt_met_para_shift[3] = new TGraphErrors(h_dt_met_para_shift[3]);
    gr_mc_met_para_shift[3] = new TGraphErrors(h_mc_met_para_shift[3]);
    gr_dt_met_para_shift[4] = new TGraphErrors(h_dt_met_para_shift[4]);
    gr_mc_met_para_shift[4] = new TGraphErrors(h_mc_met_para_shift[4]);
    gr_dt_met_para_shift[5] = new TGraphErrors(h_dt_met_para_shift[5]);
    gr_mc_met_para_shift[5] = new TGraphErrors(h_mc_met_para_shift[5]);

    h_met_para_shift_dtmc[3] = (TH1D*)h_met_para_shift_dtmc[0]->Clone("h_met_para_shift_dtmc_all_smooth");
    h_met_para_shift_dtmc[4] = (TH1D*)h_met_para_shift_dtmc[1]->Clone("h_met_para_shift_dtmc_mu_smooth");
    h_met_para_shift_dtmc[5] = (TH1D*)h_met_para_shift_dtmc[2]->Clone("h_met_para_shift_dtmc_el_smooth");
    h_met_para_shift_dtmc[3]->Smooth();
    h_met_para_shift_dtmc[4]->Smooth();
    h_met_para_shift_dtmc[5]->Smooth();

    gr_met_para_shift_dtmc[3] = new TGraphErrors(h_met_para_shift_dtmc[3]);
    gr_met_para_shift_dtmc[4] = new TGraphErrors(h_met_para_shift_dtmc[4]);
    gr_met_para_shift_dtmc[5] = new TGraphErrors(h_met_para_shift_dtmc[5]);
  
    h_ratio_met_para_sigma_dtmc[3] = (TH1D*)h_ratio_met_para_sigma_dtmc[0]->Clone("h_ratio_met_para_sigma_dtmc_all_smooth");
    h_ratio_met_para_sigma_dtmc[4] = (TH1D*)h_ratio_met_para_sigma_dtmc[1]->Clone("h_ratio_met_para_sigma_dtmc_mu_smooth");
    h_ratio_met_para_sigma_dtmc[5] = (TH1D*)h_ratio_met_para_sigma_dtmc[2]->Clone("h_ratio_met_para_sigma_dtmc_el_smooth");
    h_ratio_met_para_sigma_dtmc[3]->Smooth();
    h_ratio_met_para_sigma_dtmc[4]->Smooth();
    h_ratio_met_para_sigma_dtmc[5]->Smooth();
  
    h_ratio_met_perp_sigma_dtmc[3] = (TH1D*)h_ratio_met_perp_sigma_dtmc[0]->Clone("h_ratio_met_perp_sigma_dtmc_all_smooth");
    h_ratio_met_perp_sigma_dtmc[4] = (TH1D*)h_ratio_met_perp_sigma_dtmc[1]->Clone("h_ratio_met_perp_sigma_dtmc_mu_smooth");
    h_ratio_met_perp_sigma_dtmc[5] = (TH1D*)h_ratio_met_perp_sigma_dtmc[2]->Clone("h_ratio_met_perp_sigma_dtmc_el_smooth");
    h_ratio_met_perp_sigma_dtmc[3]->Smooth();
    h_ratio_met_perp_sigma_dtmc[4]->Smooth();
    h_ratio_met_perp_sigma_dtmc[5]->Smooth();
 
    gr_ratio_met_para_sigma_dtmc[3] = new TGraphErrors(h_ratio_met_para_sigma_dtmc[3]);
    gr_ratio_met_para_sigma_dtmc[4] = new TGraphErrors(h_ratio_met_para_sigma_dtmc[4]);
    gr_ratio_met_para_sigma_dtmc[5] = new TGraphErrors(h_ratio_met_para_sigma_dtmc[5]);
    gr_ratio_met_perp_sigma_dtmc[3] = new TGraphErrors(h_ratio_met_perp_sigma_dtmc[3]);
    gr_ratio_met_perp_sigma_dtmc[4] = new TGraphErrors(h_ratio_met_perp_sigma_dtmc[4]);
    gr_ratio_met_perp_sigma_dtmc[5] = new TGraphErrors(h_ratio_met_perp_sigma_dtmc[5]);
  }

  if (doDyJets && !isData && isDyJets){
  
    //fdyzpt = new TFile("dyjets_zpt_weight_lo_nlo.root");
    fdyzpt = new TFile("dyjets_zpt_weight_lo_nlo_sel.root");
    hdyzpt_dtmc_ratio = (TH1D*)fdyzpt->Get("hdyzpt_dtmc_ratio"); 
    hdyzpt_dtmc_lo_ratio = (TH1D*)fdyzpt->Get("hdyzpt_dtmc_lo_ratio"); 
    hdyzpt_dtmc_ratio_smooth = (TH1D*)fdyzpt->Get("hdyzpt_dtmc_ratio_smooth"); 
    hdyzpt_dtmc_lo_ratio_smooth = (TH1D*)fdyzpt->Get("hdyzpt_dtmc_lo_ratio_smooth"); 
    fcdyzpt_dtmc_ratio = (TF1*)fdyzpt->Get("fcdyzpt_dtmc_ratio"); 
    fcdyzpt_dtmc_lo_ratio = (TF1*)fdyzpt->Get("fcdyzpt_dtmc_lo_ratio"); 
    hdyzpt_mc_nlo_lo_ratio = (TH1D*)fdyzpt->Get("hdyzpt_mc_nlo_lo_ratio"); 
    fcdyzpt_mc_nlo_lo_ratio = (TF1*)fdyzpt->Get("fcdyzpt_mc_nlo_lo_ratio"); 
  }

  // correct Muon Pt
  KalmanMuonCalibrator* muCalib;

  if (correctMuonPt) {
    if (isData) muCalib = new KalmanMuonCalibrator("data/DATA_80X_13TeV.root");
    else muCalib = new KalmanMuonCalibrator("data/MC_80X_13TeV.root");
  } 

  //
  int n_interval = 5000;
  Int_t nTrueInt; 
  if (!isData) tree->SetBranchAddress( "nTrueInt", &nTrueInt );

  for (int i=0; i<(int)tree->GetEntries(); i++){
    tree->GetEntry(i);

    if ( i%n_interval == 0 ) {
      std::cout << "Event " << i << "   " << std::endl;
    }

    // correctMuonPt
    if (correctMuonPt && abs(llnunu_l1_l1_pdgId)==13&&abs(llnunu_l1_l2_pdgId)==13 && llnunu_l1_l1_pt>2. && llnunu_l1_l1_pt<200. ) {
      llnunu_l1_l1_pt = (Float_t)muCalib->getCorrectedPt(double(llnunu_l1_l1_pt), double(llnunu_l1_l1_eta), double(llnunu_l1_l1_phi), double(llnunu_l1_l1_charge)); 
      llnunu_l1_l2_pt = (Float_t)muCalib->getCorrectedPt(double(llnunu_l1_l2_pt), double(llnunu_l1_l2_eta), double(llnunu_l1_l2_phi), double(llnunu_l1_l2_charge)); 
      //llnunu_l1_l1_TuneP_pt = (Float_t)muCalib->getCorrectedPt(double(llnunu_l1_l1_TuneP_pt), double(llnunu_l1_l1_TuneP_eta), double(llnunu_l1_l1_TuneP_phi), double(llnunu_l1_l1_charge)); 
      //llnunu_l1_l2_TuneP_pt = (Float_t)muCalib->getCorrectedPt(double(llnunu_l1_l2_TuneP_pt), double(llnunu_l1_l2_TuneP_eta), double(llnunu_l1_l2_TuneP_phi), double(llnunu_l1_l2_charge)); 
      if (!isData) {
        llnunu_l1_l1_pt = (Float_t)muCalib->smear(double(llnunu_l1_l1_pt), double(llnunu_l1_l1_eta));
        llnunu_l1_l2_pt = (Float_t)muCalib->smear(double(llnunu_l1_l2_pt), double(llnunu_l1_l2_eta));
        //llnunu_l1_l1_TuneP_pt = (Float_t)muCalib->smear(double(llnunu_l1_l1_TuneP_pt), double(llnunu_l1_l1_TuneP_eta));
        //llnunu_l1_l2_TuneP_pt = (Float_t)muCalib->smear(double(llnunu_l1_l2_TuneP_pt), double(llnunu_l1_l2_TuneP_eta));
      }
      llnunu_l1_l1_ptErr = llnunu_l1_l1_pt*(Float_t)muCalib->getCorrectedError(double(llnunu_l1_l1_pt), double(llnunu_l1_l1_eta), double(llnunu_l1_l1_ptErr/llnunu_l1_l1_pt));
      llnunu_l1_l2_ptErr = llnunu_l1_l2_pt*(Float_t)muCalib->getCorrectedError(double(llnunu_l1_l2_pt), double(llnunu_l1_l2_eta), double(llnunu_l1_l2_ptErr/llnunu_l1_l2_pt));
      //llnunu_l1_l1_TuneP_ptErr = llnunu_l1_l1_TuneP_pt*(Float_t)muCalib->getCorrectedError(double(llnunu_l1_l1_TuneP_pt), double(llnunu_l1_l1_TuneP_eta), double(llnunu_l1_l1_TuneP_ptErr/llnunu_l1_l1_TuneP_pt));
      //llnunu_l1_l2_TuneP_ptErr = llnunu_l1_l2_TuneP_pt*(Float_t)muCalib->getCorrectedError(double(llnunu_l1_l2_TuneP_pt), double(llnunu_l1_l2_TuneP_eta), double(llnunu_l1_l2_TuneP_ptErr/llnunu_l1_l2_TuneP_pt));
      TLorentzVector l1v, l2v; 
      l1v.SetPtEtaPhiM(llnunu_l1_l1_pt, llnunu_l1_l1_eta, llnunu_l1_l1_phi, llnunu_l1_l1_mass);
      l2v.SetPtEtaPhiM(llnunu_l1_l2_pt, llnunu_l1_l2_eta, llnunu_l1_l2_phi, llnunu_l1_l2_mass);
      //TLorentzVector tpl1v, tpl2v;
      //tpl1v.SetPtEtaPhiM(llnunu_l1_l1_TuneP_pt, llnunu_l1_l1_TuneP_eta, llnunu_l1_l1_TuneP_phi, llnunu_l1_l1_mass);
      //tpl2v.SetPtEtaPhiM(llnunu_l1_l2_TuneP_pt, llnunu_l1_l2_TuneP_eta, llnunu_l1_l2_TuneP_phi, llnunu_l1_l2_mass);
      TLorentzVector zv = l1v+l2v;
      //TLorentzVector tpzv = tpl1v+tpl2v;

      //llnunu_l1_l1_px = (Float_t)l1v.Px();
      //llnunu_l1_l1_py = (Float_t)l1v.Py();
      //llnunu_l1_l1_pz = (Float_t)l1v.Pz();
      llnunu_l1_l1_rapidity = (Float_t)l1v.Rapidity();
      //llnunu_l1_l2_px = (Float_t)l2v.Px();
      //llnunu_l1_l2_py = (Float_t)l2v.Py();
      //llnunu_l1_l2_pz = (Float_t)l2v.Pz();
      llnunu_l1_l2_rapidity = (Float_t)l2v.Rapidity();
      //llnunu_l1_l1_TuneP_rapidity = (Float_t)tpl1v.Rapidity();
      //llnunu_l1_l2_TuneP_rapidity = (Float_t)tpl2v.Rapidity();
      llnunu_l1_pt = (Float_t)zv.Pt();
      //llnunu_l1_px = (Float_t)zv.Px();
      //llnunu_l1_py = (Float_t)zv.Py();
      //llnunu_l1_pz = (Float_t)zv.Pz();
      llnunu_l1_eta = (Float_t)zv.Eta();
      llnunu_l1_phi = (Float_t)zv.Phi();
      llnunu_l1_rapidity = (Float_t)zv.Rapidity();
      llnunu_l1_deltaPhi = (Float_t)l1v.DeltaPhi(l2v);
      llnunu_l1_deltaR = (Float_t)l1v.DeltaR(l2v);
      llnunu_l1_mt = (Float_t)zv.Mt();
      llnunu_l1_mass = (Float_t)zv.M();
      //llnunu_l1_TuneP_pt = (Float_t)tpzv.Pt();
      //llnunu_l1_TuneP_eta = (Float_t)tpzv.Eta();
      //llnunu_l1_TuneP_phi = (Float_t)tpzv.Phi();
      //llnunu_l1_TuneP_rapidity = (Float_t)tpzv.Rapidity();
      //llnunu_l1_TuneP_deltaPhi = (Float_t)tpl1v.DeltaPhi(tpl2v);
      //llnunu_l1_TuneP_deltaR = (Float_t)tpl1v.DeltaR(tpl2v);
      //llnunu_l1_TuneP_mt = (Float_t)tpzv.Mt();
      //llnunu_l1_TuneP_mass = (Float_t)tpzv.M();

      TVector2 vec_met;
      llnunu_l2_px = llnunu_l2_pt*cos(llnunu_l2_phi);
      llnunu_l2_py = llnunu_l2_pt*sin(llnunu_l2_phi);
      vec_met.Set(llnunu_l2_px, llnunu_l2_py);

      llnunu_deltaPhi = TVector2::Phi_mpi_pi(llnunu_l2_phi-llnunu_l1_phi);
      //llnunu_TuneP_deltaPhi = TVector2::Phi_mpi_pi(llnunu_l2_phi-llnunu_l1_TuneP_phi);
      //llnunu_CosdphiZMet = TMath::Cos(llnunu_deltaPhi);

      //llnunu_dPTPara = fabs(llnunu_l1_pt+llnunu_l2_pt*cos(llnunu_l2_phi-llnunu_l1_phi));
      //llnunu_dPTParaRel = llnunu_dPTPara/llnunu_l1_pt;
      //llnunu_dPTPerp = fabs(llnunu_l1_pt+llnunu_l2_pt*sin(llnunu_l2_phi-llnunu_l1_phi));
      //llnunu_dPTPerpRel = llnunu_dPTPerp/llnunu_l1_pt;
      //llnunu_metOvSqSET = llnunu_l2_pt/llnunu_l2_sumEt;

      Float_t et1 = TMath::Sqrt(llnunu_l1_mass*llnunu_l1_mass + llnunu_l1_pt*llnunu_l1_pt);
      Float_t et2 = TMath::Sqrt(llnunu_l1_mass*llnunu_l1_mass + llnunu_l2_pt*llnunu_l2_pt);
      llnunu_mt = TMath::Sqrt(2.0*llnunu_l1_mass*llnunu_l1_mass + 2.0* (et1*et2 - llnunu_l1_pt*cos(llnunu_l1_phi)*llnunu_l2_px - llnunu_l1_pt*sin(llnunu_l1_phi)*llnunu_l2_py));

      //et1 = TMath::Sqrt(llnunu_l1_TuneP_mass*llnunu_l1_TuneP_mass + llnunu_l1_TuneP_pt*llnunu_l1_TuneP_pt);
      //et2 = TMath::Sqrt(llnunu_l1_TuneP_mass*llnunu_l1_TuneP_mass + llnunu_l2_pt*llnunu_l2_pt);
      //llnunu_TuneP_mt = TMath::Sqrt(2.0*llnunu_l1_TuneP_mass*llnunu_l1_TuneP_mass + 2.0* (et1*et2 - llnunu_l1_TuneP_pt*cos(llnunu_l1_TuneP_phi)*llnunu_l2_px - llnunu_l1_TuneP_pt*sin(llnunu_l1_TuneP_phi)*llnunu_l2_py));
    }

    if (!isData) {
      for (int j=0; j<(int)pileup_tags.size(); j++){
        *(pileup_weights.at(j)) = pileup_hists.at(j)->GetBinContent(pileup_hists.at(j)->FindBin(nTrueInt));
      }
    }
    Float_t met_para = llnunu_l2_pt*cos(llnunu_l2_phi-llnunu_l1_phi);
    Float_t met_perp = llnunu_l2_pt*sin(llnunu_l2_phi-llnunu_l1_phi);
    if ( doRecoil && ((isData &&correctData)||isDyJets) ) {
      if (abs(llnunu_l1_l1_pdgId)==13&&abs(llnunu_l1_l2_pdgId)==13) {
        Int_t idd=1;
        if (doRecoilUseSmooth) idd=4;
        h_dt_met_para_shift[6] = h_dt_met_para_shift[idd];
        h_mc_met_para_shift[6] = h_mc_met_para_shift[idd];
        h_met_para_shift_dtmc[6] = h_met_para_shift_dtmc[idd];
        h_ratio_met_para_sigma_dtmc[6] = h_ratio_met_para_sigma_dtmc[idd];
        h_ratio_met_perp_sigma_dtmc[6] = h_ratio_met_perp_sigma_dtmc[idd];
        gr_dt_met_para_shift[6] = gr_dt_met_para_shift[4];
        gr_mc_met_para_shift[6] = gr_mc_met_para_shift[4];
        gr_met_para_shift_dtmc[6] = gr_met_para_shift_dtmc[4];
        gr_ratio_met_para_sigma_dtmc[6] = gr_ratio_met_para_sigma_dtmc[4];
        gr_ratio_met_perp_sigma_dtmc[6] = gr_ratio_met_perp_sigma_dtmc[4];
      }
      else if (abs(llnunu_l1_l1_pdgId)==11&&abs(llnunu_l1_l2_pdgId)==11) {
        Int_t idd=2;
        if (doRecoilUseSmooth) idd=5;
        h_dt_met_para_shift[6] = h_dt_met_para_shift[idd];
        h_mc_met_para_shift[6] = h_mc_met_para_shift[idd];
        h_met_para_shift_dtmc[6] = h_met_para_shift_dtmc[idd];
        h_ratio_met_para_sigma_dtmc[6] = h_ratio_met_para_sigma_dtmc[idd];
        h_ratio_met_perp_sigma_dtmc[6] = h_ratio_met_perp_sigma_dtmc[idd];
        gr_dt_met_para_shift[6] = gr_dt_met_para_shift[5];
        gr_mc_met_para_shift[6] = gr_mc_met_para_shift[5];
        gr_met_para_shift_dtmc[6] = gr_met_para_shift_dtmc[5];
        gr_ratio_met_para_sigma_dtmc[6] = gr_ratio_met_para_sigma_dtmc[5];
        gr_ratio_met_perp_sigma_dtmc[6] = gr_ratio_met_perp_sigma_dtmc[5];
      }
      else {
        Int_t idd=0;
        if (doRecoilUseSmooth) idd=3;
        h_dt_met_para_shift[6] = h_dt_met_para_shift[idd];
        h_mc_met_para_shift[6] = h_mc_met_para_shift[idd];
        h_met_para_shift_dtmc[6] = h_met_para_shift_dtmc[idd];
        h_ratio_met_para_sigma_dtmc[6] = h_ratio_met_para_sigma_dtmc[idd];
        h_ratio_met_perp_sigma_dtmc[6] = h_ratio_met_perp_sigma_dtmc[idd];
        gr_dt_met_para_shift[6] = gr_dt_met_para_shift[3];
        gr_mc_met_para_shift[6] = gr_mc_met_para_shift[3];
        gr_met_para_shift_dtmc[6] = gr_met_para_shift_dtmc[3];
        gr_ratio_met_para_sigma_dtmc[6] = gr_ratio_met_para_sigma_dtmc[3];
        gr_ratio_met_perp_sigma_dtmc[6] = gr_ratio_met_perp_sigma_dtmc[3];
      }

      if (isData && correctData) {
        if (doRecoilUseGraph) {
          met_para -= gr_dt_met_para_shift[6]->Eval(llnunu_l1_pt);
        }
        else {
          met_para -= h_dt_met_para_shift[6]->GetBinContent(h_dt_met_para_shift[6]->FindBin(llnunu_l1_pt));
        }
      }
      else if (isDyJets && correctData) {
        if (doRecoilUseGraph) {
          met_para -= gr_mc_met_para_shift[6]->Eval(llnunu_l1_pt);
        }
        else {
          met_para -= h_mc_met_para_shift[6]->GetBinContent(h_mc_met_para_shift[6]->FindBin(llnunu_l1_pt));
        }
      }
      else if (isDyJets && !correctData) {
        if (doRecoilUseGraph) {
          met_para += gr_met_para_shift_dtmc[6]->Eval(llnunu_l1_pt);
        }
        else {
          met_para += h_met_para_shift_dtmc[6]->GetBinContent(h_met_para_shift_dtmc[6]->FindBin(llnunu_l1_pt));
        }
      }

      if (!isData && isDyJets){
        if (doRecoilUseGraph) {
          if (correctData) {
            met_para *= gr_ratio_met_para_sigma_dtmc[6]->Eval(llnunu_l1_pt);
            met_perp *= gr_ratio_met_perp_sigma_dtmc[6]->Eval(llnunu_l1_pt);
          }
          else {
            met_para = (met_para-gr_dt_met_para_shift[6]->Eval(llnunu_l1_pt))*gr_ratio_met_para_sigma_dtmc[6]->Eval(llnunu_l1_pt) + gr_dt_met_para_shift[6]->Eval(llnunu_l1_pt);
            met_perp *= gr_ratio_met_perp_sigma_dtmc[6]->Eval(llnunu_l1_pt);
          }
        } 
        else {
          if (correctData) {
            met_para *= h_ratio_met_para_sigma_dtmc[6]->GetBinContent(h_ratio_met_para_sigma_dtmc[6]->FindBin(llnunu_l1_pt));
            met_perp *= h_ratio_met_perp_sigma_dtmc[6]->GetBinContent(h_ratio_met_perp_sigma_dtmc[6]->FindBin(llnunu_l1_pt));
          }
          else {
            met_para = (met_para-h_dt_met_para_shift[6]->GetBinContent(h_dt_met_para_shift[6]->FindBin(llnunu_l1_pt)))
                     * h_ratio_met_para_sigma_dtmc[6]->GetBinContent(h_ratio_met_para_sigma_dtmc[6]->FindBin(llnunu_l1_pt)) 
                     + h_dt_met_para_shift[6]->GetBinContent(h_dt_met_para_shift[6]->FindBin(llnunu_l1_pt));
            met_perp *= h_ratio_met_perp_sigma_dtmc[6]->GetBinContent(h_ratio_met_perp_sigma_dtmc[6]->FindBin(llnunu_l1_pt));
          }
        }
      }

      // recalculate vars

      // store some essential old vars
      llnunu_l2_pt_old = llnunu_l2_pt;
      llnunu_l2_phi_old = llnunu_l2_phi;
      llnunu_mt_old = llnunu_mt;

      // recalculate
      llnunu_l2_px = met_para*cos(llnunu_l1_phi)-met_perp*sin(llnunu_l1_phi);
      llnunu_l2_py = met_para*sin(llnunu_l1_phi)+met_perp*cos(llnunu_l1_phi);
      TVector2 vec_met;
      vec_met.Set(llnunu_l2_px, llnunu_l2_py);
      llnunu_l2_pt = vec_met.Mod();
      llnunu_l2_phi = TVector2::Phi_mpi_pi(vec_met.Phi());
      
      llnunu_deltaPhi = TVector2::Phi_mpi_pi(llnunu_l2_phi-llnunu_l1_phi);
      //llnunu_TuneP_deltaPhi = TVector2::Phi_mpi_pi(llnunu_l2_phi-llnunu_l1_TuneP_phi);
      //llnunu_CosdphiZMet = TMath::Cos(llnunu_deltaPhi);
          
      //llnunu_dPTPara = fabs(llnunu_l1_pt+llnunu_l2_pt*cos(llnunu_l2_phi-llnunu_l1_phi)); 
      //llnunu_dPTParaRel = llnunu_dPTPara/llnunu_l1_pt;
      //llnunu_dPTPerp = fabs(llnunu_l1_pt+llnunu_l2_pt*sin(llnunu_l2_phi-llnunu_l1_phi));
      //llnunu_dPTPerpRel = llnunu_dPTPerp/llnunu_l1_pt;
      //llnunu_metOvSqSET = llnunu_l2_pt/llnunu_l2_sumEt;

      Float_t et1 = TMath::Sqrt(llnunu_l1_mass*llnunu_l1_mass + llnunu_l1_pt*llnunu_l1_pt);
      Float_t et2 = TMath::Sqrt(llnunu_l1_mass*llnunu_l1_mass + llnunu_l2_pt*llnunu_l2_pt);
      llnunu_mt = TMath::Sqrt(2.0*llnunu_l1_mass*llnunu_l1_mass + 2.0* (et1*et2 - llnunu_l1_pt*cos(llnunu_l1_phi)*llnunu_l2_px - llnunu_l1_pt*sin(llnunu_l1_phi)*llnunu_l2_py));

      //et1 = TMath::Sqrt(llnunu_l1_TuneP_mass*llnunu_l1_TuneP_mass + llnunu_l1_TuneP_pt*llnunu_l1_TuneP_pt);
      //et2 = TMath::Sqrt(llnunu_l1_TuneP_mass*llnunu_l1_TuneP_mass + llnunu_l2_pt*llnunu_l2_pt);
      //llnunu_TuneP_mt = TMath::Sqrt(2.0*llnunu_l1_TuneP_mass*llnunu_l1_TuneP_mass + 2.0* (et1*et2 - llnunu_l1_TuneP_pt*cos(llnunu_l1_TuneP_phi)*llnunu_l2_px - llnunu_l1_TuneP_pt*sin(llnunu_l1_TuneP_phi)*llnunu_l2_py));

    }

    //
    ZPtWeight=1.0;
    ZPtWeight_up=1.0;
    ZPtWeight_dn=1.0;
    //PhiStarWeight=1.0;
    //PhiStarWeight_up=1.0;
    //PhiStarWeight_dn=1.0;

    if (doDyJets && !isData&&isDyJets) {
      
      // zpt weight
      Int_t zptBin=0;
      if (ngenZ>0) {
        zptBin = hdyzpt_dtmc_ratio->FindBin(genZ_pt[0]);
        if (genZ_pt[0]>1000) zptBin = hdyzpt_dtmc_ratio->FindBin(999);
      }
      else {
        zptBin = hdyzpt_dtmc_ratio->FindBin(llnunu_l1_pt);
        if (llnunu_l1_pt>1000) zptBin = hdyzpt_dtmc_ratio->FindBin(999);
      }
      if (doDyJetsUseSmooth) {
        ZPtWeight = hdyzpt_dtmc_ratio_smooth->GetBinContent(zptBin);
      }
      else if (doDyJetsUseFunction) {
        if (ngenZ>0) ZPtWeight = fcdyzpt_dtmc_ratio->Eval(genZ_pt[0]);
        else ZPtWeight = fcdyzpt_dtmc_ratio->Eval(llnunu_l1_pt);
      }
      else { 
        ZPtWeight = hdyzpt_dtmc_ratio->GetBinContent(zptBin);
      }
      ZPtWeight_up = ZPtWeight+0.5*hdyzpt_dtmc_ratio->GetBinError(zptBin);
      ZPtWeight_dn = ZPtWeight-0.5*hdyzpt_dtmc_ratio->GetBinError(zptBin);
      
      if (isDyJetsLO) {
        // for LO ZJets samples
        zptBin=0;
        if (ngenZ>0) {
          zptBin = hdyzpt_dtmc_lo_ratio->FindBin(genZ_pt[0]);
          if (genZ_pt[0]>1000) zptBin = hdyzpt_dtmc_lo_ratio->FindBin(999);
        }
        else {
          zptBin = hdyzpt_dtmc_lo_ratio->FindBin(llnunu_l1_pt);
          if (llnunu_l1_pt>1000) zptBin = hdyzpt_dtmc_lo_ratio->FindBin(999);
        }
        if (doDyJetsLOUseSmooth) {
          ZPtWeight = hdyzpt_dtmc_lo_ratio_smooth->GetBinContent(zptBin);
        }
        else if (doDyJetsLOUseFunction) {
        //if (doDyJetsLOUseFunction) {
          //if (ngenZ>0) ZPtWeight *= fcdyzpt_dtmc_lo_ratio->Eval(genZ_pt[0]);
          //else ZPtWeight *= fcdyzpt_dtmc_lo_ratio->Eval(llnunu_l1_pt);
          if (ngenZ>0) ZPtWeight *= fcdyzpt_mc_nlo_lo_ratio->Eval(genZ_pt[0]);
          else ZPtWeight *= fcdyzpt_mc_nlo_lo_ratio->Eval(llnunu_l1_pt);
        }
        else {
          //ZPtWeight *= hdyzpt_dtmc_lo_ratio->GetBinContent(zptBin);
          ZPtWeight *= hdyzpt_mc_nlo_lo_ratio->GetBinContent(zptBin);
        }

        ZPtWeight_up = ZPtWeight+0.5*hdyzpt_dtmc_lo_ratio->GetBinError(zptBin);
        ZPtWeight_dn = ZPtWeight-0.5*hdyzpt_dtmc_lo_ratio->GetBinError(zptBin);

      }
    }

    if (addZjetsNewGenWeight and !isData and isDyJets) {
      if (!isDyJetsLO) {
        ZJetsGenWeight = genWeight/ZJetsNLOSumWeights*ZJetsNLOSumEvents/(ZJetsNLOSumEvents+ZJetsLOSumEvents);
      }
      else {
        ZJetsGenWeight = genWeight/ZJetsLOSumWeights*ZJetsLOSumEvents/(ZJetsNLOSumEvents+ZJetsLOSumEvents);
      }
    }
    if (addZrapidity) {
      TLorentzVector z4vec;
      z4vec.SetPtEtaPhiM(llnunu_l1_pt,llnunu_l1_eta,llnunu_l1_phi,llnunu_l1_mass);
      llnunu_l1_rapidity = z4vec.Rapidity();
    }

    tree_out->Fill();
  }

  foutput->cd();
  tree_out->Write();
  foutput->Close();
  finput->Close();

  return 0;

}



