// Hengne Li, 2016 @ CERN, initial version.
 
// root
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TEntryList.h"
#include "TLorentzVector.h"


// std
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <time.h>

// other
#include "tdrstyle.h"
#include "PParameterReader.h"
#include "KalmanMuonCalibrator.h"



//======================================================
// ╔═╗╦  ╔═╗╔╗ ╔═╗╦    ╦  ╦╔═╗╦═╗╦╔═╗╔╗ ╦  ╔═╗╔═╗ 
// ║ ╦║  ║ ║╠╩╗╠═╣║    ╚╗╔╝╠═╣╠╦╝║╠═╣╠╩╗║  ║╣ ╚═╗ 
// ╚═╝╩═╝╚═╝╚═╝╩ ╩╩═╝   ╚╝ ╩ ╩╩╚═╩╩ ╩╚═╝╩═╝╚═╝╚═╝ 
//======================================================


// Intput output root files

std::string _file_in_name;
std::string _file_out_name;
std::string _file_config_name;

TFile* _file_in;
TFile* _file_out;


// input output root trees
TTree* _tree_in;
TTree* _tree_out;

// selected entry list
TEntryList* _selected_entries;

// output plots file
std::string _plots_out_name;

// if DYJets samples
bool _isDyJets, _isDyJetsLO;

// strings for temporary usage
char name[3000], name1[3000], name2[3000];


//======================================================
// ╔═╗╔═╗╔╗╔╔═╗╦╔═╗  ╔═╗╔═╗╦═╗╔═╗╔╦╗╔═╗╔╦╗╔═╗╦═╗╔═╗
// ║  ║ ║║║║╠╣ ║║ ╦  ╠═╝╠═╣╠╦╝╠═╣║║║║╣  ║ ║╣ ╠╦╝╚═╗
// ╚═╝╚═╝╝╚╝╚  ╩╚═╝  ╩  ╩ ╩╩╚═╩ ╩╩ ╩╚═╝ ╩ ╚═╝╩╚═╚═╝
//======================================================


//====================
// general paramters
//====================

// debug
bool _debug=false;

// starting entry number 
int _n_start = 0;

// number of entries to be run, -1 means all entries
int _n_test = -1;

// n entries interval between which print once entry number
int _n_interval = 1000;

// use slimmed tree or not
bool _useLightTree = true;

// tree selection string
std::string _selection = "(1)";

// store old branches 
bool _storeOldBranches = false;


//=========================
// add PU weights
//=========================
bool _addPUWeights = true;
// PU input files directory
std::string _PUInputDir;
// PU input tags, "puWeight<tag>" branches will be added
std::vector<std::string> _PUTags;
// PU input root files, one-to-one correspnding to the tags above
std::vector<std::string> _PUInputFileNames;
// PU weight hist name
std::string _PUWeightHistName;

// not from config file:
std::vector<TFile*> _PUFiles;
std::vector<TH1D*> _PUHists;
std::vector<Float_t*> _PUWeights;


//==========================
// recalibrate muon pt
//==========================
bool _doMuonPtRecalib = false;
// recalibrator input files
std::string _MuonPtRecalibInputForData, _MuonPtRecalibInputForMC;

// not from config file:
KalmanMuonCalibrator* _muCalib;


//========================
// Add DYJet gen reweight 
//========================
// Add DYJet gen reweight according to 
// unfolded 2015 data precision measurement results 
// will also reweight LO to NLO MC samples to gain statistics.

// default will use histograms for both NLO and LO reweight
bool _addDyZPtWeight = true;

// for NLO MC use function parametrization instead of hist.
bool _addDyZPtWeightUseFunction = true;

// for LO MC to NLO, also use function  
bool _addDyZPtWeightLOUseFunction = true;

// Input root file
std::string _DyZPtWeightInputFileName = "data/zptweight/dyjets_zpt_weight_lo_nlo_sel.root";

// add new DY jets gen-weights for mergeing NLO and LO
// to gain statistics
bool _addDyNewGenWeight = true;



// not from config file
TFile* _fdyzpt;
TH1D*  _hdyzpt_dtmc_ratio;
TF1*   _fcdyzpt_dtmc_ratio;
TH1D*  _hdyzpt_mc_nlo_lo_ratio;
TF1*   _fcdyzpt_mc_nlo_lo_ratio;



//==============================================
// Do simple version of the MET recoil tuning
//==============================================

// default use histogram
bool _doRecoil = false;

// use smoothed hist
bool _doRecoilUseSmooth = true;

// use graph from smoothed hist
bool _doRecoilUseSmoothGraph = true;

// input files
std::string _RecoilInputFileNameData_all, _RecoilInputFileNameData_mu, _RecoilInputFileNameData_el;
std::string _RecoilInputFileNameMC_all, _RecoilInputFileNameMC_mu, _RecoilInputFileNameMC_el;
std::string _RecoilInputFileNameMCLO_all, _RecoilInputFileNameMCLO_mu, _RecoilInputFileNameMCLO_el;

// not from config file
TFile* _file_dt_sigma[10];
TFile* _file_mc_sigma[10];
TH1D* _h_dt_met_para_shift[10];
TH1D* _h_mc_met_para_shift[10];
TH1D* _h_met_para_shift_dtmc[10];
TH1D* _h_dt_met_para_sigma[10];
TH1D* _h_dt_met_perp_sigma[10];
TH1D* _h_mc_met_para_sigma[10];
TH1D* _h_mc_met_perp_sigma[10];
TH1D* _h_ratio_met_para_sigma_dtmc[10];
TH1D* _h_ratio_met_perp_sigma_dtmc[10];
TGraphErrors* _gr_dt_met_para_shift[10];
TGraphErrors* _gr_mc_met_para_shift[10];
TGraphErrors* _gr_met_para_shift_dtmc[10];
TGraphErrors* _gr_ratio_met_para_sigma_dtmc[10];
TGraphErrors* _gr_ratio_met_perp_sigma_dtmc[10];


//==============================================
// eff scale
//==============================================
bool _addEffScale = true;
bool _addEffScaleOnData = false;
// Input files for:
// - el id iso eff
std::string _EffScaleInputFileName_IdIso_El = "data/eff/egammaEffi.txt_SF2D.root";
// - el tracking eff
std::string _EffScaleInputFileName_Trk_El = "data/eff/egammatracking.root";
// - mu id iso eff
std::string _EffScaleInputFileName_IdIso_Mu = "data/eff/muon80x12p9.root";
// - mu tracking eff
std::string _EffScaleInputFileName_Trk_Mu = "data/eff/muontrackingsf.root";
// - el trigger eff
std::string _EffScaleInputFileName_Trg_El = "data/eff/trigereff12p9.root";
// - mu trigger eff
std::string _EffScaleInputFileName_Trg_Mu = "data/eff/trigeff_mu.root";


// not from config file
// electron sf
TFile* _file_idiso_el;
TH2F* _h_sf_idiso_el;
// electron tracking sf
TFile* _file_trksf_el;
TH2F* _h_sf_trk_el;

// muon tracking sf
TFile* _file_trksf_mu;
TH1F* _h_sf_trk_mu;

// muon id iso sf
TFile* _file_idiso_mu;
TH2F* _h_eff_trkhpt_mu_dt;
TH2F* _h_eff_trkhpt_mu_mc;
TH2F* _h_eff_hpt_mu_dt;
TH2F* _h_eff_hpt_mu_mc;
TH2F* _h_sf_iso_mu;

// electron trigger sf
TFile* _file_trg_el;
TH2D* _h_sf_trg_el_l1;


// muon trigger sf
TFile* _file_trg_mu;
TH2D* _h_eff_trg_mu_l1_tot;
TH2D* _h_eff_trg_mu_l2_tot;
TH2D* _h_eff_trg_mu_l1_l1p;
TH2D* _h_eff_trg_mu_l2_l1p;
TH2D* _h_eff_trg_mu_l1_l1f;
TH2D* _h_eff_trg_mu_l2_l1f;
TH2D* _h_eff_trg_mu_l1_l1pl2f;
TH2D* _h_eff_trg_mu_l1_l1pl2p;
TH2D* _h_eff_trg_mu_l1_l1fl2p;
TH2D* _h_eff_trg_mu_l2_l1pl2f;
TH2D* _h_eff_trg_mu_l2_l1pl2p;
TH2D* _h_eff_trg_mu_l2_l1fl2p;
Double_t _N_eff_trg_mu_tot;
Double_t _N_eff_trg_mu_tot_err;
Double_t _N_eff_trg_mu_l1pl2f;
Double_t _N_eff_trg_mu_l1pl2f_err;
Double_t _N_eff_trg_mu_l1pl2p;
Double_t _N_eff_trg_mu_l1pl2p_err;
Double_t _N_eff_trg_mu_l1fl2p;
Double_t _N_eff_trg_mu_l1fl2p_err;
Double_t _N_eff_trg_mu_l1p;
Double_t _N_eff_trg_mu_l1p_err;
Double_t _N_eff_trg_mu_l1f;
Double_t _N_eff_trg_mu_l1f_err;
Int_t _NPtBins_eff_trg_mu;
Int_t _NEtaBins_eff_trg_mu;
TH2D* _h_eff_trg_mu_l1_tot_norm;
TH2D* _h_eff_trg_mu_l2_tot_norm;
TH2D* _h_eff_trg_mu_l1_l1p_norm;
TH2D* _h_eff_trg_mu_l1_l1f_norm;
TH2D* _h_eff_trg_mu_l2_l1p_norm;
TH2D* _h_eff_trg_mu_l2_l1f_norm;
TH2D* _h_eff_trg_mu_l1_l1pl2f_norm;
TH2D* _h_eff_trg_mu_l1_l1pl2p_norm;
TH2D* _h_eff_trg_mu_l1_l1fl2p_norm;
TH2D* _h_eff_trg_mu_l2_l1pl2f_norm;
TH2D* _h_eff_trg_mu_l2_l1pl2p_norm;
TH2D* _h_eff_trg_mu_l2_l1fl2p_norm;
TH2D* _h_eff_trg_mu_l1_l1p_norm_vs_tot;
TH2D* _h_eff_trg_mu_l1_l1f_norm_vs_tot;
TH2D* _h_eff_trg_mu_l2_l1p_norm_vs_tot;
TH2D* _h_eff_trg_mu_l2_l1f_norm_vs_tot;
TH2D* _h_eff_trg_mu_l1_l1pl2f_norm_vs_tot;
TH2D* _h_eff_trg_mu_l1_l1pl2p_norm_vs_tot;
TH2D* _h_eff_trg_mu_l1_l1fl2p_norm_vs_tot;
TH2D* _h_eff_trg_mu_l2_l1pl2f_norm_vs_tot;
TH2D* _h_eff_trg_mu_l2_l1pl2p_norm_vs_tot;
TH2D* _h_eff_trg_mu_l2_l1fl2p_norm_vs_tot;
TH2D* _h_eff_trg_mu_l1_l1pl2f_norm_vs_l1p;
TH2D* _h_eff_trg_mu_l1_l1pl2p_norm_vs_l1p;
TH2D* _h_eff_trg_mu_l1_l1fl2p_norm_vs_l1f;
TH2D* _h_eff_trg_mu_l2_l1pl2f_norm_vs_l1p;
TH2D* _h_eff_trg_mu_l2_l1pl2p_norm_vs_l1p;
TH2D* _h_eff_trg_mu_l2_l1fl2p_norm_vs_l1f;






//======================================================
// ╔╦╗╦═╗╔═╗╔═╗  ╦  ╦╔═╗╦═╗╦╔═╗╔╗ ╦  ╔═╗╔═╗
//  ║ ╠╦╝║╣ ║╣   ╚╗╔╝╠═╣╠╦╝║╠═╣╠╩╗║  ║╣ ╚═╗
//  ╩ ╩╚═╚═╝╚═╝   ╚╝ ╩ ╩╩╚═╩╩ ╩╚═╝╩═╝╚═╝╚═╝
//======================================================


// sum events, sum weights
Double_t _SumEvents, _SumWeights;

// isData
Int_t _isData;

// run, lumi, evt
UInt_t _run, _lumi;
ULong64_t _evt;

// leptons, Z, MET
Float_t _llnunu_mt;
Float_t _llnunu_l1_mass, _llnunu_l1_mt;
Float_t _llnunu_l1_pt, _llnunu_l1_phi, _llnunu_l1_eta;
Float_t _llnunu_l1_deltaPhi, _llnunu_l1_deltaR, _llnunu_l1_rapidity;
Float_t _llnunu_l2_pt, _llnunu_l2_phi;
Float_t _llnunu_l1_l1_pt, _llnunu_l1_l1_eta, _llnunu_l1_l1_phi;
Float_t _llnunu_l1_l1_rapidity, _llnunu_l1_l1_mass, _llnunu_l1_l1_ptErr;
Int_t   _llnunu_l1_l1_pdgId, _llnunu_l1_l1_charge;
Float_t _llnunu_l1_l2_pt, _llnunu_l1_l2_eta, _llnunu_l1_l2_phi;
Float_t _llnunu_l1_l2_rapidity, _llnunu_l1_l2_mass, _llnunu_l1_l2_ptErr;
Int_t   _llnunu_l1_l2_pdgId, _llnunu_l1_l2_charge;
Float_t _llnunu_l1_l1_eSCeta, _llnunu_l1_l2_eSCeta;

// MC Only
Int_t   _nTrueInt;
Float_t _genWeight;
Int_t   _ngenZ;
Float_t _genZ_pt[10];


// tree_out only

// store old branches
Float_t _llnunu_mt_old;
Float_t _llnunu_l1_mass_old, _llnunu_l1_pt_old, _llnunu_l1_phi_old, _llnunu_l1_eta_old;
Float_t _llnunu_l2_pt_old, _llnunu_l2_phi_old;
Float_t _llnunu_l1_l1_pt_old, _llnunu_l1_l1_eta_old, _llnunu_l1_l1_phi_old, _llnunu_l1_l1_ptErr_old;
Float_t _llnunu_l1_l2_pt_old, _llnunu_l1_l2_eta_old, _llnunu_l1_l2_phi_old, _llnunu_l1_l2_ptErr_old;

// zpt gen weights
Float_t _ZPtWeight, _ZPtWeight_up, _ZPtWeight_dn;
Float_t _ZJetsGenWeight;

// efficiency scale factors
Float_t _trgsf, _isosf, _idsf, _trksf, _idisotrksf;
Float_t _trgsf_err, _isosf_err, _idsf_err, _trksf_err;
Float_t _trgsf_up, _trgsf_dn, _idisotrksf_up, _idisotrksf_dn;



//======================================================
//  ╔╦╗╔═╗╔═╗╦╔╗╔╔═╗  ╔═╗╦ ╦╔╗╔╔═╗╔╦╗╦╔═╗╔╗╔╔═╗
//   ║║║╣ ╠╣ ║║║║║╣   ╠╣ ║ ║║║║║   ║ ║║ ║║║║╚═╗
//  ═╩╝╚═╝╚  ╩╝╚╝╚═╝  ╚  ╚═╝╝╚╝╚═╝ ╩ ╩╚═╝╝╚╝╚═╝
//======================================================

// read configration file
void readConfigFile();

// prepare the trees, if no entries in _tree_in, stop the program.
bool prepareTrees();

// store Old Branches
void storeOldBranches();

// prepare inputs for pu weights
void preparePUWeights();

// add more pileup weights
void addPUWeights();

// prepare inputs for muon re-calib
void prepareMuonPtRecalib();

// do muon re-calib
void doMuonPtRecalib();

// prepare inputs for addDyZPtWeight
void prepareDyZPtWeight();

// addDyZPtWeight
void addDyZPtWeight();

// prepare inputs for simple met recoil tune.
void prepareRecoil();

// do simple met recoil fine tuning
void doRecoil();

// prepare eff scale factors
void prepareEffScale();
    
// add eff scale factors
void addEffScale();

// 


//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================
//======================================================


