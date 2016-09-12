// Hengne Li, 2016 @ CERN, initial version.

#include "metcorr.h"


int main(int argc, char** argv) {

  if( argc<5 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: met kin-fit framework, actually a full analysis framework  ... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " config.file inputfile.root outputfile.root Nevts SumWeights " << std::endl ;
     exit(1) ;
  }

  // config file name
  _file_config_name = std::string((const char*)argv[1]);

  // input file name
  _file_in_name = std::string((const char*)argv[2]);

  // output file name
  _file_out_name = std::string((const char*)argv[3]);

  // sum events
  _SumEvents = atof(argv[4]);

  // sum weights
  _SumWeights = atof(argv[5]);


  // input file
  _file_in = TFile::Open(_file_in_name.c_str());

  // output file 
  _file_out = TFile::Open(_file_out_name.c_str(), "recreate");

  // special for DYJets MC samples,
  //  decide if this is the DYJets samples, and LO or NLO, 
  //  based on the output file names. So...note below:
  // ATTENTION: 
  //   Makes sure output NLO DYJets MC output file name has string
  //   "DYJets" but does not have "MGMLM", and LO DYJets one has 
  //   both strings "DYJets" and "MGMLM". 
  _isDyJets = (_file_out_name.find("DYJets")!=std::string::npos);
  _isDyJetsLO = (_file_out_name.find("DYJets")!=std::string::npos && _file_out_name.find("MGMLM")!=std::string::npos);
  
  if (_debug) std::cout << "DEBUG: isDyJets = " << _isDyJets << ", isDyJetsLO = " << _isDyJetsLO << std::endl;


  // read config file
  readConfigFile();

  // prepare the trees
  prepareTrees();

  // prepare inputs for pu weights
  if (_addPUWeights && !_isData) preparePUWeights();

  // prepare inputs for muon re-calib
  if (_doMuonPtRecalib) prepareMuonPtRecalib();


  // prepare inputs for addDyZPtWeight
  if (_addDyZPtWeight && !_isData && _isDyJets) prepareDyZPtWeight();

  // prepare inputs for simple met recoil tune.
  if (_doRecoil && !_isData && _isDyJets) prepareRecoil();


  // prepare eff scale factors
  if (_addEffScale && (!_isData || _addEffScaleOnData) ) prepareEffScale();

  // loop
  int n_pass = 0;

  // get n selected entries to loop 
  Long64_t n_selected_entries = _selected_entries->GetN();

  for (Int_t i_slist = _n_start; i_slist < n_selected_entries; i_slist++) {

    // get _tree_in's entry number from selected list
    Int_t i = _selected_entries->GetEntry(i_slist);

    // get the selected entry from tree_in
    _tree_in->GetEntry(i);

    // beak if runs more events than  you want to test
    n_pass++;
    if (_n_test>0 && n_pass>_n_test) break;

    //
    if (_debug || i%_n_interval == 0) {
      if (_debug) std::cout << "#############################################" << std::endl;
      std::cout << "##  Entry " << i << ", Run " << _run << ", Event " << _evt << std::endl;
      if (_debug) std::cout << "#############################################" << std::endl;
    }

    // _storeOldBranches
    if (_storeOldBranches) storeOldBranches();

    // add pu weights
    if (_addPUWeights && !_isData) addPUWeights();

    // do muon re-calib
    if (_doMuonPtRecalib) doMuonPtRecalib(); 

    //  addDyZPtWeight
    if (_addDyZPtWeight && !_isData && _isDyJets) addDyZPtWeight();

    // simple met recoil tune.
    if (_doRecoil && !_isData && _isDyJets) doRecoil();

    // add eff scale factors
    if (_addEffScale && (!_isData || _addEffScaleOnData) ) addEffScale();

    // fill output tree
    _tree_out->Fill(); 
  }  


  // store output tree
  _file_out->cd();
  _tree_out->Write();
  _file_out->Close();

  return 0;

}





//======================================================
// ╦═╗╔═╗╔═╗╔╦╗  ╔═╗╔═╗╔╗╔╔═╗╦╔═╗  ╔═╗╦╦  ╔═╗
// ╠╦╝║╣ ╠═╣ ║║  ║  ║ ║║║║╠╣ ║║ ╦  ╠╣ ║║  ║╣ 
// ╩╚═╚═╝╩ ╩═╩╝  ╚═╝╚═╝╝╚╝╚  ╩╚═╝  ╚  ╩╩═╝╚═╝
//======================================================

void readConfigFile() 
{

  // init paramter file
  PParameterReader parm(_file_config_name.c_str());

  // print it
  parm.Print();

  //===============================================
  // NOTE: 
  //  - See notes above in parameter definitions
  //     for the meanings and usage of them.
  //  - Please always use the same variable name
  //     in the c++ codes and in the config files,
  //     with the c++ one has "_" in front, 
  //      e.g.:   
  //       bool _AAA = parm.GetBool("AAA", kFALSE);
  //     This rule can minimize mistakes.
  //==============================================

  //======================
  // general parameters
  //======================
  _debug = parm.GetBool("debug", kFALSE);
  _n_start = parm.GetInt("n_start", 0);
  _n_test = parm.GetInt("n_test", -1);
  _n_interval = parm.GetInt("n_interval", 3000);
  _selection = parm.GetString("selection", "(1)");
  _storeOldBranches = parm.GetBool("storeOldBranches", kFALSE);


  //==========================
  // use slimmed tree or not
  //=========================
  _useLightTree = parm.GetBool("useLightTree", kTRUE);


  //=========================
  // add PU weights
  //=========================
  _addPUWeights = parm.GetBool("addPUWeights", kTRUE);
  _PUTags = parm.GetVString("PUTags");
  _PUInputDir = parm.GetString("PUInputDir", "data/pileup");
  _PUInputFileNames = parm.GetVString("PUInputFileNames");
  _PUWeightHistName = parm.GetString("PUWeightHistName", "puweight_dtmc");


  //==========================
  // muon pT recalibration
  //==========================
  _doMuonPtRecalib = parm.GetBool("doMuonPtRecalib", kFALSE);
  _MuonPtRecalibInputForData = parm.GetString("MuonPtRecalibInputForData", "data/kalman/DATA_80X_13TeV.root");
  _MuonPtRecalibInputForMC = parm.GetString("MuonPtRecalibInputForMC", "data/kalman/MC_80X_13TeV.root");


  //========================
  // Add DYJet gen reweight 
  //========================
  _addDyZPtWeight = parm.GetBool("addDyZPtWeight", kTRUE);
  _addDyZPtWeightUseFunction = parm.GetBool("addDyZPtWeightUseFunction", kTRUE);
  _addDyZPtWeightLOUseFunction = parm.GetBool("addDyZPtWeightLOUseFunction", kTRUE);
  _DyZPtWeightInputFileName = parm.GetString("DyZPtWeightInputFileName", "data/zptweight/dyjets_zpt_weight_lo_nlo_sel.root");
  _addDyNewGenWeight = parm.GetBool("addDyNewGenWeight", kTRUE);;

  //==============================================
  // Do simple version of the MET recoil tuning
  //==============================================
  _doRecoil = parm.GetBool("doRecoil", kFALSE);
  _doRecoilUseSmooth = parm.GetBool("doRecoilUseSmooth", kTRUE);
  _doRecoilUseSmoothGraph = parm.GetBool("doRecoilUseSmoothGraph", kTRUE);
  _RecoilInputFileNameData_all = parm.GetString("RecoilInputFileNameData_all", "data/recoil/SingleEMU_Run2016BCD_PromptReco_met_para_study.root"); 
  _RecoilInputFileNameData_mu = parm.GetString("RecoilInputFileNameData_mu", "data/recoil/SingleEMU_Run2016BCD_PromptReco_met_para_study_mu.root"); 
  _RecoilInputFileNameData_el = parm.GetString("RecoilInputFileNameData_el", "data/recoil/SingleEMU_Run2016BCD_PromptReco_met_para_study_el.root"); 
  _RecoilInputFileNameMC_all = parm.GetString("RecoilInputFileNameMC_all", "data/recoil/DYJetsToLL_M50_NoRecoil_met_para_study.root"); 
  _RecoilInputFileNameMC_mu = parm.GetString("RecoilInputFileNameMC_mu", "data/recoil/DYJetsToLL_M50_NoRecoil_met_para_study_mu.root"); 
  _RecoilInputFileNameMC_el = parm.GetString("RecoilInputFileNameMC_el", "data/recoil/DYJetsToLL_M50_NoRecoil_met_para_study_el.root"); 
  _RecoilInputFileNameMCLO_all = parm.GetString("RecoilInputFileNameMCLO_all", "data/recoil/DYJetsToLL_M50_MGMLM_Ext1_NoRecoil_met_para_study.root"); 
  _RecoilInputFileNameMCLO_mu = parm.GetString("RecoilInputFileNameMCLO_mu", "data/recoil/DYJetsToLL_M50_MGMLM_Ext1_NoRecoil_met_para_study_mu.root"); 
  _RecoilInputFileNameMCLO_el = parm.GetString("RecoilInputFileNameMCLO_el", "data/recoil/DYJetsToLL_M50_MGMLM_Ext1_NoRecoil_met_para_study_el.root"); 


  //==============================================
  // Add efficiency scale factors
  //==============================================  
  _addEffScale = parm.GetBool("addEffScale", kFALSE);
  _addEffScaleOnData = parm.GetBool("addEffScaleOnData", kFALSE);
  _EffScaleInputFileName_IdIso_El = parm.GetString("EffScaleInputFileName_IdIso_El", "data/eff/egammaEffi.txt_SF2D.root");
  _EffScaleInputFileName_Trk_El = parm.GetString("EffScaleInputFileName_Trk_El", "data/eff/egammatracking.root");
  _EffScaleInputFileName_IdIso_Mu = parm.GetString("EffScaleInputFileName_IdIso_Mu", "data/eff/muon80x12p9.root");
  _EffScaleInputFileName_Trk_Mu = parm.GetString("EffScaleInputFileName_Trk_Mu", "data/eff/muontrackingsf.root");
  _EffScaleInputFileName_Trg_El = parm.GetString("EffScaleInputFileName_Trg_El", "data/eff/trigereff12p9.root");
  _EffScaleInputFileName_Trg_Mu = parm.GetString("EffScaleInputFileName_Trg_Mu", "data/eff/trigeff_mu.root");

}




// prepare the trees
bool  prepareTrees() 
{

  // 1.) prepare in put tree:

  // input tree
  _tree_in = (TTree*)_file_in->Get("tree");

  // selection 
  _tree_in->Draw(">>_selected_entries", _selection.c_str(), "entrylist");
  _selected_entries = (TEntryList*)gDirectory->Get("_selected_entries");

  // isData  
  _tree_in->SetBranchAddress("isData",&_isData);

  // check if tree has events
  if (_tree_in->GetEntries()<=0) {
    std::cout << "prepareTrees():: input tree has no event pass selection. Quit." <<std::endl;
    return false;
  }

  // get isData info
  _tree_in->GetEntry(0);

  // set common branches
  _tree_in->SetBranchAddress("run", &_run);
  _tree_in->SetBranchAddress("lumi", &_lumi);
  _tree_in->SetBranchAddress("evt", &_evt);

  _tree_in->SetBranchAddress("llnunu_mt", &_llnunu_mt);
  _tree_in->SetBranchAddress("llnunu_l1_mass",&_llnunu_l1_mass);
  _tree_in->SetBranchAddress("llnunu_l1_mt", &_llnunu_l1_mt);
  _tree_in->SetBranchAddress("llnunu_l1_pt", &_llnunu_l1_pt);
  _tree_in->SetBranchAddress("llnunu_l1_phi", &_llnunu_l1_phi);
  _tree_in->SetBranchAddress("llnunu_l1_eta", &_llnunu_l1_eta);
  _tree_in->SetBranchAddress("llnunu_l1_deltaPhi", &_llnunu_l1_deltaPhi);
  _tree_in->SetBranchAddress("llnunu_l1_deltaR", &_llnunu_l1_deltaR);
  _tree_in->SetBranchAddress("llnunu_l1_rapidity", &_llnunu_l1_rapidity);

  _tree_in->SetBranchAddress("llnunu_l2_pt", &_llnunu_l2_pt);
  _tree_in->SetBranchAddress("llnunu_l2_phi", &_llnunu_l2_phi);

  _tree_in->SetBranchAddress("llnunu_l1_l1_pt", &_llnunu_l1_l1_pt);
  _tree_in->SetBranchAddress("llnunu_l1_l1_eta", &_llnunu_l1_l1_eta);
  _tree_in->SetBranchAddress("llnunu_l1_l1_phi", &_llnunu_l1_l1_phi);
  _tree_in->SetBranchAddress("llnunu_l1_l1_rapidity", &_llnunu_l1_l1_rapidity);
  _tree_in->SetBranchAddress("llnunu_l1_l1_mass", &_llnunu_l1_l1_mass);
  _tree_in->SetBranchAddress("llnunu_l1_l1_pdgId", &_llnunu_l1_l1_pdgId);
  _tree_in->SetBranchAddress("llnunu_l1_l1_charge", &_llnunu_l1_l1_charge);
  _tree_in->SetBranchAddress("llnunu_l1_l1_ptErr", &_llnunu_l1_l1_ptErr);
  _tree_in->SetBranchAddress("llnunu_l1_l1_eSCeta", &_llnunu_l1_l1_eSCeta);

  _tree_in->SetBranchAddress("llnunu_l1_l2_pt", &_llnunu_l1_l2_pt);
  _tree_in->SetBranchAddress("llnunu_l1_l2_eta", &_llnunu_l1_l2_eta);
  _tree_in->SetBranchAddress("llnunu_l1_l2_phi", &_llnunu_l1_l2_phi);
  _tree_in->SetBranchAddress("llnunu_l1_l2_rapidity", &_llnunu_l1_l2_rapidity);
  _tree_in->SetBranchAddress("llnunu_l1_l2_mass", &_llnunu_l1_l2_mass);
  _tree_in->SetBranchAddress("llnunu_l1_l2_pdgId", &_llnunu_l1_l2_pdgId);
  _tree_in->SetBranchAddress("llnunu_l1_l2_charge", &_llnunu_l1_l2_charge);
  _tree_in->SetBranchAddress("llnunu_l1_l2_ptErr", &_llnunu_l1_l2_ptErr);
  _tree_in->SetBranchAddress("llnunu_l1_l2_eSCeta", &_llnunu_l1_l2_eSCeta);


  // other branches for not light weight tree   
  if (!_useLightTree) {
    std::cout << " Warning: for not _useLightTree, to be implemented." << std::endl;
  }

  // 
  // MC only
  if (!_isData) {
    _tree_in->SetBranchAddress("nTrueInt", &_nTrueInt );
    _tree_in->SetBranchAddress("genWeight",&_genWeight);
    _tree_in->SetBranchAddress("ngenZ", &_ngenZ);
    _tree_in->SetBranchAddress("genZ_pt", _genZ_pt);
  }


  // 2.) Output tree

  // output tree
  _file_out->cd();
  _tree_out = _tree_in->CloneTree(0);

  // sum of events and weights
  _tree_out->Branch("SumEvents", &_SumEvents, "SumEvents/D");
  _tree_out->Branch("SumWeights", &_SumWeights, "SumWeights/D");

  // keep a copy of old branches
  if (_storeOldBranches) {
    _tree_out->Branch("llnunu_mt_old", &_llnunu_mt_old, "llnunu_mt_old/F");
    _tree_out->Branch("llnunu_l1_mass_old", &_llnunu_l1_mass_old, "llnunu_l1_mass_old/F");
    _tree_out->Branch("llnunu_l1_pt_old", &_llnunu_l1_pt_old, "llnunu_l1_pt_old/F");
    _tree_out->Branch("llnunu_l1_phi_old", &_llnunu_l1_phi_old, "llnunu_l1_phi_old/F");
    _tree_out->Branch("llnunu_l1_eta_old", &_llnunu_l1_eta_old, "llnunu_l1_eta_old/F");
    _tree_out->Branch("llnunu_l2_pt_old", &_llnunu_l2_pt_old, "llnunu_l2_pt_old/F");
    _tree_out->Branch("llnunu_l2_phi_old", &_llnunu_l2_phi_old, "llnunu_l2_phi_old/F");
    _tree_out->Branch("llnunu_l1_l1_pt_old", &_llnunu_l1_l1_pt_old, "llnunu_l1_l1_pt_old/F");
    _tree_out->Branch("llnunu_l1_l1_eta_old", &_llnunu_l1_l1_eta_old, "llnunu_l1_l1_eta_old/F");
    _tree_out->Branch("llnunu_l1_l1_phi_old", &_llnunu_l1_l1_phi_old, "llnunu_l1_l1_phi_old/F");
    _tree_out->Branch("llnunu_l1_l1_ptErr_old", &_llnunu_l1_l1_ptErr_old, "llnunu_l1_l1_ptErr_old/F");
    _tree_out->Branch("llnunu_l1_l2_pt_old", &_llnunu_l1_l2_pt_old, "llnunu_l1_l2_pt_old/F");
    _tree_out->Branch("llnunu_l1_l2_eta_old", &_llnunu_l1_l2_eta_old, "llnunu_l1_l2_eta_old/F");
    _tree_out->Branch("llnunu_l1_l2_phi_old", &_llnunu_l1_l2_phi_old, "llnunu_l1_l2_phi_old/F");
    _tree_out->Branch("llnunu_l1_l2_ptErr_old", &_llnunu_l1_l2_ptErr_old, "llnunu_l1_l2_ptErr_old/F");
  }


  return true;
}

// store Old Branches
void storeOldBranches()
{
  _llnunu_mt_old = _llnunu_mt;
  _llnunu_l1_mass_old = _llnunu_l1_mass;
  _llnunu_l1_pt_old = _llnunu_l1_pt;
  _llnunu_l1_phi_old = _llnunu_l1_phi;
  _llnunu_l1_eta_old = _llnunu_l1_eta;
  _llnunu_l2_pt_old = _llnunu_l2_pt;
  _llnunu_l2_phi_old = _llnunu_l2_phi;
  _llnunu_l1_l1_pt_old = _llnunu_l1_l1_pt;
  _llnunu_l1_l1_eta_old = _llnunu_l1_l1_eta;
  _llnunu_l1_l1_phi_old = _llnunu_l1_l1_phi;
  _llnunu_l1_l1_ptErr_old = _llnunu_l1_l1_ptErr;
  _llnunu_l1_l2_pt_old = _llnunu_l1_l2_pt;
  _llnunu_l1_l2_eta_old = _llnunu_l1_l2_eta;
  _llnunu_l1_l2_phi_old = _llnunu_l1_l2_phi;
  _llnunu_l1_l2_ptErr_old = _llnunu_l1_l2_ptErr;
}

// prepare inputs for pu weights
void preparePUWeights()
{
  // for each puWeight, read input weight hist and create branch
  for (int i=0; i<(int)_PUTags.size(); i++) {
    sprintf(name, "%s/%s", _PUInputDir.c_str(), _PUInputFileNames.at(i).c_str());
    _PUFiles.push_back(new TFile(name));
    _PUHists.push_back((TH1D*)_PUFiles.at(i)->Get(_PUWeightHistName.c_str()));
    _PUWeights.push_back(new Float_t(1.0));
    sprintf(name, "puWeight%s", _PUTags.at(i).c_str());
    sprintf(name1, "puWeight%s/F", _PUTags.at(i).c_str());
    _tree_out->Branch(name, _PUWeights.at(i),name1);
  }

}

// add more pileup weights
void addPUWeights() 
{
  // for each puWeight, get the weight
  for (int i=0; i<(int)_PUTags.size(); i++){
    *(_PUWeights.at(i)) = _PUHists.at(i)->GetBinContent(_PUHists.at(i)->FindBin(_nTrueInt));
  }
}


// prepare inputs for muon re-calib
void prepareMuonPtRecalib()
{
  if (_isData) _muCalib = new KalmanMuonCalibrator(_MuonPtRecalibInputForData.c_str());
  else _muCalib = new KalmanMuonCalibrator(_MuonPtRecalibInputForMC.c_str());
}


// do muon re-calib
void doMuonPtRecalib()
{
  if (abs(_llnunu_l1_l1_pdgId)==13&&abs(_llnunu_l1_l2_pdgId)==13 && _llnunu_l1_l1_pt>2. && _llnunu_l1_l1_pt<200. ) {
    _llnunu_l1_l1_pt = (Float_t)_muCalib->getCorrectedPt(double(_llnunu_l1_l1_pt), double(_llnunu_l1_l1_eta), double(_llnunu_l1_l1_phi), double(_llnunu_l1_l1_charge));
    _llnunu_l1_l2_pt = (Float_t)_muCalib->getCorrectedPt(double(_llnunu_l1_l2_pt), double(_llnunu_l1_l2_eta), double(_llnunu_l1_l2_phi), double(_llnunu_l1_l2_charge));
    if (!_isData) {
      _llnunu_l1_l1_pt = (Float_t)_muCalib->smear(double(_llnunu_l1_l1_pt), double(_llnunu_l1_l1_eta));
      _llnunu_l1_l2_pt = (Float_t)_muCalib->smear(double(_llnunu_l1_l2_pt), double(_llnunu_l1_l2_eta));
    }
    _llnunu_l1_l1_ptErr = _llnunu_l1_l1_pt*(Float_t)_muCalib->getCorrectedError(double(_llnunu_l1_l1_pt), double(_llnunu_l1_l1_eta), double(_llnunu_l1_l1_ptErr/_llnunu_l1_l1_pt));
    _llnunu_l1_l2_ptErr = _llnunu_l1_l2_pt*(Float_t)_muCalib->getCorrectedError(double(_llnunu_l1_l2_pt), double(_llnunu_l1_l2_eta), double(_llnunu_l1_l2_ptErr/_llnunu_l1_l2_pt));
    TLorentzVector l1v, l2v;
    l1v.SetPtEtaPhiM(_llnunu_l1_l1_pt, _llnunu_l1_l1_eta, _llnunu_l1_l1_phi, _llnunu_l1_l1_mass);
    l2v.SetPtEtaPhiM(_llnunu_l1_l2_pt, _llnunu_l1_l2_eta, _llnunu_l1_l2_phi, _llnunu_l1_l2_mass);
    TLorentzVector zv = l1v+l2v;
    _llnunu_l1_l1_rapidity = (Float_t)l1v.Rapidity(); 
    _llnunu_l1_l2_rapidity = (Float_t)l2v.Rapidity();
    _llnunu_l1_pt = (Float_t)zv.Pt();
    _llnunu_l1_eta = (Float_t)zv.Eta();
    _llnunu_l1_phi = (Float_t)zv.Phi();
    _llnunu_l1_rapidity = (Float_t)zv.Rapidity();
    _llnunu_l1_deltaPhi = (Float_t)l1v.DeltaPhi(l2v);
    _llnunu_l1_deltaR = (Float_t)l1v.DeltaR(l2v);
    _llnunu_l1_mt = (Float_t)zv.Mt();
    _llnunu_l1_mass = (Float_t)zv.M();

    TVector2 vec_met(_llnunu_l2_pt*cos(_llnunu_l2_phi), _llnunu_l2_pt*sin(_llnunu_l2_phi));

    Float_t et1 = TMath::Sqrt(_llnunu_l1_mass*_llnunu_l1_mass + _llnunu_l1_pt*_llnunu_l1_pt);
    Float_t et2 = TMath::Sqrt(_llnunu_l1_mass*_llnunu_l1_mass + _llnunu_l2_pt*_llnunu_l2_pt);
    _llnunu_mt = TMath::Sqrt(2.0*_llnunu_l1_mass*_llnunu_l1_mass+2.0*(et1*et2
               -_llnunu_l1_pt*cos(_llnunu_l1_phi)*_llnunu_l2_pt*cos(_llnunu_l2_phi)
               -_llnunu_l1_pt*sin(_llnunu_l1_phi)*_llnunu_l2_pt*sin(_llnunu_l2_phi)));
  }
}


// prepare inputs for addDyZPtWeight
void prepareDyZPtWeight()
{

  // needed output branches
  _tree_out->Branch("ZPtWeight", &_ZPtWeight, "ZPtWeight/F");
  _tree_out->Branch("ZPtWeight_up", &_ZPtWeight_up, "ZPtWeight_up/F");
  _tree_out->Branch("ZPtWeight_dn", &_ZPtWeight_dn, "ZPtWeight_dn/F");
  if (_addDyNewGenWeight) {
    _tree_out->Branch("ZJetsGenWeight", &_ZJetsGenWeight, "ZJetsGenWeight/F");
  }

  // get input file and building materials
  _fdyzpt = new TFile(_DyZPtWeightInputFileName.c_str());
  _hdyzpt_dtmc_ratio = (TH1D*)_fdyzpt->Get("hdyzpt_dtmc_ratio");
  _fcdyzpt_dtmc_ratio = (TF1*)_fdyzpt->Get("fcdyzpt_dtmc_ratio");
  _hdyzpt_mc_nlo_lo_ratio = (TH1D*)_fdyzpt->Get("hdyzpt_mc_nlo_lo_ratio");
  _fcdyzpt_mc_nlo_lo_ratio = (TF1*)_fdyzpt->Get("fcdyzpt_mc_nlo_lo_ratio");
  
}

// addDyZPtWeight
void addDyZPtWeight()
{
  if (_addDyZPtWeight && !_isData && _isDyJets) {
    Int_t zptBin=0;
    if (_ngenZ>0) {
      zptBin = _hdyzpt_dtmc_ratio->FindBin(_genZ_pt[0]);
      if (_genZ_pt[0]>1000) zptBin = _hdyzpt_dtmc_ratio->FindBin(999);
    }
    else {
      zptBin = _hdyzpt_dtmc_ratio->FindBin(_llnunu_l1_pt);
      if (_llnunu_l1_pt>1000) zptBin = _hdyzpt_dtmc_ratio->FindBin(999);
    }
    if (_addDyZPtWeightUseFunction) {
      if (_ngenZ>0) _ZPtWeight = _fcdyzpt_dtmc_ratio->Eval(_genZ_pt[0]);
      else _ZPtWeight = _fcdyzpt_dtmc_ratio->Eval(_llnunu_l1_pt);
    }
    else {
      _ZPtWeight = _hdyzpt_dtmc_ratio->GetBinContent(zptBin);
    }
    _ZPtWeight_up = _ZPtWeight+0.5*_hdyzpt_dtmc_ratio->GetBinError(zptBin);
    _ZPtWeight_dn = _ZPtWeight-0.5*_hdyzpt_dtmc_ratio->GetBinError(zptBin);

    if (_isDyJetsLO) {
      // for LO ZJets samples
      zptBin=0;
      if (_ngenZ>0) {
        zptBin = _hdyzpt_mc_nlo_lo_ratio->FindBin(_genZ_pt[0]);
        if (_genZ_pt[0]>1000) zptBin = _hdyzpt_mc_nlo_lo_ratio->FindBin(999);
      }
      else {
        zptBin = _hdyzpt_mc_nlo_lo_ratio->FindBin(_llnunu_l1_pt);
        if (_llnunu_l1_pt>1000) zptBin = _hdyzpt_mc_nlo_lo_ratio->FindBin(999);
      }
      if (_addDyZPtWeightLOUseFunction) {
        if (_ngenZ>0) _ZPtWeight *= _fcdyzpt_mc_nlo_lo_ratio->Eval(_genZ_pt[0]);
        else _ZPtWeight *= _fcdyzpt_mc_nlo_lo_ratio->Eval(_llnunu_l1_pt);
      }
      else {
        _ZPtWeight *= _hdyzpt_mc_nlo_lo_ratio->GetBinContent(zptBin);
      }

      _ZPtWeight_up = _ZPtWeight+0.5*_hdyzpt_dtmc_ratio->GetBinError(zptBin);
      _ZPtWeight_dn = _ZPtWeight-0.5*_hdyzpt_dtmc_ratio->GetBinError(zptBin);

    }
  }

  if (_addDyNewGenWeight && !_isData && _isDyJets) {
    Float_t ZJetsLOSumWeights(49877138);
    Float_t ZJetsNLOSumWeights(450670522117);
    Float_t ZJetsLOSumEvents(49877138);
    Float_t ZJetsNLOSumEvents(28696958);
    if (!_isDyJetsLO) {
      _ZJetsGenWeight = _genWeight/ZJetsNLOSumWeights*ZJetsNLOSumEvents/(ZJetsNLOSumEvents+ZJetsLOSumEvents);
    }
    else {
      _ZJetsGenWeight = _genWeight/ZJetsLOSumWeights*ZJetsLOSumEvents/(ZJetsNLOSumEvents+ZJetsLOSumEvents);
    }
  }

}

// prepare inputs for simple met recoil tune.
void prepareRecoil()
{
  if (_doRecoil && _isDyJets) {
    // met shift  sigma
    _file_dt_sigma[0] = new TFile(_RecoilInputFileNameData_all.c_str());
    _file_dt_sigma[1] = new TFile(_RecoilInputFileNameData_mu.c_str());
    _file_dt_sigma[2] = new TFile(_RecoilInputFileNameData_el.c_str());
    if (!_isDyJetsLO) {
      _file_mc_sigma[0] = new TFile(_RecoilInputFileNameMC_all.c_str());
      _file_mc_sigma[1] = new TFile(_RecoilInputFileNameMC_mu.c_str());
      _file_mc_sigma[2] = new TFile(_RecoilInputFileNameMC_el.c_str());
    }
    else {
      _file_mc_sigma[0] = new TFile(_RecoilInputFileNameMCLO_all.c_str());
      _file_mc_sigma[1] = new TFile(_RecoilInputFileNameMCLO_mu.c_str());
      _file_mc_sigma[2] = new TFile(_RecoilInputFileNameMCLO_el.c_str());
    }
    _h_dt_met_para_shift[0] = (TH1D*)_file_dt_sigma[0]->Get("h_met_para_vs_zpt_mean");
    _h_mc_met_para_shift[0] = (TH1D*)_file_mc_sigma[0]->Get("h_met_para_vs_zpt_mean");
    _h_met_para_shift_dtmc[0] = (TH1D*)_h_dt_met_para_shift[0]->Clone("h_met_para_shift_dtmc_all");
    _h_met_para_shift_dtmc[0]->Add(_h_mc_met_para_shift[0], -1);

    _h_dt_met_para_sigma[0] = (TH1D*)_file_dt_sigma[0]->Get("h_met_para_vs_zpt_sigma");
    _h_dt_met_perp_sigma[0] = (TH1D*)_file_dt_sigma[0]->Get("h_met_perp_vs_zpt_sigma");
    _h_mc_met_para_sigma[0] = (TH1D*)_file_mc_sigma[0]->Get("h_met_para_vs_zpt_sigma");
    _h_mc_met_perp_sigma[0] = (TH1D*)_file_mc_sigma[0]->Get("h_met_perp_vs_zpt_sigma");

    _h_ratio_met_para_sigma_dtmc[0] = (TH1D*)_h_dt_met_para_sigma[0]->Clone("h_ratio_met_para_sigma_dtmc_all");
    _h_ratio_met_perp_sigma_dtmc[0] = (TH1D*)_h_dt_met_perp_sigma[0]->Clone("h_ratio_met_perp_sigma_dtmc_all");
    _h_ratio_met_para_sigma_dtmc[0]->Divide(_h_mc_met_para_sigma[0]);
    _h_ratio_met_perp_sigma_dtmc[0]->Divide(_h_mc_met_perp_sigma[0]);

    _h_dt_met_para_shift[1] = (TH1D*)_file_dt_sigma[1]->Get("h_met_para_vs_zpt_mean");
    _h_mc_met_para_shift[1] = (TH1D*)_file_mc_sigma[1]->Get("h_met_para_vs_zpt_mean");
    _h_met_para_shift_dtmc[1] = (TH1D*)_h_dt_met_para_shift[1]->Clone("h_met_para_shift_dtmc_mu");
    _h_met_para_shift_dtmc[1]->Add(_h_mc_met_para_shift[1], -1);

    _h_dt_met_para_sigma[1] = (TH1D*)_file_dt_sigma[1]->Get("h_met_para_vs_zpt_sigma");
    _h_dt_met_perp_sigma[1] = (TH1D*)_file_dt_sigma[1]->Get("h_met_perp_vs_zpt_sigma");
    _h_mc_met_para_sigma[1] = (TH1D*)_file_mc_sigma[1]->Get("h_met_para_vs_zpt_sigma");
    _h_mc_met_perp_sigma[1] = (TH1D*)_file_mc_sigma[1]->Get("h_met_perp_vs_zpt_sigma");

    _h_ratio_met_para_sigma_dtmc[1] = (TH1D*)_h_dt_met_para_sigma[1]->Clone("h_ratio_met_para_sigma_dtmc_mu");
    _h_ratio_met_perp_sigma_dtmc[1] = (TH1D*)_h_dt_met_perp_sigma[1]->Clone("h_ratio_met_perp_sigma_dtmc_mu");
    _h_ratio_met_para_sigma_dtmc[1]->Divide(_h_mc_met_para_sigma[1]);
    _h_ratio_met_perp_sigma_dtmc[1]->Divide(_h_mc_met_perp_sigma[1]);

    _h_dt_met_para_shift[2] = (TH1D*)_file_dt_sigma[2]->Get("h_met_para_vs_zpt_mean");
    _h_mc_met_para_shift[2] = (TH1D*)_file_mc_sigma[2]->Get("h_met_para_vs_zpt_mean");
    _h_met_para_shift_dtmc[2] = (TH1D*)_h_dt_met_para_shift[2]->Clone("h_met_para_shift_dtmc_el");
    _h_met_para_shift_dtmc[2]->Add(_h_mc_met_para_shift[2], -1);

    _h_dt_met_para_sigma[2] = (TH1D*)_file_dt_sigma[2]->Get("h_met_para_vs_zpt_sigma");
    _h_dt_met_perp_sigma[2] = (TH1D*)_file_dt_sigma[2]->Get("h_met_perp_vs_zpt_sigma");
    _h_mc_met_para_sigma[2] = (TH1D*)_file_mc_sigma[2]->Get("h_met_para_vs_zpt_sigma");
    _h_mc_met_perp_sigma[2] = (TH1D*)_file_mc_sigma[2]->Get("h_met_perp_vs_zpt_sigma");

    _h_ratio_met_para_sigma_dtmc[2] = (TH1D*)_h_dt_met_para_sigma[2]->Clone("h_ratio_met_para_sigma_dtmc_el");
    _h_ratio_met_perp_sigma_dtmc[2] = (TH1D*)_h_dt_met_perp_sigma[2]->Clone("h_ratio_met_perp_sigma_dtmc_el");
    _h_ratio_met_para_sigma_dtmc[2]->Divide(_h_mc_met_para_sigma[2]);
    _h_ratio_met_perp_sigma_dtmc[2]->Divide(_h_mc_met_perp_sigma[2]);

    // smooth functions
    _h_dt_met_para_shift[0]->SetName("h_dt_met_para_shift_all");
    _h_mc_met_para_shift[0]->SetName("h_mc_met_para_shift_all");
    _h_dt_met_para_shift[1]->SetName("h_dt_met_para_shift_mu");
    _h_mc_met_para_shift[1]->SetName("h_mc_met_para_shift_mu");
    _h_dt_met_para_shift[2]->SetName("h_dt_met_para_shift_el");
    _h_mc_met_para_shift[2]->SetName("h_mc_met_para_shift_el");
    _h_dt_met_para_shift[3] = (TH1D*)_h_dt_met_para_shift[0]->Clone("h_dt_met_para_shift_all_smooth");
    _h_mc_met_para_shift[3] = (TH1D*)_h_mc_met_para_shift[0]->Clone("h_mc_met_para_shift_all_smooth");
    _h_dt_met_para_shift[4] = (TH1D*)_h_dt_met_para_shift[1]->Clone("h_dt_met_para_shift_mu_smooth");
    _h_mc_met_para_shift[4] = (TH1D*)_h_mc_met_para_shift[1]->Clone("h_mc_met_para_shift_mu_smooth");
    _h_dt_met_para_shift[5] = (TH1D*)_h_dt_met_para_shift[2]->Clone("h_dt_met_para_shift_el_smooth");
    _h_mc_met_para_shift[5] = (TH1D*)_h_mc_met_para_shift[2]->Clone("h_mc_met_para_shift_el_smooth");
    _h_dt_met_para_shift[3]->Smooth();
    _h_mc_met_para_shift[3]->Smooth();
    _h_dt_met_para_shift[4]->Smooth();
    _h_mc_met_para_shift[4]->Smooth();
    _h_dt_met_para_shift[5]->Smooth();
    _h_mc_met_para_shift[5]->Smooth();

    _gr_dt_met_para_shift[3] = new TGraphErrors(_h_dt_met_para_shift[3]);
    _gr_mc_met_para_shift[3] = new TGraphErrors(_h_mc_met_para_shift[3]);
    _gr_dt_met_para_shift[4] = new TGraphErrors(_h_dt_met_para_shift[4]);
    _gr_mc_met_para_shift[4] = new TGraphErrors(_h_mc_met_para_shift[4]);
    _gr_dt_met_para_shift[5] = new TGraphErrors(_h_dt_met_para_shift[5]);
    _gr_mc_met_para_shift[5] = new TGraphErrors(_h_mc_met_para_shift[5]);

    _h_met_para_shift_dtmc[3] = (TH1D*)_h_met_para_shift_dtmc[0]->Clone("h_met_para_shift_dtmc_all_smooth");
    _h_met_para_shift_dtmc[4] = (TH1D*)_h_met_para_shift_dtmc[1]->Clone("h_met_para_shift_dtmc_mu_smooth");
    _h_met_para_shift_dtmc[5] = (TH1D*)_h_met_para_shift_dtmc[2]->Clone("h_met_para_shift_dtmc_el_smooth");
    _h_met_para_shift_dtmc[3]->Smooth();
    _h_met_para_shift_dtmc[4]->Smooth();
    _h_met_para_shift_dtmc[5]->Smooth();

    _gr_met_para_shift_dtmc[3] = new TGraphErrors(_h_met_para_shift_dtmc[3]);
    _gr_met_para_shift_dtmc[4] = new TGraphErrors(_h_met_para_shift_dtmc[4]);
    _gr_met_para_shift_dtmc[5] = new TGraphErrors(_h_met_para_shift_dtmc[5]);

    _h_ratio_met_para_sigma_dtmc[3] = (TH1D*)_h_ratio_met_para_sigma_dtmc[0]->Clone("h_ratio_met_para_sigma_dtmc_all_smooth");
    _h_ratio_met_para_sigma_dtmc[4] = (TH1D*)_h_ratio_met_para_sigma_dtmc[1]->Clone("h_ratio_met_para_sigma_dtmc_mu_smooth");
    _h_ratio_met_para_sigma_dtmc[5] = (TH1D*)_h_ratio_met_para_sigma_dtmc[2]->Clone("h_ratio_met_para_sigma_dtmc_el_smooth");
    _h_ratio_met_para_sigma_dtmc[3]->Smooth();
    _h_ratio_met_para_sigma_dtmc[4]->Smooth();
    _h_ratio_met_para_sigma_dtmc[5]->Smooth();

    _h_ratio_met_perp_sigma_dtmc[3] = (TH1D*)_h_ratio_met_perp_sigma_dtmc[0]->Clone("h_ratio_met_perp_sigma_dtmc_all_smooth");
    _h_ratio_met_perp_sigma_dtmc[4] = (TH1D*)_h_ratio_met_perp_sigma_dtmc[1]->Clone("h_ratio_met_perp_sigma_dtmc_mu_smooth");
    _h_ratio_met_perp_sigma_dtmc[5] = (TH1D*)_h_ratio_met_perp_sigma_dtmc[2]->Clone("h_ratio_met_perp_sigma_dtmc_el_smooth");
    _h_ratio_met_perp_sigma_dtmc[3]->Smooth();
    _h_ratio_met_perp_sigma_dtmc[4]->Smooth();
    _h_ratio_met_perp_sigma_dtmc[5]->Smooth();

    _gr_ratio_met_para_sigma_dtmc[3] = new TGraphErrors(_h_ratio_met_para_sigma_dtmc[3]);
    _gr_ratio_met_para_sigma_dtmc[4] = new TGraphErrors(_h_ratio_met_para_sigma_dtmc[4]);
    _gr_ratio_met_para_sigma_dtmc[5] = new TGraphErrors(_h_ratio_met_para_sigma_dtmc[5]);
    _gr_ratio_met_perp_sigma_dtmc[3] = new TGraphErrors(_h_ratio_met_perp_sigma_dtmc[3]);
    _gr_ratio_met_perp_sigma_dtmc[4] = new TGraphErrors(_h_ratio_met_perp_sigma_dtmc[4]);
    _gr_ratio_met_perp_sigma_dtmc[5] = new TGraphErrors(_h_ratio_met_perp_sigma_dtmc[5]);
  }
  
  

}

// do simple met recoil fine tuning
void doRecoil()
{
  if ( _doRecoil && !_isData && _isDyJets) {
    Float_t met_para = _llnunu_l2_pt*cos(_llnunu_l2_phi-_llnunu_l1_phi);
    Float_t met_perp = _llnunu_l2_pt*sin(_llnunu_l2_phi-_llnunu_l1_phi);
    if (abs(_llnunu_l1_l1_pdgId)==13&&abs(_llnunu_l1_l2_pdgId)==13) {
      Int_t idd=1;
      if (_doRecoilUseSmooth) idd=4;
      _h_dt_met_para_shift[6] = _h_dt_met_para_shift[idd];
      _h_mc_met_para_shift[6] = _h_mc_met_para_shift[idd];
      _h_met_para_shift_dtmc[6] = _h_met_para_shift_dtmc[idd];
      _h_ratio_met_para_sigma_dtmc[6] = _h_ratio_met_para_sigma_dtmc[idd];
      _h_ratio_met_perp_sigma_dtmc[6] = _h_ratio_met_perp_sigma_dtmc[idd];
      _gr_dt_met_para_shift[6] = _gr_dt_met_para_shift[4];
      _gr_mc_met_para_shift[6] = _gr_mc_met_para_shift[4];
      _gr_met_para_shift_dtmc[6] = _gr_met_para_shift_dtmc[4];
      _gr_ratio_met_para_sigma_dtmc[6] = _gr_ratio_met_para_sigma_dtmc[4];
      _gr_ratio_met_perp_sigma_dtmc[6] = _gr_ratio_met_perp_sigma_dtmc[4];
    }
    else if (abs(_llnunu_l1_l1_pdgId)==11&&abs(_llnunu_l1_l2_pdgId)==11) {
      Int_t idd=2;
      if (_doRecoilUseSmooth) idd=5;
      _h_dt_met_para_shift[6] = _h_dt_met_para_shift[idd];
      _h_mc_met_para_shift[6] = _h_mc_met_para_shift[idd];
      _h_met_para_shift_dtmc[6] = _h_met_para_shift_dtmc[idd];
      _h_ratio_met_para_sigma_dtmc[6] = _h_ratio_met_para_sigma_dtmc[idd];
      _h_ratio_met_perp_sigma_dtmc[6] = _h_ratio_met_perp_sigma_dtmc[idd];
      _gr_dt_met_para_shift[6] = _gr_dt_met_para_shift[5];
      _gr_mc_met_para_shift[6] = _gr_mc_met_para_shift[5];
      _gr_met_para_shift_dtmc[6] = _gr_met_para_shift_dtmc[5];
      _gr_ratio_met_para_sigma_dtmc[6] = _gr_ratio_met_para_sigma_dtmc[5];
      _gr_ratio_met_perp_sigma_dtmc[6] = _gr_ratio_met_perp_sigma_dtmc[5];
    }
    else {
      Int_t idd=0;
      if (_doRecoilUseSmooth) idd=3;
      _h_dt_met_para_shift[6] = _h_dt_met_para_shift[idd];
      _h_mc_met_para_shift[6] = _h_mc_met_para_shift[idd];
      _h_met_para_shift_dtmc[6] = _h_met_para_shift_dtmc[idd];
      _h_ratio_met_para_sigma_dtmc[6] = _h_ratio_met_para_sigma_dtmc[idd];
      _h_ratio_met_perp_sigma_dtmc[6] = _h_ratio_met_perp_sigma_dtmc[idd];
      _gr_dt_met_para_shift[6] = _gr_dt_met_para_shift[3];
      _gr_mc_met_para_shift[6] = _gr_mc_met_para_shift[3];
      _gr_met_para_shift_dtmc[6] = _gr_met_para_shift_dtmc[3];
      _gr_ratio_met_para_sigma_dtmc[6] = _gr_ratio_met_para_sigma_dtmc[3];
      _gr_ratio_met_perp_sigma_dtmc[6] = _gr_ratio_met_perp_sigma_dtmc[3];
    }

    // peak shift
    if (_doRecoilUseSmoothGraph) {
      met_para += _gr_met_para_shift_dtmc[6]->Eval(_llnunu_l1_pt);
    }
    else {
      met_para += _h_met_para_shift_dtmc[6]->GetBinContent(_h_met_para_shift_dtmc[6]->FindBin(_llnunu_l1_pt));
    }

    // smearing
    if (_doRecoilUseSmoothGraph) {
      met_para = (met_para-_gr_dt_met_para_shift[6]->Eval(_llnunu_l1_pt))*_gr_ratio_met_para_sigma_dtmc[6]->Eval(_llnunu_l1_pt) + _gr_dt_met_para_shift[6]->Eval(_llnunu_l1_pt);
      met_perp *= _gr_ratio_met_perp_sigma_dtmc[6]->Eval(_llnunu_l1_pt);
    } 
    else {
      met_para = (met_para-_h_dt_met_para_shift[6]->GetBinContent(_h_dt_met_para_shift[6]->FindBin(_llnunu_l1_pt)))
               * _h_ratio_met_para_sigma_dtmc[6]->GetBinContent(_h_ratio_met_para_sigma_dtmc[6]->FindBin(_llnunu_l1_pt)) 
               + _h_dt_met_para_shift[6]->GetBinContent(_h_dt_met_para_shift[6]->FindBin(_llnunu_l1_pt));
      met_perp *= _h_ratio_met_perp_sigma_dtmc[6]->GetBinContent(_h_ratio_met_perp_sigma_dtmc[6]->FindBin(_llnunu_l1_pt));
    }

    // recalculate vars
    Float_t met_px = met_para*cos(_llnunu_l1_phi)-met_perp*sin(_llnunu_l1_phi);
    Float_t met_py = met_para*sin(_llnunu_l1_phi)+met_perp*cos(_llnunu_l1_phi);
    TVector2 vec_met;
    vec_met.Set(met_px, met_py);
    _llnunu_l2_pt = vec_met.Mod();
    _llnunu_l2_phi = TVector2::Phi_mpi_pi(vec_met.Phi());
    
    Float_t et1 = TMath::Sqrt(_llnunu_l1_mass*_llnunu_l1_mass + _llnunu_l1_pt*_llnunu_l1_pt);
    Float_t et2 = TMath::Sqrt(_llnunu_l1_mass*_llnunu_l1_mass + _llnunu_l2_pt*_llnunu_l2_pt);
    _llnunu_mt = TMath::Sqrt(2.0*_llnunu_l1_mass*_llnunu_l1_mass + 2.0* (et1*et2 - _llnunu_l1_pt*cos(_llnunu_l1_phi)*met_px - _llnunu_l1_pt*sin(_llnunu_l1_phi)*met_py));
  }

}



// prepareEffScale
void prepareEffScale()
{

  // needed branches 
  _tree_out->Branch("trgsf", &_trgsf, "trgsf/F");
  _tree_out->Branch("isosf", &_isosf, "isosf/F");
  _tree_out->Branch("idsf", &_idsf, "idsf/F");
  _tree_out->Branch("trksf", &_trksf, "trksf/F");
  _tree_out->Branch("trgsf_err", &_trgsf_err, "trgsf_err/F");
  _tree_out->Branch("isosf_err", &_isosf_err, "isosf_err/F");
  _tree_out->Branch("idsf_err", &_idsf_err, "idsf_err/F");
  _tree_out->Branch("trksf_err", &_trksf_err, "trksf_err/F");
  _tree_out->Branch("trgsf_up", &_trgsf_up, "trgsf_up/F");
  _tree_out->Branch("trgsf_dn", &_trgsf_dn, "trgsf_dn/F");
  _tree_out->Branch("idisotrksf", &_idisotrksf, "idisotrksf/F");
  _tree_out->Branch("idisotrksf_up", &_idisotrksf_up, "idisotrksf_up/F");
  _tree_out->Branch("idisotrksf_dn", &_idisotrksf_dn, "idisotrksf_dn/F");

  // Electron ID ISO scale factors 
  _file_idiso_el = TFile::Open(_EffScaleInputFileName_IdIso_El.c_str());
  _h_sf_idiso_el = (TH2F*)_file_idiso_el->Get("EGamma_SF2D");

  // Electron tracking scale factors
  _file_trksf_el = TFile::Open(_EffScaleInputFileName_Trk_El.c_str());
  _h_sf_trk_el = (TH2F*)_file_trksf_el->Get("EGamma_SF2D");

  // muon tracking scale factors
  _file_trksf_mu = TFile::Open(_EffScaleInputFileName_Trk_Mu.c_str());
  _h_sf_trk_mu = (TH1F*)_file_trksf_mu->Get("hist_ratio_eta");

  // muon id iso scale factors
  _file_idiso_mu = TFile::Open(_EffScaleInputFileName_IdIso_Mu.c_str());
  _h_eff_trkhpt_mu_dt = (TH2F*)_file_idiso_mu->Get("eff_trackHighPt_80Xdata_pteta");
  _h_eff_trkhpt_mu_mc = (TH2F*)_file_idiso_mu->Get("eff_trackHighPt_80Xmc_pteta");
  _h_eff_hpt_mu_dt = (TH2F*)_file_idiso_mu->Get("eff_HighPt_80Xdata_pteta");
  _h_eff_hpt_mu_mc = (TH2F*)_file_idiso_mu->Get("eff_HighPt_80Xmc_pteta");
  _h_sf_iso_mu = (TH2F*)_file_idiso_mu->Get("sf_trackerIso_80X_pteta");

  // electron trigger scale factors
  _file_trg_el = TFile::Open(_EffScaleInputFileName_Trg_El.c_str());
  _h_sf_trg_el_l1=(TH2D*)_file_trg_el->Get("ell1pteta"); 

  // muon trigger scale factors
  _file_trg_mu = TFile::Open(_EffScaleInputFileName_Trg_Mu.c_str());
  _h_eff_trg_mu_l1_tot = (TH2D*)_file_trg_mu->Get("htrg_l1_tot");
  _h_eff_trg_mu_l2_tot = (TH2D*)_file_trg_mu->Get("htrg_l2_tot");
  _h_eff_trg_mu_l1_l1p = (TH2D*)_file_trg_mu->Get("htrg_l1_l1p");
  _h_eff_trg_mu_l2_l1p = (TH2D*)_file_trg_mu->Get("htrg_l2_l1p");
  _h_eff_trg_mu_l1_l1f = (TH2D*)_file_trg_mu->Get("htrg_l1_l1f");
  _h_eff_trg_mu_l2_l1f = (TH2D*)_file_trg_mu->Get("htrg_l2_l1f");
  _h_eff_trg_mu_l1_l1pl2f = (TH2D*)_file_trg_mu->Get("htrg_l1_l1pl2f");
  _h_eff_trg_mu_l1_l1pl2p = (TH2D*)_file_trg_mu->Get("htrg_l1_l1pl2p");
  _h_eff_trg_mu_l1_l1fl2p = (TH2D*)_file_trg_mu->Get("htrg_l1_l1fl2p");
  _h_eff_trg_mu_l2_l1pl2f = (TH2D*)_file_trg_mu->Get("htrg_l2_l1pl2f");
  _h_eff_trg_mu_l2_l1pl2p = (TH2D*)_file_trg_mu->Get("htrg_l2_l1pl2p");
  _h_eff_trg_mu_l2_l1fl2p = (TH2D*)_file_trg_mu->Get("htrg_l2_l1fl2p");

  _NPtBins_eff_trg_mu = _h_eff_trg_mu_l2_tot->GetNbinsX();
  _NEtaBins_eff_trg_mu = _h_eff_trg_mu_l2_tot->GetNbinsY();
  _N_eff_trg_mu_tot = _h_eff_trg_mu_l2_tot->IntegralAndError(Int_t(1), _NPtBins_eff_trg_mu, Int_t(1), _NEtaBins_eff_trg_mu, _N_eff_trg_mu_tot_err);
  _N_eff_trg_mu_l1p = _h_eff_trg_mu_l2_l1p->IntegralAndError(Int_t(1), _NPtBins_eff_trg_mu, Int_t(1), _NEtaBins_eff_trg_mu, _N_eff_trg_mu_l1p_err);
  _N_eff_trg_mu_l1f = _h_eff_trg_mu_l2_l1f->IntegralAndError(Int_t(1), _NPtBins_eff_trg_mu, Int_t(1), _NEtaBins_eff_trg_mu, _N_eff_trg_mu_l1f_err);
  _N_eff_trg_mu_l1pl2f = _h_eff_trg_mu_l2_l1pl2f->IntegralAndError(Int_t(1), _NPtBins_eff_trg_mu, Int_t(1), _NEtaBins_eff_trg_mu, _N_eff_trg_mu_l1pl2f_err);
  _N_eff_trg_mu_l1pl2p = _h_eff_trg_mu_l2_l1pl2p->IntegralAndError(Int_t(1), _NPtBins_eff_trg_mu, Int_t(1), _NEtaBins_eff_trg_mu, _N_eff_trg_mu_l1pl2p_err);
  _N_eff_trg_mu_l1fl2p = _h_eff_trg_mu_l2_l1fl2p->IntegralAndError(Int_t(1), _NPtBins_eff_trg_mu, Int_t(1), _NEtaBins_eff_trg_mu, _N_eff_trg_mu_l1fl2p_err);

  _h_eff_trg_mu_l1_tot_norm = (TH2D*)_h_eff_trg_mu_l1_tot->Clone("htrg_l1_tot_norm");
  _h_eff_trg_mu_l2_tot_norm = (TH2D*)_h_eff_trg_mu_l2_tot->Clone("htrg_l2_tot_norm");
  _h_eff_trg_mu_l1_l1p_norm = (TH2D*)_h_eff_trg_mu_l1_l1p->Clone("htrg_l1_l1p_norm");
  _h_eff_trg_mu_l1_l1f_norm = (TH2D*)_h_eff_trg_mu_l1_l1f->Clone("htrg_l1_l1f_norm");
  _h_eff_trg_mu_l2_l1p_norm = (TH2D*)_h_eff_trg_mu_l2_l1p->Clone("htrg_l2_l1p_norm");
  _h_eff_trg_mu_l2_l1f_norm = (TH2D*)_h_eff_trg_mu_l2_l1f->Clone("htrg_l2_l1f_norm");
  _h_eff_trg_mu_l1_l1pl2f_norm = (TH2D*)_h_eff_trg_mu_l1_l1pl2f->Clone("htrg_l1_l1pl2f_norm");
  _h_eff_trg_mu_l1_l1pl2p_norm = (TH2D*)_h_eff_trg_mu_l1_l1pl2p->Clone("htrg_l1_l1pl2p_norm");
  _h_eff_trg_mu_l1_l1fl2p_norm = (TH2D*)_h_eff_trg_mu_l1_l1fl2p->Clone("htrg_l1_l1fl2p_norm");
  _h_eff_trg_mu_l2_l1pl2f_norm = (TH2D*)_h_eff_trg_mu_l2_l1pl2f->Clone("htrg_l2_l1pl2f_norm");
  _h_eff_trg_mu_l2_l1pl2p_norm = (TH2D*)_h_eff_trg_mu_l2_l1pl2p->Clone("htrg_l2_l1pl2p_norm");
  _h_eff_trg_mu_l2_l1fl2p_norm = (TH2D*)_h_eff_trg_mu_l2_l1fl2p->Clone("htrg_l2_l1fl2p_norm");


  _h_eff_trg_mu_l1_tot_norm->Scale(1./_N_eff_trg_mu_tot);
  _h_eff_trg_mu_l2_tot_norm->Scale(1./_N_eff_trg_mu_tot);
  _h_eff_trg_mu_l1_l1p_norm->Scale(1./_N_eff_trg_mu_l1p);
  _h_eff_trg_mu_l1_l1f_norm->Scale(1./_N_eff_trg_mu_l1f);
  _h_eff_trg_mu_l2_l1p_norm->Scale(1./_N_eff_trg_mu_l1p);
  _h_eff_trg_mu_l2_l1f_norm->Scale(1./_N_eff_trg_mu_l1f);
  _h_eff_trg_mu_l1_l1pl2f_norm->Scale(1./_N_eff_trg_mu_l1pl2f);
  _h_eff_trg_mu_l1_l1pl2p_norm->Scale(1./_N_eff_trg_mu_l1pl2p);
  _h_eff_trg_mu_l1_l1fl2p_norm->Scale(1./_N_eff_trg_mu_l1fl2p);
  _h_eff_trg_mu_l2_l1pl2f_norm->Scale(1./_N_eff_trg_mu_l1pl2f);
  _h_eff_trg_mu_l2_l1pl2p_norm->Scale(1./_N_eff_trg_mu_l1pl2p);
  _h_eff_trg_mu_l2_l1fl2p_norm->Scale(1./_N_eff_trg_mu_l1fl2p);

  _h_eff_trg_mu_l1_l1p_norm_vs_tot    = (TH2D*)_h_eff_trg_mu_l1_l1p_norm->Clone("htrg_l1_l1p_norm_vs_tot");
  _h_eff_trg_mu_l1_l1f_norm_vs_tot    = (TH2D*)_h_eff_trg_mu_l1_l1f_norm->Clone("htrg_l1_l1f_norm_vs_tot");
  _h_eff_trg_mu_l2_l1p_norm_vs_tot    = (TH2D*)_h_eff_trg_mu_l2_l1p_norm->Clone("htrg_l2_l1p_norm_vs_tot");
  _h_eff_trg_mu_l2_l1f_norm_vs_tot    = (TH2D*)_h_eff_trg_mu_l2_l1f_norm->Clone("htrg_l2_l1f_norm_vs_tot");
  _h_eff_trg_mu_l1_l1pl2f_norm_vs_tot = (TH2D*)_h_eff_trg_mu_l1_l1pl2f_norm->Clone("htrg_l1_l1pl2f_norm_vs_tot");
  _h_eff_trg_mu_l1_l1pl2p_norm_vs_tot = (TH2D*)_h_eff_trg_mu_l1_l1pl2p_norm->Clone("htrg_l1_l1pl2p_norm_vs_tot");
  _h_eff_trg_mu_l1_l1fl2p_norm_vs_tot = (TH2D*)_h_eff_trg_mu_l1_l1fl2p_norm->Clone("htrg_l1_l1fl2p_norm_vs_tot");
  _h_eff_trg_mu_l2_l1pl2f_norm_vs_tot = (TH2D*)_h_eff_trg_mu_l2_l1pl2f_norm->Clone("htrg_l2_l1pl2f_norm_vs_tot");
  _h_eff_trg_mu_l2_l1pl2p_norm_vs_tot = (TH2D*)_h_eff_trg_mu_l2_l1pl2p_norm->Clone("htrg_l2_l1pl2p_norm_vs_tot");
  _h_eff_trg_mu_l2_l1fl2p_norm_vs_tot = (TH2D*)_h_eff_trg_mu_l2_l1fl2p_norm->Clone("htrg_l2_l1fl2p_norm_vs_tot");
  _h_eff_trg_mu_l1_l1pl2f_norm_vs_l1p = (TH2D*)_h_eff_trg_mu_l1_l1pl2f_norm->Clone("htrg_l1_l1pl2f_norm_vs_l1p");
  _h_eff_trg_mu_l1_l1pl2p_norm_vs_l1p = (TH2D*)_h_eff_trg_mu_l1_l1pl2p_norm->Clone("htrg_l1_l1pl2p_norm_vs_l1p");
  _h_eff_trg_mu_l1_l1fl2p_norm_vs_l1f = (TH2D*)_h_eff_trg_mu_l1_l1fl2p_norm->Clone("htrg_l1_l1fl2p_norm_vs_l1f");
  _h_eff_trg_mu_l2_l1pl2f_norm_vs_l1p = (TH2D*)_h_eff_trg_mu_l2_l1pl2f_norm->Clone("htrg_l2_l1pl2f_norm_vs_l1p");
  _h_eff_trg_mu_l2_l1pl2p_norm_vs_l1p = (TH2D*)_h_eff_trg_mu_l2_l1pl2p_norm->Clone("htrg_l2_l1pl2p_norm_vs_l1p");
  _h_eff_trg_mu_l2_l1fl2p_norm_vs_l1f = (TH2D*)_h_eff_trg_mu_l2_l1fl2p_norm->Clone("htrg_l2_l1fl2p_norm_vs_l1f");

  _h_eff_trg_mu_l1_l1p_norm_vs_tot->Divide(_h_eff_trg_mu_l1_tot_norm);
  _h_eff_trg_mu_l1_l1f_norm_vs_tot->Divide(_h_eff_trg_mu_l1_tot_norm);
  _h_eff_trg_mu_l2_l1p_norm_vs_tot->Divide(_h_eff_trg_mu_l2_tot_norm);
  _h_eff_trg_mu_l2_l1f_norm_vs_tot->Divide(_h_eff_trg_mu_l2_tot_norm);
  _h_eff_trg_mu_l1_l1pl2f_norm_vs_tot->Divide(_h_eff_trg_mu_l1_tot_norm);
  _h_eff_trg_mu_l1_l1pl2p_norm_vs_tot->Divide(_h_eff_trg_mu_l1_tot_norm);
  _h_eff_trg_mu_l1_l1fl2p_norm_vs_tot->Divide(_h_eff_trg_mu_l1_tot_norm);
  _h_eff_trg_mu_l2_l1pl2f_norm_vs_tot->Divide(_h_eff_trg_mu_l2_tot_norm);
  _h_eff_trg_mu_l2_l1pl2p_norm_vs_tot->Divide(_h_eff_trg_mu_l2_tot_norm);
  _h_eff_trg_mu_l2_l1fl2p_norm_vs_tot->Divide(_h_eff_trg_mu_l2_tot_norm);
  _h_eff_trg_mu_l1_l1pl2f_norm_vs_l1p->Divide(_h_eff_trg_mu_l1_l1p_norm);
  _h_eff_trg_mu_l1_l1pl2p_norm_vs_l1p->Divide(_h_eff_trg_mu_l1_l1p_norm);
  _h_eff_trg_mu_l1_l1fl2p_norm_vs_l1f->Divide(_h_eff_trg_mu_l1_l1f_norm);
  _h_eff_trg_mu_l2_l1pl2f_norm_vs_l1p->Divide(_h_eff_trg_mu_l2_l1p_norm);
  _h_eff_trg_mu_l2_l1pl2p_norm_vs_l1p->Divide(_h_eff_trg_mu_l2_l1p_norm);
  _h_eff_trg_mu_l2_l1fl2p_norm_vs_l1f->Divide(_h_eff_trg_mu_l2_l1f_norm);

 

}


// add efficiency scale factors
void addEffScale()
{
  // muon 
  if (abs(_llnunu_l1_l1_pdgId)==13 && abs(_llnunu_l1_l2_pdgId)==13) {

    // id 
    double effdt1  = _h_eff_trkhpt_mu_dt->GetBinContent(_h_eff_trkhpt_mu_dt->FindBin(_llnunu_l1_l1_eta, _llnunu_l1_l1_pt));
    double effmc1  = _h_eff_trkhpt_mu_mc->GetBinContent(_h_eff_trkhpt_mu_mc->FindBin(_llnunu_l1_l1_eta, _llnunu_l1_l1_pt));
    double errdt1  = _h_eff_trkhpt_mu_dt->GetBinError(_h_eff_trkhpt_mu_dt->FindBin(_llnunu_l1_l1_eta, _llnunu_l1_l1_pt));
    double errmc1  = _h_eff_trkhpt_mu_mc->GetBinError(_h_eff_trkhpt_mu_mc->FindBin(_llnunu_l1_l1_eta, _llnunu_l1_l1_pt));
    double effdt1a = _h_eff_hpt_mu_dt->GetBinContent(_h_eff_hpt_mu_dt->FindBin(_llnunu_l1_l1_eta, _llnunu_l1_l1_pt));
    double effmc1a = _h_eff_hpt_mu_mc->GetBinContent(_h_eff_hpt_mu_mc->FindBin(_llnunu_l1_l1_eta, _llnunu_l1_l1_pt));
    double errdt1a = _h_eff_hpt_mu_dt->GetBinError(_h_eff_hpt_mu_dt->FindBin(_llnunu_l1_l1_eta, _llnunu_l1_l1_pt));
    double errmc1a = _h_eff_hpt_mu_mc->GetBinError(_h_eff_hpt_mu_mc->FindBin(_llnunu_l1_l1_eta, _llnunu_l1_l1_pt));
    double effdt2  = _h_eff_trkhpt_mu_dt->GetBinContent(_h_eff_trkhpt_mu_dt->FindBin(_llnunu_l1_l2_eta, _llnunu_l1_l2_pt));
    double effmc2  = _h_eff_trkhpt_mu_mc->GetBinContent(_h_eff_trkhpt_mu_mc->FindBin(_llnunu_l1_l2_eta, _llnunu_l1_l2_pt));
    double errdt2  = _h_eff_trkhpt_mu_dt->GetBinError(_h_eff_trkhpt_mu_dt->FindBin(_llnunu_l1_l2_eta, _llnunu_l1_l2_pt));
    double errmc2  = _h_eff_trkhpt_mu_mc->GetBinError(_h_eff_trkhpt_mu_mc->FindBin(_llnunu_l1_l2_eta, _llnunu_l1_l2_pt));
    double effdt2a = _h_eff_hpt_mu_dt->GetBinContent(_h_eff_hpt_mu_dt->FindBin(_llnunu_l1_l2_eta, _llnunu_l1_l2_pt));
    double effmc2a = _h_eff_hpt_mu_mc->GetBinContent(_h_eff_hpt_mu_mc->FindBin(_llnunu_l1_l2_eta, _llnunu_l1_l2_pt));
    double errdt2a = _h_eff_hpt_mu_dt->GetBinError(_h_eff_hpt_mu_dt->FindBin(_llnunu_l1_l2_eta, _llnunu_l1_l2_pt));
    double errmc2a = _h_eff_hpt_mu_mc->GetBinError(_h_eff_hpt_mu_mc->FindBin(_llnunu_l1_l2_eta, _llnunu_l1_l2_pt));
    double effdt = effdt1*effdt2a+effdt1a*effdt2-effdt1a*effdt2a;
    double effmc = effmc1*effmc2a+effmc1a*effmc2-effmc1a*effmc2a;
    if(effdt&&effmc){
      _idsf = effdt/effmc;
      _idsf_err = (  TMath::Power((effdt2-effdt2a)*errdt1a,2)
                  + TMath::Power((effdt1-effdt1a)*errdt2a,2)
                  + TMath::Power(effdt1a*errdt2,2)
                  + TMath::Power(effdt2a*errdt1,2) 
                 ) / TMath::Power(effdt,2) 
               + (  TMath::Power((effmc2-effmc2a)*errmc1a,2)
                  + TMath::Power((effmc1-effmc1a)*errmc2a,2)
                  + TMath::Power(effmc1a*errmc2,2)
                  + TMath::Power(effmc2a*errmc1,2)
                 ) / TMath::Power(effmc,2);

      _idsf_err = TMath::Power(_idsf_err,.5)*_idsf;
    }
    else {
      _idsf = 1;
      _idsf_err = 1;
    }


    // iso
    effdt1 = _h_sf_iso_mu->GetBinContent(_h_sf_iso_mu->FindBin(_llnunu_l1_l1_eta,_llnunu_l1_l1_pt));
    effdt2 = _h_sf_iso_mu->GetBinContent(_h_sf_iso_mu->FindBin(_llnunu_l1_l2_eta,_llnunu_l1_l2_pt));
    errdt1 = _h_sf_iso_mu->GetBinError(_h_sf_iso_mu->FindBin(_llnunu_l1_l1_eta,_llnunu_l1_l1_pt));
    errdt2 = _h_sf_iso_mu->GetBinError(_h_sf_iso_mu->FindBin(_llnunu_l1_l2_eta,_llnunu_l1_l2_pt));
    _isosf = effdt1*effdt2;
    _isosf_err = TMath::Power((TMath::Power(effdt1*errdt2,2) + TMath::Power(errdt1*effdt2,2)), .5);

    // tracking
    effdt1 = _h_sf_trk_mu->GetBinContent(_h_sf_trk_mu->FindBin(_llnunu_l1_l1_eta));
    effdt2 = _h_sf_trk_mu->GetBinContent(_h_sf_trk_mu->FindBin(_llnunu_l1_l2_eta));
    errdt1 = _h_sf_trk_mu->GetBinError(_h_sf_trk_mu->FindBin(_llnunu_l1_l1_eta));
    errdt2 = _h_sf_trk_mu->GetBinError(_h_sf_trk_mu->FindBin(_llnunu_l1_l2_eta));
    _trksf = effdt1*effdt2;
    _trksf_err = TMath::Power((TMath::Power(effdt1*errdt2,2) + TMath::Power(errdt1*effdt2,2)), .5);

    // id, iso, tracking combined uncertainty up/down
    double idisotrksf_err = sqrt(_idsf_err*_idsf_err+_isosf_err*_isosf_err+_trksf_err*_trksf_err);
    _idisotrksf = _idsf*_isosf*_trksf;
    _idisotrksf_up = _idisotrksf+0.5*idisotrksf_err;
    _idisotrksf_dn = _idisotrksf-0.5*idisotrksf_err;

    // trigger
    int trg_bin_l1 = _h_eff_trg_mu_l1_l1p_norm_vs_tot->FindBin(_llnunu_l1_l1_pt,fabs(_llnunu_l1_l1_eta));
    int trg_bin_l2 = _h_eff_trg_mu_l2_l1pl2f_norm_vs_l1p->FindBin(_llnunu_l1_l2_pt,fabs(_llnunu_l1_l2_eta));
    double trg_sc_l1_l1p_vs_tot = _h_eff_trg_mu_l1_l1p_norm_vs_tot->GetBinContent(trg_bin_l1);
    double trg_sc_l2_l1pl2f_vs_l1p = _h_eff_trg_mu_l2_l1pl2f_norm_vs_l1p->GetBinContent(trg_bin_l2);
    double trg_sc_l2_l1pl2p_vs_l1p = _h_eff_trg_mu_l2_l1pl2p_norm_vs_l1p->GetBinContent(trg_bin_l2);
    double trg_sc_l2_l1fl2p_vs_tot = _h_eff_trg_mu_l2_l1fl2p_norm_vs_tot->GetBinContent(trg_bin_l2);
    double trg_sc_l1_l1p_vs_tot_err = _h_eff_trg_mu_l1_l1p_norm_vs_tot->GetBinError(trg_bin_l1);
    double trg_sc_l2_l1pl2f_vs_l1p_err = _h_eff_trg_mu_l2_l1pl2f_norm_vs_l1p->GetBinError(trg_bin_l2);
    double trg_sc_l2_l1pl2p_vs_l1p_err = _h_eff_trg_mu_l2_l1pl2p_norm_vs_l1p->GetBinError(trg_bin_l2);
    double trg_sc_l2_l1fl2p_vs_tot_err = _h_eff_trg_mu_l2_l1fl2p_norm_vs_tot->GetBinError(trg_bin_l2);

    double trg_npass = _N_eff_trg_mu_l1pl2f*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2f_vs_l1p
                     + _N_eff_trg_mu_l1pl2p*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2p_vs_l1p
                     + _N_eff_trg_mu_l1fl2p*trg_sc_l2_l1fl2p_vs_tot
                     ;
    double trg_npass_err = pow(_N_eff_trg_mu_l1pl2f_err*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2f_vs_l1p,2)
                         + pow(_N_eff_trg_mu_l1pl2f*trg_sc_l1_l1p_vs_tot_err*trg_sc_l2_l1pl2f_vs_l1p,2)
                         + pow(_N_eff_trg_mu_l1pl2f*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2f_vs_l1p_err,2)
                         + pow(_N_eff_trg_mu_l1pl2p_err*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2p_vs_l1p,2)
                         + pow(_N_eff_trg_mu_l1pl2p*trg_sc_l1_l1p_vs_tot_err*trg_sc_l2_l1pl2p_vs_l1p,2)
                         + pow(_N_eff_trg_mu_l1pl2p*trg_sc_l1_l1p_vs_tot*trg_sc_l2_l1pl2p_vs_l1p_err,2)
                         + pow(_N_eff_trg_mu_l1fl2p_err*trg_sc_l2_l1fl2p_vs_tot,2)
                         + pow(_N_eff_trg_mu_l1fl2p*trg_sc_l2_l1fl2p_vs_tot_err,2)
                         ;
    trg_npass_err = sqrt(trg_npass_err);

    double trg_nfail = _N_eff_trg_mu_tot-trg_npass;
    double trg_nfail_err = sqrt(_N_eff_trg_mu_tot_err*_N_eff_trg_mu_tot_err
                                 - _N_eff_trg_mu_l1pl2f_err*_N_eff_trg_mu_l1pl2f_err
                                 - _N_eff_trg_mu_l1pl2p_err*_N_eff_trg_mu_l1pl2p_err
                                 - _N_eff_trg_mu_l1fl2p_err*_N_eff_trg_mu_l1fl2p_err);

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

    _trgsf = trg_eff;
    _trgsf_err = trg_eff_err;
    _trgsf_up = trg_eff_up;
    _trgsf_dn = trg_eff_dn;

  }
  // electron
  else if (abs(_llnunu_l1_l1_pdgId)==11 && abs(_llnunu_l1_l2_pdgId)==11) {

    // id
    double effdt1 = 1.0;
    if(_llnunu_l1_l1_pt<=200) {
      effdt1 = _h_sf_idiso_el->GetBinContent(_h_sf_idiso_el->FindBin(_llnunu_l1_l1_eSCeta,_llnunu_l1_l1_pt));
    }
    double effdt2 = 1.0;
    if(_llnunu_l1_l2_pt<=200) {
      effdt2 = _h_sf_idiso_el->GetBinContent(_h_sf_idiso_el->FindBin(_llnunu_l1_l2_eSCeta,_llnunu_l1_l2_pt));
    }
    double errdt1 = _h_sf_idiso_el->GetBinError(_h_sf_idiso_el->FindBin(_llnunu_l1_l1_eSCeta,_llnunu_l1_l1_pt));
    double errdt2 = _h_sf_idiso_el->GetBinError(_h_sf_idiso_el->FindBin(_llnunu_l1_l2_eSCeta,_llnunu_l1_l2_pt));
    _idsf = effdt1*effdt2;
    _idsf_err = TMath::Power((TMath::Power(effdt1*errdt2,2)+TMath::Power(errdt1*effdt2,2)),.5);

    // iso (iso included in id sf)
    _isosf=1;
    _isosf_err=0;

    // track
    effdt1 = _h_sf_trk_el->GetBinContent(_h_sf_trk_el->FindBin(_llnunu_l1_l1_eSCeta,100));
    effdt2 = _h_sf_trk_el->GetBinContent(_h_sf_trk_el->FindBin(_llnunu_l1_l2_eSCeta,100));
    errdt1 = _h_sf_trk_el->GetBinError(_h_sf_trk_el->FindBin(_llnunu_l1_l1_eSCeta,100));
    errdt2 = _h_sf_trk_el->GetBinError(_h_sf_trk_el->FindBin(_llnunu_l1_l2_eSCeta,100));
    _trksf = effdt1*effdt2;
    _trksf_err = TMath::Power((TMath::Power(effdt1*errdt2,2)+TMath::Power(errdt1*effdt2,2)),.5);


    // id, iso, tracking combined uncertainty up/down
    double idisotrksf_err = sqrt(_idsf_err*_idsf_err+_isosf_err*_isosf_err+_trksf_err*_trksf_err);
    _idisotrksf = _idsf*_isosf*_trksf;
    _idisotrksf_up = _idisotrksf+0.5*idisotrksf_err;
    _idisotrksf_dn = _idisotrksf-0.5*idisotrksf_err;

    // trigger
    _trgsf = _h_sf_trg_el_l1->GetBinContent(_h_sf_trg_el_l1->FindBin(_llnunu_l1_l1_pt,fabs(_llnunu_l1_l1_eta)))/100;
    _trgsf_err = _h_sf_trg_el_l1->GetBinError(_h_sf_trg_el_l1->FindBin(_llnunu_l1_l1_pt,fabs(_llnunu_l1_l1_eta)))/100;

    _trgsf_up = _trgsf+0.5*_trgsf_err;
    _trgsf_dn = _trgsf-0.5*_trgsf_err;
    if (_trgsf>=1) _trgsf = 1;
    if (_trgsf<=0) _trgsf = 0;
    if (_trgsf_up>=1) _trgsf_up=1;
    if (_trgsf_up<=0) _trgsf_up=0;
    if (_trgsf_dn>=1) _trgsf_dn=1;
    if (_trgsf_dn<=0) _trgsf_dn=0;

    _trgsf_err = fabs(_trgsf_up-_trgsf_dn);



  }

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  ,.  .:;::;@;,;;:,::...@;,@,                                                                                                                 
//  ..,,.,,.,,:.,':,:::.`.;::':                                                                                                                 
//   ,``. .`,;'.,;::::'`     +@                                                                                                                 
//  :.```.``..:;:'::.::`,@@@+;,:                                                                                                                
//  ,`. :````.,.:..,`';,,    .@++                                                                                                               
//  ```.,``.`..`.```    .:@@@@#,``                                                                                                              
//    `,```` :..,.;@@@             +                                                                                                            
//   ,,,`...,:,,:    .+@@'                                                                                                                      
//  `,,,`,``.:::;@        `:;',,,,:',                                                                                                           
//  ,,.```.`:::.':@@@@@@@@@;,  `@'::.,                                                                                                          
//  ,`,````.::,,;:,.`@+#@@@@';`` ;@.,,                                                                                                          
//  ,`...`..,,.,;,:, @ .,,@@@@@@++. ::,:                                                  `                                                     
//  ``..``.`;:`,;::. @ `..':@@@@++''+@',`;                                  `     .        `                                                    
//  .``,....;,:,:,;,`@ `.`'.;+@   @@+';':;.,                            `         `.,  ``.,.     :,                                             
//  `...```.;:,`:,:. @ ```+`;.@     @@@@@@',`,.                    ,.,,.,.    `   `,;.':':.```      ``                                          
//  :,,`.` ,,:..:`. `@ ` .'`:.@      @@#+'''+#```.               `.::::;;'.,,;` ..',,::,.```                                                    
//  :,,.,`.`,;, ;````@ ``:'.:,'       :@++';::,..``,             ``: ``;;''';.;++#,.`  ,.`  `.:'..,`` ``                                        
//  ,.,,,``  ,: '` ` ' ,,:;:,:'         @+++';:,,.``:@                  `,;;;';;::` `      .  `,  `.`:,`..                                      
//  .```  `` ```..`,:,,::;::,,;          @+#+''::,.. ``.        `          .:.    ,    `     ` `. ``..`,.``                                     
//  `   `...``` ..` .,,,,,::.,,           @+##+';:,,..```,       ``     ``..,.`     `       ` `  .` `..`,.`                                     
//  ` .``````   `,`.` ;`,.::..,`           ;@#@@@'#@:,,.```.   .@@.`    `..,,:::#@'`.`,    `  ````..``,`.,.                                     
//  ` ``` ,  :@@ .; ,.;`,`,,`,,`             @'@@@+@@@@,,.```. ,; .;@:``....,:;:;'#@;`  .  .. `.`.```..,`..                                     
//  ,    `,, :,:# ',`,;., ,`,:,`              `@###@@'@@':';,,`  :'#@@@@@,`,,:;'+@':;,` +.`.,.`.....``..` ``                                    
//  ,    `` ;`@'.:' ;`:.,  :,,,.                ,@##@@@#@@+@` `,...`,@#`:;:+,::'+'+#:.` :,,.,,.``,.``.`.` `.                                    
//  .  `,  `.``,.,;.`.`,::@+,,,,                   @@@@@#+;;#@,.:':`` '@@@@#';;''@@,`. ,;;:,`:,.```` ``````.`                                   
//   @.`,  :.``` :,,::``@` :.;,;                      @@##':.,@@@@;.`.....`,:;;++@:.;``,:::,;;.        `. `.                                    
//  `..`, ,,```` , .`:.`@``+:`,@                          @@`@. `..,,,:::::;;+++@@:::`.,,`, . ,,.,;,``      `                                   
//    , .`;`.``. .`, ,,,@. ,,`:'                            @`,@:.,::;;'''++++#@@@,;,`. .  .  ``   .:.;.``                                      
//  ```,:,. ```` .```,,,@.,:,.:@                            :@@@:'';;'''+'+++#@@#@,,,`;        `, `  ` ;``..`                                   
//  .``.,..     ```, .`:@,,,:.:@                            @::;;;'''''''''+++#@@@@,  @;`   `.` ..     `,.`````                                 
//  .` `,.` ` ```.```` `@,,':,;+                              @@@@@@@@@+''++@@@@@@@,  @@ .` .:.`.,`.` ``````` ` `                               
//  .```,`` : ,.``.`.`` @`.:.','                                       @@@@@@@@@@@#'  @'` , `.,`...```   ```                                    
//  .```.`` . ,``.`....`@.`:...;                                       @@@@@@@@@@;@  .:,  ,.`,:`....``` `.                                      
//  ``` .`` ` ,``.````` @:`;.,.@                                      ,,  ````..,:,  `, ,``,.,.,.`...``````                                     
//:#,`` . . ` ,     .`.`@:`.`.,@                                   ..`        `.,,:` , `:. ,.,...`.,,. ```                                      
//':, ``, ,```, `` , `` ';`.,.,#                           .``         ``````..,,:::   ;., ,,;:.``..,; ````                                     
//,,, ` ``.   .` ` ,`,`.,+,,,,:,@@@,;;@:::`             ..`    `` ````````.,.:;';'': ` :, . ';,   :.,:.,`````                                   
//,.,`` `.,.```  ``. :`,,@.:'::,::,..:,:.,,,` `  .,,#@..``       .`    ``...,..,;;``',  ,   @,    :.,.....``                                    
//@#:`.``.,```.``````,`,,@,;,,:,::. :: `.: .:+,`....::,..`        ``````...,,::;;; @@,,`    ;      ',..``,. `                                   
//@@@````,, .``.`.`.`,`,,@::.,,,::..``'@, ',,@,``,:;;::,.`         `..`..,,::::;: @'+`, ,   .     .  .,,`,                                      
//+++``.`.,``,`..`..`;`:,@:..:,::@+'''';::@@,#@+:::;+;;:,.`         ``..,,,:;::,';.# @ # :   .   : `. ..:                                       
// ,@`.``.,` .```` ` ,`',#:.,;:::@@@@@@@@@@@@@@@@#+';,:;,..`        ``..;'::;;:;;''  @@.:.:``     .. .```                                       
//   ,.,`.:. ,``..```:.;:':,`:;:;@@@@@@@@@@@@@@@@@@@@. ;::..`        ``.,@;;;;;;;;++ @@ ,`,.`. .`   ,`.      ,;;@@`                             
//   ,.,.`,,`. `..`  , :,;::,:::':@@@@@@@@@@@@@@@@   `.:;;:,.``      ```,:@;;;;;;;'@. @:@ `:,... ;;` `` ,`  .  ::  :,,..,:::,.                  
//   :.`.,;,,. ``.```:`:.+;,:,::@                `   .,:;;;:,.```     ``.,''';;:;;'#@@ @.:``;,.`;.,;,.   `.` ,   ``:,``  `.:+.`.,,:;'++@,,`     
//   ,...,:,:. `.`` `:,,.@.:..::@            :`` ` ``.;;#:;:,..``      ``.,@'';;;''+##@@ .,,.:.',.;..`.. , .::;@,+,,;'':;, ``  ```...  `',``...,
//   .``..,:,, `,.   :,:.@.,:::;@        `,,.````` `..;#@@::,,.```      ``,;@'''''+++@@:,.:,,`,,:::,. ; , `,@@,';;+#@@@@@+:,,:,,:'';;:;';`.:.   
//   .``.,.,,, .,`. ```..@`,,:::@        ;:,....`````..,:#;:,,.```      ``.,@+''+'+'@, ``:#@ :,.;:;,:@ ,`'@,@@@@@@@@@@@@@@@@@':;';::'@@@@: ;; ``
//   .`` :,:,.``. `` `, `@`:,;:;@         ;:,,...````...,:',,,..``       ``.,@++++#: @@:.'@ ;.,,,#,.@    :'@@@@@@@@@@@@+@@@@@@@@@@@@@@@@@@@@@;+@
//   , .`,.,`` `.`` ``.` @`::::;@         @;::,,.......,,:;',,..`         ``.;@@@+;+@@@,,,;;,,,:@';;.                 `;@@@@@@@@@@@@@@@@@@@@@#@@
//`,,` . .,.````,``  `.``@,,,::;@          @;;:::,,,,,,,,::';,..``         `..@@+#@`@@;:.,'.`,.;+''@                                .@@@@@@@@@@@
//+'@`````::`  `.`.```..`@:.,:::+           @++';;;;;::,,,,:',,.``         ``.,@@@##@@@,.,,`;`@`:+;:                                  `         
//@@@````.:,,```..`. ..``+;`.:::',,;'  `:    @@@@@@@'::,,,,,:'..``          ``.'@''+@@'.;+:;`'';#'#                                             
//     `.::,,   `` .``,``+; .`,:;,,.``..,. `,` @@@@@':,,,,,,:;:,.``          ``.;''@@#`;;'''';;';#                                              
//     `,:::, `  .`...,``:;`.`.:,+@@@@@@@@@';.   ,@@;:,,,,,:::',..````       ``.;+@@';;;;;;;;;;;'`                                              
//    `,..,::`.````...,``.@,,.,.,, `:@@@@@@@@@@@@@@@,:::::,:,:;+:,.``        ``.:@@;::::::::;;;'                                                
//    `,, ,,, ` .`.`...``:',,::,,:                `@@;:,,,,,:::;+:,``         ``:@'::::::::;;';,,,,...:,                                        
//    `....,,  `.```..`.;,+````:.:                   @:::,,,,::;+@;,.`        `..;;::::::::;;;:@@@@@';;:,:;,````  `,,',.`                       
//    `..,:,, ``.```.`.,:,:..,,`:'                    +;::,,,,::;'@'.`       ``..:;::::::::::;@@@@@+;..  `...::.  .`.```  ..  ,.`.',..`         
//    `` ,.,:```````..:;...,,.,``@                    :;;,,,,,::;@;         ```.,,::::::;:;+` `:@@@@@@@@@@@@@@@@@#''@@@@'.:.::,,`,,,..  ``.`` `,
//     ` ``,,.,.    `+,: ..;,:` @;                     '::,,,:'@.          ``...,,:,,:;;:;@                    `,'@@@@@@@@@@@@@@@@@@+@@+'';::,,,
//     ```.:,.,`....`..,`..',.,.@;                    .::,,,:;         `....,,,,,,:,,,,::@                                       `'@@@@@@@@@@@@@
//     ,`.;:,..     .``.``.+, :;,:                    ;::::.       ``...,,,::::;;+:,,,,,,                                                       
//    `. :,;   `,++`.;..``.@.:'.:,                    ;,;`      ``....,,::;;;'+@;::,.,,,,`                                                      
//    @`:` ``.` .` ,;`:,...#`,::`:                   @:       ``....,::;'''+@;;;;::,,...,;                                                      
//    `;,`  , `  , `'`,,+,.;..,:.,                  `       ``....::;;''+@':;;;;;::,,...,,'                                                     
//@@@@`..``,.;;+:.         ;.,;;.;.                      ``....,:;''++@';;;;;;:::::,,,,,...;                                     `,,:'@@@@@;;:;'
//    `.` `'`.,``@@,,,,.#'..   ;,'.  `.,``;@,         ``...,,:;'+++@+:;;;;;;;::::,,,,.........                       .`,;';:;:,,:::,;::;,::,,,::
//    ..`  :` . .@.@@@@@@        `@@@@@@;        `````..,:;''++#@+::;;;;;;;;::::,,,....```....,`             .:++''::@;,:,,,,,,,:.,:::;,:;';;;''
//+:;:```                         `            ````.,,:;'++#@@:,:;;;;;;';;;;:::,,...````````````, `.,;@+;,:,:;,,   ::@+'''::::::::::;;::;;'++@@@
//;@@'`` .,,,,,,,,,,:;',.                  `  ``..:;'+#@@;..::::;;'''''';;;:::.``````` ` ````````` @#@';;'';;:.:  ``:@@'@@@@@@@@@@@@@@@@@@.     
//@@@@.`.@@@@;:.                    `````````.,:'+@@@,`.::::;;;''+'''';::,..````````````     `````` ,@'';###@@@@@@`.:.:+@@@@                    
//    . `,,.,:::,..,,:';.:,````````..`.``.,:;++;::::,::::;;'+++++';:,..``````````````` `       `````` @@@@@@;.  ,.``..;+;@@@@                   
//    :,.                  ....,,.....:@@@@@'@@@@::::;;;'+++#':,.....```````````````.```        ``````.         ,..`@;. `;#@                    
//    ,`;@@@'':;':;;'@@@@:` `...,,.'@     :'+@@@;,:;;'+++':,,.....```       ` ``````.,````      ```````.:        @.@@@#'@;;@                    
//    .`,: :,,          .;+'#+#@@@@      @;:;@:.,:;'+#':,,,,....````          ` ```..,.```       ` ```...,       ;,:,,;@@@ @                    
//    .`., ,:@@@@@@@@@@@##+##@#@++@       @;@@@@,:+#'::,,,....`````           ````...::..``      ` ```....       .,:` `.;'@@                    
//    .`., `,:::``;..`,.`,.,.::';'@      `@;@;@;;'+'::,,,....``````         ` ```...,:#,..`   `` ````....,,       :.`,'`:';@                    
//    .`.:,`::,:``:.``,, :...;,':;@       @';;@;;';;::,,....````````   ``````````..,,:@:..``````````....,,,`      :,..`.,;:+                    
//    ```;` ;:,,`,::....`.,,.;,';;@        ;::';'';::,,,.....`````````````````....,,:;@:...```````.....,,,,:      ,,.,`.,;:'                    
//    ```,,.;::,,:',..,,.`,.,,.':;@.       ;::#:';;::,,,.......``````````````...,,,::@@:,............,,,,,,:      ,...`..:;;                    
//    ```,` ':;,:::,;.;.,,:,:..;:;',       ,::+:'';:,,,..........`````````....,,,:::;@':,,,.........,,,,,,,:.     .,...`,::;                    
//     ` ;` '::,,'::;,,.:::,:..::;:;       +,:';';::,,,,,...........``.......,,:::;;@@';::,,,,.,,.,,,,,,,,,::     ..,,`.,:,;`                   
//     . .` :,`,,':,'::.,::,;.`::,:#       +,:';';::,,,,,,.................,,,::;;'@@@';;:::,,,,,,,,,,,,::::;     ....`..,;:`                   
//     .``. ;  :,;::'::.,;::;.`::,.@       ,,:''';:::,,,,,.,...........,,,,:::;''+@@@@@+'';;:::,,,::::::::::;     ..`@ ..,:::                   
//     .`.. :. ,,;:,'::`:;:::;.:,,,@       ::::;';;:::,,,,,,,,,,,,,,,,,::::;''+@@@@@@@@@@@+'';;;;;:::;:;;;;:;     .,....,.:,;                   
//     ,`.` ,``.,,,:+::,;::,;:.,,``+       ,:,:;';;;;:::,,,,,,,,,,,,:::;'+#@@@@@:@@@@@@@@@@@@@@@++''';'';;;:;     `,..`,,,:,;                   
//     ,``  ; ` .,::'::::,,:'  `;';#       ,:,,;';;;;:::,,,,,,,,,,,,::;;;;;;;;@ @@@@@@':,:::''+#@@+'';;;;;;;'      ,,,.,,.:,+                   
//     .``` '`` ::.:':,:;,::.`,::;;@       ::,:.@'';;;::,,,,,,,,,,,,,:::::::;;'+@@@@@;:,,,,,....,,:::::;;;;;;      :..,,.,,,+                   
//     ,.`  ;.` `:,:'::;;::;.`.:;;:@       ::,,`@'';;::::,,,,,,,,,,,,,::::;:;;;#@@@@;::,,,......,,,,::::;;;;:      ,.` `..,:;                   
//     ,...`+.`  `.:;,,:::::'.`:;;:@       ::,,,:+';;;::,,,,,,,,,,,,,,,:::::::;+@@@+;:,,,.........,,::::;;;;+.,,,,,,, @ ,,.::,:,.` :            
//     ,,..`',...` .:,:;:,::':`,;:,+       .,..:`@+'';:::,,,,,,,,,,,,,,,::::;;;#@@@';::,..........,,:::;;;;+       ,,:`,,..,:`     :            
//``   :,.,`+..` ``,,,:;:.;,;;`.::,'`      ,.+;:.#''';;::,,,,,,..,.,,,,,:::;;;;@@@@';::,.........,,,::::;;;;@;;@@@@,.,.,,,.:,;@+;.`+`:,:,,,``.``
//`,   ,...`'`,``. ..,,;..;,;:,..;;+`      .,;....#++';::,,,,,.......,,,::::;;#.@@@';::,,........,,,,:::::;;:,,,,,...,,,,,,,,,  ` `      ` `   `
// :  .,... '.. `` `.``:`.',;;,..+;';;';;,;:,,,..`@+'';;:,,,,.........,,,:::;;+`@@@+;::,,.........,,,,::::;,,..````.,.,,.,,,:,.,...``.:.,.,,,`,`
//  ,` :``. @``  ` ````;..:,;:,..';;;,... ..,,,.``:#+'';;:,,,.........,,:::;;;@.@@@';::,,.........,,,,::::;````.`.`:.,.,,,,,,:.`... `,..,`.,,...
//...`.,.., @`` ` ``. `:`.,,;;:`,+;;::...``.:,:,..:'++'';::,,,........,,::::;;@ ;@@';::,,.........,,:::::::,,,,:,.,.,.,.,` ,:,.:,,,:,.,::,...``.
// `,`.,`.. @``   ` .` , `..;:;;`,',:.`..`,.:;....:.@+'';::,,,,....,,,,,:::;;':,#@@+;;:,,........,,,:::::;,  `.,,,,,,.,...@`,:';@;;,;;:,::;::,:'
//.`..`,:`` @ `     .``, `.`,`:;:...:`.,.,..,:,@...,:#+'';::,,,,...,,,,,:::;;@..:@@+;::,,.......,,,::::::;.,.`..`.`,,,.,,;..,:.`,,,..`,,..`.,,..
//..,`.,:.. @ `:;.:;;. .,;,.,`.:;,..:``..,...::,.`.: @+''';::,,,..,.,,,,::;;'@` :@@+'::,,,,.....,,,,:::::'. ..,....,.:.:,,`,,;:,;,::,:;,.,,,...,
//..`.,,,;;:@,:`,;.:;:;,;:::+,.,:;.:,..,`.`.,:::.`,;.,@+'';::,,,,,,,,,,:::;;'@``:@@+'::,,,......,,,::::;:'.,,,,`..,,,..:,,.:,;'';:;:'';,;';,:,::
//.,,,,`````@````.` ,.;::,:,''`..;'+'`.,`,.,,;:::.,'.`'#+';;:,,,,,,,,,,,:::;'@`,,@@';::,,.......,,,::::;;:,,,....,,,,,:.:,.,,;,.,,::.:`,,,..,:..
//,.,,...```;,.``.``..;::```;,'.`.';',`,.,..,:::,..+`. @+'';;:,,,,,,,,,,,::;;''.`@@';::,........,,,:::;:;,.,.,`,`,::,,,,,:.,,:..`..@,@.,.`.';,;,
//;,:.::,.:,:;`` `  .`;:...:;..';`'.,::.,::,,:,,,..:.;,,@'';;:,,,,,,,,,,,:::;;@':@#':,,,....`....,,:::::'`...,.`..`,,,,:,,;,:;:`,``,,`,;;'::,:',
//`,.'`,.:,::',,,,,;  :..:,,,,;.',,,,..``. ..:,...,'.,..,@';;;:::,,,,,,,,,:::;#;:@@';,,......`...,,,::::,           :,,::,,,:;,   .. .`` `. ...`
//...``, ....+,.,,::,`.`,,`:;;.,:,':#::@,   .,:,.:.;     `@';;:,,,,,,.,,,,,:::;, @@+::,...````...,,,:::'            :,,,:,,::,,                 
//.,,...: ``:#` .``.;,,.,;;.,`.;,..:,' `,;',:.:.;:,:      ,+;;::,,,,....,,,::::@ ,@#;:,....```...,,::,,+            ;.,,,,.,:::                 
//...`;,; `:.+;:;::..:,:.,.,+;+..`:.'@:.    ..,,``.:       @;;::,,......,,,::::@  @@';,,..`````..,,,,,:@            ;,,,,,,,:,:                 
//::.,.::`.; @.`  .'.,..:`;:+,::' ,#,#         ...:'       @;;:::,.......,,,,,:+  ;@+;:,...```...,,,,,,+            :.,.,,.:,,: `    ``         
//      .,,``#,`..,:' ``, '.,,;,`:;,;':;.::,               ;;;::,...```...,,,,:'`  @@';:,..```....,,,..:            ;,,....:,,:           ````  
//,,.   ,`: `@,``.;,`` '`. `.;';`,:'.; ;.`.`. .`;  `. ..   `';;:,..``  `..,,,,:'    @#;:,.....`..........          ,,.....,,;:,         ````    
//````  ;, ` @:;`: ;`   . ;..;'`'...#,    `  ``    ``       @;;:,.`    `.,,,,::+. `.@@+;:,,........,..``,              ``.,:':'  `              
//      ;:,.`@.:..`; @@@@.`;`@,.'.';@:     `          . @   ';;:,..` ``.,,,,,:;@, ` @@#';:,,,,.....,.```.:     ,::,`. `                         
//      ;`.;.@:@, : :+@+:.@`.@,::,':',  `  .             ` ,::;:,..`    `...,,;;   @ .@#';::,,,..,..``   `.   :.`:`;@@`..,.,,,,.,;.,..,. .'`,,;.
//  ```.::.+:+.@ `;'':.:.@;@,, @.,#,:,.,,,:..` .``.    .    ;:::,,..`   ``...,'`   `.`.@++;::,.....       .. `, `` ` `, ': .  .,`;,  `,`;   ```,
//      :,:,:;.@ ..`.` `...`.##:;,+:;,`           ,:    `` .+,:,,..``   ``...,;, .``...`@+';:,...```      `.`@`   .   .      ::.@`.           ` 
//      ,.: ;'.@:.` ,: ,,.++`..,'.@,,.,.```.`,.`.`:.        ;,,....````  `...,:.`. . .`,.@';:,,.```       `.,   , ,`              `   ` :      `
//@@@@# :`;.`;;@:::::`,,,;;,:,:;,,+:,,,`..``....`:..```.   ::,,..`````  ```...,# ``.,.:;.;@';:,.``        `..,` @ `              .   `; ;`      
//````  :..:..:@,:;,``:,`'::;:;@`''.,,+`` ````,`` .`..@@@@@::,,..`````` ```..,,:@@@,      :+;:,.``        ``..,., , `` `,        `     ` `:     
//`   ` :,`:;`:@ ````,::., :,;;:.'@.,.@````   `.,;,,:,;'''';:,,..`````` ````..,:;`,...  ```@;:,.``        ```...,`#@`.:;@ ```@  .@ .    ``      
//      .,`,;,:@:,:,,:,:.:.,`...:':,:.@`..``. ```.`.#@@@@@@;:,,...``````````...,' `  ````  ';:,..```   ``````..:                 ``.`` .```';   
//,,,',`.:....:@,:,,::,:.: .`,,,+::.:'@            ``,     .;:,,...`....``.`...,:,...`...``:;:,..`` ` `  ```.,,,,.`           `                 
//     `.:``.`,@,,,,:,::`.`,.:,:;.,`';@    .```.```.:.....``@;:,,.........`....,,'`......``,;:,..```` ` ````..,,:.::::::.,.....``.,`..`.``````  
//,`.;`..;`..`,;:...:.,:`.,:,:;;,.:.@;@..   `   `      .``.,@';:,,,,,,,,......,,,;:,.......,':,,.````` ````..,::: `  ```  ` ` .`` `````` . ``.` 
//..    `'``.`,':```+ ,`  '::+;::,,,+;+`      `  .``       .`@';::::,,,,.,....,,::,`,`..``..+:,,...````````.,::::`` `  `            ` .`..`.:`  
// .`  ``:;`.`.;:``.,`, .`+;:;:::.:,';' `..` ` ``..,`.` ```...';::,,,,,,,,,.,,,,:;;...``.   @;,,.....``..`.,::::,`.``. ` `.,:,`..`,,;'@@;:,@':,`
//       .;;: ,';`.`. . .,;::;;:.`.,':'`  ``;@@@@:,,.,.``.....@'::,,....,,,,,,,,:;,    `  ` :':,,.........,::::;`` `.,;;:.,`                    
//`     ``;';:`;:`.`,`: ::::;;:;.,.,';;:;@@@,.` ,:,,`,...`` `` @;:,,.....,,,:::::':, `,@@@@@@';:,,,.....,::::::#                                
//`      ,:,,:.::` `. ;,,;,,:;:;..,:';;.        ``:.`.:.`.,`.. :+::,......,,,,:;;+`          @'';::,,,,,:::,,::@                         `   `  
//     `.,:;,,',;`..`.',:,``.::,...;'';' @@@@@@':+'.  `         @';,,......,,,:::'`           +;;;:;:::::,:,,::' :`#.,, ` ` `:.                 
//@@@'@::..,``;,'.::;:':,,`@@;::,,,;:;:;.:``                     @;:,,......,,,,:+  .   ` `:;,@;::::::,,,,,,,:'                                 
//,,,.. `.`..`;.::::::',.,.`.,:,,.:;::;'     ,@+@          `     :+;:,,.....,,,,:@  ` ` `     :;:,,,....,,,,,:#  `                              
//.:::,:;.`.``,.;,:;,:':..,...,...:',:;' ,    ;: .,@@@@@          @';:,,...,.,,:'.       ``    @::,,....,,,,:;:            .``'@@':@`       .::#
//       , .` ,.::,:,:;:.`..`.,..,:;,.+:   ` @`  .',         `.    @';:,,.,.,,,:@ .;  , .``  . ,;:,,.....,,,:;@@@@';,+@@@;:':@@'.:;;++;@@@:,@   
//:   .@ ` .` ,+:.::,:;,,`,.`,,.`,:;``;:        .@         .       :+'::,,,.,,,;@     ,+. ;:@#+:+;:,,....,,:;;':+:@@:@@,.,@'@@@;                
//` @   `.``` `'.:,.,:':` ,` ,,`.::'``:,    :    ',`:@       ``     @;;:,,,,.,,;#   `@.#@+:#:'+@;;:,,....,::;@@+                `.` `     .``@+`
//:@@,:  `..  `# :,:::;`` ,``,.`,:::``:.  ,  `@ +,@::#@@,       @:@@,';:,,,,.,,;:+;@;@@'@;@@@@@@'':,,,,..,::'      ``   ,     `.  ` `` .  .   ` 
// @:@ ,````  `+ .`,,:'`` :``,``,:::..',@  @@@@@#@@.@'@.@;@@,@@+,@,,'@;:,,,,,.,::`  `           +;::,....,:;: ``          `` . `    `` ;, `     
//@'@@@@@.``  .' .....,```,``,..:;;;.:+,;@@` @@@@@@@@@,.              @;,,....,::'      `` ..,..';;:,,...,:'.` ,.`          `.  `. :  .,  `     
//`@@#@@@,``` .+ .```., ` ,``:`..;::,:;,:            . ,       ` ````..':,,..`.,:+`.`,`` .,`  ``,';::,,..,,' `. ,`      `` ,   #```  `   ``  `  
//;`@@.@@:    .# `````.`.`,``;`.`;:;,,',:    ``.`;.,`,``.,`.:,,``   .  @;:,.```,:'...,.;```.:`  .';;:,,.`,,# ``,` ,..   ,`         ``   `@@@@@@@
//    `  . ```.' .`.``.:.`:`.:`..;;;.::,'    `    .` .     . .: .`,::..,'::,.`..,'` ,. `';.. `;@:;;::,..`,:;``: ``    ,@@@@':@.````.`..;;,`     
// ,@  @ . ````' :....,..`,.,,`.:::;:';:;    `` ` ``  ':@@'@+::;@@;;    +;:,,,` .:;  `  `   .`;+;;;;::,.``:.   ` , ` `    . `,. . ` `    ``  `  
//      @:`,```',,,,...,,.,`,,`.;;:;;;+,'    :.    :    .;#  .;'.,;'@# ,@:;:.,,  :'@@:'`@@@@@+#@;;;;::,.``;' .`.`:@'., `,:,@@@'` ```  `.   . :``
//@++@@@@@`.` .:.``...`,,`,.,.,;,';::::,,@ . `,', :.;`  .`    , ,,  ';@@+,:'..    ,   `:@@'@:;+@:;::;;,.``':..'+@@,:;;,.                 `  `  `
//';@;@@''`.  `:. ....`..`:`,,`,,;;;;;::: . .````  .'..,.@@';@'@@@@@@@@+@:,,:;      @@@+@#@@@@'+,::,:;:.  ,, ``. `  .`  `` ` `         `` ` ` , 
//+#'''#;;`.`  ;` .``.....,.,...,';;':;',,`   ``.`, `: ``:.`:@@@+@;@@+`:,;,,,         ,`.  `,  ,:,`.,::.`  '      `   .` ` .``.,..   ```  `    .
//':,@@@@@ ,`  ;` .`,````.` ,...,;:::;''`.`     ` .@;.` .  `;            @:,``           `  `.,.',`.,;@     `.:``..,`  ``  ,   ` ``` . ``     ` 
//    '@;'`.`  '.`.......:`..`..,,:;;;+#.:       `.`,`:::@@@. `      ``` `',``          `````` :@:,.;        `   ..` `., ,`, `.`   `. ., ` , :,.
//        `.`  ', `.`...`;...`..,.;:':+@`,.  `. ` ` . `,`.,    ,'@@@@.,   `..``   `   `,  ```  @;:,, `        `.`` .   `  .            ```.```:`
//@@    :@ ``  ;,``..:.``;``..`...,:::;;';@.,.` `.,`,``, . ,``` `    ```:@@,.```        ,`,`.`.:;;`` `         .`` ` .     ` .     ,  ` .`.  .  
//.`.. `. .`` `;:`.:,:.``;``.`...`,.:,.,:.  `..,,,,``.`..,..```. ``,.`.,:@@:.....      `,.`   ':;```.          ..        ` `        `  `.  ``,  
//. ````. :  `.:;,.,.,```;``.`,.,`,`,;;.`,`....;.:;`,....;...,...,``.. ` ` `,,.`.`     `,..,;;+#``..`     ``, ,.  `` `,.``` . .`` . .`,`` ` ` .,
//  ` `.:.; ``,:;,,..;```; ```.@''+',..':; :., ;`: ..,. .,.,`` .,`.:`,.:``,,,,`@@`.    @ ,`` :;,`....`'. ` ```:... `` ,  `.` ```..`.. ,`....`   
//   ,.`.`;`.`.,',:,`:`. . ``;:`` ,.;;..`:.,`.,.. `.:.,`.,`...., , ;,,````, :, @@@:`   @..,``+.,,.,,+.,`.`  `,` ` ``..,. ,`.` .  .````  ` `     
//        :`,.:,',,.,',';#,:';.`.:;. +..`:` ,.:  .,,...`. `:;,.,  ,`..:`.`, :. @@@@`.  @ `,. ',;`.,.@@`...` . `.```.,...``,,. . ```  `. ,`,  .  
//   ``   ;. `` '.`````.,',:;: ::..;: ,;`@. .` `: :+ ;,,,:,.,`,,, .`. ,: ,  .,.;@@@@   @ `.,;::,,:;@@@@..`  ` `,`, .`.````; `.`` `. .    ``     
//```,` ..;.`,.`::`.``.``,,;:, ;``@. ,.`;; .`..'::` ;.,`.:````  ```. `.` .`,`,,,@@@@:` @ ,`.,:;,,,`@@@@@ `  `.,`.`..``.` ``  ... `   ` ``.``   `
//`  . . `.' ```,;,,::,.` ..,`., ::``., ;' ,`': ,. `  `  ` ` `  `` `    `  ` `, @#@@@   .  . ;:::,,@@@@@ . .```.``  ` ```` ````` ``` `..`   ``  
//. `  ` `.;`.  ,@.',.,,,.  ` ,, ;:``.,,,;.` ..``.  `````. `...``,``..`.  `. .. @@+@@@  @ ` .,,.,`@@@@@@   ..`.. ``,. ``,`;`   `` .` , `` ` `  `
// ``  ..` :`:` `@.':',,,:,: `., ;;,`,`.``,......``. ```,``  `  `,``...`  ` `,. @##@@@` @@  `:',`,@@@@@@.  ,````.`.,,..` `..   ` ``  ``  ` `  ,`
// . ```, ` ,,` ,@,,`.,   .:,::`,:,;. ,..:::..`.`..` ``  . ```    .``` ,,: ``,, @''+#@@  @`  :;,`@@@@@@@,  ,...,.,`    ` ` ` ..  . ` `` .` `  ``
//  ;: `. `::.,.+ ;,.. ` `` :  ,;.`:+,,,.,.,.,::,:,.`  ```   ;` . `` .., ,...`, @'#+#+@  @@  :+.@@@@@@@@;  :..`.`` `  .``  .::   , ``  ,.```,`.`
//`.`...,.,#,,, ,@,,   ',`.`@`@,`.  .``.,`,,,....,,.``.,,,`  .    .``...`.`:`., @#@@@#@@@@@@ :;@@@@@@@@@@  ,,,.. ;. ``    `,, :.`. .  :  `  ` . 
//....` ...`,`.`,. ;#: ,  ,.```,`   `.``.`:@;,...``..`.:  `.,   `     .```, `:,:@'@@@@@@@@@@@@:;.  .. .,,  .. ,.``  .`. `,  `` ``  ` .  . .  ,.,
//   ::,::'`` ..':. ,,    ,,,   .`` ,,,@#;       ` ` `..,,;#,.:.;    .', ..`,`:,+@@@@@@@@@@@#:,':..,`,. @@@  .      `` ``    ` .     .`` `` `   
//  ` `   ,+::.,;@+'#,..` :,...,,                 ``.........  `        `,+,. .  .   ` ,,.,`,   .  ,;.  ` .  .   ,   .   .@.  @ .:  .`. .   ,`. 
//
