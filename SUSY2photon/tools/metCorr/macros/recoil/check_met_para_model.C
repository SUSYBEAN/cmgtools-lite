{

bool _isDyJetsLO = false;

std::string _RecoilInputFileNameData_all, _RecoilInputFileNameData_mu, _RecoilInputFileNameData_el;
std::string _RecoilInputFileNameMC_all, _RecoilInputFileNameMC_mu, _RecoilInputFileNameMC_el;
std::string _RecoilInputFileNameMCLO_all, _RecoilInputFileNameMCLO_mu, _RecoilInputFileNameMCLO_el;

std::string _dir = "./recoil_out/";

std::string _dt_tag = "met_para_study";
//std::string _mc_tag = "RecoilNoSmooth_met_para_study";
//std::string _dt_tag = "met_para_study_fullCuts";
//std::string _mc_tag = "RecoilNoSmooth_met_para_study_fullCuts_effSf";
//std::string _dt_tag = "met_para_study";
//std::string _mc_tag = "RecoilNoPUWtNoSmooth_met_para_study";
//std::string _dt_tag = "met_para_study_fullCuts";
//std::string _mc_tag = "RecoilNoPUWtNoSmooth_met_para_study_fullCuts_effSf";
std::string _mc_tag = "met_para_study_effSf";

_RecoilInputFileNameData_all = _dir+"SingleEMU_Run2016BCD_PromptReco_"+_dt_tag+".root";
_RecoilInputFileNameData_mu = _dir+"SingleEMU_Run2016BCD_PromptReco_"+_dt_tag+"_mu.root";
_RecoilInputFileNameData_el = _dir+"SingleEMU_Run2016BCD_PromptReco_"+_dt_tag+"_el.root";
_RecoilInputFileNameMC_all = _dir+"DYJetsToLL_M50_"+_mc_tag+".root";
_RecoilInputFileNameMC_mu = _dir+"DYJetsToLL_M50_"+_mc_tag+"_mu.root";
_RecoilInputFileNameMC_el = _dir+"DYJetsToLL_M50_"+_mc_tag+"_el.root";
_RecoilInputFileNameMCLO_all = _dir+"DYJetsToLL_M50_MGMLM_Ext1_"+_mc_tag+".root";
_RecoilInputFileNameMCLO_mu = _dir+"DYJetsToLL_M50_MGMLM_Ext1_"+_mc_tag+"_mu.root";
_RecoilInputFileNameMCLO_el = _dir+"DYJetsToLL_M50_MGMLM_Ext1_"+_mc_tag+"_el.root";


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

