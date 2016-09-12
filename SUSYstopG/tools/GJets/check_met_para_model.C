{


  TFile* file_dt_sigma[10];
  TFile* file_mc_sigma[10];
  TH1D* h_dt_met_para_shift[10];
  TH1D* h_mc_met_para_shift[10];
  TH1D* h_dt_met_para_sigma[10];
  TH1D* h_dt_met_perp_sigma[10];
  TH1D* h_mc_met_para_sigma[10];
  TH1D* h_mc_met_perp_sigma[10];
  TH1D* h_ratio_met_para_sigma_dtmc[10];
  TH1D* h_ratio_met_perp_sigma_dtmc[10];
  TGraphErrors* gr_dt_met_para_shift[10];
  TGraphErrors* gr_mc_met_para_shift[10];
  TGraphErrors* gr_ratio_met_para_sigma_dtmc[10];
  TGraphErrors* gr_ratio_met_perp_sigma_dtmc[10];


  file_dt_sigma[0] = new TFile("SingleEMU_Run2016BCD_PromptReco_killdup_met_para_study.root");
  file_dt_sigma[1] = new TFile("SingleEMU_Run2016BCD_PromptReco_killdup_met_para_study_mu.root");
  file_dt_sigma[2] = new TFile("SingleEMU_Run2016BCD_PromptReco_killdup_met_para_study_el.root");
  file_mc_sigma[0] = new TFile("DYJetsToLL_M50_met_para_study.root");
  file_mc_sigma[1] = new TFile("DYJetsToLL_M50_met_para_study_mu.root");
  file_mc_sigma[2] = new TFile("DYJetsToLL_M50_met_para_study_el.root");

  h_dt_met_para_shift[0] = (TH1D*)file_dt_sigma[0]->Get("h_met_para_vs_zpt_mean");
  h_mc_met_para_shift[0] = (TH1D*)file_mc_sigma[0]->Get("h_met_para_vs_zpt_mean");

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
  gr_mc_met_para_shift[3] = new TGraphErrors(h_dt_met_para_shift[3]);
  gr_dt_met_para_shift[4] = new TGraphErrors(h_dt_met_para_shift[4]);
  gr_mc_met_para_shift[4] = new TGraphErrors(h_dt_met_para_shift[4]);
  gr_dt_met_para_shift[5] = new TGraphErrors(h_dt_met_para_shift[5]);
  gr_mc_met_para_shift[5] = new TGraphErrors(h_dt_met_para_shift[5]);

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
