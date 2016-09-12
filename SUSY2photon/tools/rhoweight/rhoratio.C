{


  TFile* fin = new TFile("CC80X_L12p9_HLTv2_ichepPU_ZPTwt_allSF_CT_tight_puWeight68075_metfilter_all_1pb.root");
  TFile* fout = new TFile("rhoratio.root", "recreate");

  gROOT->ProcessLine(".x tdrstyle.h");

  TH1D* hdt = (TH1D*)fin->Get("CC80X_L12p9_HLTv2_ichepPU_ZPTwt_allSF_CT_tight_puWeight68075_metfilter_all_1pb_rho_data_SingleEMU_Run2016BCD_PromptReco_killdup");
  TH1D* hmc = (TH1D*)((THStack*)fin->Get("CC80X_L12p9_HLTv2_ichepPU_ZPTwt_allSF_CT_tight_puWeight68075_metfilter_all_1pb_rho_stack"))->GetHistogram();
  
  TH1D* hdt_norm = (TH1D*)hdt->Clone("h_rho_dt_norm");
  TH1D* hmc_norm = (TH1D*)hmc->Clone("h_rho_mc_norm");

  hdt_norm->Scale(1./hdt_norm->Integral());
  hmc_norm->Scale(1./hmc_norm->Integral());

  hratio_dt_mc = (TH1D*)hdt_norm->Clone("h_rho_ratio_dt_mc");
  hratio_dt_mc->Divide(hmc_norm);

  hratio_dt_mc->SetTitle("Rho reweighting: Normalized data/MC, 2016 data 12.9 fb^{-1}.");
  hratio_dt_mc->GetYaxis()->SetTitle("Data/MC");

  char name[10000];
  //sprintf(name, "pol3");
  //sprintf(name, "[0]*exp(-(x-[2])*(x-[2])/[4]/[4])+[1]*exp(-(x-[3])*(x-[3])/[5]/[5])");
  sprintf(name, "gaus(0)+gaus(3)");

  TCanvas* plots = new TCanvas("plots", "plots");

  TF1* fratio_dt_mc = new TF1("f_rho_ratio_dt_mc", name, 0,50);
  fratio_dt_mc->SetParameter(0,0.5);
  fratio_dt_mc->SetParameter(1,8);
  fratio_dt_mc->SetParameter(2,10);
  fratio_dt_mc->SetParameter(3,0.5);
  fratio_dt_mc->SetParameter(4,13);
  fratio_dt_mc->SetParameter(5,20);
  fratio_dt_mc->SetParLimits(0,0,2.0);
  fratio_dt_mc->SetParLimits(1,4,10);
  fratio_dt_mc->SetParLimits(2,5,20);
  fratio_dt_mc->SetParLimits(3,0,2.0);
  fratio_dt_mc->SetParLimits(4,6,40);
  fratio_dt_mc->SetParLimits(5,10,50);

  hratio_dt_mc->Fit(fratio_dt_mc);
  fratio_dt_mc->Print("v");

  TPaveText* lumipt;
  lumipt = new TPaveText(0.2,0.66,0.8,0.8,"brNDC");
  lumipt->SetBorderSize(0);
  lumipt->SetTextAlign(12);
  lumipt->SetFillStyle(0);
  lumipt->SetTextFont(42);
  lumipt->SetTextSize(0.035);
  sprintf(name, "f(#rho) = %.2f #upoint Gaus(#rho; #mu=%.2f, #sigma=%.2f) ", fratio_dt_mc->GetParameter(0), fratio_dt_mc->GetParameter(1), fratio_dt_mc->GetParameter(2));
  lumipt->AddText(0.0,0.6, name);
  sprintf(name, "      + %.2f #upoint Gaus(#rho; #mu=%.2f, #sigma=%.2f) ", fratio_dt_mc->GetParameter(3), fratio_dt_mc->GetParameter(4), fratio_dt_mc->GetParameter(5)); 
  lumipt->AddText(0.0,0.3, name);
  lumipt->Draw();

  plots->Print("rho_rewegiht_12p9.pdf");
  plots->SaveAs("rho_rewegiht_12p9.C");


  sprintf(name, "(%.3f*exp(-0.5*pow((x-%.3f)/%.3f,2))+%.3f*exp(-0.5*pow((x-%.3f)/%.3f,2)))",
          fratio_dt_mc->GetParameter(0),
          fratio_dt_mc->GetParameter(1),
          fratio_dt_mc->GetParameter(2),
          fratio_dt_mc->GetParameter(3),
          fratio_dt_mc->GetParameter(4),
          fratio_dt_mc->GetParameter(5)
         );

  std::cout << "#############################################" << std::endl;          
  std::cout << "func = " << name << std::endl;
  std::cout << "#############################################" << std::endl;          
 
  fout->cd();
  hdt->Write();
  hmc->Write();
  hdt_norm->Write();
  hmc_norm->Write();
  hratio_dt_mc->Write();
  fratio_dt_mc->Write();


//  fout->Close(); 


}
