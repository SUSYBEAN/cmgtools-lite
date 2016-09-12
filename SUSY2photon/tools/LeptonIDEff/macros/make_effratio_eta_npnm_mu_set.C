void make_effratio_eta_npnm_mu_set(TString dir="")
{
  //gROOT->ProcessLine(".x tdrstyle.C");
  TFile* fdt = new TFile("rootfiles/"+dir+"/fiteff_npnm_data_mu_set.root");
  TFile* fmc = new TFile("rootfiles/"+dir+"/fiteff_npnm_mc_mu_set.root");
  TFile* fmcgen = new TFile("rootfiles/effhists_eta_npnm_fullmc_mu_set.root");
  TFile* fout = new TFile("rootfiles/"+dir+"/effratio_npnm_mu_set.root", "recreate");
  
  TString plotout("plots/"+dir+"/effratio_npnm_mu_set.pdf");

  TH1D* hmc_np = (TH1D*)fmc->Get("hnp");
  TH1D* hmc_nm = (TH1D*)fmc->Get("hnm");

  TH1D* hmcgen_np = (TH1D*)fmcgen->Get("hnp");
  TH1D* hmcgen_nm = (TH1D*)fmcgen->Get("hnm");

  TH1D* hdt_np = (TH1D*)fdt->Get("hnp");
  TH1D* hdt_nm = (TH1D*)fdt->Get("hnm");
  
  TH1D* hmc_eff = (TH1D*)hmc_np->Clone("hmc_eff");
  TH1D* hdt_eff = (TH1D*)hdt_np->Clone("hdt_eff");
  TH1D* hmcgen_eff = (TH1D*)hmcgen_np->Clone("hmcgen_eff");

  hmc_eff->Reset();
  hmcgen_eff->Reset();
  hdt_eff->Reset();

  int nbins = hmc_eff->GetNbinsX();

  // calculate eff and eff_err
  for (int i=0; i<=nbins; i++){
    double np, nm, np_err, nm_err, eff, eff_err;
    // mc
    np = hmc_np->GetBinContent(i);
    np_err = hmc_np->GetBinError(i);
    nm = hmc_nm->GetBinContent(i);
    nm_err = hmc_nm->GetBinError(i);
    if (np+nm!=0) eff = np/(np+nm);
    else eff = 0.;
    if (np+nm!=0) eff_err = sqrt( pow(np_err,2.)*pow(nm,2.)/pow((np+nm), 4.) + pow(nm_err,2.)*pow(np,2.)/pow((np+nm), 4.) );
    else eff_err = 0.;
    hmc_eff->SetBinContent(i, eff);
    hmc_eff->SetBinError(i, eff_err);
    // mc gen
    np = hmcgen_np->GetBinContent(i);
    np_err = hmcgen_np->GetBinError(i);
    nm = hmcgen_nm->GetBinContent(i);
    nm_err = hmcgen_nm->GetBinError(i);
    if (np+nm!=0) eff = np/(np+nm);
    else eff = 0.;
    if (np+nm!=0) eff_err = sqrt( pow(np_err,2.)*pow(nm,2.)/pow((np+nm), 4.) + pow(nm_err,2.)*pow(np,2.)/pow((np+nm), 4.) );
    else eff_err = 0.;
    hmcgen_eff->SetBinContent(i, eff);
    hmcgen_eff->SetBinError(i, eff_err);
    // dt
    np = hdt_np->GetBinContent(i);
    np_err = hdt_np->GetBinError(i);
    nm = hdt_nm->GetBinContent(i);
    nm_err = hdt_nm->GetBinError(i);
    if (np+nm!=0) eff = np/(np+nm);
    else eff = 0.;
    if (np+nm!=0) eff_err = sqrt( pow(np_err,2.)*pow(nm,2.)/pow((np+nm), 4.) + pow(nm_err,2.)*pow(np,2.)/pow((np+nm), 4.) );
    else eff_err = 0.;
    hdt_eff->SetBinContent(i, eff);
    hdt_eff->SetBinError(i, eff_err);
  }

  TH1D* hratio_dt_mc = (TH1D*)hdt_eff->Clone("hratio_dt_mc");
  hratio_dt_mc->Divide(hmc_eff);
  hratio_dt_mc->SetMarkerColor(1);
  hratio_dt_mc->SetLineColor(1);
  hratio_dt_mc->SetMarkerStyle(20);
  hratio_dt_mc->SetMarkerSize(1);

  hmc_eff->SetMarkerStyle(20);
  hmc_eff->SetLineColor(4);
  hmc_eff->SetMarkerColor(4);
  hmcgen_eff->SetMarkerStyle(20);
  hmcgen_eff->SetLineColor(3);
  hmcgen_eff->SetMarkerColor(3);
  hdt_eff->SetMarkerStyle(20);
  hdt_eff->SetLineColor(2);
  hdt_eff->SetMarkerColor(2);


  TF1* fratio_dt_mc = new TF1("fratio_dt_mc", "pol19", 25, 70);

  TCanvas* c0 = new TCanvas("c0");

  //hratio_dt_mc->Fit(fratio_dt_mc, "R");
  
  hratio_dt_mc->Draw();
  hmc_eff->Draw("same");
  hdt_eff->Draw("same");

  fratio_dt_mc->SetLineColor(kBlue);
  
  fout->cd();
  hratio_dt_mc->Write();
  hmc_eff->Write();
  hmcgen_eff->Write();
  hdt_eff->Write();

  c0->SaveAs(plotout);

}
