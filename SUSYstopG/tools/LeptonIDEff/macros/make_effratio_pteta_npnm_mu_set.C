void make_effratio_pteta_npnm_mu_set(TString dir="")
{
  //gROOT->ProcessLine(".x tdrstyle.C");
  TFile* fdt = new TFile("rootfiles/"+dir+"/fiteff_npnm_data_mu_set.root");
  TFile* fmc = new TFile("rootfiles/"+dir+"/fiteff_npnm_mc_mu_set.root");
  TFile* fmcgen = new TFile("rootfiles/effhists_pteta_npnm_fullmc_mu_set.root");
  TFile* fout = new TFile("rootfiles/"+dir+"/effratio_npnm_mu_set.root", "recreate");
  
  TString plotout("plots/"+dir+"/effratio_npnm_mu_set.pdf");

  TH2D* hmc_np = (TH2D*)fmc->Get("hnp");
  TH2D* hmc_nm = (TH2D*)fmc->Get("hnm");

  TH2D* hmcgen_np = (TH2D*)fmcgen->Get("hnp");
  TH2D* hmcgen_nm = (TH2D*)fmcgen->Get("hnm");

  TH2D* hdt_np = (TH2D*)fdt->Get("hnp");
  TH2D* hdt_nm = (TH2D*)fdt->Get("hnm");
  
  TH2D* hmc_eff = (TH2D*)hmc_np->Clone("hmc_eff");
  TH2D* hdt_eff = (TH2D*)hdt_np->Clone("hdt_eff");
  TH2D* hmcgen_eff = (TH2D*)hmcgen_np->Clone("hmcgen_eff");

  hmc_eff->Reset();
  hmcgen_eff->Reset();
  hdt_eff->Reset();

  int nbinsx = hmc_eff->GetNbinsX();
  int nbinsy = hmc_eff->GetNbinsY();
  // calculate eff and eff_err
  for (int i=0; i<=nbinsx; i++){
    for (int j=0; j<=nbinsy; j++){
    double np, nm, np_err, nm_err, eff, eff_err;
    // mc
    np = hmc_np->GetBinContent(i,j);
    np_err = hmc_np->GetBinError(i,j);
    nm = hmc_nm->GetBinContent(i,j);
    nm_err = hmc_nm->GetBinError(i,j);
    if (np+nm!=0) eff = np/(np+nm);
    else eff = 0.;
    if (np+nm!=0) eff_err = sqrt( pow(np_err,2.)*pow(nm,2.)/pow((np+nm), 4.) + pow(nm_err,2.)*pow(np,2.)/pow((np+nm), 4.) );
    else eff_err = 0.;
    hmc_eff->SetBinContent(i, j,eff);
    hmc_eff->SetBinError(i,j, eff_err);
    // mc gen
    np = hmcgen_np->GetBinContent(i,j);
    np_err = hmcgen_np->GetBinError(i,j);
    nm = hmcgen_nm->GetBinContent(i,j);
    nm_err = hmcgen_nm->GetBinError(i,j);
    if (np+nm!=0) eff = np/(np+nm);
    else eff = 0.;
    if (np+nm!=0) eff_err = sqrt( pow(np_err,2.)*pow(nm,2.)/pow((np+nm), 4.) + pow(nm_err,2.)*pow(np,2.)/pow((np+nm), 4.) );
    else eff_err = 0.;
    hmcgen_eff->SetBinContent(i, j,eff);
    hmcgen_eff->SetBinError(i,j, eff_err);
    // dt
    np = hdt_np->GetBinContent(i,j);
    np_err = hdt_np->GetBinError(i,j);
    nm = hdt_nm->GetBinContent(i,j);
    nm_err = hdt_nm->GetBinError(i,j);
    if (np+nm!=0) eff = np/(np+nm);
    else eff = 0.;
    if (np+nm!=0) eff_err = sqrt( pow(np_err,2.)*pow(nm,2.)/pow((np+nm), 4.) + pow(nm_err,2.)*pow(np,2.)/pow((np+nm), 4.) );
    else eff_err = 0.;
    hdt_eff->SetBinContent(i, j,eff);
    hdt_eff->SetBinError(i, j,eff_err);
    }
  }
  TH2D* hratio_dt_mc = (TH2D*)hdt_eff->Clone("hratio_dt_mc");
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

  
  fout->cd();
  hratio_dt_mc->Write();
  hmc_eff->Write();
  hmcgen_eff->Write();
  hdt_eff->Write();



}
