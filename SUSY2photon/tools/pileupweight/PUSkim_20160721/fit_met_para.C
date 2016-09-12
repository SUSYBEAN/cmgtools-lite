

std::string inputdir = 
  //"/datag/heli/XZZ/80X_20160721_LinksForSkim"
  "./"
  ;
std::string filename =
  //"DYJetsToLL_M50"
  "SingleEMU_Run2016BCD_PromptReco"
  ;

std::string selec = 
  "(1)"
  //"abs(llnunu_l1_l1_pdgId)==11&&abs(llnunu_l1_l2_pdgId)==11"
  //"abs(llnunu_l1_l1_pdgId)==13&&abs(llnunu_l1_l2_pdgId)==13"
  ;
std::string tag = 
  "_met_para_study"
  //"_met_para_study_el"
  //"_met_para_study_mu"
;

std::string lumiTag = 
  "CMS 13 TeV 2016 L=12.9 fb^{-1}"
  //"CMS 13 TeV Simulation for 2016 Data"
  ;

char name[1000];
std::string histname;
TFile* fin;
TFile* fout;
TPaveText* lumipt;
TCanvas* plots;
TH2D* h2d1;
TH2D* h2d2;
TH2D* h2d3;
TH2D* h2d4;
TH1D* h1d1[1000];
TH1D* h1d2[1000];
TH1D* h1d3[1000];
TH1D* h1d4[1000];
TF1* func1[1000];
TF1* func2[1000];
TF1* func3[1000];
TF1* func4[1000];

Int_t Nbins;
std::vector<Double_t> fit_min; 
std::vector<Double_t> fit_max; 
std::vector<Int_t> fit_rebin; 

int fit_slice_gaus(TH2D* h2d, TH1D** h1d);

void fit_met_para(){

  gROOT->ProcessLine(".x tdrstyle.C");
  gStyle->SetOptTitle(0);

  lumipt = new TPaveText(0.2,0.9,0.8,0.98,"brNDC");
  lumipt->SetBorderSize(0);
  lumipt->SetTextAlign(12);
  lumipt->SetFillStyle(0);
  lumipt->SetTextFont(42);
  lumipt->SetTextSize(0.03);
  lumipt->AddText(0.15,0.3, lumiTag.c_str());

  sprintf(name, "%s/%s.root", inputdir.c_str(), filename.c_str());
  fin = new TFile(name);


  sprintf(name, "%s%s.root", filename.c_str(), tag.c_str());
  fout = new TFile(name, "recreate");

  TTree* tree = (TTree*)fin->Get("tree");

  plots = new TCanvas("plots", "plots");

  sprintf(name, "%s%s.ps[", filename.c_str(), tag.c_str());
  plots->Print(name);


  // other control plots
  Double_t ZPtBins[] = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28, 30, 35, 40, 50, 60, 80, 100, 150, 250, 5000 };
  Int_t NZPtBins = sizeof(ZPtBins)/sizeof(ZPtBins[0]) - 1;
  const Int_t NMetParaBins=500;
  Double_t MetParaBins[NMetParaBins+1];
  for (int i=0; i<=NMetParaBins; i++) { MetParaBins[i] = -100.0+200.0/NMetParaBins*i; };
  const Int_t NMetPerpBins=500;
  Double_t MetPerpBins[NMetPerpBins+1];
  for (int i=0; i<=NMetPerpBins; i++) { MetPerpBins[i] = -100.0+200.0/NMetPerpBins*i; };

  h2d1 = new TH2D("h_met_para_vs_zpt", "h_met_para_vs_zpt", NZPtBins, ZPtBins, NMetParaBins, MetParaBins);
  h2d2 = new TH2D("h_met_perp_vs_zpt", "h_met_perp_vs_zpt", NZPtBins, ZPtBins, NMetPerpBins, MetPerpBins);
  h2d1->Sumw2();
  h2d2->Sumw2();
  tree->Draw("llnunu_l2_pt*cos(llnunu_l2_phi-llnunu_l1_phi):llnunu_l1_pt>>h_met_para_vs_zpt", selec.c_str(), "colz");
  tree->Draw("llnunu_l2_pt*sin(llnunu_l2_phi-llnunu_l1_phi):llnunu_l1_pt>>h_met_perp_vs_zpt", selec.c_str(), "colz");
  
  h2d1->GetXaxis()->SetTitle("P_{T}(Z) (GeV)");
  h2d1->GetYaxis()->SetTitle("MET para (GeV)");
  h2d2->GetXaxis()->SetTitle("P_{T}(Z) (GeV)");
  h2d2->GetYaxis()->SetTitle("MET para (GeV)");

  Nbins = NZPtBins;
  fit_min = { -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -25  };
  fit_max = { +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20  };
  fit_rebin = {  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4  };
  for (int ii=0; ii<Nbins; ii++) fit_rebin[ii] = 5;
  for (int ii=0; ii<Nbins; ii++) fit_min[ii] = -16;
  for (int ii=0; ii<Nbins; ii++) fit_max[ii] = 16;

  fit_min = { -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -30, -30, -30, -30  };
  fit_max = { +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +20, +30, +30, +30, +30  };

  fit_slice_gaus(h2d1, h1d1);
  fit_slice_gaus(h2d2, h1d2);



  sprintf(name, "%s%s.ps]", filename.c_str(), tag.c_str());
  plots->Print(name);


  sprintf(name, ".! ps2pdf %s%s.ps %s%s.pdf", filename.c_str(), tag.c_str(), filename.c_str(), tag.c_str());
  gROOT->ProcessLine(name);

  fout->Close();
}

int fit_slice_gaus(TH2D* h2d, TH1D** h1d){ 

  std::string hname = h2d->GetName();
  Double_t* xbins = (Double_t*)h2d->GetXaxis()->GetXbins()->GetArray();

  plots->cd();
  plots->Clear();
  plots->SetLogx(1);
  h2d->Draw("colz");
  lumipt->Draw();
  sprintf(name, "%s%s.ps", filename.c_str(), tag.c_str());
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();
  
  fout->cd();
  h2d->Write();

  sprintf(name, "%s_mean", hname.c_str());
  TH1D* h_mean = new TH1D(name, name, Nbins, xbins);
  h_mean->Sumw2();
  h_mean->GetXaxis()->SetTitle("P_{T}(Z) (GeV)");
  h_mean->GetYaxis()->SetTitle("mean (GeV)");
  h_mean->SetLineColor(2);
  h_mean->SetMarkerColor(2);
  h_mean->SetMarkerStyle(20);
  sprintf(name, "%s_sigma", hname.c_str());
  TH1D* h_sigma = new TH1D(name, name, Nbins, xbins);
  h_sigma->Sumw2();
  h_sigma->GetXaxis()->SetTitle("P_{T}(Z) (GeV)");
  h_sigma->GetYaxis()->SetTitle("sigma (GeV)");
  h_sigma->SetLineColor(4);
  h_sigma->SetMarkerColor(4);
  h_sigma->SetMarkerStyle(20);

  for (int i=0; i<Nbins; i++){
    sprintf(name, "%s_bin%d_func", hname.c_str(), i+1);
    TF1* afunc = new TF1(name, "gaus", -100,+100);
    sprintf(name, "%s_bin%d", hname.c_str(), i+1);
    TH1D* ahist = (TH1D*)h2d->ProjectionY(name, i+1, i+1, "e");
    ahist->SetTitle(name);
    ahist->Rebin(fit_rebin[i]);
    ahist->Fit(afunc, "R", "", fit_min[i], fit_max[i]);
    double mean = afunc->GetParameter(1);
    double sigma = afunc->GetParameter(2);
    ahist->Fit(afunc, "R", "", mean-1.5*sigma, mean+1.5*sigma); 
    func1[i] = afunc;
    h1d1[i] = ahist;

    h_mean->SetBinContent(i+1, afunc->GetParameter(1));
    h_mean->SetBinError(i+1, afunc->GetParError(1));
    h_sigma->SetBinContent(i+1, afunc->GetParameter(2));
    h_sigma->SetBinError(i+1, afunc->GetParError(2));
    

    plots->cd();
    plots->Clear();
    ahist->Draw();
    lumipt->Draw();
    sprintf(name, "%s%s.ps", filename.c_str(), tag.c_str());
    plots->Print(name);    
    plots->Clear();    
   
    fout->cd();
    ahist->Write();
    afunc->Write();
  }


  plots->cd();
  plots->Clear();
  plots->SetLogx(1);
  h_mean->Draw();
  lumipt->Draw();
  sprintf(name, "%s%s.ps", filename.c_str(), tag.c_str());
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();

  plots->cd();
  plots->Clear();
  plots->SetLogx(1);
  h_sigma->Draw();
  lumipt->Draw();
  sprintf(name, "%s%s.ps", filename.c_str(), tag.c_str());
  plots->Print(name);
  plots->SetLogx(0);
  plots->Clear();

  fout->cd();
  h_mean->Write();
  h_sigma->Write();


  return 0;

}
