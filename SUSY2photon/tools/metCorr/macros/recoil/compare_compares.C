{

TFile *file0 = TFile::Open("compare_dt_sigma_MzCut_vs_nothing.root");
TFile *file1 = TFile::Open("compare_dt_sigma_ZSelecLowLPt_vs_nothing.root");
TFile *file2 = TFile::Open("compare_dt_sigma_ZSelec_vs_nothing.root");
TFile *file3 = TFile::Open("compare_dt_sigma_ZSelecLowLPt_dtHLT_vs_nothing.root");

gROOT->ProcessLine(".x tdrstyle.C");

TCanvas* c0 = (TCanvas*)file0->Get("c1");
TCanvas* c1 = (TCanvas*)file1->Get("c1");
TCanvas* c2 = (TCanvas*)file2->Get("c1");
TCanvas* c3 = (TCanvas*)file3->Get("c1");


c0->SetName("c0");
c1->SetName("c1");
c2->SetName("c2");
c3->SetName("c3");


TH1D* h0 = (TH1D*)c0->FindObject("h_ratio_met_para_sigma_dtmc_all");
TH1D* h1 = (TH1D*)c1->FindObject("h_ratio_met_para_sigma_dtmc_all");
TH1D* h2 = (TH1D*)c2->FindObject("h_ratio_met_para_sigma_dtmc_all");
TH1D* h3 = (TH1D*)c3->FindObject("h_ratio_met_para_sigma_dtmc_all");


h0->SetName("compare_dt_sigma_MzCut_vs_nothing");
h1->SetName("compare_dt_sigma_ZSelecLowLPt_vs_nothing");
h2->SetName("compare_dt_sigma_ZSelec_vs_nothing");
h3->SetName("compare_dt_sigma_ZSelecLowLPt_dtHLT_vs_nothing");


h0->SetLineColor(2);
h1->SetLineColor(4);
h2->SetLineColor(6);
h3->SetLineColor(8);

h0->SetMarkerColor(2);
h1->SetMarkerColor(4);
h2->SetMarkerColor(6);
h3->SetMarkerColor(8);

TLegend* lg = new TLegend(0.5, 0.6, 0.85, 0.85);
lg->SetName("lg");
lg->AddEntry(h0, "MzCut_vs_nothing", "pl");
lg->AddEntry(h1, "ZSelecLowLPt_vs_nothing", "pl");
lg->AddEntry(h2, "ZSelec_vs_nothing", "pl");
lg->AddEntry(h3, "ZSelecLowLPt_dtHLT_vs_nothing", "pl");


TCanvas* plots = new TCanvas("plots", "plots");
h0->Draw();
h1->Draw("same");
h2->Draw("same");
h3->Draw("same");
lg->Draw();



}
