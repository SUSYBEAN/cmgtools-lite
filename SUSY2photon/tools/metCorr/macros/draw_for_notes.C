{
TFile* file1 = new TFile(
//"skim/BulkGravToZZToZlepZinv_narrow_1000_V3MetShiftBeforeJetNewSelV6NoMetLepAnyWayAllJetsBigSig1p4LepResSigProtect.root"
"skim/DYJetsToLL_M50_V3JetNewSelV6NoMetLepAnyWayAllJetsBigSig1p4LepRes.root"
);


gROOT->ProcessLine(".x tdrstyle.C");

TCanvas* plots = new TCanvas("plots", "plots");

TH1D* h1[100];
TH1D* h2[100];
h1[0] = (TH1D*)file->Get("hmta_old");
h2[0] = (TH1D*)file->Get("hmta_new");

h1[0]->GetXaxis()->SetTitle("M_{T} (GeV)");
h2[0]->GetXaxis()->SetTitle("M_{T} (GeV)");

TLegend *lg = new TLegend(0.5, 0.6, 0.85, 0.85);
lg->AddEntry(h1, "Before Fit", "pl");
lg->AddEntry(h2, "After Fit", "pl");
h1[0]->Draw();
h2[0]->Draw("same");



}
