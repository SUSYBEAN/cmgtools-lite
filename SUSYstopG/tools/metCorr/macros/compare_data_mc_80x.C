




void compare_data_mc(){


std::string tag = 
"doJetsCorrUseLepResPtErrSel8JetLepHardOnlyZpT250"
;

std::string out="compare_data_mc_"+tag+".root";

std::string fname1 = "skim2/SingleEMU_Run2015_16Dec_V4_"+tag+".root";
std::string fname2 = "skim2/DYJetsToLL_M50_V4_"+tag+".root";

char name[1000];
gROOT->ProcessLine(".x tdrstyle.C");

TFile* file1 = new TFile(fname1.c_str());
TFile* file2 = new TFile(fname2.c_str());

TTree* tree1 = (TTree*)file1->Get("tree");
TTree* tree2 = (TTree*)file2->Get("tree");


TH1D* h1d1[100];
TH1D* h1d2[100];
TH1D* h1d3[100];
TH1D* h1d4[100];
TH1D* h1d5[100];

TLegend* lg[100];

h1d1[0] = (TH1D*)file1->Get("hmet_para_new");
h1d1[1] = (TH1D*)file1->Get("hmet_perp_new");
h1d1[2] = (TH1D*)file1->Get("hdphi_new");
h1d1[3] = (TH1D*)file1->Get("hmet_new");
h1d1[4] = (TH1D*)file1->Get("hmt_new");


h1d2[0] = (TH1D*)file2->Get("hmet_para_new");
h1d2[1] = (TH1D*)file2->Get("hmet_perp_new");
h1d2[2] = (TH1D*)file2->Get("hdphi_new");
h1d2[3] = (TH1D*)file2->Get("hmet_new");
h1d2[4] = (TH1D*)file2->Get("hmt_new");


h1d1[3]->GetXaxis()->SetRangeUser(0, 200);
h1d2[3]->GetXaxis()->SetRangeUser(0, 200);
h1d1[4]->GetXaxis()->SetRangeUser(100, 500);
h1d2[4]->GetXaxis()->SetRangeUser(100, 500);


h1d1[0]->SetMarkerColor(2);
h1d1[1]->SetMarkerColor(2);
h1d1[2]->SetMarkerColor(2);
h1d1[3]->SetMarkerColor(2);
h1d1[4]->SetMarkerColor(2);

h1d1[0]->SetLineColor(2);
h1d1[1]->SetLineColor(2);
h1d1[2]->SetLineColor(2);
h1d1[3]->SetLineColor(2);
h1d1[4]->SetLineColor(2);

h1d2[0]->SetMarkerColor(4);
h1d2[1]->SetMarkerColor(4);
h1d2[2]->SetMarkerColor(4);
h1d2[3]->SetMarkerColor(4);
h1d2[4]->SetMarkerColor(4);

h1d2[0]->SetLineColor(4);
h1d2[1]->SetLineColor(4);
h1d2[2]->SetLineColor(4);
h1d2[3]->SetLineColor(4);
h1d2[4]->SetLineColor(4);

/*
h1d1[0]->SetMarkerStyle(20);
h1d1[1]->SetMarkerStyle(20);
h1d1[2]->SetMarkerStyle(20);
h1d1[3]->SetMarkerStyle(20);
h1d1[4]->SetMarkerStyle(20);
h1d1[5]->SetMarkerStyle(20);
h1d1[6]->SetMarkerStyle(20);
h1d1[7]->SetMarkerStyle(20);
*/
TCanvas* plots = new TCanvas("plots", "plots");

sprintf(name, "%s.ps[", out.c_str());
plots->Print(name);


h1d2[0]->Scale(h1d1[0]->Integral()/h1d2[0]->Integral());
h1d2[1]->Scale(h1d1[1]->Integral()/h1d2[1]->Integral());
h1d2[2]->Scale(h1d1[2]->Integral()/h1d2[2]->Integral());
h1d2[3]->Scale(h1d1[3]->Integral()/h1d2[3]->Integral());
h1d2[4]->Scale(h1d1[4]->Integral()/h1d2[4]->Integral());

lg[0] = new TLegend(0.7,0.7, 0.85,0.85);
lg[0]->AddEntry(h1d1[0], "Data ", "pl");
lg[0]->AddEntry(h1d2[0], "DY Jets", "pl");

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d1[0]->Draw("hist");
h1d2[0]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d1[1]->Draw("hist");
h1d2[1]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d1[2]->Draw("hist");
h1d2[2]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d1[3]->Draw("hist");
h1d2[3]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();


plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d1[4]->Draw("hist");
h1d2[4]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();


// log
plots->cd();
plots->Clear();
plots->SetLogy(1);
h1d1[0]->Draw("hist");
h1d2[0]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
plots->SetLogy(1);
h1d1[1]->Draw("hist");
h1d2[1]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
plots->SetLogy(1);
h1d1[2]->Draw("hist");
h1d2[2]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
plots->SetLogy(1);
h1d1[3]->Draw("hist");
h1d2[3]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();


plots->cd();
plots->Clear();
plots->SetLogy(1);
h1d1[4]->Draw("hist");
h1d2[4]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();



sprintf(name, "%s.ps]", out.c_str());
plots->Print(name);

//sprintf(name, ".! ps2pdf %s.ps %s.pdf", out.c_str(), out.c_str());
sprintf(name, ".! ps2pdf %s.ps ", out.c_str());
gROOT->ProcessLine(name);


}


