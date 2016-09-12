{
TFile* fdt = new TFile("/data/XZZ/80X_Ntuple/80X_20160606_NoHLT_Skim/SingleMuon_Run2016B_PromptReco_v2.root");
TFile* fmc = new TFile("/data/XZZ/80X_Ntuple/80X_20160606_NoHLT_Skim/DYJetsToLL_M50.root");

std::string output("studybs_2016b");
char name[1000];

sprintf(name, "%s.root", output.c_str());
TFile* fout = new TFile(name, "recreate");

gROOT->ProcessLine(".x tdrstyle.C");

gStyle->SetNdivisions(505);

TTree* tdt = (TTree*)fdt->Get("tree");
TTree* tmc = (TTree*)fmc->Get("tree");


TH1D* hdt_x = new TH1D("hdt_x", "Data Vertex X",100,0.056,0.075);
TH1D* hdt_y = new TH1D("hdt_y", "Data Vertex Y",100,0.08,0.105);
TH1D* hdt_z = new TH1D("hdt_z", "Data Vertex Z",100,-25,25);

TH1D* hmc_x = new TH1D("hmc_x", "MC Vertex X",100,0.08,0.12);
TH1D* hmc_y = new TH1D("hmc_y", "MC Vertex Y",100,0.155,0.18);
TH1D* hmc_z = new TH1D("hmc_z", "MC Vertex Z",100,-25,25);

TH1D* hmc_corr_x = new TH1D("hmc_corr_x", "Corrected MC Vertex X",100,0.056,0.075);
TH1D* hmc_corr_y = new TH1D("hmc_corr_y", "Corrected MC Vertex Y",100,0.08,0.105);
TH1D* hmc_corr_z = new TH1D("hmc_corr_z", "Corrected MC Vertex Z",100,-15,15);

hdt_x->Sumw2();
hdt_y->Sumw2();
hdt_z->Sumw2();
hmc_x->Sumw2();
hmc_y->Sumw2();
hmc_z->Sumw2();
hmc_corr_x->Sumw2();
hmc_corr_y->Sumw2();
hmc_corr_z->Sumw2();

hdt_x->SetLineColor(1);
hdt_y->SetLineColor(1);
hdt_z->SetLineColor(1);
hmc_x->SetLineColor(4);
hmc_y->SetLineColor(4);
hmc_z->SetLineColor(4);
hmc_corr_x->SetLineColor(4);
hmc_corr_y->SetLineColor(4);
hmc_corr_z->SetLineColor(4);

hdt_x->SetMarkerColor(1);
hdt_y->SetMarkerColor(1);
hdt_z->SetMarkerColor(1);
hmc_x->SetMarkerColor(4);
hmc_y->SetMarkerColor(4);
hmc_z->SetMarkerColor(4);
hmc_corr_x->SetMarkerColor(4);
hmc_corr_y->SetMarkerColor(4);
hmc_corr_z->SetMarkerColor(4);

hdt_x->SetMarkerStyle(20);
hdt_y->SetMarkerStyle(20);
hdt_z->SetMarkerStyle(20);
hmc_x->SetMarkerStyle(20);
hmc_y->SetMarkerStyle(20);
hmc_z->SetMarkerStyle(20);
hmc_corr_x->SetMarkerStyle(20);
hmc_corr_y->SetMarkerStyle(20);
hmc_corr_z->SetMarkerStyle(20);

hdt_x->GetXaxis()->SetTitle("vertex x (cm)");
hdt_y->GetXaxis()->SetTitle("vertex y (cm)");
hdt_z->GetXaxis()->SetTitle("vertex z (cm)");
hmc_x->GetXaxis()->SetTitle("vertex x (cm)");
hmc_y->GetXaxis()->SetTitle("vertex y (cm)");
hmc_z->GetXaxis()->SetTitle("vertex z (cm)");
hmc_corr_x->GetXaxis()->SetTitle("vertex x (cm)");
hmc_corr_y->GetXaxis()->SetTitle("vertex y (cm)");
hmc_corr_z->GetXaxis()->SetTitle("vertex z (cm)");

TF1* fc_dt_x = new TF1("fc_dt_x", "gaus", -15,15);
TF1* fc_dt_y = new TF1("fc_dt_y", "gaus", -15,15);
TF1* fc_dt_z = new TF1("fc_dt_z", "gaus", -15,15);
TF1* fc_mc_x = new TF1("fc_mc_x", "gaus", -15,15);
TF1* fc_mc_y = new TF1("fc_mc_y", "gaus", -15,15);
TF1* fc_mc_z = new TF1("fc_mc_z", "gaus", -15,15);

fc_dt_x->SetLineWidth(3);
fc_dt_y->SetLineWidth(3);
fc_dt_z->SetLineWidth(3);
fc_mc_x->SetLineWidth(3);
fc_mc_y->SetLineWidth(3);
fc_mc_z->SetLineWidth(3);

TCanvas* plots = new TCanvas("plots", "plots", 600,600);

sprintf(name, "%s_plots.ps[", output.c_str());
plots->Print(name);

// data
tdt->Draw("vtx_x>>hdt_x");
tdt->Draw("vtx_y>>hdt_y");
tdt->Draw("vtx_z>>hdt_z");
// mc
tmc->Draw("vtx_x>>hmc_x");
tmc->Draw("vtx_y>>hmc_y");
tmc->Draw("vtx_z>>hmc_z");

// fit data 
hdt_x->Fit(fc_dt_x, "R", "", 0.061, 0.067);
hdt_y->Fit(fc_dt_y, "R", "", 0.092, 0.098);
hdt_z->Fit(fc_dt_z, "R", "", -10, 10);
hmc_x->Fit(fc_mc_x, "R", "", 0.101, 0.109);
hmc_y->Fit(fc_mc_y, "R", "", 0.165, 0.172);
hmc_z->Fit(fc_mc_z, "R", "", -10, 10);


double dt_vtx_x_mean = fc_dt_x->GetParameter(1);
double dt_vtx_y_mean = fc_dt_y->GetParameter(1);
double dt_vtx_z_mean = fc_dt_z->GetParameter(1);
double mc_vtx_x_mean = fc_mc_x->GetParameter(1);
double mc_vtx_y_mean = fc_mc_y->GetParameter(1);
double mc_vtx_z_mean = fc_mc_z->GetParameter(1);
 

double dt_vtx_x_sigma = fc_dt_x->GetParameter(2);
double dt_vtx_y_sigma = fc_dt_y->GetParameter(2);
double dt_vtx_z_sigma = fc_dt_z->GetParameter(2);
double mc_vtx_x_sigma = fc_mc_x->GetParameter(2);
double mc_vtx_y_sigma = fc_mc_y->GetParameter(2);
double mc_vtx_z_sigma = fc_mc_z->GetParameter(2);

// correction:
// starting from equation:
//    (corr_mc_vtx_x - dt_vtx_x_mean) / (mc_vtx_x - mc_vtx_x_mean) = dt_vtx_x_sigma / mc_vtx_x_sigma
// we can have:
//   corr_mc_vtx_x = dt_vtx_x_mean + (mc_vtx_x - mc_vtx_x_mean) * dt_vtx_x_sigma / mc_vtx_x_sigma
// same applies to y and z.

// prediction:
sprintf(name, "(%f+(vtx_x-%f)*%f/%f)>>hmc_corr_x", dt_vtx_x_mean, mc_vtx_x_mean, dt_vtx_x_sigma, mc_vtx_x_sigma);
tmc->Draw(name);
sprintf(name, "(%f+(vtx_y-%f)*%f/%f)>>hmc_corr_y", dt_vtx_y_mean, mc_vtx_y_mean, dt_vtx_y_sigma, mc_vtx_y_sigma);
tmc->Draw(name);
sprintf(name, "(%f+(vtx_z-%f)*%f/%f)>>hmc_corr_z", dt_vtx_z_mean, mc_vtx_z_mean, dt_vtx_z_sigma, mc_vtx_z_sigma);
tmc->Draw(name);

// norm to data
hmc_corr_x->Scale(hdt_x->Integral()/hmc_corr_x->Integral());
hmc_corr_y->Scale(hdt_y->Integral()/hmc_corr_y->Integral());
hmc_corr_z->Scale(hdt_z->Integral()/hmc_corr_z->Integral());

// plotting
TLegend* lg_x = new TLegend(0.55,0.6,0.85,0.85);
TLegend* lg_y = new TLegend(0.2,0.6,0.5,0.85);
TLegend* lg_z = new TLegend(0.15,0.6,0.45,0.85);
lg_x->SetName("lg_x");
lg_y->SetName("lg_y");
lg_z->SetName("lg_z");

lg_x->AddEntry(hdt_x, "Data", "pl");
lg_x->AddEntry(hmc_corr_x, "Corrected MC", "l");
lg_y->AddEntry(hdt_y, "Data", "pl");
lg_y->AddEntry(hmc_corr_y, "Corrected MC", "l");
lg_z->AddEntry(hdt_z, "Data", "pl");
lg_z->AddEntry(hmc_corr_z, "Corrected MC", "l");

// dt x
plots->Clear();
hdt_x->Draw();
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();
// mc x
plots->Clear();
hmc_x->Draw();
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();
// corr x
plots->Clear();
hdt_x->Draw();
hmc_corr_x->Draw("same");
lg_x->Draw();
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();

// dt y
plots->Clear();
hdt_y->Draw();
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();
// mc y
plots->Clear();
hmc_y->Draw();
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();
// corr y
plots->Clear();
hdt_y->Draw();
hmc_corr_y->Draw("same");
lg_y->Draw();
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();

// dt z
plots->Clear();
hdt_z->Draw();
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();
// mc z
plots->Clear();
hmc_z->Draw();
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();
// corr z
plots->Clear();
hdt_z->Draw();
hmc_corr_z->Draw("same");
lg_z->Draw();
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();


// 2d plots

TH2D* hdt_xy = new TH2D("hdt_xy", "Vertex X-Y Data", 100, 0,0.15, 100,0, 0.20);
TH2D* hmc_xy = new TH2D("hmc_xy", "Vertex X-Y MC", 100, 0,0.15, 100,0, 0.20);
TH2D* hdtmc_xy = new TH2D("hdtmc_xy", "Vertex X-Y Data and MC", 100, 0,0.15, 100,0, 0.20);
hdt_xy->Sumw2();
hmc_xy->Sumw2();
hdtmc_xy->Sumw2();

tdt->Draw("vtx_y:vtx_x>>hdt_xy");
tmc->Draw("vtx_y:vtx_x>>hmc_xy");

hdt_xy->GetXaxis()->SetTitle("Vertex X (cm)");
hdt_xy->GetYaxis()->SetTitle("Vertex Y (cm)");
hmc_xy->GetXaxis()->SetTitle("Vertex X (cm)");
hmc_xy->GetYaxis()->SetTitle("Vertex Y (cm)");
hdtmc_xy->GetXaxis()->SetTitle("Vertex X (cm)");
hdtmc_xy->GetYaxis()->SetTitle("Vertex Y (cm)");

hdtmc_xy->Add(hmc_xy);
hdtmc_xy->Scale(hdt_xy->Integral()/hmc_xy->Integral());
hdtmc_xy->Add(hdt_xy);

// dt xy
plots->Clear();
hdt_xy->Draw("colz");
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();

// mc xy
plots->Clear();
hmc_xy->Draw("colz");
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();

// dt+mc xy
plots->Clear();
hdtmc_xy->Draw("colz");
sprintf(name, "%s_plots.ps", output.c_str());
plots->Print(name);
plots->Clear();




// close
sprintf(name, "%s_plots.ps]", output.c_str());
plots->Print(name);

sprintf(name, ".! ps2pdf %s_plots.ps", output.c_str());
gROOT->ProcessLine(name);


fout->cd();
hdt_x->Write();
hdt_y->Write();
hdt_z->Write();
hmc_x->Write();
hmc_y->Write();
hmc_z->Write();
hmc_corr_x->Write();
hmc_corr_y->Write();
hmc_corr_z->Write();
fc_dt_x->Write();
fc_dt_y->Write();
fc_dt_z->Write();
fc_mc_x->Write();
fc_mc_y->Write();
fc_mc_z->Write();
lg_x->Write();
lg_y->Write();
lg_z->Write();

hdt_xy->Write();
hmc_xy->Write();
hdtmc_xy->Write();

fout->Close();

// print parameters

sprintf(name, "(%f+(vtx_x-%f)*%f/%f)", dt_vtx_x_mean, mc_vtx_x_mean, dt_vtx_x_sigma, mc_vtx_x_sigma);
std::cout << "corr_vtx_x = " << name << std::endl; 
sprintf(name, "(%f+(vtx_y-%f)*%f/%f)", dt_vtx_y_mean, mc_vtx_y_mean, dt_vtx_y_sigma, mc_vtx_y_sigma);
std::cout << "corr_vtx_y = " << name << std::endl; 
sprintf(name, "(%f+(vtx_z-%f)*%f/%f)", dt_vtx_z_mean, mc_vtx_z_mean, dt_vtx_z_sigma, mc_vtx_z_sigma);
std::cout << "corr_vtx_z = " << name << std::endl; 

std::cout << std::endl;
TVector3 vec_dt_bs(dt_vtx_x_mean, dt_vtx_y_mean, dt_vtx_z_mean);
TVector3 vec_mc_bs(mc_vtx_x_mean, mc_vtx_y_mean, mc_vtx_z_mean);

std::cout << "data beam spot (x,y,z) = (" << dt_vtx_x_mean << "," << dt_vtx_y_mean << "," << dt_vtx_z_mean << ")" << std::endl;
std::cout << "MC   beam spot (x,y,z) = (" << mc_vtx_x_mean << "," << mc_vtx_y_mean << "," << mc_vtx_z_mean << ")" << std::endl;
std::cout << "data beam spot (Rxy, phi,theta,eta) = (" << vec_dt_bs.Perp() << "," << vec_dt_bs.Phi() << "," << vec_dt_bs.Theta() << "," << vec_dt_bs.Eta() << ")" << std::endl;
std::cout << "MC   beam spot (Rxy, phi,theta,eta) = (" << vec_mc_bs.Perp() << "," << vec_mc_bs.Phi() << "," << vec_mc_bs.Theta() << "," << vec_mc_bs.Eta() << ")" << std::endl;

}
