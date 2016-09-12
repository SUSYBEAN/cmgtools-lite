


TH1D*  GetROCGraph(TH1D* Hsig, TH1D* Hbkg);



bool do2016 = false;
double zptcut = -1;

void study_signalProtection(){


std::string tag = 
"V4_doJetsCorrUseLepResPtErrSel8JetLepHardOnly"
;

if (do2016) tag = "80x_"+tag;

std::string out="study_signalProtectioni_"+tag+".root";

std::string fname1 = "skim2/BulkGravToZZToZlepZinv_narrow_1000_"+tag+".root";
std::string fname2 = "skim2/DYJetsToLL_M50_"+tag+".root";

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

h1d1[0] = new TH1D("h_ut_hard_vs_zpt_sig", "sum jet pt over z pt", 1000, 0, 5);
h1d1[1] = new TH1D("h_ut_hard_vs_zpt_bkg", "sum jet pt over z pt", 1000, 0, 5);
h1d1[2] = new TH1D("h_max_jet_pt_vs_zpt_sig", "highest jet pt over zpt" , 1000,0,5);
h1d1[3] = new TH1D("h_max_jet_pt_vs_zpt_bkg", "highest jet pt over zpt" , 1000,0,5);
h1d1[4] = new TH1D("h_1st_jet_pt_vs_zpt_sig", "1st jet pt over zpt" , 1000,0,5);
h1d1[5] = new TH1D("h_1st_jet_pt_vs_zpt_bkg", "1st jet pt over zpt" , 1000,0,5);
h1d1[6] = new TH1D("h_1st_jet_para_vs_zpt_sig", "1st jet para over zpt" , 1000,0,5);
h1d1[7] = new TH1D("h_1st_jet_para_vs_zpt_bkg", "1st jet para over zpt" , 1000,0,5);

h1d1[0]->Sumw2();
h1d1[1]->Sumw2();
h1d1[2]->Sumw2();
h1d1[3]->Sumw2();
h1d1[4]->Sumw2();
h1d1[5]->Sumw2();
h1d1[6]->Sumw2();
h1d1[7]->Sumw2();
h1d1[0]->GetXaxis()->SetTitle("[Sum jets P_{T}]/P_{T}(Z) ");
h1d1[1]->GetXaxis()->SetTitle("[Sum jets P_{T}]/P_{T}(Z) ");
h1d1[2]->GetXaxis()->SetTitle("[Biggest jets P_{T}] / P_{T}(Z) ");
h1d1[3]->GetXaxis()->SetTitle("[Biggest jets P_{T}] / P_{T}(Z) ");
h1d1[4]->GetXaxis()->SetTitle("[Leading jet P_{T}] / P_{T}(Z) ");
h1d1[5]->GetXaxis()->SetTitle("[Leading jet P_{T}] / P_{T}(Z) ");
h1d1[6]->GetXaxis()->SetTitle("[Leading jet P_{para}] / P_{T}(Z) ");
h1d1[7]->GetXaxis()->SetTitle("[Leading jet P_{para}] / P_{T}(Z) ");
h1d1[0]->GetYaxis()->SetTitle("events");
h1d1[1]->GetYaxis()->SetTitle("events");
h1d1[2]->GetYaxis()->SetTitle("events");
h1d1[3]->GetYaxis()->SetTitle("events");
h1d1[4]->GetYaxis()->SetTitle("events");
h1d1[5]->GetYaxis()->SetTitle("events");
h1d1[6]->GetYaxis()->SetTitle("events");
h1d1[7]->GetYaxis()->SetTitle("events");

h1d1[0]->SetLineColor(2);
h1d1[1]->SetLineColor(4);
h1d1[2]->SetLineColor(2);
h1d1[3]->SetLineColor(4);
h1d1[4]->SetLineColor(2);
h1d1[5]->SetLineColor(4);
h1d1[6]->SetLineColor(2);
h1d1[7]->SetLineColor(4);

h1d1[0]->SetMarkerColor(2);
h1d1[1]->SetMarkerColor(4);
h1d1[2]->SetMarkerColor(2);
h1d1[3]->SetMarkerColor(4);
h1d1[4]->SetMarkerColor(2);
h1d1[5]->SetMarkerColor(4);
h1d1[6]->SetMarkerColor(2);
h1d1[7]->SetMarkerColor(4);

h1d1[0]->SetMarkerStyle(20);
h1d1[1]->SetMarkerStyle(20);
h1d1[2]->SetMarkerStyle(20);
h1d1[3]->SetMarkerStyle(20);
h1d1[4]->SetMarkerStyle(20);
h1d1[5]->SetMarkerStyle(20);
h1d1[6]->SetMarkerStyle(20);
h1d1[7]->SetMarkerStyle(20);

TCanvas* plots = new TCanvas("plots", "plots");

sprintf(name, "%s.ps[", out.c_str());
plots->Print(name);

sprintf(name, "njet_corr>0&&llnunu_l1_pt>%f", zptcut);
tree1->Draw("ut_hard_pt_old/llnunu_l1_pt>>h_ut_hard_vs_zpt_sig", name);
tree2->Draw("ut_hard_pt_old/llnunu_l1_pt>>h_ut_hard_vs_zpt_bkg", name);
tree1->Draw("Max$(jet_corr_pt_old)/llnunu_l1_pt>>h_max_jet_pt_vs_zpt_sig", name);
tree2->Draw("Max$(jet_corr_pt_old)/llnunu_l1_pt>>h_max_jet_pt_vs_zpt_bkg", name);
tree1->Draw("jet_corr_pt_old[0]/llnunu_l1_pt>>h_1st_jet_pt_vs_zpt_sig", name );
tree2->Draw("jet_corr_pt_old[0]/llnunu_l1_pt>>h_1st_jet_pt_vs_zpt_bkg", name );
tree1->Draw("abs(Min$(jet_corr_pt_old[]*cos(jet_corr_phi[]-llnunu_l1_phi[0])))/llnunu_l1_pt[0]>>h_1st_jet_para_vs_zpt_sig", name);
tree2->Draw("abs(Min$(jet_corr_pt_old[]*cos(jet_corr_phi[]-llnunu_l1_phi[0])))/llnunu_l1_pt[0]>>h_1st_jet_para_vs_zpt_bkg", name);


h1d1[0]->Scale(1./h1d1[0]->Integral());
h1d1[1]->Scale(1./h1d1[1]->Integral());
h1d1[2]->Scale(1./h1d1[2]->Integral());
h1d1[3]->Scale(1./h1d1[3]->Integral());
h1d1[4]->Scale(1./h1d1[4]->Integral());
h1d1[5]->Scale(1./h1d1[5]->Integral());
h1d1[6]->Scale(1./h1d1[6]->Integral());
h1d1[7]->Scale(1./h1d1[7]->Integral());

lg[0] = new TLegend(0.5,0.6, 0.85,0.85);
lg[0]->AddEntry(h1d1[0], "BulkG M 1000 GeV", "pl");
lg[0]->AddEntry(h1d1[1], "DY Jets", "pl");

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d1[0]->Draw("hist");
h1d1[1]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d1[2]->Draw("hist");
h1d1[3]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d1[4]->Draw("hist");
h1d1[5]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d1[6]->Draw("hist");
h1d1[7]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();


h1d2[0] = (TH1D*)h1d1[0]->GetCumulative();
h1d2[1] =  (TH1D*)h1d1[1]->GetCumulative();
h1d2[2] =  (TH1D*)h1d1[2]->GetCumulative();
h1d2[3] =  (TH1D*)h1d1[3]->GetCumulative();
h1d2[4] =  (TH1D*)h1d1[4]->GetCumulative();
h1d2[5] =  (TH1D*)h1d1[5]->GetCumulative();
h1d2[6] =  (TH1D*)h1d1[6]->GetCumulative();
h1d2[7] =  (TH1D*)h1d1[7]->GetCumulative();

h1d2[0]->GetYaxis()->SetTitle(" rejection rate ");
h1d2[1]->GetYaxis()->SetTitle(" rejection rate ");
h1d2[2]->GetYaxis()->SetTitle(" rejection rate ");
h1d2[3]->GetYaxis()->SetTitle(" rejection rate ");
h1d2[4]->GetYaxis()->SetTitle(" rejection rate ");
h1d2[5]->GetYaxis()->SetTitle(" rejection rate ");
h1d2[6]->GetYaxis()->SetTitle(" rejection rate ");
h1d2[7]->GetYaxis()->SetTitle(" rejection rate ");


plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d2[0]->Draw("hist");
h1d2[1]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d2[2]->Draw("hist");
h1d2[3]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d2[4]->Draw("hist");
h1d2[5]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();

plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d2[6]->Draw("hist");
h1d2[7]->Draw("hist same");
lg[0]->Draw();
sprintf(name, "%s.ps", out.c_str());
plots->Print(name);
plots->SetLogy(0);
plots->Clear();


h1d3[0]= (TH1D*)GetROCGraph(h1d2[0],h1d2[1]);
h1d3[1]= (TH1D*)GetROCGraph(h1d2[2],h1d2[3]);
h1d3[2]= (TH1D*)GetROCGraph(h1d2[4],h1d2[5]);
h1d3[3]= (TH1D*)GetROCGraph(h1d2[6],h1d2[7]);

h1d3[0]->SetLineColor(2);
h1d3[1]->SetLineColor(4);
h1d3[2]->SetLineColor(6);
h1d3[3]->SetLineColor(8);

h1d3[0]->GetXaxis()->SetTitle("Signal Protection Rate");
h1d3[1]->GetXaxis()->SetTitle("Signal Protection Rate");
h1d3[2]->GetXaxis()->SetTitle("Signal Protection Rate");
h1d3[3]->GetXaxis()->SetTitle("Signal Protection Rate");
h1d3[0]->GetYaxis()->SetTitle("ZJets Efficiency");
h1d3[1]->GetYaxis()->SetTitle("ZJets Efficiency");
h1d3[2]->GetYaxis()->SetTitle("ZJets Efficiency");
h1d3[3]->GetYaxis()->SetTitle("ZJets Efficiency");

lg[1] = new TLegend(0.3, 0.3, 0.7,0.7);
lg[1]->AddEntry(h1d3[0], "Sum Jets PT / ZPT", "L");
//lg[1]->AddEntry(h1d3[1], "Largest Jet PT / ZPT", "L");
//lg[1]->AddEntry(h1d3[2], "First Jet PT / ZPT", "L");
lg[1]->AddEntry(h1d3[3], "First Jet P_{para} / ZPT", "L");



plots->cd();
plots->Clear();
//plots->SetLogy(1);
h1d3[0]->Draw("hist");
//h1d3[1]->Draw("hist same");
//h1d3[2]->Draw("hist same");
h1d3[3]->Draw("hist same");
lg[1]->Draw();
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


TH1D*  GetROCGraph(TH1D* Hsig, TH1D* Hbkg) {

    TGraph* graph = new TGraph();
    int ig=0;
    for (int i=1; i< Hsig->GetNbinsX()+1; i++){
//      std::cout << " ibin " << i << std::endl;
        double i_sig = Hsig->GetBinContent(i);
        double i_bkg = Hbkg->GetBinContent(i);
        graph->SetPoint(ig,i_sig,1-i_bkg);
        ig += 1;
    }
    char xxx[1000];
    sprintf(xxx, "roc_%s", Hsig->GetName());
    TH1D* hist = new TH1D(xxx, xxx,500,0,1);
    for (int i=1; i<=500; i++) {
        hist->SetBinContent(i, graph->Eval(hist->GetBinCenter(i)));
    }

    graph->Delete();

    return hist;

}
