{


TFile* file1 = TFile::Open("/home/heli/XZZ/80X_20160818_light_Skim/DYJetsToLL_M50_NoRecoil.root");
TFile* file2 = TFile::Open("/home/heli/XZZ/80X_20160810_light/GJet_Pt_20toInf_DoubleEMEnriched/vvTreeProducer/tree.root");

std::string outtag="study_gjets";

char name[1000];

TCanvas* plots = new TCanvas("plots", "plots");

sprintf(name, "%s.pdf[", outtag.c_str());
plots->Print(name);

TTree* tree1 = (TTree*)file1->Get("tree");
TTree* tree2 = (TTree*)file2->Get("tree");

//Double_t ZPtBins[] = {0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,400,1000};
Double_t ZPtBins[] = {20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,65,70,75,80,85,90,95,100,105,110,120,130,150,150,160,170,180,190,200,220,240,260,300,400,500,1000};
Int_t NZPtBins = sizeof(ZPtBins)/sizeof(ZPtBins[0]) - 1;

TH1D* hzptbb1 = new TH1D("hzptbb1", "hzptbb1", NZPtBins, ZPtBins);
TH1D* hzptbb2 = new TH1D("hzptbb2", "hzptbb2", NZPtBins, ZPtBins);

hzptbb1->Sumw2();
hzptbb2->Sumw2();


tree1->Draw("llnunu_l1_pt>>hzptbb1", "(1)*(genWeight*ZPtWeight)");
tree2->Draw("gjet_l1_pt>>hzptbb2", "");

hzptbb1->SetMarkerColor(2);
hzptbb2->SetMarkerColor(4);
hzptbb1->SetLineColor(2);
hzptbb2->SetLineColor(4);

hzptbb1->Scale(1./hzptbb1->Integral());
hzptbb2->Scale(1./hzptbb2->Integral());

TLegend* lgzptbb = new TLegend(0.5,0.6,0.8,0.8);
lgzptbb->SetName("lgzptbb");
lgzptbb->AddEntry(hzptbb1, "ZJets", "pl");
lgzptbb->AddEntry(hzptbb2, "GJets", "pl");

plots->Clear();

hzptbb1->Draw();
hzptbb2->Draw("same");
lgzptbb->Draw();
sprintf(name, "%s.pdf", outtag.c_str());
plots->Print(name);
plots->Clear();

TH1D* hzptbbr12 = (TH1D*)hzptbb1->Clone("hzptbbr12");
hzptbbr12->Divide(hzptbb2);



TH1D* hzpt1 = new TH1D("hzpt1", "hzpg1", 200, 20, 500);
TH1D* hzpt2 = new TH1D("hzpt2", "hzpg2", 200, 20, 500);

hzpt1->Sumw2();
hzpt2->Sumw2();


tree1->Draw("llnunu_l1_pt>>hzpt1", "(1)*(genWeight*ZPtWeight)");
tree2->Draw("gjet_l1_pt>>hzpt2", "");

hzpt1->SetMarkerColor(2);
hzpt2->SetMarkerColor(4);
hzpt1->SetLineColor(2);
hzpt2->SetLineColor(4);

hzpt1->Scale(1./hzpt1->Integral());
hzpt2->Scale(1./hzpt2->Integral());

TLegend* lgzpt = new TLegend(0.5,0.6,0.8,0.8);
lgzpt->SetName("lgzpt");
lgzpt->AddEntry(hzpt1, "ZJets", "pl");
lgzpt->AddEntry(hzpt2, "GJets", "pl");

plots->Clear();

hzpt1->Draw();
hzpt2->Draw("same");
lgzpt->Draw();
sprintf(name, "%s.pdf", outtag.c_str());
plots->Print(name);
plots->Clear();

TH1D* hzptr12 = (TH1D*)hzpt1->Clone("hzptr12");
hzptr12->Divide(hzpt2);

TH1D* hzeta1 = new TH1D("hzeta1", "hzeta1", 100, -10, 10);
TH1D* hzeta2 = new TH1D("hzeta2", "hzeta2", 100, -10, 10);

hzeta1->Sumw2();
hzeta2->Sumw2();


tree1->Draw("llnunu_l1_eta>>hzeta1", "(1)*(genWeight*ZPtWeight)");
tree2->Draw("gjet_l1_eta>>hzeta2", "");

hzeta1->SetLineColor(2);
hzeta2->SetLineColor(4);


hzeta1->Scale(1./hzeta1->Integral());
hzeta2->Scale(1./hzeta2->Integral());

hzeta1->Draw();
hzeta2->Draw("same");

TH1D* hzrap1 = new TH1D("hzrap1", "hzrap1", 60, -3, 3);
TH1D* hzrap2 = new TH1D("hzrap2", "hzrap2", 60, -3, 3);

hzrap1->Sumw2();
hzrap2->Sumw2();


tree1->Draw("llnunu_l1_rapidity>>hzrap1", "(1)*(genWeight*ZPtWeight)");
tree2->Draw("gjet_l1_rapidity>>hzrap2", "");

hzrap1->SetLineColor(2);
hzrap2->SetLineColor(4);


hzrap1->Scale(1./hzrap1->Integral());
hzrap2->Scale(1./hzrap2->Integral());

hzrap1->Draw();
hzrap2->Draw("same");

TH1D* hmt1 = new TH1D("hmt1", "hmt1", 100, 0, 200);
TH1D* hmt2 = new TH1D("hmt2", "hmt2", 100, 0, 200);

hmt1->Sumw2();
hmt2->Sumw2();


tree1->Draw("llnunu_mt>>hmt1", "(1)*(genWeight*ZPtWeight)");
tree2->Draw("gjet_mt>>hmt2", "");

hmt1->SetLineColor(2);
hmt2->SetLineColor(4);


hmt1->Scale(1./hmt1->Integral());
hmt2->Scale(1./hmt2->Integral());

hmt1->Draw();
hmt2->Draw("same");

TH1D* hmet1 = new TH1D("hmet1", "hmet1", 100, 0, 200);
TH1D* hmet2 = new TH1D("hmet2", "hmet2", 100, 0, 200);

hmet1->Sumw2();
hmet2->Sumw2();


tree1->Draw("llnunu_l2_pt>>hmet1", "(1)*(genWeight*ZPtWeight)");
tree2->Draw("gjet_l2_pt>>hmet2", "");

hmet1->SetLineColor(2);
hmet2->SetLineColor(4);


hmet1->Scale(1./hmet1->Integral());
hmet2->Scale(1./hmet2->Integral());

hmet1->Draw();
hmet2->Draw("same");

TH1D* hzmass1 = new TH1D("hzmass1", "hzmass1", 100, 50, 180);
TH1D* hzmass2 = new TH1D("hzmass2", "hzmass2", 100, 50, 180);

hzmass1->Sumw2();
hzmass2->Sumw2();


tree1->Draw("llnunu_l1_mass>>hzmass1", "(1)*(genWeight*ZPtWeight)");
tree2->Draw("gjet_l1_mass>>hzmass2", "");

hzmass1->SetLineColor(2);
hzmass2->SetLineColor(4);


hzmass1->Scale(1./hzmass1->Integral());
hzmass2->Scale(1./hzmass2->Integral());

hzmass1->Draw();
hzmass2->Draw("same");



// 3D ZMass
Double_t ZMassBins[] = {50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180};
Int_t NZMassBins = sizeof(ZMassBins)/sizeof(ZMassBins[0]) - 1;
Double_t ZRapBins[] = {-3.0,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0};
Int_t NZRapBins = sizeof(ZRapBins)/sizeof(ZRapBins[0]) - 1;

TH3D* hzmass_zpt_zrap = new TH3D("hzmass_zpt_zrap", "hzmass_zpt_zrap", NZMassBins, ZMassBins, NZPtBins,ZPtBins,NZRapBins,ZRapBins);
hzmass_zpt_zrap->Sumw2();
tree1->Draw("llnunu_l1_rapidity:llnunu_l1_pt:llnunu_l1_mass>>hzmass_zpt_zrap", "(1)*(genWeight*ZPtWeight)");


// 2D zpt zrap
TH2D* hzpt_zrap1 = new TH2D("hzpt_zrap1", "hzpt_zrap1", NZPtBins,ZPtBins,NZRapBins,ZRapBins);
TH2D* hzpt_zrap2 = new TH2D("hzpt_zrap2", "hzpt_zrap2", NZPtBins,ZPtBins,NZRapBins,ZRapBins);

hzpt_zrap1->Sumw2();
hzpt_zrap2->Sumw2();
tree1->Draw("llnunu_l1_rapidity:llnunu_l1_pt>>hzpt_zrap1", "(1)*(genWeight*ZPtWeight)");
tree2->Draw("gjet_l1_rapidity:gjet_l1_pt>>hzpt_zrap2", "");

/*

TH1D* hzeta1_wt = new TH1D("hzeta1_wt", "hzeta1_wt", 100, -10, 10);
TH1D* hzeta2_wt = new TH1D("hzeta2_wt", "hzeta2_wt", 100, -10, 10);

hzeta1_wt->Sumw2();
hzeta2_wt->Sumw2();

Int_t ngenZ;
Float_t genZ_pt;
Float_t genZ_eta;
Float_t genZ_mass;

tree2->SetBranchAddress("ngenZ", &ngenZ);
tree2->SetBranchAddress("genZ_pt", &genZ_pt);
tree2->SetBranchAddress("genZ_eta", &genZ_eta);
tree2->SetBranchAddress("genZ_mass", &genZ_mass);

for (int i=0; i<(int)tree2->GetEntries(); i++){  
  tree2->GetEntry(i);
  if (ngenZ==0) continue;

  double w= hzptr12->GetBinContent(hzptr12->FindBin(genZ_pt));

  hzeta2_wt->Fill(genZ_eta, w);



}


hzeta2_wt->Scale(1./hzeta2_wt->Integral());

hzeta1->Draw();
hzeta2_wt->Draw("same");

*/

sprintf(name, "%s.root", outtag.c_str());
TFile* fout = new TFile(name, "recreate");

hzpt1->Write();
hzpt2->Write();
hzptr12->Write();
hzptbb1->Write();
hzptbb2->Write();
hzptbbr12->Write();

hzeta1->Write();
hzeta2->Write();
hzrap1->Write();
hzrap2->Write();
hmet1->Write();
hmet2->Write();

hmt1->Write();
hmt2->Write();
//hzeta2_wt->Write();

hzmass1->Write();
hzmass2->Write();
hzmass_zpt_zrap->Write();

hzpt_zrap1->Write();
hzpt_zrap2->Write();

fout->Close();

sprintf(name, "%s.pdf]", outtag.c_str());
plots->Print(name);


}
