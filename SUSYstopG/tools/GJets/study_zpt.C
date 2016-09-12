{


//TFile* file1 = TFile::Open("/home/heli/XZZ/80X_20160822/DYJetsToLL_M50_Ext/vvTreeProducer/tree.root");
//TFile* file1 = TFile::Open("/home/heli/XZZ/80X_20160822/DYJetsToLL_M50/vvTreeProducer/tree.root");
//TFile* file2 = TFile::Open("/home/heli/XZZ/80X_20160822/DYJetsToLL_M50_MGMLM_Ext1/vvTreeProducer/tree.root");
//TFile* file1 = TFile::Open("/data2/XZZ/76X_20160705/DYJetsToLL_M50_BIG/vvTreeProducer/tree.root");
TFile* file1 = TFile::Open("/home/heli/XZZ/80X_20160818_light/DYJetsToLL_M50/vvTreeProducer/tree.root");
TFile* file2 = TFile::Open("/home/heli/XZZ/80X_20160818_light/DYJetsToLL_M50_MGMLM_Ext1/vvTreeProducer/tree.root");

//std::string tag = "study_zpt_new3";
//std::string tag = "study_zpt";
std::string tag = "study_zpt_old";
char name[1000];

//std::string sel1 = "(ngenZ>0&&abs(genLep_pdgId[0])!=15)*(genWeight)";
//std::string sel2 = "(ngenZ>0&&abs(genLep_pdgId[0])!=15)*(genWeight)";
std::string sel1 = "(ngenZ>0)*(genWeight)";
std::string sel2 = "(ngenZ>0)*(genWeight)";
TTree* tree1 = (TTree*)file1->Get("tree");
TTree* tree2 = (TTree*)file2->Get("tree");

//Double_t ZPtBins[] = {0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,300,500,1000};
Double_t ZPtBins[] = {0, 1.25, 2.5, 3.75, 5, 6.25, 7.5, 8.75, 10, 11.25, 12.5, 15, 17.5, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 220, 250, 400, 1000};
Double_t ZPtBins3k[] = {0, 1.25, 2.5, 3.75, 5, 6.25, 7.5, 8.75, 10, 11.25, 12.5, 15, 17.5, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 220, 250, 400, 3000};
Int_t NZPtBins = sizeof(ZPtBins)/sizeof(ZPtBins[0]) - 1;

TH1D* hzpt1 = new TH1D("hzpt1", "hzpg1", NZPtBins, ZPtBins);
TH1D* hzpt2 = new TH1D("hzpt2", "hzpg2", NZPtBins, ZPtBins);

//TH1D* hzpt1 = new TH1D("hzpt1", "hzpg1", 100, 0, 1500);
//TH1D* hzpt2 = new TH1D("hzpt2", "hzpg2", 100, 0, 1500);

hzpt1->Sumw2();
hzpt2->Sumw2();

tree1->Draw("genZ_pt>>hzpt1", sel1.c_str());
tree2->Draw("genZ_pt>>hzpt2", sel2.c_str());

//tree1->Draw("genZ_pt>>hzpt1", "(ngenZ>0)");
//tree2->Draw("genZ_pt>>hzpt2", "(ngenZ>0)");


hzpt1->SetLineColor(2);
hzpt2->SetLineColor(4);

hzpt1->Scale(1./hzpt1->Integral(), "width");
hzpt2->Scale(1./hzpt2->Integral(), "width");

hzpt1->Draw();
hzpt2->Draw("same");

TH1D* hzptr12 = (TH1D*)hzpt1->Clone("hzptr12");
hzptr12->Divide(hzpt2);


TH1D* hzpt3k1 = new TH1D("hzpt3k1", "hzpt3k1", NZPtBins, ZPtBins3k);
TH1D* hzpt3k2 = new TH1D("hzpt3k2", "hzpt3k2", NZPtBins, ZPtBins3k);


hzpt3k1->Sumw2();
hzpt3k2->Sumw2();

tree1->Draw("genZ_pt>>hzpt3k1", sel1.c_str());
tree2->Draw("genZ_pt>>hzpt3k2", sel2.c_str());


hzpt3k1->SetLineColor(2);
hzpt3k2->SetLineColor(4);

hzpt3k1->Scale(1./hzpt3k1->Integral(), "width");
hzpt3k2->Scale(1./hzpt3k2->Integral(), "width");

hzpt3k1->Draw();
hzpt3k2->Draw("same");

TH1D* hzpt3kr12 = (TH1D*)hzpt3k1->Clone("hzpt3kr12");
hzpt3kr12->Divide(hzpt3k2);


TH1D* hzptsb1 = new TH1D("hzptsb1", "hzptsb1", 100, 0, 1500);
TH1D* hzptsb2 = new TH1D("hzptsb2", "hzptsb2", 100, 0, 1500);

hzptsb1->Sumw2();
hzptsb2->Sumw2();

tree1->Draw("genZ_pt>>hzptsb1", sel1.c_str());
tree2->Draw("genZ_pt>>hzptsb2", sel2.c_str());


hzptsb1->SetLineColor(2);
hzptsb2->SetLineColor(4);

hzptsb1->Scale(1./hzptsb1->Integral());
hzptsb2->Scale(1./hzptsb2->Integral());

hzptsb1->Draw();
hzptsb2->Draw("same");


TH1D* hzptsbr12 = (TH1D*)hzptsb1->Clone("hzptsbr12");
hzptsbr12->Divide(hzptsb2);

TH1D* hzeta1 = new TH1D("hzeta1", "hzeta1", 100, -10, 10);
TH1D* hzeta2 = new TH1D("hzeta2", "hzeta2", 100, -10, 10);

hzeta1->Sumw2();
hzeta2->Sumw2();


tree1->Draw("genZ_eta>>hzeta1", sel1.c_str());
tree2->Draw("genZ_eta>>hzeta2", sel2.c_str());

hzeta1->SetLineColor(2);
hzeta2->SetLineColor(4);


hzeta1->Scale(1./hzeta1->Integral());
hzeta2->Scale(1./hzeta2->Integral());

hzeta1->Draw();
hzeta2->Draw("same");


TH1D* hzmass1 = new TH1D("hzmass1", "hzmass1", 100, 50, 180);
TH1D* hzmass2 = new TH1D("hzmass2", "hzmass2", 100, 50,180);

hzmass1->Sumw2();
hzmass2->Sumw2();


tree1->Draw("genZ_mass>>hzmass1", sel1.c_str());
tree2->Draw("genZ_mass>>hzmass2", sel2.c_str());

hzmass1->SetLineColor(2);
hzmass2->SetLineColor(4);


hzmass1->Scale(1./hzmass1->Integral());
hzmass2->Scale(1./hzmass2->Integral());

hzmass1->Draw();
hzmass2->Draw("same");

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

sprintf(name, "%s.root", tag.c_str());
TFile* fout = new TFile(name, "recreate");

hzpt1->Write();
hzpt2->Write();
hzptr12->Write();
hzpt3k1->Write();
hzpt3k2->Write();
hzpt3kr12->Write();
hzptsb1->Write();
hzptsb2->Write();
hzptsbr12->Write();

hzeta1->Write();
hzeta2->Write();
hzmass1->Write();
hzmass2->Write();
//hzeta2_wt->Write();

fout->Close();


}
