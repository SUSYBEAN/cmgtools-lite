{


TFile* file1 = TFile::Open("80X_20160705/DYJetsToLL_M50/vvTreeProducer/tree.root");
TFile* file2 = TFile::Open("80X_20160705/DYJetsToLL_M50_MGMLM_BIG/vvTreeProducer/tree.root");


TTree* tree1 = (TTree*)file1->Get("tree");
TTree* tree2 = (TTree*)file2->Get("tree");

TH1D* hzpt1 = new TH1D("hzpt1", "hzpg1", 100, 0, 1000);
TH1D* hzpt2 = new TH1D("hzpt2", "hzpg2", 100, 0, 1000);

hzpt1->Sumw2();
hzpt2->Sumw2();


tree1->Draw("genZ_pt>>hzpt1", "ngenZ>0");
tree2->Draw("genZ_pt>>hzpt2", "ngenZ>0");

hzpt1->SetLineColor(2);
hzpt2->SetLineColor(4);

hzpt1->Scale(1./hzpt1->Integral());
hzpt2->Scale(1./hzpt2->Integral());

hzpt1->Draw();
hzpt2->Draw("same");


TH1D* hzptr12 = (TH1D*)hzpt1->Clone("hzptr12");
hzptr12->Divide(hzpt2);

TH1D* hzeta1 = new TH1D("hzeta1", "hzeta1", 100, -10, 10);
TH1D* hzeta2 = new TH1D("hzeta2", "hzeta2", 100, -10, 10);

hzeta1->Sumw2();
hzeta2->Sumw2();


tree1->Draw("genZ_eta>>hzeta1", "ngenZ>0");
tree2->Draw("genZ_eta>>hzeta2", "ngenZ>0");

hzeta1->SetLineColor(2);
hzeta2->SetLineColor(4);


hzeta1->Scale(1./hzeta1->Integral());
hzeta2->Scale(1./hzeta2->Integral());

hzeta1->Draw();
hzeta2->Draw("same");


TH1D* hzmass1 = new TH1D("hzmass1", "hzmass1", 100, -10, 10);
TH1D* hzmass2 = new TH1D("hzmass2", "hzmass2", 100, -10, 10);

hzmass1->Sumw2();
hzmass2->Sumw2();


tree1->Draw("genZ_mass>>hzmass1", "ngenZ>0");
tree2->Draw("genZ_mass>>hzmass2", "ngenZ>0");

hzmass1->SetLineColor(2);
hzmass2->SetLineColor(4);


hzmass1->Scale(1./hzmass1->Integral());
hzmass2->Scale(1./hzmass2->Integral());

hzmass1->Draw();
hzmass2->Draw("same");



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

hzeta1->Draw();
hzeta2_wt->Draw("same");


}
