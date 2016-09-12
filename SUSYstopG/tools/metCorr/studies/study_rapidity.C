{
TFile* fin = new TFile("/data2/XZZ2/80X_20160721_EffSkim_v2/DYJetsToLL_M50.root");

TTree* tree = (TTree*)fin->Get("tree");

Float_t llnunu_l1_pt, llnunu_l1_eta, llnunu_l1_phi, llnunu_l1_mass;
tree->SetBranchAddress("llnunu_l1_pt", &llnunu_l1_pt);
tree->SetBranchAddress("llnunu_l1_eta", &llnunu_l1_eta);
tree->SetBranchAddress("llnunu_l1_phi", &llnunu_l1_phi);
tree->SetBranchAddress("llnunu_l1_mass", &llnunu_l1_mass);

TH1D* h1d1 = new TH1D("hrapidity", "hrapidity", 50, -10,10);
TH1D* h1d2 = new TH1D("heta", "heta", 50, -10,10);

h1d1->Sumw2();
h1d2->Sumw2();

for (int i=0; i<(int)tree->GetEntries(); i++){
  tree->GetEntry(i);
  TLorentzVector z;
  z.SetPtEtaPhiM(llnunu_l1_pt, llnunu_l1_eta, llnunu_l1_phi, llnunu_l1_mass);
  h1d1->Fill(z.Rapidity());
  h1d2->Fill(z.PseudoRapidity());



}

h1d1->SetLineColor(2);
h1d2->SetLineColor(4);

TCanvas* plots= new TCanvas("plots", "plots");

h1d1->Draw();
h1d2->Draw("same");


}
