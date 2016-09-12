{

TFile* file1 = TFile::Open("/home/heli/XZZ/80X_20160825_light_Skim/SingleEMU_Run2016BCD_PromptReco.root");
TTree* tree = (TTree*)file1->Get("tree");

std::string sel_base = "((llnunu_l1_mass>70&&llnunu_l1_mass<110)&&abs(llnunu_l1_l1_pdgId)==13&&abs(llnunu_l1_l2_pdgId)==13&&llnunu_l1_l1_pt>50&&fabs(llnunu_l1_l1_eta)<2.4&&llnunu_l1_l2_pt>20&&fabs(llnunu_l1_l2_eta)<2.4&&(llnunu_l1_l1_highPtID>=0.991||llnunu_l1_l2_highPtID>=0.99))";

std::string sel_hlt = sel_base+"&&(HLT_MUv2)";
std::string sel_trg = sel_base+"*(trgsf)";


TH1D* h_zpt_hlt = new TH1D("h_zpt_hlt", "h_zpt_hlt", 1000,0,1000);
TH1D* h_zpt_trg = new TH1D("h_zpt_trg", "h_zpt_trg", 1000,0,1000);

h_zpt_hlt->Sumw2();
h_zpt_trg->Sumw2();


tree->Draw("llnunu_l1_pt>>h_zpt_hlt", sel_hlt.c_str());
tree->Draw("llnunu_l1_pt>>h_zpt_trg", sel_trg.c_str());

h_zpt_hlt->SetLineColor(2);
h_zpt_trg->SetLineColor(4);



}
