
{
// check if apply HLT filter before and after lepton pair selection make any difference.
// conclusion is no, they make no difference.

TFile* filep1 = new TFile("/home/heli/work/XZZ/80X_Ntuple/80X_20160606_NoHLT_Skim/SingleMuon_Run2016B_PromptReco_v2_HLT_MU.root");
TFile* filep2 = new TFile("/home/heli/work/XZZ/80X_Ntuple/80X_20160606_NoHLT_Skim/SingleElectron_Run2016B_PromptReco_v2_HLT_ELE.root");
TFile* file1 = new TFile("/home/heli/work/XZZ/80X_Ntuple/80X_20160606_NoHLT_Skim/SingleMuon_Run2016B_PromptReco_v2.root");
TFile* file2 = new TFile("/home/heli/work/XZZ/80X_Ntuple/80X_20160606_NoHLT_Skim/SingleElectron_Run2016B_PromptReco_v2.root");

gROOT->ProcessLine(".x tdrstyle.C");

std::string sel("(((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5))&&(llnunu_l1_mass>70&&llnunu_l1_mass<110))");

std::string sel_mu = sel+"&&(abs(llnunu_l1_l1_pdgId)==13)";
std::string sel_el = sel+"&&(abs(llnunu_l1_l1_pdgId)==11)";

std::string sel_mu_ptrg = sel_mu+"&&HLT_MU";
std::string sel_el_ptrg = sel_el+"&&HLT_ELE&&!HLT_MU";

std::string sel_mu_ftrg = sel_mu+"&&!HLT_MU";
std::string sel_el_ftrg = sel_el+"&&!HLT_ELE&&!HLT_MU";

TTree* treep1 = (TTree*)filep1->Get("tree");
TTree* treep2 = (TTree*)filep2->Get("tree");
TTree* tree1 = (TTree*)file1->Get("tree");
TTree* tree2 = (TTree*)file2->Get("tree");


Double_t n_mu_prehlt = (Double_t)treep1->GetEntries(sel_mu_ptrg.c_str());
Double_t n_el_prehlt = (Double_t)treep2->GetEntries(sel_el_ptrg.c_str());
Double_t n_mu = (Double_t)tree1->GetEntries(sel_mu_ptrg.c_str());
Double_t n_el = (Double_t)tree2->GetEntries(sel_el_ptrg.c_str());

/*
Double_t n_mu_prehlt = (Double_t)treep1->GetEntries(sel_mu.c_str());
Double_t n_el_prehlt = (Double_t)treep2->GetEntries(sel_el.c_str());
Double_t n_mu = (Double_t)tree1->GetEntries(sel_mu.c_str());
Double_t n_el = (Double_t)tree2->GetEntries(sel_el.c_str());
*/
std::cout << "Mu: Npre/N = " << n_mu_prehlt/n_mu << " Npre = " << n_mu_prehlt << " N = " << n_mu << std::endl;
std::cout << "El: Npre/N = " << n_el_prehlt/n_el << " Npre = " << n_el_prehlt << " N = " << n_el << std::endl;

}



