{

  TFile* fin= new TFile("/home/heli/XZZ/80X_20160825_light/SingleEMU_Run2016BCD_PromptReco/vvTreeProducer/tree.root");
  std::string sel_tight = "(abs(llnunu_l1_l1_pdgId)==11&&abs(llnunu_l1_l2_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5&&llnunu_l1_mass>70&&llnunu_l1_mass<110)";
  std::string sel_loose = "(abs(llnunu_l1_l1_pdgId)==11&&abs(llnunu_l1_l2_pdgId)==11&&llnunu_l1_l1_pt>35&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5&&llnunu_l1_mass>70&&llnunu_l1_mass<110)";
  TTree* tree = (TTree*)fin->Get("tree");

  std::string sel_passHLT = "(HLT_ELEv2)";

  std::string sel_pass;
  std::string sel_base;

  std::string sel_pass_tight = "("+sel_tight+"&&"+sel_passHLT+")";
  std::string sel_base_tight = sel_tight;


  double N_pass_tight = (double)tree->GetEntries(sel_pass.c_str());
  double N_base_tight = (double)tree->GetEntries(sel_base.c_str());

  std::cout << " N_pass_tight = " << N_pass_tight << "; N_base_tight = " << N_base_tight << ";  eff = " << N_pass_tight/N_base_tight << std::endl;


  std::string sel_pass_loose = "("+sel_loose+"&&"+sel_passHLT+")";
  std::string sel_base_loose = sel_loose;


  double N_pass_loose = (double)tree->GetEntries(sel_pass_loose.c_str());
  double N_base_loose = (double)tree->GetEntries(sel_base_loose.c_str());

  std::cout << " N_pass_loose = " << N_pass_loose << "; N_base_loose = " << N_base_loose << ";  eff = " << N_pass_loose/N_base_loose << std::endl

}
