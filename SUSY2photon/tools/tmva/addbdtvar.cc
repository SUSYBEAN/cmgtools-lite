#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "config.h"

// author Hengne Li 

int main(int argc, char** argv) {

  if( argc<4 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: add mva variables to the ntuple... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " inputfile.root outputfile.root config_file " << std::endl ;
     exit(1) ;
  }

  // input file name
  std::string inputfile((const char*)argv[1]); 

  // output file name
  std::string outputfile((const char*)argv[2]);

  //config file 
  config cfg(std::string((const char*)argv[3]));

  // cuts
  std::string cuts(cfg.getString("cuts"));

  // tags 
  std::vector<std::string> tags = cfg.getStringArray("tags");  
 
  // weight files
  std::vector<std::string> weight_file_names;
  for (int itg=0; itg<(int)tags.size(); itg++) {
    std::string a_weight_file_name(cfg.getString(tags.at(itg)));
    std::cout << "tag = " << tags.at(itg) << " weight file = " << a_weight_file_name << std::endl;
    weight_file_names.push_back(a_weight_file_name);
  }
  
 

  // initialize

  // root files
  TFile* finput = new TFile(inputfile.c_str());
  TFile* foutput = new TFile(outputfile.c_str(), "recreate");

  // tree
  TTree* tree = (TTree*)finput->Get("tree");

  // out_tree
  TTree* tree_out = tree->CopyTree(cuts.c_str());

  // define variables
  Float_t zpt, met;
  Float_t llnunu_mt, cosDPhiZMET, dPTParaRel;

  // bdt variables
  Float_t bdt_var;

  // book bdt branches
  std::vector< TBranch* > bdt_var_branches;

  // set branch address
  tree_out->SetBranchAddress( "llnunu_l1_pt", &zpt );
  tree_out->SetBranchAddress( "met_pt", &met );
  
  // TMVA readers
  std::vector< TMVA::Reader* > tmva_readers;

  for (int itg=0; itg<(int)tags.size(); itg++) {

    // define a reader
    TMVA::Reader* a_reader = new TMVA::Reader("!Color:!Silent");

    // add variables
    a_reader->AddVariable("llnunu_l1_pt", &zpt);
    a_reader->AddVariable("met_pt", &met);

    a_reader->AddSpectator( "llnunu_mt",  &llnunu_mt );
    a_reader->AddSpectator( "cosDPhiZMET := cos(llnunu_deltaPhi)", &cosDPhiZMET );
    a_reader->AddSpectator( "dPTParaRel := abs(llnunu_l1_pt+met_pt*cos(llnunu_deltaPhi))/llnunu_l1_pt", &dPTParaRel );

    // Define the MVA type and read the weights file
    a_reader->BookMVA("BDT",weight_file_names.at(itg));

    // push back
    tmva_readers.push_back(a_reader);

    // add new branch to the output tree for this BDT tag
    TString bdt_var_name = "bdt_"+tags.at(itg);
    TBranch *b_bdt_var = tree_out->Branch(bdt_var_name, &bdt_var, bdt_var_name+TString("/F"));


    // push back
    bdt_var_branches.push_back(b_bdt_var);    
     
  }

  //  

  // loop over events and calculate bdt vars
  for (int i=0; i<(int)tree_out->GetEntries(); i++){ 
    tree_out->GetEntry(i);

    // loop over bdt vars
    for (int itg=0; itg<(int)tags.size(); itg++){
      bdt_var = tmva_readers.at(itg)->EvaluateMVA("BDT");
      bdt_var_branches.at(itg)->Fill();
    }

  }

  tree_out->Write();
  foutput->Close();
  finput->Close();

  // delete readers
  for (int itg=0; itg<(int)tags.size(); itg++){
    delete tmva_readers.at(itg);
  }

  return 0;

}



