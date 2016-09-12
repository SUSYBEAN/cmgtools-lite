#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TVector2.h"
#include "TEntryList.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <utility>

// Hengne Li @ CERN, 2016
 
int main(int argc, char** argv) {

  if( argc<3 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: very fast algorithm to remove duplicated events.... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " inputfile.root outputfile.root  " << std::endl ;
     exit(1) ;
  }

  // input file name
  std::string inputfile((const char*)argv[1]); 
  // output file name
  std::string outputfile((const char*)argv[2]);

  // initialize
  // root files
  TFile* finput = new TFile(inputfile.c_str());
  TFile* foutput = new TFile(outputfile.c_str(), "recreate");

  // tree
  TTree* tree = (TTree*)finput->Get("tree");


  UInt_t          run;
  ULong64_t       evt;
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("evt",&evt);

  // out_tree
  TTree* tree_out = tree->CloneTree(0);;

  TTree* tree_out_duplic;
  tree_out_duplic = tree->CloneTree(0);
  tree_out_duplic->SetName("tree_duplic");
  
  // duplication removal
  int ntotal=0;
  int nduplic=0;
  std::cout << " - remove duplicated entries" << std::endl;
  std::set< std::pair<UInt_t, ULong64_t> > evt_set; 
  std::set< std::pair<UInt_t, ULong64_t> >::iterator it_evt_set; 
  for (int i=0; i<(int)tree->GetEntries(); i++) {  
    tree->GetEntry(i);

    if ( i%50000 == 0 ) {
      std::cout << "Event " << i << "   " << std::endl;
    }

    std::pair<UInt_t, ULong64_t> evt_pair = std::make_pair(run, evt);
    it_evt_set = evt_set.find(evt_pair);
    if (it_evt_set == evt_set.end()){
      evt_set.insert(evt_pair);
      tree_out->Fill();
      ntotal++;
    }
    else {
      tree_out_duplic->Fill();
      nduplic++;
    }
  }

  std::cout << " - N total entries = " << ntotal << "; N duplicate = " << nduplic << std::endl;

  foutput->cd();
  tree_out->Write();
  tree_out_duplic->Write();
  foutput->Close();
  finput->Close();

  return 0;

}



