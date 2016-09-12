#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

int main(int argc, char** argv) {

  if( argc<6 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: skimming... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " inputfile.root outputfile.root Nevts SumWeights cut_string " << std::endl ;
     exit(1) ;
  }

  // input file name
  std::string inputfile((const char*)argv[1]); 
  // output file name
  std::string outputfile((const char*)argv[2]);


  // nevts 
  double SumEvents = atof(argv[3]);  
  // sum weights 
  double SumWeights = atof(argv[4]);  
  // cuts
  std::string Cuts((const char*)argv[5]);

  // initialize
  // root files
  TFile* finput = new TFile(inputfile.c_str());
  TFile* foutput = new TFile(outputfile.c_str(), "recreate");


  // tree
  TTree* tree = (TTree*)finput->Get("tree");

  // out_tree
  //TTree* tree_out = tree->CopyTree(Cuts.c_str());
  TTree* tree_out = tree->CloneTree(0);;

  TBranch *b_SumEvents=tree_out->Branch("SumEvents",&SumEvents,"SumEvents/D");
  TBranch *b_SumWeights=tree_out->Branch("SumWeights",&SumWeights,"SumWeights/D");


  char name[3000];
  char name1[3000];

  // pileup file tags
  std::vector<std::string> pileup_tags = {
   "62525" ,
   "61000", "61100", "61200", "61300", "61400", "61500", "61600", "61700", "61800", "61900", "62000", "62100", "62200", "62300", "62400", "62500", "62600", "62700", "62800", "62900", "63000",
     "57500", "58000", "58500", "59000", "59500", "60000", "60500", "61000", "61500", "62000", "62500", "63000", "63500", "64000", "64500", "65000", "65500", "66000", "66500", "67000", "67500", "68000", "68500", "69000", "69500", "70000", "70500", "71000", "71500", "72000", "72500", "73000", "73500", "74000", "74500", "75000", "75500", "76000", "76500", "77000", "77500"
  //    "69452"
  //  "69100", "69200", "69300", "69400", "69500", "69600", "69700", "69800", "69900", "70000", "70100", "70200", "70300", "70400", "70500"
    };
  // pileup inputs
  std::vector<TFile*> pileup_files;
  std::vector<TH1D*> pileup_hists;
  std::vector<Double_t*> pileup_weights;
  std::vector<TBranch*> pileup_branches;

  for (int i=0; i<(int)pileup_tags.size(); i++) {
    sprintf(name, "results/pileup_MC_80x_271036-275783_%s.root",pileup_tags.at(i).c_str());
    pileup_files.push_back(new TFile(name));
    pileup_hists.push_back((TH1D*)pileup_files.at(i)->Get("puweight_dtmc"));
    pileup_weights.push_back(new Double_t(1.0));
    sprintf(name, "puWeight%s", pileup_tags.at(i).c_str());
    sprintf(name1, "puWeight%s/D", pileup_tags.at(i).c_str());
    pileup_branches.push_back(tree_out->Branch(name,pileup_weights.at(i),name1));
    

  }


  //
  int n_interval = 5000;
  Int_t nTrueInt; 
  tree->SetBranchAddress( "nTrueInt", &nTrueInt );

  for (int i=0; i<(int)tree->GetEntries(); i++){
  //for (int i=0; i<(int)1000; i++){
    tree->GetEntry(i);

    if ( i%n_interval == 0 ) {
      std::cout << "Event " << i << "   " << std::endl;
    }

    for (int j=0; j<(int)pileup_tags.size(); j++){
      *(pileup_weights.at(j)) = pileup_hists.at(j)->GetBinContent(pileup_hists.at(j)->FindBin(nTrueInt));
    }

    tree_out->Fill();
  }

  foutput->cd();
  tree_out->Write();
  foutput->Close();
  finput->Close();

  return 0;

}



