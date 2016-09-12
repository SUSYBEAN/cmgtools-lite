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
      "697", "713"
      //"575", "580", "585", "590", "595", "600", "605", "610", "615", "620", "625", "630", "633", "635", "640", "645", "646", "650", "655", "660", "665", "670", "675", "680", "685", "690", "691", "692", "693", "694", "695", "696", "697", "698", "699", "700", "701", "702", "703", "704", "705", "710", "713", "715", "720", "725", "730", "735", "740", "745", "750", "755", "760", "765"
      //"633", "646", "713"
      //"575", "580", "585", "590", "595", "600", "605", "610", "615", "620", "625", "630", "635", "640", "645", "650", "655", "660", "665", "670", "675", "680", "685", "690", "695", "700", "705", "710", "713", "715", "720", "725", "730", "735", "740", "745", "750", "755", "760", "765"
      //"660", "665", "670", "675", "680", "685", "690", "695", "700", "705", "710", "713", "715", "720", "725", "730", "735", "740", "745", "750", "755", "760", "765"
       };
  // pileup inputs
  std::vector<TFile*> pileup_files;
  std::vector<TH1D*> pileup_hists;
  std::vector<Double_t*> pileup_weights;
  std::vector<TBranch*> pileup_branches;

  for (int i=0; i<(int)pileup_tags.size(); i++) {
    sprintf(name, "pileup_MC_80x_271036-274443_%s00.root",pileup_tags.at(i).c_str());
    pileup_files.push_back(new TFile(name));
    pileup_hists.push_back((TH1D*)pileup_files.at(i)->Get("puweight"));
    pileup_weights.push_back(new Double_t(1.0));
    sprintf(name, "puWeight%s", pileup_tags.at(i).c_str());
    sprintf(name1, "puWeight%s/D", pileup_tags.at(i).c_str());
    pileup_branches.push_back(tree_out->Branch(name,pileup_weights.at(i),name1));
    

  }


  //
  Int_t nTrueInt; 
  tree->SetBranchAddress( "nTrueInt", &nTrueInt );

  for (int i=0; i<(int)tree->GetEntries(); i++){
  //for (int i=0; i<(int)1000; i++){
    tree->GetEntry(i);

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



