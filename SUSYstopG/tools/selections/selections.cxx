#define selections_cxx
#include <TROOT.h>
#include <TTree.h>
#include "selections.h"

#include <iostream>
//#include <iomanip>
#include <fstream>
#include <sstream>
//#include <algorithm>
//#include <vector>
#include <cstdio> // for printf

using namespace std;

void Selections::skimming()
{
  //infile = new TFile("./76X/BulkGravToZZToZlepZinv_narrow_800/vvTreeProducer/tree.root");
  TTree* t1 = (TTree*)infile->Get("tree");
  
  //  outfile = new TFile("fout.root", "recreate");
  readCut();
  //  cout<<"double check: "<<cutChain<<endl;
  TTree* t2 = t1->CopyTree(cutChain);

  TBranch *b_SumWeights=t2->Branch("SumWeights",&SumWeights,"SumWeights/D");
  TBranch *b_SumEvents=t2->Branch("SumEvents",&SumEvents,"SumEvents/D");

  for (int ii=0; ii<(int)t2->GetEntries(); ii++){
    b_SumWeights->Fill();
    b_SumEvents->Fill();
  }
  
  t2->Write();
  outfile->ls();
  outfile->Close();
  infile->Close();
  return;
  
}


void Selections::readCut(string config)
{
  printf("Using customized cuts...\n");
  std::ifstream ifs( config );
  std::string line;
  TString channel;
  
  while (std::getline(ifs, line)){
    if ( !line.empty() ){
      while ( line[0] == ' ' ) line.erase(0,1);
    }
    if ( !line.empty() && line[0] != '#' ){
      stringstream is(line);
      TString cut;
      is >> cut;
      if (cutChain.IsNull())  cutChain=cut;
      else{
	//if (cut.Contains("pdgId")) cutChain.Append("||"+cut);
	//else
	cutChain.Append("&&"+cut);
      }
    }
  }
  cout<<"^_^ You are using the following cuts: \n"<<cutChain<<endl;
  return;
}
