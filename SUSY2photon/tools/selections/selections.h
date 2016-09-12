#ifndef selections_h
#define selections_h

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TChain.h>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>      // std::ifstream
#include <sys/stat.h>   // mkdir()

class Selections {
 public :
  TFile  *infile; 
  TFile  *outfile;

  Selections(TString indir="", TString outname="");
  virtual ~Selections();
  
  virtual void readCut(std::string config="config.txt");
  virtual void skimming();
  virtual void mymkdir(TString dirname);
  virtual void SetSumWeight(TString indir);
  
 private:
  TString cutChain; //="llnunu_l1_pt>200";
  Double_t SumWeights=-999.0;
  Double_t SumEvents=-999.0;
};

Selections::Selections(TString indir, TString outname){
  TString fname;
  if (indir.IsNull()) fname="./76X/BulkGravToZZToZlepZinv_narrow_800/vvTreeProducer/tree.root";
  else fname=indir+"/vvTreeProducer/tree.root";

  infile = new TFile(fname);
  std::cout << "[Info] You are openning: " << fname << std::endl;

  TString outdir="AnalysisRegion";
  if(gSystem->AccessPathName(outdir)) mymkdir(outdir);
  
  if (outname.IsNull()) outname=outdir+"/fout.root";
  else if (!outname.Contains(".root")) outname=outdir+"/"+outname+".root";
  else outname=outdir+"/"+outname;
  outfile = new TFile(outname, "recreate");
  std::cout << "[Info] You are saving reduced tree to " << outname << std::endl;

  SetSumWeight(indir);
}

Selections::~Selections(){
  std::cout << "Object is being deleted" << std::endl;
}

void Selections::mymkdir(TString dirname){
  const int dir_err = mkdir(dirname.Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (-1 == dir_err)
    {
      printf("Error creating directory!\n");
      exit(1);
    }
}


void Selections::SetSumWeight(TString indir){
  //--- To read numbers in the format of /skimAnalyzerCount/SkimReport.txt
  indir+="/skimAnalyzerCount/SkimReport.txt";

  ifstream in(indir.Data());
  if (!in) return;
  string line;
  while ( getline(in,line) ) {
    if ( !line.empty() ) { // skip leading blanks (in case there are some in front of a comment)
      while ( line[0] == ' ' ) line.erase(0,1);
    }
    if ( !line.empty() ){
      stringstream is(line);
      TString str1, str2;
      Double_t num;
      is >> str1 >> str2 >> num;
      str1+=str2;
      if (str1.EqualTo("SumWeights",TString::ECaseCompare::kIgnoreCase)) SumWeights=num;
      if (str1.EqualTo("AllEvents",TString::ECaseCompare::kIgnoreCase)) SumEvents=num;
    }
  }
  cout<<fixed<<SumWeights<<" "<<SumEvents<<endl;
  return;
}
#endif
