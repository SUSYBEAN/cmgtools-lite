#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "config.hpp"

int main(int argc, char** argv) {

  if( argc<2 ) {
     std::cout << " storehist:     " << std::endl ;
     std::cout << " Author : Hengne Li @ LPSC 2010 " << std::endl ;
     std::cout << " Functionality: Store hists for track match efficiency."  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: storehist storehist.steer" << std::endl ;
     exit(1) ;
  }

  std::cout << " storehist:     " << std::endl ;
  std::cout << " Author : Hengne Li @ LPSC 2010 " << std::endl ;

  // readme:
  // The propose of this program is to fill Mee histograms in 
  // bins of an efficiency dependence variable. An efficiency 
  // dependence varialbe can be elecPT, SET, Lumi, or whatever 
  // you want.
  // This is done in parallel for at lease two histgrams sets, 
  // i.e. hnume/hdeno or htrk0/htrk1/htrk2. For each of the 
  // histogram sets, e.g. htrk2, it contains an array of 
  // histograms each of them corresponding to a bin in this 
  // efficiency dependence variable.
  // The user has to tell the program what is the BaseSelection 
  // to be applied to all the histograms; the BaseSelectionVec 
  // to be applied deferently for each histogram set; and the 
  // binning of this given efficiency dependence variable.

  // open steering file names
  config  steer(std::string((const char*)argv[1]));
  
  // input file name
  std::vector<std::string> inputfile = steer.getStringArray("inputfile");
  // output file name
  std::string outputfile = steer.getString("outputfile");
  // Electron Basis, if true, it will draw twice for each electron
  bool ElectronBasis = steer.getBool("ElectronBasis");

  // Dimension of efficiency dependece vars: e.g. elecPt 1D , SET-Lumi 2D
  int DepVarDimension = steer.getInt("DepVarDimension");
  if (DepVarDimension>3) {
    std::cout << " Error:: DepVarDimension maximum to be 3 dimintion, no more, modify me if needed! " << std::endl;
    abort();
  }
  if (DepVarDimension==0) {
    std::cout << " Warning:: DepVarDimension should at least be 1 -- set it to be 1 " << std::endl;
    DepVarDimension = 1;
  }
   
  // Dependence Variable is of electrons?  
  //  For SET or Lumi these Event based variable, it should be set to false, binning
  //   will be written as "set>0&&set<20" ;
  //  For elecPt these electron based variable, it should be set to true, because the binning should be written as 
  //   "elecPt[0]>0&&elecPt[0]<20".
  bool ElectronDepVar1(false),ElectronDepVar2(false),ElectronDepVar3(false);
  ElectronDepVar1 = steer.getBool("ElectronDepVar");
  ElectronDepVar2 = steer.getBool("ElectronDepVar2");
  ElectronDepVar3 = steer.getBool("ElectronDepVar3");
  // ZTree name
  std::string ZTreeName = steer.getString("ZTreeName");
  // mee varaible names for each HistNames: e.g. mee
  std::vector<std::string> MeeVarNames = steer.getStringArray("MeeVarNames");
  // Mee Variable bins
  std::vector<double> MeeVarBins = steer.getDoubleArray("MeeVarBins");
  // histogram names : e.g. hnum hdeno or htrk0, htrk1, htrk2
  std::vector<std::string> HistNames = steer.getStringArray("HistNames");
  // If MeeVarNames is not an array of string one-to-one corresponding to HistNames,
  // all HistNames will use the first one in "MeeVarNames"; 
  // If does not have the MeeVarNames, it will look for "MeeVarName"
  std::string MeeVarName;
  if (MeeVarNames.size()!=HistNames.size()){
    if (MeeVarNames.size()==0) {
      MeeVarName = steer.getString("MeeVarName");
    }
    else {
      MeeVarName = MeeVarNames.at(0);
    }
    MeeVarNames.clear();
    for (int i=0; i<HistNames.size(); i++){
      MeeVarNames.push_back(MeeVarName);
    }
  }
  
  // efficiency dependence variable name
  std::vector<std::string> DepVarNames1, DepVarNames2, DepVarNames3, DepVarNames1lepII, DepVarNames2lepII, DepVarNames3lepII;
  DepVarNames1 = steer.getStringArray("DepVarNames");
  DepVarNames1lepII = steer.getStringArray("DepVarNameslepII");
  // only if DepVarDimension above 1, the following will be read 
  if (DepVarDimension>1) {
    DepVarNames2 = steer.getStringArray("DepVarNames2");
    DepVarNames2lepII = steer.getStringArray("DepVarNames2lepII");}
  if (DepVarDimension>2) {
    DepVarNames3 = steer.getStringArray("DepVarNames3");
    DepVarNames3lepII = steer.getStringArray("DepVarNames3lepII");}

  // dependence variable bins
  std::vector<double> DepVarBins1, DepVarBins2, DepVarBins3;
  DepVarBins1 = steer.getDoubleArray("DepVarBins");
  // only if DepVarDimension above 1, the following will be read 
  if (DepVarDimension>1) DepVarBins2 = steer.getDoubleArray("DepVarBins2");
  if (DepVarDimension>2) DepVarBins3 = steer.getDoubleArray("DepVarBins3");

  // base selection string to be applied everywhere 
  std::string BaseSelection = steer.getString("BaseSelection");
  // base selection array for the corresponding histnames
  std::vector<std::string> BaseSelectionVec = steer.getStringArray("BaseSelectionVec");
  // if ElectronBasis, we will need the BaseSelectionVec2 for the 2nd time fill
  std::vector<std::string> BaseSelectionVec2;
  if (ElectronBasis) BaseSelectionVec2 = steer.getStringArray("BaseSelectionVec2");
 
  // smear Mee
  bool SmearMee = steer.getBool("SmearMee");
  // Mee smearing resolution
  std::vector<double> SmearMeeResolutions;
  double SmearMeeResolution = 0;
  if (SmearMee) {
    SmearMeeResolutions = steer.getDoubleArray("SmearMeeResolutions"); 
    if (SmearMeeResolutions.size()!=HistNames.size()) {
      if (SmearMeeResolutions.size()==0){
        SmearMeeResolution = steer.getDouble("SmearMeeResolution");
      }
      else {
        SmearMeeResolution = SmearMeeResolutions.at(0);
      } 
      SmearMeeResolutions.clear();         
      for (int i=0; i<HistNames.size(); i++){
        SmearMeeResolutions.push_back(SmearMeeResolution);
      }
    }
  }
  

  // Apply weight
  bool UseWeight = steer.getBool("UseWeight");
  // weight variable name
  std::string WeightVarName = "";
  if (UseWeight) WeightVarName = steer.getString("WeightVarName");
  // 
  // debug
  std::cout << "inputfile = "; for (int i=0; i<inputfile.size(); i++) std::cout << inputfile.at(i) << " "; std::cout << std::endl;
  std::cout << "outputfile = " << outputfile << std::endl;
  std::cout << "ElectronBasis = " << ElectronBasis << std::endl;
  std::cout << "ZTreeName = " << ZTreeName << std::endl;
  std::cout << "MeeVarName = " << MeeVarName << std::endl;
  std::cout << "MeeVarBins = "; for (int i=0; i<MeeVarBins.size(); i++) std::cout << MeeVarBins.at(i) << " "; std::cout << std::endl;
  std::cout << "HistNames = "; for (int i=0; i<HistNames.size(); i++) std::cout << HistNames.at(i) << " "; std::cout << std::endl;
  std::cout << "DepVarNames = "; for (int i=0; i<HistNames.size(); i++) std::cout << DepVarNames1.at(i) << " "; std::cout << std::endl;
  std::cout << "DepVarBins = "; for (int i=0; i<DepVarBins1.size(); i++) std::cout << DepVarBins1.at(i) << " "; std::cout << std::endl;
  std::cout << "BaseSelection = " << BaseSelection << std::endl;
  std::cout << "BaseSelectionVec = "; for (int i=0; i<BaseSelectionVec.size(); i++) std::cout << BaseSelectionVec.at(i) << " "; std::cout << std::endl;
  // 

  // initialize
  // root files
  //TFile* finput = new TFile(inputfile.c_str());
  TFile* foutput = new TFile(outputfile.c_str(), "recreate");

  // tree
  //TTree* tree = (TTree*)finput->Get(ZTreeName.c_str());

  // use TChain for several input files
  TChain* tree = new TChain(ZTreeName.c_str());
  // read in files
  for (int i=0; i<inputfile.size(); i++){
    tree->Add(inputfile.at(i).c_str());
  }

  // Mee binning
  int NBinsMee = MeeVarBins.size()-1;
  double BinsMee[ (const int)MeeVarBins.size() ];
  for (int i=0; i<(int)MeeVarBins.size(); i++){
    BinsMee[i] = MeeVarBins.at(i);
  }
 

  // dependence histogram sets
  int Nbins1 = DepVarBins1.size()-1; 
  int Nbins2(1), Nbins3(1);
  if (DepVarDimension>1) Nbins2 = DepVarBins2.size()-1;
  if (DepVarDimension>2) Nbins3 = DepVarBins3.size()-1;

  // Hists in bins
  int Nhists = HistNames.size();
  std::vector< std::vector< std::vector< std::vector<TH1D*> > > > Hists;
  // average all bins, an overall reference
  //std::vector< TH1D* > RefHists ; 

  for (int i=0; i<Nhists; i++){
    std::vector < std::vector< std::vector<TH1D*> > > jhistvec;
    for (int j=0; j<Nbins1; j++){
      std::vector< std::vector<TH1D*> > khistvec;
      for (int k=0; k<Nbins2; k++){    
        std::vector<TH1D*> lhistvec;
        for (int l=0; l<Nbins3; l++){
          char hname[100];
          sprintf(hname, "%sbin%d_%d_%d", HistNames.at(i).c_str(), j, k, l);
          TH1D* hist = new TH1D(hname, hname, NBinsMee, BinsMee);
          hist->Sumw2();
          lhistvec.push_back(hist);
        }
        khistvec.push_back(lhistvec);
      }
      jhistvec.push_back(khistvec);
    }
    Hists.push_back(jhistvec);
    // overall reference hist
    //char hname[100];
    //sprintf(hname, "%sref", HistNames.at(i).c_str());
    //TH1D* hist = new TH1D(hname, hname, NBinsMee, BinsMee);
    //RefHists.push_back(hist);
  }



  // Smear Mee Resolution Tree.
  // The method to smear the Mee distribution is to generate 
  // another Tree with random numbers and add it as a Frind of
  // the Z tree.
  //std::vector<TTree*> randtree_vec;
  if (SmearMee) {
    double resolution;
    for (int i=0; i<Nhists; i++){
      char tname[50];
      sprintf(tname, "randtree_%d", i);
      TTree* randtree = new TTree(tname, tname);
      randtree->Branch("resolution", &resolution, "resolution/D");
      TRandom3 rand(1016);
      // fill the randtree
      for (int j=0; j<(int)tree->GetEntries(); j++){
        resolution = 1.+rand.Gaus(0., 1.)*SmearMeeResolutions.at(i);
        randtree->Fill();
      }
      //randtree_vec.push_back(randtree);
      // add friend
      tree->AddFriend(randtree);
    }
  }

  // fill histograms
  for (int i=0; i<Nhists; i++){

    // draw and sele string
    char draw[500], sele[2000], seletemp[2000];

    // first draw overal reference hists
/*
    //print draw string
    sprintf(draw, "%s>>%s", MeeVarNames.at(i).c_str(), RefHists.at(i)->GetName());

    // if sear Mee, add Friend smear factors tree
    if (SmearMee) {
      sprintf(draw, "(%s)*(randtree_%d.resolution)>>%s", MeeVarNames.at(i).c_str(), i, RefHists.at(i)->GetName());
    }

    // print selection string
    sprintf(sele, "(%s)&&(%s)", BaseSelection.c_str(), BaseSelectionVec.at(i).c_str());

    // if use weight, times it.
    if (UseWeight) {
      sprintf(sele, "(%s)*%s" , sele, WeightVarName.c_str());
    }

    // draw
    tree->Draw(draw, sele);

    // electron basis, draw the other electron
    if (ElectronBasis) {

      // print draw string
      sprintf(draw, "%s>>+%s", MeeVarNames.at(i).c_str(), RefHists.at(i)->GetName());

      // if sear Mee, add Friend smear factors tree
      if (SmearMee) {
        sprintf(draw, "(%s)*(randtree_%d.resolution)>>%s", MeeVarNames.at(i).c_str(), i, RefHists.at(i)->GetName());
      }

      // print selection string
      sprintf(sele, "(%s)&&(%s)", BaseSelection.c_str(), BaseSelectionVec2.at(i).c_str());

      // if use weight, times it.
      if (UseWeight) {
        sprintf(sele, "(%s)*%s" , sele, WeightVarName.c_str());
      }

      // draw
      tree->Draw(draw, sele);

    }

*/
    // draw hists in bins
    for (int j=0; j<Nbins1; j++){
      for (int k=0; k<Nbins2; k++){
        for (int l=0; l<Nbins3; l++){
          
          // print draw string
          sprintf(draw, "%s>>%s", MeeVarNames.at(i).c_str(), Hists.at(i).at(j).at(k).at(l)->GetName());    
          // if sear Mee, add Friend smear factors tree
          if (SmearMee) {
            sprintf(draw, "(%s)*(randtree_%d.resolution)>>%s", MeeVarNames.at(i).c_str(), i, Hists.at(i).at(j).at(k).at(l)->GetName());
          } 
          
          // print binning of dependence variables
          char str_depvar1[100], str_depvar2[100], str_depvar3[100]; 

          // dimension 1
          sprintf(str_depvar1, "%s>%e&&%s<%e", DepVarNames1.at(i).c_str(), DepVarBins1.at(j), DepVarNames1.at(i).c_str(), DepVarBins1.at(j+1));

          // dimension 2
          if (DepVarDimension>1) {
            sprintf(str_depvar2, "%s>%e&&%s<%e", DepVarNames2.at(i).c_str(), DepVarBins2.at(k), DepVarNames2.at(i).c_str(), DepVarBins2.at(k+1));
          }
          else {
            sprintf(str_depvar2, "1");
          }
          
          // dimension 3
          if (DepVarDimension>2) {
            sprintf(str_depvar3, "%s>%e&&%s<%e", DepVarNames3.at(i).c_str(), DepVarBins3.at(l), DepVarNames3.at(i).c_str(), DepVarBins3.at(l+1));
          }
          else {
            sprintf(str_depvar3, "1");
          }

          // print selection string
          sprintf(sele, "(%s)&&(%s)&&(%s&&%s&&%s)", BaseSelection.c_str(), BaseSelectionVec.at(i).c_str(), str_depvar1, str_depvar2, str_depvar3);

          // if use weight, times it.
          if (UseWeight) {
	    sprintf(seletemp, "(%s)*%s" , sele, WeightVarName.c_str());
	    sprintf(sele, "%s" , seletemp);
          } 

          // debug
          std::cout << " -- draw = "  << draw << std::endl;
          std::cout << " -- sele = "  << sele << std::endl;

          // draw     
          tree->Draw(draw, sele);

          // electron basis, draw the other electron
          if (ElectronBasis) {

            // print draw string
            sprintf(draw, "%s>>+%s", MeeVarNames.at(i).c_str(), Hists.at(i).at(j).at(k).at(l)->GetName());
            // if sear Mee, add Friend smear factors tree
            if (SmearMee) {
              sprintf(draw, "(%s)*(randtree_%d.resolution)>>+%s", MeeVarNames.at(i).c_str(), i, Hists.at(i).at(j).at(k).at(l)->GetName());
            }
  
            // print binning of dependence variables

            // dimension 1
            if (ElectronDepVar1) sprintf(str_depvar1, "%s>%e&&%s<%e", DepVarNames1lepII.at(i).c_str(), DepVarBins1.at(j), DepVarNames1lepII.at(i).c_str(), DepVarBins1.at(j+1));
            else sprintf(str_depvar1, "%s>%e&&%s<%e", DepVarNames1.at(i).c_str(), DepVarBins1.at(j), DepVarNames1.at(i).c_str(), DepVarBins1.at(j+1));


            // dimension 2
            if (DepVarDimension>1) {
              if (ElectronDepVar2) sprintf(str_depvar2, "%s>%e&&%s<%e", DepVarNames2lepII.at(i).c_str(), DepVarBins2.at(k), DepVarNames2lepII.at(i).c_str(), DepVarBins2.at(k+1));
              else sprintf(str_depvar2, "%s>%e&&%s<%e", DepVarNames2.at(i).c_str(), DepVarBins2.at(k), DepVarNames2.at(i).c_str(), DepVarBins2.at(k+1));
            }
            else {
              sprintf(str_depvar2, "1");
            }
   
            // dimension 3
            if (DepVarDimension>2) {
              if (ElectronDepVar3) sprintf(str_depvar3, "%s>%e&&%s<%e",DepVarNames3lepII.at(i).c_str(), DepVarBins3.at(l), DepVarNames3lepII.at(i).c_str(),DepVarBins3.at(l+1));
              else  sprintf(str_depvar3, "%s>%e&&%s<%e",DepVarNames3.at(i).c_str(), DepVarBins3.at(l), DepVarNames3.at(i).c_str(),DepVarBins3.at(l+1));
            }
            else {
              sprintf(str_depvar3, "1");
            }

            // print selection string
            sprintf(sele, "(%s)&&(%s)&&(%s&&%s&&%s)", BaseSelection.c_str(), BaseSelectionVec2.at(i).c_str(), str_depvar1, str_depvar2, str_depvar3);

            // if use weight, times it.
            if (UseWeight) {
              sprintf(seletemp, "(%s)*%s" , sele, WeightVarName.c_str());
              sprintf(sele, "%s" , seletemp);
            }

            // debug
            std::cout << " -- draw = "  << draw << std::endl;
            std::cout << " -- sele = "  << sele << std::endl;           
            
            // draw 
            tree->Draw(draw, sele);

          } // if ElectronBasis
        } // for l Depvar3
      } // for k DepVar2
    } // for j DepVar1
  } // for i Hists

  foutput->cd();
  for (int i=0; i<Nhists; i++){
    //RefHists.at(i)->Write();
    for (int j=0; j<Nbins1; j++){
      for (int k=0; k<Nbins2; k++){
        for (int l=0; l<Nbins3; l++){
          Hists.at(i).at(j).at(k).at(l)->Write();
        }
      }
    }
  }

  foutput->Close();
  return 0;
}
