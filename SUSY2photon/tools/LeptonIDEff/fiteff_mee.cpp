#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TMinuit.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#include "config.hpp"

// global vars
TH1D* hdata;
TH1D* hbkgd;
TH1D* hpmcs;
TH1D* hsig;
int eff_bin_min, eff_bin_max;
int fit_bin_min, fit_bin_max;
double GetChi2A(Double_t* par);
double GetChi2B(Double_t* par);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void PrintPlotsFillSignalHists(TCanvas* plots, std::string OutPlotFile, double frac, double frac_err, int sighist_bin, int stat);

int main(int argc, char** argv) {

  if( argc<2 ) {
     std::cout << " fiteff:     " << std::endl ;
     std::cout << " Author : Hengne Li @ LPSC 2010 " << std::endl ;
     std::cout << " Functionality: fit for track match efficiency."  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: fiteff fiteff.config" << std::endl ;
     exit(1) ;
  }

  std::cout << " fiteff:     " << std::endl ;
  std::cout << " Author : Hengne Li @ LPSC 2010 " << std::endl ;

  // readme:

  // open steering file names
  config  steer(std::string((const char*)argv[1]));
  
  // output root file name
  std::string OutRootFile = steer.getString("OutRootFile");
  // output ps plot file name
  std::string OutPlotFile = steer.getString("OutPlotFile");
  // data file name
  std::string DATAFile = steer.getString("DATAFile");
  // PMCS file name
  std::string PMCSFile = steer.getString("PMCSFile");
  // Background file name
  std::string BKGDFile = steer.getString("BKGDFile");
  // DATA hist names
  std::vector<std::string> DATAHistNames = steer.getStringArray("DATAHistNames");
  // PMCS hist names
  std::vector<std::string> PMCSHistNames = steer.getStringArray("PMCSHistNames");
  // BKGD hist names
  std::vector<std::string> BKGDHistNames = steer.getStringArray("BKGDHistNames");
  // Names of Background subtracted Data Signal hist as a function of dependence variable
  std::vector<std::string> SignalHistNames = steer.getStringArray("SignalHistNames");
  // dependence variable bins
  std::vector<double> DepVarBins = steer.getDoubleArray("DepVarBins");
  // fit range in Mee histograms
  std::vector<double> FitRange = steer.getDoubleArray("FitRange");
  // efficiency calculation range in Mee: only use number of events in this Mee range 
  // to fill the output histograms for later on efficiency calculation 
  std::vector<double> EffCalcRange = steer.getDoubleArray("EffCalcRange");

  // initialize
  // root files
  TFile* fdata = new TFile(DATAFile.c_str());
  TFile* fpmcs = new TFile(PMCSFile.c_str());
  TFile* fbkgd = new TFile(BKGDFile.c_str());
  TFile* foutput = new TFile(OutRootFile.c_str(), "recreate");

  // fit range
  double fit_min = FitRange.at(0);
  double fit_max = FitRange.at(1);

  // efficiency calculation range
  double eff_min = EffCalcRange.at(0);
  double eff_max = EffCalcRange.at(1);

  // canvas to print out plots
  TCanvas* plots = new TCanvas("plots", "plots");  
  plots->SetFillColor(10);
  plots->SetFrameFillColor(10);
  plots->SetObjectStat(false);

  // Dep variable binning
  const int Nbins = DepVarBins.size()-1;
  double Bins[ (const int)DepVarBins.size() ];
  for (int i=0; i<(int)DepVarBins.size(); i++){
    Bins[i] = DepVarBins.at(i);
  }
 
  // dependence histogram sets
  const int Nhists = DATAHistNames.size();
  std::vector< std::vector<TH1D*> > DATAHists;
  std::vector< std::vector<TH1D*> > PMCSHists;
  std::vector< std::vector<TH1D*> > BKGDHists;
  for (int i=0; i<Nhists; i++){
    std::vector<TH1D*> DATAHistVec;
    std::vector<TH1D*> PMCSHistVec;
    std::vector<TH1D*> BKGDHistVec;
    for (int j=0; j<Nbins; j++){
      char hname[100];
      TH1D* hist;
      sprintf(hname, "%sbin%d", DATAHistNames.at(i).c_str(), j);
      hist = (TH1D*)fdata->Get(hname);
      hist->Sumw2();
      DATAHistVec.push_back(hist);
      sprintf(hname, "%sbin%d", PMCSHistNames.at(i).c_str(), j);
      hist = (TH1D*)fpmcs->Get(hname);
      hist->Sumw2();
      PMCSHistVec.push_back(hist);
      sprintf(hname, "%sbin%d", BKGDHistNames.at(i).c_str(), j);
      hist = (TH1D*)fbkgd->Get(hname);
      hist->Sumw2();
      BKGDHistVec.push_back(hist);
    }
    DATAHists.push_back(DATAHistVec);
    PMCSHists.push_back(PMCSHistVec);
    BKGDHists.push_back(BKGDHistVec);
  }

  // background subtracted histograms as a function of dependence variable
  std::vector<TH1D*> SignalHists;
  for (int i=0; i<Nhists; i++){
      char hname[100];
      sprintf(hname, "%s", SignalHistNames.at(i).c_str());
      TH1D* hist = new TH1D(hname, hname, Nbins, Bins);
      hist->Sumw2();
      SignalHists.push_back(hist);
  }

  // fit
  for (int i=0; i<Nbins; i++){
    for (int j=0; j<Nhists; j++){
      hdata = DATAHists.at(j).at(i); 
      hpmcs = PMCSHists.at(j).at(i); 
      hbkgd = BKGDHists.at(j).at(i); 
      hsig = SignalHists.at(j);

      // fit bin range 
      fit_bin_min = hdata->FindBin(fit_min);      
      fit_bin_max = hdata->FindBin(fit_max);      
 
      // efficiency calculation bin range
      eff_bin_min = hdata->FindBin(eff_min);
      eff_bin_max = hdata->FindBin(eff_max);

      // initialize Minuit
      TMinuit* gMinuit = new TMinuit(10);
      gMinuit->SetFCN(fcn);
      Double_t arglist[10];
      Int_t ierflg = 0;
      arglist[0] = 5000;
      arglist[1] = 0.01;
      gMinuit->mnparm(0, "Frac", 0.5, 0.001, 0.0, 1.0, ierflg);

      // fit
      gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
 
      // Print results
      Double_t amin,edm,errdef;
      Int_t nvpar,nparx,icstat;
      gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
      gMinuit->mnprin(3, amin);

      double frac, frac_err;
      gMinuit->GetParameter(0, frac, frac_err);

      // plot hists to ps file and fill Signal Hists
      int plotstat = 0;
      if (i==0&&j==0) plotstat = 0;
      else if (i==Nbins-1 && j==Nhists-1) plotstat = 1;
      else plotstat = 2;
      
      int sighist_bin = i+1;

      PrintPlotsFillSignalHists(plots, OutPlotFile, frac, frac_err, sighist_bin, plotstat);
      gMinuit->Delete();
    }
  }

  foutput->cd();
  for (int i=0; i<Nhists; i++){
    SignalHists.at(i)->Write();
  }
  foutput->Close();
  return 0;
}


double GetChi2A(Double_t* par){

  double frac = par[0];

  double Ndata = hdata->Integral(fit_bin_min, fit_bin_max);
  double Npmcs = hpmcs->Integral(fit_bin_min, fit_bin_max);
  double Nbkgd = hbkgd->Integral(fit_bin_min, fit_bin_max);

  TH1D* hdt = (TH1D*)hdata->Clone("hdt");
  TH1D* hmc = (TH1D*)hpmcs->Clone("hmc");

  hdt->GetXaxis()->SetRange(fit_bin_min, fit_bin_max);
  hmc->GetXaxis()->SetRange(fit_bin_min, fit_bin_max);

  hmc->Scale(frac*Ndata/Npmcs);
  hmc->Add(hbkgd, (1.-frac)*Ndata/Nbkgd);

  double chi2 = hdt->Chi2Test(hmc, "UW P CHI2");
  
  hdt->Delete();
  hmc->Delete();

  return chi2;
}

double GetChi2B(Double_t* par){

  double S = par[0];
  double B = par[1];

  double Ndata = hdata->Integral(fit_bin_min, fit_bin_max);
  double Npmcs = hpmcs->Integral(fit_bin_min, fit_bin_max);
  double Nbkgd = hbkgd->Integral(fit_bin_min, fit_bin_max);


  TH1D* hdt = (TH1D*)hdata->Clone("hdt");
  TH1D* hmc = (TH1D*)hpmcs->Clone("hmc");
  hdt->GetXaxis()->SetRange(fit_bin_min, fit_bin_max);
  hmc->GetXaxis()->SetRange(fit_bin_min, fit_bin_max);

  hmc->Scale(S/Npmcs);
  hmc->Add(hbkgd, B/Nbkgd);

  double chi2 = hdt->Chi2Test(hmc, "UW P CHI2");

  hdt->Delete();
  hmc->Delete();

  return chi2;

}


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
   f = GetChi2A(par);
}

void PrintPlotsFillSignalHists(TCanvas* plots, std::string OutPlotFile, double frac, double frac_err, int sighist_bin, int stat=0){

  // N raw in fit range
  double Nraw_data = hdata->Integral(fit_bin_min, fit_bin_max);
  double Nraw_pmcs = hpmcs->Integral(fit_bin_min, fit_bin_max);
  double Nraw_bkgd = hbkgd->Integral(fit_bin_min, fit_bin_max);

  // N eff in efficiency calculation range
  double Neff_data = hdata->Integral(eff_bin_min, eff_bin_max);
  double Neff_pmcs = hpmcs->Integral(eff_bin_min, eff_bin_max);
  double Neff_bkgd = hbkgd->Integral(eff_bin_min, eff_bin_max);

  // fractions of Neff/Nraw 
  double feff_data = 0.;
  double feff_pmcs = 0.;
  double feff_bkgd = 0.;
  if (Nraw_data>1e-10) feff_data = Neff_data/Nraw_data;
  if (Nraw_pmcs>1e-10) feff_pmcs = Neff_pmcs/Nraw_pmcs;
  if (Nraw_bkgd>1e-10) feff_bkgd = Neff_bkgd/Nraw_bkgd;

  // N signal and N background and their errors in fit range
  double Nsig = frac*Nraw_data;
  double Nsig_err = sqrt( Nsig + pow(frac_err*Nraw_data,2) );
  double Nbkg = (1.-frac)*Nraw_data;
  double Nbkg_err = sqrt( Nbkg + pow(frac_err*Nraw_data,2) );

  // N signal and N background and their errors in efficiency calculation range.
  // Calculated using Neff_data, Nraw_data and feff_bkgd.
  double Nsig_eff = Neff_data - feff_bkgd*Nraw_data + feff_bkgd*frac*Nraw_data;
  double Nsig_eff_err = sqrt( Nsig_eff + pow(frac_err*feff_bkgd*Nraw_data,2) );
  double Nbkg_eff = feff_bkgd*(1.-frac)*Nraw_data;
  double Nbkg_eff_err = sqrt( Nbkg_eff + pow(frac_err*feff_bkgd*Nraw_data,2) );

  // clone hists
  TH1D* hdt = (TH1D*)hdata->Clone("hdt");
  TH1D* hmc = (TH1D*)hpmcs->Clone("hmc");
  TH1D* hbk = (TH1D*)hbkgd->Clone("hbk");

  // compose and scale hists
  hmc->Scale(Nsig/Nraw_pmcs);
  hmc->Add(hbkgd, Nbkg/Nraw_bkgd);
  hbk->Scale(Nbkg/Nraw_bkgd);

  // get chi2 between data and S+B
  double chi2ndf = hdt->Chi2Test(hmc, "UW CHI2/NDF");

  // check
  if (Nsig<1 || std::isnan(chi2ndf)){
    hdt->Reset();
    hmc->Reset();
    hbk->Reset();
    Nsig_eff = 0.;
    Nsig_eff_err = 0.;
  }

  // set output background subtracted signal hist bins
  hsig->SetBinContent(sighist_bin, Nsig_eff);
  hsig->SetBinError(sighist_bin, Nsig_eff_err); 

  // goto canvas
  plots->cd();

  // not show stat
  hdt->SetStats(false);
  hmc->SetStats(false);
  hbk->SetStats(false);

  // make ups
  hdt->SetLineColor(kRed);
  hmc->SetLineColor(kBlue);
  hbk->SetLineColor(kGreen+3);
  hbk->SetFillColor(kGreen);

  hdt->SetMarkerColor(kRed);
  hmc->SetMarkerColor(kBlue);
  hbk->SetMarkerColor(kGreen+3);

  hdt->GetXaxis()->SetTitle("M(ee), GeV");
  hdt->GetXaxis()->SetTitleSize(0.04);

  // painting texts
  char txt_dt[100], txt_frac[100], txt_sg[100], txt_bk[100], txt_chi[100];
  char txt_n_eff[100], txt_s_eff[100], txt_b_eff[100];
  char txt_eff_n[100], txt_eff_s[100], txt_eff_b[100];
  sprintf(txt_dt, "Data: %d", (int)Nraw_data);
  sprintf(txt_frac, "Frac: %f+-%f", frac, frac_err);
  sprintf(txt_sg, "S_fit: %d+-%d", (int)Nsig, (int)Nsig_err);
  sprintf(txt_bk, "B_fit: %d+-%d", (int)Nbkg, (int)Nbkg_err);
  sprintf(txt_chi, "Chi2/Ndf: %f", chi2ndf);
  sprintf(txt_n_eff, "N_eff: %d", (int)Neff_data);
  sprintf(txt_s_eff, "S_eff: %d+-%d", (int)Nsig_eff, (int)Nsig_eff_err);
  sprintf(txt_b_eff, "B_eff: %d+-%d", (int)Nbkg_eff, (int)Nbkg_eff_err);

  sprintf(txt_eff_n, "Eff_N: %f", feff_data);
  sprintf(txt_eff_s, "Eff_S: %f", feff_pmcs);
  sprintf(txt_eff_b, "Eff_B: %f", feff_bkgd);
  // legend
  TLegend* lg = new TLegend(0.6,0.6,0.9,0.904531);
  lg->AddEntry(hdt, "DATA", "pl"); 
  lg->AddEntry(hmc, "Sig+Bkg", "pl"); 
  lg->AddEntry(hbk, "Bkg", "pl"); 
  lg->AddEntry(hdt, txt_dt, "");
  lg->AddEntry(hdt, txt_frac, "");
  lg->AddEntry(hdt, txt_sg, "");
  lg->AddEntry(hdt, txt_bk, "");
  lg->AddEntry(hdt, txt_chi, "");
  lg->AddEntry(hdt, txt_n_eff, "");
  lg->AddEntry(hdt, txt_s_eff, "");
  lg->AddEntry(hdt, txt_b_eff, "");
  lg->AddEntry(hdt, txt_eff_n, "");
  lg->AddEntry(hdt, txt_eff_s, "");
  lg->AddEntry(hdt, txt_eff_b, "");

  // draw
  hdt->Draw("Ep");
  hbk->Draw("b LF2 sames");
  hmc->Draw("sames");
  lg->Draw("same");

  // print plots to ps file
  if (stat==0) {
    std::string plotfile = OutPlotFile + "[";
    plots->Print(plotfile.c_str());
  } 

  std::cout << " --  Printing plot " << hdata->GetName() << std::endl;
  plots->Print(OutPlotFile.c_str());  

  if (stat==1) {
    std::string plotfile = OutPlotFile + "]";
    plots->Print(plotfile.c_str());
  }

  // clean up hists
  hdt->Delete();
  hmc->Delete();
  hbk->Delete();
  lg->Delete();
}

