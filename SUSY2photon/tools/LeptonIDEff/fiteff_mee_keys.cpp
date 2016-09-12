#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
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
TH1* hsig;

TF1* fsig;
TF1* fbkg;
TF1* fall;

int eff_bin_min, eff_bin_max;
int fit_bin_min, fit_bin_max;
double fit_min, fit_max;
double eff_min, eff_max;
void KeysPar(const int nbin, double* par);
Double_t KeysPdf(Double_t* x, Double_t* par);
Double_t FBKG(Double_t *x, Double_t *par);
Double_t FALL(Double_t *x, Double_t *par);
void PrintPlotsFillSignalHists(TCanvas* plots, std::string OutPlotFile, int sighist_bin, int stat);

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
  // allow signal fitting shift
  bool sgshift = steer.getBool("SignalShift");
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
  std::cout << "DepVarDimension = " << DepVarDimension << std::endl;
  // dependence variable bins
  std::vector<double> DepVarBins1, DepVarBins2, DepVarBins3;
  DepVarBins1 = steer.getDoubleArray("DepVarBins");
  // only if DepVarDimension above 1, the following will be read 
  if (DepVarDimension>1) DepVarBins2 = steer.getDoubleArray("DepVarBins2");
  if (DepVarDimension>2) DepVarBins3 = steer.getDoubleArray("DepVarBins3");
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
  fit_min = FitRange.at(0);
  fit_max = FitRange.at(1);

  // efficiency calculation range
  eff_min = EffCalcRange.at(0);
  eff_max = EffCalcRange.at(1);

  // canvas to print out plots
  TCanvas* plots = new TCanvas("plots", "plots");  
  plots->SetFillColor(10);
  plots->SetFrameFillColor(10);
  plots->SetObjectStat(false);

  // Dep variable binning
  int Nbins1 = DepVarBins1.size()-1;
  int Nbins2(1), Nbins3(1);
  if (DepVarDimension>1) Nbins2 = DepVarBins2.size()-1;
  if (DepVarDimension>2) Nbins3 = DepVarBins3.size()-1;

  double Bins1[ (const int)(Nbins1+1) ];
  double Bins2[ (const int)(Nbins2+1) ];
  double Bins3[ (const int)(Nbins3+1) ];
  for (int i=0; i<(int)DepVarBins1.size(); i++){
    Bins1[i] = DepVarBins1.at(i);
  }
  if (DepVarDimension>1) {
    for (int i=0; i<(int)DepVarBins2.size(); i++){
      Bins2[i] = DepVarBins2.at(i);
    }
  }
  if (DepVarDimension>2) {
    for (int i=0; i<(int)DepVarBins3.size(); i++){
      Bins3[i] = DepVarBins3.at(i);
    }
  } 

  // dependence histogram sets
  int Nhists = DATAHistNames.size();
  std::vector< std::vector< std::vector< std::vector<TH1D*> > > > DATAHists;
  std::vector< std::vector< std::vector< std::vector<TH1D*> > > > PMCSHists;
  std::vector< std::vector< std::vector< std::vector<TH1D*> > > > BKGDHists;
  for (int i=0; i<Nhists; i++){
    std::vector< std::vector< std::vector<TH1D*> > > jDATAHistVec;
    std::vector< std::vector< std::vector<TH1D*> > > jPMCSHistVec;
    std::vector< std::vector< std::vector<TH1D*> > > jBKGDHistVec;
    for (int j=0; j<Nbins1; j++){
      std::vector< std::vector<TH1D*> > kDATAHistVec;
      std::vector< std::vector<TH1D*> > kPMCSHistVec;
      std::vector< std::vector<TH1D*> > kBKGDHistVec;
      for (int k=0; k<Nbins2; k++){
        std::vector<TH1D*> lDATAHistVec;
        std::vector<TH1D*> lPMCSHistVec;
        std::vector<TH1D*> lBKGDHistVec;
        for (int l=0; l<Nbins3; l++){
          char hname[100];
          TH1D* hist;
          sprintf(hname, "%sbin%d_%d_%d", DATAHistNames.at(i).c_str(), j, k, l);
          hist = (TH1D*)fdata->Get(hname);
          hist->Sumw2();
          lDATAHistVec.push_back(hist);
          sprintf(hname, "%sbin%d_%d_%d", PMCSHistNames.at(i).c_str(), j, k, l);
          hist = (TH1D*)fpmcs->Get(hname);
          hist->Sumw2();
          lPMCSHistVec.push_back(hist);
          sprintf(hname, "%sbin%d_%d_%d", BKGDHistNames.at(i).c_str(), j, k, l);
          hist = (TH1D*)fbkgd->Get(hname);
          hist->Sumw2();
          lBKGDHistVec.push_back(hist);
        }  
        kDATAHistVec.push_back(lDATAHistVec);
        kPMCSHistVec.push_back(lPMCSHistVec);
        kBKGDHistVec.push_back(lBKGDHistVec);
      }
      jDATAHistVec.push_back(kDATAHistVec);
      jPMCSHistVec.push_back(kPMCSHistVec);
      jBKGDHistVec.push_back(kBKGDHistVec);
    }
    DATAHists.push_back(jDATAHistVec);
    PMCSHists.push_back(jPMCSHistVec);
    BKGDHists.push_back(jBKGDHistVec);
  }

  // background subtracted histograms as a function of dependence variable
  std::vector<TH1*> SignalHists;
  for (int i=0; i<Nhists; i++){
      char hname[100];
      TH1* hist;
      sprintf(hname, "%s", SignalHistNames.at(i).c_str());
      if (DepVarDimension==1) {
        hist = new TH1D(hname, hname, Nbins1, Bins1);
      }
      else if (DepVarDimension==2) {
        hist = new TH2D(hname, hname, Nbins1, Bins1, Nbins2, Bins2);
      }
      else if (DepVarDimension==3) {
        hist = new TH3D(hname, hname, Nbins1, Bins1, Nbins2, Bins2, Nbins3, Bins3);
      }
      else{
        std::cout << "Error:: DepVarDimension can be 1 to 3, no more. Modify me if needed!" << std::endl;
        abort();
      }
      hist->Sumw2();
      SignalHists.push_back(hist);
  }

  // fit
  for (int i=0; i<Nhists; i++){
    for (int j=0; j<Nbins1; j++){
      for (int k=0; k<Nbins2; k++){
        for (int l=0; l<Nbins3; l++){
          hdata = DATAHists.at(i).at(j).at(k).at(l); 
          hpmcs = PMCSHists.at(i).at(j).at(k).at(l); 
          hbkgd = BKGDHists.at(i).at(j).at(k).at(l); 
          hsig = SignalHists.at(i);

          // fit bin range 
          fit_bin_min = hdata->FindBin(fit_min);      
          fit_bin_max = hdata->FindBin(fit_max);      
 
          // efficiency calculation bin range
          eff_bin_min = hdata->FindBin(eff_min);
          eff_bin_max = hdata->FindBin(eff_max);

          // number of data
          double Ndata_fit = hdata->Integral(fit_bin_min, fit_bin_max);
          double Npmcs_fit = hpmcs->Integral(fit_bin_min, fit_bin_max);
          double Nbkgd_fit = hbkgd->Integral(fit_bin_min, fit_bin_max);

          // initialize fitting functions

          // initialize KeysPdf for sig
          const int nbin = hpmcs->GetNbinsX();
          double sig_par[5+3*nbin];

          // get KeysPdf parameters
          if (Npmcs_fit>0) {
            KeysPar(nbin, sig_par);
          }
          else {
            for (int ik=0; ik<5+3*nbin; ik++){
              sig_par[ik] =  0.;
            }
          }

          // sig
          fsig = new TF1("fsig", KeysPdf, fit_min, fit_max, 5+3*nbin);
  
          for (int ik=0; ik<5+3*nbin; ik++){
            fsig->FixParameter(ik, sig_par[ik]);
            //std::cout << "sig_par[" << ik << "]= " << sig_par[ik] << std::endl;
          }

          // normalize fsig so as the integral is 1.
          std::cout << " -- Normalize fbkg :: " << std::endl;
          double norm_fsig = fsig->Integral(fit_min, fit_max);
          if (norm_fsig>0){
            fsig->SetParameter(0, fsig->GetParameter(0)/norm_fsig);
          }

          //   bkg
          fbkg = new TF1("fbkg", FBKG, fit_min, fit_max, 11);

          //   all
          fall = new TF1("fall", FALL, fit_min, fit_max, 2+5+3*nbin+11);

          // first fit to determine background shape
          std::cout << " -- Fit to BKG :: " << std::endl;
          for (int ik=0; ik<11; ik++){
            fbkg->SetParameter(ik, 0.);
          }
        
          fbkg->SetParameter(1, 90.);
          fbkg->SetParameter(2, 10.);

          fbkg->SetParLimits(1, 10, 500);
          fbkg->SetParLimits(2, 5, 500);
      
          if ( Nbkgd_fit>0 ){
            hbkgd->Fit(fbkg, "R N Q");
          }

          // normalize fbkg so as the integral is 1.
          std::cout << " -- Normalize fbkg :: " << std::endl;
          double norm_fbkg = fbkg->Integral(fit_min, fit_max);
          if (norm_fbkg>0){
            fbkg->SetParameter(0, fbkg->GetParameter(0)/norm_fbkg);
            for ( int ik=3; ik<11; ik++){
              fbkg->SetParameter(ik, fbkg->GetParameter(ik)/norm_fbkg);
            }
          }
          // fix parameters
          std::cout << " -- Set Pars fall :: " << std::endl;
          // par[0] : Ndata_fit 
          //fall->FixParameter(0, Ndata_fit);
          // release total
          fall->SetParameter(0, Ndata_fit);
          fall->SetParLimits(0, Ndata_fit*0.99, Ndata_fit*1.01);

          // par[1] : signal fraction initial value 0.5
          fall->SetParameter(1, 0.99);
	  fall->SetParLimits(1,0,1);
          // par[2] to par[2+5+3*nbin-1] : sig_par[0] to sig_par[5+3*nbin-1], i.e. sigal pars except normalization
          for ( int ik=2; ik<2+5+3*nbin; ik++ ){
            fall->FixParameter(ik, fsig->GetParameter(ik-2));
            //std::cout << "par[" << ik << "]= " << fall->GetParameter(ik) << std::endl;
          }

          // par[2+5+3*nbin] to par[2+5+3*nbin+3-1] : background pars
          for ( int ik=2+5+3*nbin; ik< 2+5+3*nbin+11; ik++ ){
            fall->FixParameter(ik, fbkg->GetParameter(ik-(2+5+3*nbin)));
            //std::cout << "par[" << ik << "]= " << fall->GetParameter(ik) << std::endl;
          }

          //   release frac
          fall->ReleaseParameter(1);
          //   release shift
	  if (sgshift)
	  fall->ReleaseParameter(3+2);
          //   release sigma
          fall->ReleaseParameter(4+2);


          // set par limits
          fall->SetParLimits(1, 0., 1.);
	  if (sgshift)
	  fall->SetParLimits(3+2, -1, 1);
          fall->SetParLimits(4+2, 0, 2);

          // fix frac to be 1: for FullMC fit because no background actually
          //fall->FixParameter(1, 1);

          // fix shift for the stability test
          //fall->FixParameter(3+2, -0.66);

          // set function drawing precision
          fall->SetNpx(3*nbin);
          fbkg->SetNpx(3*nbin);
          fsig->SetNpx(3*nbin);

          // do the fit to data 1st time
          std::cout << " -- Fit to Data :: " << std::endl;
          if (Ndata_fit>0 && Npmcs_fit>0 && Nbkgd_fit>0){
            hdata->Fit(fall, "R N Q");
          }
          std::cout << " -- Fit to Data :: Done" << std::endl;

          // if sigma (par 4+2) hit limits, shift err (par 3+2) set to 10 GeV (a big number)
	  if (fall->GetParameter(4+2)>=3&&sgshift) fall->SetParError(3+2, 10);

          // comment: should not fit 2nd time with shift and sigma fixed, because that error should be considered.
          // fix shift and sigma; 
          //fall->FixParameter(3+2, fall->GetParameter(3+2));
          //fall->FixParameter(4+2, fall->GetParameter(4+2));

          // do the fit to data 2nd time
          //hdata->Fit(fall, "R N Q");

          // plot hists to ps file and fill Signal Hists
          int plotstat = 0;
          if (i==0&&j==0&&k==0&&l==0) plotstat = 0;
          else if (i==Nhists-1 && j==Nbins1-1 && k==Nbins2-1 && l==Nbins3-1) plotstat = 1;
          else plotstat = 2;
      
          int sighist_bin(0);
          if (DepVarDimension==1) sighist_bin = (int)hsig->GetBin(j+1);
          if (DepVarDimension==2) sighist_bin = (int)hsig->GetBin(j+1, k+1);
          if (DepVarDimension==3) sighist_bin = (int)hsig->GetBin(j+1, k+1, l+1);
      
          std::cout << " -- Print plots and Fill Signal Hists :: " << std::endl;
          PrintPlotsFillSignalHists(plots, OutPlotFile, sighist_bin, plotstat);

          // delete fsig, fbkg, fall
          std::cout << " -- Delete funtions :: " << std::endl;
          fsig->Delete();
          fbkg->Delete();
          fall->Delete();

        } // for l
      } // for k
    } // for j
  } // for i

  foutput->cd();
  for (int i=0; i<Nhists; i++){
    SignalHists.at(i)->Write();
  }
  foutput->Close();
  return 0;
}


void KeysPar(const int nbin, double* par){

  double N = hpmcs->Integral();

  double t[nbin];
  double n[nbin];
  double w[nbin];

  for ( int i=0; i<nbin; i++ ){
    t[i] = hpmcs->GetBinCenter(i);
    n[i] = hpmcs->GetBinContent(i);
    w[i] = hpmcs->GetBinWidth(i);
  }
 
  par[0] = N;
  par[1] = N;
  par[2] = (double)nbin;
  par[3] = 0.0; // initial shift
  par[4] = 0.1; // initial sigma

  for ( int i=0; i<nbin; i++ ){
    par[i+5] = t[i];
    par[i+5+nbin] = n[i];
    par[i+5+2*nbin] = w[i];
  }

  return ;
}

Double_t KeysPdf(Double_t* x, Double_t* par){

  double N = par[1];
  const int nbin = (const int)par[2];
  double delta = par[3];
  double sigma = par[4];
  double X = x[0]-delta;

  double t[nbin];
  double n[nbin];
  double w[nbin];

  for ( int i=0; i<nbin; i++ ){
    t[i] = par[i+5];
    n[i] = par[i+5+nbin];
    w[i] = par[i+5+2*nbin];
  }

  double h[nbin];

  for ( int i=0; i<nbin; i++ ){
    if (n[i]>1e-30) {
      //h[i] = pow(4./3., 0.2) * pow(N, 0.3) * sqrt(w[i]/2./sqrt(3.)) * pow(n[i], -0.5);
      //h[i] = pow(4./3., 0.2) * pow(N, 0.5) * sqrt(w[i]/2./sqrt(3.)) * pow(n[i], -0.5);
      h[i] = 0.5*pow(4./3., 0.2) * pow(N, 0.5) * sqrt(w[i]/2./sqrt(3.)) * pow(n[i], -0.5);
    }
    else {
      h[i]=0;
    }
  }

  double var(0);

  for ( int i=0; i<nbin; i++ ){
    if ( h[i]>0 ){
      //var +=  n[i]/h[i] * exp(-(X-t[i])*(X-t[i])/(2.*h[i]*h[i]));
      //var += n[i]/sqrt(h[i]*h[i]+sigma*sigma) * exp(-(X-t[i])*(X-t[i])/(2.*(h[i]*h[i]+sigma*sigma)));
      var += n[i]/sqrt((h[i]+sigma)*(h[i]+sigma)) * exp(-(X-t[i])*(X-t[i])/(2.*(h[i]+sigma)*(h[i]+sigma)));
    }
  }

  if (N>0) {
    return par[0]*var/N/sqrt(2*TMath::Pi()) ;
  }
  else {
    return 0;
  }
  

}

Double_t FBKG(Double_t *x, Double_t *par){

  Double_t bkg = 0;

  if (par[2]>0.) {
    bkg += par[0] * exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*par[2]*par[2]));
    //bkg += par[0] * TMath::Gaus(x[0], par[1], par[2]);
    //bkg += par[0] * TMath::Landau(x[0], par[1], par[2]);
  }

  for (int i=3; i<11; i++){
    bkg += par[i]*pow(x[0], i-3);
  }

  return bkg;

}

Double_t FALL(Double_t *x, Double_t *par){

  Double_t sigPar[5+3*(int)par[2+2]]; 
  Double_t bkgPar[11];

  sigPar[0] = par[0]*par[1]*par[2];
  for (int i=1; i<5+3*(int)par[2+2]; i++){
    sigPar[i] = par[i+2];
  } 

  bkgPar[0] = par[0]*(1.-par[1])*par[2+5+3*(int)par[2+2]];
  for (int i=1; i<3; i++){
    bkgPar[i] = par[i+2+5+3*(int)par[2+2]];
  }
  for (int i=3; i<11; i++){
    bkgPar[i] = par[0]*(1.-par[1])*par[i+2+5+3*(int)par[2+2]];
  }

  Double_t sig = KeysPdf(x, sigPar);

  Double_t bkg = FBKG(x, bkgPar);

  return sig+bkg;
}

void PrintPlotsFillSignalHists(TCanvas* plots, std::string OutPlotFile, int sighist_bin, int stat=0){

  // N raw in fit range
  double Nraw_data = hdata->Integral(fit_bin_min, fit_bin_max);
  double Nraw_bkgd = hbkgd->Integral(fit_bin_min, fit_bin_max);

  // N eff in efficiency calculation range
  double Neff_data = hdata->Integral(eff_bin_min, eff_bin_max);
  double Neff_bkgd = hbkgd->Integral(eff_bin_min, eff_bin_max);

  // fractions of Neff/Nraw 
  double feff_data = 0.;
  if (Nraw_data>1e-10) feff_data = Neff_data/Nraw_data;
  double feff_bkgd = 0.;
  double feff_bkgd_err = 0.;
  if (Nraw_bkgd>1e-10) {
    feff_bkgd = Neff_bkgd/Nraw_bkgd;
    feff_bkgd_err = sqrt( Neff_bkgd*(Nraw_bkgd+Neff_bkgd)/pow(Nraw_bkgd,3) );
  }

  // Sig fraction: par[1]
  double frac = fall->GetParameter(1);
  double frac_err = fall->GetParError(1);

  // Shift : par[3-1+2]
  double shift = fall->GetParameter(3+2);
  double shift_err = fall->GetParError(3+2);

  // Sigma : par[4-1+2]
  double sigma = fall->GetParameter(4+2);
  double sigma_err = fall->GetParError(4+2);

  // N signal and N background and their errors in fit range
  double Nsig = Nraw_data*frac;
  double Nsig_err = sqrt( Nsig + pow(frac_err*Nraw_data,2) );
  double Nbkg = Nraw_data*(1.-frac);
  double Nbkg_err = sqrt( Nbkg + pow(frac_err*Nraw_data,2) );

  // N signal and N background and their errors in efficiency calculation range.
  // Calculated using Neff_data, Nraw_data and feff_bkgd.
  double Nsig_eff = Neff_data - feff_bkgd*Nraw_data + feff_bkgd*frac*Nraw_data;
  double Nbkg_eff = feff_bkgd*(1.-frac)*Nraw_data;
  //double Nsig_eff_err = sqrt( Nsig_eff + pow(frac_err*feff_bkgd*Nraw_data,2) );
  //double Nbkg_eff_err = sqrt( Nbkg_eff + pow(frac_err*feff_bkgd*Nraw_data,2) );
  // should use Nraw_data for error calculation, as the 2 lines below
  //double Nsig_eff_err = sqrt( Nsig_eff + pow(frac_err*Nraw_data,2) );
  //double Nbkg_eff_err = sqrt( Nbkg_eff + pow(frac_err*Nraw_data,2) );
  // update again
  //double Nsig_eff_err = sqrt(Nraw_data + Nbkg_eff + pow(Nraw_data*frac_err,2));
  //double Nbkg_eff_err = sqrt(Nraw_data + Nsig_eff + pow(Nraw_data*frac_err,2));
  //double Nsig_eff_err = sqrt(Neff_data + Nbkg_eff);
  //double Nbkg_eff_err = sqrt(Neff_data + Nsig_eff);
  double Nsig_eff_err = sqrt(Neff_data + pow(feff_bkgd*(1.0-frac),2)*Nraw_data + pow(feff_bkgd*Nraw_data*frac_err,2) + pow(Nraw_data*(1.0-frac)*feff_bkgd_err,2));  
  double Nbkg_eff_err = sqrt(pow(feff_bkgd*(1.0-frac),2)*Nraw_data + pow(feff_bkgd*Nraw_data*frac_err,2) + pow(Nraw_data*(1.0-frac)*feff_bkgd_err,2));  

  // set back fitted normalizations and shift to fsig and fbkg
  fsig->SetParameter(0, Nsig*fsig->GetParameter(0));
  fsig->SetParError(0, Nsig_err*fsig->GetParameter(0));
  fsig->SetParameter(3, shift);
  fsig->SetParError(3, shift_err);
  fsig->SetParameter(4, sigma);
  fsig->SetParError(4, sigma_err);

  fbkg->SetParameter(0, Nbkg*fbkg->GetParameter(0));
  fbkg->SetParError(0, Nbkg_err*fbkg->GetParameter(0));
  for (int i=3; i<11; i++){
    fbkg->SetParameter(i, Nbkg*fbkg->GetParameter(i));
  }


  // clone hists
  TH1D* hdt = (TH1D*)hdata->Clone("hdt");
  TH1D* hbk = (TH1D*)hbkgd->Clone("hbk");

  // compose and scale hists
  hbk->Scale(Nbkg/Nraw_bkgd);

  // get chi2 between data and fall
  double chi2ndf = fall->GetChisquare()/fall->GetNDF();

  // check
  if (Nsig<=0.){
    hdt->Reset();
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
  hbk->SetStats(false);

  // make ups
  hdt->SetLineColor(kBlack);
  hbk->SetLineColor(kGreen+3);
  hbk->SetFillColor(kGreen);

  hdt->SetMarkerColor(kBlack);
  hbk->SetMarkerColor(kGreen+3);

  hdt->SetMarkerStyle(20);
  hbk->SetMarkerStyle(20);

  hdt->SetMarkerSize(0.6);
  hbk->SetMarkerSize(0.6);

  hdt->GetXaxis()->SetTitle("M(ll), GeV");
  hdt->GetXaxis()->SetTitleSize(0.04);

  fsig->SetLineColor(kMagenta);
  fall->SetLineColor(kBlue);
  fbkg->SetLineColor(kGreen+1);
 
  fall->SetLineWidth(1);
  fsig->SetLineWidth(1);
  fbkg->SetLineWidth(1);

  // painting texts
  char txt_dt[100], txt_sg[100], txt_bk[100], txt_chi[100];
  char txt_frac[100], txt_shift[100], txt_sigma[100];
  char txt_n_eff[100], txt_s_eff[100], txt_b_eff[100];
  char txt_eff_n[100], txt_eff_b[100];
  sprintf(txt_dt, "Data: %d", (int)Nraw_data);
  sprintf(txt_frac, "Frac: %f+-%f", frac, frac_err);
  sprintf(txt_shift, "Shift: %f+-%f", shift, shift_err);
  sprintf(txt_sigma, "Sigma: %f+-%f", sigma, sigma_err);
  sprintf(txt_sg, "S_fit: %d+-%d", (int)Nsig, (int)Nsig_err);
  sprintf(txt_bk, "B_fit: %d+-%d", (int)Nbkg, (int)Nbkg_err);
  sprintf(txt_chi, "Chi2/Ndf: %f", chi2ndf);
  sprintf(txt_n_eff, "N_eff: %d", (int)Neff_data);
  sprintf(txt_s_eff, "S_eff: %d+-%d", (int)Nsig_eff, (int)Nsig_eff_err);
  sprintf(txt_b_eff, "B_eff: %d+-%d", (int)Nbkg_eff, (int)Nbkg_eff_err);
  sprintf(txt_eff_n, "Eff_N: %f", feff_data);
  sprintf(txt_eff_b, "Eff_B: %f", feff_bkgd);

  // legend
  TLegend* lg = new TLegend(0.6,0.6,0.99,0.99);
  lg->AddEntry(hdt, "DATA", "pl"); 
  lg->AddEntry(hbk, "Bkg", "pl"); 
  lg->AddEntry(fall, "Fit to Data", "l"); 
  lg->AddEntry(fsig, "Fit to Sig", "l"); 
  lg->AddEntry(fbkg, "Fit to Bkg", "l"); 
  lg->AddEntry(hdt, txt_dt, "");
  lg->AddEntry(hdt, txt_frac, "");
  lg->AddEntry(hdt, txt_shift, "");
  lg->AddEntry(hdt, txt_sigma, "");
  lg->AddEntry(hdt, txt_sg, "");
  lg->AddEntry(hdt, txt_bk, "");
  lg->AddEntry(hdt, txt_chi, "");
  lg->AddEntry(hdt, txt_n_eff, "");
  lg->AddEntry(hdt, txt_s_eff, "");
  lg->AddEntry(hdt, txt_b_eff, "");
  lg->AddEntry(hdt, txt_eff_n, "");
  lg->AddEntry(hdt, txt_eff_b, "");

  // draw
  hdt->Draw();
  hbk->Draw("b LF2 sames");
  fall->Draw("same");
  fsig->Draw("same");
  fbkg->Draw("same");
  lg->Draw("same");

  // print plots to ps file
  std::string plotfile;
  if (stat==0) {
    plotfile = OutPlotFile + ".ps[";
    plots->Print(plotfile.c_str());
  } 

  plotfile = OutPlotFile + ".ps";
  plots->Print(plotfile.c_str());  

  if (stat==1) {
    plotfile = OutPlotFile + ".ps]";
    plots->Print(plotfile.c_str());
    std::string makepdf = ".! ps2pdf " + OutPlotFile + ".ps " + OutPlotFile + ".pdf";
    gROOT->ProcessLine(makepdf.c_str());
  }

  // write to root file
  char plotname[100];
  sprintf(plotname,"plot_%s", hdata->GetName()); 
  plots->Write(plotname);

  // clean up hists
  hdt->Delete();
  hbk->Delete();
  lg->Delete();
}

