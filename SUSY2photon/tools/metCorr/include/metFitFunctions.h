#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/MnUserParameters.h"
#include "TComplex.h"
#include "TH1D.h"
#include "math.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <time.h>

bool debug=false;

// 
class MetChi2Fcn : public ROOT::Minuit2::FCNBase {

public:
  MetChi2Fcn(void) {};
  ~MetChi2Fcn() {}

  void InitFunction(const double& z_pt, const double& z_phi, 
              const double& met_pt, const double& met_phi,
              const std::vector<double>& jets_pt, const std::vector<double>& jets_phi,
              const std::vector<double>& jets_sigma, std::string option)
  {
     _z_pt = z_pt;
     _z_phi = z_phi;
     _met_pt = met_pt;
     _met_phi = met_phi;
     _jets_pt = jets_pt;
     _jets_phi = jets_phi;
     _jets_sigma = jets_sigma;
     _njets = (int)jets_pt.size();

     _option = option;
  }

  virtual double operator()(const std::vector<double>& par) const
  {
     return getChi2(par);
  }


  double getChi2(const std::vector<double>& par) const
  {
    double sigma_scale=1.4;//1.4;//1.7; //1.4; //1.25886; //1.33924; // 3.; // 1.11754;
    //double sigma_met_para = 15.3016 + 0.0262152*_z_pt;
    //double sigma_met_perp = 15.2906 + 0.025963 *_z_pt;
    // 1.11754 = 16.3326/14.6148 
    //double sigma_met_para = 1.11754*1.11754*(13.1218 + 0.0457877*_z_pt);
    //double sigma_met_perp = 1.11754*1.11754*(13.2456 + 0.0189742*_z_pt);
    //double sigma_met_para = (13.2456 + 0.0189742*_z_pt);
    //double sigma_met_perp = (13.2456 + 0.0189742*_z_pt);
    //double sigma_met_para = 1.11754*(13.2456 + 0.0189742*_z_pt);
    //double sigma_met_perp = 1.11754*(13.2456 + 0.0189742*_z_pt);
    double sigma_met_para = sigma_scale*(13.2456 + 0.0189742*_z_pt);
    double sigma_met_perp = sigma_scale*(13.2456 + 0.0189742*_z_pt);
    double new_met_para = _met_pt*cos(_met_phi-_z_phi);
    double new_met_perp = _met_pt*sin(_met_phi-_z_phi);
    double chi2 = 0;
    for (int i=0; i<(int)_njets; i++){
      chi2 += pow((par.at(i)-1.0)/_jets_sigma.at(i), 2);
      new_met_para += (1.0-par.at(i))*_jets_pt.at(i)*cos(_jets_phi.at(i)-_z_phi);
      new_met_perp += (1.0-par.at(i))*_jets_pt.at(i)*sin(_jets_phi.at(i)-_z_phi);
    }

    if (_option=="useMetShift") {
      chi2 += pow((new_met_para + par.at(_njets))/sigma_met_para, 2); 
    }
    else {
      chi2 += pow(new_met_para/sigma_met_para, 2); 
    }
    chi2 += pow(new_met_perp/sigma_met_perp, 2); 


    return chi2;
 
  }

  virtual double Up() const {return 1.0; }

 // double GetMetPara() const { return _metPara; }
 // std::vector<double> GetJetsPara() const { return _jetsPara; }
 // std::vector<double> GetJetsSigma() const { return _jetsSigma; }

  
private:
  double _z_pt;
  double _z_phi;
  double _met_pt;
  double _met_phi;
  std::vector<double> _jets_pt;
  std::vector<double> _jets_phi;
  std::vector<double> _jets_sigma;
  
  int _njets;   
  std::string _option;

};


// some functions

int find_ut_jet(std::vector<int>& index, double z_pt, double z_phi, 
                int njet, float* jet_pt, float* jet_phi, int* jet_id, 
                float* jet_chargedEmEnergyFraction, float* jet_muonEnergyFraction) {
 
  int idx=-100;
  if (z_pt<10.0) return idx;
  double sum_jet_para = 0;
  for (int j=0; j<(int)index.size(); j++) {
    sum_jet_para += jet_pt[index.at(j)]*cos(jet_phi[index.at(j)]-z_phi);
  }
  double min_dpara=100000000000;
  for (int j=0; j<njet; j++){
    bool repeated = false;
    for (int k=0; k<(int)index.size(); k++){
      if (j==index.at(k)) repeated = true;
    }
    if (repeated) continue;
    bool jetLepVeto = (jet_chargedEmEnergyFraction[j]<0.8&&jet_muonEnergyFraction[j]<0.8);
    bool jetTightId = (jet_id[j]>=3);
    if (!jetLepVeto) continue;
    if (!jetTightId) continue;
    double jet_para = jet_pt[j]*cos(jet_phi[j]-z_phi);
    if (jet_para>0) continue;
    double jet_perp = jet_pt[j]*sin(jet_phi[j]-z_phi);
    if (fabs(jet_perp)>3.*fabs(jet_para)) continue;
    double dpara = jet_para+sum_jet_para+z_pt;
    if (fabs(dpara)>fabs(sum_jet_para+z_pt)) continue;
    if (fabs(dpara) < min_dpara) {
      min_dpara = dpara;
      idx = j;
    }
  }

  return idx;

}


int find_ut_jet_v3(std::vector<int>& index, double z_pt, double z_phi,
                double met, double met_phi,
                int njet, float* jet_pt, float* jet_phi, int* jet_id,
                float* jet_chargedEmEnergyFraction, float* jet_muonEnergyFraction) {

  int idx=-100;
  if (z_pt<10.0) return idx;
  double sum_jet_para = 0;
  for (int j=0; j<(int)index.size(); j++) {
    sum_jet_para += jet_pt[index.at(j)]*cos(jet_phi[index.at(j)]-z_phi);
  }
  double min_dpara=100000000000;
  for (int j=0; j<njet; j++){
    bool repeated = false;
    for (int k=0; k<(int)index.size(); k++){
      if (j==index.at(k)) repeated = true;
    }
    if (repeated) continue;
    //if (idx==j) continue; 
    bool jetLepVeto = (jet_chargedEmEnergyFraction[j]<0.8&&jet_muonEnergyFraction[j]<0.8);
    bool jetTightId = (jet_id[j]>=3);
    if (!jetLepVeto) continue;
    if (!jetTightId) continue;
    double jet_para = (double)jet_pt[j]*cos(jet_phi[j]-z_phi);
    if (jet_para>0) continue;
    double jet_perp = (double)jet_pt[j]*sin(jet_phi[j]-z_phi);
    if (fabs(jet_perp)>3.*fabs(jet_para)) continue;
    double met_para = met*cos(met_phi-z_phi);
    double dpara = jet_para+sum_jet_para+met_para+z_pt;
    if (fabs(dpara)>fabs(sum_jet_para+met_para+z_pt)) continue;
    if (fabs(dpara) < min_dpara) {
      min_dpara = fabs(dpara);
      idx = j;
    }
  }
  return idx;
}

int find_ut_jet_v4(std::vector<int>& index, double z_pt, double z_phi,
                double met, double met_phi,
                int njet, float* jet_pt, float* jet_phi, int* jet_id,
                float* jet_chargedEmEnergyFraction, float* jet_muonEnergyFraction) {

  int idx=-100;
  if (z_pt<10.0) return idx;
  double sum_jet_para = 0;
  int n_idx=(int)index.size();
  for (int j=0; j<n_idx; j++) {
    sum_jet_para += jet_pt[index.at(j)]*cos(jet_phi[index.at(j)]-z_phi);
  }
  double min_dpara=100000000000;
  for (int j=0; j<njet; j++){
    bool repeated = false;
    for (int k=0; k<n_idx; k++){
      if (j==index.at(k)) repeated = true;
    }
    if (repeated) continue;
    //if (idx==j) continue; 
    bool jetLepVeto = (jet_chargedEmEnergyFraction[j]<0.8&&jet_muonEnergyFraction[j]<0.8);
    bool jetTightId = (jet_id[j]>=3);
    if (!jetLepVeto) continue;
    if (!jetTightId) continue;
    double jet_para = (double)jet_pt[j]*cos(jet_phi[j]-z_phi);
    //if (jet_para>0) continue;
    double jet_perp = (double)jet_pt[j]*sin(jet_phi[j]-z_phi);
    //if (fabs(jet_perp)>3.*fabs(jet_para)) continue;
    double met_para = met*cos(met_phi-z_phi);
    double dpara = jet_para+sum_jet_para+met_para+z_pt;
    if (fabs(dpara)>fabs(sum_jet_para+met_para+z_pt)) continue;
    if (fabs(dpara) < min_dpara) {
      min_dpara = fabs(dpara);
      idx = j;
    }
  }
  return idx;
}


std::vector<int> find_ut_jet_v5(double z_pt, double z_phi,
                double met_pt, double met_phi,
                int njet, float* jet_pt, float* jet_phi, int* jet_id,
                float* jet_chargedEmEnergyFraction, float* jet_muonEnergyFraction) {

  std::vector<int> index;
  std::vector<int> init_index;
  std::vector<double> jets_para;
  std::vector<double> jets_perp;
  for (int j=0; j<njet; j++) {
    bool jetLepVeto = (jet_chargedEmEnergyFraction[j]<0.8&&jet_muonEnergyFraction[j]<0.8);
    bool jetTightId = (jet_id[j]>=3);
    if (!jetTightId) continue;
    double jet_para = jet_pt[j]*cos(jet_phi[j]-z_phi);
    double jet_perp = jet_pt[j]*sin(jet_phi[j]-z_phi);
    if (!jetLepVeto) continue;
    //if (jetLepVeto) z_pt -= jet_para;
    init_index.push_back(j);
    jets_para.push_back(jet_para);
    jets_perp.push_back(jet_perp);
  }    

  // sort jets_para from small to big
  int njets_sel = (int)jets_para.size();

  // return empty index vector if no jets pass id
  if (njets_sel==0) return index;

  if (njets_sel>1) {
    for (int i=0; i<njets_sel; i++){
      for (int j=i+1; j<njets_sel; j++){
        if (jets_para.at(i)>jets_para.at(j)) {
          double coucou = jets_para.at(i);
          jets_para.at(i) = jets_para.at(j);
          jets_para.at(j) = coucou;
          coucou = jets_perp.at(i);
          jets_perp.at(i) = jets_perp.at(j);
          jets_perp.at(j) = coucou;
          int icoucou = init_index.at(i);
          init_index.at(i) = init_index.at(j);
          init_index.at(j) = icoucou;
        }
      }
    }
  }
  // selecting the group of jets satisfy
  //  sum_jets_para + met_para + z_pt => 0
  std::vector<double> sum_paras;
  sum_paras.push_back(0);
  double sum_para_i = 0;
  double min_dpara = 1e30;
  int idx_first = 0;
  int idx_last = 0;
  double met_para = met_pt*cos(met_phi-z_phi);
  double target = -(z_pt+met_para);
  if (debug) {
    std::cout << "### jet selection v5, njets pass id = " << njets_sel << "###" << std::endl;
    std::cout << " z_pt + met_para = " << z_pt+met_para << std::endl;
  }

  for(int i = 0; i < (int)jets_para.size(); i++){
    sum_para_i += jets_para.at(i);
    double dpara = sum_para_i - target;
    for (int j=0; j<(int)sum_paras.size(); j++){
      double sum_para_j = sum_paras.at(j);
      double diff =  sum_para_i - sum_para_j - target;
      if (fabs(diff) < fabs(min_dpara)){
        min_dpara = diff;
        idx_last = i;
        idx_first = j;
      }
    }
    if (debug) {
      std::cout << " jet " << i << " : jet_para = " << jets_para.at(i) << ", sum_para_i = " << sum_para_i << ", dpara_i = " << dpara
              << ", idx_last = " << idx_last << ", idx_first = " << idx_first << ", min_dpara = " << min_dpara << std::endl;
    }
    sum_paras.push_back(sum_para_i);
  }

  for (int i=idx_first; i<=idx_last; i++) {
    index.push_back(init_index.at(i));  
  }
  if (debug) {
    std::cout << " min_dpara = " << min_dpara << ", min_sum_para_i = " << sum_paras.at(idx_last) << ", min_sum_para_j = " << sum_paras.at(idx_first) << ", "
            << " z_pt + met_para + min_dpara = " << -target+ sum_paras.at(idx_first)-sum_paras.at(idx_last)
            << std::endl;
  }
  return index;
}


std::vector<int> find_ut_jet_v6(double z_pt, double z_phi,
                double met_pt, double met_phi,
                int njet, float* jet_pt, float* jet_phi, int* jet_id,
                float* jet_chargedEmEnergyFraction, float* jet_muonEnergyFraction) {
  // remove lepton veto
  std::vector<int> index;
  std::vector<int> jets_index;
  std::vector<int> lepton_index; 
  for (int j=0; j<njet; j++) {
    bool jetLepVeto = (jet_chargedEmEnergyFraction[j]<0.2&&jet_muonEnergyFraction[j]<0.2);
    bool jetTightId = (jet_id[j]>=3);
    // reserve vetoed leptons 
    if (!jetLepVeto) lepton_index.push_back(j);
    if (!jetTightId) continue;
    if (!jetLepVeto) continue;
    jets_index.push_back(j);
  }    
 
  // keep only two leading  lepton jets
  if (lepton_index.size()>=2) {
    // sort according to pt
    for (int i=0; i<(int)lepton_index.size(); i++){
      for (int j=i+1; j<(int)lepton_index.size(); j++) {
        int idx_i = lepton_index.at(i);
        int idx_j = lepton_index.at(j);
        if (jet_pt[idx_i]<jet_pt[idx_j]) {
          lepton_index.at(i) = idx_j;
          lepton_index.at(j) = idx_i;
        }
      }   
    }
    // add the reststo jets, and erase them 
    for (int i=(int)lepton_index.size()-1; i>=2; i--){
      jets_index.push_back(lepton_index.at(i));
      lepton_index.pop_back();
    }
    // then sort again from small to big according to jet_para
    int idx_i = lepton_index.at(0);
    int idx_j = lepton_index.at(1);
    double jet_para_i = jet_pt[idx_i]*cos(jet_phi[idx_i]-z_phi);
    double jet_para_j = jet_pt[idx_j]*cos(jet_phi[idx_j]-z_phi);
    if (jet_para_i>jet_para_j) {
      lepton_index.at(0) = idx_j;
      lepton_index.at(1) = idx_i;
    }
    
  }

  // sort jets according to para, from negative to positive
  int njets_sel = (int)jets_index.size();
  if (njets_sel>1) {
    for (int i=0; i<njets_sel; i++){
      for (int j=i+1; j<njets_sel; j++){
        int idx_i = jets_index.at(i);
        int idx_j = jets_index.at(j);
        double jet_para_i = jet_pt[idx_i]*cos(jet_phi[idx_i]-z_phi); 
        double jet_para_j = jet_pt[idx_j]*cos(jet_phi[idx_j]-z_phi); 
        if (jet_para_i>jet_para_j) {
          jets_index.at(i) = idx_j;
          jets_index.at(j) = idx_i;
        }
      }
    }
  }

  // selecting the group of jets satisfy
  //  sum_jets_para + met_para + z_pt => 0
  std::vector<double> sum_paras;
  sum_paras.push_back(0);
  double sum_para_i = 0;
  double min_dpara = 1e30;
  int idx_first = 0;
  int idx_last = 0;
  double met_para = 0.0;// met_pt*cos(met_phi-z_phi);
  double target = -(z_pt+met_para);
  if (debug) {
    std::cout << "### jet selection v6, njets pass id = " << njets_sel << "###" << std::endl;
    std::cout << " z_pt + met_para = " << z_pt+met_para << std::endl;
  }
  
  for(int i = 0; i < (int)jets_index.size(); i++){
    int idx_i = jets_index.at(i);
    double jet_para_i = jet_pt[idx_i]*cos(jet_phi[idx_i]-z_phi);
    sum_para_i += jet_para_i;
    double dpara = sum_para_i - target;
    for (int j=0; j<(int)sum_paras.size(); j++){
      double sum_para_j = sum_paras.at(j);
      double diff =  sum_para_i - sum_para_j - target;
      if (fabs(diff) < fabs(min_dpara)){
        min_dpara = diff;
        idx_last = i;
        idx_first = j;
      }
    }
    if (debug) {
      std::cout << " jet " << i << " : jet_para = " << jet_para_i << ", sum_para_i = " << sum_para_i << ", dpara_i = " << dpara
              << ", idx_last = " << idx_last << ", idx_first = " << idx_first << ", min_dpara = " << min_dpara << std::endl;
    }
    sum_paras.push_back(sum_para_i);
  }

  if (debug) {
    std::cout << " min_dpara = " << min_dpara << ", min_sum_para_i = " << sum_paras.at(idx_last) << ", min_sum_para_j = " << sum_paras.at(idx_first) << ", "
            << " z_pt + met_para + min_dpara = " << -target+ sum_paras.at(idx_first)-sum_paras.at(idx_last)
            << std::endl;
  }
  // put selected jets to index vector
  for (int i=idx_first; i<=idx_last; i++) {
    if ((int)jets_index.size()<=idx_last) break;
    if ((int)jets_index.size()<=idx_first) break;
    index.push_back(jets_index.at(i));  
  }
  // also put the lepton jets in index
  //for (int i=0; i<(int)lepton_index.size(); i++){
  //  index.push_back(lepton_index.at(i));
  //}

  return index;
}


std::vector<int> find_ut_jet_v7(double z_pt, double z_phi,
                double met_pt, double met_phi,
                int njet, float* jet_pt, float* jet_phi, int* jet_id,
                float* jet_chargedEmEnergyFraction, float* jet_muonEnergyFraction) {
  // take all tigh id jets except leptons
  std::vector<int> index;
  std::vector<int> jets_index;
  std::vector<int> lepton_index; 
  for (int j=0; j<njet; j++) {
    bool jetLepVeto = (jet_chargedEmEnergyFraction[j]<0.2&&jet_muonEnergyFraction[j]<0.2);
    bool jetTightId = (jet_id[j]>=3);
    // reserve vetoed leptons 
    if (!jetLepVeto) lepton_index.push_back(j);
    if (!jetTightId) continue;
    if (!jetLepVeto) continue;
    jets_index.push_back(j);
  }    
 
  // keep only two leading  lepton jets
  if (lepton_index.size()>=2) {
    // sort according to pt
    for (int i=0; i<(int)lepton_index.size(); i++){
      for (int j=i+1; j<(int)lepton_index.size(); j++) {
        int idx_i = lepton_index.at(i);
        int idx_j = lepton_index.at(j);
        if (jet_pt[idx_i]<jet_pt[idx_j]) {
          lepton_index.at(i) = idx_j;
          lepton_index.at(j) = idx_i;
        }
      }   
    }
    // add the reststo jets, and erase them 
    for (int i=(int)lepton_index.size()-1; i>=2; i--){
      jets_index.push_back(lepton_index.at(i));
      lepton_index.pop_back();
    }
    // then sort again from small to big according to jet_para
    int idx_i = lepton_index.at(0);
    int idx_j = lepton_index.at(1);
    double jet_para_i = jet_pt[idx_i]*cos(jet_phi[idx_i]-z_phi);
    double jet_para_j = jet_pt[idx_j]*cos(jet_phi[idx_j]-z_phi);
    if (jet_para_i>jet_para_j) {
      lepton_index.at(0) = idx_j;
      lepton_index.at(1) = idx_i;
    }
    
  }

  // sort jets according to para, from negative to positive
  int njets_sel = (int)jets_index.size();
  if (njets_sel>1) {
    for (int i=0; i<njets_sel; i++){
      for (int j=i+1; j<njets_sel; j++){
        int idx_i = jets_index.at(i);
        int idx_j = jets_index.at(j);
        double jet_para_i = jet_pt[idx_i]*cos(jet_phi[idx_i]-z_phi); 
        double jet_para_j = jet_pt[idx_j]*cos(jet_phi[idx_j]-z_phi); 
        if (jet_para_i>jet_para_j) {
          jets_index.at(i) = idx_j;
          jets_index.at(j) = idx_i;
        }
      }
    }
  }

  index = jets_index;
  return index;
}



std::vector<int> find_ut_jet_v8(double z_pt, double z_phi,
                double met_pt, double met_phi,
                int njet, float* jet_pt, float* jet_phi, int* jet_id,
                float* jet_chargedEmEnergyFraction, float* jet_muonEnergyFraction) {
  // take all tigh id jets including lepton jets, will be splitted to lepton and jet later
  std::vector<int> index;
  std::vector<int> jets_index;
  std::vector<int> lepton_index; 
  for (int j=0; j<njet; j++) {
    bool jetTightId = (jet_id[j]>=3);
    bool jetLooseId = (jet_id[j]>=1);
    if (!jetTightId) continue;
    //if (!jetLooseId) continue;
    jets_index.push_back(j);
  }    
 
  // sort jets according to para, from negative to positive
  int njets_sel = (int)jets_index.size();
  if (njets_sel>1) {
    for (int i=0; i<njets_sel; i++){
      for (int j=i+1; j<njets_sel; j++){
        int idx_i = jets_index.at(i);
        int idx_j = jets_index.at(j);
        double jet_para_i = jet_pt[idx_i]*cos(jet_phi[idx_i]-z_phi); 
        double jet_para_j = jet_pt[idx_j]*cos(jet_phi[idx_j]-z_phi); 
        if (jet_para_i>jet_para_j) {
          jets_index.at(i) = idx_j;
          jets_index.at(j) = idx_i;
        }
      }
    }
  }

  index = jets_index;
  return index;
}







