#!/usr/bin/env python

import ROOT
import os

### Please keep the keys same in tex_dic={}  and cuts={}
class SetCuts:
    def __init__(self):
         
        self.Tex_dic={'SR': 'Region A', 'CRb': 'Region B','CRc': 'Region C','CRd': 'Region D'}
        self.met_pt_in=raw_input("[info] 'SetCuts' -> please give the MET cut: (enter to use default MET>100GeV)\n")
        self.met_pt='100' if self.met_pt_in=='' else self.met_pt_in
        # define a cutflow for signal region
        self.cutflow=("abs(llnunu_l1_mass-91.1876)<20.0", "llnunu_l1_pt>100.0",
                      "met_pt>"+self.met_pt, "llnunu_deltaPhi>2.5", "llnunu_l1_deltaR<0.6")
        

    # Cuts used for alpha-method to estimate non-resonant bkgs
    def alphaCuts(self, isIn=True, isll=True):
        astr = "({0}&&{1}&&{2}&&{3}&&{4})"
        cuts = astr.format(*self.cutflow)
        if not isIn:
            cuts=ROOT.TString(cuts).ReplaceAll("91.1876)<20.0","91.1876)<=55.0&&abs(llnunu_l1_mass-91.1876)>=25.0")
        #if not isll: cuts=ROOT.TString(cuts).ReplaceAll("llnunu","elmununu")
        
        return cuts

    def abcdCuts(self,channel, whichRegion="", isPreSelect=False):
        if whichRegion=="": whichRegion=raw_input("Please choose a benchmarck Region (SR or VR): \n")
        if whichRegion=='SR':
            preSelection='(abs(llnunu_l1_mass-91.1876)<20.0&&llnunu_l1_pt>100.0&&met_pt>'+self.met_pt+')'
        elif whichRegion=='VR':
            preSelection='(abs(llnunu_l1_mass-91.1876)<20.0&&llnunu_l1_pt>100.0&&met_pt<'+self.met_pt+')'
        else:
            print "I do not understand your benchmark Region, should be either 'SR' for MET>"+self.met_pt+"GeV or 'VR' for MET<"+self.met_pt+"GeV\n"
            sys.exit(0)
            
        pdgID={'el':'11','mu':'13' }
        if isPreSelect:
            cuts='('+preSelection+'&&abs(llnunu_l1_l1_pdgId)=='+pdgID[channel]+')'
        else:
            dR_ll={'el':'0.6','mu':'0.6'}
            dPhi_ZMET='2.5'
            var1='llnunu_deltaPhi'
            var2='llnunu_l1_deltaR'
            
            cuts_a='('+preSelection+'&&abs(llnunu_l1_l1_pdgId)=='+pdgID[channel]+'&&'+var1+'>'+dPhi_ZMET+'&&'+var2+'<'+dR_ll[channel]+')'
            cuts_c='('+preSelection+'&&abs(llnunu_l1_l1_pdgId)=='+pdgID[channel]+'&&'+var1+'>'+dPhi_ZMET+'&&'+var2+'>'+dR_ll[channel]+')'
            cuts_b='('+preSelection+'&&abs(llnunu_l1_l1_pdgId)=='+pdgID[channel]+'&&'+var1+'<'+dPhi_ZMET+'&&'+var2+'<'+dR_ll[channel]+')'
            cuts_d='('+preSelection+'&&abs(llnunu_l1_l1_pdgId)=='+pdgID[channel]+'&&'+var1+'<'+dPhi_ZMET+'&&'+var2+'>'+dR_ll[channel]+')'
            cuts={'SR':cuts_a, 'CRb':cuts_b, 'CRc':cuts_c, 'CRd':cuts_d}
        
        #cuts_test='1'
        return cuts
   
