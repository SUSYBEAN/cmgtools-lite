import ROOT
import os, types
from math import *
from PhysicsTools.HeppyCore.utils.deltar import *

class JetResolution:
    def __init__(self,globalTag="",jetFlavour="",jerPath="" ,debug=False):
        """Create resolution and scalefactor objects to read the payloads from the text dumps from global tag under
            CMGTools/XZZ2l2nu/data/jer  (provided from https://github.com/cms-jet/JRDatabase/tree/master/textFiles/Fall15_25nsV2_MC).
           It will jet PtResolution and resolution scale factors for jets.
           Receipt from 
           https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution
           """
        if globalTag=="" and jetFlavour=="" and jerPath=="": raise RuntimeError, 'Configuration not supported'
            
        self.globalTag = globalTag # "Summer15_25nsV6_DATA" or "Summer15_25nsV6_MC"
        self.jetFlavour = jetFlavour # "AK4PFchs" 
        self.jerPath = jerPath #"%s/src/CMGTools/XZZ2l2nu/data/jer" % os.environ['CMSSW_BASE'];
        path = os.path.expandvars(jerPath) 
        self.debug = debug
        if self.debug:
            print '[Debug] I am JetResolution.py: ',path
            print "[Debug] %s/%s_PtResolution_%s.txt" % (path,globalTag,jetFlavour)
            print "[Debug] %s/%s_SF_%s.txt" % (path,globalTag,jetFlavour)

        # Step1 : Construct a FactorizedJetResObject object
        self.JetResObject = ROOT.JME.JetResolutionObject("%s/%s_PtResolution_%s.txt" % (path,globalTag,jetFlavour))
        self.JetSFObject = ROOT.JME.JetResolutionObject("%s/%s_SF_%s.txt" % (path,globalTag,jetFlavour))

        if self.debug : 
            self.JetResObject.dump() 
            self.JetSFObject.dump()

        # Step2 : initialize JetParameters object for pt and eta (rho if it is sf)
        self.resPar = ROOT.JME.JetParameters()
        self.sfPar = ROOT.JME.JetParameters()
       

    def getResolution(self,jet,rho,resolution=None):
        if not resolution: resolution = self.JetResObject
        if resolution != self.JetResObject : raise RuntimeError, 'Configuration not supported'
        self.resPar.setJetPt(jet.pt())
        self.resPar.setJetEta(jet.eta())
        self.resPar.setRho(rho)

        res = 1.0
        if resolution.getRecord(self.resPar):
            res = resolution.evaluateFormula(resolution.getRecord(self.resPar), self.resPar)
        else:
            if self.debug: print '[Info] [JetResolution.py] no resolution record for [jet (pT, eta), rho] = [(%.2f, %.2f), %.2f]' % (jet.pt(), jet.eta(),rho)
        return res

    def getScaleFactor(self,jet,variation='norminal',scalefactor=None):

        Variation={'norminal':0,'down':1,'up':2} 

        if not scalefactor: scalefactor = self.JetSFObject
        if scalefactor != self.JetSFObject : raise RuntimeError, 'Configuration not supported'
        self.sfPar.setJetEta(jet.eta())
        res_sf = 1.0
        if scalefactor.getRecord(self.sfPar):
            res_sf = scalefactor.getRecord(self.sfPar).getParametersValues().at(Variation[variation])
        else:
            if self.debug: print '[Info] [JetResolution.py] no scale factor record for jet(pT, eta) = (%.2f, %.2f)' % (jet.pt(), jet.eta()) 
        return res_sf

