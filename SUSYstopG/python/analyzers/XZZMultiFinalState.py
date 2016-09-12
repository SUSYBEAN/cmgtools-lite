from CMGTools.XZZ2l2nu.analyzers.XZZEventInterpretationBase import *
from CMGTools.XZZ2l2nu.tools.Pair import Pair
import itertools
import ROOT
from array import array

class XZZMultiFinalState( XZZEventInterpretationBase ):
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(XZZMultiFinalState,self).__init__(cfg_ana, cfg_comp, looperName)
#        if "LLNuNu" in self.processTypes:
#            self.eeff=ROOT.TF1('trigeff', '50*[0]*(1+TMath::Erf((x-[1])/(sqrt(2)*[2])))*(1. + TMath::Erf((x-[3])/(sqrt(2)*[4])))')
#            self.ep1=array('d',[4.88895e-01,1.20166e+02,1.17300e+01,1.52207e+02,3.01742e+01])
#            self.ep2=array('d',[4.91171e-01,1.15409e+02,1.03195e+01,1.50434e+02,2.85255e+01])
#            self.ep1e=array('d',[1.41647e-03,2.48648,1.55638,9.45006e-01,1.69037])
#            self.ep2e=array('d',[6.27821e-04,8.85865e-01,6.79409e-01,3.08905e-01,5.33998e-01])
#            self.etrgsf=lambda x:self.eeff.EvalPar(array('d',[x]),self.ep1)/self.eeff.EvalPar(array('d',[x]),self.ep2) if self.eeff.EvalPar(array('d',[x]),self.ep2) else 1
#            #msff=ROOT.gROOT.GetFunction("pol5")
#            #mp=array('d',[0.98395,-4.05067e-05,1.74509e-06,-1.1423e-08,2.67736e-11,-2.00434e-14])
#            self.mtrgsf=lambda x:0.983879
#            self.mtrgsfe=lambda x:0.00068

#    def gettrigersf(self,zpt,pdgid):
#        if abs(pdgid)==11:return self.etrgsf(zpt)
#        if abs(pdgid)==13:return self.mtrgsf(zpt)

#    def gettrigersfUpper(self,zpt,pdgid):
#        if abs(pdgid)==11:
#            ep1=array('d',[self.ep1[0]+self.ep1e[0],self.ep1[1]-self.ep1e[1],self.ep1[2]-self.ep1e[2],self.ep1[3]-self.ep1e[3],self.ep1[4]-self.ep1e[4]])
#            ep2=array('d',[self.ep2[0]-self.ep2e[0],self.ep2[1]+self.ep2e[1],self.ep2[2]+self.ep2e[2],self.ep2[3]+self.ep2e[3],self.ep2[4]+self.ep2e[4]])
#            if zpt<self.ep1[1]:ep1[2]=self.ep1[2]+self.ep1e[2]
#            if zpt<self.ep1[3]:ep1[4]=self.ep1[4]+self.ep1e[4]
#            if zpt<self.ep2[1]:ep2[2]=self.ep2[2]-self.ep2e[2]
#            if zpt<self.ep2[3]:ep2[4]=self.ep2[4]-self.ep2e[4]
#            return self.eeff.EvalPar(array('d',[zpt]),ep1)/self.eeff.EvalPar(array('d',[zpt]),ep2) if self.eeff.EvalPar(array('d',[zpt]),ep2) else 1
#        if abs(pdgid)==13:return self.mtrgsf(zpt)+self.mtrgsfe(zpt)

#    def gettrigersfLower(self,zpt,pdgid):
#        if abs(pdgid)==11:
#            ep1=array('d',[self.ep1[0]-self.ep1e[0],self.ep1[1]+self.ep1e[1],self.ep1[2]+self.ep1e[2],self.ep1[3]+self.ep1e[3],self.ep1[4]+self.ep1e[4]])
#            ep2=array('d',[self.ep2[0]+self.ep2e[0],self.ep2[1]-self.ep2e[1],self.ep2[2]-self.ep2e[2],self.ep2[3]-self.ep2e[3],self.ep2[4]-self.ep2e[4]])
#            if zpt<self.ep1[1]:ep1[2]=self.ep1[2]-self.ep1e[2]
#            if zpt<self.ep1[3]:ep1[4]=self.ep1[4]-self.ep1e[4]
#            if zpt<self.ep2[1]:ep2[2]=self.ep2[2]+self.ep2e[2]
#            if zpt<self.ep2[3]:ep2[4]=self.ep2[4]+self.ep2e[4]
#            return self.eeff.EvalPar(array('d',[zpt]),ep1)/self.eeff.EvalPar(array('d',[zpt]),ep2) if self.eeff.EvalPar(array('d',[zpt]),ep2) else 1
#        if abs(pdgid)==13:return self.mtrgsf(zpt)-self.mtrgsfe(zpt)

    def process(self, event):
        super(XZZMultiFinalState,self).process(event)
        self.counters.counter('events').inc('all events')

        LLNuNu=[]
        ElMuNuNu=[]
        PhotonJets=[]

        # do LL+MET combination
        if ("LLNuNu" in self.processTypes) and hasattr(event, "LL") and len(event.LL)>0:
            # Take the Z->ll nearest to the Z mass and the highest pt jets
            bestZ = min(event.LL,key = lambda x: abs(x.M()-91.1876))
            VV=Pair(bestZ,event.met)
            if self.selectPairLLNuNu(VV):
                selected = {'pair':VV}
                LLNuNu.append(selected)
#                if self.cfg_comp.isMC:
#                    event.trgsf = self.gettrigersf(bestZ.pt(),bestZ.leg1.pdgId())
#                    event.trgsfUp = self.gettrigersfUpper(bestZ.pt(),bestZ.leg1.pdgId())
#                    event.trgsfLo = self.gettrigersfLower(bestZ.pt(),bestZ.leg1.pdgId())


        if ("LLNuNu" in self.processTypes) and len(LLNuNu)>0: 
            self.counters.counter('events').inc('pass llNuNu events')
            setattr(event,'LLNuNu'+self.cfg_ana.suffix,LLNuNu)

        # do ElMu+MET combination
        if ("ElMuNuNu" in self.processTypes) and hasattr(event, 'ElMu') and len(event.ElMu)>0:
            # Take the elmu nearest to the Z mass and the highest pt
            bestPair = min(event.ElMu,key = lambda x: abs(x.M()-91.1876))
            fakeVV=Pair(bestPair,event.met)
            if self.selectPairElMuNuNu(fakeVV):
                selected = {'pair':fakeVV}
                ElMuNuNu.append(selected)
                # if self.cfg_comp.isMC:
                #     event.trgsf = self.gettrigersf(bestPair.pt(),bestPair.leg1.pdgId())
                #     event.trgsfUp = self.gettrigersfUpper(bestPair.pt(),bestPair.leg1.pdgId())
                #     event.trgsfLo = self.gettrigersfLower(bestPair.pt(),bestPair.leg1.pdgId())
        if ("ElMuNuNu" in self.processTypes) and len(ElMuNuNu)>0 :
            self.counters.counter('events').inc('pass elmuNuNu events')
            setattr(event,'ElMuNuNu'+self.cfg_ana.suffix,ElMuNuNu)

        #
        # do Photon+jets combination
        if ("PhotonJets" in self.processTypes) and hasattr(event, 'selectedPhotons') and len(event.selectedPhotons)>0 :
            for g in event.selectedPhotons: 
                gjet = Pair(g,event.met)
                if self.selectPhotonJets(gjet):
                    selected = {'pair':gjet}
                    PhotonJets.append(selected)
        if ("PhotonJets" in self.processTypes) and len(PhotonJets)>0:
            self.counters.counter('events').inc('pass PhotonJets events')
            setattr(event,'PhotonJets'+self.cfg_ana.suffix,PhotonJets)
            
        




