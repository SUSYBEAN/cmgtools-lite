import ROOT
import random
import math
from  itertools import combinations
from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import *
import PhysicsTools.HeppyCore.framework.config as cfg
from CMGTools.XZZ2l2nu.tools.Pair import *


class XZZLeptonEffTree( Analyzer ):

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(XZZLeptonEffTree,self).__init__(cfg_ana,cfg_comp,looperName)
        self.genfilter=getattr(self.cfg_ana,"genfilter",False)
        self.eithercharge=getattr(self.cfg_ana,"eithercharge",False)

    def checkgen(self,lep):
        try:
            mom=lep.physObj.genParticle()
            while mom.pdgId()!=23:
                mom=mom.mother()
            return True
        except:
            return False
        
    def process(self, event):
        self.readCollections( event.input )
        event.llpair=[]
        if len(event.selectedElectrons)>2 and (event.selectedElectrons[0].pdgId()==-event.selectedElectrons[1].pdgId() or self.eithercharge):
            event.llpair.append(Pair(event.selectedElectrons[0],event.selectedElectrons[1],23))
        if len(event.selectedMuons)>2 and (event.selectedMuons[0].pdgId()==-event.selectedMuons[1].pdgId() or self.eithercharge):
            event.llpair.append(Pair(event.selectedMuons[0],event.selectedMuons[1],23))
        if not event.llpair: return False
        if self.genfilter and self.cfg_comp.isMC:
            for i in event.llpair:
                i.leg1.xdaughter=self.checkgen(i.leg1)
                i.leg2.xdaughter=self.checkgen(i.leg2)
        return True




        
            

        


                
                
