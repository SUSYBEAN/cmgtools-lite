from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import *

class XZZMultTrgEff(Analyzer):
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(XZZMultTrgEff, self).__init__(cfg_ana, cfg_comp, looperName)
        self.HLTlist=cfg_ana.HLTlist if hasattr(cfg_ana,'HLTlist') else []

    def declareHandles(self):
        super(XZZMultTrgEff, self).declareHandles()
        self.handles['trgresults']=AutoHandle(("TriggerResults","","HLT"),"edm::TriggerResults")
        self.handles['selectedtrg']=AutoHandle('selectedPatTrigger','std::vector<pat::TriggerObjectStandAlone>')

    def process(self,event):
        self.readCollections( event.input )
        Zll = event.LLNuNu[0]['pair'].leg1
        names = event.input.object().triggerNames(self.handles['trgresults'].product())
        trgoblist=self.handles['selectedtrg'].product()
        Zll.leg1.triggerbit=0
        Zll.leg2.triggerbit=0
        leptrgobs1=[i for i in trgoblist if deltaR(Zll.leg1.eta(),Zll.leg1.phi(),i.eta(),i.phi())<min(.3,deltaR(Zll.leg2.eta(),Zll.leg2.phi(),i.eta(),i.phi()))]
        leptrgobs2=[i for i in trgoblist if deltaR(Zll.leg2.eta(),Zll.leg2.phi(),i.eta(),i.phi())<min(.3,deltaR(Zll.leg1.eta(),Zll.leg1.phi(),i.eta(),i.phi()))]
        for itg in leptrgobs1:
            itg.unpackPathNames(names)
            pNames=list(itg.pathNames(True))
            for n,i in enumerate(self.HLTlist):
                if [pN for pN in pNames if i in pN]: Zll.leg1.triggerbit|=1<<n
        for itg in leptrgobs2:
            itg.unpackPathNames(names)
            pNames=list(itg.pathNames(True))
            for n,i in enumerate(self.HLTlist):
                if [pN for pN in pNames if i in pN]: Zll.leg2.triggerbit|=1<<n
        return True


