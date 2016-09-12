import ROOT
from  itertools import combinations
from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import *
from CMGTools.XZZ2l2nu.tools.Pair import *

class XZZMuonEffTree( Analyzer ):

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(XZZMuonEffTree,self).__init__(cfg_ana,cfg_comp,looperName)
        self.checktag=getattr(cfg_ana,'checktag',False)
        self.pfbkg=getattr(cfg_ana,'pfbkg',False)
        self.muHLT=getattr(cfg_ana,'muHLT',[])

    def declareHandles(self): 
        super(XZZMuonEffTree, self).declareHandles()
        if self.checktag:
            self.handles['trgresults']=AutoHandle(("TriggerResults","","HLT"),"edm::TriggerResults")
            self.handles['selectedtrg']=AutoHandle('selectedPatTrigger','std::vector<pat::TriggerObjectStandAlone>')
        if self.pfbkg:
            self.handles['pfcandidate']=AutoHandle("packedPFCandidates","std::vector<pat::PackedCandidate>")

    def istag(self,mu):
        itg=1
        if not mu.muonID("POG_ID_Tight"): itg-=1
#        if (mu.physObj.pfIsolationR04().sumChargedHadronPt + max(0., mu.physObj.pfIsolationR04().sumNeutralHadronEt + mu.physObj.pfIsolationR04().sumPhotonEt - 0.5*mu.physObj.pfIsolationR04().sumPUPt))/mu.pt()>.2: itg-=2
        if mu.absIsoWithFSR(0.4)/mu.pt()>.2: itg-=2
        return itg

    def trgbit(self,mu,names):
        mu.triggerbit=0
        leptrgobs=[i for i in self.handles['selectedtrg'].product() if deltaR(mu.eta(),mu.phi(),i.eta(),i.phi())<0.3]
        for itg in leptrgobs:
            itg.unpackPathNames(names)
            pNames=list(itg.pathNames(True))
            for n,i in enumerate(self.muHLT):
                if [pN for pN in pNames if i in pN]: mu.triggerbit|=1<<n

    def process(self, event):
        self.readCollections( event.input )
        event.llpair=[]
        if self.checktag:
            names = event.input.object().triggerNames(self.handles['trgresults'].product())
            for i in event.selectedMuons: 
                i.istag=self.istag(i)
                self.trgbit(i,names)
            
        if self.pfbkg:
            event.selectedhadrons=[i for i in self.handles['pfcandidate'].product() if i.pt()>20 and abs(i.pdgId())==211 and abs(i.eta())<2.4]
            if self.checktag: event.selectedMuonstagged=[i for i in event.selectedMuons if i.istag==1]
            if not (event.selectedMuons and event.selectedhadrons): return False
            for i in event.selectedMuons:event.llpair.extend([Pair(i,j,23) for j in event.selectedhadrons])
        else:
            event.zllcandidates=[i for i in event.selectedMuons if i.track().isNonnull()]
            if len(event.zllcandidates)<2: return False
            for l1,l2 in combinations(event.zllcandidates,2):
                if l1.pdgId() == -l2.pdgId():
                    pair = Pair(l1,l2,23)
                    event.llpair.append(pair)
            if not event.llpair: return False
            event.llpair=[min(event.llpair,key = lambda x: abs(x.M()-91.118))]
        return True




        
            

        


                
                
