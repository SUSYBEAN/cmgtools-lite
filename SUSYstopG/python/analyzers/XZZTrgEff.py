from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from ROOT import TFile,TH1F
class TrgEffHists(object):
    def __init__(self, filename):
        self.file=TFile(filename,'recreate')
        self.dRmunoHLT=TH1F('dRmunoHLT','deltaR with mu no HLT',6,0,1.2)
        self.dRelnoHLT=TH1F('dRelnoHLT','deltaR with ele no HLT',6,0,1.2)
        self.dRmuHLT=TH1F('dRmuHLT','deltaR with muon HLT',6,0,1.2)
        self.dRmuHLTmatch=TH1F('dRmuHLTmatch','deltaR with muon HLT matching',6,0,1.2)
        self.dRelHLT=TH1F('dRelHLT','deltaR with ele HLT',6,0,1.2)
        self.dRelHLTmatch=TH1F('dRelHLTmatch','deltaR with ele HLT matching',6,0,1.2)
    def write(self):
        self.file.Write()

def deltR(l1,l2):
    pi=3.1415926536
    de=l1.eta()-l2.eta()
    dp=l1.phi()-l2.phi()
    while dp>pi: dp-=2*pi
    while dp<-pi: dp+=2*pi
    return (de**2+dp**2)**.5

class TrgObj(object):
    def __init__(self, pt=-100, eta=-100, phi=-100):
        self.pt=pt
        self.eta=eta
        self.phi=phi

class XZZTrgEff(Analyzer):
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(XZZTrgEff, self).__init__(cfg_ana, cfg_comp, looperName)
        self.trgeff=TrgEffHists('/'.join([self.dirName,'trgeff.root']))
        self.eleHLT=cfg_ana.eleHLT if hasattr(cfg_ana,'eleHLT') else ''
        self.muHLT=cfg_ana.muHLT if hasattr(cfg_ana,'muHLT') else ''
    def declareHandles(self):
        super(XZZTrgEff, self).declareHandles()
        self.handles['trgresults']=AutoHandle(("TriggerResults","","HLT"),"edm::TriggerResults")
        self.handles['selectedtrg']=AutoHandle('selectedPatTrigger','std::vector<pat::TriggerObjectStandAlone>')
    def beginLoop(self,setup):
        super(XZZTrgEff,self).beginLoop(setup)
        self.counters.addCounter('events')
        count = self.counters.counter('events')
        count.register('events with mu no HLT')
        count.register('events with ele no HLT')
        count.register('events pass mu HLT')
        count.register('events pass mu HLT and matching')
        count.register('events pass ele HLT')
        count.register('events pass ele HLT and matching')
    def process(self,event):
        self.readCollections( event.input )
        Zll = event.LLNuNu[0]['pair'].leg1
        deltr = deltR(Zll.leg1,Zll.leg2)
        names = event.input.object().triggerNames(self.handles['trgresults'].product())
        eleob=[]
        muob=[]
        Zll.leg1.triggerob=TrgObj()
        Zll.leg2.triggerob=TrgObj()
        for i in self.handles['selectedtrg'].product():
            i.unpackPathNames(names)
            pNames=list(i.pathNames())
            if self.eleHLT and [pN for pN in pNames if self.eleHLT in pN]:eleob.append(i)
            if self.muHLT and [pN for pN in pNames if self.muHLT in pN]:muob.append(i)
        if abs(Zll.leg1.pdgId())==11 and self.eleHLT:
            self.trgeff.dRelnoHLT.Fill(deltr)
            self.counters.counter('events').inc('events with ele no HLT')
            if eleob:
                self.trgeff.dRelHLT.Fill(deltr)
                self.counters.counter('events').inc('events pass ele HLT')
                lmatch1=[i for i in eleob if deltR(Zll.leg1,i)<.3 and deltR(Zll.leg1,i)<=deltR(Zll.leg2,i)]
                lmatch2=[i for i in eleob if deltR(Zll.leg2,i)<.3 and deltR(Zll.leg2,i)<deltR(Zll.leg1,i)]
                if lmatch1:
                    m1=min(lmatch1,key=lambda x:deltR(Zll.leg1,x))
                    Zll.leg1.triggerob=TrgObj(m1.pt(),m1.eta(),m1.phi())
                if lmatch2:
                    m2=min(lmatch2,key=lambda x:deltR(Zll.leg2,x))
                    Zll.leg2.triggerob=TrgObj(m2.pt(),m2.eta(),m2.phi())
                if lmatch1+lmatch2:
                    self.trgeff.dRelHLTmatch.Fill(deltr)
                    self.counters.counter('events').inc('events pass ele HLT and matching')
        if abs(Zll.leg1.pdgId())==13 and self.muHLT:
            self.trgeff.dRmunoHLT.Fill(deltr)
            self.counters.counter('events').inc('events with mu no HLT')
            if muob:
                self.trgeff.dRmuHLT.Fill(deltr)
                self.counters.counter('events').inc('events pass mu HLT')
                lmatch1=[i for i in muob if deltR(Zll.leg1,i)<.3 and deltR(Zll.leg1,i)<=deltR(Zll.leg2,i)]
                lmatch2=[i for i in muob if deltR(Zll.leg2,i)<.3 and deltR(Zll.leg2,i)<deltR(Zll.leg1,i)]
                if lmatch1:
                    m1=min(lmatch1,key=lambda x:deltR(Zll.leg1,x))
                    Zll.leg1.triggerob=TrgObj(m1.pt(),m1.eta(),m1.phi())
                if lmatch2:
                    m2=min(lmatch2,key=lambda x:deltR(Zll.leg2,x))
                    Zll.leg2.triggerob=TrgObj(m2.pt(),m2.eta(),m2.phi())
                if lmatch1+lmatch2:
                    self.trgeff.dRmuHLTmatch.Fill(deltr)
                    self.counters.counter('events').inc('events pass mu HLT and matching')
        return True
    def write(self,setup):
        super(XZZTrgEff, self).write(setup)
        self.trgeff.write()

