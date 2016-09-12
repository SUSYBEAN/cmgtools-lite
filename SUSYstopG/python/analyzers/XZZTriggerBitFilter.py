import ROOT

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.statistics.counter import Counter, Counters        

class XZZTriggerBitFilter( Analyzer ):
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(XZZTriggerBitFilter,self).__init__(cfg_ana,cfg_comp,looperName)
        self.processName = getattr(self.cfg_ana,"processName","HLT")
        triggers = cfg_comp.triggers
        self.autoAccept = True if len(triggers) == 0 else False
        vetoTriggers = cfg_comp.vetoTriggers if hasattr(cfg_comp, 'vetoTriggers') else []
        trigVec = ROOT.vector(ROOT.string)()
        for t in triggers: trigVec.push_back(t)
        self.mainFilter = ROOT.heppy.TriggerBitChecker(trigVec)
        if len(vetoTriggers):
            vetoVec = ROOT.vector(ROOT.string)()
            for t in vetoTriggers: vetoVec.push_back(t)
            self.vetoFilter = ROOT.heppy.TriggerBitChecker(vetoVec)
        else:
            self.vetoFilter = None 
        
    def declareHandles(self):
        super(XZZTriggerBitFilter, self).declareHandles()
        self.handles['TriggerResults'] = AutoHandle( ('TriggerResults','',self.processName), 'edm::TriggerResults' )

    def beginLoop(self, setup):
        super(XZZTriggerBitFilter,self).beginLoop(setup)

        self.counters.addCounter('TriggerReport')
        self.count = self.counters.counter('TriggerReport')
        self.count.register('All Events')
        self.count.register('Pass Trigger Events')

    def process(self, event):
        if self.autoAccept: return True
        self.readCollections( event.input )
        self.count.inc('All Events')

        if not self.mainFilter.check(event.input.object(), self.handles['TriggerResults'].product()):
            return False
        if self.vetoFilter != None and self.vetoFilter.check(event.input.object(), self.handles['TriggerResults'].product()):
            return False
       
        self.count.inc('Pass Trigger Events')
 
        return True

