from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.statistics.counter import Counter, Counters
import PhysicsTools.HeppyCore.framework.config as cfg
from DataFormats.FWLite import Handle
from ROOT.gen import WeightsInfo

class XZZLHEWeightAnalyzer( Analyzer ):
    """Read the WeightsInfo objects of the LHE branch and store them
       in event.LHE_weights list.

       If the WeightsInfo.id is a string, replace it with an integer.

       So far the only allowed string format is "mg_reweight_X",
       which gets stored as str(10000+int(X))

       If w.id is an unknown string or anything but a string or int,
       a RuntimeError is raised.
    """
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(XZZLHEWeightAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)

    def declareHandles(self):
        super(XZZLHEWeightAnalyzer, self).declareHandles()
        self.mchandles['GenInfo'] = AutoHandle( ('generator','',''), 'GenEventInfoProduct' )
        self.mchandles['LHEweights'] = AutoHandle('externalLHEProducer',
                                                  'LHEEventProduct',
                                                  mayFail=True,
                                                  fallbackLabel='source',
                                                  lazy=False )

    def beginLoop(self, setup):
        super(XZZLHEWeightAnalyzer,self).beginLoop(setup)
        self.counters.addCounter('LHEWeightReport')
        self.count = self.counters.counter('LHEWeightReport')
        self.count.register('All Events')
        if self.cfg_comp.isMC:
            self.count.register('Sum Weights')
            self.count.register('SumLHEOrigWeights')

    def process(self, event):
        self.readCollections( event.input )

        # if not MC, nothing to do
        if not self.cfg_comp.isMC:
            return True

        self.count.inc('All Events')
        self.count.inc('Sum Weights', self.mchandles['GenInfo'].product().weight())

        # Add LHE weight info
        event.LHE_weights = []
        event.LHE_originalWeight = 1.0

        if self.mchandles['LHEweights'].isValid():

            event.LHE_originalWeight = self.mchandles['LHEweights'].product().originalXWGTUP()
            self.count.inc('SumLHEOrigWeights', self.mchandles['LHEweights'].product().originalXWGTUP())

            for w in self.mchandles['LHEweights'].product().weights():
                # Check if id is string or int and convert to int if it's a string
                try:
                    int(w.id)
                    event.LHE_weights.append(w)
                    # register counter if not exist
                    if not ( str(w.id) in self.count.dico.keys() ):
                        self.count.register(str(w.id))
                    self.count.inc(str(w.id), w.wgt)
                       
                except ValueError:
                    if not type(w.id) == str:
                        raise RuntimeError('Non int or string type for LHE weight id')

                    newweight = WeightsInfo()
                    newweight.wgt = w.wgt
                    if w.id.startswith('mg_reweight'):
                        newid = str(10000 + int(w.id.rsplit('_',1)[1]))
                        newweight.id = newid

                    else: raise RuntimeError('Unknown string id in LHE weights')
                    event.LHE_weights.append(newweight)

                    # register counter if not exist
                    if not ( str(newid) in self.count.dico.keys() ):
                        self.count.register(str(w.id))
                    self.count.inc(str(newid), w.wgt)

        return True

setattr(XZZLHEWeightAnalyzer,"defaultConfig",
    cfg.Analyzer(XZZLHEWeightAnalyzer,
    )
)
