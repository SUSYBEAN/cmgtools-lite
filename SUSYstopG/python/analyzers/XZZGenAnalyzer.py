from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsutils.genutils import *
from PhysicsTools.HeppyCore.statistics.counter import Counter, Counters
from CMGTools.XZZ2l2nu.tools.Pair import Pair
import math

class XZZGenAnalyzer( Analyzer ):
    """ Only select X->ZZ->2l2nu events

       """

    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(XZZGenAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)

        self.filter = getattr(self.cfg_ana, 'filter', "None")
 
    def declareHandles(self):
        super(XZZGenAnalyzer, self).declareHandles()
        self.mchandles['prunedGenParticles'] = AutoHandle( 'prunedGenParticles', 'std::vector<reco::GenParticle>' )
        self.mchandles['packedGenParticles'] = AutoHandle( 'packedGenParticles', 'std::vector<pat::PackedGenParticle>' )
        self.mchandles['LHEinfo'] = AutoHandle('externalLHEProducer',
                                                  'LHEEventProduct',
                                                  mayFail=True,
                                                  fallbackLabel='source',
                                                  lazy=False )
               
    def beginLoop(self,setup):
        super(XZZGenAnalyzer,self).beginLoop(setup)

        self.counters.addCounter('XZZGenReport')
        self.count = self.counters.counter('XZZGenReport')
        self.count.register('XZZ2l2nu Events')
        self.count.register('XZZ2l2j Events')
        self.count.register('emu2nu Events')
        self.count.register('pass events')
        

    def makeMCInfo(self, event, bosonID = 23):
        verbose = getattr(self.cfg_ana, 'verbose', False)
        pruned = self.mchandles['prunedGenParticles'].product()
        packed = self.mchandles['packedGenParticles'].product()

        event.genParticles = [ p for p in packed ]
        selectedGenParticles = []
        event.genZBosons = []
        event.genLeptons = []
        event.mother_genLeptons = []
        event.genElectrons = []
        event.genMuons = []
        event.genTaus = []
	event.genLeptonsFsr = []
	event.genElectronsFsr = []
	event.genMuonsFsr = []
	event.genTausFsr = []
        event.genNeutrinos = []
        event.genJets = []
        event.genXZZ = []

        #event.genZBosons = [ p for p in pruned if (abs(p.pdgId()) == bosonID) and p.numberOfDaughters() > 0 ]
        
        #for zboson in event.genZBosons:
        #    selectedGenParticles.append(zboson)
        #    for i in xrange( zboson.numberOfDaughters() ):
        #        dau = zboson.daughter(i)
        #        dauid = dau.pdgId()
        #        if abs(dauid) in [11,13]:
        #            event.genLeptons.append(dau)
        #            selectedGenParticles.append(dau)
        #            if abs(dauid)==11: event.genElectrons.append(dau)
        #            if abs(dauid)==13: event.genMuons.append(dau)
        #        elif abs(dauid) in [12,14,16]:
        #            event.genNeutrinos.append(dau)
        #            selectedGenParticles.append(dau)
        #        elif abs(dauid) in range(7):
        #            event.genJets.append(dau)

        event.genMuons = [ p for p in pruned if abs(p.pdgId())==13 and  p.status()==1 and p.isPromptFinalState() and p.fromHardProcessFinalState() ]
        event.genElectrons = [ p for p in pruned if abs(p.pdgId())==11 and  p.status()==1 and p.isPromptFinalState() and p.fromHardProcessFinalState() ]
        #event.genTaus = [ p for p in pruned if abs(p.pdgId())==15 and p.status()==2  and p.isPromptDecayed() ]
        event.genTaus = [ p for p in pruned if abs(p.pdgId())==15 and p.status()==2  and p.isPromptDecayed() and p.fromHardProcessDecayed() ]
        event.genNeutrinos = [ p for p in pruned if (abs(p.pdgId()) in [12,14,16]) and p.isPromptFinalState() and p.fromHardProcessFinalState() ]
        event.genMuonsFsr = [ p for p in pruned if abs(p.pdgId())==13 and  p.status()==1 and p.isPromptFinalState() and not p.fromHardProcessFinalState() ]
        event.genElectronsFsr = [ p for p in pruned if abs(p.pdgId())==11 and  p.status()==1 and p.isPromptFinalState() and not p.fromHardProcessFinalState() ]
        event.genTausFsr = [ p for p in pruned if abs(p.pdgId())==15 and p.status()==2  and p.isPromptDecayed() and not p.fromHardProcessDecayed() ]

        event.genLeptons = event.genMuons + event.genElectrons + event.genTaus
        #event.genLeptonsFsr = event.genMuonsFsr + event.genElectronsFsr
        event.genLeptonsFsr = event.genMuonsFsr + event.genElectronsFsr + event.genTausFsr

        if len(event.genLeptonsFsr)==0:
	    # for no fsr case
            if len(event.genLeptons)==2:
                Z = event.genLeptons[0].p4()+event.genLeptons[1].p4()
                event.genZBosons.append(Z)
        elif len(event.genLeptonsFsr)>=2 and len(event.genLeptons)==2:
            # for fsr case, find the mother Fsr leptons and create the Z
            fsr_mothers = [ p.mother(0) for p in event.genLeptonsFsr ]
            event.mother_genLeptons = event.genLeptons
            for ip, p in enumerate(event.mother_genLeptons):
                if p.mother(0) in fsr_mothers:
                    event.mother_genLeptons[ip] = p.mother(0)
            Z = event.mother_genLeptons[0].p4()+event.mother_genLeptons[1].p4()
            event.genZBosons.append(Z)

        # neutrinos to create another Z
        if len(event.genNeutrinos)==2:
            Z = event.genNeutrinos[0].p4()+event.genNeutrinos[1].p4()
            event.genZBosons.append(Z)
                    
            
            

        if len(event.genZBosons)>=2 and len(event.genLeptons)>=2 and len(event.genNeutrinos)>=2:
            event.genIsXZZ2l2nu = True
        else: 
            event.genIsXZZ2l2nu = False
 
        if len(event.genZBosons)>=2 and len(event.genLeptons)>=2 and len(event.genJets)>=2:
            event.genIsXZZ2l2j = True
        else:
            event.genIsXZZ2l2j = False
            
        if len(event.genZBosons)>=2 and len(event.genElectrons)>=1 and len(event.genMuons)>=1 and len(event.genNeutrinos)>=2:
            event.genIsXZZemu2nu = True
        else:
            event.genIsXZZemu2nu = False

        if event.genIsXZZemu2nu:
            self.count.inc('emu2nu Events')

        if event.genIsXZZ2l2nu:
            self.count.inc('XZZ2l2nu Events') 

        if event.genIsXZZ2l2j:
            self.count.inc('XZZ2l2j Events')

        
        if len(event.genZBosons)>=2 :
            genX = event.genZBosons[0]+event.genZBosons[1]
            event.genXZZ.append(genX)
            

        if self.cfg_ana.verbose and event.genIsXZZ2l2nu:
            print "N gen Z bosons: "+str(len(event.genZBosons))
            print "N gen Leptons: "+str(len(event.genLeptons))
            print "N gen neutrinos: "+str(len(event.genNeutrinos))

        # get n partons info
        event.lheNb = 0
        event.lheNj = 0
        lheEvent = self.mchandles['LHEinfo'].product().hepeup();
        lheParticles = lheEvent.PUP;
        for idxParticle in range(len(lheParticles)):
            idx = abs(lheEvent.IDUP[idxParticle])
            status = lheEvent.ISTUP[idxParticle]
            if status == 1 and idx==5:  event.lheNb += 1
            if status == 1 and ((idx >= 1 and idx <= 6) or idx == 21) : event.lheNj += 1
              

    def process(self, event):
        self.readCollections( event.input )

        

        # if not MC, nothing to do
        if not self.cfg_comp.isMC: 
            return True
        # do MC level analysis
        self.makeMCInfo(event)
        
        # no gen filter
        if self.filter == "None":
            self.count.inc('pass events')
            return True
        # filter
        if self.filter=="ZZ2mu":
            if len(event.genMuons)>=2: 
                self.count.inc('pass events')
                return True
        elif self.filter=="ZZ2el":
            if len(event.genElectrons)>=2:
                self.count.inc('pass events')
                return True
        elif self.filter=="ZZ2l2nu":
            if event.genIsXZZ2l2nu:
                self.count.inc('pass events')
                return True

        # false if use filter but not pass any
        return False


import PhysicsTools.HeppyCore.framework.config as cfg
setattr(XZZGenAnalyzer,"defaultConfig",
    cfg.Analyzer(XZZGenAnalyzer,
        # Print out debug information
        verbose = False,
    )
)
