import math, os, random, ROOT
from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsobjects.PhysicsObjects import Jet #include definition for jetID etc.
from PhysicsTools.HeppyCore.utils.deltar import deltaR2, deltaPhi, matchObjectCollection, matchObjectCollection2, bestMatch,matchObjectCollection3
#from PhysicsTools.Heppy.physicsutils.JetReCalibrator import JetReCalibrator
import PhysicsTools.HeppyCore.framework.config as cfg

from PhysicsTools.Heppy.physicsutils.QGLikelihoodCalculator import QGLikelihoodCalculator
from PhysicsTools.HeppyCore.statistics.counter import Counter, Counters

from CMGTools.XZZ2l2nu.analyzers.JetResolution import JetResolution
from CMGTools.XZZ2l2nu.analyzers.JetReCalibrator import JetReCalibrator

import copy

def XZZmatchObjectCollection( objects, matchCollection,
                              deltaR2Max, dPtMaxFactor, 
                              filter = lambda x,y : True):
    ''' @ Apr-6-16 Adds genJets matching for JER according to 
    https://github.com/cms-sw/cmssw/blob/518995a62433825453d7ccf9ec222d782086cad9/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
    
    used to produce the sync table:
    https://twiki.cern.ch/twiki/bin/view/CMS/JERCReference
    '''
    pairs = {}
    if len(objects)==0:
        return pairs
    if len(matchCollection)==0:
        return dict( zip(objects, [None]*len(objects)) )
    for object in objects:
        bm, dr2 = bestMatch( object, [mob for mob in matchCollection if filter(object,mob)])
        dPt = abs(bm.pt() - object.pt())
        if dr2<deltaR2Max and dPt < dPtMaxFactor * object.pt() * object.resolution:
            pairs[object] = bm
        else:
            pairs[object] = None
    return pairs

def cleanNearestJetOnly(jets,leptons,deltaR):
    dr2 = deltaR**2
    good = [ True for j in jets ]
    for l in leptons:
        ibest, d2m = -1, dr2
        for i,j in enumerate(jets):
            d2i = deltaR2(l.eta(),l.phi(), j.eta(),j.phi())
            if d2i < d2m:
                ibest, d2m = i, d2i
        if ibest != -1: good[ibest] = False
    return [ j for (i,j) in enumerate(jets) if good[i] == True ] 

def cleanJetsAndLeptons(jets,leptons,deltaR,arbitration):
    dr2 = deltaR**2
    goodjet = [ True for j in jets ]
    goodlep = [ True for l in leptons ]
    for il, l in enumerate(leptons):
        ibest, d2m = -1, dr2
        for i,j in enumerate(jets):
            d2i = deltaR2(l.eta(),l.phi(), j.eta(),j.phi())
            if d2i < dr2:
                choice = arbitration(j,l)
                if choice == j:
                   # if the two match, and we prefer the jet, then drop the lepton and be done
                   goodlep[il] = False
                   break 
                elif choice == (j,l) or choice == (l,j):
                   # asked to keep both, so we don't consider this match
                   continue
            if d2i < d2m:
                ibest, d2m = i, d2i
        # this lepton has been killed by a jet, then we clean the jet that best matches it
        if not goodlep[il]: continue 
        if ibest != -1: goodjet[ibest] = False
    return ( [ j for (i ,j) in enumerate(jets)    if goodjet[i ] == True ], 
             [ l for (il,l) in enumerate(leptons) if goodlep[il] == True ] )



# def shiftJERfactor(JERShift, aeta):
#         factor = 1.079 + JERShift*0.026
#         if   aeta > 3.2: factor = 1.056 + JERShift * 0.191
#         elif aeta > 2.8: factor = 1.395 + JERShift * 0.063
#         elif aeta > 2.3: factor = 1.254 + JERShift * 0.062
#         elif aeta > 1.7: factor = 1.208 + JERShift * 0.046
#         elif aeta > 1.1: factor = 1.121 + JERShift * 0.029
#         elif aeta > 0.5: factor = 1.099 + JERShift * 0.028
#         return factor 


class XZZJetAnalyzer( Analyzer ):
    """Taken from RootTools.JetAnalyzer, simplified, modified, added corrections    
       Please sync your jet output and MET output with: (Apr-2016)
       https://twiki.cern.ch/twiki/bin/view/CMS/JERCReference#Reference_table_for_MC
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.debug = getattr(cfg_ana, 'debug', False)

        super(XZZJetAnalyzer,self).__init__(cfg_ana, cfg_comp, looperName)
        mcGT   = cfg_ana.mcGT   if hasattr(cfg_ana,'mcGT')   else "PHYS14_25_V2"
        dataGT = cfg_ana.dataGT if hasattr(cfg_ana,'dataGT') else "GR_70_V2_AN1"

        mcGT_jer   = cfg_ana.mcGT_jer   if hasattr(cfg_ana,'mcGT_jer')   else "Summer15_25nsV6_MC"
        dataGT_jer = cfg_ana.dataGT_jer if hasattr(cfg_ana,'dataGT_jer') else "Summer15_25nsV6_DATA"
        
        self.shiftJEC = self.cfg_ana.shiftJEC if hasattr(self.cfg_ana, 'shiftJEC') else 0
        self.recalibrateJets = self.cfg_ana.recalibrateJets
        self.addJECShifts = self.cfg_ana.addJECShifts if hasattr(self.cfg_ana, 'addJECShifts') else 0
        if   self.recalibrateJets == "MC"  : self.recalibrateJets =     self.cfg_comp.isMC
        elif self.recalibrateJets == "Data": self.recalibrateJets = not self.cfg_comp.isMC
        elif self.recalibrateJets not in [True,False]: raise RuntimeError, "recalibrateJets must be any of { True, False, 'MC', 'Data' }, while it is %r " % self.recalibrateJets
       
        calculateSeparateCorrections = getattr(cfg_ana,"calculateSeparateCorrections", False);
        calculateType1METCorrection  = getattr(cfg_ana,"calculateType1METCorrection",  False);
        self.doJEC = self.recalibrateJets or (self.shiftJEC != 0) or self.addJECShifts or calculateSeparateCorrections or calculateType1METCorrection
        if self.doJEC:
            if self.debug: print "[Debug] I am doing Jet energy calibration :D"
            doResidual = getattr(cfg_ana, 'applyL2L3Residual', 'Data')
            if   doResidual == "MC":   doResidual = self.cfg_comp.isMC
            elif doResidual == "Data": doResidual = not self.cfg_comp.isMC
            elif doResidual not in [True,False]: raise RuntimeError, "If specified, applyL2L3Residual must be any of { True, False, 'MC', 'Data'(default)}"
            GT = getattr(cfg_comp, 'jecGT', mcGT if self.cfg_comp.isMC else dataGT)
            GT_jer= getattr(cfg_comp, 'jerGT', mcGT_jer if self.cfg_comp.isMC else dataGT_jer)
            # Now take care of the optional arguments
            kwargs = { 'calculateSeparateCorrections':calculateSeparateCorrections,
                       'calculateType1METCorrection' :calculateType1METCorrection, }
            if kwargs['calculateType1METCorrection']: kwargs['type1METParams'] = cfg_ana.type1METParams
            # instantiate the jet re-calibrator
            if self.debug: print "[Debug] check input for jetReCalibrator: tell me if it is MC -- %r; tell me the doResidual is %r\n" % (self.cfg_comp.isMC, doResidual)
            
            self.jetReCalibrator = JetReCalibrator(GT, cfg_ana.recalibrationType, doResidual, cfg_ana.jecPath, **kwargs)

        self.smearJets = getattr(self.cfg_ana, 'smearJets', False)
        if self.smearJets:  
            self.jetResolution = JetResolution(GT_jer, cfg_ana.recalibrationType, cfg_ana.jerPath)
        self.doPuId = getattr(self.cfg_ana, 'doPuId', True)
        self.jetLepDR = getattr(self.cfg_ana, 'jetLepDR', 0.4)
        self.jetLepArbitration = getattr(self.cfg_ana, 'jetLepArbitration', lambda jet,lepton: lepton) 
        self.lepPtMin = getattr(self.cfg_ana, 'minLepPt', -1)
        self.lepSelCut = getattr(self.cfg_ana, 'lepSelCut', lambda lep : True)
        self.jetGammaDR =  getattr(self.cfg_ana, 'jetGammaDR', 0.4)
        if(self.cfg_ana.doQG):
            qgdefname="{CMSSW_BASE}/src/PhysicsTools/Heppy/data/pdfQG_AK4chs_13TeV_v2b.root"
            self.qglcalc = QGLikelihoodCalculator(getattr(self.cfg_ana,"QGpath",qgdefname).format(CMSSW_BASE= os.environ['CMSSW_BASE']))
        if not hasattr(self.cfg_ana ,"collectionPostFix"):self.cfg_ana.collectionPostFix=""

    def declareHandles(self):
        super(XZZJetAnalyzer, self).declareHandles()
        self.handles['jets']   = AutoHandle( self.cfg_ana.jetCol, 'std::vector<pat::Jet>' )
        #self.handles['jets_raw']   = AutoHandle( self.cfg_ana.jetCol, 'std::vector<pat::Jet>' )
        self.handles['genJet'] = AutoHandle( self.cfg_ana.genJetCol, 'std::vector<reco::GenJet>' )
        self.shiftJER = self.cfg_ana.shiftJER if hasattr(self.cfg_ana, 'shiftJER') else 0
        self.addJERShifts = self.cfg_ana.addJERShifts if hasattr(self.cfg_ana, 'addJERShifts') else 0
        self.handles['rho'] = AutoHandle( self.cfg_ana.rho[0], 'double' )
        self.handles['rho_jer'] = AutoHandle( self.cfg_ana.rho[1] if self.cfg_ana.rho[1]!='' else self.cfg_ana.rho[0], 'double' )
        
    def beginLoop(self, setup):
        super(XZZJetAnalyzer,self).beginLoop(setup)
        self.counters.addCounter('jets')
        count = self.counters.counter('jets')
        count.register('All Events')
    #    count.register('veto b jest')

    def process(self, event):
        self.readCollections( event.input )
        rho  = float(self.handles['rho'].product()[0])
        rho_jer = float(self.handles['rho_jer'].product()[0])
        self.rho = rho

        ## Read jets, if necessary recalibrate and shift MET
        #self.jets_raw = map(Jet, self.handles['jets_raw'].product())
        if self.cfg_ana.copyJetsByValue: 
          import ROOT
          allJets = map(lambda j:Jet(ROOT.pat.Jet(ROOT.edm.Ptr(ROOT.pat.Jet)(ROOT.edm.ProductID(),j,0))), self.handles['jets'].product()) #copy-by-value is safe if JetAnalyzer is ran more than once
        else: 
          allJets = map(Jet, self.handles['jets'].product()) 
        self.jets_raw = copy.deepcopy(allJets)
        #set dummy MC flavour for all jets in case we want to ntuplize discarded jets later
        for jet in allJets:
            jet.mcFlavour = 0

        self.deltaMetFromJEC = [0.,0.]
        self.type1METCorr    = [0.,0.,0.]
        self.type1METCorrUp  = [0.,0.,0.]
        self.type1METCorrDown = [0.,0.,0.]
        self.sumRawJetsforT1 = [0.,0.,0.]

#        print "before. rho",self.rho,self.cfg_ana.collectionPostFix,'allJets len ',len(allJets),'pt', [j.pt() for j in allJets]

        if self.doJEC:
            if not self.recalibrateJets:  # check point that things won't change
                jetsBefore = [ (j.pt(),j.eta(),j.phi(),j.rawFactor()) for j in allJets ]
            #print '[Debug] I am event = ', event.input.eventAuxiliary().id().event()
                
            self.jetReCalibrator.correctAll(allJets, rho, delta=self.shiftJEC, 
                                            addCorr=True, addShifts=self.addJECShifts,
                                            metShift=self.deltaMetFromJEC, type1METCorr=self.type1METCorr,
                                            type1METCorrUp=self.type1METCorrUp, type1METCorrDown=self.type1METCorrDown,
                                            sumRawJetsforT1=self.sumRawJetsforT1)
            
            if not self.recalibrateJets: 
                jetsAfter = [ (j.pt(),j.eta(),j.phi(),j.rawFactor()) for j in allJets ]
                if len(jetsBefore) != len(jetsAfter): 
                    print "ERROR: I had to recompute jet corrections, and they rejected some of the jets:\nold = %s\n new = %s\n" % (jetsBefore,jetsAfter)
                else:
                    for (told, tnew) in zip(jetsBefore,jetsAfter):
                        if (deltaR2(told[1],told[2],tnew[1],tnew[2])) > 0.0001:
                            print "ERROR: I had to recompute jet corrections, and one jet direction moved: old = %s, new = %s\n" % (told, tnew)
                        elif abs(told[0]-tnew[0])/(told[0]+tnew[0]) > 0.5e-3 or abs(told[3]-tnew[3]) > 0.5e-3:
                            print "ERROR: I had to recompute jet corrections, and one jet pt or corr changed: old = %s, new = %s\n" % (told, tnew)
        self.allJetsUsedForMET = allJets
        self.sumJetsInT1 = {"rawP4forT1": self.sumRawJetsforT1,
                            "type1METCorr":self.type1METCorr, 
                            "corrP4forT1":[x + y for x,y in zip(self.sumRawJetsforT1,self.type1METCorr)]}

#        print "after. rho",self.rho,self.cfg_ana.collectionPostFix,'allJets len ',len(allJets),'pt', [j.pt() for j in allJets]

        if self.cfg_comp.isMC:
            self.genJets = [ x for x in self.handles['genJet'].product() ]
        
        # ==> following matching example from PhysicsTools/Heppy/python/analyzers/examples/JetAnalyzer.py
        if self.cfg_comp.isMC:
            for jet in allJets:
                if self.smearJets:
                    jet.resolution = self.jetResolution.getResolution(jet,rho_jer)
                    jet.jerfactor = self.jetResolution.getScaleFactor(jet)
                else:
                    jet.resolution = 1.0
                    jet.jerfactor = 1.0

                if self.genJets:
                    # Use DeltaR = 0.2 matching jet and genJet from 
                    # https://github.com/blinkseb/JMEReferenceSample/blob/master/test/createReferenceSample.py
                    pairs = XZZmatchObjectCollection( [jet], self.genJets, 0.2*0.2, 3)
                    if pairs[jet] is None:
                        pass
                    else:
                        jet.matchedGenJet = pairs[jet]
                        jet.matchedGenJetIdx = self.genJets.index(pairs[jet]) 
                else: pass
                if self.debug: print '[Debug] I am event = ', event.input.eventAuxiliary().id().event()
                if self.smearJets:
                    self.jerCorrection(jet, rho_jer)

	##Sort Jets by pT 
        allJets.sort(key = lambda j : j.pt(), reverse = True)
        
        leptons = []
        if hasattr(event, 'selectedLeptons'):
            leptons = [ l for l in event.selectedLeptons if l.pt() > self.lepPtMin and self.lepSelCut(l) ]
        if self.cfg_ana.cleanJetsFromTaus and hasattr(event, 'selectedTaus'):
            leptons = leptons[:] + event.selectedTaus
        if self.cfg_ana.cleanJetsFromIsoTracks and hasattr(event, 'selectedIsoCleanTrack'):
            leptons = leptons[:] + event.selectedIsoCleanTrack

        self.counters.counter('jets').inc('All Events')
	## Apply jet selection
        self.jets = []
        self.jetsFailId = []
        self.jetsAllNoID = []
        self.jetsIdOnly = []
        for jet in allJets:
            #Check if lepton and jet have overlapping PF candidates 
            leps_with_overlaps = []
            if getattr(self.cfg_ana, 'checkLeptonPFOverlap', True):
              for i in range(jet.numberOfSourceCandidatePtrs()):
                p1 = jet.sourceCandidatePtr(i) #Ptr<Candidate> p1
                for lep in leptons:
                    for j in range(lep.numberOfSourceCandidatePtrs()):
                        p2 = lep.sourceCandidatePtr(j)
                        has_overlaps = p1.key() == p2.key() and p1.refCore().id().productIndex() == p2.refCore().id().productIndex() and p1.refCore().id().processIndex() == p2.refCore().id().processIndex()
                        if has_overlaps:
                            leps_with_overlaps += [lep]
            if len(leps_with_overlaps)>0:
                for lep in leps_with_overlaps:
                    lep.jetOverlap = jet
            if self.testJetNoID( jet ): 
                self.jetsAllNoID.append(jet) 
                if(self.cfg_ana.doQG):
                    jet.qgl_calc =  self.qglcalc.computeQGLikelihood
                    jet.qgl_rho =  rho
                if self.testJetID( jet ):
                    self.jets.append(jet)
                    self.jetsIdOnly.append(jet)
                else:
                    self.jetsFailId.append(jet)
            elif self.testJetID (jet ):
                self.jetsIdOnly.append(jet)

        jetsEtaCut = [j for j in self.jets if abs(j.eta()) <  self.cfg_ana.jetEta ]
        self.cleanJetsAll, cleanLeptons = cleanJetsAndLeptons(jetsEtaCut, leptons, self.jetLepDR, self.jetLepArbitration)

        self.cleanJets    = [j for j in self.cleanJetsAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        self.cleanJetsFwd = [j for j in self.cleanJetsAll if abs(j.eta()) >= self.cfg_ana.jetEtaCentral ]
        self.discardedJets = [j for j in self.jets if j not in self.cleanJetsAll]
        if hasattr(event, 'selectedLeptons') and self.cfg_ana.cleanSelectedLeptons:
            event.discardedLeptons = [ l for l in leptons if l not in cleanLeptons ]
            event.selectedLeptons  = [ l for l in event.selectedLeptons if l not in event.discardedLeptons ]
        for lep in leptons:
            if hasattr(lep, "jetOverlap"):
                if lep.jetOverlap in self.cleanJetsAll:
                    #print "overlap reco", lep.p4().pt(), lep.p4().eta(), lep.p4().phi(), lep.jetOverlap.p4().pt(), lep.jetOverlap.p4().eta(), lep.jetOverlap.p4().phi()
                    lep.jetOverlapIdx = self.cleanJetsAll.index(lep.jetOverlap)
                elif lep.jetOverlap in self.discardedJets:
                    #print "overlap discarded", lep.p4().pt(), lep.p4().eta(), lep.p4().phi(), lep.jetOverlap.p4().pt(), lep.jetOverlap.p4().eta(), lep.jetOverlap.p4().phi()
                    lep.jetOverlapIdx = 1000 + self.discardedJets.index(lep.jetOverlap)

        ## First cleaning, then Jet Id
        self.noIdCleanJetsAll, cleanLeptons = cleanJetsAndLeptons(self.jetsAllNoID, leptons, self.jetLepDR, self.jetLepArbitration)
        self.noIdCleanJets = [j for j in self.noIdCleanJetsAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        self.noIdCleanJetsFwd = [j for j in self.noIdCleanJetsAll if abs(j.eta()) >=  self.cfg_ana.jetEtaCentral ]
        self.noIdDiscardedJets = [j for j in self.jetsAllNoID if j not in self.noIdCleanJetsAll]

        ## Clean Jets from photons (first cleaning, then Jet Id)
        photons = []
        if hasattr(event, 'selectedPhotons'):
            if self.cfg_ana.cleanJetsFromFirstPhoton:
                photons = event.selectedPhotons[:1]
            else:
                photons = [ g for g in event.selectedPhotons ] 

        self.gamma_cleanJetsAll = cleanNearestJetOnly(self.cleanJetsAll, photons, self.jetGammaDR)
        self.gamma_cleanJets    = [j for j in self.gamma_cleanJetsAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        self.gamma_cleanJetsFwd = [j for j in self.gamma_cleanJetsAll if abs(j.eta()) >= self.cfg_ana.jetEtaCentral ]
        self.gamma_noIdCleanJetsAll = cleanNearestJetOnly(self.noIdCleanJetsAll, photons, self.jetGammaDR)
        self.gamma_noIdCleanJets    = [j for j in self.gamma_noIdCleanJetsAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        self.gamma_noIdCleanJetsFwd = [j for j in self.gamma_noIdCleanJetsAll if abs(j.eta()) >= self.cfg_ana.jetEtaCentral ]
        ###

        if self.cfg_ana.alwaysCleanPhotons:
            self.cleanJets = self.gamma_cleanJets
            self.cleanJetsAll = self.gamma_cleanJetsAll
            self.cleanJetsFwd = self.gamma_cleanJetsFwd
            #
            self.noIdCleanJets = self.gamma_noIdCleanJets
            self.noIdCleanJetsAll = self.gamma_noIdCleanJetsAll
            self.noIdCleanJetsFwd = self.gamma_noIdCleanJetsFwd

        ## Jet Id, after jet/lepton cleaning
        self.cleanJetsFailIdAll = []
        for jet in self.noIdCleanJetsAll:
            if not self.testJetID( jet ):
                self.cleanJetsFailIdAll.append(jet)
        
        self.cleanJetsFailId = [j for j in self.cleanJetsFailIdAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        
        ## Jet Id, after jet/photon cleaning
        self.gamma_cleanJetsFailIdAll = []
        for jet in self.gamma_noIdCleanJetsAll:
            if not self.testJetID( jet ):
                self.gamma_cleanJetsFailIdAll.append(jet)

        self.gamma_cleanJetsFailId = [j for j in self.gamma_cleanJetsFailIdAll if abs(j.eta()) <  self.cfg_ana.jetEtaCentral ]
        
        ## Associate jets to leptons
        incleptons = event.selectedLeptons if hasattr(event, 'selectedLeptons') else []
        if hasattr(event, 'inclusiveLeptons'): 
            incleptons = event.inclusiveLeptons
        jlpairs = matchObjectCollection(incleptons, allJets, self.jetLepDR**2)

        for jet in allJets:
            jet.leptons = [l for l in jlpairs if jlpairs[l] == jet ]
        for lep in incleptons:
            jet = jlpairs[lep]
            if jet is None:
                setattr(lep,"jet"+self.cfg_ana.collectionPostFix,lep)
            else:
                setattr(lep,"jet"+self.cfg_ana.collectionPostFix,jet)
        ## Associate jets to taus 
        taus = getattr(event,'selectedTaus',[])
        jtaupairs = matchObjectCollection( taus, allJets, self.jetLepDR**2)

        for jet in allJets:
            jet.taus = [l for l in jtaupairs if jtaupairs[l] == jet ]
        for tau in taus:
            setattr(tau,"jet"+self.cfg_ana.collectionPostFix,jtaupairs[tau])

        #MC stuff
        if self.cfg_comp.isMC:
            self.deltaMetFromJetSmearing = [0, 0]
            for j in self.cleanJetsAll:
                if hasattr(j, 'deltaMetFromJetSmearing'):
                    self.deltaMetFromJetSmearing[0] += j.deltaMetFromJetSmearing[0]
                    self.deltaMetFromJetSmearing[1] += j.deltaMetFromJetSmearing[1]

            self.cleanGenJets = cleanNearestJetOnly(self.genJets, leptons, self.jetLepDR)
            
            if self.cfg_ana.cleanGenJetsFromPhoton:
                self.cleanGenJets = cleanNearestJetOnly(self.cleanGenJets, photons, self.jetLepDR)

	    if getattr(self.cfg_ana, 'attachNeutrinos', True) and hasattr(self.cfg_ana,"genNuSelection") :
		jetNus=[x for x in event.genParticles if abs(x.pdgId()) in [12,14,16] and self.cfg_ana.genNuSelection(x) ]
	 	pairs= matchObjectCollection (jetNus, self.genJets, 0.4**2)
                
		for (nu,genJet) in pairs.iteritems() :
		     if genJet is not None :
			if not hasattr(genJet,"nu") :
				genJet.nu=nu.p4()
			else :
				genJet.nu+=nu.p4()
			
            
            if self.cfg_ana.do_mc_match:
                self.jetFlavour(event)

        if hasattr(event,"jets"+self.cfg_ana.collectionPostFix): raise RuntimeError, "Event already contains a jet collection with the following postfix: "+self.cfg_ana.collectionPostFix
        setattr(event,"rho"                    +self.cfg_ana.collectionPostFix, self.rho                    ) 
        setattr(event,"deltaMetFromJEC"        +self.cfg_ana.collectionPostFix, self.deltaMetFromJEC        ) 
        setattr(event,"type1METCorr"           +self.cfg_ana.collectionPostFix, self.type1METCorr           ) 
        setattr(event,"type1METCorrUp"         +self.cfg_ana.collectionPostFix, self.type1METCorrUp         ) 
        setattr(event,"type1METCorrDown"       +self.cfg_ana.collectionPostFix, self.type1METCorrDown       )
        setattr(event,"sumJetsInT1"            +self.cfg_ana.collectionPostFix, self.sumJetsInT1            )
        setattr(event,"allJetsUsedForMET"      +self.cfg_ana.collectionPostFix, self.allJetsUsedForMET      ) 
        setattr(event,"jets_raw"               +self.cfg_ana.collectionPostFix, self.jets_raw               )
        setattr(event,"jets"                   +self.cfg_ana.collectionPostFix, self.jets                   ) 
        setattr(event,"jetsFailId"             +self.cfg_ana.collectionPostFix, self.jetsFailId             ) 
        setattr(event,"jetsAllNoID"            +self.cfg_ana.collectionPostFix, self.jetsAllNoID            ) 
        setattr(event,"jetsIdOnly"             +self.cfg_ana.collectionPostFix, self.jetsIdOnly             ) 
        setattr(event,"cleanJetsAll"           +self.cfg_ana.collectionPostFix, self.cleanJetsAll           ) 
        setattr(event,"cleanJets"              +self.cfg_ana.collectionPostFix, self.cleanJets              ) 
        setattr(event,"cleanJetsFwd"           +self.cfg_ana.collectionPostFix, self.cleanJetsFwd           ) 
        setattr(event,"cleanJetsFailIdAll"           +self.cfg_ana.collectionPostFix, self.cleanJetsFailIdAll           ) 
        setattr(event,"cleanJetsFailId"              +self.cfg_ana.collectionPostFix, self.cleanJetsFailId              ) 
        setattr(event,"discardedJets"          +self.cfg_ana.collectionPostFix, self.discardedJets          ) 
        setattr(event,"gamma_cleanJetsAll"     +self.cfg_ana.collectionPostFix, self.gamma_cleanJetsAll     ) 
        setattr(event,"gamma_cleanJets"        +self.cfg_ana.collectionPostFix, self.gamma_cleanJets        ) 
        setattr(event,"gamma_cleanJetsFwd"     +self.cfg_ana.collectionPostFix, self.gamma_cleanJetsFwd     ) 
        setattr(event,"gamma_cleanJetsFailIdAll"     +self.cfg_ana.collectionPostFix, self.gamma_cleanJetsFailIdAll     ) 
        setattr(event,"gamma_cleanJetsFailId"        +self.cfg_ana.collectionPostFix, self.gamma_cleanJetsFailId        ) 


        if self.cfg_comp.isMC:
            setattr(event,"deltaMetFromJetSmearing"+self.cfg_ana.collectionPostFix, self.deltaMetFromJetSmearing) 
            setattr(event,"cleanGenJets"           +self.cfg_ana.collectionPostFix, self.cleanGenJets           )
            setattr(event,"genJets"                +self.cfg_ana.collectionPostFix, self.genJets                )
            if self.cfg_ana.do_mc_match:
                setattr(event,"bqObjects"              +self.cfg_ana.collectionPostFix, self.bqObjects              )
                setattr(event,"cqObjects"              +self.cfg_ana.collectionPostFix, self.cqObjects              )
                setattr(event,"partons"                +self.cfg_ana.collectionPostFix, self.partons                )
                setattr(event,"heaviestQCDFlavour"     +self.cfg_ana.collectionPostFix, self.heaviestQCDFlavour     )

 
        return True

        

    def testJetID(self, jet):
        jet.puJetIdPassed = jet.puJetId() 
        jet.pfJetIdPassed = jet.jetID('POG_PFID_Loose') if jet.isPFJet() else False
        if self.cfg_ana.relaxJetId:
            return True
        else:
            return jet.pfJetIdPassed and (jet.puJetIdPassed or not(self.doPuId)) 
        
    def testJetNoID( self, jet ):
        # 2 is loose pile-up jet id
        return jet.pt() > self.cfg_ana.jetPt and \
               abs( jet.eta() ) < self.cfg_ana.jetEta;

    def jetFlavour(self,event):
        def isFlavour(x,f):
            id = abs(x.pdgId())
            if id > 999: return (id/1000)%10 == f
            if id >  99: return  (id/100)%10 == f
            return id % 100 == f



        self.bqObjects = [ p for p in event.genParticles if (p.status() == 2 and isFlavour(p,5)) ]
        self.cqObjects = [ p for p in event.genParticles if (p.status() == 2 and isFlavour(p,4)) ]

        self.partons   = [ p for p in event.genParticles if ((p.status() == 23 or p.status() == 3) and abs(p.pdgId())>0 and (abs(p.pdgId()) in [1,2,3,4,5,21]) ) ]
        match = matchObjectCollection2(self.cleanJetsAll,
                                       self.partons,
                                       deltaRMax = 0.3)

        for jet in self.cleanJetsAll:
            parton = match[jet]
            jet.partonId = (parton.pdgId() if parton != None else 0)
            jet.partonMotherId = (parton.mother(0).pdgId() if parton != None and parton.numberOfMothers()>0 else 0)
        
        for jet in self.jets:
           (bmatch, dr) = bestMatch(jet, self.bqObjects)
           if dr < 0.4:
               jet.mcFlavour = 5
           else:
               (cmatch, dr) = bestMatch(jet, self.cqObjects) 
               if dr < 0.4:
                   jet.mcFlavour = 4
               else:
                   jet.mcFlavour = 0

        self.heaviestQCDFlavour = 5 if len(self.bqObjects) else (4 if len(self.cqObjects) else 1);
    
    def jerCorrection(self, jet, rho_jer):
        ''' @ Apr-6-16 Adds JER correction 
        according to the recommended 'hybrid' method at
        https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        
        algorithm following:
        https://github.com/cms-sw/cmssw/blob/518995a62433825453d7ccf9ec222d782086cad9/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
        '''
        res = jet.resolution
        factor = jet.jerfactor
        ptScale = 1.0
        jetpt = jet.pt()
        #jetpt=jet.energy()
        if hasattr(jet, 'matchedGenJet'):
            genpt = jet.matchedGenJet.pt()
            #genpt=jet.matchedGenJet.energy()
            ptScale = 1 + (factor - 1)*(jetpt - genpt)/jetpt
            if self.debug: print '[Debug] I am scaling jet: jet(pT, eta, rho, rho_jer) = (%.4f, %.4f, %.4f, %.4f)' % (jet.pt(), jet.eta(), self.rho, rho_jer) 
        elif factor > 1.0:
            #random.seed(37428479) # for syncing to the reference table
            ran = ROOT.TRandom(37428479)
            #ptScale = 1 + random.gauss(0, res*math.sqrt(factor**2 - 1))/jetpt
            ptScale = 1 + ran.Gaus(0, res*math.sqrt(factor**2 - 1))/jetpt
            if self.debug: print '[Debug] I am smearing jet: jet(pT, eta, rho, rho_jer) = (%.4f, %.4f, %.4f, %.4f) ' % (jet.pt(), jet.eta(), self.rho, rho_jer) 
        else:
            if self.debug: print '[Info] Impossible to smear this jet -- check if genJet (%r), and factor>1.0 (%r) '% (hasattr(jet, 'matchedGenJet'), (factor > 1.0))

        # prepare the variable to compute corrected MET after JER correction:
        jet.deltaMetFromJetSmearing = [ -(ptScale-1)*jet.px(), -(ptScale-1)*jet.py() ]

        if ptScale!=0: 
            jet.scaleEnergy(ptScale)            
        else:
            pass
        setattr(jet, "corrJER", ptScale )
        
        if self.debug and factor > 1.0 and not hasattr(jet, 'matchedGenJet'):
            print '[Debug] I am scaling jet: jer_sf = %.4f, jer_jetpt = %.4f' % (jet.jerfactor, jet.pt() )
            print '[Debug] I am smearing jet: jer_sf = %.4f, res = %.4f, jer_jetpt = %.4f' % (jet.jerfactor, jet.resolution, jet.pt() )
            
setattr(XZZJetAnalyzer,"defaultConfig", cfg.Analyzer(
    debug=False,    
    class_object = XZZJetAnalyzer,
    jetCol = 'slimmedJets',
    copyJetsByValue = False,      #Whether or not to copy the input jets or to work with references (should be 'True' if JetAnalyzer is run more than once)
    genJetCol = 'slimmedGenJets',
    rho = ('fixedGridRhoFastjetAll','',''),
    jetPt = 25.,
    jetEta = 4.7,
    jetEtaCentral = 2.4,
    jetLepDR = 0.4,
    jetLepArbitration = (lambda jet,lepton : lepton), # you can decide which to keep in case of overlaps; e.g. if the jet is b-tagged you might want to keep the jet
    cleanSelectedLeptons = True, #Whether to clean 'selectedLeptons' after disambiguation. Treat with care (= 'False') if running Jetanalyzer more than once
    minLepPt = 10,
    lepSelCut = lambda lep : True,
    relaxJetId = False,  
    doPuId = False, # Not commissioned in 7.0.X
    doQG = False, 
    checkLeptonPFOverlap = True,
    recalibrateJets = False,
    applyL2L3Residual = 'Data', # if recalibrateJets, apply L2L3Residual to Data only
    recalibrationType = "AK4PFchs",
    shiftJEC = 0, # set to +1 or -1 to apply +/-1 sigma shift to the nominal jet energies
    addJECShifts = False, # if true, add  "corr", "corrJECUp", and "corrJECDown" for each jet (requires uncertainties to be available!)
    smearJets =False,
    shiftJER = 0, # set to +1 or -1 to get +/-1 sigma shifts    
    jecPath = "",
    jerPath = "",
    calculateSeparateCorrections = False,
    calculateType1METCorrection  = False,
    type1METParams = { 'jetPtThreshold':15., 'skipEMfractionThreshold':0.9, 'skipMuons':True },
    addJERShifts = 0, # add +/-1 sigma shifts to jets, intended to be used with shiftJER=0
    cleanJetsFromFirstPhoton = False,
    cleanJetsFromTaus = False,
    cleanJetsFromIsoTracks = False,
    alwaysCleanPhotons = False,
    do_mc_match=True,
    cleanGenJetsFromPhoton = False,
    attachNeutrinos = True,
    genNuSelection = lambda nu : True, #FIXME: add here check for ispromptfinalstate
    collectionPostFix = ""
    )
)
