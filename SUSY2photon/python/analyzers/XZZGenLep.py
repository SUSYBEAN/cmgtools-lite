from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsobjects.Electron import Electron
from PhysicsTools.Heppy.physicsobjects.Muon import Muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR
import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.HeppyCore.utils.deltar import * 
from PhysicsTools.Heppy.physicsutils.genutils import *
import ROOT

from ROOT import heppy

class XZZGenLep( Analyzer ):
    #----------------------------------------
    # DECLARATION OF HANDLES OF LEPTONS STUFF   
    #----------------------------------------
    def declareHandles(self):
        super(XZZGenLep, self).declareHandles()

        #leptons
        self.handles['muons'] = AutoHandle(self.cfg_ana.muons,"std::vector<pat::Muon>")            
        self.handles['electrons'] = AutoHandle(self.cfg_ana.electrons,"std::vector<pat::Electron>")            
        # packedCandidates
        self.handles['packedCandidates'] = AutoHandle( self.cfg_ana.packedCandidates, 'std::vector<pat::PackedCandidate>') 
        #rho for muons
        self.handles['rhoMu'] = AutoHandle( self.cfg_ana.rhoMuon, 'double')
        #rho for electrons
        self.handles['rhoEle'] = AutoHandle( self.cfg_ana.rhoElectron, 'double')
        # effective area
        self.ele_effectiveAreas = getattr(self.cfg_ana, 'ele_effectiveAreas', "Spring15_25ns_v1")
        self.mu_effectiveAreas  = getattr(self.cfg_ana, 'mu_effectiveAreas',  "Spring15_25ns_v1")

        self.miniIsolationPUCorr = self.cfg_ana.miniIsolationPUCorr
        if self.miniIsolationPUCorr == "weights":
            self.IsolationComputer = heppy.IsolationComputer(0.4)
        else:
            self.IsolationComputer = heppy.IsolationComputer()

    def checkgen(self,lep):
        try:
            mom=lep.physObj.genParticle()
            while mom.pdgId()!=23:
                mom=mom.mother()
            return True
        except:
            return False
        
    #------------------
    # MAKE LEPTON LISTS
    #------------------
    def makeLeptons(self, event):
        self.IsolationComputer.setPackedCandidates(self.handles['packedCandidates'].product())
        event.selectedLeptons = self.makeAllMuons(event)+self.makeAllElectrons(event)
        event.selectedLeptons.sort(key = lambda l : l.pt(), reverse = True)
        if self.cfg_comp.isMC:
            for i in event.selectedLeptons: i.hasgen=i.physObj.genParticle().__nonzero__()
    def makeAllMuons(self, event):
        """
               make a list of all muons, and apply basic corrections to them
        """
        allmuons = map( Muon, self.handles['muons'].product() )

        # Attach EAs for isolation:
        for mu in allmuons:
          mu.rho = float(self.handles['rhoMu'].product()[0])
          if self.mu_effectiveAreas == "Spring15_25ns_v1":
              aeta = abs(mu.eta())
              if   aeta < 0.800: mu.EffectiveArea03 = 0.0735
              elif aeta < 1.300: mu.EffectiveArea03 = 0.0619
              elif aeta < 2.000: mu.EffectiveArea03 = 0.0465
              elif aeta < 2.200: mu.EffectiveArea03 = 0.0433
              elif aeta < 2.500: mu.EffectiveArea03 = 0.0577
              else:              mu.EffectiveArea03 = 0.0000
          else: 
              raise RuntimeError,  "Unsupported value for mu_effectiveAreas, can use Spring15_25ns_v1"
 
        # calculate miniIso
        for mu in allmuons:
            self.attachMiniIsolation(mu)
            mu.xdaughter=self.checkgen(mu)

        # Attach the vertex to them, for dxy/dz calculation
        for mu in allmuons:
            mu.associatedVertex = event.goodVertices[0] if len(event.goodVertices)>0 else event.vertices[0]

        # define muon id
        for mu in allmuons:
            mu.highPtID = mu.physObj.isHighPtMuon(mu.associatedVertex)
            mu.trackerHighPtID = mu.physObj.isTrackerMuon() \
                         and mu.physObj.track().isNonnull() \
                         and mu.physObj.innerTrack().isNonnull() \
                         and mu.physObj.numberOfMatchedStations()>1 \
                         and (mu.physObj.muonBestTrack().ptError()/mu.physObj.muonBestTrack().pt())<0.3 \
                         and abs(mu.physObj.muonBestTrack().dxy(mu.associatedVertex.position()))<0.2 \
                         and abs(mu.physObj.muonBestTrack().dz(mu.associatedVertex.position()))<0.5 \
                         and mu.physObj.innerTrack().hitPattern().numberOfValidPixelHits()>0 \
                         and mu.physObj.track().hitPattern().trackerLayersWithMeasurement()>5 
        for mu in allmuons:
            mu.trackerIso=mu.physObj.isolationR03().sumPt
            minc=[i for i in allmuons if i != mu and (i.highPtID or i.trackerHighPtID ) and deltaR(mu.eta(),mu.phi(),i.eta(),i.phi())<0.3]
            mu.nminc=len(minc)
            for i in minc:
                if i.physObj.innerTrack().isNonnull():mu.trackerIso-=i.physObj.innerTrack().pt()
            mu.trackerIso/=mu.pt()
        return allmuons


    def makeAllElectrons(self, event):
        """
               make a list of all electrons, and apply basic corrections to them
        """
        allelectrons = map( Electron, self.handles['electrons'].product() )

        # fill EA for rho-corrected isolation
        for ele in allelectrons:
          ele.rho = float(self.handles['rhoEle'].product()[0])
          if self.ele_effectiveAreas == "Spring15_25ns_v1":
              SCEta = abs(ele.superCluster().eta())
              ## ----- https://github.com/ikrav/cmssw/blob/egm_id_747_v2/RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt
              if   SCEta < 1.000: ele.EffectiveArea03 = 0.1752
              elif SCEta < 1.479: ele.EffectiveArea03 = 0.1862
              elif SCEta < 2.000: ele.EffectiveArea03 = 0.1411
              elif SCEta < 2.200: ele.EffectiveArea03 = 0.1534
              elif SCEta < 2.300: ele.EffectiveArea03 = 0.1903
              elif SCEta < 2.400: ele.EffectiveArea03 = 0.2243
              elif SCEta < 5.000: ele.EffectiveArea03 = 0.2687
              else:               ele.EffectiveArea03 = 0.0000
          else: 
              raise RuntimeError,  "Unsupported value for ele_effectiveAreas: can only use Data2012 (rho: ?), Phys14_v1 and Spring15_v1 (rho: fixedGridRhoFastjetAll)"

        # calculate miniIso and relisoea
        for ele in allelectrons:
            self.attachMiniIsolation(ele)
#            ele.relIsoea03=ele.absIsoWithFSR(0.3)/ele.pt()
#            if abs(ele.physObj.superCluster().eta())<1.479:
#                ele.looseiso=True if ele.relIsoea03<0.0893 else False
#            else:
#                ele.looseiso=True if ele.relIsoea03<0.121 else False
            ele.xdaughter=self.checkgen(ele)
        # Attach the vertex
        for ele in allelectrons:
            ele.associatedVertex = event.goodVertices[0] if len(event.goodVertices)>0 else event.vertices[0]

        # define electron ID and loose id no iso
        for ele in allelectrons:
            ele.loose_nonISO=ele.electronID("POG_Cuts_ID_full5x5_SPRING15_25ns_v1_ConvVetoDxyDz_Loose")
            ele.heepV60_noISO_EB = ele.pt()>35.0 \
                         and abs(ele.superCluster().eta())<1.4442 \
                         and abs(ele.deltaEtaSeedClusterTrackAtVtx())<0.004 \
                         and abs(ele.deltaPhiSuperClusterTrackAtVtx())<0.06 \
                         and (ele.full5x5_e5x5()>0 and ( ele.full5x5_e2x5Max()>ele.full5x5_e5x5()*0.94 or ele.full5x5_e1x5()>ele.full5x5_e5x5()*0.83)) \
                         and (ele.hadronicOverEm()*ele.superCluster().energy()<1.0+0.05*ele.superCluster().energy()) \
                         and abs(ele.gsfTrack().dxy(ele.associatedVertex.position()))<0.02 \
                         and (ele.gsfTrack().hitPattern().numberOfHits(ROOT.reco.HitPattern.MISSING_INNER_HITS)<=1) \
                         and ele.ecalDriven() 

            ele.heepV60_noISO_EE = ele.pt()>35.0 \
                         and abs(ele.superCluster().eta())>1.566 and abs(ele.superCluster().eta())<2.5 \
                         and abs(ele.deltaEtaSeedClusterTrackAtVtx())<0.006 \
                         and abs(ele.deltaPhiSuperClusterTrackAtVtx())<0.06 \
                         and ele.full5x5_sigmaIetaIeta()<0.03 \
                         and (ele.hadronicOverEm()*ele.superCluster().energy()<5.0+0.05*ele.superCluster().energy()) \
                         and abs(ele.gsfTrack().dxy(ele.associatedVertex.position()))<0.05 \
                         and (ele.gsfTrack().hitPattern().numberOfHits(ROOT.reco.HitPattern.MISSING_INNER_HITS)<=1) \
                         and ele.ecalDriven() 

            ele.heepV60_noISO = ele.heepV60_noISO_EB or ele.heepV60_noISO_EE
        return allelectrons 


    def attachMiniIsolation(self, lep):
        lep.miniIsoR = 10.0/min(max(lep.pt(), 50.0),200.0)
        what = "mu" if (abs(lep.pdgId()) == 13) else ("eleB" if lep.isEB() else "eleE")

        if what=="mu":
            lep.miniAbsIsoChargedHad = self.IsolationComputer.chargedHadAbsIso(lep.physObj, lep.miniIsoR, 0.0001, 0.0, self.IsolationComputer.selfVetoNone);
        else:
            lep.miniAbsIsoChargedHad = self.IsolationComputer.chargedHadAbsIso(lep.physObj, lep.miniIsoR, {"mu":0.0001,"eleB":0,"eleE":0.015}[what], 0.0, self.IsolationComputer.selfVetoNone);

        if self.miniIsolationPUCorr == None: puCorr = self.cfg_ana.mu_isoCorr if what=="mu" else self.cfg_ana.ele_isoCorr
        else: puCorr = self.miniIsolationPUCorr

        if puCorr == "weights":
            if what == "mu":
                lep.miniAbsIsoNeutral = self.IsolationComputer.neutralAbsIsoWeighted(lep.physObj, lep.miniIsoR, 0.01, 0.5,self.IsolationComputer.selfVetoNone);
            else:
                lep.miniAbsIsoNeutral = ( self.IsolationComputer.photonAbsIsoWeighted(lep.physObj, lep.miniIsoR, 0.08 if what=="eleE" else 0.0, 0.0, self.IsolationComputer.selfVetoNone) + self.IsolationComputer.neutralHadAbsIsoWeighted(lep.physObj, lep.miniIsoR, 0.0, 0.0, self.IsolationComputer.selfVetoNone) )
        else:
            if what == "mu":
                #lep.miniAbsIsoNeutral = self.IsolationComputer.neutralAbsIsoRaw(lep.physObj, lep.miniIsoR, 0.01, 0.5);
                lep.miniAbsIsoPho  = self.IsolationComputer.photonAbsIsoRaw(lep.physObj, lep.miniIsoR, 0.01, 0.5,self.IsolationComputer.selfVetoNone)
                lep.miniAbsIsoNeutralHad = self.IsolationComputer.neutralHadAbsIsoRaw(lep.physObj, lep.miniIsoR, 0.01, 0.5,self.IsolationComputer.selfVetoNone)
            else:
                lep.miniAbsIsoPho  = self.IsolationComputer.photonAbsIsoRaw(lep.physObj, lep.miniIsoR, 0.08 if what == "eleE" else 0.0, 0.0, self.IsolationComputer.selfVetoNone)
                lep.miniAbsIsoNeutralHad = self.IsolationComputer.neutralHadAbsIsoRaw(lep.physObj, lep.miniIsoR, 0.0, 0.0, self.IsolationComputer.selfVetoNone)
            # only calculate rhoArea
            lep.miniAbsIsoNeutral = lep.miniAbsIsoPho + lep.miniAbsIsoNeutralHad
            lep.miniAbsIsoNeutral = max(0.0, lep.miniAbsIsoNeutral - lep.rho * lep.EffectiveArea03 * (lep.miniIsoR/0.3)**2)

        if lep.pt()<5.0:
            lep.miniAbsIso = 9999.0
            lep.miniRelIso = 9999.0
        else:
            lep.miniAbsIso = lep.miniAbsIsoChargedHad + lep.miniAbsIsoNeutral
            lep.miniRelIso = lep.miniAbsIso/lep.pt()

    def process(self, event):
        self.readCollections( event.input )
        self.makeLeptons(event)
        if not event.selectedLeptons:
            return False
        else:
            return True

