#!/bin/env python
from math import *
import ROOT
#from CMGTools.TTHAnalysis.signedSip import *
from PhysicsTools.Heppy.analyzers.core.autovars import *
from PhysicsTools.HeppyCore.utils.deltar import deltaR

objectFloat = NTupleObjectType("builtInType", variables = [
    NTupleVariable("",    lambda x : x),
])
objectInt = NTupleObjectType("builtInType", variables = [
    NTupleVariable("",    lambda x : x,int),
])

twoVectorType = NTupleObjectType("twoVector", variables = [
    NTupleVariable("pt",    lambda x : x.pt()),
    NTupleVariable("phi",   lambda x : x.phi()),
])

threeVectorType = NTupleObjectType("threeVector", variables = [
    NTupleVariable("x",    lambda x : x.x()),
    NTupleVariable("y",    lambda x : x.y()),
    NTupleVariable("z",    lambda x : x.z()),
])

fourVectorType = NTupleObjectType("fourVector", variables = [
    NTupleVariable("pt",    lambda x : x.pt()),
    #NTupleVariable("px",    lambda x : x.px()),
    #NTupleVariable("py",    lambda x : x.py()),
    #NTupleVariable("pz",    lambda x : x.pz()),
    NTupleVariable("eta",   lambda x : x.eta()),
    NTupleVariable("rapidity",   lambda x : x.rapidity()),
    NTupleVariable("phi",   lambda x : x.phi()),
    NTupleVariable("mass",  lambda x : x.mass()),
    NTupleVariable("p4",    lambda x : x, "TLorentzVector", default=ROOT.reco.Particle.LorentzVector(0.,0.,0.,0.), filler = lambda vector, obj: vector.SetPtEtaPhiM(obj.pt(), obj.eta(), obj.phi(), obj.mass())),
    #               ^^^^------- Note: p4 normally is not saved unless 'saveTLorentzVectors' is enabled in the tree producer
])
tlorentzFourVectorType = NTupleObjectType("tlorentzFourVectorType", variables = [
    NTupleVariable("pt",    lambda x : x.Pt()),
    #NTupleVariable("px",    lambda x : x.Px()),
    #NTupleVariable("py",    lambda x : x.Py()),
    #NTupleVariable("pz",    lambda x : x.Pz()),
    NTupleVariable("eta",   lambda x : x.Eta()),
    NTupleVariable("rapidity",   lambda x : x.Rapidity()),
    NTupleVariable("phi",   lambda x : x.Phi()),
    NTupleVariable("energy",  lambda x : x.E()),
    NTupleVariable("mass",  lambda x : x.M()),
    NTupleVariable("mt",  lambda x : x.Mt()),
])
particleType = NTupleObjectType("particle", baseObjectTypes = [ fourVectorType ], variables = [
    NTupleVariable("pdgId",   lambda x : x.pdgId(), int),
])

weightsInfoType = NTupleObjectType("WeightsInfo", mcOnly=True, variables = [
    NTupleVariable("id",   lambda x : x.id, int),
    NTupleVariable("wgt",   lambda x : x.wgt),
])

##------------------------------------------  
## LEPTON
##------------------------------------------  

### BASIC VERSION WITH ONLY MAIN LEPTON ID CRITERIA
leptonType = NTupleObjectType("lepton", baseObjectTypes = [ particleType ], variables = [
    NTupleVariable("charge",   lambda x : x.charge(), int),
    NTupleVariable("hasgen",  lambda x : getattr(x,"hasgen",-1), int, help="has gen particle"),
    NTupleVariable("ptErr",  lambda x : x.ptErr() if abs(x.pdgId())==13 else -999, help="pt Error"),
    #NTupleVariable("TuneP_pt",  lambda x : x.TuneP_pt() if abs(x.pdgId())==13 else x.pt(), help="TuneP Pt for muon"),
    #NTupleVariable("TuneP_ptErr",  lambda x : x.physObj.tunePMuonBestTrack().ptError() if abs(x.pdgId())==13 else -999, help="TuneP Pt error for muon"),
    #NTupleVariable("TuneP_eta",  lambda x : x.TuneP_eta() if abs(x.pdgId())==13 else x.eta(), help="TuneP eta for muon"),
    #NTupleVariable("TuneP_phi",  lambda x : x.TuneP_phi() if abs(x.pdgId())==13 else x.phi(), help="TuneP phi for muon"),
    #NTupleVariable("TuneP_type",  lambda x : x.physObj.tunePMuonBestTrackType() if abs(x.pdgId())==13 else -999, help="TuneP type for muon, https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_6_3_patch2/doc/html/df/de3/classreco_1_1Muon.html#afceb985a23ee1d456e4dc91391f2e7fe"),
    # Impact parameter
    NTupleVariable("dxy",   lambda x : x.dxy(), help="d_{xy} with respect to PV, in cm (with sign)"),
    NTupleVariable("dz",    lambda x : x.dz() , help="d_{z} with respect to PV, in cm (with sign)"),
    NTupleVariable("edxy",  lambda x : x.edB() if hasattr(x,'edB') else -999, help="#sigma(d_{xy}) with respect to PV, in cm"),
    NTupleVariable("edz",   lambda x : x.edz() if hasattr(x,'edz') else -999, help="#sigma(d_{z}) with respect to PV, in cm"),
    NTupleVariable("ip3d",  lambda x : x.ip3D() if hasattr(x,'ip3D') else -999, help="d_{3d} with respect to PV, in cm (absolute value)"),
    # Conversion rejection
    #NTupleVariable("miniRelIso",  lambda x : x.miniRelIso if hasattr(x,'miniRelIso') else  -999, help="PF Rel miniRel, pile-up corrected"),
    NTupleVariable("muonincore03",  lambda x : getattr(x,"nminc",-999), help="number of other muons in core0.3 "),
    NTupleVariable("trackerIso",  lambda x : getattr(x,"trackerIso",-999), help="muon tracker isolation"),
    # id and iso
    #NTupleVariable("heepV60_noISO",  lambda x : x.heepV60_noISO if hasattr(x,'heepV60_noISO') else  -999, help="heepV60_noISO"),
    NTupleVariable("highPtID",  lambda x : x.highPtID if hasattr(x,'highPtID') else  -999, help="highPtID"),
    NTupleVariable("trackerHighPtID",  lambda x : x.trackerHighPtID if hasattr(x,'trackerHighPtID') else  -999, help="trackerHighPtID"),
    NTupleVariable("isfromX", lambda x : x.xdaughter if hasattr(x,'xdaughter') else -100, help="is from X"),
    # Extra electron ID working points
    NTupleVariable("looseelectron",   lambda x : x.physObj.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose") if abs(x.pdgId())==11 else -999, int, help="electron POG Loose id"),
    NTupleVariable("eSCeta",   lambda x : x.physObj.superCluster().eta() if abs(x.pdgId())==11 else -999, float, help="electron SC eta"),
    NTupleVariable("looseelectronnoiso",   lambda x : x.loose_nonISO if abs(x.pdgId())==11 else -999, int, help="electron POG Loose id no iso"),
    NTupleVariable("looseisoelectron",   lambda x : x.looseiso if abs(x.pdgId())==11 else -999, int, help="electron POG Loose id default iso"),
    NTupleVariable("electronrelIsoea03",   lambda x : x.relIsoea03 if abs(x.pdgId())==11 else -999, float, help="electron relisoea03"),
    #NTupleVariable("lepsf",  lambda x : getattr(x,'lepsf',1), help="lepton sf"),
    #NTupleVariable("lepsfUp",  lambda x : getattr(x,'lepsfUp',1), help="lepton sf upper"),
    #NTupleVariable("lepsfLo",  lambda x : getattr(x,'lepsfLo',1), help="lepton sf lower"),
    # PF candidates info
    #NTupleVariable("npf",  lambda x : x.numberOfSourceCandidatePtrs(), help="n pfparticles"),
    #NTupleVariable("pf0_pt",  lambda x : x.sourceCandidatePtr(0).pt() if x.numberOfSourceCandidatePtrs()>0 else -999 , help="lepton pf0 pt"),
    #NTupleVariable("pf0_eta",  lambda x : x.sourceCandidatePtr(0).eta() if x.numberOfSourceCandidatePtrs()>0 else -999 , help="lepton pf0 eta"),
    #NTupleVariable("pf0_phi",  lambda x : x.sourceCandidatePtr(0).phi() if x.numberOfSourceCandidatePtrs()>0 else -999 , help="lepton pf0 phi"),
    #NTupleVariable("pf1_pt",  lambda x : x.sourceCandidatePtr(1).pt() if x.numberOfSourceCandidatePtrs()>1 else -999 , help="lepton pf1 pt"),
    #NTupleVariable("pf1_eta",  lambda x : x.sourceCandidatePtr(1).eta() if x.numberOfSourceCandidatePtrs()>1 else -999 , help="lepton pf1 eta"),
    #NTupleVariable("pf1_phi",  lambda x : x.sourceCandidatePtr(1).phi() if x.numberOfSourceCandidatePtrs()>1 else -999 , help="lepton pf1 phi"),
    NTupleVariable("trigerob_HLTbit", lambda x : x.triggerbit if hasattr(x,'triggerbit')  else -100, int, help="Electron matched HLT object path bit"),
])

### EXTENDED VERSION WITH INDIVIUAL DISCRIMINATING VARIABLES
leptonTypeExtra = NTupleObjectType("leptonExtra", baseObjectTypes = [ leptonType ], variables = [
    # Extra isolation variables
    NTupleVariable("chargedHadRelIso03",  lambda x : x.chargedHadronIsoR(0.3)/x.pt(), help="PF Rel Iso, R=0.3, charged hadrons only"),
    NTupleVariable("chargedHadRelIso04",  lambda x : x.chargedHadronIsoR(0.4)/x.pt(), help="PF Rel Iso, R=0.4, charged hadrons only"),
    # Extra muon ID working points
    NTupleVariable("softMuonId", lambda x : x.muonID("POG_ID_Soft") if abs(x.pdgId())==13 else 1, int, help="Muon POG Soft id"),
    NTupleVariable("pfMuonId",   lambda x : x.muonID("POG_ID_Loose") if abs(x.pdgId())==13 else 1, int, help="Muon POG Loose id"),
    NTupleVariable("tightMuonId",   lambda x : x.muonID("POG_ID_Tight") if abs(x.pdgId())==13 else 1, int, help="Muon POG Tight id"),
    # Extra electron ID working points
    NTupleVariable("eleCutId2012_full5x5",     lambda x : (1*x.electronID("POG_Cuts_ID_2012_full5x5_Veto") + 1*x.electronID("POG_Cuts_ID_2012_full5x5_Loose") + 1*x.electronID("POG_Cuts_ID_2012_full5x5_Medium") + 1*x.electronID("POG_Cuts_ID_2012_full5x5_Tight")) if abs(x.pdgId()) == 11 else -1, int, help="Electron cut-based id (POG 2012, full5x5 shapes): 0=none, 1=veto, 2=loose, 3=medium, 4=tight"),
    # Extra tracker-related variables
    NTupleVariable("trackerLayers", lambda x : (x.track() if abs(x.pdgId())==13 else x.gsfTrack()).hitPattern().trackerLayersWithMeasurement(), int, help="Tracker Layers"),
    NTupleVariable("pixelLayers", lambda x : (x.track() if abs(x.pdgId())==13 else x.gsfTrack()).hitPattern().pixelLayersWithMeasurement(), int, help="Pixel Layers"),
    NTupleVariable("trackerHits", lambda x : (x.track() if abs(x.pdgId())==13 else x.gsfTrack()).hitPattern().numberOfValidTrackerHits(), int, help="Tracker hits"),
    NTupleVariable("lostOuterHits",    lambda x : (x.gsfTrack() if abs(x.pdgId())==11 else x.innerTrack()).hitPattern().numberOfLostHits(ROOT.reco.HitPattern.MISSING_OUTER_HITS), int, help="Number of lost hits on inner track"),
    NTupleVariable("innerTrackValidHitFraction", lambda x : (x.gsfTrack() if abs(x.pdgId())==11 else x.innerTrack()).validFraction(), help="fraction of valid hits on inner track"), 
    NTupleVariable("innerTrackChi2",      lambda x : (x.gsfTrack() if abs(x.pdgId())==11 else x.innerTrack()).normalizedChi2(), help="Inner track normalized chi2"), 
    # Extra muon ID variables
    NTupleVariable("nStations",    lambda lepton : lepton.numberOfMatchedStations() if abs(lepton.pdgId()) == 13 else 4, help="Number of matched muons stations (4 for electrons)"),
    NTupleVariable("caloCompatibility",      lambda lepton : lepton.caloCompatibility() if abs(lepton.pdgId()) == 13 else 0, help="Calorimetric compatibility"), 
    NTupleVariable("globalTrackChi2",      lambda lepton : lepton.globalTrack().normalizedChi2() if abs(lepton.pdgId()) == 13 and lepton.globalTrack().isNonnull() else 0, help="Global track normalized chi2"), 
    NTupleVariable("trkKink",      lambda lepton : lepton.combinedQuality().trkKink if abs(lepton.pdgId()) == 13 else 0, help="Tracker kink-finder"), 
    NTupleVariable("segmentCompatibility", lambda lepton : lepton.segmentCompatibility() if abs(lepton.pdgId()) == 13 else 0, help="Segment-based compatibility"), 
    NTupleVariable("chi2LocalPosition",    lambda lepton : lepton.combinedQuality().chi2LocalPosition if abs(lepton.pdgId()) == 13 else 0, help="Tracker-Muon matching in position"), 
    NTupleVariable("chi2LocalMomentum",    lambda lepton : lepton.combinedQuality().chi2LocalMomentum if abs(lepton.pdgId()) == 13 else 0, help="Tracker-Muon matching in momentum"), 
    NTupleVariable("glbTrackProbability",  lambda lepton : lepton.combinedQuality().glbTrackProbability if abs(lepton.pdgId()) == 13 else 0, help="Global track pseudo-probability"), 
    # Extra electron ID variables
    NTupleVariable("sigmaIEtaIEta",  lambda x : x.full5x5_sigmaIetaIeta() if abs(x.pdgId())==11 else 0, help="Electron sigma(ieta ieta), with full5x5 cluster shapes"),
    NTupleVariable("dEtaScTrkIn",    lambda x : x.deltaEtaSuperClusterTrackAtVtx() if abs(x.pdgId())==11 else 0, help="Electron deltaEtaSuperClusterTrackAtVtx (without absolute value!)"),
    NTupleVariable("dPhiScTrkIn",    lambda x : x.deltaPhiSuperClusterTrackAtVtx() if abs(x.pdgId())==11 else 0, help="Electron deltaPhiSuperClusterTrackAtVtx (without absolute value!)"),
    NTupleVariable("hadronicOverEm", lambda x : x.hadronicOverEm() if abs(x.pdgId())==11 else 0, help="Electron hadronicOverEm"),
    NTupleVariable("eInvMinusPInv",  lambda x : ((1.0/x.ecalEnergy() - x.eSuperClusterOverP()/x.ecalEnergy()) if x.ecalEnergy()>0. else 9e9) if abs(x.pdgId())==11 else 0, help="Electron 1/E - 1/p  (without absolute value!)"),
#    new version used by EGM in Spring15, 7_4_14:
    NTupleVariable("eInvMinusPInv_tkMom", lambda x: ((1.0/x.ecalEnergy()) - (1.0 / x.trackMomentumAtVtx().R() ) if (x.ecalEnergy()>0. and x.trackMomentumAtVtx().R()>0.) else 9e9) if abs(x.pdgId())==11 else 0, help="Electron 1/E - 1/p_tk_vtx  (without absolute value!)"),
    NTupleVariable("etaSc", lambda x : x.superCluster().eta() if abs(x.pdgId())==11 else -100, help="Electron supercluster pseudorapidity"),
#    NTupleVariable("trigerob_pt", lambda x : x.triggerob.pt if hasattr(x,'triggerob') else -100, help="Electron matched HLT object pt"),
#    NTupleVariable("trigerob_eta", lambda x : x.triggerob.eta if hasattr(x,'triggerob') else -100, help="Electron matched HLT object eta"),
#    NTupleVariable("trigerob_phi", lambda x : x.triggerob.phi if hasattr(x,'triggerob') else -100, help="Electron matched HLT object phi"),
#    NTupleVariable("trigerob_deltaR", lambda x : x.triggerob.dR if hasattr(x,'triggerob') else -100, help="Electron matched HLT object phi"),
])
 

##------------------------------------------  
## TAU
##------------------------------------------  

tauType = NTupleObjectType("tau",  baseObjectTypes = [ particleType ], variables = [
    NTupleVariable("charge",   lambda x : x.charge(), int),
    NTupleVariable("decayMode",   lambda x : x.decayMode(), int),
    NTupleVariable("idDecayMode",   lambda x : x.idDecayMode, int),
    NTupleVariable("idDecayModeNewDMs",   lambda x : x.idDecayModeNewDMs, int),
    NTupleVariable("dxy",   lambda x : x.dxy(), help="d_{xy} of lead track with respect to PV, in cm (with sign)"),
    NTupleVariable("dz",    lambda x : x.dz() , help="d_{z} of lead track with respect to PV, in cm (with sign)"),
    NTupleVariable("idMVA", lambda x : x.idMVA, int, help="1,2,3,4,5,6 if the tau passes the very loose to very very tight WP of the MVA3oldDMwLT discriminator"),
    NTupleVariable("idMVANewDM", lambda x : x.idMVANewDM, int, help="1,2,3,4,5,6 if the tau passes the very loose to very very tight WP of the MVA3newDMwLT discriminator"),
    NTupleVariable("idCI3hit", lambda x : x.idCI3hit, int, help="1,2,3 if the tau passes the loose, medium, tight WP of the By<X>CombinedIsolationDBSumPtCorr3Hits discriminator"),
    NTupleVariable("idAntiMu", lambda x : x.idAntiMu, int, help="1,2 if the tau passes the loose/tight WP of the againstMuon<X>3 discriminator"),
    NTupleVariable("idAntiE", lambda x : x.idAntiE, int, help="1,2,3,4,5 if the tau passes the v loose, loose, medium, tight, v tight WP of the againstElectron<X>MVA5 discriminator"),
    NTupleVariable("isoCI3hit",  lambda x : x.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"), help="byCombinedIsolationDeltaBetaCorrRaw3Hits raw output discriminator"),
    # MC-match info
    NTupleVariable("mcMatchId",  lambda x : getattr(x, 'mcMatchId', -99), int, mcOnly=True, help="Match to source from hard scatter (pdgId of heaviest particle in chain, 25 for H, 6 for t, 23/24 for W/Z), zero if non-prompt or fake"),
])

##------------------------------------------  
##  ISOTRACK
##------------------------------------------  

isoTrackType = NTupleObjectType("isoTrack",  baseObjectTypes = [ particleType ], variables = [
    NTupleVariable("charge",   lambda x : x.charge(), int),
    NTupleVariable("dz",    lambda x : x.dz() , help="d_{z} of lead track with respect to PV, in cm (with sign)"),
    NTupleVariable("absIso",  lambda x : x.absIso, float, mcOnly=False, help="abs charged iso with condition for isolation such that Min(0.2*pt, 8 GeV)"),
    NTupleVariable("relIsoAn04",  lambda x : x.relIsoAn04 if hasattr(x,'relIsoAn04') else  -999, help="PF Activity Annulus, pile-up corrected"),
    NTupleVariable("mcMatchId",  lambda x : getattr(x, 'mcMatchId', -99), int, mcOnly=True, help="Match to source from hard scatter (pdgId of heaviest particle in chain, 25 for H, 6 for t, 23/24 for W/Z), zero if non-prompt or fake"),
])


##------------------------------------------  
## PHOTON
##------------------------------------------  

photonType = NTupleObjectType("gamma", baseObjectTypes = [ particleType ], variables = [
    NTupleVariable("idCutBased", lambda x : x.idCutBased, int, help="1,2,3 if the gamma passes the POG_SPRING15_25ns_(Loose,Medium,Tight)"),
    NTupleVariable("hOverE",  lambda x : x.hOVERe(), float, help="hoverE for photons"),
    NTupleVariable("r9",  lambda x : x.full5x5_r9(), float, help="r9 for photons"),
    NTupleVariable("sigmaIetaIeta",  lambda x : x.full5x5_sigmaIetaIeta(), float, help="sigmaIetaIeta for photons"),
   # NTupleVariable("chHadIso04",  lambda x : x.chargedHadronIso(), float, help="chargedHadronIsolation for photons (PAT method, deltaR = 0.4)"),
    NTupleVariable("chHadIso", lambda x : x.chargedHadronIso(x.isoCorr), float, help="chargedHadronIsolation for photons with footprint removal and pile-up correction"),
    NTupleVariable("phIso", lambda x : x.photonIso(x.isoCorr), float, help="gammaIsolation for photons with footprint removal and pile-up correction"),
    NTupleVariable("neuHadIso", lambda x : x.neutralHadronIso(x.isoCorr), float, help="neutralHadronIsolation for photons with footprint removal and pile-up correction"),
    #NTupleVariable("chHadIso", lambda x : x.chargedHadronIso(), float, help="chargedHadronIsolation for photons with footprint removal"),
    #NTupleVariable("phIso", lambda x : x.photonIso(), float, help="gammaIsolation for photons with footprint removal"),
    #NTupleVariable("neuHadIso", lambda x : x.neutralHadronIso(), float, help="neutralHadronIsolation for photons with footprint removal"),
    #NTupleVariable("relIso", lambda x : x.ftprRelIso03 if hasattr(x,'ftprRelIso03') else x.relIso, float, help="relativeIsolation for photons with footprint removal and pile-up correction"),
    NTupleVariable("mcMatchId",  lambda x : getattr(x, 'mcMatchId', -99), int, mcOnly=True, help="Match to source from hard scatter (pdgId of heaviest particle in chain, 25 for H, 6 for t, 23/24 for W/Z), zero if non-prompt or fake"),
    NTupleVariable("mcPt",   lambda x : x.mcGamma.pt() if getattr(x,"mcGamma",None) else 0., mcOnly=True, help="p_{T} of associated gen photon"),
    NTupleVariable("mcEta",   lambda x : x.mcGamma.eta() if getattr(x,"mcGamma",None) else 0., mcOnly=True, help="eta of associated gen photon"),
    NTupleVariable("mcPhi",   lambda x : x.mcGamma.phi() if getattr(x,"mcGamma",None) else 0., mcOnly=True, help="phi of associated gen photon"),
])

##------------------------------------------  
## JET
##------------------------------------------  

jetType = NTupleObjectType("jet",  baseObjectTypes = [ fourVectorType ], variables = [
    NTupleVariable("puId", lambda x : getattr(x, 'puJetIdPassed', -99), int,     mcOnly=False, help="puId (full MVA, loose WP, 5.3.X training on AK5PFchs: the only thing that is available now)"),
    NTupleVariable("rawPt",  lambda x : x.pt() * x.rawFactor(), help="p_{T} before JEC"),
    NTupleVariable("hadronFlavour", lambda x : x.hadronFlavour(), int,     mcOnly=True, help="hadron flavour (ghost matching to B/C hadrons)"),
    NTupleVariable("mcMatchId",  lambda x : getattr(x, 'matchedGenJetIdx', -99), int, mcOnly=True, help="Match to source from hard scatter (pdgId of heaviest particle in chain, 25 for H, 6 for t, 23/24 for W/Z), zero if non-prompt or fake"),
    #NTupleVariable("corr_JECUp",  lambda x : getattr(x, 'corrJECUp', -99), float, mcOnly=True, help="JEC correction factor up"),
    #NTupleVariable("corr_JECDown",  lambda x : getattr(x, 'corrJECDown', -99), float, mcOnly=True, help="JEC correction factor down"),
    #NTupleVariable("corr",  lambda x : getattr(x, 'corr', -99), float, mcOnly=True, help="JEC correction factor"),
])

genJetType = NTupleObjectType("genJetType",  baseObjectTypes = [ fourVectorType ], variables = [
    NTupleVariable("emEnergy", lambda x : x.emEnergy(), float, mcOnly=True, help="xxx"),
    NTupleVariable("hadEnergy", lambda x : x.hadEnergy(), float, mcOnly=True, help="xxx"),
    NTupleVariable("invisibleEnergy", lambda x : x.invisibleEnergy(), float, mcOnly=True, help="xxx"),
    NTupleVariable("auxiliaryEnergy", lambda x : x.auxiliaryEnergy(), float, mcOnly=True, help="xxx"),
    NTupleVariable("nGenConstituents", lambda x : x.numberOfDaughters(), int, mcOnly=True, help="xxx"),

])
jetTypeExtra = NTupleObjectType("jetExtra",  baseObjectTypes = [ jetType ], variables = [
    NTupleVariable("area",   lambda x : x.jetArea(), help="Catchment area of jet"),
    # QG variables:
    NTupleVariable("qgl",   lambda x :x.qgl() , float, mcOnly=False,help="QG Likelihood"),
    NTupleVariable("ptd",   lambda x : getattr(x.computeQGvars(),'ptd', 0), float, mcOnly=False,help="QG input variable: ptD"),
    NTupleVariable("axis2",   lambda x : getattr(x.computeQGvars(),'axis2', 0) , float, mcOnly=False,help="QG input variable: axis2"),
    NTupleVariable("mult",   lambda x : getattr(x.computeQGvars(),'mult', 0) , int, mcOnly=False,help="QG input variable: total multiplicity"),
    NTupleVariable("partonId", lambda x : getattr(x,'partonId', 0), int,     mcOnly=True, help="parton flavour (manually matching to status 23 particles)"),
    NTupleVariable("partonMotherId", lambda x : getattr(x,'partonMotherId', 0), int,     mcOnly=True, help="parton flavour (manually matching to status 23 particles)"),
    NTupleVariable("nLeptons",   lambda x : len(x.leptons) if  hasattr(x,'leptons') else  0 , float, mcOnly=False,help="Number of associated leptons"),
])

      
##------------------------------------------  
## MET
##------------------------------------------  
  
metType = NTupleObjectType("met", baseObjectTypes = [ fourVectorType ], variables = [
    NTupleVariable("sumEt", lambda x : x.sumEt() ),
    NTupleVariable("rawPt",  lambda x : x.uncorPt() ),
    NTupleVariable("rawPhi", lambda x : x.uncorPhi() ),
    NTupleVariable("rawSumEt", lambda x : x.uncorSumEt() ),
    NTupleVariable("genPt",  lambda x : x.genMET().pt() if x.genMET() else 0 , mcOnly=True ),
    NTupleVariable("genPhi", lambda x : x.genMET().phi() if x.genMET() else 0, mcOnly=True ),
    NTupleVariable("genEta", lambda x : x.genMET().eta() if x.genMET() else 0, mcOnly=True ),
    #NTupleVariable("metSig", lambda x : x.metSignificance() ),
    #NTupleSubObject("metNoJet",  lambda x : x.metNoJet, fourVectorType),
])

##------------------------------------------  
## GENPARTICLE
##------------------------------------------  

genParticleType = NTupleObjectType("genParticle", baseObjectTypes = [ particleType ], mcOnly=True, variables = [
    NTupleVariable("charge",   lambda x : x.threeCharge()/3.0, float),
    NTupleVariable("status",   lambda x : x.status(),int),
])
genParticleWithMotherInfo = NTupleObjectType("genParticleWithMotherInfo", baseObjectTypes = [ genParticleType ], mcOnly=True, variables = [
    NTupleVariable("isLastCopyBeforeFSR", lambda x : x.isLastCopyBeforeFSR(), int, help="isLastCopyBeforeFSR"),
    NTupleVariable("motherId", lambda x : x.mother(0).pdgId() if x.mother(0) else 0, int, help="pdgId of the mother of the particle"),
    NTupleVariable("motherStatus", lambda x : x.mother(0).status() if x.mother(0) else 0, int, help="status of the mother of the particle"),
    NTupleVariable("motherisLastCopyBeforeFSR", lambda x : x.mother(0).isLastCopyBeforeFSR() if x.mother(0) else 0, int, help="isLastCopyBeforeFSR of the mother of the particle"),
    NTupleVariable("motherPt", lambda x : x.mother(0).pt() if x.mother(0) else 0, float, help="pt of the mother of the particle"),
    NTupleVariable("motherEta", lambda x : x.mother(0).eta() if x.mother(0) else 0, float, help="eta of the mother of the particle"),
    NTupleVariable("motherPhi", lambda x : x.mother(0).phi() if x.mother(0) else 0, float, help="phi of the mother of the particle"),
    NTupleVariable("grandmotherId", lambda x : x.mother(0).mother(0).pdgId() if x.mother(0) and x.mother(0).mother(0) else 0, int, help="pdgId of the grandmother of the particle"),
    NTupleVariable("grandmotherStatus", lambda x : x.mother(0).mother(0).status() if x.mother(0) and x.mother(0).mother(0) else 0, int, help="status of the grandmother of the particle"),
    NTupleVariable("grandmotherisLastCopyBeforeFSR", lambda x : x.mother(0).mother(0).isLastCopyBeforeFSR() if x.mother(0) and x.mother(0).mother(0) else 0, int, help="isLastCopyBeforeFSR of the grandmother of the particle"),
    NTupleVariable("grandmotherPt", lambda x : x.mother(0).mother(0).pt() if x.mother(0) and x.mother(0).mother(0) else 0, float, help="pt of the grandmother of the particle"),
    NTupleVariable("grandmotherEta", lambda x : x.mother(0).mother(0).eta() if x.mother(0) and x.mother(0).mother(0) else 0, float, help="eta of the grandmother of the particle"),
    NTupleVariable("grandmotherPhi", lambda x : x.mother(0).mother(0).phi() if x.mother(0) and x.mother(0).mother(0) else 0, float, help="phi of the grandmother of the particle"),
])

genParticleWithMotherId = NTupleObjectType("genParticleWithMotherId", baseObjectTypes = [ genParticleType ], mcOnly=True, variables = [
    NTupleVariable("motherId", lambda x : x.mother(0).pdgId() if x.mother(0) else 0, int, help="pdgId of the mother of the particle"),
    NTupleVariable("grandmotherId", lambda x : x.mother(0).mother(0).pdgId() if x.mother(0) and x.mother(0).mother(0) else 0, int, help="pdgId of the grandmother of the particle")
])
genParticleWithAncestryType = NTupleObjectType("genParticleWithAncestry", baseObjectTypes = [ genParticleType ], mcOnly=True, variables = [
    NTupleVariable("motherId", lambda x : x.motherId, int, help="pdgId of the mother of the particle"),
    NTupleVariable("grandmotherId", lambda x : x.grandmotherId, int, help="pdgId of the grandmother of the particle"),
    NTupleVariable("sourceId", lambda x : x.sourceId, int, help="origin of the particle (heaviest ancestor): 6=t, 25=h, 23/24=W/Z"),
])
genParticleWithLinksType = NTupleObjectType("genParticleWithLinks", baseObjectTypes = [ genParticleWithAncestryType ], mcOnly=True, variables = [
    NTupleVariable("motherIndex", lambda x : x.motherIndex, int, help="index of the mother in the generatorSummary")
])

