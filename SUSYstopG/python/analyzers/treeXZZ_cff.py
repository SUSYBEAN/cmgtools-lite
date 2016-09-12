from CMGTools.XZZ2l2nu.analyzers.AutoFillTreeProducer  import * 
from CMGTools.XZZ2l2nu.analyzers.XZZTypes  import * 
from CMGTools.XZZ2l2nu.analyzers.Skimmer  import * 
import PhysicsTools.HeppyCore.framework.config as cfg

vvSkimmer = cfg.Analyzer(
    Skimmer,
    name='vvSkimmer',
    #required = ['LLNuNu', 'ElMuNuNu']
    required = ['LLNuNu']
)

leptonSkimmer = cfg.Analyzer(
    Skimmer,
    name='leptonSkimmer',
    required = ['inclusiveLeptons']
)

fullTreeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='fullTreeProducer',
     vectorTree = True,
     saveTLorentzVectors = True,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
     globalVariables = [
         NTupleVariable("nLL",lambda ev: len(ev.LL) , int),      
         NTupleVariable("nElMu",lambda ev: len(ev.ElMu) , int),       
         NTupleVariable("nLLNuNu",lambda ev: len(ev.LLNuNu) , int),       
         NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"), 
         NTupleVariable("triggersf",  lambda x : getattr(x,'trgsf',1), help="singleelectron/muon trigger sf"),
     ],
     globalObjects =  {
         "met" : NTupleObject("met", metType, help="PF E_{T}^{miss}, with 76X type 1 corrections"),
         "met_miniAod" : NTupleObject("met_74x", metType, help="PF E_{T}^{miss}, with 74X type 1 corrections stored in miniAOD"),
         "metNoJet" : NTupleObject("metNoJet", fourVectorType, help="PF E_{T}^{miss}, with out jet"),
         #"met_JEC" : NTupleObject("met_JEC", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC added"),
         #"met_JECUp" : NTupleObject("met_JECUp", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC up"),
         #"met_JECDown" : NTupleObject("met_JECDown", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC down"),
         #"met_JECJER" : NTupleObject("met_JECJER", metType, mcOnly=True, help="PF E_{T}^{miss}, after type 1 corrections with JEC+JER jets"),
         "sumJetsInT1": NTupleObject("sumJetsInT1", jet4metType, help="sum of Delta(px, py, et) of jets used for type1METCorr"),
     },
     collections = {
         #"LL"  : NTupleCollection("Zll",LLType,5, help="Z to ll"),
         "ElMu" : NTupleCollection("elmu",LLType,5, help="electron - muon pair for non-resonant bkg"),
         "selectedLeptons" : NTupleCollection("lep",leptonType,10, help="selected leptons"),
         "genLeptons" : NTupleCollection("genLep", genParticleType, 10, help="Generated leptons (e/mu) from W/Z decays"),
         "LLNuNu"     : NTupleCollection("llnunu",LLNuNuType ,5, help="VV candidate with di-lepton and MET"),
         "ElMuNuNu"   : NTupleCollection("elmununu",LLNuNuType ,5, help="Fake VV candidate with el-mu and MET"),
         "jets"       : NTupleCollection("jet_corr",JetType,15, help="all jets with new JEC for 76X applied"),
         #"jets_raw"   : NTupleCollection("jet",jetType,15, help="all jets from miniAOD"),
     }
)

vvTreeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='vvTreeProducer',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
     globalVariables = [
         NTupleVariable("nLL",lambda ev: len(ev.LL) , int),      
         #NTupleVariable("LHEweight_original", lambda ev: ev.LHE_originalWeight if  hasattr(ev,'LHE_originalWeight') else  0, mcOnly=True, help="original LHE weight"), 
         #NTupleVariable("nElMu",lambda ev: len(ev.ElMu) , int),       
         NTupleVariable("nLLNuNu",lambda ev: len(ev.LLNuNu) , int),       
         NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"), 
         #NTupleVariable("nVertAll",  lambda ev: len(ev.vertices), int, help="Number of good vertices"), 
         NTupleVariable("vtx_x",  lambda ev: ev.goodVertices[0].x(), float, help="primary vertex x"), 
         NTupleVariable("vtx_y",  lambda ev: ev.goodVertices[0].y(), float, help="primary vertex y"), 
         NTupleVariable("vtx_z",  lambda ev: ev.goodVertices[0].z(), float, help="primary vertex z"), 
         #NTupleVariable("triggersf",  lambda x : getattr(x,'trgsf',1), float, mcOnly=True, help="singleelectron/muon trigger sf"),
         #NTupleVariable("triggersfUp",  lambda x : getattr(x,'trgsfUp',1),float, mcOnly=True, help="singleelectron/muon trigger sf upper"),
         #NTupleVariable("triggersfLo",  lambda x : getattr(x,'trgsfLo',1), float, mcOnly=True, help="singleelectron/muon trigger sf lower"),
         NTupleVariable("rho", lambda ev: ev.rho, float),
         NTupleVariable("lheNb", lambda ev: ev.lheNb, int, mcOnly=True),
         NTupleVariable("lheNj", lambda ev: ev.lheNj, int, mcOnly=True),
     ],

     globalObjects =  {
         #"met" : NTupleObject("met", metType, help="PF E_{T}^{miss}, after type 1 corrections"),
         #"met_miniAod" : NTupleObject("met_raw", metType, help="PF E_{T}^{miss}, after type 1 corrections in miniAOD"),
         #"met_JEC" : NTupleObject("met_JEC", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC added"),
         #"met_JECUp" : NTupleObject("met_JECUp", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC up"),
         #"met_JECDown" : NTupleObject("met_JECDown", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC down"),
         #"met_JECJER" : NTupleObject("met_JECJER", metType, mcOnly=True, help="PF E_{T}^{miss}, after type 1 corrections with JEC+JER jets"),
     },

     collections = {
         #"LL"  : NTupleCollection("Zll",LLType,5, help="Z to ll"),
         #"ElMu"  : NTupleCollection("elmu",LLType,5, help="electron - muon pair for non-resonant bkg"),
         "selectedLeptons" : NTupleCollection("lep",leptonType,100, help="selected leptons"),
         "genLeptons" : NTupleCollection("genLep", genParticleType, 100, mcOnly=True, help="Generated leptons (e/mu) from W/Z decays"),
         "genZBosons" : NTupleCollection("genZ", tlorentzFourVectorType, 100, mcOnly=True, help="Generated V bosons"),
         "LLNuNu"     : NTupleCollection("llnunu",LLNuNuType ,5, help="VV candidate with di-lepton and MET"),
         "genXZZ" : NTupleCollection("genX", tlorentzFourVectorType, 10, mcOnly=True, help="Generated X->ZZ"),
	 "jets"       : NTupleCollection("jet",JetTypeExtra,300, help="all jets in miniaod"),
         #"LHE_weights"    : NTupleCollection("LHEweight",  weightsInfoType, 1000, mcOnly=True, help="LHE weight info"),
         "genJets"       : NTupleCollection("genjet",genJetType,300, mcOnly=True, help="genJets in miniaod"),
         #"selectedPhotons"       : NTupleCollection("photon",photonType,100, help="selected photons in miniaod"),
         #"cleanJetsAll"     : NTupleCollection("cleanJetsAll",jetType,100, help="cleaned jets"),
         #"vertices"   : NTupleCollection("allvtx", threeVectorType, 300, help="all vertecies"),
         #"goodVertices"   : NTupleCollection("goodvtx", threeVectorType, 300, help="good vertecies"),
     }
)


lepeffTreeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='lepeffTreeProducer',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
     globalVariables = [
        NTupleVariable("nMuon",lambda ev: len(ev.selectedMuons) if hasattr(ev,"selectedMuons") else 0, int),
        NTupleVariable("nElectron",lambda ev: len(ev.selectedElectrons) if hasattr(ev,"selectedElectrons") else 0, int),
        ],
     collections = {
         "llpair"  : NTupleCollection("llpair",llpairType ,15, help="lepton eff study"),
     }
)

leptonTreeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='leptonTreeProducer',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
#     PDFWeights = PDFWeights,
#     globalVariables = susyMultilepton_globalVariables,
#     globalObjects = susyMultilepton_globalObjects,
     collections = {
            "genleps"          : NTupleCollection("gen",     genParticleWithLinksType, 10, help="Generated leptons (e/mu) from W/Z decays"),                                                                                                
            "inclusiveLeptons" : NTupleCollection("l",    leptonTypeExtra, 10, help="Inclusive Leptons"),                                                                                                
     }
)


jetTreeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='jetTreeProducer',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
     globalVariables = [
         NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"), 
     ],
     globalObjects =  {
         "met" : NTupleObject("met", metType, help="PF E_{T}^{miss}, after full type 1 corrections"),
         "met_miniAod" : NTupleObject("met_raw", metType, help="PF E_{T}^{miss}, after type 1 corrections in miniAOD"),
         "met_JEC" : NTupleObject("met_JEC", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC added"),
         "met_JECUp" : NTupleObject("met_JECUp", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC up"),
         "met_JECDown" : NTupleObject("met_JECDown", metType, help="PF E_{T}^{miss}, after type 1 corrections with new 76X JEC down"),
         "met_JECJER" : NTupleObject("met_JECJER", metType, help="PF E_{T}^{miss}, after type 1 corrections with JEC+JER jets"),
     },

     collections = {
         "selectedLeptons" : NTupleCollection("lep",leptonType,10, help="selected leptons"),
         "genLeptons" : NTupleCollection("genLep", genParticleType, 10, help="Generated leptons (e/mu) from W/Z decays"),
         "genZBosons" : NTupleCollection("genZ", genParticleType, 10, help="Generated V bosons"),
         "genLeptons" : NTupleCollection("genLep", genParticleType, 10, help="Generated leptons (e/mu) from W/Z decays"),
         "genZBosons" : NTupleCollection("genZ", genParticleType, 10, help="Generated V bosons"),
         "jets"  : NTupleCollection("jet_corr",JetType,15, help="all jets with new JEC for 76X applied"),
         "jets_raw"  : NTupleCollection("jet",JetType,15, help="all jets from miniAOD"),
     }
)
