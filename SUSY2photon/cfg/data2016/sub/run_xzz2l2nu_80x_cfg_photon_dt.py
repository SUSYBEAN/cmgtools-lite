##########################################################
##      configuration for XZZ2l2nu 
##########################################################

import CMGTools.XZZ2l2nu.fwlite.Config as cfg
from CMGTools.XZZ2l2nu.fwlite.Config import printComps
from CMGTools.XZZ2l2nu.RootTools import *
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption


#Load all common analyzers
from CMGTools.XZZ2l2nu.analyzers.coreXZZ_cff import *

#-------- SAMPLES AND TRIGGERS -----------
from CMGTools.XZZ2l2nu.samples.loadSamples80x import *
selectedComponents = mcSamples+dataSamples

triggerFlagsAna.triggerBits ={
    "ISOMU":triggers_1mu_iso,
    "MU":triggers_1mu_noniso,
    "MUv2":triggers_1mu_noniso_v2,
    "MU50":triggers_1mu_noniso_M50,
    "ISOELE":triggers_1e,
    "ELE":triggers_1e_noniso,
    "ELEv2":triggers_1e_noniso_v2,
    "ELE115":triggers_1e_noniso_E115,
    "MUMU": triggers_mumu,
    "MUMUNOISO":triggers_mumu_noniso,
    "ELEL": triggers_ee,
    "HT800":triggers_HT800,
    "HT900":triggers_HT900,
    "JJ":triggers_dijet_fat,
    "MET90":triggers_met90_mht90+triggers_metNoMu90_mhtNoMu90,
    "MET120":triggers_metNoMu120_mhtNoMu120,
    "PHOTONHZZ": triggers_photon_unbias,
    "ALLPHOTON": triggers_all_photons
}

#-------- Analyzer
from CMGTools.XZZ2l2nu.analyzers.treeXZZ_cff import *

multiStateAna.processTypes = ['PhotonJets']
multiStateAna.selectPhotonJets = (lambda x: x.leg1.pt()>20.0 and x.leg2.pt()>-0.0)
vvSkimmer.required = ['PhotonJets']

vvTreeProducer.globalVariables = [
         NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"), 
         NTupleVariable("nVertAll",  lambda ev: len(ev.vertices), int, help="Number of good vertices"), 
         NTupleVariable("rho", lambda ev: ev.rho , float),
     ]
vvTreeProducer.globalObjects =  {  }

vvTreeProducer.collections = {
	 "jets"       : NTupleCollection("jet",JetType,100, help="all jets in miniaod"),
         "selectedPhotons"       : NTupleCollection("photon",photonType,100, help="selected photons in miniaod"),
         "PhotonJets"     : NTupleCollection("gjet",PhotonJetType ,100, help="photon and MET"),
     }




#-------- SEQUENCE
coreSequence = [
    skimAnalyzer,
    genAna,
    jsonAna,
    triggerAna,
    pileUpAna,
    vertexAna,
    #lepAna,
    photonAna, 
    jetAna,
    metAna,
    multiStateAna,
    eventFlagsAna,
    triggerFlagsAna,
]
    
#sequence = cfg.Sequence(coreSequence)
sequence = cfg.Sequence(coreSequence+[vvSkimmer,vvTreeProducer])
#sequence = cfg.Sequence(coreSequence+[vvSkimmer,multtrg,vvTreeProducer])
#sequence = cfg.Sequence(coreSequence+[vvSkimmer,fullTreeProducer])
 

#-------- HOW TO RUN
test = 1
if test==1:
    # test a single component, using a single thread.
    #selectedComponents = dataSamples
    #selectedComponents = mcSamples
    selectedComponents = SinglePhoton
    #selectedComponents = [SinglePhoton_Run2016D_PromptReco_v2]
    #selectedComponents = [GJet_Pt_20toInf_DoubleEMEnriched]
    #selectedComponents = [GJet_Pt_20to40_DoubleEMEnriched, GJet_Pt_40toInf_DoubleEMEnriched]
    #selectedComponents = [SingleMuon_Run2015D_Promptv4,SingleElectron_Run2015D_Promptv4]
    #selectedComponents = [SingleMuon_Run2015C_25ns_16Dec]
    #selectedComponents = [SingleMuon_Run2016B_PromptReco_v2] 
    #selectedComponents = SingleMuon+SingleElectron
    #selectedComponents = [SingleMuon_Run2016B_PromptReco_v2,SingleElectron_Run2016B_PromptReco_v2] 
    #selectedComponents = [SingleMuon_Run2016D_PromptReco_v2,SingleElectron_Run2016D_PromptReco_v2] 
    #selectedComponents = [MuonEG_Run2015D_16Dec] #MuEG
    #selectedComponents = [RSGravToZZToZZinv_narrow_800]
    #selectedComponents = [DYJetsToLL_M50]
    #selectedComponents = [DYJetsToLL_M50_MGMLM_Ext1]
    #selectedComponents = [BulkGravToZZToZlepZinv_narrow_600] 
    #selectedComponents = signalSamples
    #selectedComponents = [TTTo2L2Nu]
    #selectedComponents = [BulkGravToZZ_narrow_800]
    #selectedComponents = [BulkGravToZZToZlepZhad_narrow_800]
    for c in selectedComponents:
        #c.files = c.files[3:10]
        c.splitFactor = (len(c.files)/5 if len(c.files)>5 else 1)
        #c.splitFactor = 7
        #c.triggers=triggers_1mu_noniso
        #c.triggers=triggers_1e_noniso

## output histogram
outputService=[]
from PhysicsTools.HeppyCore.framework.services.tfile import TFileService
output_service = cfg.Service(
    TFileService,
    'outputfile',
    name="outputfile",
    fname='vvTreeProducer/tree.root',
    option='recreate'
    )
outputService.append(output_service)

from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
from CMGTools.TTHAnalysis.tools.EOSEventsWithDownload import EOSEventsWithDownload
event_class = EOSEventsWithDownload
event_class = Events
if getHeppyOption("nofetch"):
    event_class = Events
config = cfg.Config( components = selectedComponents,
                     sequence = sequence,
                     services = [],
                     events_class = event_class)




