##########################################################
##      configuration for XZZ2l2nu 
##########################################################

import CMGTools.XZZ2l2nu.fwlite.Config as cfg
from CMGTools.XZZ2l2nu.fwlite.Config import printComps
from CMGTools.XZZ2l2nu.RootTools import *
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption


#Load all common analyzers
from CMGTools.XZZ2l2nu.analyzers.coreXZZ_cff import *
from CMGTools.XZZ2l2nu.analyzers.XZZMuonEffTree import *

#-------- SAMPLES AND TRIGGERS -----------
from CMGTools.XZZ2l2nu.samples.loadSamples76x import *

#-------- Analyzer
from CMGTools.XZZ2l2nu.analyzers.treeXZZ_cff import *

#-------- SEQUENCE
#sequence = cfg.Sequence(coreSequence+[vvSkimmer,vvTreeProducer])
meffAna = cfg.Analyzer(
    XZZMuonEffTree,
    name='mEffTree',
    genfilter=True,
    pfbkg=False,
    eithercharge=False,
    checktag=True,
    muHLT="HLT_IsoMu20_v"
    )

lepAna.applyIso=False
lepAna.applyID=False
vertexAna.keepFailingEvents=True
leptonType.variables.extend([
    # Extra muon ID working points
    NTupleVariable("softMuonId", lambda x : x.muonID("POG_ID_Soft") if abs(x.pdgId())==13 else -100, int, help="Muon POG Soft id"),
    NTupleVariable("pfMuonId",   lambda x : x.muonID("POG_ID_Loose") if abs(x.pdgId())==13 else -100, int, help="Muon POG Loose id"),
    NTupleVariable("MediumMuonId",   lambda x : x.muonID("POG_ID_Medium") if abs(x.pdgId())==13 else -100, int, help="Muon POG Medium id"),
    NTupleVariable("tightMuonId",   lambda x : x.muonID("POG_ID_Tight") if abs(x.pdgId())==13 else -100, int, help="Muon POG Tight id"),
    NTupleVariable("isTag",   lambda x : x.istag if abs(x.pdgId())==13 else -100, int, help="Muon POG Tight id"),
]) 
sequence = [
    jsonAna,
    pileUpAna,
    vertexAna,
    lepAna,
    meffAna,
    lepeffTreeProducer
]
    

#-------- HOW TO RUN
test = 1
if test==1:
    # test a single component, using a single thread.
    selectedComponents = [DYJetsToLL_M50]
    #selectedComponents = mcSamples
    #selectedComponents = [SingleMuon_Run2015D_Promptv4,SingleElectron_Run2015D_Promptv4]
    #[SingleElectron_Run2015D_Promptv4,SingleElectron_Run2015D_05Oct]
    #selectedComponents = [RSGravToZZToZZinv_narrow_800]
    #selectedComponents = [DYJetsToLL_M50]
    #selectedComponents = [BulkGravToZZ_narrow_800]
    #selectedComponents = [BulkGravToZZToZlepZhad_narrow_800]
    for c in selectedComponents:
        #c.files = c.files[:1]
        c.splitFactor = (len(c.files)/10 if len(c.files)>10 else 1)
        #c.splitFactor = 1
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




