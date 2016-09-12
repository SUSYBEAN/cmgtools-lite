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
from CMGTools.XZZ2l2nu.samples.loadSamples76x import *
selectedComponents = mcSamples+dataSamples


#-------- Analyzer
from CMGTools.XZZ2l2nu.analyzers.treeXZZ_cff import *


vvTreeProducer.globalVariables = [  NTupleVariable("LHEweight_original", lambda ev: ev.LHE_originalWeight if  hasattr(ev,'LHE_originalWeight') else  0, mcOnly=True, help="original LHE weight"), ]
vvTreeProducer.globalObjects =  { }
vvTreeProducer.collections = { "LHE_weights"    : NTupleCollection("LHEweight",  weightsInfoType, 1000, help="LHE weight info"),}



#-------- SEQUENCE
coreSequence = [
    skimAnalyzer,
    lheWeightAna,
]
    
sequence = cfg.Sequence(coreSequence+[vvTreeProducer])

 

#-------- HOW TO RUN
test = 1
if test==1:
    # test a single component, using a single thread.
    #selectedComponents = signalSamples
    selectedComponents = [BulkGravToZZToZlepZinv_narrow_800]
    for c in selectedComponents:
        c.files = c.files[0]
        #c.splitFactor = (len(c.files)/10 if len(c.files)>10 else 1)
        c.splitFactor = 1
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




