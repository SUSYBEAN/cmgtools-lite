##########################################################
##      configuration for XZZ2l2nu 
##########################################################

import CMGTools.XZZ2l2nu.fwlite.Config as cfg
from CMGTools.XZZ2l2nu.fwlite.Config import printComps
from CMGTools.XZZ2l2nu.RootTools import *
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption


#Load all common analyzers
from CMGTools.XZZ2l2nu.analyzers.coreXZZ_cff import *
from CMGTools.XZZ2l2nu.analyzers.XZZTrgEff import *
#-------- SAMPLES AND TRIGGERS -----------
from CMGTools.XZZ2l2nu.samples.loadSamples80x import *

#-------- Analyzer
from CMGTools.XZZ2l2nu.analyzers.treeXZZ_cff import *

coreSequence = [
    skimAnalyzer,
    genAna,
]

vvTreeProducer.globalVariables = [
         NTupleVariable("lheNb", lambda ev: ev.lheNb, int),
         NTupleVariable("lheNj", lambda ev: ev.lheNj, int),
]
vvTreeProducer.globalObjects =  {  }
vvTreeProducer.collections = {
         "genLeptons" : NTupleCollection("genLep", genParticleType, 100, help="Generated leptons (e/mu) from W/Z decays"),
         "genLeptonsFsr" : NTupleCollection("genLepFsr", genParticleType, 100, help="Generated leptons (e/mu) from W/Z decays"),
         "mother_genLeptons" : NTupleCollection("genLepMother", genParticleType, 100, help="Generated leptons (e/mu) from W/Z decays"),
         "genMuons" : NTupleCollection("genMuon", genParticleWithMotherInfo, 100, help="Generated leptons (e/mu) from W/Z decays"),
         "genElectrons" : NTupleCollection("genElec", genParticleWithMotherInfo, 100, help="Generated leptons (e/mu) from W/Z decays"),
         "genTaus" : NTupleCollection("genTau", genParticleType, 100, help="Generated leptons (e/mu) from W/Z decays"),
         "genMuonsFsr" : NTupleCollection("genMuonFsr", genParticleWithMotherInfo, 100, help="Generated leptons (e/mu) from W/Z decays"),
         "genElectronsFsr" : NTupleCollection("genElecFsr", genParticleWithMotherInfo, 100, help="Generated leptons (e/mu) from W/Z decays"),
         "genTausFsr" : NTupleCollection("genTauFsr", genParticleWithMotherInfo, 100, help="Generated leptons (e/mu) from W/Z decays"),
         "genNeutrinos" : NTupleCollection("genNeu", genParticleType, 300, help="Generated leptons (e/mu) from W/Z decays"),
         "genZBosons" : NTupleCollection("genZ", tlorentzFourVectorType, 100, help="Generated V bosons"),
     }
    
sequence = cfg.Sequence(coreSequence+[vvTreeProducer])
 

#-------- HOW TO RUN
test = 1
if test==1:
    #selectedComponents = [BulkGravToZZToZlepZinv_narrow_1000]
    #selectedComponents = [DY1JetsToLL_M50_MGMLM]
    #selectedComponents = [DY1JetsToLL_M50_MGMLM, DY2JetsToLL_M50_MGMLM, DY3JetsToLL_M50_MGMLM, DY4JetsToLL_M50_MGMLM, DYBJetsToLL_M50_MGMLM]
    #selectedComponents = [DYJetsToLL_M50, DYJetsToLL_M50_MGMLM_Ext1, DY1JetsToLL_M50_MGMLM, DY2JetsToLL_M50_MGMLM, DY3JetsToLL_M50_MGMLM, DY4JetsToLL_M50_MGMLM, DYBJetsToLL_M50_MGMLM]
    #selectedComponents = [DYJetsToLL_M50_Ext]
    selectedComponents = [DYJetsToLL_M50]
    #selectedComponents = [DYJetsToLL_M50_MGMLM_Ext1]
    for c in selectedComponents:
        c.files = c.files[1]
        #c.splitFactor = (len(c.files)/30 if len(c.files)>30 else 1)
        c.splitFactor = 1


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




