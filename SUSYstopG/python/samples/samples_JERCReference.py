import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()

JERCRef_MC=kreator.makeMCComponent("JERCRef_MC","/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/MINIAODSIM","CMS",".*root",1)

#/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/02837459-03C2-E511-8EA2-002590A887AC.root
JERCRef_MC_eos=kreator.makeMCComponentFromEOS('JERCRef_MC_eos', '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/', '/store/mc/RunIIFall15MiniAODv2/%s', '.*root', 1)

#json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
# update json HN: https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2577.html
json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'

# July17 for run<=251561 with MINIAOD reprocessed with correct MET filters
run_range = (251244, 251562)
label = "_runs%s_%s"%(run_range[0], run_range[1])

JERCRef_data=kreator.makeDataComponent("DoubleMuon_Run2015D_16Dec"    , "/DoubleMuon/Run2015D-16Dec2015-v1/MINIAOD"    , "CMS", ".*root", json)
