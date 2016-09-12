import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()

#json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
# update json HN: https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2577.html
json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'

# July17 for run<=251561 with MINIAOD reprocessed with correct MET filters
run_range = (251244, 251562)
label = "_runs%s_%s"%(run_range[0], run_range[1])

JetHT_Run2015B_17Jul2015          = kreator.makeDataComponent("JetHT_Run2015B_17Jul2015"         , "/JetHT/Run2015B-17Jul2015-v1/MINIAOD"         , "CMS", ".*root", json, run_range)
HTMHT_Run2015B_17Jul2015          = kreator.makeDataComponent("HTMHT_Run2015B_17Jul2015"         , "/HTMHT/Run2015B-17Jul2015-v1/MINIAOD"         , "CMS", ".*root", json, run_range)
MET_Run2015B_17Jul2015            = kreator.makeDataComponent("MET_Run2015B_17Jul2015"           , "/MET/Run2015B-17Jul2015-v1/MINIAOD"           , "CMS", ".*root", json, run_range)
SingleElectron_Run2015B_17Jul2015 = kreator.makeDataComponent("SingleElectron_Run2015B_17Jul2015", "/SingleElectron/Run2015B-17Jul2015-v1/MINIAOD", "CMS", ".*root", json, run_range)
SingleMuon_Run2015B_17Jul2015     = kreator.makeDataComponent("SingleMuon_Run2015B_17Jul2015"    , "/SingleMuon/Run2015B-17Jul2015-v1/MINIAOD"    , "CMS", ".*root", json, run_range)
SinglePhoton_Run2015B_17Jul2015   = kreator.makeDataComponent("SinglePhoton_Run2015B_17Jul2015"  , "/SinglePhoton/Run2015B-17Jul2015-v1/MINIAOD"  , "CMS", ".*root", json, run_range)
DoubleEG_Run2015B_17Jul2015       = kreator.makeDataComponent("DoubleEG_Run2015B_17Jul2015"      , "/DoubleEG/Run2015B-17Jul2015-v1/MINIAOD"      , "CMS", ".*root", json, run_range)
MuonEG_Run2015B_17Jul2015         = kreator.makeDataComponent("MuonEG_Run2015B_17Jul2015"        , "/MuonEG/Run2015B-17Jul2015-v1/MINIAOD"        , "CMS", ".*root", json, run_range)
DoubleMuon_Run2015B_17Jul2015     = kreator.makeDataComponent("DoubleMuon_Run2015B_17Jul2015"    , "/DoubleMuon/Run2015B-17Jul2015-v1/MINIAOD"    , "CMS", ".*root", json, run_range)


# PromptReco-v1 for run > 251561
run_range = (251643, 251883)
label = "_runs%s_%s"%(run_range[0], run_range[1])

JetHT_Run2015B_PromptReco          = kreator.makeDataComponent("JetHT_Run2015B_PromptReco"         , "/JetHT/Run2015B-PromptReco-v1/MINIAOD"         , "CMS", ".*root", json, run_range)
HTMHT_Run2015B_PromptReco          = kreator.makeDataComponent("HTMHT_Run2015B_PromptReco"         , "/HTMHT/Run2015B-PromptReco-v1/MINIAOD"         , "CMS", ".*root", json, run_range)
MET_Run2015B_PromptReco            = kreator.makeDataComponent("MET_Run2015B_PromptReco"           , "/MET/Run2015B-PromptReco-v1/MINIAOD"           , "CMS", ".*root", json, run_range)
SingleElectron_Run2015B_PromptReco = kreator.makeDataComponent("SingleElectron_Run2015B_PromptReco", "/SingleElectron/Run2015B-PromptReco-v1/MINIAOD", "CMS", ".*root", json, run_range)
SingleMuon_Run2015B_PromptReco     = kreator.makeDataComponent("SingleMuon_Run2015B_PromptReco"    , "/SingleMuon/Run2015B-PromptReco-v1/MINIAOD"    , "CMS", ".*root", json, run_range)
SinglePhoton_Run2015B_PromptReco   = kreator.makeDataComponent("SinglePhoton_Run2015B_PromptReco"  , "/SinglePhoton/Run2015B-PromptReco-v1/MINIAOD"  , "CMS", ".*root", json, run_range)
DoubleEG_Run2015B_PromptReco       = kreator.makeDataComponent("DoubleEG_Run2015B_PromptReco"      , "/DoubleEG/Run2015B-PromptReco-v1/MINIAOD"      , "CMS", ".*root", json, run_range)
MuonEG_Run2015B_PromptReco         = kreator.makeDataComponent("MuonEG_Run2015B_PromptReco"        , "/MuonEG/Run2015B-PromptReco-v1/MINIAOD"        , "CMS", ".*root", json, run_range)
DoubleMuon_Run2015B_PromptReco     = kreator.makeDataComponent("DoubleMuon_Run2015B_PromptReco"    , "/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD"    , "CMS", ".*root", json, run_range)



### ----------------------------- 17July re-reco ----------------------------------------
# https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDataReprocessing747reMiniAod2015B

Jet_Run2015B_17Jul            = kreator.makeDataComponent("Jet_Run2015B_17Jul"           , "/Jet/Run2015B-17Jul2015-v1/MINIAOD"           , "CMS", ".*root", json)
JetHT_Run2015B_17Jul          = kreator.makeDataComponent("JetHT_Run2015B_17Jul"         , "/JetHT/Run2015B-17Jul2015-v1/MINIAOD"         , "CMS", ".*root", json)
HTMHT_Run2015B_17Jul          = kreator.makeDataComponent("HTMHT_Run2015B_17Jul"         , "/HTMHT/Run2015B-17Jul2015-v1/MINIAOD"         , "CMS", ".*root", json)
MET_Run2015B_17Jul            = kreator.makeDataComponent("MET_Run2015B_17Jul"           , "/MET/Run2015B-17Jul2015-v1/MINIAOD"           , "CMS", ".*root", json)
SingleElectron_Run2015B_17Jul = kreator.makeDataComponent("SingleElectron_Run2015B_17Jul", "/SingleElectron/Run2015B-17Jul2015-v1/MINIAOD", "CMS", ".*root", json)
SingleMu_Run2015B_17Jul       = kreator.makeDataComponent("SingleMu_Run2015B_17Jul"      , "/SingleMu/Run2015B-17Jul2015-v1/MINIAOD"      , "CMS", ".*root", json)
SingleMuon_Run2015B_17Jul     = kreator.makeDataComponent("SingleMuon_Run2015B_17Jul"    , "/SingleMuon/Run2015B-17Jul2015-v1/MINIAOD"    , "CMS", ".*root", json)
SinglePhoton_Run2015B_17Jul   = kreator.makeDataComponent("SinglePhoton_Run2015B_17Jul"  , "/SinglePhoton/Run2015B-17Jul2015-v1/MINIAOD"  , "CMS", ".*root", json)
EGamma_Run2015B_17Jul         = kreator.makeDataComponent("EGamma_Run2015B_17Jul"        , "/EGamma/Run2015B-17Jul2015-v1/MINIAOD"        , "CMS", ".*root", json)
DoubleEG_Run2015B_17Jul       = kreator.makeDataComponent("DoubleEG_Run2015B_17Jul"      , "/DoubleEG/Run2015B-17Jul2015-v1/MINIAOD"      , "CMS", ".*root", json)
MuonEG_Run2015B_17Jul         = kreator.makeDataComponent("MuonEG_Run2015B_17Jul"        , "/MuonEG/Run2015B-17Jul2015-v1/MINIAOD"        , "CMS", ".*root", json)
DoubleMuon_Run2015B_17Jul     = kreator.makeDataComponent("DoubleMuon_Run2015B_17Jul"    , "/DoubleMuon/Run2015B-17Jul2015-v1/MINIAOD"    , "CMS", ".*root", json)

minBias_Run2015B_17Jul  = kreator.makeDataComponent("minBias_Run2015B_17Jul" , "/MinimumBias/Run2015B-17Jul2015-v1/MINIAOD", "CMS", ".*root", json)
zeroBias_Run2015B_17Jul = kreator.makeDataComponent("zeroBias_Run2015B_17Jul", "/ZeroBias/Run2015B-17Jul2015-v1/MINIAOD"   , "CMS", ".*root", json)



### ----------------------------- Run2015C ----------------------------------------

Jet_Run2015C            = kreator.makeDataComponent("Jet_Run2015C"           , "/Jet/Run2015C-PromptReco-v1/MINIAOD"           , "CMS", ".*root", json)
JetHT_Run2015C          = kreator.makeDataComponent("JetHT_Run2015C"         , "/JetHT/Run2015C-PromptReco-v1/MINIAOD"         , "CMS", ".*root", json)
HTMHT_Run2015C          = kreator.makeDataComponent("HTMHT_Run2015C"         , "/HTMHT/Run2015C-PromptReco-v1/MINIAOD"         , "CMS", ".*root", json)
MET_Run2015C            = kreator.makeDataComponent("MET_Run2015C"           , "/MET/Run2015C-PromptReco-v1/MINIAOD"           , "CMS", ".*root", json)
SingleElectron_Run2015C = kreator.makeDataComponent("SingleElectron_Run2015C", "/SingleElectron/Run2015C-PromptReco-v1/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2015C     = kreator.makeDataComponent("SingleMuon_Run2015C"    , "/SingleMuon/Run2015C-PromptReco-v1/MINIAOD"    , "CMS", ".*root", json)
SinglePhoton_Run2015C   = kreator.makeDataComponent("SinglePhoton_Run2015C"  , "/SinglePhoton/Run2015C-PromptReco-v1/MINIAOD"  , "CMS", ".*root", json)
EGamma_Run2015C         = kreator.makeDataComponent("EGamma_Run2015C"        , "/EGamma/Run2015C-PromptReco-v1/MINIAOD"        , "CMS", ".*root", json)
DoubleEG_Run2015C       = kreator.makeDataComponent("DoubleEG_Run2015C"      , "/DoubleEG/Run2015C-PromptReco-v1/MINIAOD"      , "CMS", ".*root", json)
MuonEG_Run2015C         = kreator.makeDataComponent("MuonEG_Run2015C"        , "/MuonEG/Run2015C-PromptReco-v1/MINIAOD"        , "CMS", ".*root", json)
DoubleMuon_Run2015C     = kreator.makeDataComponent("DoubleMuon_Run2015C"    , "/DoubleMuon/Run2015C-PromptReco-v1/MINIAOD"    , "CMS", ".*root", json)

minBias_Run2015C  = kreator.makeDataComponent("minBias_Run2015C" , "/MinimumBias/Run2015C-PromptReco-v1/MINIAOD", "CMS", ".*root", json)
zeroBias_Run2015C = kreator.makeDataComponent("zeroBias_Run2015C", "/ZeroBias/Run2015C-PromptReco-v1/MINIAOD"   , "CMS", ".*root", json)


### ----------------------------- Run2015D miniAODv1 ----------------------------------------

#Jet_Run2015D            = kreator.makeDataComponent("Jet_Run2015D"           , "/Jet/Run2015D-PromptReco-v3/MINIAOD"           , "CMS", ".*root", json)
JetHT_Run2015D          = kreator.makeDataComponent("JetHT_Run2015D"         , "/JetHT/Run2015D-PromptReco-v3/MINIAOD"         , "CMS", ".*root", json)
HTMHT_Run2015D          = kreator.makeDataComponent("HTMHT_Run2015D"         , "/HTMHT/Run2015D-PromptReco-v3/MINIAOD"         , "CMS", ".*root", json)
MET_Run2015D            = kreator.makeDataComponent("MET_Run2015D"           , "/MET/Run2015D-PromptReco-v3/MINIAOD"           , "CMS", ".*root", json)
SingleElectron_Run2015D = kreator.makeDataComponent("SingleElectron_Run2015D", "/SingleElectron/Run2015D-PromptReco-v3/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2015D     = kreator.makeDataComponent("SingleMuon_Run2015D"    , "/SingleMuon/Run2015D-PromptReco-v3/MINIAOD"    , "CMS", ".*root", json)
SinglePhoton_Run2015D   = kreator.makeDataComponent("SinglePhoton_Run2015D"  , "/SinglePhoton/Run2015D-PromptReco-v3/MINIAOD"  , "CMS", ".*root", json)
#EGamma_Run2015D         = kreator.makeDataComponent("EGamma_Run2015D"        , "/EGamma/Run2015D-PromptReco-v3/MINIAOD"        , "CMS", ".*root", json)
DoubleEG_Run2015D       = kreator.makeDataComponent("DoubleEG_Run2015D"      , "/DoubleEG/Run2015D-PromptReco-v3/MINIAOD"      , "CMS", ".*root", json)
MuonEG_Run2015D         = kreator.makeDataComponent("MuonEG_Run2015D"        , "/MuonEG/Run2015D-PromptReco-v3/MINIAOD"        , "CMS", ".*root", json)
DoubleMuon_Run2015D     = kreator.makeDataComponent("DoubleMuon_Run2015D"    , "/DoubleMuon/Run2015D-PromptReco-v3/MINIAOD"    , "CMS", ".*root", json)


### ----------------------------- Run2015D miniAODv2 ----------------------------------------
###https://hypernews.cern.ch/HyperNews/CMS/get/physics-announcements/3561.html

JetHT_Run2015D_Promptv4          = kreator.makeDataComponent("JetHT_Run2015D_v4"         , "/JetHT/Run2015D-PromptReco-v4/MINIAOD"         , "CMS", ".*root", json)
HTMHT_Run2015D_Promptv4          = kreator.makeDataComponent("HTMHT_Run2015D_v4"         , "/HTMHT/Run2015D-PromptReco-v4/MINIAOD"         , "CMS", ".*root", json)
MET_Run2015D_Promptv4            = kreator.makeDataComponent("MET_Run2015D_v4"           , "/MET/Run2015D-PromptReco-v4/MINIAOD"           , "CMS", ".*root", json)
SingleElectron_Run2015D_Promptv4 = kreator.makeDataComponent("SingleElectron_Run2015D_v4", "/SingleElectron/Run2015D-PromptReco-v4/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2015D_Promptv4     = kreator.makeDataComponent("SingleMuon_Run2015D_v4"    , "/SingleMuon/Run2015D-PromptReco-v4/MINIAOD"    , "CMS", ".*root", json)
SinglePhoton_Run2015D_Promptv4   = kreator.makeDataComponent("SinglePhoton_Run2015D_v4"  , "/SinglePhoton/Run2015D-PromptReco-v4/MINIAOD"  , "CMS", ".*root", json)
DoubleEG_Run2015D_Promptv4       = kreator.makeDataComponent("DoubleEG_Run2015D_v4"      , "/DoubleEG/Run2015D-PromptReco-v4/MINIAOD"      , "CMS", ".*root", json)
MuonEG_Run2015D_Promptv4         = kreator.makeDataComponent("MuonEG_Run2015D_v4"        , "/MuonEG/Run2015D-PromptReco-v4/MINIAOD"        , "CMS", ".*root", json)
DoubleMuon_Run2015D_Promptv4     = kreator.makeDataComponent("DoubleMuon_Run2015D_v4"    , "/DoubleMuon/Run2015D-PromptReco-v4/MINIAOD"    , "CMS", ".*root", json)



### ----------------------------- Run2015D-05Oct2015 ----------------------------------------
## https://hypernews.cern.ch/HyperNews/CMS/get/datasets/4154.html

JetHT_Run2015D_05Oct          = kreator.makeDataComponent("JetHT_Run2015D_05Oct"         , "/JetHT/Run2015D-05Oct2015-v1/MINIAOD"         , "CMS", ".*root", json)
HTMHT_Run2015D_05Oct          = kreator.makeDataComponent("HTMHT_Run2015D_05Oct"         , "/HTMHT/Run2015D-05Oct2015-v1/MINIAOD"         , "CMS", ".*root", json)
MET_Run2015D_05Oct            = kreator.makeDataComponent("MET_Run2015D_05Oct"           , "/MET/Run2015D-05Oct2015-v1/MINIAOD"           , "CMS", ".*root", json)
SingleElectron_Run2015D_05Oct = kreator.makeDataComponent("SingleElectron_Run2015D_05Oct", "/SingleElectron/Run2015D-05Oct2015-v1/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2015D_05Oct     = kreator.makeDataComponent("SingleMuon_Run2015D_05Oct"    , "/SingleMuon/Run2015D-05Oct2015-v1/MINIAOD"    , "CMS", ".*root", json)
SinglePhoton_Run2015D_05Oct   = kreator.makeDataComponent("SinglePhoton_Run2015D_05Oct"  , "/SinglePhoton/Run2015D-05Oct2015-v1/MINIAOD"  , "CMS", ".*root", json)
DoubleEG_Run2015D_05Oct       = kreator.makeDataComponent("DoubleEG_Run2015D_05Oct"      , "/DoubleEG/Run2015D-05Oct2015-v1/MINIAOD"      , "CMS", ".*root", json)
MuonEG_Run2015D_05Oct         = kreator.makeDataComponent("MuonEG_Run2015D_05Oct"        , "/MuonEG/Run2015D-05Oct2015-v2/MINIAOD"        , "CMS", ".*root", json)
DoubleMuon_Run2015D_05Oct     = kreator.makeDataComponent("DoubleMuon_Run2015D_05Oct"    , "/DoubleMuon/Run2015D-05Oct2015-v1/MINIAOD"    , "CMS", ".*root", json)



### ----------------------------- Run2015B-05Oct2015 ----------------------------------------
## https://hypernews.cern.ch/HyperNews/CMS/get/datasets/4153.html

JetHT_Run2015B_05Oct          = kreator.makeDataComponent("JetHT_Run2015B_05Oct"         , "/JetHT/Run2015B-05Oct2015-v1/MINIAOD"         , "CMS", ".*root", json)
HTMHT_Run2015B_05Oct          = kreator.makeDataComponent("HTMHT_Run2015B_05Oct"         , "/HTMHT/Run2015B-05Oct2015-v1/MINIAOD"         , "CMS", ".*root", json)
MET_Run2015B_05Oct            = kreator.makeDataComponent("MET_Run2015B_05Oct"           , "/MET/Run2015B-05Oct2015-v1/MINIAOD"           , "CMS", ".*root", json)
SingleElectron_Run2015B_05Oct = kreator.makeDataComponent("SingleElectron_Run2015B_05Oct", "/SingleElectron/Run2015B-05Oct2015-v1/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2015B_05Oct     = kreator.makeDataComponent("SingleMuon_Run2015B_05Oct"    , "/SingleMuon/Run2015B-05Oct2015-v1/MINIAOD"    , "CMS", ".*root", json)
SinglePhoton_Run2015B_05Oct   = kreator.makeDataComponent("SinglePhoton_Run2015B_05Oct"  , "/SinglePhoton/Run2015B-05Oct2015-v1/MINIAOD"  , "CMS", ".*root", json)
DoubleEG_Run2015B_05Oct       = kreator.makeDataComponent("DoubleEG_Run2015B_05Oct"      , "/DoubleEG/Run2015B-05Oct2015-v1/MINIAOD"      , "CMS", ".*root", json)
MuonEG_Run2015B_05Oct         = kreator.makeDataComponent("MuonEG_Run2015B_05Oct"        , "/MuonEG/Run2015B-05Oct2015-v1/MINIAOD"        , "CMS", ".*root", json)
DoubleMuon_Run2015B_05Oct     = kreator.makeDataComponent("DoubleMuon_Run2015B_05Oct"    , "/DoubleMuon/Run2015B-05Oct2015-v1/MINIAOD"    , "CMS", ".*root", json)



### ----------------------------- Run2015C_25ns-16Dec2015  ----------------------------------------
## https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset+dataset%3D%2F*%2F*2015C_25ns*16Dec2015*%2F*+status%3D*

JetHT_Run2015C_25ns_16Dec          = kreator.makeDataComponent("JetHT_Run2015C_25ns_16Dec"         , "/JetHT/Run2015C_25ns-16Dec2015-v1/MINIAOD"         , "CMS", ".*root", json)
HTMHT_Run2015C_25ns_16Dec          = kreator.makeDataComponent("HTMHT_Run2015C_25ns_16Dec"         , "/HTMHT/Run2015C_25ns-16Dec2015-v1/MINIAOD"         , "CMS", ".*root", json)
MET_Run2015C_25ns_16Dec            = kreator.makeDataComponent("MET_Run2015C_25ns_16Dec"           , "/MET/Run2015C_25ns-16Dec2015-v1/MINIAOD"           , "CMS", ".*root", json)
SingleElectron_Run2015C_25ns_16Dec = kreator.makeDataComponent("SingleElectron_Run2015C_25ns_16Dec", "/SingleElectron/Run2015C_25ns-16Dec2015-v1/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2015C_25ns_16Dec     = kreator.makeDataComponent("SingleMuon_Run2015C_25ns_16Dec"    , "/SingleMuon/Run2015C_25ns-16Dec2015-v1/MINIAOD"    , "CMS", ".*root", json)
SinglePhoton_Run2015C_25ns_16Dec   = kreator.makeDataComponent("SinglePhoton_Run2015C_25ns_16Dec"  , "/SinglePhoton/Run2015C_25ns-16Dec2015-v1/MINIAOD"  , "CMS", ".*root", json)
DoubleEG_Run2015C_25ns_16Dec       = kreator.makeDataComponent("DoubleEG_Run2015C_25ns_16Dec"      , "/DoubleEG/Run2015C_25ns-16Dec2015-v1/MINIAOD"      , "CMS", ".*root", json)
MuonEG_Run2015C_25ns_16Dec         = kreator.makeDataComponent("MuonEG_Run2015C_25ns_16Dec"        , "/MuonEG/Run2015C_25ns-16Dec2015-v1/MINIAOD"        , "CMS", ".*root", json)
DoubleMuon_Run2015C_25ns_16Dec     = kreator.makeDataComponent("DoubleMuon_Run2015C_25ns_16Dec"    , "/DoubleMuon/Run2015C_25ns-16Dec2015-v1/MINIAOD"    , "CMS", ".*root", json)


### ----------------------------- Run2015D-16Dec2015 ----------------------------------------
## https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset+dataset%3D%2F*%2F*2015D*16Dec2015*%2F*+status%3D*

JetHT_Run2015D_16Dec          = kreator.makeDataComponent("JetHT_Run2015D_16Dec"         , "/JetHT/Run2015D-16Dec2015-v1/MINIAOD"         , "CMS", ".*root", json)
HTMHT_Run2015D_16Dec          = kreator.makeDataComponent("HTMHT_Run2015D_16Dec"         , "/HTMHT/Run2015D-16Dec2015-v1/MINIAOD"         , "CMS", ".*root", json)
MET_Run2015D_16Dec            = kreator.makeDataComponent("MET_Run2015D_16Dec"           , "/MET/Run2015D-16Dec2015-v1/MINIAOD"           , "CMS", ".*root", json)
SingleElectron_Run2015D_16Dec = kreator.makeDataComponent("SingleElectron_Run2015D_16Dec", "/SingleElectron/Run2015D-16Dec2015-v1/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2015D_16Dec     = kreator.makeDataComponent("SingleMuon_Run2015D_16Dec"    , "/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD"    , "CMS", ".*root", json)
SinglePhoton_Run2015D_16Dec   = kreator.makeDataComponent("SinglePhoton_Run2015D_16Dec"  , "/SinglePhoton/Run2015D-16Dec2015-v1/MINIAOD"  , "CMS", ".*root", json)
DoubleEG_Run2015D_16Dec       = kreator.makeDataComponent("DoubleEG_Run2015D_16Dec"      , "/DoubleEG/Run2015D-16Dec2015-v2/MINIAOD"      , "CMS", ".*root", json)
MuonEG_Run2015D_16Dec         = kreator.makeDataComponent("MuonEG_Run2015D_16Dec"        , "/MuonEG/Run2015D-16Dec2015-v1/MINIAOD"        , "CMS", ".*root", json)
DoubleMuon_Run2015D_16Dec     = kreator.makeDataComponent("DoubleMuon_Run2015D_16Dec"    , "/DoubleMuon/Run2015D-16Dec2015-v1/MINIAOD"    , "CMS", ".*root", json)
ZeroBias_Run2015D_16Dec       = kreator.makeDataComponent("ZeroBias_Run2015D_16Dec"       , "/ZeroBias/Run2015D-16Dec2015-v1/MINIAOD"     , "CMS", ".*root", json, useAAA=True)
