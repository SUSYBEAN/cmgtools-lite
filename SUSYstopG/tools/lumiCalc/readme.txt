#!/bin/sh

#####
2016-07-14

Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt

 brilcalc lumi -b "STABLE BEAMS" --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt 

 brilcalc lumi -b "STABLE BEAMS" --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt

 
#####
2016-06-17
brilcalc lumi -b "STABLE BEAMS" --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt

2.596909574 fb-1

####

2016:

brilcalc lumi -b "STABLE BEAMS" --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt


json:
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt

#http://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html

export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH

pip install --install-option="--prefix=$HOME/.local" brilws

pip show brilws

brilcalc -h

briltag -h

brilcalc lumi -h

brilcalc beam -h

# calculate int lumi for 2015D
brilcalc lumi --normtag /afs/cern.ch/user/c/cmsbril/public/normtag_json/OfflineNormtagV1.json -i Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2_2015D.txt

# the output is in lumi_2015D.txt, L=2153 pb-1 


# calculate int lumi for re-reco 2015C 25ns + 2015D
# rereco golden json: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt
brilcalc lumi --normtag /afs/cern.ch/user/c/cmsbril/public/normtag_json/OfflineNormtagV1.json -i Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt
# the output is in lumi_rereoco_2015CD25ns.txt, L=2169 pb-1

OfflineNormtagV2.json

update:
https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM 

 3/3/2016 new normtag for 3.8T 50 ns and 25 ns pp data, currently covering data included in the 25ns Muon JSON and 50ns Muon JSON. The uncertainty on the integrated luminosity is 2.7% if this normtag is used, i.e.:.

brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json -i [yourjson] 

brilcalc lumi --normtag  /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json  -i Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt >lumi_rereoco_2015CD25ns_20160303.txt 

2.3183 fb-1

Then update JSON:
brilcalc lumi --normtag  /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json  -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt > lumi_rereoco_2015CD25ns_20160303.txt
This is recommended:
https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2015Analysis
2.318 fb-1

also try v3:
brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/OfflineNormtagV3.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt 

v3 gives : 2.263483265176 fb-1
not yet understood if we can use it.



