#!/bin/sh


pileupCalc.py -i Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt --inputLumiJSON pileup_latest_20160617.txt --calcMode true  --minBiasXsec 69735 --maxPileupBin 100 --numPileupBins 100  pileup_DATA_80x_2016B_Run_271036-274443_mb_69735.root

./mcpileup80x.py --tag "2016B_Run_271036-274443_mb_69735" 
