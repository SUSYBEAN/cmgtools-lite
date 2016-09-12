#!/bin/sh


pileupCalc.py -i Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt --inputLumiJSON pileup_latest_20160705.txt --calcMode true  --minBiasXsec 69452 --maxPileupBin 100 --numPileupBins 100  pileup_DATA_80x_2016B_Run_271036-275125_mb_69452.root

./mcpileup80x.py --tag "2016B_Run_271036-275125_mb_69452" 
