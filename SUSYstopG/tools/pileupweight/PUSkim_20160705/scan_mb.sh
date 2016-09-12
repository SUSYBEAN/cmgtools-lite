#!/bin/sh

mbs="69100 69200 69300 69400 69500 69600 69700 69800 69900 70000 70100 70200 70300 70400 70500"
for mb in $mbs;
do

pileupCalc.py -i Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt --inputLumiJSON pileup_latest_20160705.txt --calcMode true  --minBiasXsec ${mb} --maxPileupBin 100 --numPileupBins 100  pileup_DATA_80x_271036-275125_${mb}.root

  ./mcpileup80x.py --tag "271036-275125_${mb}"
done
