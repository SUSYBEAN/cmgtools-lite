#!/bin/sh

./vvMakeLimitPlot.py -b 1 -o SigmaGZZPlot --minX 500 --maxX 4000 --minY 0.001 --maxY 100 --titleY "95% C.L. limit on #sigma(BulkG#rightarrowZZ) (pb)" --scale 25 --drawXsec 1  results/higgsCombinexzz2l2nu.root
#./vvMakeLimitPlot.py -b 1 -o limitExpSig0 --minX 500 --maxX 4000 --minY 0.1 --maxY 10000000 --titleY "#sigma_{95%}/#sigma_{th.}" --drawXsec 0  results/higgsCombinexzz2l2nu_zptcut100_metcut200.root
#./vvMakeLimitPlot.py -b 1 -o SigmaGZZPlot --minX 500 --maxX 4000 --minY 0.001 --maxY 100 --titleY '95% C.L. limit on #sigma(BulkG#rightarrowZZ) (pb)' --scale 25 --drawXsec 1 higgsCombine2l2nuLimit.Asymptotic.root
#./vvMakeLimitPlot.py -b 1 -o limitExpSig0 --minX 500 --maxX 4000 --minY 0.1 --maxY 10000000 --titleY "#sigma_{95%}/#sigma_{th.}" --drawXsec 0  higgsCombine2l2nuLimit.Asymptotic.root 
#./vvMakeLimitPlot.py -b 1 -o limitExpSig1 --minX 500 --maxX 3500 --minY 1 --maxY 10000000 --titleY="#sigma_{95%}/#sigma_{th.}" higgsCombine2l2nuLimitExpSig1.Asymptotic.root 
#./vvMakeLimitPlot.py -b 0 -o SigmaPlotExpSig1 --minX 500 --maxX 3500 --minY 0.0001 --maxY 1 --titleY '95% limit on #sigma(X#rightarrowZZ#rightarrow ll#nu#nu) (pb)' higgsCombine2l2nuLimitXsecExpSig1.Asymptotic.root
#./vvMakeLimitPlot.py -b 0 -o SigmaPlot --minX 500 --maxX 3500 --minY 0.0001 --maxY 1 --titleY '95% limit on #sigma(X#rightarrowZZ#rightarrow ll#nu#nu) (pb)' higgsCombine2l2nuLimitXsec.Asymptotic.root
#./vvMakeLimitPlot.py -b 0 -o SigmaGZZPlot --minX 500 --maxX 4500 --minY 0.005 --maxY 50 --titleY '95% limit on #sigma(BulkG#rightarrowZZ) (pb)' --scale 25 higgsCombine2l2nuLimitXsec.Asymptotic.root

