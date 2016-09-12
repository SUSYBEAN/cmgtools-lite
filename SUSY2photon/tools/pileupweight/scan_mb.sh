#!/bin/sh

mbs="697 713"
#mbs="691 692 693 694 695 696 697 698 699 700 701 702 703 704 705"
#mbs="575 580 585 590 595 600 605 610 615 620 625 630 635 640 645 650 655"
#mbs="660 665 670 675 680 685 690 695 700 705 710 715 720 725 730 735 740 745 750 755 760 765"
for mb in $mbs;
do

pileupCalc.py -i Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt --inputLumiJSON pileup_latest_20160617.txt --calcMode true  --minBiasXsec ${mb}00 --maxPileupBin 100 --numPileupBins 100  pileup_DATA_80x_271036-274443_${mb}00.root

#  pileupCalc.py -i Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec ${mb}00 --maxPileupBin 50 --numPileupBins 50  pileup_DATA_80x_${mb}.root

  ./mcpileup80x.py --tag "271036-274443_${mb}00"
done
