#!/bin/sh
# --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json
lumiCalc2.py lumibylsXing \
--xingMinLum 0.1 \
-b stable \
--normtag normtag_DATACERT.json  \
-i json_DCSONLY.txt \
-o lumi_DCSONLY.csv

estimatePileup_makeJSON.py --csvInput lumi_DCSONLY.csv pileup_JSON.txt
