#!/bin/sh
mbs="62525 "
#mbs="62127 "
#mbs="61000 61100 61200 61300 61400 61500 61600 61700 61800 61900 62000 62100 62200 62300 62400 62500 62600 62700 62800 62900 63000 "
#mbs="59000 "
#mbs="57500 58000 58500 59000 59500 60000 60500 61000 61500 62000 62500 63000 63500 64000 64500 65000 65500 66000 66500 67000 67500 68000 68500 69000 69500 70000 70500 71000 71500 72000 72500 73000 73500 74000 74500 75000 75500 76000 76500 77000 77500 " 
#mbs="69100 69200 69300 69400 69500 69600 69700 69800 69900 70000 70100 70200 70300 70400 70500"

out="results"
mkdir -p $out

for mb in $mbs;
do
  echo "Scanning mb xsec = $mb ub :"
  echo "  - calculating data pileup profile ... "
  pileupCalc.py -i \
    Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt \
    --inputLumiJSON \
    pileup_JSON2016_972016_v3.txt \
    --calcMode true  \
    --minBiasXsec ${mb} \
    --maxPileupBin 100 \
    --numPileupBins 100 \
    ${out}/pileup_DATA_80x_271036-275783_${mb}.root &>  ${out}/pileup_DATA_80x_271036-275783_${mb}.log &

done

#for mb in $mbs;
#do
#  echo "Scanning mb xsec = $mb ub :"
#  echo "  - calculating Data/MC pileup reweighting ... "
#  ./mcpileup80x.py -l -b -q \
#    --input_tag="$out/pileup_DATA_80x_271036-275783_${mb}" \
#    --output_tag="$out/pileup_MC_80x_271036-275783_${mb}" 
#done
