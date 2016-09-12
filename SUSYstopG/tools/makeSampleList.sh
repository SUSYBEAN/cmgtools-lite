#!/bin/sh


for data in `das_client.py --query="dataset=/BulkGravToZZ_narrow_M-*/RunIISpring15MiniAODv2-74X*/MINIAODSIM" --limit=0`;
do
   id=`echo ${data} | sed "s/_13TeV-madgraph\/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1\/MINIAODSIM//g" | sed "s/\///g"` 
   echo "${id} = kreator.makeMCComponent(\"${id}\",\"${data}\", \"CMS\", \".*root\",1.0)"
done

