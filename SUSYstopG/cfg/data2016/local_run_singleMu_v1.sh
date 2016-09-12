#!/bin/sh

tag="SingleMuon_Run2016B_PromptReco_Chunk"

njob="0"
for dir in $tag*; 
do
  echo $dir
  cd $dir
  rm -rf Loop* 
  python $CMSSW_BASE/src/PhysicsTools/HeppyCore/python/framework/looper.py pycfg.py config.pck --options=options.json &> local.log && mv Loop/* . &
  cd ../
  njob=$(( njob + 1 ))
  if [ "$njob" -eq "6" ]; then
    wait
    njob="0"
  fi
done

