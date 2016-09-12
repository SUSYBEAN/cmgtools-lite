#!/bin/sh

tag="SingleElectron_Run2016B_PromptReco_v2_Chunk"

njob="0"
for dir in $tag*; 
do
  echo $dir
  cd $dir
  rm -rf Loop*
  log=`grep processed log.txt`; 
  if [ -z "$log" ]; then 
    echo "run the job now" ; 
    python $CMSSW_BASE/src/PhysicsTools/HeppyCore/python/framework/looper.py pycfg.py config.pck --options=options.json &> local.log && mv Loop/* . &
  else 
    echo "processed before"; 
  fi
  cd ../
  njob=$(( njob + 1 ))
  if [ "$njob" -eq "6" ]; then
    wait
    njob="0"
  fi
done

