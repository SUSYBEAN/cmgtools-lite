#!/bin/sh

tag=`date +%F-%H-%M-%S`
mkdir -p old_${tag}
mv Loop* old_${tag}/


python $CMSSW_BASE/src/PhysicsTools/HeppyCore/python/framework/looper.py pycfg.py config.pck --options=options.json &> local.log 
cp -rf Loop/* . 


