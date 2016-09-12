#!/bin/sh

heppy_batch.py -o mc1 run_xzz2l2nu_80x_cfg_loose_mc.py -b 'bsub -q sssss < ./batchScript.sh'
#heppy_batch.py -o mc1 run_xzz2l2nu_80x_cfg_photon_mc.py -b 'bsub -q ssssss < ./batchScript.sh'

#cd mc1 
#heppy_check.py * -b "bsub -q 2nd"
#cd -
