#!/bin/sh

heppy_batch.py -o LSF_mc run_xzz2l2nu_New_cfg_mc.py -b 'bsub -q 1nd < ./batchScript.sh'

