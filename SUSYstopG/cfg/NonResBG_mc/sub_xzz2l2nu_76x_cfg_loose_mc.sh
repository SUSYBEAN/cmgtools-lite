#!/bin/sh

heppy_batch.py -o LSF1 run_76x_cfg_nonResBG_mc.py -b 'bsub -q 1nd < ./batchScript.sh'

