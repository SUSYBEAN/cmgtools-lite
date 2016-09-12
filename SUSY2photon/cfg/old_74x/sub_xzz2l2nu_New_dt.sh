#!/bin/sh

heppy_batch.py -o LSF_dt run_xzz2l2nu_New_cfg_dt.py -b 'bsub -q 1nd < ./batchScript.sh'

