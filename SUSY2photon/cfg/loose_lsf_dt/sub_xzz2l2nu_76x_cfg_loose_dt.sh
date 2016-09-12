#!/bin/sh

heppy_batch.py -o LSF run_xzz2l2nu_76x_cfg_loose_dt.py -b 'bsub -q 8nh < ./batchScript.sh'

