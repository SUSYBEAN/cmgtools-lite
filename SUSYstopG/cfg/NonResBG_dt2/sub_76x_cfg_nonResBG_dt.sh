#!/bin/sh

heppy_batch.py -o LSF run_76x_cfg_nonResBG_dt.py -b 'bsub -q 8nh < ./batchScript.sh'

