#!/bin/sh


g++ killdup.cc -o killdup.exe `root-config --cflags` `root-config --libs`

#inputs
infile=/home/heli/XZZ/80X_20160825/SingleEMU_Run2016BCD_PromptReco/vvTreeProducer/tree.root
outfile=/home/heli/XZZ/80X_20160825/SingleEMU_Run2016BCD_PromptReco/vvTreeProducer/tree_killdup.root


echo -- Command: ./killdup.exe $infile $outfile
./killdup.exe $infile $outfile

