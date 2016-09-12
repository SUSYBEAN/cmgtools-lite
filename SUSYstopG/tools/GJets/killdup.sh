#!/bin/sh


g++ killdup.cc -o killdup.exe `root-config --cflags` `root-config --libs`

#inputs
infile=/data/XZZ/80X_20160721_SkimV4/SingleEMU_Run2016BCD_PromptReco.root
outfile=/data2/XZZ2/80X_20160721_SkimV4_EffSkim/SingleEMU_Run2016BCD_PromptReco_killdup.root
#infile=/data2/XZZ2/80X_20160724_Skim/SingleEMU_Run2016BCD_PromptReco.root
#outfile=/data2/XZZ2/80X_20160724_Skim/SingleEMU_Run2016BCD_PromptReco_killdup.root
#infile=/data/XZZ/80X_20160721_Skim/SingleEMU_Run2016BCD_PromptReco.root
#outfile=/data/XZZ/80X_20160721_Skim/SingleEMU_Run2016BCD_PromptReco_killdup.root


#infile=/data2/XZZ/80X_20160721/SingleEMU_Run2016BCD_PromptReco/vvTreeProducer/tree.root
#outfile=/data2/XZZ/80X_20160721/SingleEMU_Run2016BCD_PromptReco/vvTreeProducer/tree_killdup.root


echo -- Command: ./killdup.exe $infile $outfile
./killdup.exe $infile $outfile

