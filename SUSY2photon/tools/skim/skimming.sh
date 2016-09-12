#!/bin/sh


# selection
selection="(1)"
#selection="((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5)&&(llnunu_l1_pt>100.0&&abs(llnunu_l1_mass-91.1876)<20.0&&met_pt>0.0))"
#selection="((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5)&&(llnunu_l1_pt>0.0&&abs(llnunu_l1_mass-91.1876)<20.0&&met_pt>0.0))"


# compile
g++ skimming.cc -o skimming.exe `root-config --cflags` `root-config --libs`

#inputs
inputdir=/datag/heli/XZZ/80X_20160716
outputdir=/dataf/heli/XZZ/80X_20160716_Skim
#inputdir=/datab/heli/XZZ/80X_20160716
#outputdir=/dataf/heli/XZZ/80X_20160716_Skim
#inputdir=/data/XZZ/80X_Ntuple/80X_20160621
#outputdir=/home/heli/work/XZZ/80X_Ntuple/80X_20160621_Skim
#inputdir=/data/XZZ/80X_Ntuple/80X_20160618
#outputdir=/home/heli/work/XZZ/80X_Ntuple/80X_20160618_Skim
#outputdir=/data/XZZ/80X_Ntuple/80X_20160618_Skim
#inputdir=/data/XZZ/80X_Ntuple/80X_20160606_NoHLT
#inputdir=/afs/cern.ch/user/h/heli/work/private/cms/2tev/XZZ2l2nu_80X/CMSSW_8_0_10/src/CMGTools/XZZ2l2nu/cfg/mc80x/sub/mc1
#outputdir=/data/XZZ/80X_Ntuple/80X_20160606_NoHLT_Skim
#outputdir=/home/heli/work/XZZ/80X_Ntuple/80X_20160606_NoHLT_Skim
#outputdir=/home/heli/work/XZZ/80X_Ntuple/
mkdir -p ${outputdir}

#for infile in $inputdir/*/vvTreeProducer/tree.root ; 
for infile in $inputdir/Single*/vvTreeProducer/tree.root ; 
do
  echo "+++ skimming $infile +++"
  outfile="${outputdir}/${infile/$inputdir\//}"
  outfile="${outfile/\/vvTreeProducer\/tree/}"

  inSkimFile=${infile/vvTreeProducer\/tree.root/skimAnalyzerCount\/SkimReport.txt}

  #echo $inSkimFile
  AllEvents=`grep "All Events" ${inSkimFile} | awk {'print $3'}`
  SumWeights=`grep "Sum Weights" ${inSkimFile} | awk {'print $3'}`

  if [ -z $AllEvents ]; then
    echo Fail to get All Events from file ${inSkimFile}
    continue
  fi
  if [ -z $SumWeights ]; then
    SumWeights=$AllEvents
  fi

  echo -- Input file: $infile
  echo -- Output file: $outfile
  echo -- AllEvents: $AllEvents , SumWeights: $SumWeights
  echo -- Selection: $selection
  echo -- Command:   ./skimming.exe $infile $outfile $AllEvents $SumWeights $selection
  ./skimming.exe $infile $outfile $AllEvents $SumWeights $selection &> $outfile.log &

done

