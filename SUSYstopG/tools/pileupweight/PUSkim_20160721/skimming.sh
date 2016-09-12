#!/bin/sh


# selection
selection="(1)"
#selection="((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5)&&(llnunu_l1_pt>100.0&&abs(llnunu_l1_mass-91.1876)<20.0&&met_pt>0.0))"
#selection="((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5)&&(llnunu_l1_pt>0.0&&abs(llnunu_l1_mass-91.1876)<20.0&&met_pt>0.0))"


# compile
g++ skimming.cc -o skimming.exe `root-config --cflags` `root-config --libs`

#inputs
#inputdir=/data2/XZZ2/80X_20160724
#outputdir=/data2/XZZ2/80X_20160724_Skim
#inputdir=80X_20160721_LinksForSkim_v2
#outputdir=/data/XZZ/80X_20160721_Skim
inputdir=/datag/heli/XZZ/80X_20160721_LinksForSkim
outputdir=/home/heli/XZZ/80X_20160721_SkimV2
#outputdir=/dataf/heli/XZZ/80X_20160721_Skim
mkdir -p ${outputdir}

njob="0"

#for infile in $inputdir/*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/Bulk*/vvTreeProducer/tree.root ; 
#for infile in $(ls $inputdir/Single*/vvTreeProducer/tree.root ); 
#for infile in $(ls $inputdir/Single*/vvTreeProducer/tree.root ); 
#for infile in $inputdir/DYJetsToLL_M50/vvTreeProducer/tree.root ; 
for infile in $inputdir/*/vvTreeProducer/tree.root ;
do
  echo "+++ skimming $infile +++"
  outfile="${outputdir}/${infile/$inputdir\//}"
  outfile="${outfile/\/vvTreeProducer\/tree/}"

  # if do PUscan using Zjets:
  #outfile="${outfile/\/vvTreeProducer\/tree/_PUScanV2}"

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
  echo -- Command: ./skimming.exe $infile $outfile $AllEvents $SumWeights $selection

  ./skimming.exe $infile $outfile $AllEvents $SumWeights $selection &> ${outfile}.skim.log &

  njob=$(( njob + 1 ))
  if [ "$njob" -eq "100" ]; then
   # wait
    njob="0"
  fi

done

