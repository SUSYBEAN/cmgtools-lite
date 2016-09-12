#!/bin/sh


# selection
selection="(1)"
#selection="((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5)&&(llnunu_l1_pt>100.0&&abs(llnunu_l1_mass-91.1876)<20.0&&met_pt>0.0))"
#selection="((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5)&&(llnunu_l1_pt>0.0&&abs(llnunu_l1_mass-91.1876)<20.0&&met_pt>0.0))"


# compile
g++ -c KalmanMuonCalibrator.cc -o KalmanMuonCalibrator.o `root-config --cflags` -I.
g++ -c skimming.cc -o skimming.o `root-config --cflags` -I. 
g++ KalmanMuonCalibrator.o skimming.o -o skimming.exe `root-config --cflags` `root-config --libs` -lMinuit2 -I. -L.

#inputs
inputdir=/home/heli/XZZ/80X_20160818_light
outputdir=/home/heli/XZZ/80X_20160818_light_Skim
#outputdir=/home/heli/XZZ/80X_20160818_light_Skim
#inputdir=/data2/XZZ2/80X_20160721_LinksForSkimV3
#outputdir=/data/XZZ/80X_20160721_SkimRecoilOnlyMC
#outputdir=/data/XZZ/80X_20160721_SkimV3
#outputdir=/data/XZZ/80X_20160721_SkimV4
#inputdir=/data2/XZZ2/80X_20160724
#outputdir=/data2/XZZ2/80X_20160724_Skim
#inputdir=80X_20160721_LinksForSkim_v2
#outputdir=/data/XZZ/80X_20160721_Skim
#inputdir=/datag/heli/XZZ/80X_20160721_LinksForSkim
#outputdir=/home/heli/XZZ/80X_20160721_SkimV2
#outputdir=/dataf/heli/XZZ/80X_20160721_Skim
mkdir -p ${outputdir}

njob="0"

#for infile in $inputdir/*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/Bulk*/vvTreeProducer/tree.root ; 
#for infile in $(ls $inputdir/Single*/vvTreeProducer/tree.root ); 
#for infile in $(ls $inputdir/Single*/vvTreeProducer/tree.root ); 
#for infile in $inputdir/*/vvTreeProducer/tree.root ;
#for infile in $inputdir/DYJetsToLL_M50_MGMLM/vvTreeProducer/tree.root ; 
#for infile in $inputdir/DYJetsToLL_M50_MGMLM_Ext1/vvTreeProducer/tree.root ; 
#for infile in $inputdir/Sing*/vvTreeProducer/tree.root ; 
#for infile in $(ls $inputdir/*/vvTreeProducer/tree.root | grep -v Single | grep -v DYJets ); 
#for infile in $inputdir/DYJetsToLL_M50/vvTreeProducer/tree.root ; 
for infile in $inputdir/DYJetsToLL*/vvTreeProducer/tree.root ; 
do
  echo "+++ skimming $infile +++"
  outfile="${outputdir}/${infile/$inputdir\//}"

  # options for outputs
  #outfile="${outfile/\/vvTreeProducer\/tree/}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_NoRecoil}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_RecoilNoSmooth}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_RecoilSmooth}"
  outfile="${outfile/\/vvTreeProducer\/tree/_RecoilGraph}"

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

