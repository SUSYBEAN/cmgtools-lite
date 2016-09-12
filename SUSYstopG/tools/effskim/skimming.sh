#!/bin/sh


# selection
selection="(1)"
#selection="((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5)&&(llnunu_l1_pt>100.0&&abs(llnunu_l1_mass-91.1876)<20.0&&met_pt>0.0))"
#selection="((abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.1&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4)||(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5)&&(llnunu_l1_pt>0.0&&abs(llnunu_l1_mass-91.1876)<20.0&&met_pt>0.0))"


# compile
g++ skimming.cc -o skimming.exe `root-config --cflags` `root-config --libs`


alias eos='/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'

#inputs
inputdir=/home/heli/XZZ/80X_20160818_light_Skim
outputdir=/home/heli/XZZ/80X_20160818_light_Skim_EffSkim
#inputdir=/data/XZZ/80X_20160721_SkimRecoilOnlyMC
#outputdir=/home/heli/XZZ/80X_20160721_SkimV2_EffSkim
#outputdir=/data2/XZZ2/80X_20160721_SkimRecoilOnlyMC_EffSkim
#outputdir=/data2/XZZ2/80X_20160721_SkimRecoilOnlyMC_EffSkim
#inputdir=/data/XZZ/80X_20160721_SkimV4
#outputdir=/data2/XZZ2/80X_20160721_SkimV4_EffSkim
#inputdir=/data/XZZ/80X_20160721_SkimV3
#outputdir=/data2/XZZ2/80X_20160721_SkimV3_EffSkim
#inputdir=/data2/XZZ2/80X_20160721_SkimV2
#outputdir=/data/XZZ/80X_20160721_SkimV2_EffSkim
#inputdir=/eos/cms/store/caf/user/heli/XZZ/80X_20160721_Skim
#outputdir=/eos/cms/store/caf/user/heli/XZZ/80X_20160721_EffSkim_v3
#inputdir=/home/heli/XZZ/80X_20160721_SkimV2
#outputdir=/datab/yanchu/XZZ/80X_20160721_EffSkim_v2
#inputdir=/datab/yanchu/XZZ/80X_20160721_EffSkim_v2
#outputdir=/home/heli/XZZ/80X_20160721_EffSkim_v2
#outputdir=/dataf/heli/XZZ/80X_20160721_EffSkim_v2
#inputdir=/data/XZZ/80X_20160721_Skim
#outputdir=/data2/XZZ2/80X_20160721_EffSkim_v3
#inputdir=/dataf/heli/XZZ/80X_20160705_L9p17_Skim
#outputdir=/datag/heli/XZZ/80X_20160705_L9p17_EffSkim
#inputdir=/dataf/heli/XZZ/80X_20160705_L6p26_Skim
#outputdir=/datad/heli/XZZ/80X_20160705_L6p26_Skim_ZpTLepTSkim
#inputdir=/data/XZZ/76X_Ntuple/76X_20160514
#inputdir=/data/XZZ/80X_Ntuple/80X_20160705_Skim
#inputdir=/data/mewu/76X_new
#inputdir=/afs/cern.ch/work/m/mewu/public/76X_new
#outputdir=AnalysisRegion
#outputdir=/data2/yanchu/80X_Ntuple/80X_20160705_ZPTSkim

mkdir -p ${outputdir}
#eos mkdir -p ${outputdir}

njob="0"

#for infile in $(ls $inputdir/*.root|grep -v 2016B|grep -v DYJetsToLL|grep -v BulkGrav);
#for infile in $(eos ls $inputdir/*.root | grep root |grep  Run2016 | grep killdup );
#for infile in $(ls $inputdir/*.root|grep -v Run2016);
#for infile in $(ls $inputdir/*.root | grep root |grep  Run2016 | grep killdup );
#for infile in $(ls $inputdir/*.root | grep root |grep  Run2016BCD );
#for infile in $inputdir/DYJetsToLL_M50_MGMLM_Ext1.root;
#for infile in $inputdir/DYJetsToLL_M50.root;
#for infile in $inputdir/DYJetsToLL_M50_RecoilSmooth.root;
#for infile in $inputdir/DYJetsToLL_M50_Recoil*.root;
#for infile in  $inputdir/*.root ;
#for infile in $inputdir/DYJetsToLL_M50_MGMLM_Ext1_RecoilSmooth.root;
#for infile in $inputdir/DYJetsToLL_M50_MGMLM_Ext1_RecoilGraph.root;
#for infile in $inputdir/DYJetsToLL_M50_MGMLM_Ext1_RecoilNoSmooth.root;
#for infile in $inputdir/DYJetsToLL_M50_RecoilGraph.root;
#for infile in $inputdir/SingleEMU_Run2016BCD_PromptReco.root;
#for infile in $(ls $inputdir/*.root|grep -v Single | grep -v DYJets);
for infile in $inputdir/DYJetsToLL_*_RecoilGraph.root;
do
#  echo "+++ skimming $infile +++"
  outfile="${outputdir}/${infile/$inputdir\//}"
  outfile="${outfile/\/vvTreeProducer\/tree/}"

  #infile="root://eoscms/"${inputdir}/${infile} 
  #outfile="root://eoscms/"${outfile}

  echo -- Command: ./skimming.exe $infile $outfile 

  
  ./skimming.exe $infile $outfile > ${outfile}.log &

  njob=$(( njob + 1 ))
  if [ "$njob" -eq "100" ]; then
    wait
    njob="0"
  fi
done

