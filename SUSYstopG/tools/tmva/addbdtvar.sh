#!/bin/sh



# compile
g++ addbdtvar.cc -o addbdtvar.exe `root-config --cflags` -I.  `root-config --libs` -lMLP -lMinuit -lTreePlayer -lTMVA -lTMVAGui -lXMLIO

#inputs
configfile=addbdtvar.cfg
inputdir=76X_JEC_Skim
outputdir=76X_JEC_Skim_BDT
mkdir -p ${outputdir}


#for infile in $inputdir/*.root ; 
for infile in $inputdir/Sing*.root ; 
do
  echo "+++ skimming $infile +++"
  outfile="${outputdir}/${infile/$inputdir\//}"

  echo $outfile
  echo -- Input file: $infile
  echo -- Output file: $outfile
  echo -- Congfit file: $configfile
  echo -- Command: ./addbdtvar.exe $infile $outfile $configfile
 
  tag="${infile/$inputdir\//}"
  ./addbdtvar.exe $infile $outfile $configfile &> addbdtvar_$tag.log &

done

