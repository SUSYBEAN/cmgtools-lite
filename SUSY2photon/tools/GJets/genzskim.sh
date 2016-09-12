#!/bin/sh



# compile
g++ genzskim.cc -o genzskim.exe `root-config --cflags` `root-config --libs`

#inputs
inputdir=/home/heli/XZZ/80X_20160818_light
inputdirgen=/home/heli/XZZ/80X_20160822
#inputdirgen=/data2/XZZ2/80X_20160822
outputdir=/home/heli/XZZ/80X_20160818_lightgen

mkdir -p $outputdir

for dd in ${inputdir}/DYJets*/vvTreeProducer;
do 
  infile="${dd}/tree.root";
  gg="${dd/$inputdir/$inputdirgen}"
  genfile="${gg}/tree.root";
  oo="${dd/$inputdir/$outputdir}";
  outfile="${oo}/tree_gen.root";
  echo mkdir -p $oo;
  mkdir -p $oo ;
  echo "./genzskim.exe $infile $genfile $outfile &> ${outfile}.log & "
  ./genzskim.exe $infile $genfile $outfile &> ${outfile}.log &
done
