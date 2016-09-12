#!/bin/sh



# compile
g++ preskim_gjets.cc -o preskim_gjets.exe `root-config --cflags` `root-config --libs`

#inputs
inputdir=/data2/XZZ2/80X_20160810_LinksForLight
outputdir=/home/heli/XZZ/80X_20160810
for dd in ${inputdir}/*/vvTreeProducer;
do 
  infile="${dd}/tree.root";
  oo="${dd/$inputdir/$outputdir}";
  outfile="${oo}/tree_light.root";
  echo mkdir -p $oo;
  mkdir -p $oo ;
  echo "./preskim_gjets.exe $infile $outfile &> ${outfile}.log & "
  ./preskim_gjets.exe $infile $outfile &> ${outfile}.log &
done
