#!/bin/sh



# compile
g++ preskim.cc -o preskim.exe `root-config --cflags` `root-config --libs`

#inputs
inputdir=/data2/XZZ2/80X_20160818
outputdir=/home/heli/XZZ/80X_20160818
for dd in ${inputdir}/*/vvTreeProducer;
do 
  infile="${dd}/tree.root";
  oo="${dd/$inputdir/$outputdir}";
  outfile="${oo}/tree_light.root";
  echo mkdir -p $oo;
  mkdir -p $oo ;
  echo "./preskim.exe $infile $outfile &> ${outfile}.log & "
  ./preskim.exe $infile $outfile &> ${outfile}.log &
done
