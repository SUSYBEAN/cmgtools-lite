#!/bin/sh

mkdir -p ../not_finished

for dd in `ls`; 
do 
  echo $dd;
  n1=`rootls ${dd}/*/*.root | grep root | wc -l`;  
  n2=`ls -l ${dd}/* | grep pck |wc -l`;
  if [ "$n1" -eq "4" -a  "$n2" -eq "16" ] ; then
    echo "job is done with ${n1} root files and ${n2} pck files."; 
  else
    mv ${dd} ../not_finished/ 
  fi;  

done
