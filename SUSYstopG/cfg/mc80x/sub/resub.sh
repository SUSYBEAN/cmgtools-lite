#!/bin/sh

dir="mc1"
queue="2nd"

# need to check your jobs with the expected n root files and n pck files
# to verify if the job is finished sucessfully.
n_root_files="3"
n_pck_files="14"

cd $dir
list=`ls -1`
for dd in $list;
do
  echo "### check job $dd"
  n1=`rootls ${dd}/*/*.root | grep root | wc -l`;
  n2=`ls -l ${dd}/* | grep pck |wc -l`;
  if [ "$n1" -eq "$n_root_files" -a  "$n2" -eq "$n_pck_files" ];
  then
    echo "job is done with ${n1} root files and ${n2} pck files.";
  else
    rr=`bjobs -o "name" -J ${dd} | grep -v JOB_NAME | grep ${dd}`;
    if [ "$rr" == "$dd" ]; 
    then
      echo "  job $dd is still running";
    else 
      echo " submit $dd";
      cd $dd ;
      echo "  bsub -q $queue -J $dd  < batchScript.sh" ;
      bsub -q $queue -J $dd  < batchScript.sh ;
      cd ../ ;
    fi;
  fi;

done 

cd ../

