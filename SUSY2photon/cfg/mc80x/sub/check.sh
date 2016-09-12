#!/bin/sh

# ask Hengne.Li@cern.ch for questions.

# Warning!! 
# always test before really run this script, 
# it has potential to delete your good and finished job folders!!!!

# privide chunks running directory, 
# and directory to copy jobs out

dir=mc1
#out=/home/heli/80X_20160825_Chunks
out=/data2/XZZ2/80X_20160825_Chunks

mkdir -p $out

# need to check your jobs with the expected n root files and n pck files
# to verify if the job is finished sucessfully.
n_root_files="3"
n_pck_files="14"


if [ ! -e "$out" ]; then
  echo "ERROR:: Do not exist output directory $out, exist... "
  exit 0
fi

#echo "#### first do a overall sync"
#rsync -var $dir $out/

echo "#### go to check each directory "
cd $dir

list=`ls -1`

for job in $list;
do
  echo "## check $job "
  n1=`rootls ${job}/*/*.root | grep root | wc -l`;
  n2=`ls -l ${job}/* | grep pck |wc -l`;

  if [ "$n1" -eq "$n_root_files" -a  "$n2" -eq "$n_pck_files" ];
  then
    echo "- job is done correctly with ${n1} root files and ${n2} pck files."
    echo "  - copy out and delete .. "
    echo " > rsync -var $job $out/$dir/"
    rsync -var $job $out/$dir/
    echo " > rm -rf $job"
    rm -rf $job
  else
    echo "- job is not finished or has problem to be resubmitted .. "
  fi

done

cd ../


