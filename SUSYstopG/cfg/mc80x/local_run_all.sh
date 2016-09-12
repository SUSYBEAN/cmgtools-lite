#!/bin/sh

listdd="
ZZTo2L2Nu_Chunk0
ZZTo2L2Nu_Chunk1
ZZTo2L2Nu_Chunk10
ZZTo2L2Nu_Chunk11
ZZTo2L2Nu_Chunk12
ZZTo2L2Nu_Chunk13
ZZTo2L2Nu_Chunk14
ZZTo2L2Nu_Chunk15
ZZTo2L2Nu_Chunk16
ZZTo2L2Nu_Chunk17
ZZTo2L2Nu_Chunk18
ZZTo2L2Nu_Chunk19
ZZTo2L2Nu_Chunk2
ZZTo2L2Nu_Chunk20
ZZTo2L2Nu_Chunk3
ZZTo2L2Nu_Chunk4
ZZTo2L2Nu_Chunk5
ZZTo2L2Nu_Chunk6
ZZTo2L2Nu_Chunk7
ZZTo2L2Nu_Chunk8
ZZTo2L2Nu_Chunk9
"
script="local_run_one_job.sh"

njob="0"

for dd in $listdd;
do
  cd $dd
  cp ../$script .
  echo "submit $dd"
  ./$script &> log &
  cd ../


  njob=$(( njob + 1 ))
  if [ "$njob" -eq "8" ]; then
    wait
    njob="0"
  fi

done 


