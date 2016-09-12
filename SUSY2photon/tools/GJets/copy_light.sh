#!/bin/sh

indir=/data2/XZZ2/80X_20160818
outdir=/home/heli/XZZ/80X_20160818_light

for dd in  ${indir}/* ; 
do 
  echo $dd;  
  ddo=${dd/$indir/$outdir}
  echo "mkdir -p ${ddo} " 
  mkdir -p ${ddo}  

  for dd2 in ${dd}/* ; 
  do 
    if [[ ${dd2} != *"vvTree"* ]]; then
      #echo $dd2;
      ddo2="${dd2/$indir/$outdir} "
      echo "cp -rp $dd2 $ddo2 " 
      cp -rp $dd2 $ddo2 
    fi
  done; 

  ttin="${ddo/_light/}/vvTreeProducer"
   
  echo "cp -rp $ttin $ddo/"
  cp -rp $ttin $ddo/
  echo "mv $ddo/vvTreeProducer/tree_light.root $ddo/vvTreeProducer/tree.root"
  mv $ddo/vvTreeProducer/tree_light.root $ddo/vvTreeProducer/tree.root
done
