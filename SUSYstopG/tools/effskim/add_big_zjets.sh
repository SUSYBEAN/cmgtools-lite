#!/bin/sh

dir=/home/heli/XZZ/80X_20160818_light_Skim_EffSkim
fout=DYJetsToLL_M50_BIG_RecoilGraph.root
f1=DYJetsToLL_M50_RecoilGraph.root
f2=DYJetsToLL_M50_MGMLM_Ext1_RecoilGraph.root

rm  $dir/$fout
#hadd -f9 $dir/$fout $dir/$f1 $dir/$f2

hadd  $dir/$fout $dir/$f1 $dir/$f2
