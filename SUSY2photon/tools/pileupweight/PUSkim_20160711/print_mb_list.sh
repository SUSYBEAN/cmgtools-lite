#!/bin/sh

#cd results
#echo `ls pileup_MC_80x*00.root | sed "s/pileup_MC_80x_271036-275783_/\"/g" | sed "s/.root/\",/g" `
#echo `ls pileup_MC_80x*00.root | sed "s/pileup_MC_80x_271036-275783_//g" | sed "s/.root/,/g" `
#cd ../


mbs=" 61000 61100 61200 61300 61400 61500 61600 61700 61800 61900 62000 62100 62200 62300 62400 62500 62600 62700 62800 62900 63000 "
echo `echo $mbs | sed "s/ /\", \"/g"`
echo `echo $mbs | sed "s/ /, /g"`

