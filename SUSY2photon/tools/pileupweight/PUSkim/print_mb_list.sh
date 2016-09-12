#!/bin/sh


echo `ls pileup_MC_80x*.root | sed "s/pileup_MC_80x_/\"/g" | sed "s/.root/\",/g" `

echo `ls pileup_MC_80x*.root | sed "s/pileup_MC_80x_//g" | sed "s/.root/,/g" `


