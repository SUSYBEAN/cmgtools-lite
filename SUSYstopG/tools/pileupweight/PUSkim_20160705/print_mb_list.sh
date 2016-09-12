#!/bin/sh


echo `ls pileup_MC_80x*00.root | sed "s/pileup_MC_80x_271036-275125_/\"/g" | sed "s/.root/\",/g" `

echo `ls pileup_MC_80x*00.root | sed "s/pileup_MC_80x_271036-275125_//g" | sed "s/.root/,/g" `


