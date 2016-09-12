#!/bin/sh


echo `ls pileup_MC_80x*00.root | sed "s/pileup_MC_80x_271036-274443_/\"/g" | sed "s/00.root/\",/g" `

echo `ls pileup_MC_80x*00.root | sed "s/pileup_MC_80x_271036-274443_//g" | sed "s/00.root/,/g" `


