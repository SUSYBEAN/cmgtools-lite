#!/bin/env python

import os

cardtag = '2l2nu_both'
#cardtag = '2l2nusig1pb_both'
masses = [600,
800,1000,1200,1400,
         1600,1800,
         2000,2500,3000,
         3500,4000,
        #4500
         ]
#tag='2l2nuPvalue'
#tag='2l2nuSignif'
#tag='2l2nuLimitXsecExpSig1'
tag='2l2nuLimit'
#tag='2l2nuLimitExpSig1'
options = '-M Asymptotic --run blind --rMax 10000  '
#options = '-M Asymptotic --run blind --rMax 0.2  '
#options = '-M ProfileLikelihood --significance --pvalue  '
#options = '-M ProfileLikelihood --significance --pvalue -t -1 --expectSignal 1 '
#options = '-M ProfileLikelihood --significance -t -1 --expectSignal 1 '
#options = '-M Asymptotic --rMax=0.1'
#options = '-M Asymptotic --rMax=10000'
#options = '-M Asymptotic -t -1 --expectSignal 1 --rMax=0.2'
#options = '-M Asymptotic -t -1 --expectSignal 1 '
#options = '-M Asymptotic -t -1 --expectSignal 0 '
#options = '-M Asymptotic -t -1 --expectSignal 1 --rMax=10000'
#options = '-M Asymptotic -t -1 --expectSignal 0 --rMax=10000'
#options = '-M Asymptotic -t -1 --expectSignal 1 --minimizerTolerance 0.0001'
for mass in masses:
    card=cardtag+'_m'+str(mass)+'.txt'
    cmd='combine -n '+tag+' -m {mass} {options} {card}\n'.format(mass=mass, options=options, card=card)
    print '### Run combine mass point '+str(mass)+' ###'
    print 'command: '+cmd
    os.system(cmd)

# merge root file
cmd = 'rm -rf higgsCombine'+tag+'.Asymptotic.root'
os.system(cmd)
cmd = 'hadd  higgsCombine'+tag+'.Asymptotic.root  higgsCombine'+tag+'.Asymptotic.mH*.root'
os.system(cmd)

