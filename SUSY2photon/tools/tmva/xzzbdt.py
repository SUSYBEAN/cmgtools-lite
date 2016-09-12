#!/usr/bin/env python

import os

scale = 1.0

inputdir = '76X_JEC_Skim'

tags = [
 "m600", 
# "m800", "m1000",
# "m1200",
# "m1400", "m1600", "m1800", "m2000",
# "m2500", "m3000", "m3500", "m4000",
]

files = {
"m600":"BulkGravToZZToZlepZinv_narrow_600.root",
"m800":"BulkGravToZZToZlepZinv_narrow_800.root",
"m1000":"BulkGravToZZToZlepZinv_narrow_1000.root",
"m1200":"BulkGravToZZToZlepZinv_narrow_1200.root",
"m1400":"BulkGravToZZToZlepZinv_narrow_1400.root",
"m1600":"BulkGravToZZToZlepZinv_narrow_1600.root",
"m1800":"BulkGravToZZToZlepZinv_narrow_1800.root",
"m2000":"BulkGravToZZToZlepZinv_narrow_2000.root",
"m2500":"BulkGravToZZToZlepZinv_narrow_2500.root",
"m3000":"BulkGravToZZToZlepZinv_narrow_3000.root",
"m3500":"BulkGravToZZToZlepZinv_narrow_3500.root",
"m4000":"BulkGravToZZToZlepZinv_narrow_4000.root",
"m4500":"BulkGravToZZToZlepZinv_narrow_4500.root",
}

sigXsec = {
'm600'  : 3.44631e-04,
'm800'  : 6.31859e-05,
'm1000' : 1.68661e-05,
'm1200' : 5.59677e-06,
'm1400' : 2.13168e-06,
'm1600' : 8.97713e-07,
'm1800' : 4.06090e-07,
'm2000' : 1.94415e-07,
'm2500' : 3.63496e-08,
'm3000' : 7.95422e-09,
'm3500' : 1.95002e-09,
'm4000' : 5.03747e-10,
'm4500' : 1.0,
}

# compile
cmd = "g++ xzzbdt.cc -o xzzbdt.exe  `root-config --cflags` `root-config --libs` -lMLP -lMinuit -lTreePlayer -lTMVA -lTMVAGui -lXMLIO"
os.system(cmd)

# run
for tag in tags:
    print "running "+tag
    cmd = "./xzzbdt.exe "+tag+" "+inputdir+" "+files[tag]+" "+str(sigXsec[tag])+" "+str(scale)+" &> xzzbdt_"+tag+".log &" 
    print " command "+cmd
    os.system(cmd)






