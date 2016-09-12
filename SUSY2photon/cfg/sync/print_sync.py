#!/usr/bin/env python
import os, sys, re, optparse,pickle,shutil


channel='mu'
#channel='el'
tag='_v2'
dir='sync_'+channel+tag
samples = {
           'M-0800':'BulkGravToZZToZlepZhad_narrow_800',
           'M-1000':'BulkGravToZZToZlepZhad_narrow_1000',
           'M-1200':'BulkGravToZZToZlepZhad_narrow_1200',
           'M-1400':'BulkGravToZZToZlepZhad_narrow_1400',
           'M-1600':'BulkGravToZZToZlepZhad_narrow_1600',
           'M-1800':'BulkGravToZZToZlepZhad_narrow_1800',
           'M-2000':'BulkGravToZZToZlepZhad_narrow_2000',
           'M-2500':'BulkGravToZZToZlepZhad_narrow_2500',
           'M-3000':'BulkGravToZZToZlepZhad_narrow_3000',
           'M-3500':'BulkGravToZZToZlepZhad_narrow_3500',
           'M-4000':'BulkGravToZZToZlepZhad_narrow_4000',
           'M-4500':'BulkGravToZZToZlepZhad_narrow_4500',
          }

class counter:
    def __init__(self, name):
        self.name = name
        
counters = []

sampleNames = sorted(samples.keys())
 
for name in sampleNames:
    sample = samples[name]
    directory = dir+'/'+sample
    count = counter(name)
    pkcount = pickle.load(open(directory+'/skimAnalyzerCount/SkimReport.pck'))
    count.AllEvents = pkcount['All Events'][1] 
    pkcount = pickle.load(open(directory+'/XZZGenAnalyzer/XZZGenReport.pck'))
    count.GenZlep = pkcount['pass events'][1] 
    pkcount = pickle.load(open(directory+'/TriggerBitFilter/TriggerReport.pck'))
    count.HLT = pkcount['Pass Triger Events'][1]
    pkcount = pickle.load(open(directory+'/VertexAnalyzer/GoodVertex.pck'))
    count.Vertex = pkcount['Events With Good Vertex'][1]
    pkcount = pickle.load(open(directory+'/leptonAnalyzer/events.pck'))
    count.LeptonKin = pkcount['pass 1mu kin events' if channel=='mu' else 'pass 1el kin events'][1]
    count.LeptonID = pkcount['pass 1mu kin+id events' if channel=='mu' else 'pass 1el kin+id events'][1]
    count.LeptonIso = pkcount['pass 1mu kin+id+iso events' if channel=='mu' else 'pass 1el kin+id+iso events'][1]
    pkcount = pickle.load(open(directory+'/leptonicVMaker/events.pck'))
    count.Zdaughter = pkcount['pass mu events' if channel=='mu' else 'pass el events'][1]
    count.Zpeak = pkcount['pass events'][1]
   
    counters.append(count)

outfile = open(dir+'_print.txt', 'w')

line = '*'*12
for count in counters: line += '*'*12
line += '*'
print line
outfile.write(line+'\n')
line = '*'+' '*11
for count in counters: line += '*'+count.name.center(11)
line += '*'
print line
outfile.write(line+'\n')
line = '*'*12
for count in counters: line += '*'*12
line += '*'
print line
outfile.write(line+'\n')
line = '* '+'AllEvts'.ljust(10)
for count in counters: line += '*'+str(count.AllEvents).center(11)
line += '*'
print line
outfile.write(line+'\n')
line = '* '+'GenZlep'.ljust(10)
for count in counters: line += '*'+str(count.GenZlep).center(11)
line += '*'
print line
outfile.write(line+'\n')
line = '* '+'HLT'.ljust(10)
for count in counters: line += '*'+str(count.HLT).center(11)
line += '*'
print line
outfile.write(line+'\n')
line = '* '+'Vertex'.ljust(10)
for count in counters: line += '*'+str(count.Vertex).center(11)
line += '*'
print line
outfile.write(line+'\n')
line = '* '+'LeptonKin'.ljust(10)
for count in counters: line += '*'+str(count.LeptonKin).center(11)
line += '*'
print line
outfile.write(line+'\n')
line = '* '+'LeptonID'.ljust(10)
for count in counters: line += '*'+str(count.LeptonID).center(11)
line += '*'
print line
outfile.write(line+'\n')
line = '* '+'LeptonIso'.ljust(10)
for count in counters: line += '*'+str(count.LeptonIso).center(11)
line += '*'
print line
outfile.write(line+'\n')
line = '* '+'Zdaughter'.ljust(10)
for count in counters: line += '*'+str(count.Zdaughter).center(11)
line += '*'
print line
outfile.write(line+'\n')
line = '* '+'Zpeak'.ljust(10)
for count in counters: line += '*'+str(count.Zpeak).center(11)
line += '*'
print line
outfile.write(line+'\n')
line = '*'*12
for count in counters: line += '*'*12
line += '*'
print line
outfile.write(line+'\n')

outfile.close()



