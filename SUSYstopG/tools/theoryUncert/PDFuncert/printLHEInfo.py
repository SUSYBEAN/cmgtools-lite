#!/usr/bin/env python

from DataFormats.FWLite import Handle, Events, Runs
import sys
import os


samples = {
'WWTo2L2Nu':'/store/mc/RunIIFall15MiniAODv2/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000/5E230663-7BBC-E511-962F-18A905C5A684.root',
'WWToLNuQQ':'/store/mc/RunIIFall15MiniAODv2/WWToLNuQQ_13TeV-powheg/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/007D46BD-84D1-E511-8BD5-002590A88802.root',
'WZTo1L1Nu2Q':'/store/mc/RunIIFall15MiniAODv2/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/006191B4-83CA-E511-8231-0CC47A78A30E.root',
'WZTo2L2Q':'/store/mc/RunIIFall15MiniAODv2/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0045FB1F-C2CD-E511-9ADE-0026B937D207.root',
'WZTo3LNu':'/store/mc/RunIIFall15MiniAODv2/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000/022EC2EB-90B8-E511-AED0-0026B937D37D.root',
'ZZTo2L2Nu':'/store/mc/RunIIFall15MiniAODv2/ZZTo2L2Nu_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/029968AE-7AB9-E511-AB91-44A8423DE2C0.root',
'ZZTo2L2Q':'/store/mc/RunIIFall15MiniAODv2/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/004BED0E-B7C9-E511-AE1A-00145EDD7301.root',
'ZZTo4L':'/store/mc/RunIIFall15MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/00AE5257-01BA-E511-AA8F-002590596484.root',
'DYJetsToLL_M50':'/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/06532BBC-05C8-E511-A60A-F46D043B3CE5.root',
'TTTo2L2Nu':'/store/mc/RunIIFall15MiniAODv2/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/02A468DA-E8B9-E511-942C-0022195E688C.root',
'BulkGravToZZToZlepZinv_narrow_600':'/store/mc/RunIIFall15MiniAODv2/BulkGravToZZToZlepZinv_narrow_M-600_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000/3ABE89FB-F5D6-E511-A603-001E67E6920C.root',


}

LHAID_file = open('LHAPDFID_6.1.txt')
LHAID = {}

for line in LHAID_file.readlines():
    if not line.rstrip(): continue
    LHAID[line.split()[0]] = line.split()[1]

#print LHAID


for key in samples.keys():
    print '#################################'
    print ' LHE Info for dataset: '
    print '    '+key
    print '#################################'
#    print EDMWeightInfo.getWeightIDs(samples[key])

    edm_file_name = samples[key]

    if "/store/" in edm_file_name :
        #edm_file_name = "/".join(["root://cmsxrootd.fnal.gov/",
        edm_file_name = "/".join(["root://cms-xrd-global.cern.ch/",
            edm_file_name])
    elif not os.path.isfile(edm_file_name) :
        raise FileNotFoundException("File %s was not found." % edm_file_name)
    runs = Runs(edm_file_name)
    runInfoHandle = Handle("LHERunInfoProduct")
    run = runs.__iter__().next()
    run.getByLabel('externalLHEProducer', runInfoHandle)
    lheStuff = runInfoHandle.product()

    # id number of PDF used
    pdfidx = lheStuff.heprup().PDFSUP.first
    pdfname = LHAID[str(pdfidx)] if pdfidx!=-1 else ' undefined, please check the LHE header.'

    print '## Id number of PDF used: '+str(pdfidx)
    print '## PDF set name: '+pdfname
    print '####'

    lines = []
    for i, h in enumerate(lheStuff.headers_begin()) :
        if i == lheStuff.headers_size() :
            break
        hlines = []
        isWeights = False
        for line in h.lines() :
            hlines.append(line)
            if 'weightgroup' in line :
                isWeights = True
        if isWeights :
            lines.extend(hlines)
    #        break

    print ''.join(lines).rstrip("<")

    print '#################################'
