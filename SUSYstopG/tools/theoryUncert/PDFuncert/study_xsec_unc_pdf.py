#!/usr/bin/env python

import ROOT 
import pprint,pickle
import math

indir='LHEWeights'

#bkg_list=['DYJetsToLL_M50_BIG']
bkg_list=['DYJetsToLL_M50_BIG','ZZTo4L', 'ZZTo2L2Q', 'ZZTo2L2Nu', 'TTTo2L2Nu', 'WZTo2L2Q', 'WZTo3LNu', 'WZTo1L1Nu2Q', 'WWTo2L2Nu']

sig_list=['BulkGravToZZToZlepZinv_narrow_1000',
'BulkGravToZZToZlepZinv_narrow_1200',
'BulkGravToZZToZlepZinv_narrow_1400',
'BulkGravToZZToZlepZinv_narrow_1600',
'BulkGravToZZToZlepZinv_narrow_1800',
'BulkGravToZZToZlepZinv_narrow_2000',
'BulkGravToZZToZlepZinv_narrow_2500',
'BulkGravToZZToZlepZinv_narrow_3000',
'BulkGravToZZToZlepZinv_narrow_3500',
'BulkGravToZZToZlepZinv_narrow_4000',
'BulkGravToZZToZlepZinv_narrow_600',
'BulkGravToZZToZlepZinv_narrow_800']

samples = bkg_list+sig_list


#masses = []

masses = [
         600,
         800,
          1000,
         1200,
          1400,
         1600,1800,
          2000,
         2500,
          3000,
         3500,
         4000,]


#file= open('DYJetsToLL_M50_LHEWeightReport.pck')

#data = pickle.load(file)

#pprint.pprint(data)

# pdf weights id for bkg: 2001-2100,  as up/down 2101,2102
wtIds_pdf_bkg = range(2001,2101)
wtIds_as_bkg = [2101,2102]

# pdf weights id for sig: 111-210, no as up/down
wtIds_pdf_sig = range(111,211)

# qcd weights id
wtIds_qcd_bkg = [1002,1003,1004,1005,1007,1009]
wtIds_qcd_sig = [2,3,4,5,7,9] 

sample='DYJetsToLL_M50_BIG'


fout = ROOT.TFile("study_xsec_unc_pdf.root", "recreate")



hunc_pdf = {}
all_unc_pdf = {}
all_unc_as = {}
all_unc_pdf_as = {}

hunc_qcd = {}
all_var_qcd = {}
all_unc_qcd = {}

#for sample in [sample]:
for sample in samples:

    # open pickle file
    file= open(indir+'/'+sample+'/LHEWeightAnalyzer/LHEWeightReport.pck')

    data = pickle.load(file)

    print 'PDF for '+sample
    wtIds_pdf = wtIds_pdf_bkg
    if 'BulkGrav' in sample:
        wtIds_pdf = wtIds_pdf_sig

    hunc_pdf[sample] = ROOT.TH1F(sample+'_hunc_pdf', sample+'_hunc_pdf', 1000,0,2)

    for idx in wtIds_pdf:
        print ' - PDF '+str(idx)
        # this is n_i/n_0
        n = float(data[str(idx)][1])/float(data['SumLHEOrigWeights'][1])
        print '  ni/n0 = '+str(n)
        hunc_pdf[sample].Fill(n)

    pdf_unc = hunc_pdf[sample].GetRMS()
    print '  pdf_unc = '+str(pdf_unc)
    hunc_pdf[sample].Write()

    all_unc_pdf[sample] = pdf_unc
    #
    pdf_as_unc = pdf_unc

    if not ('BulkGrav' in sample):
        print ' - Alpha_s '
        n1 = float(data[str(wtIds_as_bkg[0])][1])/float(data['SumLHEOrigWeights'][1])
        print '  n1/n0 = '+str(n1)
        n2 = float(data[str(wtIds_as_bkg[1])][1])/float(data['SumLHEOrigWeights'][1])
        print '  n2/n0 = '+str(n2)
        as_unc = abs(n1-n2)/2.0*1.5  # scale by 1.5
        print '  as_unc = '+str(as_unc)
        pdf_as_unc = math.sqrt(pdf_unc**2+as_unc**2)
        all_unc_as[sample] = as_unc

    all_unc_pdf_as[sample] = pdf_as_unc


    print ' pdf_as_unc = '+str(pdf_as_unc)

    # QCD
    print 'QCD for '+sample

    wtIds_qcd = wtIds_qcd_bkg
    if 'BulkGrav' in sample:
        wtIds_qcd = wtIds_qcd_sig

    hunc_qcd[sample] = ROOT.TH1F(sample+'_hunc_qcd', sample+'_hunc_qcd', 1000,0,2)

    all_var_qcd[sample] = []
    for idx in wtIds_qcd:
        print ' - QCD '+str(idx)
        # n is actually ni/n0
        n = float(data[str(idx)][1])/float(data['SumLHEOrigWeights'][1])
        print '  ni/n0 = '+str(n)
        hunc_qcd[sample].Fill(n)
        all_var_qcd[sample].append(n)

    qcd_up = max(all_var_qcd[sample])
    qcd_dn = min(all_var_qcd[sample])
    qcd_unc = abs(qcd_up-qcd_dn)/2.0
    print '  qcd_unc = '+str(qcd_unc)
    hunc_qcd[sample].Write()

    all_unc_qcd[sample] = qcd_unc

print '#######################'
print 'All PDF uncertainties: '
print 'all_unc_pdf =',all_unc_pdf
print 'all_unc_as = ',all_unc_as
print 'all_unc_pdf_as = ',all_unc_pdf_as
print '#######################'
print 'All QCD uncertainties: '
print 'all_unc_qcd = ',all_unc_qcd
print 'all_var_qcd = ',all_var_qcd
print '#######################'

fout.Close()

   

# print to tables

# PDF
# background
print '\\begin{table}[htdp]'
print '\\caption{PDF uncertainties on cross-sections for background.}'
print '\\begin{center}'
print '\\begin{footnotesize}'
print '\\begin{tabular}{c c c c }'
print '\\hline'
print ' MC sample & PDF & $\\alpha_{s}$  &  PDF+$\\alpha_{s}$  \\\\'
print '\\hline'
for bkg in bkg_list:
    print '{bkg} & {bpdf:.3f} & {bas:.3f} & {bpdfas:.3f} \\\\'.format(bkg=bkg,bpdf=all_unc_pdf[bkg],bas=all_unc_as[bkg],bpdfas=all_unc_pdf_as[bkg])
print '\\hline'
print '\\end{tabular}'
print '\\end{footnotesize}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'

# signal
print '\\begin{table}[htdp]'
print '\\caption{PDF uncertainties on cross-sections for signal.}'
print '\\begin{center}'
print '\\begin{footnotesize}'
print '\\begin{tabular}{c c }'
print '\\hline'
print ' Mass (GeV) & PDF unc.  \\\\'
print '\\hline'
for m in masses:
    bkg='BulkGravToZZToZlepZinv_narrow_'+str(m)
    print '{m:n} & {bpdf:.3f} \\\\'.format(m=m,bpdf=all_unc_pdf[bkg])
print '\\hline'
print '\\end{tabular}'
print '\\end{footnotesize}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'


# QCD
# background
print '\\begin{table}[htdp]'
print '\\caption{QCD uncertainties on cross-sections for background.}'
print '\\begin{center}'
print '\\begin{footnotesize}'
print '\\begin{tabular}{c c }'
print '\\hline'
print ' MC sample & Uncertainty  \\\\'
print '\\hline'
for bkg in bkg_list:
    print '{bkg} & {bpdf:.3f}  \\\\'.format(bkg=bkg,bpdf=all_unc_qcd[bkg])
print '\\hline'
print '\\end{tabular}'
print '\\end{footnotesize}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'

# signal
print '\\begin{table}[htdp]'
print '\\caption{QCD uncertainties on cross-sections for signal.}'
print '\\begin{center}'
print '\\begin{footnotesize}'
print '\\begin{tabular}{c c }'
print '\\hline'
print ' Mass (GeV) & Uncertainty  \\\\'
print '\\hline'
for m in masses:
    bkg='BulkGravToZZToZlepZinv_narrow_'+str(m)
    print '{m:n} & {bpdf:.3f} \\\\'.format(m=m,bpdf=all_unc_qcd[bkg])
print '\\hline'
print '\\end{tabular}'
print '\\end{footnotesize}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'


 
