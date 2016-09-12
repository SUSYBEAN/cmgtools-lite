#!/usr/bin/env python


both = {'ZZTo2L2Q': 0.17645240703300102, 'ZZTo2L2Nu': 0.029099583374124938, 'WZTo1L1Nu2Q': 0.16083390677305714, 'TTTo2L2Nu': 0.1881922369141552, 'ZZTo4L': 0.0828129581759191, 'WZTo2L2Q': 0.09708491526598456, 'WZTo3LNu': 0.07734907737311292, 'WWTo2L2Nu': 0.06378467333531118, 'BulkGravToZZToZlepZinv_narrow_1400': 0.0002000872144263255, 'BulkGravToZZToZlepZinv_narrow_800': 0.0005639329969880369, 'BulkGravToZZToZlepZinv_narrow_1000': 0.00037762536808189084, 'BulkGravToZZToZlepZinv_narrow_1800': 0.00014161995638645175, 'BulkGravToZZToZlepZinv_narrow_600': 0.0008258787621844843, 'BulkGravToZZToZlepZinv_narrow_3000': 3.4652330997031466e-05, 'BulkGravToZZToZlepZinv_narrow_4000': 2.9022145294099744e-05, 'BulkGravToZZToZlepZinv_narrow_1600': 0.00015593582501316483, 'BulkGravToZZToZlepZinv_narrow_1200': 0.0002923457488859804, 'BulkGravToZZToZlepZinv_narrow_2000': 0.0001411442281681885, 'BulkGravToZZToZlepZinv_narrow_2500': 6.359653884641103e-05, 'BulkGravToZZToZlepZinv_narrow_3500': 9.552707982196651e-06,'DYJetsToLL_M50_BIG': 0.16835596126929742}

el = {'ZZTo2L2Q': 0.19494506285226298, 'ZZTo2L2Nu': 0.028697211211152707, 'TTTo2L2Nu': 0.2529110636664725, 'ZZTo4L': 0.07897210620459105, 'WZTo2L2Q': 0.0972134053458823, 'WZTo3LNu': 0.07092410153249873, 'WWTo2L2Nu': 0.06796417021790524, 'BulkGravToZZToZlepZinv_narrow_1400': 0.00032368309482416757, 'BulkGravToZZToZlepZinv_narrow_800': 0.000577405900962058, 'BulkGravToZZToZlepZinv_narrow_1000': 0.00046901404265786706, 'BulkGravToZZToZlepZinv_narrow_1800': 0.00022242336489330938, 'BulkGravToZZToZlepZinv_narrow_600': 0.0007708363078096081, 'BulkGravToZZToZlepZinv_narrow_3000': 2.447441547931417e-05, 'BulkGravToZZToZlepZinv_narrow_4000': 4.4695888886614465e-05, 'BulkGravToZZToZlepZinv_narrow_1600': 0.00026821143532501646, 'BulkGravToZZToZlepZinv_narrow_1200': 0.00036149277951325054, 'BulkGravToZZToZlepZinv_narrow_2000': 0.00016444283780803204, 'BulkGravToZZToZlepZinv_narrow_2500': 0.00014447808547696228, 'BulkGravToZZToZlepZinv_narrow_3500': 1.282535273355867e-05, 'DYJetsToLL_M50_BIG': 0.22678024307397066}

mu = {'ZZTo2L2Q': 0.1747936023652248, 'ZZTo2L2Nu': 0.029328695149584583, 'WZTo1L1Nu2Q': 0.16083390677305714, 'TTTo2L2Nu': 0.1787936054637918, 'ZZTo4L': 0.08518630521906595, 'WZTo2L2Q': 0.10376126622292886, 'WZTo3LNu': 0.08090111784812537, 'WWTo2L2Nu': 0.06252982919033717, 'BulkGravToZZToZlepZinv_narrow_1400': 9.567742143024027e-05, 'BulkGravToZZToZlepZinv_narrow_800': 0.0005532200908990315, 'BulkGravToZZToZlepZinv_narrow_1000': 0.00030389078665504776, 'BulkGravToZZToZlepZinv_narrow_1800': 7.253361036724248e-05, 'BulkGravToZZToZlepZinv_narrow_600': 0.0008665025432698092, 'BulkGravToZZToZlepZinv_narrow_3000': 4.332945692869794e-05, 'BulkGravToZZToZlepZinv_narrow_4000': 1.6378560443064405e-05, 'BulkGravToZZToZlepZinv_narrow_1600': 6.128422461987304e-05, 'BulkGravToZZToZlepZinv_narrow_1200': 0.00023397581377881949, 'BulkGravToZZToZlepZinv_narrow_2000': 0.00012062604564688773, 'BulkGravToZZToZlepZinv_narrow_2500': 6.3318670388556875e-06, 'BulkGravToZZToZlepZinv_narrow_3500': 6.870940128911318e-06, 'DYJetsToLL_M50_BIG': 0.13133808808157216}



el.update({'WZTo1L1Nu2Q':0.0})

bkg_list=['DYJetsToLL_M50_BIG','ZZTo4L', 'ZZTo2L2Q', 'ZZTo2L2Nu', 'TTTo2L2Nu', 'WZTo2L2Q', 'WZTo3LNu', 'WZTo1L1Nu2Q', 'WWTo2L2Nu']
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

# QCD
# background
print '\\begin{table}[htdp]'
print '\\caption{QCD uncertainties on acceptance for background.}'
print '\\begin{center}'
print '\\begin{footnotesize}'
print '\\begin{tabular}{c c c c}'
print '\\hline'
print ' MC sample & inc. & el. & mu.  \\\\'
print '\\hline'
for bkg in bkg_list:
    print '{m} & {bpdf:.3f} & {epdf:.3f} & {mpdf:.3f} \\\\'.format(m=bkg,bpdf=both[bkg],epdf=el[bkg],mpdf=mu[bkg])
print '\\hline'
print '\\end{tabular}'
print '\\end{footnotesize}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'


# signal
print '\\begin{table}[htdp]'
print '\\caption{QCD uncertainties on acceptance for signal.}'
print '\\begin{center}'
print '\\begin{footnotesize}'
print '\\begin{tabular}{c c c c}'
print '\\hline'
print ' Mass (GeV) &  inc. & el. & mu.   \\\\'
print '\\hline'
for m in masses:
    bkg='BulkGravToZZToZlepZinv_narrow_'+str(m)
    print '{m} & {bpdf:.3f} & {epdf:.3f} & {mpdf:.3f} \\\\'.format(m=m,bpdf=both[bkg],epdf=el[bkg],mpdf=mu[bkg])
print '\\hline'
print '\\end{tabular}'
print '\\end{footnotesize}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'





