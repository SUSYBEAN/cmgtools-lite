#!/usr/bin/env python


both_pdf =  {'ZZTo2L2Q': 0.05471847950257513, 'ZZTo2L2Nu': 0.007637328337378304, 'WZTo1L1Nu2Q': 0.015811933150282276, 'TTTo2L2Nu': 0.015694194669118255, 'ZZTo4L': 0.020907014013051912, 'WZTo2L2Q': 0.07309134977731079, 'WZTo3LNu': 0.011325756743461571, 'WWTo2L2Nu': 0.008107971127675091,'BulkGravToZZToZlepZinv_narrow_1400': 0.004002804519443842, 'BulkGravToZZToZlepZinv_narrow_800': 0.006464697164413783, 'BulkGravToZZToZlepZinv_narrow_1000': 0.005307383499276305, 'BulkGravToZZToZlepZinv_narrow_1800': 0.0035957969659051307, 'BulkGravToZZToZlepZinv_narrow_600': 0.007143295578790101, 'BulkGravToZZToZlepZinv_narrow_3000': 0.0028518766926849152, 'BulkGravToZZToZlepZinv_narrow_4000': 0.004679912139220447, 'BulkGravToZZToZlepZinv_narrow_1600': 0.003107094984354692, 'BulkGravToZZToZlepZinv_narrow_1200': 0.005037572293931746, 'BulkGravToZZToZlepZinv_narrow_2000': 0.0038913125791678045, 'BulkGravToZZToZlepZinv_narrow_2500': 0.002519711699268357, 'BulkGravToZZToZlepZinv_narrow_3500': 0.0018048176852852834, 'DYJetsToLL_M50_BIG': 0.03061736872190807}

el_pdf =  {'ZZTo2L2Q': 0.06384851971523786, 'ZZTo2L2Nu': 0.008581577336354406, 'TTTo2L2Nu': 0.024263970179630657, 'ZZTo4L': 0.022357114693074075, 'WZTo2L2Q': 0.07050115350648561, 'WZTo3LNu': 0.01206019813110783, 'WWTo2L2Nu': 0.014244853624046205,'BulkGravToZZToZlepZinv_narrow_1400': 0.006013099263873355, 'BulkGravToZZToZlepZinv_narrow_800': 0.006233230177175719, 'BulkGravToZZToZlepZinv_narrow_1000': 0.0064994989023878435, 'BulkGravToZZToZlepZinv_narrow_1800': 0.005285926581221518, 'BulkGravToZZToZlepZinv_narrow_600': 0.0063658416418692024, 'BulkGravToZZToZlepZinv_narrow_3000': 0.0033379178791895273, 'BulkGravToZZToZlepZinv_narrow_4000': 0.0019478475742998218, 'BulkGravToZZToZlepZinv_narrow_1600': 0.005228883676442283, 'BulkGravToZZToZlepZinv_narrow_1200': 0.005873542566174066, 'BulkGravToZZToZlepZinv_narrow_2000': 0.004687713715847046, 'BulkGravToZZToZlepZinv_narrow_2500': 0.004804464556928468, 'BulkGravToZZToZlepZinv_narrow_3500': 0.001683766863406105, 'DYJetsToLL_M50_BIG': 0.03383722654588172}

mu_pdf =  {'ZZTo2L2Q': 0.07334424155824634, 'ZZTo2L2Nu': 0.007153241012682882, 'WZTo1L1Nu2Q': 0.015811933150282276, 'TTTo2L2Nu': 0.015135767594320913, 'ZZTo4L': 0.023559278009013876, 'WZTo2L2Q': 0.09330043292354906, 'WZTo3LNu': 0.011051750376789517, 'WWTo2L2Nu': 0.007298783555030545, 'BulkGravToZZToZlepZinv_narrow_1400': 0.002398045406088539, 'BulkGravToZZToZlepZinv_narrow_800': 0.006658479961393173, 'BulkGravToZZToZlepZinv_narrow_1000': 0.004350047864235286, 'BulkGravToZZToZlepZinv_narrow_1800': 0.0021997497078086576, 'BulkGravToZZToZlepZinv_narrow_600': 0.007721186129519206, 'BulkGravToZZToZlepZinv_narrow_3000': 0.0025390341100096406, 'BulkGravToZZToZlepZinv_narrow_4000': 0.008593956775935562, 'BulkGravToZZToZlepZinv_narrow_1600': 0.0013402973218664647, 'BulkGravToZZToZlepZinv_narrow_1200': 0.004350845850794779, 'BulkGravToZZToZlepZinv_narrow_2000': 0.0032408956547417133, 'BulkGravToZZToZlepZinv_narrow_2500': 0.0009060679735222341, 'BulkGravToZZToZlepZinv_narrow_3500': 0.0020928815185550444, 'DYJetsToLL_M50_BIG': 0.03298307584767903}


both_as =  {'ZZTo2L2Q': 0.007582453791299182, 'ZZTo2L2Nu': 0.005190205299191009, 'WZTo1L1Nu2Q': 0.003040920968416544, 'TTTo2L2Nu': 0.011217431482971163, 'ZZTo4L': 0.008855109470409295, 'WZTo2L2Q': 0.0013038648196066205, 'WZTo3LNu': 0.00860639664270249, 'WWTo2L2Nu': 0.00033665617569661466, 'DYJetsToLL_M50_BIG': 0.008484491447177778}

el_as =  {'ZZTo2L2Q': 0.028754021319167905, 'ZZTo2L2Nu': 0.005676042367939754, 'TTTo2L2Nu': 0.016407150012736488, 'ZZTo4L': 0.010202184108802026, 'WZTo2L2Q': 0.019197218700559876, 'WZTo3LNu': 0.008855107539317009, 'WWTo2L2Nu': 0.0012378849221417743, 'DYJetsToLL_M50_BIG': 0.005906984957651212}

mu_as =  {'ZZTo2L2Q': 0.018376538018274546, 'ZZTo2L2Nu': 0.004913568385289863, 'WZTo1L1Nu2Q': 0.003040920968416544, 'TTTo2L2Nu': 0.010463767452888573, 'ZZTo4L': 0.008022722376662517, 'WZTo2L2Q': 0.006958492539193617, 'WZTo3LNu': 0.00846889711067883, 'WWTo2L2Nu': 0.0008093933814286025, 'DYJetsToLL_M50_BIG': 0.010117610510392505}


both_pdf_as =  {'ZZTo2L2Q': 0.055241339634108454, 'ZZTo2L2Nu': 0.009234014034030492, 'WZTo1L1Nu2Q': 0.016101690292796934, 'TTTo2L2Nu': 0.019290892031923575, 'ZZTo4L': 0.022704981803007052, 'WZTo2L2Q': 0.07310297856952887, 'WZTo3LNu': 0.014224725972179077, 'WWTo2L2Nu': 0.00811495737437034, 'DYJetsToLL_M50_BIG': 0.03177121122290599}

el_pdf_as =  {'ZZTo2L2Q': 0.07002447580560872, 'ZZTo2L2Nu': 0.010288873920039973, 'TTTo2L2Nu': 0.029290524413510332, 'ZZTo4L': 0.024574888361682565, 'WZTo2L2Q': 0.07306809051550597, 'WZTo3LNu': 0.014961995471675753, 'WWTo2L2Nu': 0.014298538871191289, 'DYJetsToLL_M50_BIG': 0.03434895008012967}

mu_pdf_as =  {'ZZTo2L2Q': 0.07561134120812485, 'ZZTo2L2Nu': 0.008678249320136432, 'WZTo1L1Nu2Q': 0.016101690292796934, 'TTTo2L2Nu': 0.01840059482672844, 'ZZTo4L': 0.02488782141206832, 'WZTo2L2Q': 0.09355956072010861, 'WZTo3LNu': 0.013923483926881643, 'WWTo2L2Nu': 0.007343524973000687, 'DYJetsToLL_M50_BIG': 0.03449999035961826}


el_pdf.update({'WZTo1L1Nu2Q':0.0})
el_as.update({'WZTo1L1Nu2Q':0.0})
el_pdf_as.update({'WZTo1L1Nu2Q':0.0})

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



# background
print '\\begin{table}[htdp]'
print '\\caption{PDF uncertainties on acceptance for background.}'
print '\\begin{center}'
print '\\begin{footnotesize}'
print '\\begin{tabular}{c c c c}'
print '\\hline'
print ' MC sample & inc. PDF & inc. $\\alpha_s$  & inc. PDF+$\\alpha_s$ & el. PDF & el. $\\alpha_s$ & el. PDF+$\\alpha_s$ & mu. PDF & mu. $\\alpha_s$ & mu. PDF+$\\alpha_s$ \\\\'
print '\\hline'
for bkg in bkg_list:
    print '{bkg} & {bpdf:.3f} & {bas:.3f} & {bpdfas:.3f} & {epdf:.3f} & {eas:.3f} & {epdfas:.3f} & {mpdf:.3f} & {mas:.3f} & {mpdfas:.3f} \\\\'.format(bkg=bkg,bpdf=both_pdf[bkg],bas=both_as[bkg],bpdfas=both_pdf_as[bkg],epdf=el_pdf[bkg],eas=el_as[bkg],epdfas=el_pdf_as[bkg],mpdf=mu_pdf[bkg],mas=mu_as[bkg],mpdfas=mu_pdf_as[bkg])


print '\\hline'
print '\\end{tabular}'
print '\\end{footnotesize}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'


# signal
print '\\begin{table}[htdp]'
print '\\caption{PDF uncertainties on acceptance for signal.}'
print '\\begin{center}'
print '\\begin{footnotesize}'
print '\\begin{tabular}{c c c c}'
print '\\hline'
print ' Mass (GeV) &  inc. & el. & mu.   \\\\'
print '\\hline'
for m in masses:
    bkg='BulkGravToZZToZlepZinv_narrow_'+str(m)
    print '{m:n} & {bpdf:.3f} & {epdf:.3f} & {mpdf:.3f} \\\\'.format(m=m,bpdf=both_pdf[bkg],epdf=el_pdf[bkg],mpdf=mu_pdf[bkg])
print '\\hline'
print '\\end{tabular}'
print '\\end{footnotesize}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'





