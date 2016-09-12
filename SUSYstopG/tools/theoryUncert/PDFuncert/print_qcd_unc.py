#!/usr/bin/env python 






both_qcd=  {'BulkGravToZZToZlepZinv_narrow_1400': 0.12241041020755528, 'ZZTo2L2Q': 0.17157885032755216, 'ZZTo2L2Nu': 0.03141690984403456, 'BulkGravToZZToZlepZinv_narrow_2500': 0.15616350132214624, 'BulkGravToZZToZlepZinv_narrow_600': 0.07372138846600179, 'WZTo1L1Nu2Q': 0.198487478233489, 'BulkGravToZZToZlepZinv_narrow_1000': 0.10344219928826914, 'TTTo2L2Nu': 0.3291266056332909, 'BulkGravToZZToZlepZinv_narrow_1800': 0.1365688640950583, 'BulkGravToZZToZlepZinv_narrow_1600': 0.12996139222356107, 'BulkGravToZZToZlepZinv_narrow_4000': 0.18939128275760397, 'BulkGravToZZToZlepZinv_narrow_800': 0.0906588342979816, 'BulkGravToZZToZlepZinv_narrow_3000': 0.16810681221646895, 'DYJetsToLL_M50_BIG': 0.15227216082174544, 'BulkGravToZZToZlepZinv_narrow_1200': 0.11373246952873717, 'WZTo2L2Q': 0.10789990823283352, 'ZZTo4L': 0.03851996217279202, 'WZTo3LNu': 0.07163078271584955, 'BulkGravToZZToZlepZinv_narrow_2000': 0.14262460074740602, 'WWTo2L2Nu': 0.07304338531167504, 'BulkGravToZZToZlepZinv_narrow_3500': 0.17918095785546106}

el_qcd=  {'BulkGravToZZToZlepZinv_narrow_1400': 0.12228375181504014, 'ZZTo2L2Q': 0.1771876951591218, 'ZZTo2L2Nu': 0.03101337748889854, 'BulkGravToZZToZlepZinv_narrow_2500': 0.15607966523395778, 'BulkGravToZZToZlepZinv_narrow_600': 0.07377717832681624, 'TTTo2L2Nu': 0.39575384856553414, 'BulkGravToZZToZlepZinv_narrow_1000': 0.10334910010684883, 'BulkGravToZZToZlepZinv_narrow_2000': 0.14260049526700014, 'BulkGravToZZToZlepZinv_narrow_1800': 0.13648573151830712, 'BulkGravToZZToZlepZinv_narrow_1600': 0.12984612071306645, 'BulkGravToZZToZlepZinv_narrow_4000': 0.18940935949646492, 'BulkGravToZZToZlepZinv_narrow_800': 0.09064542537011572, 'BulkGravToZZToZlepZinv_narrow_3000': 0.16811731899655274, 'DYJetsToLL_M50_BIG': 0.2103700961600778, 'BulkGravToZZToZlepZinv_narrow_1200': 0.11366201357823441, 'WZTo2L2Q': 0.11077151492917875, 'ZZTo4L': 0.033090632925140195, 'WZTo3LNu': 0.06520597118156746, 'WWTo2L2Nu': 0.07724615775045607, 'BulkGravToZZToZlepZinv_narrow_3500': 0.1791771475579833}

mu_qcd= {'BulkGravToZZToZlepZinv_narrow_1400': 0.12251740711150694, 'ZZTo2L2Q': 0.17774603416439355, 'ZZTo2L2Nu': 0.031646682235985546, 'BulkGravToZZToZlepZinv_narrow_2500': 0.156235984159667, 'BulkGravToZZToZlepZinv_narrow_600': 0.0736802130657011, 'WZTo1L1Nu2Q': 0.198487478233489, 'BulkGravToZZToZlepZinv_narrow_1000': 0.10351731394776553, 'TTTo2L2Nu': 0.31945082916610457, 'BulkGravToZZToZlepZinv_narrow_1800': 0.13663994186354655, 'BulkGravToZZToZlepZinv_narrow_1600': 0.13005856945482275, 'BulkGravToZZToZlepZinv_narrow_4000': 0.18937670074148955, 'BulkGravToZZToZlepZinv_narrow_800': 0.09066949633382265, 'BulkGravToZZToZlepZinv_narrow_3000': 0.16809785471885746, 'DYJetsToLL_M50_BIG': 0.12010918996649439, 'BulkGravToZZToZlepZinv_narrow_1200': 0.11379194437849888, 'WZTo2L2Q': 0.10657392759907497, 'ZZTo4L': 0.04187486431751758, 'WZTo3LNu': 0.0751827323542959, 'BulkGravToZZToZlepZinv_narrow_2000': 0.14264582950973947, 'WWTo2L2Nu': 0.07178155295825811, 'BulkGravToZZToZlepZinv_narrow_3500': 0.17918408020269477}

el_qcd.update({'WZTo1L1Nu2Q':0.0})

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


for bkg in bkg_list:
    print '{m} & {bpdf:.3f} & {epdf:.3f} & {mpdf:.3f} \\\\'.format(m=bkg,bpdf=both_qcd[bkg],epdf=el_qcd[bkg],mpdf=mu_qcd[bkg])


for m in masses:
    bkg='BulkGravToZZToZlepZinv_narrow_'+str(m)
    print '{m} & {bpdf:.3f} & {epdf:.3f} & {mpdf:.3f} \\\\'.format(m=m,bpdf=both_qcd[bkg],epdf=el_qcd[bkg],mpdf=mu_qcd[bkg])


