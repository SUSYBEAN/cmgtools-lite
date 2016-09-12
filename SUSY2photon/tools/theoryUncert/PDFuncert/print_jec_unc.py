#!/usr/bin/env python 


both={'BulkGravToZZToZlepZinv_narrow_1400': 6.574431447858853e-05, 'ZZTo2L2Q': 0.04710649107317172, 'ZZTo2L2Nu': 0.004105121240172038, 'BulkGravToZZToZlepZinv_narrow_2500': 3.2142438173343724e-05, 'BulkGravToZZToZlepZinv_narrow_600': 0.0028029433476314813, 'WZTo1L1Nu2Q': 0.0, 'BulkGravToZZToZlepZinv_narrow_1000': 0.00038950878010390877, 'TTTo2L2Nu': 0.042511649722035305, 'BulkGravToZZToZlepZinv_narrow_1800': 4.616878752694564e-06, 'BulkGravToZZToZlepZinv_narrow_1600': 1.3069322067904565e-05, 'BulkGravToZZToZlepZinv_narrow_4000': 3.8490792946044294e-05, 'BulkGravToZZToZlepZinv_narrow_800': 0.0008754539673569983, 'BulkGravToZZToZlepZinv_narrow_3000': 5.897763867529493e-05, 'DYJetsToLL_M50_BIG': 0.04641909529342553, 'BulkGravToZZToZlepZinv_narrow_1200': 0.00013928784638400638, 'WZTo2L2Q': 0.004812611467655992, 'ZZTo4L': 0.0006725611844433654, 'WZTo3LNu': 0.01107334128223747, 'BulkGravToZZToZlepZinv_narrow_2000': 4.016196239570169e-05, 'WWTo2L2Nu': 0.03378082672156557, 'BulkGravToZZToZlepZinv_narrow_3500': 1.1925461753492517e-05}

el={'BulkGravToZZToZlepZinv_narrow_1400': 5.2669474322475374e-05, 'ZZTo2L2Q': 0.04046269303105626, 'ZZTo2L2Nu': 0.0013832113840990257, 'BulkGravToZZToZlepZinv_narrow_2500': 0.0, 'BulkGravToZZToZlepZinv_narrow_600': 0.0028319207934106583, 'TTTo2L2Nu': 0.023062257455513546, 'BulkGravToZZToZlepZinv_narrow_1000': 0.000580822915617285, 'BulkGravToZZToZlepZinv_narrow_2000': 0.0, 'BulkGravToZZToZlepZinv_narrow_1800': 0.0, 'BulkGravToZZToZlepZinv_narrow_1600': 2.6916120993147885e-07, 'BulkGravToZZToZlepZinv_narrow_4000': 0.0, 'BulkGravToZZToZlepZinv_narrow_800': 0.0009125520903262574, 'BulkGravToZZToZlepZinv_narrow_3000': 0.0, 'DYJetsToLL_M50_BIG': 0.0, 'BulkGravToZZToZlepZinv_narrow_1200': 0.0002388542472250843, 'WZTo2L2Q': 0.0008325364319831108, 'ZZTo4L': 0.004528765310877203, 'WZTo3LNu': 0.007242858015256537, 'WWTo2L2Nu': 0.07150152595130499, 'BulkGravToZZToZlepZinv_narrow_3500': 0.0}

mu={'BulkGravToZZToZlepZinv_narrow_1400': 7.678951557849034e-05, 'ZZTo2L2Q': 0.049080092404494846, 'ZZTo2L2Nu': 0.005654983922408152, 'BulkGravToZZToZlepZinv_narrow_2500': 5.9932083109226486e-05, 'BulkGravToZZToZlepZinv_narrow_600': 0.002781556702955379, 'WZTo1L1Nu2Q': 0.0, 'BulkGravToZZToZlepZinv_narrow_1000': 0.00023515193261691714, 'TTTo2L2Nu': 0.04533613958100913, 'BulkGravToZZToZlepZinv_narrow_1800': 8.56427761758427e-06, 'BulkGravToZZToZlepZinv_narrow_1600': 2.3860229435568314e-05, 'BulkGravToZZToZlepZinv_narrow_4000': 6.954027412137131e-05, 'BulkGravToZZToZlepZinv_narrow_800': 0.0008459555972725741, 'BulkGravToZZToZlepZinv_narrow_3000': 0.00010925870035188856, 'DYJetsToLL_M50_BIG': 0.07583043129933081, 'BulkGravToZZToZlepZinv_narrow_1200': 5.5239633428960655e-05, 'WZTo2L2Q': 0.006650433626562591, 'ZZTo4L': 0.0038865754058347024, 'WZTo3LNu': 0.013191019533596526, 'BulkGravToZZToZlepZinv_narrow_2000': 7.553104611396355e-05, 'WWTo2L2Nu': 0.022455636177723837, 'BulkGravToZZToZlepZinv_narrow_3500': 2.1697778146367863e-05}


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


for bkg in bkg_list:
    print '{m} & {bpdf:.3f} & {epdf:.3f} & {mpdf:.3f} \\\\'.format(m=bkg,bpdf=both[bkg],epdf=el[bkg],mpdf=mu[bkg])


for m in masses:
    bkg='BulkGravToZZToZlepZinv_narrow_'+str(m)
    print '{m} & {bpdf:.3f} & {epdf:.3f} & {mpdf:.3f} \\\\'.format(m=m,bpdf=both[bkg],epdf=el[bkg],mpdf=mu[bkg])


