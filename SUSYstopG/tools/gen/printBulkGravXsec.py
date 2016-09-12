#!/bin/env python


def latex_float(f):
    float_str = "{0:.3e}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))
    else:
        return float_str


#total xsec:
# https://github.com/acarvalh/Cross_sections_CMS/blob/master/WED/bulk_KKgrav_LHC13.txt
#if you want the cross section (CX) using a ktilda' you do
#CX(ktilda',pb) = CX(ktilda=0.1,pb) * (ktilda' / 0.1)^2

ktilda=0.5

masses = [600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500]
BulkGXsec = {200:5.649417401247585,
260:1.9429371431011526,
300:1.0684873854214705,
350:0.5514532256075042,
400:0.30646520928197674,
450:0.1796483415002634,
500:0.1102129890727813,
550:0.07045280534342299,
600:0.046310769831789535,
650:0.03127154730329527,
700:0.02161133022878741,
750:0.01524141818968972,
800:0.010953371285280107,
850:0.007973233764684683,
900:0.005886101099004809,
950:0.004386360642407879,
1000:0.0033153186581553304,
1050:0.0025337224464699252,
1100:0.0019499014923576293,
1150:0.0015106660056920108,
1200:0.0011796819730577413,
1250:0.0009290485361395383,
1300:0.0007342676766258561,
1350:0.0005857353481100774,
1400:0.0004688386674106958,
1450:0.0003780591771340791,
1500:0.00030704401954085035,
1600:2.03E-04,
1700:1.37E-04,
1800:9.36E-05,
1900:6.49E-05,
2000:0.00004542880185027153,
2100:3.21E-05,
2200:2.28E-05,
2300:1.64E-05,
2400:1.19E-05,
2500:8.674571731526816e-6,
2600:6.34E-06,
2700:4.67E-06,
2800:3.46E-06,
2900:2.58E-06,
3000:1.9201182127764607e-6,
3100:1.44E-06,
3200:1.09E-06,
3300:8.19E-07,
3400:6.23E-07,
3500:4.74E-07,
3600:3.60E-07,
3700:2.74E-07,
3800:2.09E-07,
3900:1.61E-07,
4000:1.23E-07,
4100:9.45E-08,
4200:7.27E-08,
4300:5.59E-08,
4400:4.33E-08,
4500:3.34E-08,
}

# apply ktilda scale
#CX(ktilda',pb) = CX(ktilda=0.1,pb) * (ktilda' / 0.1)^2
for mass,value in BulkGXsec.iteritems():
    BulkGXsec[mass] = BulkGXsec[mass]*(ktilda/0.1)**2

# BulkGrav to ZZ branch ratio:
# https://github.com/acarvalh/Cross_sections_CMS/blob/master/WED/bulk_KKgrav_decay.txt
ZZBr = {200:0.316483,
300:0.355847,
400:0.308107,
500:0.230339,
600:0.184219,
700:0.158352,
800:0.142802,
900:0.132777,
1000:0.125936,
1100:0.121055,
1200:0.117445,
1300:0.114697,
1400:0.112554,
1500:0.110851,
1600:0.109472,
1700:0.108341,
1800:0.107401,
1900:0.106611,
2000:0.10594,
2100:0.105366,
2200:0.104871,
2300:0.10444,
2400:0.104063,
2500:0.103732,
2600:0.103439,
2700:0.103178,
2800:0.102945,
2900:0.102737,
3000:0.102549,
3100:0.102379,
3200:0.102225,
3300:0.102086,
3400:0.101958,
3500:0.101841,
3600:0.101734,
3700:0.101636,
3800:0.101546,
3900:0.101462,
4000:0.101384}

BrZll = (3.363+3.366+3.370)*0.01 # PDG value
BrZinv = 0.200 #PDG value

#calculate BulkG->ZZ->2l2nu Xsec
BulkGZZ2l2nuXsec = {}
for mass in masses:
    if mass in BulkGXsec.keys() and mass in ZZBr.keys():
        BulkGZZ2l2nuXsec[mass] = BulkGXsec[mass]*ZZBr[mass]*BrZll*BrZinv*2.0


#calculate BulkG->ZZ Xsec
BulkGZZXsec = {}
for mass in masses:
    if mass in BulkGXsec.keys() and mass in ZZBr.keys():
        BulkGZZXsec[mass] = BulkGXsec[mass]*ZZBr[mass]

# print BulkG->ZZ->2l2nu Xsec
print "BulkG->ZZ->2l2nu Xsec:"
for mass in masses:
    if mass in BulkGZZ2l2nuXsec.keys():
        print mass,"{:.5e}".format(BulkGZZ2l2nuXsec[mass])
    else:
        print mass,"--"
print 

# print BulkG->ZZ->2l2nu Xsec dictionary codes
print "BulkG->ZZ->2l2nu Xsec Dictionary codes:"
print "BulkGZZ2l2nuXsec = {"
for mass in masses:
    if mass in BulkGZZ2l2nuXsec.keys():
        print str(mass)+":{:.5e}".format(BulkGZZ2l2nuXsec[mass])+","
print "}"
print

# another version of print codes
print "BulkG->ZZ->2l2nu Xsec Dictionary codes, another version:"
print "BulkGZZ2l2nuXsec = {"
for mass in masses:
    if mass in BulkGZZ2l2nuXsec.keys():
        print '\'m'+str(mass)+'\''+":{:.5e}".format(BulkGZZ2l2nuXsec[mass])+","
print "}"
print

#print latex table
print "BulkG->ZZ->2l2nu Xsec Latex table:"
print '\\begin{table}[htdp]'
print '\\caption{BulkGrav$\\rightarrow$ZZ$\\rightarrow ll\\nu\\nu$ cross-sections with $\\tilde{k}='+str(ktilda)+'$.}'
print '\\begin{center}'
print '\\begin{tabular}{c c}'
print '\\hline'
print 'mass points (GeV) & cross-section (pb) \\\\'
print '\\hline'

for mass in masses:
    if mass in BulkGZZ2l2nuXsec.keys() :
        print str(mass)+" & "+latex_float(BulkGZZ2l2nuXsec[mass])+" \\\\"
    else:
        print str(mass)+" & "+" --- "+" \\\\"

print '\\hline'
print '\\end{tabular}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'
print 

# print BulkG->ZZ Xsec
print "BulkG->ZZ Xsec:"
for mass in masses:
    if mass in BulkGZZXsec.keys():
        print mass,"{:.5e}".format(BulkGZZXsec[mass])
    else:
        print mass,"--"
print

# print BulkG->ZZ Xsec dictionary codes
print "BulkG->ZZ Xsec Dictionary codes:"
print "BulkGZZXsec = {"
for mass in masses:
    if mass in BulkGZZXsec.keys():
        print str(mass)+":{:.5e}".format(BulkGZZXsec[mass])+","
print "}"
print

# another version, print BulkG->ZZ Xsec dictionary codes
print "BulkG->ZZ Xsec Dictionary codes, another version:"
print "BulkGZZXsec = {"
for mass in masses:
    if mass in BulkGZZXsec.keys():
        print '\'m'+str(mass)+'\''+":{:.5e}".format(BulkGZZXsec[mass])+","
print "}"
print


# print BulkG->ZZ Xsec Latex table:
print "BulkG->ZZ Xsec Latex table:"
print '\\begin{table}[htdp]'
print '\\caption{BulkGrav$\\rightarrow$ZZ cross-sections with $\\tilde{k}='+str(ktilda)+'$.}'
print '\\begin{center}'
print '\\begin{tabular}{c c}'
print '\\hline'
print 'mass points (GeV) & cross-section (pb) \\\\'
print '\\hline'

for mass in masses:
    if mass in BulkGZZXsec.keys() :
        print str(mass)+" & "+latex_float(BulkGZZXsec[mass])+" \\\\"
    else:
        print str(mass)+" & "+" --- "+" \\\\"

print '\\hline'
print '\\end{tabular}'
print '\\end{center}'
print '\\label{default}'
print '\\end{table}'


