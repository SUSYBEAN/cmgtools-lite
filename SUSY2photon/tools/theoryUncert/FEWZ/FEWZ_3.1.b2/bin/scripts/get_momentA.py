"""
Handle normalization of Collin-Soper angles and moments
Used by "finish.sh"
Usage: python $0 <outputfile>
   or  python $0 <outputfile> <newoutputfile> 
Note: <outputfile> should be the full file name and path of the file can be added in front
      <newoutputfile> is set to be the file name of <outputfile> by default (without any path in front) if not given
      <outputfile> will be overwritten if <outputfile> is in the current directory and <newoutputfile> is not present 
Fact: A common-separated-values (CSV) file named <newoutputfile>.csv will be output for histograms if 'GenerateCSV' is set to be true
      'HistNamesXxx' gives a list of key words appear in the histogram titles in the output file
      'NewHistNamesXX' gives the titles for the normalized histograms for the new output file
      Making all the lists mentioned above empty will generate an output file identical to the input file
"""

import sys
import math
import os

### whether to output part of the file which does not require normalization
KeeptheRest = True

### whether to generate CSV file for histograms
GenerateCSV = False
#GenerateCSV = True

### Data Type 1 w/ histogram bin by bin division
HistNamesNum=['A_FB numer','A_0 numer','A_1 numer','A_2 numer','A_3 numer','A_4 numer']
HistNamesDen=['A_FB denom','A_i denom','A_i denom','A_i denom','A_i denom','A_i denom']
NewHistNamesND=['A_FB vs Qll','A0 vs ZpT','A1 vs ZpT','A2 vs ZpT','A3 vs ZpT','A4 vs ZpT']

### Data Type 2 w/ histogram bin by overall factor division
HistNamesRaw=['----   phi CS         ----', \
              '----   costh CS       ----']
HistNamesFac=['----   phi CS norm    ----', \
              '----   cth CS norm    ----',]
NewHistNamesRF=['phi CS','costh CS']

# Trick:
#   making all the lists empty and set 'GenerateCSV' to be true will simply leave the output file intact and generate a CSV file for histograms
#   i.e.
#      GenerateCSV = True
#      HistNamesNum = []
#      HistNamesDen = []
#      NewHistNamesND = []
#      HistNamesRaw = []
#      HistNamesFac = []
#      NewHistNamesRF = []

### import necessary functions and possibly overide above variable definitions
from defs import *

### Define the division is operated with the error calculated at the same time
#def myisnan(fltnum):
#    return str(fltnum) == 'nan'
#
#def myiszero(fltnum):
#    return fltnum == 0 and not myisnan(fltnum)
#
#def myfloat(numstr):
#    try:
#        return float(numstr)
#    except ValueError:
#        if numstr == 'undefined':
#            return 'undefined'
#        return float('nan')
#
#def normdivision(num=0,den=0,numerr=0,denerr=0):
#    try:
#        return num/den, math.sqrt((numerr/den)**2+(denerr*num/den**2)**2)
#    except ZeroDivisionError:
#        if myiszero(num):
#            return 0, 'undefined'
#        return 'undefined', 'undefined'
### Imported form defs.py already

try:
    if len(sys.argv) < 2:
        raise Exception
    else:
        infile = sys.argv[1]
        infileloclist = pathseparator(infile)
        tmpoutfile =  infileloclist[1]+'_tmp'
        while os.path.isfile(tmpoutfile):
            tmpoutfile = tmpoutfile + '_'
        if len(sys.argv) == 2:
            outfile = infileloclist[1]
        else:
            outfile = sys.argv[2]
except Exception:
    print('Missing arguments.')
    raise

try:
    fin = open(infile, 'r')
except IOError:
    print("Can't open or find file.")
    raise

try:
    lenHn = len(HistNamesNum)
    lenHr = len(HistNamesRaw)
    if lenHn != len(HistNamesNum) or lenHr != len(HistNamesFac):
        raise Exception
except Exception:
    print('Error initializing histogram names.')
    raise

try:
    fout = open(tmpoutfile, 'w')
except IOError:
    print("Can't open temporary file for writing.")
    raise

if GenerateCSV:
    CSVfile = outfile + '.csv'
    try:
        fCSV = open(CSVfile, 'w')
    except IOError:
        print("Can't open csv file for writing.")
        raise

try:
    HistNumValues = [[] for i in range(lenHn)]
    HistDenValues = [[] for i in range(lenHn)]
    HistRawValues = [[] for i in range(lenHr)]
    HistFacValues = [[] for i in range(lenHr)]
    nextline = fin.readline()
    spacebuffer = ''
    HistMatchFlag = False
    tmpLineName = '----'
    tmpHistValues = []
    while len(nextline) > 0: # length = 0 if end of file reached
        tmpLineName = nextline
        for HNid in range(lenHn):
            if tmpLineName.find(HistNamesNum[HNid]) != -1:
                if not HistMatchFlag:
                    HistMatchFlag = True
                    nextline = fin.readline()
                    data = nextline.split()
                    while len(nextline) > 0 and len(data) != 3:
                        nextline = fin.readline()
                        data = nextline.split()
                    while len(nextline) > 0 and len(data) == 3:
                        tmpHistValues.append([data[0]]+map(myfloat,data[1:]))
                        nextline = fin.readline()
                        data = nextline.split()
                HistNumValues[HNid] = tmpHistValues
            elif tmpLineName.find(HistNamesDen[HNid]) != -1:
                if not HistMatchFlag:
                    HistMatchFlag = True
                    nextline = fin.readline()
                    data = nextline.split()
                    while len(nextline) > 0 and len(data) != 3:
                        nextline = fin.readline()
                        data = nextline.split()
                    while len(nextline) > 0 and len(data) == 3:
                        tmpHistValues.append([data[0]]+map(myfloat,data[1:]))
                        nextline = fin.readline()
                        data = nextline.split()
                HistDenValues[HNid] = tmpHistValues
        for HNid in range(lenHr):
            if tmpLineName.find(HistNamesRaw[HNid]) != -1:
                if not HistMatchFlag:
                    HistMatchFlag = True
                    nextline = fin.readline()
                    data = nextline.split()
                    while len(nextline) > 0 and len(data) != 3:
                        nextline = fin.readline()
                        data = nextline.split()
                    while len(nextline) > 0 and len(data) == 3:
                        tmpHistValues.append([data[0]]+map(myfloat,data[1:]))
                        nextline = fin.readline()
                        data = nextline.split()
                HistRawValues[HNid] = tmpHistValues
            elif tmpLineName.find(HistNamesFac[HNid]) != -1:
                if not HistMatchFlag:
                    HistMatchFlag = True
                    nextline = fin.readline()
                    data = nextline.split()
                    while len(nextline) > 0 and len(data) != 3:
                        nextline = fin.readline()
                        data = nextline.split()
                    while len(nextline) > 0 and len(data) == 3:
                        tmpHistValues.append([data[0]]+map(myfloat,data[1:]))
                        nextline = fin.readline()
                        data = nextline.split()
                HistFacValues[HNid] = tmpHistValues
        if HistMatchFlag:
            HistMatchFlag = False
            tmpHistValues = []  
            while nextline.isspace():
                nextline = fin.readline()
        elif KeeptheRest:
            if nextline.isspace():
                spacebuffer = spacebuffer + nextline
            else:
                fout.write(spacebuffer + nextline)
                spacebuffer = ''
            if GenerateCSV:
                if nextline.find('----') != -1:
                    fCSV.write(nextline.replace('----','').strip()+'\n')
                else:
                    data = nextline.split()
                    try:
                        ### Must at least one number and all be number to be outputed to CSV
                        if len(map(float,data[1:])) > 0:
                            fCSV.write(data[0])
                            if len(data) > 1:
                                for n in range(len(data)-1):
                                    fCSV.write(',' + data[n+1])
                            fCSV.write('\n')
                    except ValueError:
                        pass
            nextline = fin.readline()
        else:
            nextline = fin.readline()
    fin.close()
except IOError:
    print('Error reading/writing file')
    raise

### Check if histograms are empty and remove it
while [] in HistNumValues:
    #print('Error specifying the histogram names of data type 1.')
    EmptyCellIndex = HistNumValues.index([])
    HistNumValues.pop(EmptyCellIndex)
    HistDenValues.pop(EmptyCellIndex)
    NewHistNamesND.pop(EmptyCellIndex)
while [] in HistDenValues:
    #print('Error specifying the histogram names of data type 1.')
    EmptyCellIndex = HistDenValues.index([])
    HistNumValues.pop(EmptyCellIndex)
    HistDenValues.pop(EmptyCellIndex)
    NewHistNamesND.pop(EmptyCellIndex)
while [] in HistRawValues:
    #print('Error specifying the histogram names of data type 2.')
    EmptyCellIndex = HistRawValues.index([])
    HistRawValues.pop(EmptyCellIndex)
    HistFacValues.pop(EmptyCellIndex)
    NewHistNamesRF.pop(EmptyCellIndex)
while [] in HistFacValues:
    #print('Error specifying the histogram names of data type 2.')
    EmptyCellIndex = HistFacValues.index([])
    HistRawValues.pop(EmptyCellIndex)
    HistFacValues.pop(EmptyCellIndex)
    NewHistNamesRF.pop(EmptyCellIndex)
lenHn = len(HistNumValues)
lenHr = len(HistRawValues)
HistValues = [[] for i in range(lenHn+lenHr)]

try:
    for HNid in range(lenHn):
	n = len(HistNumValues[HNid])
        if n != len(HistDenValues[HNid]):
            raise Exception
        else:
            for i in range(n):
                if HistNumValues[HNid][i][0] != HistDenValues[HNid][i][0]:
                    raise Exception
        HistValues[HNid] = [[HistNumValues[HNid][i][0]]+ \
                            list(normdivision(HistNumValues[HNid][i][1],HistDenValues[HNid][i][1], \
                                            HistNumValues[HNid][i][2],HistDenValues[HNid][i][2])) for i in range(n)]
except Exception:
    print('Error: Check inconsistency in data type 1')
    raise

try:
    for HNid in range(lenHn,lenHn+lenHr):
	HNid2 = HNid-lenHn
	n = len(HistRawValues[HNid2])
        n2 = len(HistFacValues[HNid2])
        if n2 > 1:
            if (n % n2) != 0:
                raise Exception
            else:
                for i in range(0, n, n2):
                    for j in range(i+1, i+n2):
                        if HistRawValues[HNid2][j-1][0] != HistRawValues[HNid2][j][0]:
                            raise Exception
        #if 0 in [HistFacValues[HNid2][j][1] for j in range(n2)]:
        #    for i in range(n):
        #        if myiszero(HistRawValues[HNid2][i][1]) and myiszero(HistFacValues[HNid2][i % n2][1]):
        #            HistValues[HNid].append([HistRawValues[HNid2][i][0], 0, 0])
        #        else:
        #            HistValues[HNid].append([HistRawValues[HNid2][i][0]]+ \
        #                                    list(normdivision(HistRawValues[HNid2][i][1], \
        #                                                    HistFacValues[HNid2][i % n2][1], \
        #                                                    HistRawValues[HNid2][i][2], \
        #                                                    HistFacValues[HNid2][i % n2][2])))
        #else:
        #    HistValues[HNid] = [[HistRawValues[HNid2][i][0]]+ \
        #                        list(normdivision(HistRawValues[HNid2][i][1],HistFacValues[HNid2][i % n2][1], \
        #                                        HistRawValues[HNid2][i][2],HistFacValues[HNid2][i % n2][2])) for i in range(n)]
        HistValues[HNid] = [[HistRawValues[HNid2][i][0]]+ \
                            list(normdivision(HistRawValues[HNid2][i][1],HistFacValues[HNid2][i % n2][1], \
                                            HistRawValues[HNid2][i][2],HistFacValues[HNid2][i % n2][2])) for i in range(n)]
except Exception:
    print('Error: Check inconsistency in data type 2')
    raise

NewHistNames = NewHistNamesND + NewHistNamesRF
for HNid in range(lenHn+lenHr):
    fout.write('\n  ----   ' + (NewHistNames[HNid]+'               ')[:15] + '----\n\n')
    n = len(HistValues[HNid])
    for i in range(n):
        fout.write( HistValues[HNid][i][0].rjust(9) \
                   + flt2str(HistValues[HNid][i][1]).rjust(18) \
                   + flt2str(HistValues[HNid][i][2]).rjust(18) + '\n' )
        #if type(HistValues[HNid][i][1]) is float and type(HistValues[HNid][i][2]) is float:
        #    fout.write( HistValues[HNid][i][0].rjust(9) \
        #           + ('%g' % HistValues[HNid][i][1]).rjust(18) \
        #           + ('%g\n' % HistValues[HNid][i][2]).rjust(19) )
        #else:
        #    fout.write( HistValues[HNid][i][0].rjust(9) \
        #           + str(HistValues[HNid][i][1]).rjust(18) \
        #           + str(HistValues[HNid][i][2]).rjust(18) \
        #           + '\n' )

try:
    fout.close()
    os.rename(tmpoutfile,outfile) # rename/move the temporary output file to the specified output file
except OSError:
    print('Unable to generate the output file')

if GenerateCSV:
    for HNid in range(lenHn+lenHr):
        fCSV.write(NewHistNames[HNid] + '\n')
        n = len(HistValues[HNid])
        for i in range(n):
            if type(HistValues[HNid][i][1]) is float and type(HistValues[HNid][i][2]) is float:
                fCSV.write( HistValues[HNid][i][0] \
                       + (',%.12g,' % HistValues[HNid][i][1]) \
                       + ('%.12g' % HistValues[HNid][i][2]) \
                       + '\n' )
            else:
                fCSV.write( HistValues[HNid][i][0] \
                       + ','+str(HistValues[HNid][i][1])+','\
                       + str(HistValues[HNid][i][2]) \
                       + '\n' )
    fCSV.close()
