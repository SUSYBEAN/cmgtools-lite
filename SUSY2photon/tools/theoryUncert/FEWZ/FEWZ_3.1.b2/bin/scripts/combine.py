"""
Perform the +,-,*,%,.,or 'A' operation between two output files
Used by "finish.sh"
Usage: python $0 <outputfile1> <op> <outputfile2> <newfile>
Note: <outputfile#> should be the full file name and path of the file can be added in front
      two output files must be identical in formats (line to line correspondence) 
"""

import sys
import os
from defs import *

###### Begine Pre-defined Variables and Functions ######
PLUSOP = 1 # Add two results
MINUSOP = 2 # Subtract one from the other
TIMESOP = 3 # Multiply one by the other
DIVOP = 4 # Divide one by the other
AVGOP = 5 # Average one with the other
ASYMOP = 6 # Asymmetry (difference over sum)
NORMOP = 7 # Normalization (normalize all histogram bins of file1 by the total cross-section of file2)

def CheckOperator(opchar):
    if opchar == '+':
        return PLUSOP
    elif opchar == '-':
        return MINUSOP
    elif opchar == '*':
        return TIMESOP
    elif opchar == '/':
        return DIVOP
    elif opchar == '.':
        return AVGOP
    elif opchar == 'A':
        return ASYMOP
    elif opchar == 'N':
        return NORMOP
    else:
        print('Invalid operator, defaulting to +')
        return PLUSOP
def CombineSigma(sigmanum1,sigmanum2,errornum1,errornum2,opnum):
    if opnum == PLUSOP:
        return sigmanum1 + sigmanum2
    elif opnum == MINUSOP:
        return sigmanum1 - sigmanum2
    elif opnum == TIMESOP:
        return sigmanum1*sigmanum2
    elif opnum == DIVOP:
        if sigmanum2 != 0:
            return sigmanum1/sigmanum2
        else:
            return 'undefined'
    elif opnum == AVGOP:
        #return (sigmanum1 + sigmanum2)/2
        if errornum1 == 0 and errornum2 == 0:
            return (sigmanum1 + sigmanum2)/2
        elif errornum1 == 0:
            return sigmanum1
        elif errornum2 == 0:
            return sigmanum2
        else:
            opmwt = mysquare(errornum2)/(mysquare(errornum1)+mysquare(errornum2)) 
            return opmwt*sigmanum1+(1-opmwt)*sigmanum2
    elif opnum == ASYMOP:
        denom = sigmanum1 + sigmanum2
        if denom != 0:
            return (sigmanum1 - sigmanum2)/denom
        else:
            return 'undefined'
    else:
        return 'undefined'
def CombineError(sigmanum1,sigmanum2,errornum1,errornum2,opnum):
    if opnum == PLUSOP or opnum == MINUSOP:
        return mysqrt(mysquare(errornum1) + mysquare(errornum2))
    elif opnum == TIMESOP:
        return mysqrt(mysquare(sigmanum1*errornum2) + mysquare(sigmanum2*errornum1))
    elif opnum == DIVOP:
        if sigmanum2 != 0:
            return mysqrt(mysquare(errornum1/sigmanum2) + mysquare(errornum2*sigmanum1/mysquare(sigmanum2)))
        else:
            return 'undefined'
    elif opnum == AVGOP:
        #return mysqrt(mysquare(errornum1) + mysquare(errornum2))/2
        if errornum1 == 0 or errornum2 == 0:
            return 0
        else:
            return errornum1*errornum2/mysqrt( mysquare(errornum1)+mysquare(errornum2) )
    elif opnum == ASYMOP:
        denom = sigmanum1 + sigmanum2
        if denom != 0:
            return mysqrt(mysquare(2*sigmanum2*errornum1) + mysquare(2*sigmanum1*errornum2))/mysquare(denom)
        else:
            return 'undefined'
    else:
        return 'undefined'
###### End Pre-defined Variables and Functions ######

try:
    if len(sys.argv) < 4:
        raise Exception
    else:
        outfile1 = sys.argv[1]
        op = CheckOperator(sys.argv[2])
        outfile2 = sys.argv[3]
        newfile = sys.argv[4]
except Exception:
    print('Missing or invalid arguments.')
    raise


try:
    fin1 = open(outfile1, 'r')
    fin2 = open(outfile2, 'r')
except IOError:
    print("Can't open or find file.")
    raise

try:
    tmpfile = newfile
    while os.path.isfile(tmpfile):
        tmpfile = tmpfile + '_'
    fout = open(tmpfile, 'w')
except IOError:
    print("Can't open new file for writing.")
    raise

try:
    nextline1 = fin1.readline()
    nextline2 = fin2.readline()
except IOError:
    print("Error reading file.")

histonamebuf = []
if op == NORMOP:
    histonorm = []
    normseq = 0
    file2done = False
if op == AVGOP:
    NoHistoSkip = True
else:
    NoHistoSkip = False
while len(nextline1) > 0 and len(nextline2) > 0: # length = 0 if end of file reached
    splitline1 = nextline1.split() # parses line of first file into pieces
    splitline2 = nextline2.split()
    if nextline1[:6] == ' Sigma' or nextline2[:6] == ' Sigma': 
        # this is cross section line, replace w/ total
        while nextline1[:6] != ' Sigma':
            nextline1 = fin1.readline()
            splitline1 = nextline1.split()
        while nextline2[:6] != ' Sigma':
            nextline2 = fin2.readline()
            splitline2 = nextline2.split()
        # allow dis-aligned files to match
        sigma1 = float(splitline1[3])
        sigma2 = float(splitline2[3])
        # fout.write(nextline1.replace(splitline1[3], '%f' % sigma))
        # hold the output in buffer until error is read
        sigmaline = nextline1
        originalsig = splitline1[3]
        writebuffer = ''
        nextline1 = ''
        while nextline1[:6] != ' Error' and nextline2[:6] != ' Error':
            writebuffer = writebuffer + nextline1
            nextline1 = fin1.readline()
            nextline2 = fin2.readline()
        splitline1 = nextline1.split()
        splitline2 = nextline2.split()
        continue
    elif nextline1[:6] == ' Error' or nextline2[:6] == ' Error':
        # this is error line, replace w/ total
        while nextline1[:6] != ' Error':
            nextline1 = fin1.readline()
            splitline1 = nextline1.split()
        while nextline2[:6] != ' Error':
            nextline2 = fin2.readline()
            splitline2 = nextline2.split()
        # allow dis-aligned files to match
        error1 = float(splitline1[3])
        error2 = float(splitline2[3])
        # assumes we have already encountered sigma line
        if op == NORMOP:
            sigma = CombineSigma(sigma1,sigma2,error1,error2,DIVOP)
            error = CombineError(sigma1,sigma2,error1,error2,DIVOP)
        else:
            sigma = CombineSigma(sigma1,sigma2,error1,error2,op)
            error = CombineError(sigma1,sigma2,error1,error2,op)
        # release buffered output
        fout.write( sigmaline.replace(originalsig, flt2str(sigma)) + writebuffer )
        fout.write( nextline1.replace(splitline1[3], flt2str(error)) )
	# store sigma2 for special renormalized histogram procedure
	if op == NORMOP:
	    histonorm = histonorm + [[sigma2,error2]]
    elif nextline1[:6] == ' chi^2' or nextline2[:6] == ' chi^2': # ignore chi^2/it line
        pass
    elif op == NORMOP and \
         ( len(splitline1) >0 and splitline1[0] == '----' ) :
        normseq = 0 # reset it whenever we move to the next histogram
        file2done = True
        fout.write(nextline1)
    elif op == NORMOP and \
         ( len(splitline1) > 0 and splitline1[0][-1].isdigit() ) :
        sigma1 = float(splitline1[1])
        error1 = float(splitline1[2])
        sigma = CombineSigma(sigma1,histonorm[normseq][0],error1,histonorm[normseq][1],DIVOP)
        error = CombineError(sigma1,histonorm[normseq][0],error1,histonorm[normseq][1],DIVOP)
        if normseq < len(histonorm)-1:
	    normseq = normseq + 1
        else:
	    normseq = 0
        fout.write( splitline1[0].rjust(9) )
        fout.write( flt2str(sigma).rjust(18) )
        fout.write( flt2str(error).rjust(18) + '\n' )
    elif op != NORMOP and (\
           ( len(splitline1) >0 and splitline1[0] == '----' ) or \
           ( len(splitline2) >0 and splitline2[0] == '----' ) \
         ):
        # first make sure file 1 reaches a histogram name
        while len(nextline1)>0 and not ( len(splitline1) >0 and splitline1[0] == '----' ):
            fout.write(nextline1)
            nextline1 = fin1.readline()
            splitline1 = nextline1.split()
        # align histogram title in case two files have different number of histograms assuming histograms appear in the same order
        if nextline1 != nextline2: # misalignment found
            if len(histonamebuf) != 0: # file 2 complete search has been done once and histogram names are stored in buffer
                if nextline1 in histonamebuf: # match found, simply skip lines until reach the matched histogram
                    while nextline1 != nextline2 and len(nextline2) > 0:
                        if NoHistoSkip:
                            fout.write(nextline2)
                        nextline2 = fin2.readline()
                    splitline2 == nextline2.split()
                    fout.write(nextline1)
                else: # no match, skip lines on file 1 to the histogram that is a match
                    while len(nextline1) > 0 and ( not nextline1 in histonamebuf ):
                        if NoHistoSkip:
                            fout.write(nextline1)
                        nextline1 = fin1.readline()
                        splitline1 = nextline1.split()
                        while len(nextline1)>0 and not ( len(splitline1) >0 and splitline1[0] == '----' ):
                            if NoHistoSkip:
                                fout.write(nextline1)
                            nextline1 = fin1.readline()
                            splitline1 = nextline1.split()
                    if len(nextline1) != 0:
                        fout.write(nextline1)
                    else: # no match found and file 1 ends, break
                        break
            else:
                fincrtpos = fin2.tell() # record file 2 current position
                histonamebuf.append(nextline2) # save the current line in case that match is found
                while nextline1 != nextline2 and len(nextline2) > 0: # search file 2 and try to find a match
                    nextline2 = fin2.readline()
                    splitline2 = nextline2.split()
                    if len(splitline2) > 0 and splitline2[0] == '----':
                        histonamebuf.append(nextline2)
                if len(nextline2) != 0: # match found, clear histogram name buffer
                    if NoHistoSkip:
                        fin2.seek(fincrtpos) # go back, complication due to that search did not output skipped histograms 
                        nextline2 = histonamebuf[0]
                        while nextline1 != nextline2:
                            fout.write(nextline2)
                            nextline2 = fin2.readline()
                    histonamebuf = []
                    splitline2 = nextline2.split()
                    fout.write(nextline1)
                else: # no match, go back to original position of file 2 and skip lines on file 1 to the histogram that is a match
                    fin2.seek(fincrtpos)
                    nextline2 = histonamebuf[0]
                    splitline2 = nextline2.split()
                    while len(nextline1) > 0 and ( not nextline1 in histonamebuf ):
                        if NoHistoSkip:
                            fout.write(nextline1)
                        nextline1 = fin1.readline()
                        splitline1 = nextline1.split()
                        while len(nextline1)>0 and not ( len(splitline1) >0 and splitline1[0] == '----' ):
                            if NoHistoSkip:
                                fout.write(nextline1)
                            nextline1 = fin1.readline()
                            splitline1 = nextline1.split()
                    if len(nextline1) != 0:
                        fout.write(nextline1)
                    else: # no match found and file 1 ends, break
                        break
        else:
            fout.write(nextline1)
    elif op != NORMOP and ( \
	  ( len(splitline1) > 0 and splitline1[0][-1].isdigit() ) or \
          ( len(splitline2) > 0 and splitline2[0][-1].isdigit() ) \
        ):
        # line starts with number, process bin (no short-circuit would throw error!)
        while not ( len(splitline1) > 0 and splitline1[0][-1].isdigit() ):
            nextline1 = fin1.readline()
            splitline1 = nextline1.split()
        while not ( len(splitline2) > 0 and splitline2[0][-1].isdigit() ):
            nextline2 = fin2.readline()
            splitline2 = nextline2.split()
        # allow dis-aligned files to match
        sigma1 = float(splitline1[1])
        sigma2 = float(splitline2[1])
        error1 = float(splitline1[2])
        error2 = float(splitline2[2])
        sigma = CombineSigma(sigma1,sigma2,error1,error2,op)
        error = CombineError(sigma1,sigma2,error1,error2,op)
        fout.write( splitline1[0].rjust(9) )
        fout.write( flt2str(sigma).rjust(18) )
        fout.write( flt2str(error).rjust(18) + '\n' )
    else: # nothing to change, regurgitate line into output
        fout.write(nextline1)
    nextline1 = fin1.readline()
    if op != NORMOP or (op == NORMOP and (not file2done)):
        nextline2 = fin2.readline()
fin1.close()
fin2.close()
fout.close()
if tmpfile != newfile:
    os.rename(tmpfile,newfile) # rename/move the temporary output file to the specified output file
 
