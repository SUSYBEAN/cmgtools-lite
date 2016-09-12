"""
Handle PDF errors of the output file
A new outputfile
Used by "finish.sh"
Usage: python $0 <outputfile>
   or  python $0 <outputfile> <pdfoutputfile>
   or  python $0 <outputfile> <pdfoutputfile> <newoutputfile>
Note: <outputfile> should be the full file name
      <pdfoutputfile> and <newoutputfile> are set by default rule if they are not given in the argument list
      <pdfoutputfile> is by default set as <order>.pdf.<filename_extension> given <outputfile> in the form of <order>.<filename_extension>
                      (otherwise error comes up), and it is assumed to exist in the same directory as <outputfile>
      <newoutputfile> is by default set to be the file name of <outputfile> (without any path in front) if not given
      <outputfile> will be overwritten if <outputfile> is in the current directory (./) and <newoutputfile> is not specified
"""

import sys
import math
from defs import *

try:
    if len(sys.argv) < 2:
        raise Exception
    else:
        fname = sys.argv[1]
        fnameloclist = pathseparator(fname)
        if len(sys.argv) >= 3:
            pdffname = sys.argv[2]
        else:
            try:
                fsplit = fnameloclist[1].split('.',1)
                pdffname = fnameloclist[0] + fsplit[0] + '.pdf.' + fsplit[1]
            except IndexError:
                print("Filename does not contain order prefix; cannot tell the default PDF file name")
                raise
        if len(sys.argv) >= 4:
            newfname = sys.argv[3]
        else:
            newfname = fnameloclist[1]
except Exception:
    print('Missing arguments.')
    raise

try:
    fin = open(fname, 'r')
except IOError:
    print("Can't open or find file.")
    raise

try:
    finpdf = open(pdffname, 'r')
except IOError:
    print("No PDF error file; output contains no PDF errors.")
    fin.close()
    sys.exit()

try:
    nextline = fin.readline()
    pdfline = finpdf.readline()
except IOError:
    print("Error reading file.")

### Two methods for PDF error calculation:
# 1. HEPDATA method:
#    there is numalfs>=1 central pdf value(s) and pdf errors are calculated through the variation of entries of the pdf eigenvector(s) w.r.t the central one(s)
#    1) sym=0 means pdf eigenvector entries come in pairs and each representing plus and minus PDF errors
#    2) sym=1 menas pdf eigenvector entries varied the PDF both positively and negatively i.e. the PDF error is symmetric (+- the same deviation)
#    3) numalfs>=1 works for both sym=0 and sym=1 cases, and specifies the number of alphas variation(s), of which each variation has its own central pdf and pdf eigenvector;
#       the central pdf result of the central alphas is given in the central output file, the pdf output file first lists the pdf eigenvector result for the central alphas;
#       the central pdf results of the alphas variations are listed after, each followed by their own pdf eigenvector results.
#       note: numalfs=1 could mean that alphas error is already included in entries of the only pdf eigenvector
# 2. Monte Carlo method
#    1) sym=-1: there is no central pdf value, and the mean and standard deviation of all pdf eigenvector entries give separately the pdf "central" value and pdf errors;
#           the central output file only gives the result for the first pdf eigenvector entry;
#    2) numalfs should always be 1 since all variations (either due to pdf eigenvector or alphas changes) will be averaged equally even if alphas variation are supported

outlines = []
sym = 0 # whether pdf error is symmetric: 0=asymmetric, 1=symmetric; defaulting to 0
numalfs = 1 # number of alphas values (for PDF sets with separate alphaS not included in PDF errors); defaulting to 1

while len(nextline) > 0: # length = 0 if end of file reached
    splitline = nextline.split() # parses line of first file into pieces
    if nextline[:4] == ' chi': # skip chi^2 line
        pass
    elif nextline[:8] == ' PDF set': # PDF set line, grab set
        #finpdf.readline() # skip PDF line <--- Don't need it anymore due to while loop after each "elif" for the alignment
        pdfset = splitline[3]
        if pdfset.find('NNPDF') != -1:
            sym = -1 # sym=-1 signals that we should use Monte Carlo method for PDF error calculation
        elif pdfset.find('ABKM09N') != -1:
            sym = 1
	elif pdfset.find('MSTW2008N') != -1 and len(splitline) > 5 and splitline[5] == 'alphaS':
            numalfs = 5
        outlines.append(nextline)
    elif len(splitline) > 3 and splitline[0] == 'bin' and splitline[1] == 'weight' and splitline[2] == 'numerical' and splitline[3] == 'error':
        #the histogram notation line
        if len(splitline) > 4:
            tmpmark = (nextline.find('numerical error')+len('numerical error'))
            strleft = nextline[tmpmark:-1]
        else:
            tmpmark = -1
            strleft = ''
        if numalfs == 1:
            outlines.append( nextline[:tmpmark] + '+ pdf error'.rjust(18) + '- pdf error'.rjust(18) + strleft + '\n' )
        else:
            outlines.append( nextline[:tmpmark] + '+ aS+pdf error'.rjust(18) + '- aS+pdf error'.rjust(18) + \
                             '+ onlypdf err.'.rjust(18) + '- onlypdf err.\n'.rjust(19) + strleft + '\n' )
    elif nextline[:6] == ' Sigma': # this is cross section line, replace w/ total
        nextline = fin.readline() ### read Error line for statistical error; NOTE: starting here nextline and splitline are not related any more
        ### the following two shouldn't be changed unless using Monte Carlo method (sym==-1)
        ### use them to avoid string -> float -> string conversion for central output file
        sigmastr = splitline[3]
        staterrstr = nextline.split()[3]
        ### align pdf variation output file with central value output file
        pdfsplitline = pdfline.split()
        while len(pdfsplitline) == 0 or pdfline[:6] != ' Sigma': 
            pdfline = finpdf.readline()
            pdfsplitline = pdfline.split()
        sigma = myfloat(sigmastr)
        pluserr = 0
        minuserr = 0 ### only used in the case of asymmetric pdf error (sym==0)
        sigmalist = [sigma] ### only used in the case of Monte Carlo method (sym==-1)
        staterrlist = [myfloat(staterrstr)] ### only used in the case of Monte Carlo method (sym==-1)
        ### SPECIAL part to count number of pdf eigenvectors given 'Sigma' appears first
        numpdfs = 0
        tmpstorage = []
        while pdfline[:6] == ' Sigma':
            numpdfs = numpdfs+1
            tmpstorage.append(pdfline)
            tmpstorage.append(finpdf.readline()) ### Also store statistical error line for sym==-1 case
            finpdf.readline()
            pdfline = finpdf.readline()
        # Note: sym == -1: Monte Carlo method central value is regarded also as one eigenvector entry
        #       the true number of pdf eigenvector entries is numpdfs+1 but in the pdf output file there are only numpdfs of them
        if numalfs > 1: # there are multiple chunks each corresponding a different alphaS
            numpdfs = (numpdfs+1)/numalfs-1 # all chunks have a central value in the beginning except for the first one which corresponds to the central alphas value
        if sym == 0: # asymmetric pdf erros come in pairs
            numpdfs = numpdfs/2
        ### use the temporary storage for lines to process pdf error for 'Sigma'
        if sigma =='undefined':
            pluserr = 'undefined'
            minuserr = 'undefined'
        for i in range(numpdfs):
            tmpline = tmpstorage.pop(0)
            sig1 = myfloat(tmpline.split()[3])
            tmpline = tmpstorage.pop(0)
            if sym == 0: # get opposite pdf set if asymmetric errors, otherwise go on
                tmpline = tmpstorage.pop(0)
                sig2 = myfloat(tmpline.split()[3])
                tmpline = tmpstorage.pop(0)
                if pluserr != 'undefined': # pluserr and minuserr must be either both 'undefined' or neither
                    if sig1 == 'undefined' or sig2 == 'undefined':
                        pluserr = 'undefined'
                        minuserr = 'undefined'
                    else:
                        pluserr = pluserr + mysquare(mymax(sig1-sigma,sig2-sigma)) # avoid nan and inf raising OverflowError by using **2
	                minuserr = minuserr + mysquare(mymax(sigma-sig1,sigma-sig2)) # avoid nan and inf raising OverflowError by using **2
            elif sym == 1:
                if pluserr != 'undefined': # pluserr and minuserr must be either both 'undefined' or neither
                    if sig1 == 'undefined':
                        pluserr = 'undefined'
                    else:
                        pluserr = pluserr + mysquare(sig1-sigma) # avoid nan and inf raising OverflowError by using **2
            elif sym == -1:
                sigmalist.append(sig1)
                staterrlist.append(myfloat(tmpline.split()[3]))
        if sym == 0:
            minuserr = mysqrt(minuserr)
        if sym >= 0:
            pluserr = mysqrt(pluserr)
        elif sym == -1:
            sigma,pluserr = my_avg_and_stddev(sigmalist)
            sigmastr = flt2str(sigma)
            staterr = my_quadavg(staterrlist)
            staterrstr = flt2str(staterr)
        outlines.append(' Sigma (pb)                  =    ' + sigmastr + '\n')
        outlines.append(' Integration Error (pb)      =    ' + staterrstr + '\n')
        ### handle alphaS error if it's separate from other pdf eigenvectors, i.e. when numalfs > 1
        if numalfs > 1 and sym >= 0:
            # when sigmas is 'undefined', pervious pluserr and minuserr are both 'undefined'
            pdfpluserr = pluserr
            pdfminuserr = minuserr
            alfpluserr = pluserr
            alfminuserr = minuserr
            if pluserr != 'undefined' and minuserr != 'undefined' :
                for j in range(1,numalfs): # the rest of chunks, corresponding to different non-central alphaS'
                    if alfpluserr == 'undefined': # alfpluserr and alfminuserr must be either both 'undefined' or neither
                        break
                    pluserr = 0
                    minuserr = 0
                    tmpline = tmpstorage.pop(0)
                    alfsigma = myfloat(tmpline.split()[3]) # read central value for the chunk
                    tmpstorage.pop(0)
                    if alfsigma == 'undefined':
                        alfpluserr = 'undefined'
                        alfminuserr = 'undefined'
                        break
                    for i in range(numpdfs):
                        tmpline = tmpstorage.pop(0)
                        sig1 = myfloat(tmpline.split()[3])
                        tmpstorage.pop(0)
                        if sym == 0: # get opposite pdf set if asymmetric errors, otherwise go on
                            tmpline = tmpstorage.pop(0)
                            sig2 = myfloat(tmpline.split()[3])
                            tmpstorage.pop(0)
                            if pluserr != 'undefined': # pluserr and minuserr must be either both 'undefined' or neither
                                if sig1 == 'undefined' or sig2 == 'undefined':
                                    pluserr = 'undefined'
                                    minuserr = 'undefined'
                                else:
                                    pluserr = pluserr + mysquare(mymax(sig1-alfsigma,sig2-alfsigma)) # avoid nan and inf raising OverflowError by using **2
                                    minuserr = minuserr + mysquare(mymax(alfsigma-sig1,alfsigma-sig2)) # avoid nan and inf raising OverflowError by using **2
                            else:
                                break
                        else:
                            if pluserr != 'undefined':
                                if sig1 == 'undefined':
                                    pluserr = 'undefined'
                                else:
                                    pluserr = pluserr + mysquare(sig1-alfsigma) # avoid nan and inf raising OverflowError by using **2
                            else:
                                break
                    if pluserr == 'undefined': # pluserr and minuserr must be either both 'undefined' or neither
                        alfpluserr = 'undefined'
                        alfminuserr = 'undefined'
                        break
                    else:
                        pluserr = mysqrt(pluserr)
                        alfpluserr = mymax(alfpluserr,pluserr+alfsigma-sigma)
                        if sym == 0:
                            minuserr = mysqrt(minuserr)
                            alfminuserr = mymax(alfminuserr,sigma-(alfsigma-minuserr))
            if sym == 0:
                outlines.append(' PDF Error (pb)              +   ' + flt2str(alfpluserr).rjust(11) \
                              + flt2str(pdfpluserr).rjust(12) + ' (w/o alphaS error)\n')
                outlines.append('                             -   ' + flt2str(alfminuserr).rjust(11) \
                              + flt2str(pdfminuserr).rjust(12) + ' (w/o alphaS error)\n')
            else:
                outlines.append(' PDF Error (pb)              =   ' + flt2str(alfpluserr).rjust(11) \
                              + flt2str(pdfpluserr).rjust(12) + ' (w/o alphaS error)\n')
        else:
            if sym == 0:
                outlines.append(' PDF Error (pb)              +   ' + flt2str(pluserr).rjust(11) + '\n')
                outlines.append('                             -   ' + flt2str(minuserr).rjust(11) + '\n')
            else:
                outlines.append(' PDF Error (pb)              =   ' + flt2str(pluserr).rjust(11) + '\n')
    elif len(splitline) > 0 and splitline[0][-1].isdigit():
    # line starts with number, process bin (no short-circuit would throw error!)
        ### the following two shouldn't be changed unless using Monte Carlo method (sym==-1)
        ### use them to avoid string -> float -> string conversion for central output file
        sigmastr = splitline[1]
        staterrstr = splitline[2]
        if len(splitline) > 3:
            strleft = nextline[ (nextline.find(staterrstr)+len(staterrstr)) : -1 ]
            strleft = strleft[ (strleft.find(splitline[3].rjust(18))) :]
        else:
            strleft = ''
        ### align pdf variation output file with central value output file
        pdfsplitline = pdfline.split()
        while len(pdfsplitline) == 0 or pdfsplitline[0] != splitline[0]:
            pdfline = finpdf.readline()
            pdfsplitline = pdfline.split()
        # test if this bin has 'undefined' value
        # 'undefined' will normally only appear in histograms during some normalizatin done by getMomentA.py
        sigma = myfloat(sigmastr)
        tmpline = splitline[0].rjust(9)
        if sigma == 'undefined': ### skip pdf error calculation because sigma is undefined
            if sym == 0:
                for i in range((2*numpdfs+1)*numalfs-1):
                    pdfline = finpdf.readline() ### to skip all
                tmpline = tmpline + sigmastr.rjust(18) + str('undefined').rjust(18) + str('undefined').rjust(18) + str('undefined').rjust(18)
                if numalfs > 1:
                    tmpline = tmpline + str('undefined').rjust(18) + str('undefined').rjust(18)
            else:
                for i in range((numpdfs+1)*numalfs-1):
                    pdfline = finpdf.readline() ### to skip all
                tmpline =  tmpline + sigmastr.rjust(18) + str('undefined').rjust(18) + str('undefined').rjust(18)
                if numalfs > 1:
                    tmpline = tmpline + str('undefined').rjust(18)
        else:
            # proceed to normal procedure
            # calculate the pdf error
            pluserr = 0
            minuserr = 0 ### only used in the case of asymmetric pdf error (sym==0)
            sigmalist = [sigma] ### only used in the case of Monte Carlo method (sym==-1)
            staterrlist = [myfloat(staterrstr)] ### only used in the case of Monte Carlo method (sym==-1)
            pdfundefflag = False
            try:
                for i in range(numpdfs):
                    if sym == 0: # must read line first to guarantee skipping 2 lines even when error occurs
                        pdfline2 = finpdf.readline()
                    sig1 = float(pdfline.split()[1]) # tested for ValueError
                    if sym == 0:
                        sig2 = float(pdfline2.split()[1]) # tested for ValueError
                        pluserr = pluserr + mysquare(mymax(sig1-sigma,sig2-sigma)) # avoid nan and inf raising OverflowError by using **2
                        minuserr = minuserr + mysquare(mymax(sigma-sig1,sigma-sig2)) # avoid nan and inf raising OverflowError by using **2
                    elif sym == 1:
                        pluserr = pluserr + mysquare(sig1-sigma) # avoid nan and inf raising OverflowError by using **2
                    elif sym == -1:
                        sigmalist.append(sig1)
                        staterrlist.append(myfloat(pdfline.split()[2])) # statistical error shouldn't be checked for 'undefined'
                    pdfline = finpdf.readline()
#                pluserr = mysqrt(pluserr)
                if sym == 0:
                    minuserr = mysqrt(minuserr)
                if sym >= 0:
                    pluserr = mysqrt(pluserr)
                elif sym == -1:
                    sigma,pluserr = my_avg_and_stddev(sigmalist)
                    sigmastr = flt2str(sigma)
                    staterr = my_quadavg(staterrlist)
                    staterrstr = flt2str(staterr)
            except ValueError:
                pdfundefflag = True
                if sym == -1:
                    sigmastr = 'undefined'
                    staterrstr = 'undefined'
                for ii in range(i+1,numpdfs):
                    pdfline = finpdf.readline() # skip other pdf lines
                    if sym == 0:
                        pdfline = finpdf.readline()
                pdfline = finpdf.readline() # the beginning of next chunk/bin
            tmpline = tmpline + sigmastr.rjust(18) + staterrstr.rjust(18) #+ strleft 
            # handle alphaS error separately if it's not included in the pdf eigenvectors
            # one exception: skip them all if 'undefined' pdf error is found
            if pdfundefflag:
                # some pre-processing to deal with 'undefined'
                if sym == 0:
                    tmpline = tmpline + 'undefined'.rjust(18) + 'undefined'.rjust(18)
                    if numalfs > 1:
                        for i in range((2*numpdfs+1)*(numalfs-1)):
                            pdfline = finpdf.readline() ### to skip all the rest
                        tmpline = tmpline + 'undefined'.rjust(18) + 'undefined'.rjust(18)
                else:
                    tmpline = tmpline + 'undefined'.rjust(18)
                    if numalfs > 1:
                        for i in range((numpdfs+1)*(numalfs-1)):
                            pdfline = finpdf.readline() ### to skip all the rest
                        tmpline = tmpline + 'undefined'.rjust(18)
            elif numalfs > 1 and sym >= 0:
                # handle alphaS error is calculated separately
                pdfpluserr = pluserr
                pdfminuserr = minuserr
                alfpluserr = pluserr
                alfminuserr = minuserr
                for j in range(1,numalfs):
	            pluserr = 0
	            minuserr = 0
	            try:
                        i = -1 # initialize i in case that there is previous ValueError or alfsigma throws ValueError
                        if pdfundefflag:
                            raise ValueError
                        alfsigma = float(pdfline.split()[1]) # read central value for the chunk
	                pdfline = finpdf.readline()
	                for i in range(numpdfs):
	                    if sym == 0: # must read line first to guarantee skipping 2 lines even when error occurs
	                         pdfline2 = finpdf.readline()
	                    sig1 = float(pdfline.split()[1]) # tested for ValueError
	                    if sym == 0:
	                        sig2 = float(pdfline2.split()[1]) # tested for ValueError
	                        pluserr = pluserr + mysquare(mymax(sig1-alfsigma,sig2-alfsigma)) # avoid nan and inf raising OverflowError by using **2
	                        minuserr = minuserr + mysquare(mymax(alfsigma-sig1,alfsigma-sig2)) # avoid nan and inf raising OverflowError by using **2
	                    else:
	                        pluserr = pluserr + mysquare(sig1-alfsigma) # avoid nan and inf raising OverflowError by using **2
	                    pdfline = finpdf.readline()
                        pluserr = mysqrt(pluserr)
                        alfpluserr = mymax(alfpluserr,pluserr+alfsigma-sigma)
                        if sym == 0:
                            minuserr = mysqrt(minuserr)
                            alfminuserr = mymax(alfminuserr,sigma-(alfsigma-minuserr))
	            except ValueError:
	                pdfundefflag = True
	                for ii in range(i+1,numpdfs):
	                    pdfline = finpdf.readline() # skip other pdf lines
	                    if sym == 0:
	                        pdfline = finpdf.readline()
                        pdfline = finpdf.readline()
		# deal with new 'undefined' pdf error from alphaS error calculation
                if pdfundefflag:
                    if sym == 0:
                        tmpline = tmpline + 'undefined'.rjust(18) + 'undefined'.rjust(18) \
                                          + flt2str(pdfpluserr).rjust(18) + flt2str(pdfminuserr).rjust(18)
                    else:
                        tmpline = tmpline + 'undefined'.rjust(18) + flt2str(pdfpluserr).rjust(18)
                else:
                    if sym == 0:
                        tmpline = tmpline + flt2str(alfpluserr).rjust(18) + flt2str(alfminuserr).rjust(18) \
                                          + flt2str(pdfpluserr).rjust(18) + flt2str(pdfminuserr).rjust(18)
                    else:
                        tmpline = tmpline + flt2str(alfpluserr).rjust(18) + flt2str(pdfpluserr).rjust(18)
            else:
                # normal procedure (alphaS error included in pdf errors)
                if sym == 0:
                    tmpline = tmpline + flt2str(pluserr).rjust(18) + flt2str(minuserr).rjust(18)
                else:
                    tmpline = tmpline + flt2str(pluserr).rjust(18)
	# output to file
        outlines.append(tmpline + strleft + '\n')
    else: # nothing to change, regurgitate line into output
        outlines.append(nextline)
        #pdfline = finpdf.readline() # skip PDF lines <--- Don't need it anymore due to while loop after each "elif" for the alignment
    nextline = fin.readline()

try:
    fin.close()
    fout = open(newfname, 'w')
    fout.writelines(outlines)
    fout.close()
except IOError:
    print("Can't open or find file.")
    raise
