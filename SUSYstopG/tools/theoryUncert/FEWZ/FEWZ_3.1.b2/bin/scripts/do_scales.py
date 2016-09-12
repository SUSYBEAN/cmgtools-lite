"""
Perform the +,-,*,%,.,or 'A' operation between two output files
Used by "finish.sh"
Usage: python $0 <outputfile> 
   or  python $0 <outputfile> <new_outputfile>
   or  python $0 <outputfile> <large_sclae_outputfile> <small_scale_outputfile>
   or  python $0 <outputfile> <large_sclae_outputfile> <small_scale_outputfile> <new_outputfile>
Note: <outputfile> should be the full file name and path of the file can be added in front
      three output files must be identical in formats (line to line correspondence) 
"""

import sys
import os
from defs import *

try:
    if len(sys.argv) < 2:
        raise Exception
    else:
        fname = sys.argv[1]
        fnameloclist = pathseparator(fname)
        if len(sys.argv) < 4:
            try:
                fsplit = fnameloclist[1].split('.',1)
                pfname = fnameloclist[0] + fsplit[0] + '.p_' + fsplit[1]
                mfname = fnameloclist[0] + fsplit[0] + '.m_' + fsplit[1]
            except IndexError:
                print("Filename does not contain order prefix; cannot tell the default scale variation file name")
                raise
        if len(sys.argv) >= 4:
            pfname = sys.argv[2]
            mfname = sys.argv[3]
        if len(sys.argv) == 3:
            newfname = sys.argv[2]
        elif len(sys.argv) >= 5:
            newfname = sys.argv[4]
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
    finp = open(pfname, 'r')
    finm = open(mfname, 'r')
except IOError:
    print("No scale variation file; output contains no scale variation.")
    fin.close()
    sys.exit()

try:
    nextline = fin.readline()
    pnextline = finp.readline()
    mnextline = finm.readline()
except IOError:
    print("Error reading file.")

outlines = []
while len(nextline) > 0: # length = 0 if end of file reached
    splitline = nextline.split() # parses line of first file into pieces
    psplitline = pnextline.split()
    msplitline = mnextline.split()
    if nextline[:4] == ' chi': # skip chi^2 line
        pass
    elif nextline.find('Factorization scale') != -1 or nextline.find('Renormalization scale') != -1:
        if 'dynamical' in splitline:
            if splitline[-1] == 'dynamical':
                outlines.append( nextline[:-1].rstrip() + '   x' + psplitline[-1] + '   x' + msplitline[-1] + '\n' )
            else:
                outlines.append( nextline[:(nextline.find('dynamical')+len('dynamical'))] + ' x' + splitline[-1] + \
                                 '   x' + psplitline[-1] + '   x' + msplitline[-1] + '\n' )
        else:
            outlines.append( nextline[:-1].rstrip() + '   ' + psplitline[-1] + '   ' + msplitline[-1] + '\n' )
    elif nextline[:6] == ' Sigma' and pnextline[:6] == ' Sigma' and mnextline[:6] == ' Sigma':
        sigma = myfloat(splitline[3])
        psigma = myfloat(psplitline[3])
        msigma = myfloat(msplitline[3])
        outlines.append(nextline) ### output the central Sigma line
        nextline = fin.readline() ### read Error line for statistical error
        pnextline = finp.readline() ### skip Error line
        mnextline = finm.readline() ### skip Error line
        outlines.append(nextline) ### output the central Error line
        if sigma == 'undefined':
            outlines.append(' Scale Error (pb)            +   ' + 'undefined'.rjust(11) + '\n')
            outlines.append('                             -   ' + 'undefined'.rjust(11) + '\n')
        else:
            if psigma != 'undefined':
                psigma = psigma - sigma
            outlines.append(' Scale Error (pb)            +   ' + flt2str(psigma).rjust(11) + '\n')
            if msigma != 'undefined':
                msigma = sigma - msigma
            outlines.append('                             -   ' + flt2str(msigma).rjust(11) + '\n')
    elif len(splitline) > 3 and splitline[0] == 'bin' and splitline[1] == 'weight' and splitline[2] == 'numerical' and splitline[3] == 'error':
        outlines.append( nextline[:-1] + '+ scale var.'.rjust(18) + '- scale var.\n'.rjust(19) )
    elif ( len(splitline) > 0 and splitline[0][-1].isdigit() ) and \
         ( len(psplitline) > 0 and psplitline[0][-1].isdigit() ) and ( len(msplitline) > 0 and msplitline[0][-1].isdigit() ):
        # line starts with number, process bin (no short-circuit would throw error!)
        sigma = myfloat(splitline[1])
        psigma = myfloat(psplitline[1])
        msigma = myfloat(msplitline[1])
        tmpline = ( nextline[:-1] ) # delete return
        if sigma == 'undefined':
            outlines.append( tmpline + 'undefined'.rjust(18) + 'undefined\n'.rjust(19) )
        else:
            if psigma != 'undefined':
                psigma = psigma - sigma
            tmpline = tmpline + flt2str(psigma).rjust(18)
            if msigma != 'undefined':
                msigma = sigma - msigma
            tmpline = tmpline + flt2str(msigma).rjust(18) + '\n'
            outlines.append( tmpline )
    else: # nothing related to scale variation, regurgitate line into output
        outlines.append( nextline )
    nextline = fin.readline()
    pnextline = finp.readline()
    mnextline = finm.readline()

try:
    fin.close()
    finp.close()
    finm.close()
    fout = open(newfname, 'w')
    fout.writelines(outlines)
    fout.close()
except IOError:
    print("Can't open or find file.")
    raise 
