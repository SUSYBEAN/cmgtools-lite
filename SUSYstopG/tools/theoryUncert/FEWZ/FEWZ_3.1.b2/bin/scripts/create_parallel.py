"""
Create the subdirectory structure for all sectors in the run directory and
Generate different input files for different sectors in the corresponding subdirectories
Used by "condor_run.sh" and "local_run.sh"
Usage: python $0 w/z <inputfile> <outputdir> [ <initsect> <lastsect> <sectstep> ]
Note: <outputdir> should be a directory name without a trailing '/' and a preceding path
      By default subdirectory for all sectors of the corresponding perturbation order are set up in the run directory
      <initsect> <lastsect> <sectstep> are optional and for advanced users only and only sets up subdirectories in the run directory
                                       for sectors given by the loop: FOR sector_num FROM <initsect> TO <lastsect> WITH STEP <sectstep>
      It is always safe to ignore the optional arguments and set up the subdirectory structure all at once, i.e. using $0 <inputfile> <outputdir>
"""

import sys
import os
from defs import pathseparator
from defs import getsects
from defs import getscalevar

#def modify_inputs(inputfile, newfiles, sectlist, listsize):
#    """Parses an input file, and writes a new input file with sector i."""
#    try:
#        nextline = inputfile.readline()
#    except IOError:
#        print('Error reading input file.')
#    try:
#        while len(nextline) > 0: # length = 0 at end of file
#            if nextline[:13] == '\'All Sectors?': # ensure this line does one sector
#                for i in range(listsize):
#                    newfiles[i].write('\'All Sectors?		      = \' 0\n')
#            elif nextline[:12] == '\'NNLO sector': # do one sector at a time
#                for i in range(listsize):
#                    newfiles[i].write('\'NNLO sector                  = \' %i\n' \
#                                      % sectlist[i])
#            else:
#                for i in range(listsize):
#                    newfiles[i].write(nextline)
#            nextline = inputfile.readline()
#    except IOError:
#        print('Error reading or writing file.')

def modify_scales(inputfile, newfile, muf, mur):
    """Parses an input file, and writes a new input file with scale specified in muf and mur."""
    try:
        nextline = inputfile.readline()
    except IOError:
        print('Error reading input file.')
    try:
        while len(nextline) > 0: # length = 0 at end of file
            if nextline.find('Factorization scale') != -1:
                newfile.write('\'Factorization scale  (GeV)    = \' ' + ('%f' % muf)   + 'd0\n')
            elif nextline.find('Renormalization scale') != -1:
                newfile.write('\'Renormalization scale  (GeV)    = \' ' + ('%f' % mur)   + 'd0\n')
            elif nextline.find('Turn off PDF error') != -1:
                newfile.write('\'Turn off PDF error (1=Yes, 0=No)    = \' 1\n')
            else:
                newfile.write(nextline)
            nextline = inputfile.readline()
    except IOError:
        print('Error reading or writing file.')


try:
    if len(sys.argv) < 4:
        raise Exception
    else:
        dirlocation = sys.argv[3]
        dirname = pathseparator(dirlocation)[1]
        sectors = getsects(sys.argv[2], sys.argv[1].upper())
        scalevars = getscalevar(sys.argv[2]) # see whether input file contain special scale variation instruction
	if len(sys.argv) > 6:
	    initsect = int(sys.argv[4])
	    lastsect = int(sys.argv[5])
	    sectstep = int(sys.argv[6])
	    if sectstep != 0 and initsect > 0 and lastsect > 0 and initsect <= sectors and lastsect <= sectors:
	        sectlist = range(initsect, lastsect+sectstep, sectstep)
	    else:
		sectlist = range(1, sectors+1)
	else:
	    sectlist = range(1, sectors+1)
	listsize = len(sectlist)
except IOError:
    print('Error reading input file.')
    raise
except Exception:
    print('Missing arguments.')
    raise

try:
    if not os.path.isdir(dirlocation):
	os.mkdir(dirlocation)
    for i in sectlist:
        subdirloc = dirlocation + '/' + dirname + '%i' % (i-1)
	if not os.path.isdir(subdirloc):
            os.mkdir(subdirloc)
            if scalevars[2]:
                os.mkdir(subdirloc + '/pscale')
                os.mkdir(subdirloc + '/mscale')
except OSError:
    print("Error creating directories.")

try:
#  ## write new input file for each sector
#    newfiles = [open(dirlocation + '/' + dirname + '%i/' % (i-1) + sys.argv[2], 'w') \
#                for i in sectlist]
#    inputfile = open(sys.argv[2], 'r')
#    modify_inputs(inputfile, newfiles, sectlist, listsize)
#    for i in range(listsize):
#	newfiles[i].close()
#    inputfile.close()
    ## write new input files for the scale variation plus and minus
    if scalevars[2]:
        newfile = open(dirlocation+'/p_'+sys.argv[2], 'w')
        inputfile = open(sys.argv[2], 'r')
        modify_scales(inputfile, newfile, scalevars[0][0], scalevars[0][1])
        newfile.close()
        inputfile.close()
        newfile = open(dirlocation+'/m_'+sys.argv[2], 'w')
        inputfile = open(sys.argv[2], 'r')
        modify_scales(inputfile, newfile, scalevars[1][0], scalevars[1][1])
        newfile.close()
        inputfile.close()
    ## write out number of sectors (in file "sectors") for later in case Condor is not used
    jobdescfile = open(dirlocation + '/job_desc', 'w')
    jobdescfile.write('Queue %i\n' % sectors)
    jobdescfile.close()
except IOError:
#    print('Error creating input files.')
    raise

