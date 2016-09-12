"""
Find if files specifying bin boundaries need to be copied, and alter histogram file to point to correct file in run directory structure.
Usage: python $0 <histofile> <outputdir>
"""

import sys
import os
import shutil

try:
    if len(sys.argv) < 2:
        raise Exception
    else:
        histofilen = sys.argv[1]
        outdir = sys.argv[2]
        newfilen = outdir + '/' + histofilen
except Exception:
    print('Missing arguments.')
    raise

try:
    newfile = open(newfilen, 'w')
except OSError:
    print('Error creating new histogram file.')
    raise

try:
    histofile = open(histofilen, 'r')
except OSError:
    print('Error reading histogram file.')
    raise

try:
    nextline = histofile.readline()
    while len(nextline) > 0:
        splitline = nextline.split('\'')
        if (len(splitline) > 4): # will read as 5 splits if two Fortran strings
            binfilen = splitline[3]
# if this script is called, bin file will be one directory above where executable is running
            nextline = nextline.replace(binfilen, '../' + binfilen)
            shutil.copy(binfilen, outdir + '/' + binfilen)
        newfile.write(nextline)
        nextline = histofile.readline()
    newfile.close()
except IOError:
    print('Error reading histogram file.')
    raise
