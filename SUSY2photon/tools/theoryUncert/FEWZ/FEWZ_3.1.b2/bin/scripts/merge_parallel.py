"""
Merge output files of all integration sectors into one file in the run directory
Sector subdirectory structure is expected (the one created by "create_parallel.py")
Output files in sector subdirectory are expected to be named as: <order>.sect#.<filename_extension> and <order>.sect#.pdf.<filename_extension>
Both output file and PDF output file will be merged and they will be named as <order>.<filename_extension> and <order>.pdf.<filename_extension>
Used by "finish.sh"
Usage: python $0 <outdir> <outputfile> [ <initsect> <lastsect> <sectstep> ]
Note: <outdir> should have no trailing '/' and path can be added in front
      <outputfile> should be the expected full name of the merged file, i.e. in the form of <order>.<filename_extension>
      PDF output files will be neglected if not present
      <initsect> <lastsect> <sectstep> are optional and for advanced users only and only sets up subdirectories in the run directory
                                       for sectors given by the loop: FOR sector_num FROM <initsect> TO <lastsect> WITH STEP <sectstep>
"""

import sys
import math
from defs import *

try:
    if len(sys.argv) < 3:
        raise Exception
    else:
        outdirloc = sys.argv[1]
        outdirname = pathseparator(outdirloc)[1]
        outnamesplit = sys.argv[2].split('.',1)
        sectors = getsects(outdirloc + '/job_desc')
        if len(sys.argv) > 5:
            initsect = int(sys.argv[3])
            lastsect = int(sys.argv[4])
            sectstep = int(sys.argv[5])
            if sectstep != 0 and initsect > 0 and lastsect > 0 and initsect <= sectors and lastsect <= sectors:
                sectlist = range(initsect, lastsect+sectstep, sectstep)
            else:
		print("Sector loop out of range, default to looping through all sectors")
                sectlist = range(1, sectors+1)
        else:
            sectlist = range(1, sectors+1)
        listsize = len(sectlist)
        outfile1 = outnamesplit[1]
        outfile2 = 'pdf.' + outnamesplit[1]
        outfile3 = 'p_' + outnamesplit[1] 
        outfile4 = 'm_' + outnamesplit[1]
        if sectors == 1:
            print("Only 1 sector exists, no merging necessary, no file is generated")
            sys.exit()
except Exception:
    print('Missing arguments.')
    raise

for outfile in [outfile1, outfile3, outfile4, outfile2]:

    try:
        fin = [open(outdirloc + '/' + outdirname + '%i/' % (i-1) + \
                    outnamesplit[0] + '.sect%(sect)i.' % \
                    {'sect': i} + outfile, 'r') \
               for i in sectlist]
    except IOError:
        if outfile == outfile2:
            print("PDF error sector file missing; no PDF error file" \
                  " will be created.")
            sys.exit()
        elif  outfile == outfile3 or  outfile == outfile4:
            continue # no scale variatio files
        else:
            print("Can't open or find file.")
            raise

    try:
        fout = open(outdirloc + '/' + outnamesplit[0]  + '.' + outfile, 'w')
    except IOError:
        print("Can't open new file in run directory for writing, write merged file in current directory.")
        try:
            fout = open(outnamesplit[0] + '.' + outfile, 'w')
        except IOError:
            print("Can't open new file in current directory for writing.")
            raise

    try:
        nextline = [fin[i].readline() for i in range(listsize)]
    except IOError:
        print("Error reading file.")

    badsects = []
    while len(nextline[0]) > 0: # length = 0 if end of file reached
        splitline = nextline[0].split() # parses line of first file into pieces
        if nextline[0][:6] == ' Sigma': # this is cross section line, replace w/ total
            sigmas = [float(nextline[i].split()[3]) for i in range(listsize)]
            try:
                sigma = list_add(sigmas)
            except SectError, secerr:
                sigma = secerr.retnval
                badsects = badsects + [sectlist[j-1] for j in secerr.excpinfo]
            fout.write(nextline[0].replace(nextline[0].split()[3], flt2str(sigma)))
        elif nextline[0][:6] == ' Error': # this is error line, replace w/ total
            errors = [float(nextline[i].split()[3]) for i in range(listsize)]
            try:
                error = quad_add(errors)
            except SectError, secerr:
                error = secerr.retnval
                badsects = badsects + [sectlist[j-1] for j in secerr.excpinfo]
            fout.write(nextline[0].replace(nextline[0].split()[3], flt2str(error)))
        elif nextline[0][:6] == ' chi^2': # ignore chi^2/it line
            pass
        elif len(splitline) > 0 and splitline[0][-1].isdigit():
        # line starts with number, process bin (no short-circuit would throw error!)
            sigmas = [float(nextline[i].split()[1]) for i in range(listsize)]
            errors = [float(nextline[i].split()[2]) for i in range(listsize)]
            try:
                sigma = list_add(sigmas)
            except SectError, secerr:
                sigma = secerr.retnval
                badsects = badsects + [sectlist[j-1] for j in secerr.excpinfo]
            try:
                error = quad_add(errors)
            except SectError, secerr:
                error = secerr.retnval
                badsects = badsects + [sectlist[j-1] for j in secerr.excpinfo]
            fout.write(splitline[0].rjust(9) + flt2str(sigma).rjust(18) \
                                         + flt2str(error).rjust(18) + '\n')
        else: # nothing to change, regurgitate line into output
            fout.write(nextline[0])
        nextline = [fin[i].readline() for i in range(listsize)]

    for i in range(listsize):
        fin[i].close()

# print warning message if necessary
if len(badsects) != 0:
    sys.stdout.write('Warning, bad numerical results found in sector: ')
    for i in set(badsects):
        sys.stdout.write(str(i)+' ')
    sys.stdout.write('\n')
