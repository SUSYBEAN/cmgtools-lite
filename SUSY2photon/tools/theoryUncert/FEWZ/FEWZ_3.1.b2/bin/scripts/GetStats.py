"""
Collect numerical integration statistical information of all sectors
Subdirectory structure created by "create_parallel.py" in the run directory are expected
A file named <outputfile>.stats will be generated, recording numerical information from all sectors
Usage: python $0 <rundir> <inputfiles> <outputfile>
Advanced Usage: python $0 <rundir> <inputfile> [<order>.][pdf.]<outputfile> [<whichoccurrence>] [<histogramname> <histogrambin>]
Note: <rundir> must have no trailing '/' and path can be added in front
      arguments in [] are optional
      [<order>] is the order prefix, which must be consistent with the setting in input file and can be neglected
      [pdf.] determines whether to test the output file or the PDF output file
      [<whichoccurrence>] determines to test the numerical result for its <whichoccurrence>-th occurrence in the output file
                          and is useful for the pdf output file which has different results of the same variable for different PDF eigenvectors
      [histogramname] is the keyword appearing in the title of the histogram in the output file
      [histogrambin] is the bin of the histogram, of which one want to collect numerical information
"""

import sys
import math
import os
from datetime import timedelta
from defs import getsects
from defs import pathseparator

try:
    if len(sys.argv) < 4:
        raise Exception
    else:
        outdir = sys.argv[1]
        outdirname = pathseparator(outdir)[1]
        infile = sys.argv[2]
        outfile = sys.argv[3]
        if len(sys.argv) >= 5:
            WhichOccurrence = int(sys.argv[4]) 
        else:
            WhichOccurrence = 1 ### Default to first occurrence
        if len(sys.argv) >= 7:
            checktotal = False
            HistName = sys.argv[5]
            HistBin = float(sys.argv[6])
        else:
            checktotal = True ### Default to only check the total cross section
except Exception:
    print('Missing arguments.')
    print('Usage: python GetStats.py <run_dir> <input_file> <output_file>\n Note: no order prefix for output filename')
    raise

# read information from inputfile
inputfileinfo = ''
try:
    #sectors = getsects(infile)
    tempfile = open(infile, 'r')
    templine = tempfile.readline()
    while len(templine) > 0:
        templine = tempfile.readline()
        if templine.find('Relative accuracy') != -1:
	    inputfileinfo = inputfileinfo + templine
	    #racc = float((templine[39:].replace('d','e')).replace('D','E'))
            racc = float((templine.split()[-1].replace('d','e')).replace('D','E'))
        elif templine.find('Absolute accuracy') != -1:
            inputfileinfo = inputfileinfo + templine
	    #aacc = float((templine[39:].replace('d','e')).replace('D','E'))
            aacc = float((templine.split()[-1].replace('d','e')).replace('D','E'))
        elif templine.find('Maximum number of evaluations') != -1:
            inputfileinfo = inputfileinfo + templine
	    #maxeval = int(templine[39:])
	    maxeval = int(templine.split()[-1])
        elif templine.find('W production') != -1:
            boson = 'W'
        elif templine.find('Z production') != -1:
            boson = 'Z'
        elif templine.find('QCD Perturb. Order') != -1:
            inputfileinfo = inputfileinfo + templine
	    #podr = int(templine[41:])
	    podr = int(templine.split()[-1])
	    sectors = getsects(podr,boson)
	    if podr == 0:
	        oprefix = 'LO'
	    elif podr == 1:
	        oprefix = 'NLO'
	    elif podr == 2:
	        oprefix = 'NNLO'
	    else:
	        print("Unrecognized specification of perturbation order in input file.")
	        raise
            ofsplit = outfile.split('.',1)
            if ofsplit[0] != oprefix: # outfile name given has no or wrong order prefix
                if ofsplit[0] == 'LO' or ofsplit[0] == 'NLO' or ofsplit[0] == 'NNLO':
                    print("Output file name inconsistent with input file perturbation order setting.")
                    raise
            else: # outfile name given has order prefix
                #outfile = '.'.join(ofsplit[1:])
                outfile = ofsplit[1]
	else:
           pass
    tempfile.close()
    aacc = aacc / math.sqrt(sectors)
except IOError:
    print("Can't read input file.")
    raise

fin = []
sectlist = []
for i in range(sectors):
    try:
	sectfilename = outdir + '/' + outdirname + '%i/' % i + \
                    oprefix + '.sect%(sect)i.' % \
                    {'sect': i+1} + outfile
        fin.append(open(sectfilename, 'r')) 
	sectlist.append(i+1)
    except IOError:
	pass #print("Can't open, find or read file: " + sectfilename)
if len(sectlist) == 0:
    raise
else:
    goodsects = len(sectlist)

try:
    nextline = [fin[i].readline() for i in range(goodsects)]
    sigmas = []
    errors = []
    runtimes = []
    actevals = []
    sigmaocr = 0
    #errorocr = 0
    while len(nextline[0]) > 0 and (0 in map(len,[sigmas,errors,runtimes,actevals]) ): # length = 0 if end of file reached
        if nextline[0][:6] == ' Total':
            runtimes = [timedelta(seconds=float(nextline[i].split()[5])).__str__() for i in range(goodsects)]
            #runtimes = [nextline[i].split()[5] for i in range(goodsects)]
        elif nextline[0][:7] == ' Actual':
            actevals = [int(nextline[i].split()[5]) for i in range(goodsects)]
        else:
            if checktotal: # read total cross section
                if nextline[0][:6] == ' Sigma':
                    sigmaocr = sigmaocr + 1
                    if sigmaocr >= WhichOccurrence:
                        sigmas = [float(nextline[i].split()[3]) for i in range(goodsects)]
                elif nextline[0][:6] == ' Error':
                    #errorocr = errorocr + 1
                    #if errorocr == WhichOccurrence:
                    if sigmaocr >= WhichOccurrence:
                        errors = errors + [float(nextline[i].split()[3]) for i in range(goodsects)]
            else: # read certain bin of histograms
                if nextline[0].find(HistName) != -1:
                    nextline = [fin[i].readline() for i in range(goodsects)]
                    while len(nextline[0].split()) == 0:
                        nextline = [fin[i].readline() for i in range(goodsects)]
                    binline = nextline[0].split()
                    if float(binline[0]) == HistBin:
                        sigmaocr = sigmaocr + 1
                    while len(binline) != 0 and ( float(binline[0]) != HistBin or sigmaocr < WhichOccurrence ) :
                        nextline = [fin[i].readline() for i in range(goodsects)]
			binline = nextline[0].split()
                        if float(binline[0]) == HistBin:
                            sigmaocr = sigmaocr + 1
                    if len(binline) != 0:
                        sigmas = [float(nextline[i].split()[1]) for i in range(goodsects)]
                        errors = [float(nextline[i].split()[2]) for i in range(goodsects)]
### Backdoor tool to sum N bins down from the designated one as the total
			nbindown = 1 # N=1
                        for toolcount in range(nbindown-1):
			    nextline = [fin[i].readline() for i in range(goodsects)]
                            for i in range(goodsects):
			        sigmas[i] = sigmas[i]+float(nextline[i].split()[1])
				errors[i] = math.sqrt((errors[i])**2+(float(nextline[i].split()[2]))**2)
### End of the Backdoor tool
        nextline = [fin[i].readline() for i in range(goodsects)]
except IOError:
    print("Error reading file for sector " + '%i' % sectlist[i]) 
except Exception:
    print('Error getting statistics')
    raise

if len(sigmas) == 0 or len(errors) == 0:
    print('No matched data found!')
    raise Exception

localmode = True
for i in sectlist:
    if os.path.isfile(outdir + '/' + outdirname + '%i/' % (i-1) + 'condor_log.log'):
	localmode = False
	break

if not localmode:
    for i in range(goodsects):
        try:
	    #print(outdir + '/' + outdirname + '%i/' % (i-1) + 'condor_log.log')
	    clogin = open(outdir + '/' + outdirname + '%i/' % (sectlist[i]-1) + \
                        'condor_log.log', 'r')
            nextline = clogin.readline()
            parsed = nextline.split()
#----------------- Expected Format in the Condor log file ------------------#
#        (1) Normal termination (return value #)                            #
#                Usr # ##:##:##, Sys # ##:##:##  -  Run Remote Usage        #
#                Usr # ##:##:##, Sys 0 ##:##:##  -  Run Local Usage         #
#                Usr # ##:##:##, Sys 0 ##:##:##  -  Total Remote Usage      #
#                Usr # ##:##:##, Sys 0 ##:##:##  -  Total Local Usage       #
#        #########  -  Run Bytes Sent By Job                                #
#        #########  -  Run Bytes Received By Job                            #
#        ##########  -  Total Bytes Sent By Job                             #
#        ##########  -  Total Bytes Received By Job                         #
#---------------------------------------------------------------------------#
            while not((len(parsed) > 2) and (parsed[1] == 'Normal') and (parsed[2] == 'termination')) \
                  and (len(nextline) > 0): # end-of-file, we did not find stats
                nextline = clogin.readline()
                parsed = nextline.split()       
            if (len(nextline) > 0): # loop terminated normally, we found stats
                clogin.readline()
                clogin.readline() # clear next two lines, we want third
                nextline = clogin.readline()
                parsed = nextline.split()
                if parsed[1] != '0':
                    if parsed[1] == '1':
                        runtimes[i] = parsed[1] + ' day, ' + (parsed[2])[:-1]
                    else:
                        runtimes[i] = parsed[1] + ' days, ' + (parsed[2])[:-1]
                else:
                    runtimes[i] = (parsed[2])[:-1]
            else:
                print("sector %i did not terminate normally / is still running, runtime is estimated from clock time" % (i+1))
	    clogin.close()
        except IOError:
            print("Error reading Condor log for sector %i, runtime is estimated from clock time." % sectlist[i])
        except Exception:
            print("Error in Condor log for sector %i, runtime is estimated from clock time." % sectlist[i])
                

try:
    fout = open(oprefix + '.' + outfile + '.stats', 'w')
except IOError:
    print("Can't open new stats file for writing.")
    raise

fout.write(inputfileinfo)
fout.write('\n')
fout.write('Sector'.rjust(8) + 'Sigma'.rjust(14) + '+/-Error'.rjust(13) \
           + 'Run Time (sec)'.rjust(18) + 'Evaluations'.rjust(15) \
           + 'Done?\n'.rjust(12))
fout.write('------'.rjust(8) + '-----'.rjust(14) + '--------'.rjust(13) \
           + '--------------'.rjust(18) + '-----------'.rjust(15) \
           + '-----\n'.rjust(12))
for i in range(goodsects):
    ifdone = (actevals[i] >= maxeval or errors[i] <= aacc \
              or errors[i] <= racc*math.fabs(sigmas[i])/100.0) 
    fout.write(str(sectlist[i]).rjust(8) + ('%f' % sigmas[i]).rjust(14) \
                                 + ('%f' % errors[i]).rjust(13) \
                                 + runtimes[i].rjust(18) \
                                 + ('%i' % actevals[i]).rjust(15) \
                                 + (str(ifdone) + '\n').rjust(12))
fout.close()
