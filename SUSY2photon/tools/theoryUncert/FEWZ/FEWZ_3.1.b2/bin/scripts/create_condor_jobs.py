"""
Produce condor submission file "job_desc"
Run "condor_submit job_desc" to submit
Used by "condor_run.sh"
Usage: python $0 z/w <outputdir> <inputfile> <histofile> <outputfile_extension> <pdfdir> [<jobdescfile> [<initsect> [<lastsect> <sectstep>]]]
Note: <outputdir> should be a directory name without a trailing '/' and a preceding path
      <pdfdir> should be a directory name or path without trailing '/'
"""

import sys
from defs import getsects
from defs import getscalevar

jobname = 'job_desc'
condor_out = 'condor_output.out'
condor_err = 'condor_error.err'
condor_log = 'condor_log.log'

try:
    if len(sys.argv) < 7:
        raise Exception

    else:
        boson = (sys.argv[1]).upper()
        outputdir = sys.argv[2]
        inputfile = sys.argv[3]
        sectors = getsects(sys.argv[3], boson)
        scalevars = getscalevar(sys.argv[3]) # see whether input file contain special scale variation instruction
        histofile = sys.argv[4]
        outputfile = sys.argv[5]
        pdfdir = sys.argv[6]
        if len(sys.argv) > 7:
            jobname = sys.argv[7]

        if len(sys.argv) > 8:
            initsect = int(sys.argv[8])
            if len(sys.argv) > 10:
               lastsect = int(sys.argv[9])
               sectstep = int(sys.argv[10])
            else:
               lastsect = initsect
               sectstep = 1
        else:
            initsect = 1
            lastsect = sectors
            sectstep = 1
        sectlist = range(initsect, lastsect+sectstep, sectstep)

except Exception:
    print('Missing arguments.')
    raise

try:
    job_file = open(outputdir + '/' + jobname, 'w')

except IOError:
    print('Error creating job file.')
    raise

try:
    if (boson == 'W'):
        job_file.write('Executable = condor_fewzw\n')
    else:
        if (boson != 'Z'):
            print('Warning: unrecognized parameter; defaulting to neutral current\n')
        job_file.write('Executable = condor_fewzz\n')
    
    job_file.write('Environment = CUBACORES=0\n')
    job_file.write('Universe = standard\n')
    job_file.write('Requirements = Memory >= 128 && Disk >= 5000\n')
    job_file.write('Rank = kflops\n')
    job_file.write('Getenv = True\n')
    job_file.write('Notification = Error\n')
    job_file.write('Output = ' + condor_out + '\n')
    job_file.write('Log = ' + condor_log + '\n')
    job_file.write('Error = ' + condor_err + '\n')
#    job_file.write('+WantCheckpoint = False\n')

    job_file.write('\n')
    if len(sys.argv) > 8:
        for i in sectlist:
            job_file.write('arguments = -i ../' + inputfile + ' -h ../' + histofile \
                           + ' -o ' + outputfile + ' -p ../../' + pdfdir + ' -s %i\n' % (i-1))
            job_file.write('initialdir = ' + outputdir + '%i\n' % (i-1))
            job_file.write('queue\n')
    else:
        job_file.write('Arguments = -i ../' + inputfile + ' -h ../' + histofile \
                       + ' -o ' + outputfile + ' -p ../../' + pdfdir + ' -s $(Process)\n')
        job_file.write('Initialdir = ' + outputdir + '$(Process)\n')
        job_file.write('Queue %i\n' % sectors)

    if scalevars[2]:
        # for scale variation plus
        job_file.write('\n')
        for i in sectlist:
            job_file.write('arguments = -i ../../p_' + inputfile + ' -h ../../' + histofile \
                           + ' -o p_' + outputfile + ' -p ../../../' + pdfdir + ' -l ..' + ' -s %i\t\n' % (i-1))
            job_file.write('initialdir = ' + outputdir + '%i/pscale\t\n' % (i-1))
            job_file.write('queue\t\n')
        # for scale variation minus
        job_file.write('\n')
        for i in sectlist:
            job_file.write('arguments = -i ../../m_' + inputfile + ' -h ../../' + histofile \
                           + ' -o m_' + outputfile + ' -p ../../../' + pdfdir + ' -l ..' + ' -s %i\t\n' % (i-1))
            job_file.write('initialdir = ' + outputdir + '%i/mscale\t\n' % (i-1))
            job_file.write('queue\t\n')

    job_file.close()

except IOError:
    print('Error writing to job file.')
    raise






