### Header file for python scripts

import math

### funtion to analyze file name
def pathseparator(filenamestr):
    filesplit = filenamestr.split('/')
    filepath = ''
    if len(filesplit) > 1:
       filepath = '/'.join(filesplit[0:-1]) + '/'
    filename = filesplit[-1]
    return [filepath,filename]

### general function to adapt to different versions of Python
def myisnan(fltnum):
    return str(fltnum) == 'nan'

def myisinf(fltnum):
    return str(fltnum) == 'inf'

def myiszero(fltnum):
    return fltnum == 0 and not myisnan(fltnum)

def myfloat(numstr):
    try:
        return float(numstr)
    except ValueError:
        if numstr == 'undefined':
            return numstr
        return float('nan')

### to avoid OverflowError of x**2 by using x*x in case x=inf
def mysquare(fltnum):
    if fltnum == 'undefined':
        return fltnum
    return fltnum*fltnum

def mysqrt(fltnum):
    try:
        if fltnum == 'undefined':
            return fltnum
        return math.sqrt(fltnum)
    except:
        if myisnan(fltnum) or myisinf(fltnum):
            return fltnum
	else:
	    raise

### convert float number to string
def flt2str(fltnum):
    if fltnum == 'undefined':
        return 'undefined'
    else:
        return ( '%g' % fltnum )

### NNLO sector information
NNLOsectorsW = 154
NNLOsectorsZ = 127

### for merging script error instance
class SectError(Exception):
    def __init__(self, retnval, excpinfo):
        self.retnval = retnval
	self.excpinfo = excpinfo
    def __str__(self):
        return repr(retnval) + repr(self.excpinfo)

### for merging scripts: result (it will throw SectError defined above when non-number appears)
def list_add(ls):
    """list_add adds the elements of a list together and returns the result."""
    try:
        if 'undefined' in ls:
            raise
        tot = sum(ls)
        if myisnan(tot) or myisinf(tot):
            raise
        return tot
    except:
        badmem = []
        for i in range(len(ls)):
            if myisnan(ls[i]) or myisinf(ls[i]):
                badmem.append(i+1)
        raise SectError(tot,badmem)

### for merging scripts: error (it will throw SectError defined above when non-number appears)
def quad_add(ls):
    """quad_add adds the elements of a list in quadrature."""
    try:
        if 'undefined' in ls:
            raise
        tot = sum(map(mysquare,ls))
        tot = math.sqrt(tot)
        if myisnan(tot) or myisinf(tot):
            raise
        return tot
    except:
        badmem = []
        for i in range(len(ls)):
            if myisnan(ls[i]) or myisinf(ls[i]):
                badmem.append(i+1)
        raise SectError(tot,badmem)

### for pdf asymmetric errors
def mymax(num1, num2):
    if num1 == 'undefined' or num2 == 'undefined':
            return 'undefined'
    return max(num1, num2, 0)

### for calculating sigma and PDF error when using Monte Carlo method to obtain PDF error
def my_avg_and_stddev(ls):
    if 'undefined' in ls:
        return 'undefined','undefined'
    num = len(ls)
    if num == 0:
        print "return 0 and 0 as average and standard deviation for the empty list"
        return 0,0
    avg = sum(ls)/num
    stddev = 0.0
    if num > 1:
        for i in range(num):
            stddev = stddev + mysquare(ls[i]-avg)
        stddev = stddev/(num-1)
        stddev = mysqrt(stddev)
    return avg,stddev

### for calculating statistical error when using Monte Carlo method to obtain PDF error
def my_quadavg(ls):
    if 'undefined' in ls:
        return 'undefined','undefined'
    num = len(ls)
    if num == 0:
        print "return 0 as quadratic average for the empty list"
        return 0
    if num == 1:
        return ls[0]
    quadavg = sum(map(mysquare,ls))/num
    quadavg = mysqrt(quadavg)
    return quadavg 

### for get A moment script
def normdivision(num=0,den=0,numerr=0,denerr=0):
    try:
        if num == 'undefined' or den == 'undefined':
            return 'undefined','undefined'
        if numerr == 'undefined' or denerr == 'undefined':
            return num/den,'undefined'
        diverr = mysquare(numerr/den)+mysquare(denerr*num/mysquare(den))
	diverr = mysqrt(diverr)
        return num/den, diverr
    except ZeroDivisionError:
        if myiszero(num):
            return 0, 'undefined'
        return 'undefined', 'undefined'

## spit out sector number according to given perturbation order
def getsects(namestr, bosonchar='Z'):
    # namestr coule be:
    #    1. name string of the input file
    #    2. name string of job description file / condor submission file
    #    3. string "LO", "NLO" or "NNLO"
    #    4. integer 0, 1 or 2 representing LO, NLO or NNLO
    # bosonchar could be:
    #    1. 'W' or 'w' representing W boson
    #    2. Anything else defaulting to Z boson
    #    3. Can be neglected when namestr is a file
    # getsects(namestr, bosonchar) returns the corresponding number of sectors
    bosonid = bosonchar
    try:
	### pretend it is a file and try to open it
        fin = open(str(namestr), 'r')
        #for nextline in fin:
        #    if nextline[:5] == 'Queue':
        #        return int(nextline.split()[1])
        #    if nextline[:5] == '\'Perturb.':
        #        break
        nextline = fin.readline()
        while nextline[:5] != 'Queue' and nextline[:13] != '\'QCD Perturb.' and len(nextline) != 0:
	    if nextline[:13] == '\'W production':
                bosonid='W'
            elif nextline[:13] == '\'Z production':
                bosonid='Z'
            nextline = fin.readline()
        if nextline[:13] == '\'QCD Perturb.': ### read the next line for possible EW order
            nextline2 = fin.readline() 
        fin.close()
	if len(nextline) == 0:
	    raise IOError
	else:
            if nextline[:5] == 'Queue':
                return int(nextline.split()[1]) 
            elif nextline[:13] == '\'QCD Perturb.':
                order = int(nextline.split()[-1])
                if nextline2[:12] == '\'EW Perturb.':
                    order = max(order,int(nextline2.split()[-1]))
    except IOError:
	if type(namestr) is int:
	    order = namestr
	elif namestr == 'NNLO':
	    order = 2
	elif namestr == 'NLO':
	    order = 1
	elif namestr == 'LO':
	    order = 0
	else:
	    print("illegal input of getsects()")
	    raise
    if order <= 1:
        return 1
    else:
        if bosonid == 'W' or bosonid == 'w':
            return NNLOsectorsW
        else:
            return NNLOsectorsZ

## read the factorization and normalization scale rescale factor for scale variation
def getscalevar(infile):
    # infile is name string of the input file
    # getscalevar(infile) returns scale (mu_r and mu_f) variation factor of the input file in format:
    fin = open(infile, 'r')
    nextline = fin.readline()
    mufpfct = 1
    mufmfct = 1
    murpfct = 1
    murmfct = 1
    while ( mufpfct == 1 or murpfct == 1 ) and len(nextline) != 0:
        while ( nextline.find('Factorization scale') == -1 and nextline.find('Renormalization scale') == -1 ) and len(nextline) != 0:
            nextline = fin.readline()  
        if len(nextline) == 0:
            break
        splitline = nextline.split("'")
        splitline = splitline[2].split()
        if nextline.find('Factorization scale') != -1:
            muf = float(splitline[0].replace('d','e').replace('D','e'))
        elif nextline.find('Renormalization scale') != -1:
            mur = float(splitline[0].replace('d','e').replace('D','e'))
        if len(splitline) > 1:
            if len(splitline) == 2:
                pfct = float(splitline[1].replace('d','e').replace('D','e'))
                if pfct != 0:
                    mfct = 1/pfct
            elif len(splitline) > 2:
                pfct = float(splitline[1].replace('d','e').replace('D','e'))
                mfct = float(splitline[2].replace('d','e').replace('D','e'))
            if nextline.find('Factorization scale') != -1 and pfct > 0 and mfct > 0:
                muf = float(splitline[0].replace('d','e').replace('D','e'))
                mufpfct = pfct
                mufmfct = mfct
            elif nextline.find('Renormalization scale') != -1 and pfct > 0 and mfct > 0:
                mur = float(splitline[0].replace('d','e').replace('D','e'))
                murpfct = pfct
                murmfct = mfct
        nextline = fin.readline()
    if mufpfct == 1 and murpfct == 1 and mufmfct == 1 and murmfct == 1:
        return [[muf,mur],[muf,mur],False]
    else:
        if muf == 0:
            muf = -1
        if mur == 0:
            mur = -1
        return [[muf*mufpfct,mur*murpfct],[muf*mufmfct,mur*murmfct],True]
