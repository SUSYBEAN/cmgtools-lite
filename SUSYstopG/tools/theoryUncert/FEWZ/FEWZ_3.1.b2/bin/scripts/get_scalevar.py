"""
Call python getscalevar() function in "defs.py" and print on screen whether there are scale variation job running
Usage: python $0 <inputfile>
"""

import sys
from defs import getscalevar

scalevars = getscalevar(sys.argv[1])
if scalevars[2]:
    print 3
else:
    print 1
