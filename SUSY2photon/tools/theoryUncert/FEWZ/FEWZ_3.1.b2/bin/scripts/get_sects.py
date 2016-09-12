"""
Call python getsects() function in "defs.py" and print the result on screen
Usage: python $0 <inputfile> <boson>
"""

import sys
from defs import getsects

print getsects(sys.argv[1], sys.argv[2].upper())

