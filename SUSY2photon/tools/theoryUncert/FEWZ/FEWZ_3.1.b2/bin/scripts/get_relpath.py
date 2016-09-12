#!/bin/python

import sys

#import os.path ### older version (<2.5) didn't have os.path.relpath, directly use the source code
#print os.path.relpath(sys.argv[1], sys.argv[2])

from posixpath import curdir, sep, pardir, join, abspath, commonprefix
def myrelpath(path, start=curdir):
    """Return a relative version of a path"""
    if not path:
        raise ValueError("no path specified")
    start_list = abspath(start).split(sep)
    path_list = abspath(path).split(sep)
    # Work out how much of the filepath is shared by start and path.
    i = len(commonprefix([start_list, path_list]))
    rel_list = [pardir] * (len(start_list)-i) + path_list[i:]
    if not rel_list:
        return curdir
    return join(*rel_list)
print myrelpath(sys.argv[1])
