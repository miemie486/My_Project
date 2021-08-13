#! /usr/bin/python
"""
put the last column of a data file in the first
"""

import sys
import os
import shutil
import fnmatch

def rewind(fname, pos):
  "for a given file, place the last column in the front"
  try:
    fr = open(fname, 'rt')
    lines = fr.readlines()
    fr.close()
    shutil.copy(fname, fname+'.BAK~')
    fw = open(fname, 'wt')
    for line in lines:
      if not '#' in line:
        lst = line.rsplit(None, pos)
        if len(lst) == pos + 1:
          Lmbd = lst.pop(1)
          line = Lmbd + '  ' + '  '.join(lst) + '\n'
      #~ print line
      fw.write(line)
    fw.close()
  except IOError:
    pass
    print '1'

if len(sys.argv) > 2:
  try:
    pos = int(sys.argv[1])
  except ValueError:
    pos = -1
    print "counting from the last, which column?"
  if pos != -1:
    for fname in sys.argv[2:]:
      rewind(fname, pos)
