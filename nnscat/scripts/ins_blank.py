#!  /usr/bin/python

"""
If there is a comment line immediately above 'subroutine', this program insert
a blank in between, to separate the comments and the definition of the
subroutine
"""

import sys
import os
import shutil

def insert_blank(fname):
  fr = open(fname, "rt")
  lines = fr.readlines()
  fr.close()
  shutil.copy(fname, fname+'.BAK~')
  fw = open(fname, 'wt')
  last_is_comment = False
  for line in lines:
    if last_is_comment:
      if line.lstrip().startswith('subroutine') or line.lstrip().startswith('function'):
        fw.write('\n')
    if line.lstrip().startswith('!'):
      last_is_comment = True
    else:
      last_is_comment = False
    fw.write(line)
  fw.close()


if len(sys.argv) > 1:
  for fname in sys.argv[1:]:
    insert_blank(fname)
