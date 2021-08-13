#! /usr/bin/python
"""
Assuming the first column is Tlab, calculate the kcm and put it in the new file
as the first column. The old file is named *.OLD_inskcm~
"""
import math, sys, shutil

def extract_Tlab_and_rest(fname):
  " Extrac the first element as a number and the rest as a string "
  try:
    f = open(fname, "rt")
    lines = f.readlines()
    f.close()
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
  rslt = []
  comments = []
  for line in lines:
    if not '#' in line:
      li = line.split(None, 1)
      if len(li) == 2:
        try:
          rslt.append([float(li[0]), li[1]])
        except ValueError:
          pass
      else:
        comments.append(line)
    else:
      comments.append(line)
  return (rslt, comments)

PC_mN  = 939.0

def kcm_from_Tlab(Tlab):
  " Calculate the corresponding kcm for the given Tlab "
  if Tlab < 0.001:
    kcm = 0.0
  else:
    kcm = math.sqrt(0.5*(PC_mN*Tlab))
  return kcm

def write_new_file(fname):
  " Read the old file and create the new file "
  try:
    f = open(fname, "rt")
    firstline = f.readline()
    f.close()
    if firstline == '# kcm inserted\n':
      print fname + ' appears to have been modified.'
      return
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    return
  source_arrays, comments = extract_Tlab_and_rest(fname)
  try:
    shutil.copy(fname, fname+'.OLD_inskcm~')
  except shutil.Error as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
  new_lines = ['# kcm inserted\n'] + comments
  for array in source_arrays:
    new_lines.append("%10.3f" % kcm_from_Tlab(array[0]) + '    ' + "%10.3f" % array[0] + '      ' + array[1])
  try:
    f = open(fname, "wt")
    f.writelines(new_lines)
    f.close()
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)

# print sys.argv[0], sys.argv[1]
for fname in sys.argv[1:]:
  write_new_file(fname)

