#! /usr/bin/python
"""
Compute dibaryonic energy-dependent 1s0 from C0's that are fitted to two
energies:
C0(Ecm) = y/(Ecm + Delta) - CgA
"""

import math, time

PC_mN  = 939.0
PC_fpi = 92.4
PC_gA  = 1.29
PI_NE  = math.pi
CgA    = PC_gA*PC_gA*PC_mN/(8.0*PI_NE*PI_NE*PC_fpi*PC_fpi)

def extract_C0_Lambda(fname):
  f = open(fname, "rt")
  CLmbd = []
  for line in f:
    if line[0] != '#':
      li = line.split()
      try:
        CLmbd.append([float(li[0]), float(li[1])])
      except ValueError:
        pass
  f.close()
  return CLmbd

tlab1  = 2.0
tlab2  = 15.0
fname1 = "C0_1s0_02MeV.arx"
fname2 = "C0_1s0_15MeV.arx"

# tlab1  = 5.0
# tlab2  = 15.0
# fname1 = "C0_1s0_05MeV.arx"
# fname2 = "C0_1s0_15MeV.arx"

# tlab1  = 0.1
# tlab2  = 0.2
# fname1 = "C0_1s0_0.1MeV.arx"
# fname2 = "C0_1s0_0.2MeV.arx"

# tlab1  = 0.1
# tlab2  = 2.0
# fname1 = "C0_1s0_0.1MeV.arx"
# fname2 = "C0_1s0_02MeV_68deg.arx"

# tlab1  = 0.1
# tlab2  = 5.0
# fname1 = "C0_1s0_0.1MeV.arx"
# fname2 = "C0_1s0_05MeV_66deg.arx"

CLmbd1 = extract_C0_Lambda(fname1)
CLmbd2 = extract_C0_Lambda(fname2)
k1     = math.sqrt(PC_mN*tlab1*0.5)
k2     = math.sqrt(PC_mN*tlab2*0.5)
#~ ecm1   = 2.0*(math.sqrt(PC_mN**2+k1**2) - PC_mN)
#~ ecm2   = 2.0*(math.sqrt(PC_mN**2+k2**2) - PC_mN)
ecm1 = k1**2/PC_mN
ecm2 = k2**2/PC_mN
DyLmbd = []
ii = 0
for elem1 in CLmbd1:
  elem2 = CLmbd2[ii]
  C1    = elem1[1] + CgA
  C2    = elem2[1] + CgA
  ii += 1
  if math.fabs((elem1[0] - elem2[0])/elem1[0]) < 1.0E-4:
    Delta = -(C1*ecm1 - C2*ecm2)/(C1 - C2)
    y     = C1*(ecm1 + Delta)
    DyLmbd.append([elem1[0], Delta, y])

comments = '# '+time.ctime()+'\n'
comments += 'LBW LO parameters'
comments += "# Lambda  Delta  y\n"
print comments+'#   '+fname1+',   '+fname2
print '\n'.join(["%8.2f    %20.17E   %20.17E" % (e[0], e[1], e[2]) for e in DyLmbd])

