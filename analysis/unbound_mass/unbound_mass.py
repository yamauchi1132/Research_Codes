import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

def calc_unbound_mass(p, mass_total):
  mass_unbound = 0.

  for i in range(len(p)):
    vx2 = p[i].velx * p[i].velx
    vy2 = p[i].vely * p[i].vely
    vz2 = p[i].velz * p[i].velz
    v2 = vx2 + vy2 + vz2

    ene = 0.5*v2 + p[i].pot + p[i].uene
    if(ene > 0):
      mass_unbound += p[i].mass

  mass_bound = mass_total - mass_unbound

  return mass_bound, mass_unbound

if __name__ == '__main__':
  args = sys.argv

  if(len(args) < 3):
    sys.stderr.write('Error : no input file\n')
    exit()

  data1 = np.loadtxt(args[1])
  data2 = np.loadtxt(args[2])

  p1 = [Particle() for i in range(len(data1))]
  p2 = [Particle() for i in range(len(data2))]

  readfile(data1, p1)
  readfile(data2, p2)

  mass_total = 0.
  for i in range(len(p1)):
    mass_total += p1[i].mass
  
  m1_bound, m1_unbound = calc_unbound_mass(p1, mass_total)
  m2_bound, m2_unbound = calc_unbound_mass(p2, mass_total)

  sys.stderr.write('\nTotal Mass = %e\n\n' %mass_total)

  sys.stderr.write('Rp = 1.0\n')
  sys.stderr.write('Bound Mass = %e, Unbound Mass = %e\n\n' %(m1_bound, m1_unbound))

  sys.stderr.write('Rp = 1.2\n')
  sys.stderr.write('Bound Mass = %e, Unbound Mass = %e\n\n' %(m2_bound, m2_unbound))