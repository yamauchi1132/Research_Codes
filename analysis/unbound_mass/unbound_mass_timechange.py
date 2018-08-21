import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

Msun = 1.989e+33

def calc_unbound_mass(p, mass_total, vel_cg):
  mass_unbound = 0.

  for i in range(len(p)):
    vx = p[i].velx - vel_cg[0]
    vy = p[i].vely - vel_cg[1]
    vz = p[i].velz - vel_cg[2]
    vx2 = vx * vx
    vy2 = vy * vy
    vz2 = vz * vz
    v2 = vx2 + vy2 + vz2

    ene = 0.5*v2 + p[i].pot + p[i].uene
    if(ene > 0):
      mass_unbound += p[i].mass

  mass_bound = mass_total - mass_unbound

  return mass_bound, mass_unbound

if __name__ == '__main__':
  args = sys.argv

  if(len(args) < 2):
    sys.stderr.write('Error : no input file\n')
    exit()

  result_bound = []
  result_unbound = []
  for i in range(1,len(args)):
    data = np.loadtxt(args[i])
    p = [Particle() for i in range(len(data))]
    readfile(data, p)
  
    pos_cg = np.array([0.,0.,0.])
    vel_cg = np.array([0.,0.,0.])
    pos_cg, vel_cg = calc_center_of_gravity(p)

    mass_total = 0.
    for i in range(len(p)):
      mass_total += p[i].mass

    m_bound, m_unbound = calc_unbound_mass(p, mass_total, vel_cg)

    result_bound.append(m_bound)
    result_unbound.append(m_unbound)

  sys.stderr.write('\nTotal Mass = %.1lf Msun\n\n' %(mass_total/Msun))

  for i in range(len(result_bound)):
    sys.stderr.write(args[i+1])
    sys.stderr.write('\nBound Mass = %.3lf Msun, Unbound Mass = %.3lf Msun\n\n' %(result_bound[i]/Msun, result_unbound[i]/Msun))
