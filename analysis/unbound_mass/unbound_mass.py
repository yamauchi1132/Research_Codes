import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
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

  if(len(args) < 3):
    sys.stderr.write('Error : no input file\n')
    exit()

  data1 = np.loadtxt(args[1])
  data2 = np.loadtxt(args[2])

  p1 = [Particle() for i in range(len(data1))]
  p2 = [Particle() for i in range(len(data2))]

  readfile(data1, p1)
  readfile(data2, p2)

  pos1_cg = np.array([0.,0.,0.])
  vel1_cg = np.array([0.,0.,0.])
  pos2_cg = np.array([0.,0.,0.])
  vel2_cg = np.array([0.,0.,0.])

  pos1_cg, vel1_cg = calc_center_of_gravity(p1)
  pos2_cg, vel2_cg = calc_center_of_gravity(p2)

  mass_total = 0.
  for i in range(len(p1)):
    mass_total += p1[i].mass
  
  m1_bound, m1_unbound = calc_unbound_mass(p1, mass_total, vel1_cg)
  m2_bound, m2_unbound = calc_unbound_mass(p2, mass_total, vel2_cg)

  sys.stderr.write('\nTotal Mass = %.1lf Msun\n\n' %(mass_total/Msun))

  sys.stderr.write('Rp = 1.0\n')
  sys.stderr.write('Bound Mass = %.3lf Msun, Unbound Mass = %.3lf Msun\n\n' %(m1_bound/Msun, m1_unbound/Msun))

  sys.stderr.write('Rp = 1.2\n')
  sys.stderr.write('Bound Mass = %.3lf Msun, Unbound Mass = %.3lf Msun\n\n' %(m2_bound/Msun, m2_unbound/Msun))