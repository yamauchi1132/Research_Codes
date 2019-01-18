import sys, os
sys.path.append(os.pardir)
sys.path.append("../../")
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *
'''
dirname = "data_thesis/snap_unbound_10.0Msun_4.0Rsun_pori1.5_rp1.0R_vinf1.00e+06/sph_t%04d.dat"
start = 100
end = 1100
step = 100
'''


dirname = "data_thesis/snap_unbound_10.0Msun_4.0Rsun_pori1.5_rp1.2R_vinf1.00e+06/sph_t%04d.dat"
start = 600
end = 1600
step = 100

Msun = 1.989e+33

def calc_unbound_mass(p, mass_total, vel_cg):
  mass_unbound = 0.

  for i in range(len(p)):
    vx = p[i].vel[0] - vel_cg[0]
    vy = p[i].vel[1] - vel_cg[1]
    vz = p[i].vel[2] - vel_cg[2]
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

  result_bound = []
  result_unbound = []
  time = []
  for t in range(start,end+step,step):
    data = np.loadtxt(dirname % t)
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
    time.append(t*1e+4)

  # f = open('t_unbound_1.0.data', 'w')
  f = open('t_unbound_1.2.data', 'w')
  for i in range(len(time)):
    f.write("%e %e\n"%(time[i]-(start*1e+4), result_unbound[i]))

  f.close()

  #plot(result_unbound, time, result2_unbound, time2)
  