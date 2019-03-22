import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

dirname = "data/sph_t0000.dat"

def calc_r_v_ratio(p1, p2, pos1_cg, vel1_cg, pos2_cg, vel2_cg):
  r_x = []
  v_x = []

  for i in range(len(p1)):
    vx = p1[i].vel[0] - vel1_cg[0]
    vy = p1[i].vel[1] - vel1_cg[1]
    vz = p1[i].vel[2] - vel1_cg[2]
    vx2 = vx * vx
    vy2 = vy * vy
    vz2 = vz * vz
    v2 = vx2 + vy2 + vz2
    ene = 0.5*v2 + p1[i].pot + p1[i].uene

    if(ene < 0):
      r_x.append(p1[i].pos[0]-pos1_cg[0])
      v_x.append(vx)

  return r_x, v_x

if __name__ == '__main__':
  args = sys.argv
  data = np.loadtxt(dirname)

  p = [Particle() for i in range(len(data))]
  readfile(data, p)
  
  p.sort(key=operator.attrgetter("p_id"))
  p1 = [Particle() for i in range(len(data)//2)]
  p2 = [Particle() for i in range(len(data)//2)]
  
  for i in range(len(p)):
    if i < len(p)//2:
      p1[i] = p[i]
    else:
      p2[i-(len(p))] = p[i]

  pos1_cg = np.array([0.,0.,0.])
  vel1_cg = np.array([0.,0.,0.])
  pos2_cg = np.array([0.,0.,0.])
  vel2_cg = np.array([0.,0.,0.])

  pos1_cg, vel1_cg = calc_center_of_gravity(p1)
  pos2_cg, vel2_cg = calc_center_of_gravity(p2)

  r_x, v_x_ratio = calc_r_v_ratio(p1, p2, pos1_cg, vel1_cg, pos2_cg, vel2_cg)

  #plot(r1, omega1, r2, omega2)
  f = open('r_contraction', 'w')
  # f = open('axis_1.3.data', 'w')
  for i in range(len(r_x)):
    f.write("%e %e\n"%(r_x[i], v_x_ratio[i]))
  f.close()