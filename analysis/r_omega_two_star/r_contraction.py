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
  r_y = []
  v_y = []
  r_z = []
  v_z = []
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
      r_y.append(p1[i].pos[1]-pos1_cg[1])
      v_y.append(vy)
      r_z.append(p1[i].pos[2]-pos1_cg[2])
      v_z.append(vz)

  for i in range(len(p2)):
    vx = p2[i].vel[0] - vel2_cg[0]
    vy = p2[i].vel[1] - vel2_cg[1]
    vz = p2[i].vel[2] - vel2_cg[2]
    vx2 = vx * vx
    vy2 = vy * vy
    vz2 = vz * vz
    v2 = vx2 + vy2 + vz2
    ene = 0.5*v2 + p2[i].pot + p2[i].uene

    if(ene < 0):
      r_x.append(p2[i].pos[0]-pos2_cg[0])
      v_x.append(vx)
      r_y.append(p2[i].pos[1]-pos2_cg[1])
      v_y.append(vy)
      r_z.append(p2[i].pos[2]-pos2_cg[2])
      v_z.append(vz)

  return r_x, v_x, r_y, v_y, r_z, v_z

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

  r_x, v_x, r_y, v_y, r_z, v_z = calc_r_v_ratio(p1, p2, pos1_cg, vel1_cg, pos2_cg, vel2_cg)

  #plot(r1, omega1, r2, omega2)
  fx = open('r_contraction_x.data', 'w')
  fy = open('r_contraction_y.data', 'w')
  fz = open('r_contraction_z.data', 'w')
  # f = open('axis_1.3.data', 'w')
  for i in range(len(r_x)//2):
    fx.write("%e %e %e %e\n"%(r_x[i], v_x[i], r_x[i+(len(r_x)//2)], v_x[i+(len(r_x)//2)]))
  for i in range(len(r_y)//2):
    fy.write("%e %e %e %e\n"%(r_y[i], v_y[i], r_y[i+(len(r_y)//2)], v_y[i+(len(r_y)//2)]))
  for i in range(len(r_z)//2):
    fz.write("%e %e %e %e\n"%(r_z[i], v_z[i], r_z[i+(len(r_z)//2)], v_z[i+(len(r_z)//2)]))

  fx.close()  
  fy.close()
  fz.close()