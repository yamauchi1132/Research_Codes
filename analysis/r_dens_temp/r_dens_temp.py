import sys, os
sys.path.append(os.pardir)
sys.path.append("../")
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

n = 1.5
max_r = 1e+12 #rsun = 695700e+5;

def calc_r_cg(p, pos_cg):
  dx = p.pos[0] - pos_cg[0]
  dy = p.pos[1] - pos_cg[1]
  dz = p.pos[2] - pos_cg[2]
  r_cg = math.sqrt(dx*dx + dy*dy + dz*dz)
  p.r_cg = r_cg

def calc_r_dens_temp(p, pos_cg, vel_cg):
  r = []
  dens = []
  shell = 1e+03
  ration = 2
  inside = 0.
  outside = inside + shell
  mass = 0.
  temp = []
  for i in range(len(p)):
    vx = p[i].vel[0] - vel_cg[0]
    vy = p[i].vel[1] - vel_cg[1]
    vz = p[i].vel[2] - vel_cg[2]
    vx2 = vx * vx
    vy2 = vy * vy
    vz2 = vz * vz
    v2 = vx2 + vy2 + vz2
    ene = 0.5*v2 + p[i].pot + p[i].uene

    if ene < 0:
      if outside >= p[i].r_cg:
        mass += p[i].mass
      else:
        dv = 4 * np.pi * (outside*outside*outside - inside*inside*inside)/ 3.0
        row = mass / dv
        r.append(outside)
        dens.append(row)

        shell *= ration
        inside = outside
        outside += shell
        mass = 0.

  for i in range(len(r)):
    left = 0
    right = len(p)-1
    while(1):
      mid = (left + right) // 2
      if(p[mid].r_cg == r[i]):
        if(dens[i] == 0):
          temp.append(0)
        else:
          pres = dens[i] * p[mid].uene * (1 / n)
          temp_1 = (pres/dens[i]) * (1/8.315e+7)
          temp.append(temp_1)
      if(p[mid].r_cg < r[i]):
        left = mid + 1
      if(p[mid].r_cg > r[i]):
        right = mid - 1
      if(left >= right):
        if(dens[i] == 0):
          temp.append(0)
        else:
          pres = dens[i] * p[mid].uene * (1 / n)
          temp_1 = (pres/dens[i]) * (1/8.315e+7)
          temp.append(temp_1)
        break

  return r, dens, temp

def plot(r, temp):
  fig = plt.figure()

  plt.plot(r, temp)
  plt.xscale('log')
  plt.yscale('log')
  plt.show()
  plt.close()

if __name__ == '__main__':
  args = sys.argv
  data = np.loadtxt("data/sph_t1100.dat")

  p = [Particle() for i in range(len(data))]
  readfile(data, p)

  pos_cg = np.array([0.,0.,0.])
  vel_cg = np.array([0.,0.,0.])
  pos_cg, vel_cg = calc_center_of_gravity(p)

  for i in range(len(data)):
    calc_r_cg(p[i], pos_cg)
  
  p.sort(key=operator.attrgetter("r_cg"))

  r, dens, temp = calc_r_dens_temp(p, pos_cg, vel_cg)

  f = open('r_dens_temp_1.0_1100.data', 'w')
  # f = open('axis_1.3.data', 'w')
  for i in range(len(r)):
    f.write("%e %e %e\n"%(r[i], dens[i], temp[i]))
  f.close()
  # plot(r, temp)
