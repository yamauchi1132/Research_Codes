import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
from common import *

def calc_r_and_dens(p, pos_cg, vel_cg):
  r = []
  dens = []
  for i in range(len(p)):
    vx = p[i].velx - vel_cg[0]
    vy = p[i].vely - vel_cg[1]
    vz = p[i].velz - vel_cg[2]
    vx2 = vx * vx
    vy2 = vy * vy
    vz2 = vz * vz
    v2 = vx2 + vy2 + vz2
    ene = 0.5*v2 + 0.5*p[i].pot + p[i].uene

    if(ene < 0):
      x = p[i].posx - pos_cg[0]
      y = p[i].posy - pos_cg[1]
      z = p[i].posz - pos_cg[2]
      x2 = x * x
      y2 = y * y
      z2 = z * z
      r_1 = math.sqrt(x2 + y2 + z2)
      r.append(r_1)
      dens.append(p[i].dens)

  return r, dens

def plot(r1, dens1, r2, dens2):
  fig = plt.figure()
  mpl.rcParams['agg.path.chunksize'] = 10000
  plt.plot(r1, dens1, '.', ms=0.5, label=1.0)
  plt.plot(r2, dens2, '.', ms=0.5, label=1.2)

  plt.xscale('log')
  plt.yscale('log')

  plt.xlabel(r'$r\,[cm]$', fontsize=18)
  plt.ylabel(r'$\rho\,[g\,cm^{-3}]$', fontsize=18)

  plt.legend(fontsize=18)
  plt.tick_params(labelsize=18)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  
  plt.show()
  # plt.savefig("r_dens.png", dpi=600)
  plt.close()

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
 
  r1, dens1 = calc_r_and_dens(p1, pos1_cg, vel1_cg)
  r2, dens2 = calc_r_and_dens(p2, pos2_cg, vel2_cg)

  plot(r1, dens1, r2, dens2)
