import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
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
    ene = 0.5*v2 + p[i].pot + p[i].uene

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

def plot(r_list, dens_list, args):
  fig = plt.figure()
  mpl.rcParams['agg.path.chunksize'] = 10000

  for i in range(len(r_list)):
    plt.plot(r_list[i], dens_list[i], '.', ms=0.5, label=args[i+1])

  plt.xscale('log')
  plt.yscale('log')

  plt.xlabel(r'$r\,[cm]$')
  plt.ylabel(r'$\rho\,[g\,cm^{-3}]$')

  plt.legend(fontsize=10)
  plt.tick_params(labelsize=10)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  
  plt.show()
  # plt.savefig("r_dens.png", dpi=600)
  plt.close()

if __name__ == '__main__':
  args = sys.argv

  if(len(args) < 2):
    sys.stderr.write('Error : no input file\n')
    exit()

  r_list = []
  dens_list = []
  for i in range(1, len(args)):
    data = np.loadtxt(args[i])
    p = [Particle() for i in range(len(data))]
    readfile(data, p)

    pos_cg = np.array([0.,0.,0.])
    vel_cg = np.array([0.,0.,0.])
    pos_cg, vel_cg = calc_center_of_gravity(p)

    r, dens = calc_r_and_dens(p, pos_cg, vel_cg)  

    r_list.append(r)
    dens_list.append(dens)

  plot(r_list, dens_list, args)
