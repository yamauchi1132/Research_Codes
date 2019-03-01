import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

def calc_r_and_vel_phi(p, pos_cg, vel_cg):
  r = []
  omega = []
  for i in range(len(p)):
    vx = p[i].vel[0] - vel_cg[0]
    vy = p[i].vel[1] - vel_cg[1]
    vz = p[i].vel[2] - vel_cg[2]
    vx2 = vx * vx
    vy2 = vy * vy
    vz2 = vz * vz
    v2 = vx2 + vy2 + vz2
    ene = 0.5*v2 + p[i].pot + p[i].uene

    if(ene < 0):
      x = p[i].pos[0] - pos_cg[0]
      y = p[i].pos[1] - pos_cg[1]
      z = p[i].pos[2] - pos_cg[2]
      x2 = x * x
      y2 = y * y
      z2 = z * z
      r_1 = math.sqrt(x2 + y2 + z2)
      r.append(r_1)

      v = x*vy - y*vx
      w = v / (x2+y2)
      omega.append(w)

  return r, omega

def plot(r1, omega1):
  fig = plt.figure()
  mpl.rcParams['agg.path.chunksize'] = 10000
  plt.plot(r1, omega1, '.', ms=0.5, label=1.0)

  plt.xscale('log')
  plt.yscale('log')

  plt.xlabel(r'$r\,[cm]$', fontsize=18)
  plt.ylabel(r'$\omega\,[rad\,s^{-1}]$', fontsize=18)

  plt.legend(fontsize=18)
  plt.tick_params(labelsize=18)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  
  plt.savefig("r_omega.pdf", transparent=True, dpi=300, bbox_inches = 'tight', pad_inches = 0)
  plt.show()
  plt.close()

def plot2(r1, omega1, r2, omega2):
  fig = plt.figure()
  mpl.rcParams['agg.path.chunksize'] = 10000
  plt.plot(r1, omega1, '.', ms=0.5, label=1.0)
  plt.plot(r2, omega2, '.', ms=0.5, label=1.2)

  plt.xscale('log')
  plt.yscale('log')

  plt.xlabel(r'$r\,[cm]$', fontsize=18)
  plt.ylabel(r'$\omega\,[rad\,s^{-1}]$', fontsize=18)

  plt.legend(fontsize=18)
  plt.tick_params(labelsize=18)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  
  plt.show()
  # plt.savefig("r_omega.png", dpi=600)
  plt.close()

if __name__ == '__main__':
  args = sys.argv

  if(len(args) < 3):
    sys.stderr.write('Error : no input file\n')
    exit()

  data1 = np.loadtxt(args[1])
  #data2 = np.loadtxt(args[2])

  p1 = [Particle() for i in range(len(data1))]
  #p2 = [Particle() for i in range(len(data2))]

  readfile(data1, p1)
  #readfile(data2, p2)

  pos1_cg = np.array([0.,0.,0.])
  vel1_cg = np.array([0.,0.,0.])
  #pos2_cg = np.array([0.,0.,0.])
  #vel2_cg = np.array([0.,0.,0.])

  pos1_cg, vel1_cg = calc_center_of_gravity(p1)
  #pos2_cg, vel2_cg = calc_center_of_gravity(p2)

  r1, omega1 = calc_r_and_vel_phi(p1, pos1_cg, vel1_cg)
  #r2, omega2 = calc_r_and_vel_phi(p2, pos2_cg, vel2_cg)

  plot(r1, omega1)
  #plot2(r1, omega1, r2, omega2)
