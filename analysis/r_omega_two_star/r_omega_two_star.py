import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *
import matplotlib.ticker as ptick ##これが必要！

dirname = "data/sph_t0000.dat"

def calc_r_omega(p1, p2, pos1_cg, vel1_cg, pos2_cg, vel2_cg):
  r1 = []
  r2 = []
  omega1 = []
  omega2 = []
  for i in range(0, len(p1), 100):
    vx = p1[i].vel[0] - vel1_cg[0]
    vy = p1[i].vel[1] - vel1_cg[1]
    vz = p1[i].vel[2] - vel1_cg[2]
    vx2 = vx * vx
    vy2 = vy * vy
    vz2 = vz * vz
    v2 = vx2 + vy2 + vz2
    ene = 0.5*v2 + p1[i].pot + p1[i].uene

    if(ene < 0):
      x = p1[i].pos[0] - pos1_cg[0]
      y = p1[i].pos[1] - pos1_cg[1]
      x2 = x * x
      y2 = y * y
      r_2 = x2 + y2

      l_phi = x*vy - y*vx
      omega = l_phi / (r_2)

      r1.append(math.sqrt(r_2))
      omega1.append(omega)

  for i in range(0, len(p2), 100):
    vx = p2[i].vel[0] - vel2_cg[0]
    vy = p2[i].vel[1] - vel2_cg[1]
    vz = p2[i].vel[2] - vel2_cg[2]
    vx2 = vx * vx
    vy2 = vy * vy
    vz2 = vz * vz
    v2 = vx2 + vy2 + vz2
    ene = 0.5*v2 + p2[i].pot + p2[i].uene

    if(ene < 0):
      x = p2[i].pos[0] - pos2_cg[0]
      y = p2[i].pos[1] - pos2_cg[1]
      x2 = x * x
      y2 = y * y
      r_2 = x2 + y2

      l_phi = x*vy - y*vx
      omega = l_phi / (r_2)

      r2.append(math.sqrt(r_2))
      omega2.append(omega)

  return r1, omega1, r2, omega2

def plot(r1, omega1, r2, omega2):
  fsize = 14
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  ax.plot(r1, omega1, linewidth=0, marker='o', color='red', ms=1)
  ax.plot(r2, omega2, linewidth=0, marker='o', color='b', ms=1)

  # plt.xscale('log')
  # plt.yscale('log')

  plt.xlabel(r'$r\,[cm]$', fontsize=fsize)
  plt.ylabel(r'$\omega\,[rad\,s^{-1}]$', fontsize=fsize)

  plt.legend(fontsize=fsize)
  plt.tick_params(labelsize=fsize)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  ax.xaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True)) 

  # plt.savefig("r_omega.pdf", transparent=True, dpi=300, bbox_inches = 'tight', pad_inches = 0)
  plt.show()
  plt.close()

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

  r1, omega1, r2, omega2 = calc_r_omega(p1, p2, pos1_cg, vel1_cg, pos2_cg, vel2_cg)

  plot(r1, omega1, r2, omega2)

