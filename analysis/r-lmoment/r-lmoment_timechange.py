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
  lmoment = []
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

      l_phi = x*vy - y*vx
      #l_phi = l_phi / (r_1 * r_1)
      lmoment.append(l_phi)

  return r, lmoment

def plot(r_list, lmoment_list, args):
  fig = plt.figure()
  mpl.rcParams['agg.path.chunksize'] = 10000

  for i in range(len(r_list)):
    plt.plot(r_list[i], lmoment_list[i], '.', ms=0.5)
  
  plt.xscale('log')
  plt.yscale('log')

  plt.xlabel(r'$radius\,[cm]$')
  plt.ylabel(r'$L$')

  plt.legend(fontsize=10)
  plt.tick_params(labelsize=10)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  
  plt.show()
  # plt.savefig("r_lmoment.png", dpi=600)
  plt.close()

def plot2(r_list, lmoment_list, args):
  fig = plt.figure()
  mpl.rcParams['agg.path.chunksize'] = 10000

  plt.plot(r_list[0], lmoment_list[0], '.', ms=0.5, label=r'$1.0r_{p}$')
  plt.plot(r_list[1], lmoment_list[1], '.', ms=0.5, label=r'$1.2r_{p}$')

  plt.xscale('log')
  plt.yscale('log')

  plt.xlabel(r'$radius\,[cm]$', fontsize=18)
  plt.ylabel(r'$L$', fontsize=18)

  plt.legend(fontsize=14)
  plt.tick_params(labelsize=18)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  
  # plt.show()
  # plt.savefig("L_0.png", dpi=600)
  plt.savefig("L_400.png", dpi=600)
  plt.close()

if __name__ == '__main__':
  args = sys.argv

  if(len(args) < 2):
    sys.stderr.write('Error : no input file\n')
    exit()

  r_list = []
  lmoment_list = []
  for i in range(1, len(args)):
    data = np.loadtxt(args[i])
    p = [Particle() for i in range(len(data))]
    readfile(data, p)

    pos_cg = np.array([0.,0.,0.])
    vel_cg = np.array([0.,0.,0.])
    pos_cg, vel_cg = calc_center_of_gravity(p)

    r, lmoment = calc_r_and_vel_phi(p, pos_cg, vel_cg)

    r_list.append(r)
    lmoment_list.append(lmoment)

  # plot(r_list, lmoment_list, args)
  plot2(r_list, lmoment_list, args)
