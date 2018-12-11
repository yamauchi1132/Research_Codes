import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *
from make_slide import *

make_slide = 1

def calc_r_cg(p, pos_cg):
  dx = p.pos[0] - pos_cg[0]
  dy = p.pos[1] - pos_cg[1]
  dz = p.pos[2] - pos_cg[2]
  r_cg = math.sqrt(dx*dx + dy*dy + dz*dz)
  p.r_cg = r_cg

def calc_r_and_dens(p, pos_cg, vel_cg):
  r = []
  dens = []
  shell = 1e+03
  ration = 2
  inside = 0.
  outside = inside + shell
  mass = 0.

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

  return r, dens

def plot(r_list, dens_list):
  fig = plt.figure()
  
  for i in range(len(r_list)):
    plt.plot(r_list[i], dens_list[i])

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
  for file in range(1, len(args)):
    data = np.loadtxt(args[file])
    p = [Particle() for i in range(len(data))]
    readfile(data, p)

    pos_cg = np.array([0.,0.,0.])
    vel_cg = np.array([0.,0.,0.])
    pos_cg, vel_cg = calc_center_of_gravity(p)

    for i in range(len(data)):
      calc_r_cg(p[i], pos_cg)
    p.sort(key=operator.attrgetter("r_cg"))

    if make_slide == 0:
      r, dens = calc_r_and_dens(p, pos_cg, vel_cg)
    else :
      if file == 1:
        r, dens = calc_r_and_dens_pori(p, pos_cg, vel_cg)
      else:
        r, dens = calc_r_and_dens(p, pos_cg, vel_cg)

    r_list.append(r)
    dens_list.append(dens)

  if make_slide == 0:
    plot(r_list, dens_list)
  else:
    plot2(r_list, dens_list)