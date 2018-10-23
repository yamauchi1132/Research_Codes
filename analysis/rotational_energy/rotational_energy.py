import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

def calc_rotational_energy(p, pos_cg, vel_cg):
  sum_rotational_energy = 0.

  for i in range(len(p)):
    vx = p[i].vel[0] - vel_cg[0]
    vy = p[i].vel[1] - vel_cg[1]
    vz = p[i].vel[2] - vel_cg[2]

    x = p[i].pos[0] - pos_cg[0]
    y = p[i].pos[1] - pos_cg[1]
    z = p[i].pos[2] - pos_cg[2]
    x2 = x * x
    y2 = y * y
    z2 = z * z
    r_1 = math.sqrt(x2 + y2 + z2)

    l_phi = x*vy - y*vx
    omga = l_phi / (r_1 * r_1)

if __name__ == '__main__':
  args = sys.argv

  if(len(args) < 2):
    sys.stderr.write('Error : no input file\n')
    exit()

  for file in range(1, len(args)):
    data = np.loadtxt(args[file])
    p = [Particle() for i in range(len(data))]
    readfile(data, p)

    pos_cg = np.array([0.,0.,0.])
    vel_cg = np.array([0.,0.,0.])
    pos_cg, vel_cg = calc_center_of_gravity(p)

    calc_rotational_energy(p, pos_cg, vel_cg)