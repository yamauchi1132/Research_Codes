import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

def calc_r_and_vel_phi(p):
  r = []
  vel_phi = []
  for i in range(len(p)):
    x2 = p[i].posx * p[i].posx
    y2 = p[i].posy * p[i].posy
    z2 = p[i].posz * p[i].posz
    r_1 = math.sqrt(x2 + y2 + z2)
    r.append(r_1)

  return r, vel_phi

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

  r1, vel1_phi = calc_r_and_vel_phi(p1)
