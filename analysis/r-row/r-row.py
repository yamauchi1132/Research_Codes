import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

def calc_r_and_dens(p):
  r = []
  dens = []
  for i in range(len(p)):
    x2 = p[i].posx * p[i].posx
    y2 = p[i].posy * p[i].posy
    z2 = p[i].posz * p[i].posz
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
  # plt.savefig("dens_r.png", dpi=600)
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

  r1, dens1 = calc_r_and_dens(p1)
  r2, dens2 = calc_r_and_dens(p2)

  plot(r1, dens1, r2, dens2)
