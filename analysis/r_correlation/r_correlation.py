import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *
import matplotlib.ticker as ptick ##これが必要！

data = np.loadtxt("data/sph_t0000.dat")
data1 = np.loadtxt("data/sph_t1100.dat")

def plot(r_ini, r_end):
  r_ini_1 = []
  r_ini_2 = []

  r_end_1 = []
  r_end_2 = []

  for i in range(0, len(r_ini)//2, 100):
    r_ini_1.append(r_ini[i])
    r_end_1.append(r_end[i])

  for i in range(len(r_ini)//2, len(r_ini), 100):
    r_ini_2.append(r_ini[i])
    r_end_2.append(r_end[i])

  fsize =14
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  ax.plot(r_ini_1, r_end_1, linewidth=0, marker='o', color='red', ms=1, label=r"Star 1")
  ax.plot(r_ini_2, r_end_2, linewidth=0, marker='o', color='b', ms=1, alpha=0.3, label=r"Star 2")
  plt.legend(fontsize=fsize)

  plt.tick_params(labelsize=fsize)
  plt.rcParams['font.family'] ='sans-serif'
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  ax.xaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True)) 

  # plt.xscale("log")
  plt.yscale("log")
  plt.ylim(1e+12, 1e+13)

  plt.xlabel(r"$r_{\rm{ini}}\ [cm]$", fontsize=fsize)
  plt.ylabel(r"$r_{\rm{remnant}}\ [cm]$", fontsize=fsize)
  plt.tight_layout()
  plt.savefig("mixing.pdf", transparent=True, dpi=300, bbox_inches = 'tight', pad_inches = 0)
  plt.show()
  plt.close()

def calc_r(p):
  r = []
  for i in range(len(p)):
    r2 = p[i].pos[0]*p[i].pos[0] + p[i].pos[1]*p[i].pos[1] + p[i].pos[2]*p[i].pos[2]
    r.append(math.sqrt(r2))

  return r

if __name__ == '__main__':
  p_ini = [Particle() for i in range(len(data))]
  readfile(data, p_ini)

  p_end = [Particle() for i in range(len(data1))]
  readfile(data1, p_end)

  p_ini.sort(key=operator.attrgetter("p_id"))
  p_end.sort(key=operator.attrgetter("p_id"))

  r_ini = calc_r(p_ini)
  r_end = calc_r(p_end)

  f = open('mixing_1.0_1100.data', 'w')
  # f = open('axis_1.3.data', 'w')
  for i in range(0, len(r_ini)//2, 100):
    f.write("%e %e %e %e\n"%(r_ini[i], r_end[i], r_ini[i+(len(r_ini)//2)], r_end[i+(len(r_ini)//2)]))
  f.close()

  # plot(r_ini, r_end)
