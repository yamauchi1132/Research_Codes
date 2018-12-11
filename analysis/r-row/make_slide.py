import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt

def calc_r_and_dens_pori(p, pos_cg, vel_cg):
  r = []
  dens = []
  
  for i in range(0, len(p), 100):
    x = p[i].pos[0]
    y = p[i].pos[1]
    z = p[i].pos[2]
    x2 = x * x
    y2 = y * y
    z2 = z * z
    r_1 = math.sqrt(x2 + y2 + z2)
    r.append(r_1)
    dens.append(p[i].dens)
    
    if(i+100 >= len(p)):
      r.append(r_1)
      dens.append(1e-10)
      '''
      dr = r_1 - math.sqrt(p[i-100].pos[0]**2+p[i-100].pos[1]**2+p[i-100].pos[2]**2)
      drow = p[i].dens - p[i-100].dens
      a = drow / dr
      row_next = 1e-9
      drow_row_next = row_next - p[i].dens
      r_next = (drow_row_next / a) + r_1
      r.append(1.5e+12)
      dens.append(row_next)
      '''

  return r, dens

def plot2(r_list, dens_list):
  fig = plt.figure()
  # mpl.rcParams['agg.path.chunksize'] = 10000

  plt.plot(r_list[0], dens_list[0])
  plt.plot(r_list[1], dens_list[1])

  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(1e+11)

  plt.xlabel(r'$r\,[cm]$', fontsize=18)
  plt.ylabel(r'$\rho\,[g\,cm^{-3}]$', fontsize=18)

  plt.legend(fontsize=12, loc='upper right')
  plt.tick_params(labelsize=18)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  
  # plt.show()
  plt.savefig("r_dens.png", dpi=600)
  plt.close()
