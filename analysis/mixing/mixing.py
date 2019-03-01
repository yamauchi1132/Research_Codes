import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

box_size = 1e+12
dirname = "data/sph_t1100.dat"

def plot_xy(p, pos_cg):
  limit = 0.8e+10
  x1 = []
  y1 = []
  for i in range(len(p)//2):
    x = p[i].pos[0] - pos_cg[0]
    y = p[i].pos[1] - pos_cg[1]
    z = p[i].pos[2] - pos_cg[2]
    if(abs(z) < limit):
      x1.append(x)
      y1.append(y)

  x2 = []
  y2 = []
  for i in range(len(p)//2, len(p)):
    x = p[i].pos[0] - pos_cg[0]
    y = p[i].pos[1] - pos_cg[1]
    z = p[i].pos[2] - pos_cg[2]
    if(abs(z) < limit):
      x2.append(x)
      y2.append(y)

  fsize = 14
  box_size = 1.5e+12
  fig = plt.figure()

  plt.tick_params(labelsize=fsize)
  plt.scatter(x1, y1, s=0.1, color = 'red')
  plt.scatter(x2, y2, s=0.01, color='b')
  plt.rcParams['font.family'] ='sans-serif'
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.xlim(-box_size, box_size)
  plt.ylim(-box_size, box_size)
  plt.rcParams['xtick.major.width'] = 1#x軸主目盛り線の線幅
  plt.rcParams['ytick.major.width'] = 1#y軸主目盛り線の線幅
  plt.locator_params(axis='y',nbins=8)
  plt.locator_params(axis='x',nbins=8)#y軸，6個以内．
  plt.xlabel(r"$x\ [cm]$", fontsize=fsize)
  plt.ylabel(r"$y\ [cm]$", fontsize=fsize)
  plt.tight_layout()
  plt.gca().set_aspect('equal', adjustable='box')

  # plt.savefig("mixing_xy.png", transparent=True, dpi=300, bbox_inches = 'tight', pad_inches = 0)
  plt.show()
  plt.close()

def plot_xz(p, pos_cg):
  limit = 1e+10
  x1 = []
  z1 = []
  for i in range(len(p)//2):
    x = p[i].pos[0] - pos_cg[0]
    y = p[i].pos[1] - pos_cg[1]
    z = p[i].pos[2] - pos_cg[2]
    if(abs(y) < limit):
      x1.append(x)
      z1.append(z)

  x2 = []
  z2 = []
  for i in range(len(p)//2, len(p)):
    x = p[i].pos[0] - pos_cg[0]
    y = p[i].pos[1] - pos_cg[1]
    z = p[i].pos[2] - pos_cg[2]
    if(abs(y) < limit):
      x2.append(x)
      z2.append(z)

  fsize = 14
  box_size = 1.5e+12
  fig = plt.figure()

  plt.tick_params(labelsize=fsize)
  plt.scatter(x1, z1, s=0.01, color = 'red')
  plt.scatter(x2, z2, s=0.01, color='b')
  plt.rcParams['font.family'] ='sans-serif'
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.xlim(-box_size, box_size)
  plt.ylim(-box_size, box_size)
  plt.rcParams['xtick.major.width'] = 1#x軸主目盛り線の線幅
  plt.rcParams['ytick.major.width'] = 1#y軸主目盛り線の線幅
  plt.locator_params(axis='y',nbins=8)
  plt.locator_params(axis='x',nbins=8)#y軸，6個以内．
  plt.xlabel(r"$x\ [cm]$", fontsize=fsize)
  plt.ylabel(r"$z\ [cm]$", fontsize=fsize)
  plt.tight_layout()
  plt.gca().set_aspect('equal', adjustable='box')

  # plt.savefig("mixing_xz.png", transparent=True, dpi=300, bbox_inches = 'tight', pad_inches = 0)
  plt.show()
  plt.close()
if __name__ == '__main__':
  data = np.loadtxt(dirname)

  p = [Particle() for i in range(len(data))]
  readfile(data, p)
  p.sort(key=operator.attrgetter("p_id"))

  pos_cg = np.array([0., 0., 0.])
  vel_cg = np.array([0., 0., 0.])
  pos_cg, vel_cg = calc_center_of_gravity(p)

  # plot_xy(p, pos_cg)
  plot_xz(p, pos_cg)