import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

start, end, step = 100, 1100, 100
dirname = '/lustre/jh180021l/k26002/sph_yamauchi/snap_unbound/snap_unbound_10.0Msun_4.0Rsun_pori1.5_rp1.0R_vinf1.00e+06/sph_t%04d.dat'

'''
start, end, step = 100, 800, 100
dirname = '../../run1/snap_unbound_10.0Msun_4.0Rsun_pori1.5_rp1.0R_vinf1.00e+0\
6/sph_t%04d.dat'
#start, end, step = 600, 1200, 100
#dirname = '../../run2/snap_unbound_10.0Msun_4.0Rsun_pori1.5_rp1.2R_vinf1.00e+06/sph_t%04d.dat'
'''
max_r = 1e+12 #rsun = 695700e+5;

def calc_inertia_tensor(p, pos_cg, vel_cg):
  I = np.zeros((3,3))
  ene_neg_r = []
  #count = 0
  row_sum = 0.
  mass = p[0].mass
  for i in range(len(p)):
    vx = p[i].vel[0] - vel_cg[0]
    vy = p[i].vel[1] - vel_cg[1]
    vz = p[i].vel[2] - vel_cg[2]
    vx2 = vx * vx
    vy2 = vy * vy
    vz2 = vz * vz
    v2 = vx2 + vy2 + vz2
    ene = 0.5*v2 + p[i].pot + p[i].uene

    x = p[i].pos[0] - pos_cg[0]
    y = p[i].pos[1] - pos_cg[1]
    z = p[i].pos[2] - pos_cg[2]
    x2 = x * x
    y2 = y * y
    z2 = z * z
    r2 = x2 + y2 + z2
    r = math.sqrt(r2)

    if(r < max_r and ene < 0):
      #count = count + 1
      #dens = p[i].dens
      #row_sum += dens

      I[0,0] += (r2-x2)
      I[1,1] += (r2-y2)
      I[2,2] += (r2-z2)

      i_01 = -(x * y)
      i_02 = -(x * z)
      i_12 = -(y * z)

      I[0,1] += i_01
      I[1,0] += i_01
      I[0,2] += i_02
      I[2,0] += i_02
      I[1,2] += i_12
      I[2,1] += i_12
  #print(I)
  #calculation of eigenvalue and Diagonal matrix
  #print(I)
  l, P = np.linalg.eig(I)
  #print(l)
  #print(P)
  #DI = np.linalg.inv(P) @ I @ P
  #print(DI)
  #print("\n")
  ######## CHECK ##############
  #print(DI)
  #print(np.diag(l))
  #print(np.dot(np.dot(np.linalg.inv(P),I),P))
  #check PDP^{-1} = I
  #print(np.dot(np.dot(P,np.diag(l)), np.linalg.inv(P))) 
  ##############################
  l_sort = np.sort(l)
  
  return l_sort, mass

def calc_ratio_of_xy_to_z(l, mass):
  x = math.sqrt((l[2]+l[1]-l[0]) / 2)
  y = math.sqrt((l[2]+l[0]-l[1]) / 2)
  z = math.sqrt((l[0]+l[1]-l[2]) / 2)
  #print("%e %e %e"%(x,y,z))
  z_x = z / x
  z_y = z / y

  return z_x, z_y

def plot(time, z_x, z_y):
  fig = plt.figure()

  plt.plot(time, z_x, label='z_x')
  plt.plot(time, z_y, label='z_y')

  #plt.xscale('log')
  #plt.yscale('log')

  '''
  plt.xlabel(r'$radius\,[cm]$', fontsize=18)
  plt.ylabel(r'$L$', fontsize=18)

  plt.legend(fontsize=14)
  plt.tick_params(labelsize=18)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  '''
  plt.show()
  plt.close()
  
if __name__ == '__main__':
  args = sys.argv

  time_list = []
  z_x_list = []
  z_y_list = []

  for time in range(start, end+step, step):
    data = np.loadtxt(dirname % time)
    p = [Particle() for i in range(len(data))]
    readfile(data, p)

    pos_cg = np.array([0.,0.,0.])
    vel_cg = np.array([0.,0.,0.])
    pos_cg, vel_cg = calc_center_of_gravity(p)

    l, mass = calc_inertia_tensor(p, pos_cg, vel_cg)
    time_list.append(time*1e+04)

    z_x, z_y = calc_ratio_of_xy_to_z(l, mass)
    z_x_list.append(z_x)
    z_y_list.append(z_y)

  f = open('axis_1.0.data', 'w')
  # f = open('axis_1.3.data', 'w')
  for i in range(len(time_list)):
    f.write("%e %e %e\n"%(time_list[i]-(start*1e+4), z_x_list[i], z_y_list[i]))
  f.close()
  # plot(time_list, z_x_list, z_y_list)
