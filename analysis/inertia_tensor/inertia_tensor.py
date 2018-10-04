import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

start, end, step = 600, 800, 100
def calc_inertia_tensor(p, pos_cg, vel_cg):
  I = np.zeros((3,3))

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

    if(ene < 0):
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

  #calculation of eigenvalue and Diagonal matrix
  #print(I)
  l, P = np.linalg.eig(I)
  DI = np.linalg.inv(P) @ I @ P
  ######## CHECK ##############
  #print(DI)
  #print(np.diag(l))
  #print(np.dot(np.dot(np.linalg.inv(P),I),P))
  #check PDP^{-1} = I
  #print(np.dot(np.dot(P,np.diag(l)), np.linalg.inv(P))) 
  ##############################
  l_sort = np.sort(l)
  return l_sort

def plot(time, l1, l2, l3):
  fig = plt.figure()

  plt.plot(time, l1, label='l1')
  plt.plot(time, l2, label='l2')
  plt.plot(time, l3, label='l3')

  #plt.xscale('log')
  plt.yscale('log')

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

  if(len(args) < 2):
    sys.stderr.write('Error : no input file\n')
    exit()

  time_list = []
  l1_list = []
  l2_list = []
  l3_list = []
  file_num = 1
  for time in range(start, end+step, step):
    data = np.loadtxt(args[file_num])
    file_num += 1
    p = [Particle() for i in range(len(data))]
    readfile(data, p)

    pos_cg = np.array([0.,0.,0.])
    vel_cg = np.array([0.,0.,0.])
    pos_cg, vel_cg = calc_center_of_gravity(p)

    l = calc_inertia_tensor(p, pos_cg, vel_cg)
    time_list.append(time)
    l1_list.append(l[0])
    l2_list.append(l[1])
    l3_list.append(l[2])
    
  plot(time_list, l1_list, l2_list, l3_list)