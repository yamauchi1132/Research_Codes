#Code revies is hear : https://qiita.com/NatsukiLab/items/83e6d31aa4a49acc528c
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

########################################
#if 0 open show window, 
#if 1 save png file, 
show_or_save = 1
########################################

########################################
#set animation axis x-y or x-z
#if 1, set this axis. if 0, not set.
x_y = 1 #0 or 1
x_z = 0 #0 or 1
#######################################

norm = 1e+10
box_size = 8.0
Nnode = 1

def plot_dens(x, y, dens):
  fig = plt.figure()
  plt.rcParams["font.size"] = 18
  plt.scatter(x,y, c=dens, s=0.05, cmap="jet")
  clb = plt.colorbar()
  clb.set_label(label=r'$\rho\,[g\,cm^{-3}]$')
  plt.clim(0,10)
  plt.xlim(-box_size, box_size)
  plt.ylim(-box_size, box_size)
  plt.xlabel(r'$x\,[10^{10}cm]$')
  plt.ylabel(r'$y\,[10^{10}cm]$')
  plt.gca().set_aspect('equal', adjustable='box')
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()

  if(show_or_save == 0):
    plt.show()

  if(show_or_save == 1):
    plt.savefig("final.png", dpi=1000)
    
  plt.close()

if __name__ == '__main__':
  if(x_y==x_z):
    print("Error : choose axis x-y or x-z")
    exit(0)
  if(show_or_save == 0 or show_or_save == 1):
    for n in range(Nnode):
      tag = []
      x = []
      y = []
      z = []
      dens = []
      # data = np.loadtxt("dat/sph_t%04d.dat" % (n));
      data = np.loadtxt("final.dat");


      tag.extend(data[:,0])
      x.extend(data[:,3])
      y.extend(data[:,4])
      z.extend(data[:,5])
      dens.extend(data[:,15])

      for i in range(len(x)):
        x[i] = x[i]/norm
        y[i] = y[i]/norm
        z[i] = z[i]/norm
        
      x_tag = [[],[]]
      y_tag = [[],[]]
      z_tag = [[],[]]

      for i in range(len(x)):
        if tag[i] // (len(x)/2) == 0:
          x_tag[0].append(x[i])
          y_tag[0].append(y[i])
          z_tag[0].append(z[i])
        else:
          x_tag[1].append(x[i])
          y_tag[1].append(y[i])
          z_tag[1].append(z[i])

      # plot_3D(x_tag,y_tag,z_tag)
      # plot_2D(x_tag,y_tag,z_tag)
      plot_dens(x, y, dens)
      del x, y, z, dens, x_tag, y_tag, z_tag, data