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
#if 2 save animation into mp4 format
show_or_save = 2
########################################

########################################
#set animation axis x-y or x-z
#if 1, set this axis. if 0, not set.
x_y = 0 #0 or 1
x_z = 1 #0 or 1
#######################################

box_size = 6.0
Nnode = 501

def plot_3D(x_tag, y_tag, z_tag):
  fig = plt.figure()
  plt.rcParams["font.size"] = 15
  clr = ["orange", "gray"]
  ax = Axes3D(fig)
  for tag in range(2):
    ax.scatter(x_tag[tag], y_tag[tag], z_tag[tag], s=1.0, c=clr[tag], edgecolor=clr[tag], alpha=0.1)
  ax.set_aspect('equal')
  ax.set_xlim3d(-box_size, box_size)
  ax.set_ylim3d(-box_size, box_size)
  ax.set_zlim3d(-box_size, box_size)
  ax.view_init(9, 45)

  if(show_or_save == 0):
    plt.show()

  if(show_or_save == 1):
    plt.savefig("./png_data3D/%04d.png" % (n), dpi=600)

  plt.close()

def plot_2D(x_tag, y_tag, z_tag):
  fig = plt.figure()
  plt.rcParams["font.size"] = 15
  clr = ["orange", "gray"]
  for tag in range(2):
    plt.scatter(x_tag[tag],y_tag[tag], c=clr[tag], s=1, edgecolor=clr[tag], alpha=0.1)
  # for tag in range(2):
    # plt.scatter(x_tag[tag],z_tag[tag], c=clr[tag], s=1, edgecolor=clr[tag], alpha=0.1)
  plt.xlim(-box_size, box_size)
  plt.ylim(-box_size, box_size)
  plt.xlabel(r'$x\,[10^{11}cm]$')
  plt.ylabel(r'$y\,[10^{11}cm]$')
  plt.gca().set_aspect('equal', adjustable='box')

  if(show_or_save == 0):
    plt.show()

  if(show_or_save == 1):
    plt.savefig("./png_data2D/%04d.png" % (n), dpi=600)
    # plt.savefig("./png_data2D_2/%04d.png" % (n), dpi=600)

  plt.close()

def plot_dens(x, y, dens):
  fig = plt.figure()
  plt.rcParams["font.size"] = 14
  plt.scatter(x,y, c=dens, s=0.001, cmap="jet")
  clb = plt.colorbar()
  clb.set_label(label=r'$\rho\,[g\,cm^{-3}]$')
  plt.clim(0,10)
  plt.xlim(-box_size, box_size)
  plt.ylim(-box_size, box_size)
  plt.xlabel(r'$x\,[10^{11}cm]$')
  plt.ylabel(r'$y\,[10^{11}cm]$')
  plt.gca().set_aspect('equal', adjustable='box')
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()

  if(show_or_save == 0):
    plt.show()

  if(show_or_save == 1):
    plt.savefig("./png_data2D_dens/%04d.png" % (n), dpi=1000)
    
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
      data = np.loadtxt("dat/sph_t%04d.dat" % (n));

      tag.extend(data[:,0])
      x.extend(data[:,3])
      y.extend(data[:,4])
      z.extend(data[:,5])
      dens.extend(data[:,15])

      for i in range(len(x)):
        x[i] = x[i]/1e+11
        y[i] = y[i]/1e+11
        z[i] = z[i]/1e+11

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

  if(show_or_save == 2):
    fig = plt.figure()
    plt.rcParams["font.size"] = 14
    ims = []
    for n in range(Nnode):
      tag = []
      x = []
      y = []
      z = []
      dens = []
      data = np.loadtxt("dat/sph_t%04d.dat" % (n));

      tag.extend(data[:,0])
      x.extend(data[:,3])
      y.extend(data[:,4])
      z.extend(data[:,5])
      dens.extend(data[:,15])

      for i in range(len(x)):
        x[i] = x[i]/1e+11
        y[i] = y[i]/1e+11
        z[i] = z[i]/1e+11

      if(x_y == 1):
        im = plt.scatter(x, y, c=dens, s=0.001, cmap="jet")
      if(x_z == 1):
        im = plt.scatter(x, z, c=dens, s=0.001, cmap="jet")
      if(n==0):
        clb = plt.colorbar()
        clb.set_label(label=r'$\rho\,[g\,cm^{-3}]$')
        plt.clim(0,10)
      plt.xlim(-box_size, box_size)
      plt.ylim(-box_size, box_size)
      if(x_y == 1):
        plt.xlabel(r'$x\,[10^{11}cm]$')
        plt.ylabel(r'$y\,[10^{11}cm]$')
      if(x_z == 1):
        plt.xlabel(r'$x\,[10^{11}cm]$')
        plt.ylabel(r'$z\,[10^{11}cm]$')
      plt.gca().set_aspect('equal', adjustable='box')
      mpl.rcParams['axes.xmargin'] = 0
      mpl.rcParams['axes.ymargin'] = 0
      plt.tight_layout()
      ims.append([im])

    ani = animation.ArtistAnimation(fig, ims, interval=100)
    if(x_y == 1):
      ani.save('x_y_dens.mp4', writer="ffmpeg", dpi=600)
    if(x_z == 1):
      ani.save('x_z_dens.mp4', writer="ffmpeg", dpi=600)
    plt.close()