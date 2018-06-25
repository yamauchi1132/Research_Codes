import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib import gridspec

box_size = 6.0
Nnode = 501
ims = []
fig = plt.figure(figsize=(12,5))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.2)
# gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1]) 
plt.rcParams["font.size"] = 13
if __name__ == '__main__':
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

   
    ax1 = fig.add_subplot(1,2,1, adjustable='box', aspect=1)
    ax2 = fig.add_subplot(1,2,2,sharey=ax1, sharex=ax1)
    # ax1 = plt.subplot(gs[0])
    # ax2 = plt.subplot(gs[1])
    im1 = ax1.scatter(x, y, c=dens, s=0.001, cmap="jet")
    im2 = ax2.scatter(x, z, c=dens, s=0.001, cmap="jet")

    ax1.set_xlim(-box_size, box_size)
    ax1.set_ylim(-box_size, box_size)
    ax2.set_xlim(-box_size, box_size)
    ax2.set_ylim(-box_size, box_size)

    ax1.set_xlabel(r'$x\,[10^{11}cm]$')
    ax2.set_xlabel(r'$x\,[10^{11}cm]$')
    ax1.set_ylabel(r'$y\,[10^{11}cm]$')
    ax2.set_ylabel(r'$z\,[10^{11}cm]$')
    # ax2.tick_params(labelleft="off", left="off")


    if(n == 0):
      clb1 = plt.colorbar(im1, ax=ax1)
      clb2 = plt.colorbar(im2, ax=ax2)
      clb1.set_label(label=r'$\rho\,[g\,cm^{-3}]$')
      clb2.set_label(label=r'$\rho\,[g\,cm^{-3}]$')
      im1.set_clim(vmin=0, vmax=10)
      im2.set_clim(vmin=0, vmax=10)

    # plt.gca().set_aspect('equal', adjustable='box')
    mpl.rcParams['axes.xmargin'] = 0
    mpl.rcParams['axes.ymargin'] = 0
    plt.tight_layout()

    ims.append([im1]+[im2])

  ani = animation.ArtistAnimation(fig, ims, interval=100)
  ani.save('x_yAndx_z.mp4', writer="ffmpeg", dpi=400)
  plt.close()


