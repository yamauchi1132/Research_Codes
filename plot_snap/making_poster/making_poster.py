#Code revies is hear : https://qiita.com/NatsukiLab/items/83e6d31aa4a49acc528c
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp

box_size = 5.5
Nnode = 9
# Nnode = 3
dat_num = [0, 21, 41, 96, 153, 178, 191, 212, 299]
# dat_num = [0,191, 299]


tag = [[], [], [], [], [], [], [], [], []]
x = [[], [], [], [], [], [], [], [], []]
y = [[], [], [], [], [], [], [], [], []]
dens = [[], [], [], [], [], [], [], [], []]

'''
tag = [[], [], []]
x = [[], [], []]
y = [[], [], []]
dens = [[], [], []]
'''
for n in range(Nnode):
  data = np.loadtxt("sph_t%04d.dat" % (dat_num[n]))
  tag[n].append(data[:,0])
  x[n].append(data[:,3])
  y[n].append(data[:,4])
  dens[n].append(data[:,15])

  for i in range(len(x[n])):
    x[n][i] = x[n][i]/1e+11
    y[n][i] = y[n][i]/1e+11

ax = [[], [], [], [], [], [], [], [], []]
plot = [[], [], [], [], [], [], [], [], []]
side = 3
ver = Nnode // side
fig = plt.figure()
# fig = plt.figure(figsize=(16.2, 8.1)) #figure size is 1920*999
# plt.subplots_adjust(left=0.075, bottom=0.07, right=0.95, top=0.93, wspace=0, hspace=0.2)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0.2)
plt.rcParams["font.size"] = 15

for n in range(Nnode):
  ax[n] = fig.add_subplot(ver, side, n+1)
  plot[n] = ax[n].scatter(x[n],y[n], c=dens[n], s=0.001, cmap="jet")
  ax[n].set_xlim(-box_size, box_size)
  ax[n].set_ylim(-box_size, box_size)

#erase y axis memory other than bottom end
  if(n!=0 and n!=3 and n!=6):
    ax[n].tick_params(labelleft="off",left="off") # y軸の削除  
  
#erase x axis memory other than left side end
  if(n!=6 and n!= 7 and n!= 8):
    ax[n].tick_params(labelbottom="off",bottom="off") # x軸の削除

#set the ylabel on bottom end
  if(n==0 or n==3 or n==6):
    ax[n].set_ylabel(r'$y\,[10^{11}cm]$')

#set the xlabel on left side end
  if(n==6 or n==7 or n==8):
    ax[n].set_xlabel(r'$x\,[10^{11}cm]$')

#set the ratio 1
  ax[n].set_aspect('equal', adjustable='box')

#set the colorbar on right side end
  if(n == 2 or n == 5 or n == 8):
    clb = plt.colorbar(plot[n],ax=ax[n])
    clb.set_label(label=r'$\rho\,[g\,cm^{-3}]$')
    plot[n].set_clim(vmin=0, vmax=10)

#plt.savefig("image.png", dpi=600)
#plt.savefig("imageGTC.png", dpi=600)
plt.show()
plt.close()
del x, y, dens, data, ax