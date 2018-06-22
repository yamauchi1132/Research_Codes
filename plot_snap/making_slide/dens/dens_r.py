import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
import math

if __name__ == '__main__':
  r_real = []
  dens_real = []
  
  r_0 = []
  dens_0 = []
  
  r_end = []
  dens_end = []

  data_real = np.loadtxt("real.data")
  data_0 = np.loadtxt("sph_t0000.dat")
  data_end = np.loadtxt("final.dat")

  r_real.extend(data_real[:,4]/1e+10)
  dens_real.extend(data_real[:,5])

  for i in range(len(data_0)):
    r_0.append(math.sqrt(data_0[i,3]**2+data_0[i,4]**2+data_0[i,5]**2)/1e+10)
    dens_0.append(data_0[i,15])


  for i in range(len(data_end)):
    r_end.append(math.sqrt(data_end[i,3]**2+data_end[i,4]**2+data_end[i,5]**2)/1e+10)
    dens_end.append(data_end[i,15])

  fig = plt.figure()
  plt.plot(r_real, dens_real, label="polytropic sphere")
  plt.plot(r_0, dens_0, '.', label="before damping")
  plt.plot(r_end, dens_end, '.', label="after damping")
  
  # plt.xscale('log')
  # plt.yscale('log')
  plt.xlabel(r'$r\,[10^{10}cm]$', fontsize=18)
  plt.ylabel(r'$\rho\,[g\,cm^{-3}]$', fontsize=18)
  plt.tick_params(labelsize=18)
  plt.legend(fontsize=18)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  
  # plt.show()
  plt.savefig("dens_r.png", dpi=600)
  plt.close()
