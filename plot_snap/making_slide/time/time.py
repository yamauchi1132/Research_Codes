import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp

if __name__ == '__main__':
  p_num = []
  dens_t = []
  hydro_t = []
  grav_t = []
  all_t = []

  p_numg = []
  dens_tg = []
  hydro_tg = []
  grav_tg = []
  all_tg = []

  data_c = np.loadtxt("time_cpu.data")
  data_g = np.loadtxt("time_gpu.data")

  p_num.extend(data_c[:,0])
  dens_t.extend(data_c[:,1])
  hydro_t.extend(data_c[:,2])
  grav_t.extend(data_c[:,3])
  all_t.extend(data_c[:,4])

  p_numg.extend(data_g[:,0])
  dens_tg.extend(data_g[:,1])
  hydro_tg.extend(data_g[:,2])
  grav_tg.extend(data_g[:,3])
  all_tg.extend(data_g[:,4])

  fig = plt.figure()
  x = p_num
  y = all_t
  y2 = all_tg
  plt.plot(x, y, label="CPU(core i5 4core)")
  plt.plot(x, y2, label="GPU(GTX 1080)")

  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel('number of particle', fontsize=18)
  plt.ylabel('time (s)', fontsize=18)
  plt.tick_params(labelsize=18)
  plt.legend(fontsize=18)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  
  #plt.show()
  plt.savefig("time.png", dpi=600)
  plt.close()