import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp

fsize = 14
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

  data_c = np.loadtxt("time_cpu_two_body.data")
  data_g = np.loadtxt("time_gpu_two_body.data")

  p_num.extend(data_c[:,0])
  dens_t.extend(data_c[:,2])
  hydro_t.extend(data_c[:,3])
  grav_t.extend(data_c[:,4])
  all_t.extend(data_c[:,1])

  p_numg.extend(data_g[:,0])
  dens_tg.extend(data_g[:,2])
  hydro_tg.extend(data_g[:,3])
  grav_tg.extend(data_g[:,4])
  all_tg.extend(data_g[:,1])

  fig = plt.figure()
  x = p_num
  y = all_t
  y2 = all_tg
  plt.plot(x, y, label="CPU(core i5 4core+SIMD)", marker='o', ms=10, linewidth=2.0, color='black')
  plt.plot(x, y2, label="GPU(GTX 1080)", marker='D', ms=10, linewidth=2.0, color='black', linestyle='dashed')

  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel('Number of particles', fontsize=fsize)
  plt.ylabel('Time (s)', fontsize=fsize)
  plt.tick_params(labelsize=fsize)
  plt.legend(fontsize=fsize)
  mpl.rcParams['axes.xmargin'] = 0
  mpl.rcParams['axes.ymargin'] = 0
  plt.tight_layout()
  
  # plt.show()
  plt.savefig("time_two_body.pdf", transparent=True, dpi=300, bbox_inches = 'tight', pad_inches = 0)
  plt.close()