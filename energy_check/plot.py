#This is included into energy_check.py.
import matplotlib as mpl
import matplotlib.pyplot as plt

def visualization(time_data, data_name):
  time = time_data[:,0]
  ene = time_data[:,1]

  fig = plt.figure()
  x = time
  y = ene - time_data[0,1]

  plt.plot(x, y, label='time.log', linewidth = 0.5)
  plt.title(data_name)

  plt.show()
  plt.close()