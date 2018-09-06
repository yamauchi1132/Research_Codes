import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
import sys, os
sys.path.append(os.pardir)

if __name__ == '__main__':
	args = sys.argv
	if(len(args) < 2):
		sys.stderr.write('Error : no input file\n')
		exit()

	data = np.loadtxt(args[1])
	mass = []
	result1 = []
	result2 = []
	result3 = []
	result4 = []

	mass.extend(data[0, :])
	result1.extend(data[1, :])
	result2.extend(data[2, :])
	# result3.extend(data[3, :])
	# result4.extend(data[4. :])

	fig = plt.figure()
	plt.plot(mass, result1, label="pori1.5, vel10")
	plt.plot(mass, result2, label="pori2.5, vel10")

	plt.xlabel('Mass[Msun]', fontsize=18)
	plt.ylabel('Rp(R1+R2)', fontsize=18)
	plt.tick_params(labelsize=18)
	plt.legend(fontsize=18)
	mpl.rcParams['axes.xmargin'] = 0
	mpl.rcParams['axes.ymargin'] = 0
	plt.tight_layout()
	plt.show()
  