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
	result3.extend(data[3, :])
	result4.extend(data[4, :])

	fig = plt.figure()
	plt.plot(mass, result1, color='blue', linestyle='solid', marker='o')
	plt.plot(mass, result2, color='red', linestyle='solid', marker='v')
	plt.plot(mass, result3, color='blue', linestyle='dashed', marker='o')
	plt.plot(mass, result4, color='red', linestyle='dashed', marker='v')
	# label=r'$n=1.5,v_{\infty}=10 km\,s^{-1}$'
	plt.xlabel(r'$Mass\,[M_{\odot}]$', fontsize=18)
	plt.ylabel(r'$\frac{r_{\rm{p}}}{R_1+R_2}$', fontsize=18)
	plt.xscale("log")
	plt.tick_params(labelsize=18)
	plt.legend(fontsize=12)
	mpl.rcParams['axes.xmargin'] = 0
	mpl.rcParams['axes.ymargin'] = 0
	plt.tight_layout()
	
	plt.savefig("result.png", transparent=True, dpi=300, bbox_inches = 'tight', pad_inches = 0)
	plt.show()
  