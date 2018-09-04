import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

if __name__ == '__main__':
	args = sys.argv
	if(len(args) < 2):
		sys.stderr.write('Error : no input file\n')
		exit()

	for i in range(1, len(args)):
		data = np.loadtxt(args[i])
		p = [Particle() for i in range(len(data))]
		readfile(data, p)



