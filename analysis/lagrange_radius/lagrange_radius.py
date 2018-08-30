import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common2 import *

dirname = "../data/sph_t%04d.dat"
start = 600
end = 800
step = 100

Msun = 1.989e+33
Rsun = 695700e+5;

def calc_r_cg(p, pos_cg):
	dx = p.pos[0] - pos_cg[0]
	dy = p.pos[1] - pos_cg[1]
	dz = p.pos[2] - pos_cg[2]
	r_cg = math.sqrt(dx*dx + dy*dy + dz*dz)
	p.r_cg = r_cg

def calc_lagrange_radius(p, lagrange_r, mass_total, lfile):
	mass_ratio = []
	for i in range(1, 10):
		m = mass_total * 0.1 * i
		mass_ratio.append(m)

	mass = 0.
	j = 0
	for i in range(lfile):
		mass += p[i].mass
		if(mass > mass_ratio[j]):
			lagrange_r[j].append(p[i-1].r_cg)
			# print("%e %e\n" %(mass, p[i-1].r_cg))
			j += 1
		if(j == 9):
			break

def plot(t, lagrange_r):
	fig = plt.figure()

	for i in range(9):
		for j in range(len(lagrange_r[i])):
			lagrange_r[i][j] = lagrange_r[i][j] / Rsun


	for i in range(9):
		plt.plot(t, lagrange_r[i], label="%d per"%((i+1)*10))

	plt.xlabel('time [s]', fontsize=12)
	plt.ylabel('radius [Rsun]', fontsize=12)
	plt.yscale('log')
	plt.tick_params(labelsize=12)
	plt.legend(fontsize=12)
	mpl.rcParams['axes.xmargin'] = 0
	mpl.rcParams['axes.ymargin'] = 0
	plt.tight_layout()

	plt.show()
	plt.close()
  
if __name__ == '__main__':
	args = sys.argv

	pre_data = np.loadtxt(dirname % start)
	lfile = len(pre_data)
	mass_total = 0.
	for i in range(lfile):
		mass_total += pre_data[i,2]

	sys.stderr.write('Total mass = %.1lf Msun\n'%(mass_total/Msun))

	t = []
	lagrange_r = [[] for i in range(9)]

	for time in range(start, end+step, step):
		data = np.loadtxt(dirname % time)
		p = [Particle() for i in range(lfile)]
		readfile(data, p, lfile)
		
		pos_cg = np.array([0.,0.,0.])
		pos_cg = calc_center_of_gravity(p)

		for i in range(lfile):
			calc_r_cg(p[i], pos_cg)

		p.sort(key=operator.attrgetter("r_cg"))	

		calc_lagrange_radius(p, lagrange_r, mass_total, lfile)
		t.append(time*1e+04)

	plot(t, lagrange_r)
