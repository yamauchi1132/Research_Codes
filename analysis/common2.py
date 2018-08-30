#This is included into energy_check.py.
import numpy as np

class Particle:
	def __init__(self):
		self.p_id = 0.
		self.mass = 0.
		self.pos = np.array([0.,0.,0.])
		self.vel = np.array([0.,0.,0.])
		self.acc = np.array([0.,0.,0.])
		self.dens = 0.
		self.r_cg = 0.


def readfile(data, p, lfile):
	for i in range(lfile):
		p[i].p_id = data[i,0] 
		p[i].mass = data[i,2]
		p[i].pos[0] = data[i,3]
		p[i].pos[1] = data[i,4]
		p[i].pos[2] = data[i,5]
		p[i].vel[0] = data[i,6]
		p[i].vel[1] = data[i,7]
		p[i].vel[2] = data[i,8]
		p[i].acc[0] = data[i,9]
		p[i].acc[1] = data[i,10]
		p[i].acc[2] = data[i,11]
		p[i].dens = data[i,15]

def calc_center_of_gravity(p):
	row_sum = 0.
	pos_cg = np.array([0.,0.,0.])

	for i in range(len(p)):
		row_sum += p[i].dens
		for j in range(3):
			pos_cg[j] += p[i].dens * p[i].pos[j]

	for j in range(3):
		pos_cg[j] = pos_cg[j] / row_sum

	return pos_cg