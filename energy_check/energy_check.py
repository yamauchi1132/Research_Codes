import sys
import numpy as np
import operator
import pandas
from user_define import Particle

G = 6.67259e-8

def readfile(data, p):
	for i in range(len(data)):
		p[i].p_id = data[i,0] 
		p[i].istar = data[i,1]
		p[i].mass = data[i,2]
		p[i].posx = data[i,3]
		p[i].posy = data[i,4]
		p[i].posz = data[i,5]
		p[i].velx = data[i,6]
		p[i].vely = data[i,7]
		p[i].velz = data[i,8]
		p[i].accx = data[i,9]
		p[i].accy = data[i,10]
		p[i].accz = data[i,11]
		p[i].uene = data[i,12]
		p[i].dalph = data[i,13]
		p[i].alphu = data[i,14]
		p[i].dens = data[i,15]
		p[i].ksr = data[i,16]
		p[i].np = data[i,17]
		p[i].vsnd = data[i,18]
		p[i].pres = data[i,19]
		p[i].emp = data[i,20]
		p[i].divv = data[i,21]
		p[i].rotv = data[i,22]
		p[i].bswt = data[i,23]
		p[i].pot = data[i,24]
		p[i].abar = data[i,25]
		p[i].zbar = data[i,26]
		p[i].enuc = data[i,27]
		p[i].vsmx = data[i,28]
		p[i].udot = data[i,29]
		p[i].dnuc = data[i,30]
		for j in range(18):
			p[i].cmps[j] = data[i, 31+j]
		
def calc_point_energy(p):
	row1_sum = 0.
	row2_sum = 0.
	pos_sum = np.array([0.,0.,0.,0.]).reshape(2,2)
	vel_sum = np.array([0.,0.,0.,0.]).reshape(2,2)

	for i in range(int(len(p)//2)):
		row1_sum += p[i].dens
		pos_sum[0][0] += p[i].dens * p[i].posx
		pos_sum[0][1] += p[i].dens * p[i].posy
		vel_sum[0][0] += p[i].dens * p[i].velx
		vel_sum[0][1] += p[i].dens * p[i].vely

	for i in range(int(len(p)//2), len(p)):
		row2_sum += p[i].dens
		pos_sum[1][0] += p[i].dens * p[i].posx
		pos_sum[1][1] += p[i].dens * p[i].posy
		vel_sum[1][0] += p[i].dens * p[i].velx
		vel_sum[1][1] += p[i].dens * p[i].vely
	
	point_pos1_x = pos_sum[0][0] / row1_sum;
	point_pos1_y = pos_sum[0][1] / row1_sum;
	point_pos2_x = pos_sum[1][0] / row2_sum;
	point_pos2_y = pos_sum[1][1] / row2_sum;

	point_vel1_x = vel_sum[0][0] / row1_sum;
	point_vel1_y = vel_sum[0][1] / row1_sum;
	point_vel2_x = vel_sum[1][0] / row2_sum;
	point_vel2_y = vel_sum[1][1] / row2_sum;


def calc_energy_difference(p_s, p_f):
	ene_s = calc_energy(p_s)

if __name__ == '__main__':
	args = sys.argv

	if(len(args) < 4):
		sys.stderr.write('Error : no input file\n')
		exit()

	data_s = np.loadtxt(args[1])
	data_f = np.loadtxt(args[2])
	data_time = np.loadtxt(args[3], usecols=(2,4))

	p_s = [Particle() for i in range(len(data_s))]
	p_f = [Particle() for i in range(len(data_f))]

	readfile(data_s, p_s)
	readfile(data_f, p_f)

	p_s.sort(key=operator.attrgetter("p_id"))
	p_f.sort(key=operator.attrgetter("p_id"))
	
	'''
	for i in range(int(len(p_s)//2)):
		print('%d %e %e %e %e' % (p_s[i].p_id, p_s[i].mass, p_s[i].posx, p_s[i].posy, p_s[i].dens))
		print('%d %e %e %e %e' % (p_f[i].p_id, p_f[i].mass, p_f[i].posx, p_f[i].posy, p_s[i].dens))
		pass
	'''
	ene_diff = calc_energy_difference(p_s, p_f)
