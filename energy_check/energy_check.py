import sys
import numpy as np
import operator
import math
from user_define import Particle
from user_define import readfile

rsun = 695700e+5
G = 6.67259e-8

####### you need to set #########
R1 = 1.0 * rsun
R2 = 1.0 * rsun
# R1 = 4.0 * rsun
# R2 = 4.0 * rsun
# R1 = 17.2 * rsun
# R2 = 17.2 * rsun
#################################

def calc_energy(p):
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

	point_m1 = 0.
	point_m2 = 0.
	for i in range(int(len(p)//2)):
		dx = p[i].posx - point_pos1_x
		dy = p[i].posy - point_pos1_y
		r_2 = dx*dx + dy*dy
		r = math.sqrt(r_2)
		if(r < 2*R1):
			point_m1 += p[i].mass

	for i in range(int(len(p)//2), len(p)):
		dx = p[i].posx - point_pos2_x
		dy = p[i].posy - point_pos2_y
		r_2 = dx*dx + dy*dy
		r = math.sqrt(r_2)
		if(r < 2*R2):
			point_m2 += p[i].mass

	vel1_2 = point_vel1_x*point_vel1_x + point_vel1_y*point_vel1_y
	vel2_2 = point_vel2_x*point_vel2_x + point_vel2_y*point_vel2_y
	k_e = 0.5*point_m1*vel1_2 + 0.5*point_m2*vel2_2

	dx = point_pos1_x - point_pos2_x
	dy = point_pos1_y - point_pos2_y
	p_e = (G*point_m1*point_m2) / math.sqrt(dx*dx + dy*dy)

	ene = k_e - p_e
	return ene

def calc_energy_difference(p_s, p_f):
	ene_s = calc_energy(p_s)
	ene_f = calc_energy(p_f)
	ene_diff = abs(ene_f - ene_s)

	sys.stderr.write('\nTotal Energy Of Start(E_s) : %e\n' %ene_s)
	sys.stderr.write('Total Energy Of End(E_e) : %e\n' %ene_f)
	sys.stderr.write('Energy Differece(dE = |E_s - E_e|) : %e\n' %ene_diff)

	return ene_diff

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
