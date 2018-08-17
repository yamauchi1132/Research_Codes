import sys
import numpy as np
import operator
import math
from common import *

## visualization(if view = 0, no plot. if view = 1, plot) ##
view = 1
############################################################
if(view == 1):
	import matplotlib as mpl
	import matplotlib.pyplot as plt
	from plot import visualization

msun = 1.989e+33
rsun = 695700e+5
G = 6.67259e-8

####### you need to set #########
#R1 = 1.0 * rsun
#R2 = 1.0 * rsun
#R1 = 4.0 * rsun
#R2 = 4.0 * rsun
R1 = 17.2 * rsun
R2 = 17.2 * rsun
#################################

def calc_energy(p):
	row1_sum = 0.
	row2_sum = 0.
	pos_sum = np.array([0.,0.,0.,0.,0.,0.]).reshape(2,3)
	vel_sum = np.array([0.,0.,0.,0.,0.,0.]).reshape(2,3)

	for i in range(int(len(p)//2)):
		row1_sum += p[i].dens
		pos_sum[0][0] += p[i].dens * p[i].posx
		pos_sum[0][1] += p[i].dens * p[i].posy
		pos_sum[0][2] += p[i].dens * p[i].posz
		vel_sum[0][0] += p[i].dens * p[i].velx
		vel_sum[0][1] += p[i].dens * p[i].vely
		vel_sum[0][2] += p[i].dens * p[i].velz

	for i in range(int(len(p)//2), len(p)):
		row2_sum += p[i].dens
		pos_sum[1][0] += p[i].dens * p[i].posx
		pos_sum[1][1] += p[i].dens * p[i].posy
		pos_sum[1][2] += p[i].dens * p[i].posz
		vel_sum[1][0] += p[i].dens * p[i].velx
		vel_sum[1][1] += p[i].dens * p[i].vely
		vel_sum[1][2] += p[i].dens * p[i].velz
	
	point_pos1_x = pos_sum[0][0] / row1_sum
	point_pos1_y = pos_sum[0][1] / row1_sum
	point_pos1_z = pos_sum[0][2] / row1_sum
	point_pos2_x = pos_sum[1][0] / row2_sum
	point_pos2_y = pos_sum[1][1] / row2_sum
	point_pos2_z = pos_sum[1][2] / row2_sum

	point_vel1_x = vel_sum[0][0] / row1_sum
	point_vel1_y = vel_sum[0][1] / row1_sum
	point_vel1_z = vel_sum[0][2] / row1_sum
	point_vel2_x = vel_sum[1][0] / row2_sum
	point_vel2_y = vel_sum[1][1] / row2_sum
	point_vel2_z = vel_sum[1][2] / row2_sum

	point_m1 = 0.
	point_m2 = 0.
	for i in range(int(len(p)//2)):
		dx = p[i].posx - point_pos1_x
		dy = p[i].posy - point_pos1_y
		dz = p[i].posz - point_pos1_z
		r_2 = dx*dx + dy*dy + dz*dz
		r = math.sqrt(r_2)
		if(r < 2*R1):
			point_m1 += p[i].mass

	for i in range(int(len(p)//2), len(p)):
		dx = p[i].posx - point_pos2_x
		dy = p[i].posy - point_pos2_y
		dz = p[i].posz - point_pos2_z
		r_2 = dx*dx + dy*dy + dz*dz
		r = math.sqrt(r_2)
		if(r < 2*R2):
			point_m2 += p[i].mass

	vel1_2 = point_vel1_x*point_vel1_x + point_vel1_y*point_vel1_y + point_vel1_z*point_vel1_z
	vel2_2 = point_vel2_x*point_vel2_x + point_vel2_y*point_vel2_y + point_vel2_z*point_vel2_z
	k_e = 0.5*point_m1*vel1_2 + 0.5*point_m2*vel2_2

	dx = point_pos1_x - point_pos2_x
	dy = point_pos1_y - point_pos2_y
	dz = point_pos1_z - point_pos2_z
	p_e = (G*point_m1*point_m2) / math.sqrt(dx*dx + dy*dy + dz*dz)

	ene = k_e - p_e
	return ene

def calc_energy_difference(p_s, p_f):
	ene_s = calc_energy(p_s)
	ene_f = calc_energy(p_f)
	ene_diff = ene_f - ene_s

	sys.stderr.write('\nTotal Energy Of Start(E_s) : %e\n' %ene_s)
	sys.stderr.write('Total Energy Of End(E_e) : %e\n' %ene_f)
	sys.stderr.write('Energy Differece(E_d = E_e - E_s) : %e\n' %ene_diff)

	return ene_diff

def check_energy_error(time_data, max_time):
	ene_init = time_data[0,1]
	max_ene_error = 0.
	time = 0.
	for i in range(len(time_data)):
		if(time_data[i,0] >= max_time):
			ene_error_end = time_data[i,1] - ene_init
			time_end = time_data[i,0]
			if(abs(ene_error_end) >= abs(max_ene_error)):
				max_ene_error = ene_error_end
				time_max = time_data[i,0]
			break

		ene_error = time_data[i,1] - ene_init
		if(abs(ene_error) > abs(max_ene_error)):
			max_ene_error = ene_error
			time_max = time_data[i,0]

	return max_ene_error, time_max, ene_error_end, time_end

if __name__ == '__main__':
	args = sys.argv

	if(len(args) < 4):
		sys.stderr.write('Error : no input file\n')
		exit()

	sys.stderr.write('R1 = %.1lfRsun, R2 = %.1lfRsun\n' %(R1/rsun, R2/rsun))

	data_s = np.loadtxt(args[1])
	data_f = np.loadtxt(args[2])
	time_data = np.loadtxt(args[3], usecols=(2,4))

	if(len(args) < 5):
		max_timestep = time_data[len(time_data)-1,0]
	elif(float(args[4]) > time_data[len(time_data)-1,0]):
		sys.stderr.write('Error : max timestep is larger than final time of time.log\n')
		exit()
	else:
		max_timestep = float(args[4])

	p_s = [Particle() for i in range(len(data_s))]
	p_f = [Particle() for i in range(len(data_f))]

	readfile(data_s, p_s)
	readfile(data_f, p_f)

	p_s.sort(key=operator.attrgetter("p_id"))
	p_f.sort(key=operator.attrgetter("p_id"))

	ene_diff = calc_energy_difference(p_s, p_f)

	max_ene_error, max_time, ene_error_end, end_time = check_energy_error(time_data, max_timestep)
	ene_error_rate = abs(ene_error_end / ene_diff) * 100

	sys.stderr.write('\nMax Energy Error(E_m) : %e, Time : %lf\n' %(max_ene_error, max_time))
	sys.stderr.write('Energy Error Of End(E_e) : %e, Time : %lf\n\n' %(ene_error_end, end_time))
	sys.stderr.write('Energy Error Rate(E_r = |E_e / E_d| * 100): %.3lf per\n\n' %ene_error_rate)

	if(view == 1):
		visualization(time_data, args[3])
