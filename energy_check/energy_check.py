import sys
import numpy as np
import operator
import math
from user_define import Particle
from user_define import readfile

msun = 1.989e+33
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

def check_energy_error(time_data, max_time):
	ene_init = time_data[0,1]
	max_ene_error = 0.
	time = 0.
	for i in range(len(time_data)):
		if(time_data[i,0] >= max_time):
			ene_error_end = abs(time_data[i,1] - ene_init)
			time_end = time_data[i,0]
			if(ene_error_end >= max_ene_error):
				max_ene_error = ene_error_end
				time_max = time_data[i,0]
			break

		ene_error = abs(time_data[i,1] - ene_init)
		if(ene_error > max_ene_error):
			max_ene_error = ene_error
			time_max = time_data[i,0]

	return max_ene_error, time_max, ene_error_end, time_end

def calc_energy_difference(p_s, p_f):
	ene_s = calc_energy(p_s)
	ene_f = calc_energy(p_f)
	ene_diff = abs(ene_f - ene_s)

	sys.stderr.write('\nTotal Energy Of Start(E_s) : %e\n' %ene_s)
	sys.stderr.write('Total Energy Of End(E_e) : %e\n' %ene_f)
	sys.stderr.write('Energy Differece(E_d = |E_s - E_e|) : %e\n' %ene_diff)

	return ene_diff

if __name__ == '__main__':
	args = sys.argv

	if(len(args) < 4):
		sys.stderr.write('Error : no input file\n')
		exit()

	data_s = np.loadtxt(args[1])
	data_f = np.loadtxt(args[2])
	time_data = np.loadtxt(args[3], usecols=(2,4))

	sys.stderr.write('R1 = %.1lfRsun, R2 = %.1lfRsun\n' %(R1/rsun, R2/rsun))

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
	ene_error_rate = (ene_error_end / ene_diff) * 100

	sys.stderr.write('\nMax Energy Error(|E_m|) : %e, Time : %lf\n' %(max_ene_error, max_time))
	sys.stderr.write('Energy Error Of End(|E_e|) : %e, Time : %lf\n\n' %(ene_error_end, end_time))
	sys.stderr.write('Energy Error Rate(E_r = |E_e / E_d|): %.1lf per\n\n' %ene_error_rate)
