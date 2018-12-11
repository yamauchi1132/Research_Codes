import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

msun = 1.989e+33
rsun = 695700e+5
G = 6.67259e-8

########set parameter #######
#vinf = 1e+06
#M1 = 1.0 * msun
#M2 = 1.0 * msun
############################
def calc_kinetic_energy(p):
  point_m1 = 0.
  point_m2 = 0.

  for i in range(int(len(p)//2)):
    point_m1 += p[i].mass

  for i in range(int(len(p)//2), len(p)):
    point_m2 += p[i].mass

  row_sum = np.array([0.,0.])
  pos_sum = np.array([0.,0.,0.,0.,0.,0.]).reshape(2,3)
  vel_sum = np.array([0.,0.,0.,0.,0.,0.]).reshape(2,3)

  for i in range(int(len(p)//2)):
    row_sum[0] += p[i].dens
    for j in range(3):
      pos_sum[0][j] += p[i].dens * p[i].pos[j]
      vel_sum[0][j] += p[i].dens * p[i].vel[j]

  for i in range(int(len(p)//2), len(p)):
    row_sum[1] += p[i].dens
    for j in range(3):
      pos_sum[1][j] += p[i].dens * p[i].pos[j]
      vel_sum[1][j] += p[i].dens * p[i].vel[j]

  point_pos = np.array([0.,0.,0.,0.,0.,0.]).reshape(2,3)
  point_vel = np.array([0.,0.,0.,0.,0.,0.]).reshape(2,3)

  for i in range(2):
    for j in range(3):
      point_pos[i][j] = pos_sum[i][j] / row_sum[i]
      point_vel[i][j] = vel_sum[i][j] / row_sum[i]

  vel1_2 = point_vel[0][0]*point_vel[0][0] + point_vel[0][1]*point_vel[0][1] + point_vel[0][2]*point_vel[0][2]
  vel2_2 = point_vel[1][0]*point_vel[1][0] + point_vel[1][1]*point_vel[1][1] + point_vel[1][2]*point_vel[1][2]
  k_e = 0.5*point_m1*vel1_2 + 0.5*point_m2*vel2_2

  dx = point_pos[0][0] - point_pos[1][0]
  dy = point_pos[0][1] - point_pos[1][1]
  dz = point_pos[0][2] - point_pos[1][2]
  p_e = (G*point_m1*point_m2) / math.sqrt(dx*dx + dy*dy + dz*dz)

  ene = k_e - p_e
  return ene

def calc_energy_ratio(p, total_energy):
  kinetic_ene = calc_kinetic_energy(p)
  binding_ene = kinetic_ene - total_energy
  #k_e_inf = 0.5 * ((M1*M2)/(M1+M2)) * vinf * vinf
  energy_ratio = binding_ene / kinetic_ene

  m = 0.
  for i in range(len(p)//2):
    m += p[i].mass
  r_eff = (2 * G * m * m) / binding_ene
  print("%e"%r_eff)
  return energy_ratio

if __name__ == '__main__':
  args = sys.argv

  if(len(args) < 3):
    sys.stderr.write("Error : no input file\n")
    exit()

  data = np.loadtxt(args[1])
  time_data = np.loadtxt(args[2], usecols=(2,4))

  p = [Particle() for i in range(len(data))]
  total_energy = time_data[0,1]

  readfile(data, p)

  p.sort(key=operator.attrgetter("p_id"))

  energy_ratio = calc_energy_ratio(p, total_energy)
  sys.stderr.write('Energy ratio(kinetic_ene/binding_ene) : %e\n' %energy_ratio)