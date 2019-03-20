#This is included by each program
import numpy as np

class Particle:
  def __init__(self):
    self.p_id = 0.
    self.istar = 0.
    self.mass = 0.
    self.pos = np.array([0.,0.,0.])
    self.vel = np.array([0.,0.,0.])
    self.acc = np.array([0.,0.,0.])
    self.uene = 0.
    self.dalph = 0.
    self.alphu = 0.
    self.dens = 0.
    self.ksr = 0.
    self.np = 0.
    self.vsnd = .0
    self.pres = 0.
    self.emp = 0.
    self.divv = 0.
    self.rotv = 0.
    self.bswt = 0.
    self.pot = 0.
    self.abar = 0.
    self.zbar = 0.
    self.enuc = 0.
    self.vsmx = 0.
    self.udot = 0.
    self.dnuc = 0.
    self.cmps = [0. for i in range(18)]
    self.r_cg = 0.
    self.temp = 0.

def readfile(data, p):
  for i in range(len(data)):
    p[i].p_id = data[i,0] 
    p[i].istar = data[i,1]
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

def calc_center_of_gravity(p):
  row_sum = 0.
  pos_cg = np.array([0.,0.,0.])
  vel_cg = np.array([0.,0.,0.])

  for i in range(len(p)):
    row_sum += p[i].dens
    pos_cg[0] += p[i].dens * p[i].pos[0]
    pos_cg[1] += p[i].dens * p[i].pos[1]
    pos_cg[2] += p[i].dens * p[i].pos[2]

    vel_cg[0] += p[i].dens * p[i].vel[0]
    vel_cg[1] += p[i].dens * p[i].vel[1]
    vel_cg[2] += p[i].dens * p[i].vel[2]
  
  for i in range(3):
    pos_cg[i] = pos_cg[i] / row_sum
    vel_cg[i] = vel_cg[i] / row_sum

  return pos_cg, vel_cg
