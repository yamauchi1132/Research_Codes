#This is included into energy_check.py.

class Particle:
  p_id = 0.
  istar = 0.
  mass = 0.
  posx = 0.
  posy = 0.
  posz = 0.
  velx = 0.
  vely = 0.
  velz = 0.
  accx = 0.
  accy = 0.
  accz = 0.
  uene = 0.
  dalph = 0.
  alphu = 0.
  dens = 0.
  ksr = 0.
  np = 0.
  vsnd = .0
  pres = 0.
  emp = 0.
  divv = 0.
  rotv = 0.
  bswt = 0.
  pot = 0.
  abar = 0.
  zbar = 0.
  enuc = 0.
  vsmx = 0.
  udot = 0.
  dnuc = 0.
  cmps = [0. for i in range(18)]

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
