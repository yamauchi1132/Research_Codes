import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

box_size = 1e+12

def plot(x_tag, y_tag, z_tag, pos_cg):
	fig = plt.figure()
	clr = ["red", "blue"]

	for tag in range(2):
		for i in range(len(x_tag[tag])):
			x_tag[tag][i] = x_tag[tag][i] - pos_cg[0]
			y_tag[tag][i] = y_tag[tag][i] - pos_cg[1]
			z_tag[tag][i] = z_tag[tag][i] - pos_cg[2]
	plt.scatter(x_tag[0],y_tag[0], c=clr[0], s=1, edgecolor=clr[0], alpha=0.1)
	plt.scatter(x_tag[1],y_tag[1], c=clr[1], s=1, edgecolor=clr[1], alpha=0.1)

	plt.xlim(-box_size, box_size)
	plt.ylim(-box_size, box_size)
	plt.xlabel(r'$x$')
	plt.ylabel(r'$y$')
	plt.gca().set_aspect('equal', adjustable='box')
	plt.show()
	plt.close()

def plot2(x_tag, y_tag, z_tag, pos_cg):
  fig = plt.figure()
  clr = ["red", "blue"]

  for tag in range(2):
    for i in range(len(x_tag[tag])):
      x_tag[tag][i] = x_tag[tag][i] - pos_cg[0]
      y_tag[tag][i] = y_tag[tag][i] - pos_cg[1]
      z_tag[tag][i] = z_tag[tag][i] - pos_cg[2]

  ax = [[],[]]
  for i in range(2):
    ax[i] = fig.add_subplot(1, 2, i+1, aspect=1.0)
    ax[i].scatter(x_tag[i],y_tag[i], c=clr[i], s=1, edgecolor=clr[i], alpha=0.1)
    ax[i].set_xlim(-box_size, box_size)
    ax[i].set_ylim(-box_size, box_size)
    ax[i].set_xlabel(r'$x$')
    ax[i].set_ylabel(r'$y$')
  plt.show()

def calc_cg(x, y, z, dens):
	row_sum = 0.
	pos_cg = np.array([0., 0., 0.])

	for i in range(len(x)):
		row_sum += dens[i]
		pos_cg[0] += dens[i] * x[i]
		pos_cg[1] += dens[i] * y[i]
		pos_cg[2] += dens[i] * z[i]

	for i in range(3):
		pos_cg[i] = pos_cg[i] / row_sum

	return pos_cg


if __name__ == '__main__':
	args = sys.argv
	if(len(args) < 2):
		sys.stderr.write('Error : no input file\n')
		exit()

	for i in range(1, len(args)):
		data = np.loadtxt(args[i])
		tag = []
		x = []
		y = []
		z = []
		dens = []

		tag.extend(data[:,0])
		x.extend(data[:,3])
		y.extend(data[:,4])
		z.extend(data[:,5])
		dens.extend(data[:,15])

		x_tag = [[],[]]
		y_tag = [[],[]]
		z_tag = [[],[]]

		for i in range(len(x)):
			if tag[i] // (len(x)/2) == 0:
				x_tag[0].append(x[i])
				y_tag[0].append(y[i])
				z_tag[0].append(z[i])
			else:
				x_tag[1].append(x[i])
				y_tag[1].append(y[i])
				z_tag[1].append(z[i])

		pos_cg = np.array([0., 0., 0.])
		pos_cg = calc_cg(x, y, z, dens)
		#plot(x_tag, y_tag, z_tag, pos_cg)
		plot2(x_tag, y_tag, z_tag, pos_cg)


