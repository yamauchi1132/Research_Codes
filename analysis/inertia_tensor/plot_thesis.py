import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common import *

fsize = 14
data = np.loadtxt('result1.data')
l1 = []
l2 = []
l3 = []
time = []

l1.extend(data[:,0])
l2.extend(data[:,1])
l3.extend(data[:,2])
time.extend(data[:,3])

fig = plt.figure()

plt.plot(time, l1, label=r'$\hat{I}_1$',linewidth=2.0)
plt.plot(time, l2, label=r'$\hat{I}_2$',linewidth=2.0)
plt.plot(time, l3, label=r'$\hat{I}_3$',linewidth=2.0)

#plt.xscale('log')
#plt.yscale('log')

plt.xlabel(r'$t\ [s]$', fontsize=fsize)
plt.ylabel(r'$\rm{Inertia\ tensor}$', fontsize=fsize)
plt.xlim(0,1.2e+07)
plt.legend(fontsize=fsize)
plt.tick_params(labelsize=fsize)
mpl.rcParams['axes.xmargin'] = 0
mpl.rcParams['axes.ymargin'] = 0
plt.tight_layout()

plt.savefig("inertia_tensor.pdf", transparent=True, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()
plt.close()
