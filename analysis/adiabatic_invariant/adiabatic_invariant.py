import sys, os
sys.path.append(os.pardir)
import numpy as np
import operator
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from common2 import *

dirname = "../data/sph_t%04d.dat"
start = 600
end = 800
step = 1

Msun = 1.989e+33