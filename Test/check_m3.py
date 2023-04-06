import numpy as np
import math
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate


def s2k(E, k, m):
    return (E-math.sqrt(k**2 + m**2))**2 - k**2


fig,ax = plt.subplots(figsize=(5,5))

s3 = 9.1
E3 = math.sqrt(s3)

kini = 0.0
kfin = 10.0

delk = 0.01

number = int( abs(kfin - kini)/delk)

m = 1

s2kres = []
kres = []
for i in range(number):
    k = kini + i*delk
    kres.append(k)
    s2kres.append(s2k(E3,k,m))

ax.scatter(kres,s2kres)

plt.show()

