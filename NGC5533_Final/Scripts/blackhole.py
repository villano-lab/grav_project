#Imports
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

#Constants
Mbh = 2.7e9                                                 #mass of the central black hole in (solar mass)
G = 4.300e-6                                                #gravitational constant (kpc/solar mass*(km/s)^2)

#equation for orbital velocity
def vbh(r, G, Mbh):
    return np.sqrt((G*Mbh)/r)