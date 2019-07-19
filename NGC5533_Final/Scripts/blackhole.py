#Imports
#import numpy as np

#Constants
#Mbh = 2.7e9                                                 #mass of the central black hole in (solar mass)
#G = 4.300e-6                                                #gravitational constant (kpc/solar mass*(km/s)^2)

#equation for orbital velocity
def vbh(r):
    return np.sqrt((G*Mbh)/r)