#Imports
#import numpy as np

#Constants
#G = 4.300e-6                                                #gravitational constant (kpc/solar mass*(km/s)^2)
#rho00 = 0.31e9                                              #central surface density (solar mass/kpc^3)

#using a velocity equation from the paper "Dark halo properties from rotation curves" by Jimenez et.al.
vcdm = lambda r: np.sqrt(4*np.pi*G*rho00*(rc**2)*(1-((rc/r)*np.arctan(r/rc))))       #eq 9 from Jimenez paper #rho00
