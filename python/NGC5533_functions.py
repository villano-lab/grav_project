import numpy as np
#---------
G = 4.30091e-6  #gravitational constant (kpc/solar mass*(km/s)^2)

def bh(r,M=2.7e9):

  #M in solar masses
  #r in kpc

  return np.sqrt(G*M/r)
