#Imports
#import numpy as np
#import sympy as sym
#from scipy.integrate import quad
#from scipy.integrate import dblquad
#from scipy.special import ellipe
#from scipy.special import ellipk

#Constants
#G = 4.300e-6                                                #gravitational constant (kpc/solar mass*(km/s)^2)
#h = 8.9                                                     #radial scale-length (kpc)
#rho00 = 0.31e9                                              #prefactor that will cancel
#epsdisk = 5.0                                               #from Noordermeer's paper
#absmag = -22.02                                             #absolute magnitude 
#magsun = 4.42                                               #absolute magnitude of the sun
#z0 = 0.2*h                                                  #half-thickness (kpc)
#R = 4*h                                                     #cut-off radius (kpc)
#d = 0.2*h                                                   #cut-off length upper limits (kpc)
#L0 = np.power(10, (0.4*(magsun-absmag)))                    #Absolute Magnitude to luminosity

###########################################Intermediate Functions#####################################################

#Define Initial Functions
def x(r,u,xi):
    return ((r**2)+(u**2)+(xi**2))/(2*r*u)
def px(r,u,xi):
    return x(r,u,xi)-(np.sqrt((x(r,u,xi)**2)-1))

#Density Piecewise Function
def rho0(r, R, h, d):                                      
    condlist = [r <= R, (r > R) & (r <= (R+d)), r > (R+d)]
    funclist = [lambda r: rho00*np.exp(-r/h), lambda r: rho00*np.exp(-R/h)*(1-((r-R)/d)), lambda r: 0]
    return np.piecewise(r, condlist, funclist)
#Partial Derivative of rho(u,xi)
def durho0(r, R, h, d):                                  
    condlist = [r <= R, (r > R) & (r <= (R+d)), r > (R+d)]
    funclist = [lambda r: -(1/h)*rho00*np.exp(-r/h), lambda r: -(1/d)*rho00*np.exp(-R/h), lambda r: 0]
    return np.piecewise(r, condlist, funclist)

#Disk Density Distribution
rho_rz = lambda r,z: rho0(r, R, h, d)*(np.power(np.cosh(z/z0), (-2)))
drho_rz = lambda r,z: durho0(r, R, h, d)*(np.power(np.cosh(z/z0), -2))
#Complete Elliptic Integral
K = lambda r,u,xi: ellipk(px(r,u,xi)) - ellipe(px(r,u,xi))
#Inner Function (3D)
f = lambda r,u,z: u*drho_rz(u, z)*(2*K(r,u,z))/(np.pi*np.sqrt(r*u*px(r,u,z)))
#Integrate Function
f3 = lambda z,r,u: f(r,u,z)
intf = lambda r,u: quad(f3, 0, np.inf, args=(r,u,))[0]

#integrate outer function
intf3 = lambda u,r: intf(r,u)
intintf = lambda r: quad(intf3, 0.1, 125, args=(r,))[0]

#mass of disk
rho_rz_r = lambda z,r: rho_rz(r,z)*r
Mdblintrho = dblquad(rho_rz_r,0,125,-125,125)

#epsdisk = Mdblintrho/L0
pref = epsdisk*(L0/Mdblintrho)

#multiplying by epsylon
F = lambda r: 4*np.pi*G*intintf(r)*pref


######################################################################################################################
#Disk Velocity
rd = np.linspace(0.1, 125, num=100)
Fv = np.vectorize(F)                                  #save to file, open file, calculate for requested h if not already hdf5, h5py
vd = np.sqrt(-rd*Fv(rd))                            