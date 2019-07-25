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

######################################## Save vd for given h #########################################################

#Savedata name values
#r0id = "varray_"+str(r0[0])+"-"+str(r0[len(r0)-1])+"_"+str(len(r0))+".hdf5"
#location = "Inputs/"+r0id
hval = "h"+str(h)

#Load savedata if available, else create file
#try:
#    saved = h5.File(location,"w")
#except OSError:
#    saved = h5.File(location,"r")
    
try:
    grp = saved.create_group("disk")
except ValueError:
    grp = saved["disk"]

try:   
    #Open
    vdd = grp[hval]
    #Read
    vdr = vdd[:]
except KeyError:                                                                                                    #calculate vd and save it for hval
    #Equations
    def x0(r,u,xi):
        return ((r**2)+(u**2)+(xi**2))/(2*r*u)
    def px(r,u,xi):                                                                                                 #Define Initial Functions
        return x0(r,u,xi)-(np.sqrt((x0(r,u,xi)**2)-1))        
    def rho0(r, R, h, d):                                                                                           #Density Piecewise Function                                   
        condlist = [r <= R, (r > R) & (r <= (R+d)), r > (R+d)]
        funclist = [lambda r: rho00*np.exp(-r/h), lambda r: rho00*np.exp(-R/h)*(1-((r-R)/d)), lambda r: 0]
        return np.piecewise(r, condlist, funclist)
    def durho0(r, R, h, d):                                                                                         #Partial Derivative of rho(u,xi)                      
        condlist = [r <= R, (r > R) & (r <= (R+d)), r > (R+d)]
        funclist = [lambda r: -(1/h)*rho00*np.exp(-r/h), lambda r: -(1/d)*rho00*np.exp(-R/h), lambda r: 0]
        return np.piecewise(r, condlist, funclist)
    rho_rz = lambda r,z: rho0(r, R, h, d)*(np.power(sym.cosh(z/z0).evalf(), (-2)))                                           #Disk Density Distribution
    drho_rz = lambda r,z: durho0(r, R, h, d)*(np.power(np.cosh(z/z0), -2))
    K = lambda r,u,xi: ellipk(px(r,u,xi)) - ellipe(px(r,u,xi))                                                      #Complete Elliptic Integral
    f = lambda z,r,u: u*drho_rz(u, z)*(2*K(r,u,z))/(np.pi*np.sqrt(r*u*px(r,u,z)))                                   #Inner Function (3D)                                                                                 
    intf = lambda u,r: quad(f, 0, np.inf, args=(r,u,))[0]                                                           #Integrate Function
    intintf = lambda r: nquad(intf, [[0.1, 125]], args=(r,),opts=[options,options])[0]                              #integrate outer function
    rho_rz_r = lambda z,r: rho_rz(r,z)*r                                                                            #mass of disk
    Mintrho = lambda r: quad(rho_rz_r, -125, 125, args=(r,))[0]
    Mdblintrho = quad(Mintrho, 0, 125)[0]                                                                           #epsdisk = Mdblintrho/L0
    pref = epsdisk*(L0/Mdblintrho)                                                                                  #multiplying by epsylon
    F = lambda r: 4*np.pi*G*intintf(r)*pref
    
    #Final Function
    #rd = np.linspace(0.1, 125, num=100)                                                                             
    #Fv = np.vectorize(F)                                                                                            #Disk Velocity
    vd = lambda r: np.sqrt(-r*F(r))
    
    #Calculate
    vdr = np.zeros(len(r0))
    for i,n in enumerate(r0):
        vdr[i] = vd(n)
    vdr[np.isnan(vdr)] = 0                                                                                          #Set all nan to 0
    
    #Save
    grp.create_dataset(hval,data=vdr)                                                                                #save to file, open file, calculate for requested h if not already hdf5, h5py