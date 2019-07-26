print("========================= Running Scripts/bulge.py ==========================")
print()

####################################### Setup. May be done in parent code. ############################################
"""
#Imports
print("Importing libraries...")
import numpy as np
import h5py as h5
from scipy.integrate import quad
from scipy.integrate import dblquad
import scipy.optimize as so
import scipy.special as ss

#Constants
print("Setting constants...")
n = 2.7                                     #concentration parameter that describes the curvature of the profile in a radius-magnitude plot, n=4 is de Vaucoileurs profile
re = 2.6                                    #1kpc
G = 4.300e-6                                #gravitational constant (kpc/solar mass*(km/s)^2)
ups = 2.8                                   #mass-to-light ratio (from Rotation Curves of Sersic Bulges paper)
q = 0.33                                    #intrinsic axis ratio
i = 45*(np.pi/180)                          #inclination angle
L = 3.27e10                                 #luminosity

#Variables
print("Setting plotting arrays...")
xd = np.linspace(1, 10, 100)                        #x from/to and line smoothness
r0 = np.linspace(0.1, 10, 100)                       #radius over which plotting will be done

#Savedata naming conventions
print("Setting name values...")
startval = str(r0[0])
endval = str(r0[len(r0)-1])
r0id = "varray_"+startval.replace('.','_')+"-"+endval.replace('.','_')+"_"+str(len(r0))+".hdf5"
location = "Inputs/"+r0id

#Load savedata if available, else create file
print("Setting up savedata file '"+r0id+"'...")
saved = h5.File(location,'a')
"""
######################################## Save vb for given n ##########################################

nval = "n"+str(n)

try:
    grp = saved.create_group("bulge")
    print("Created group 'bulge'.")
except ValueError:
    grp = saved["bulge"]
    print("Loaded group 'bulge'.")

try:     
    #Open
    vbd = grp[nval]
    #Read
    vbr = vbd[:]
    print("Loaded savedata.")
except KeyError:
    #Gamma Function
    print("Defining gamma function...")
    f1 = lambda x: ss.gammainc(2*n,x)*ss.gamma(2*n)-0.5*ss.gamma(2*n)
    root = so.brentq(f1,0,500000,rtol=0.000001,maxiter=100) #come within 1% of exact root within 100 iterations
    ra = re/root**n

    #Functions
    print("Defining velocity function...")
    def vb(r):
        fb = lambda x,m: ((np.exp(-np.power(x/ra, (1/n))))*(np.power(x/ra, ((1/n)-1))))/(np.sqrt((x**2)-(m**2)));
        g = lambda m: quad(fb, m, np.inf,args=(m,))[0]
        I0 = (L*(root**(2*n)))/(((re**2)*2*np.pi*n)*ss.gamma(2*n))
        C = (4*G*q*ups*I0)/(ra*np.float(n))*(np.sqrt((np.sin(i)**2)+(1/(q**2))*(np.cos(i)**2)))
        e2 = 1-(q**2)
        H = lambda m,r: C*g(m)*(m**2)/(np.sqrt((r**2)-((m**2)*(e2))))
        y = lambda r: quad(H,0,r,args=(r,))[0]
        return y(r)**0.5
    
    #Calculate
    print("Calculating... This may take a few minutes.")
    vbr = np.zeros(len(r0))
    for i,n in enumerate(r0):
        vbr[i] = vb(n)
    
    #Save
    print("Saving dataset '"+nval+"'...")
    grp.create_dataset(nval,data=vbr)

#Close file. Comment out if opening file elsewhere.
#print("Closing file '"+r0id+"'...")
#saved.close()    
    
print()
print("Finished running script for bulge velocity.")
print(":)")
