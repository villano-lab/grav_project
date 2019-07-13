#Imports
#import numpy as np
#from scipy.integrate import quad
#from scipy.integrate import dblquad
#import scipy.optimize as so
#import scipy.special as ss

#Constants
n = 2.7                                     #concentration parameter that describes the curvature of the profile in a radius-magnitude plot, n=4 is de Vaucoileurs profile
re = 2.6                                    #1kpc
G = 4.300e-6                                #gravitational constant (kpc/solar mass*(km/s)^2)
ups = 2.8                                   #mass-to-light ratio (from Rotation Curves of Sersic Bulges paper)
q = 0.33                                    #intrinsic axis ratio
i = 45*(np.pi/180)                          #inclination angle
L = 3.27e10                                 #luminosity

#Variables
x = np.linspace(1, 10, 100)                        #x from/to and line smoothness

#Gamma Function
f = lambda x: ss.gammainc(2*n,x)*ss.gamma(2*n)-0.5*ss.gamma(2*n)
root = so.brentq(f,0,500000,rtol=0.000001,maxiter=100) #come within 1% of exact root within 100 iterations
r0 = re/root**n


########################################Intermediate Equations#########################################

f = lambda x,m: ((np.exp(-np.power(x/r0, (1/n))))*(np.power(x/r0, ((1/n)-1))))/(np.sqrt((x**2)-(m**2)));
g = lambda m: quad(f, m, np.inf,args=(m,))[0]
I0 = (L*(root**(2*n)))/(((re**2)*2*np.pi*n)*ss.gamma(2*n))
C = (4*G*q*ups*I0)/(r0*np.float(n))*(np.sqrt((np.sin(i)**2)+(1/(q**2))*(np.cos(i)**2)))
e2 = 1-(q**2)

#integrate outer function
h = lambda m,r: C*g(m)*(m**2)/(np.sqrt((r**2)-((m**2)*(e2))))

y = np.zeros(np.shape(x))
for j,r in enumerate(x):
    y[j] = quad(h, 0, r,args=(r,))[0]

#######################################################################################################

#Final Equation
vb = np.sqrt(y)

######################################## Save vb for given n ##########################################

nval = "n"+str(n)

try:
    saved = h5.File("inputs.hdf5","w")
except OSError:
    saved = h5.File("inputs.hdf5","r")
    
try:
    grp = f.create_group("bulge")
except RuntimeError:
    grp = f["bulge"]

try:                                                
    dset = grp.create_dataset(nval,vb,dtype='a')
except ValueError:
    dset = grp[nval]