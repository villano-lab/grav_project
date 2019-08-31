import numpy as np
import scipy.special as ss
import scipy.optimize as so
import scipy.integrate as si

#---------Definitely Constant---------
G = 4.30091e-6    #gravitational constant (kpc/solar mass*(km/s)^2)
rhocrit = 9.3e-18 #critical density of the Universe (kg/km^3)

#---------Measured Directly-----------
L = 3.27e10                              #luminosity
absmag = -22.02                          #absolute magnitude
magsun = 4.42                            #absolute magnitude of the sun
L0 = np.power(10, (0.4*(magsun-absmag))) #Absolute Magnitude to luminosity

#---------Measured Indirectly---------
ups = 2.8                         #mass-to-light ratio (from Rotation Curves of Sersic Bulges paper)
q = 0.33                          #intrinsic axis ratio
i = 45*(np.pi/180)                #inclination angle
h = 8.9                           #radial scale-length (kpc)
z0 = 0.2*h                        #half-thickness (kpc)
R = 4*h                           #cut-off radius (kpc)
d = 0.2*h                         #cut-off length upper limits (kpc)
h_rc = 1.4                        #core radius (kpc)
h_rho00 = 0.31e9                  #central surface density (solar mass/kpc^3)
c = 1e-12                         #concentration parameter
Mvir = 1e11*((c/(11.7))**(-40/3)) #virial mass (in solar mass) solved from eq(5)

#---------Definitely Variable---------
#n = 2.7

#---------Uncertain-------------------
re = 2.6                                                 #1kpc
d_rho00 = 0.31e9                                         #prefactor that will cancel #central density
epsdisk = 5.0                                            #from Noordermeer's paper
rs = (1/c)*(((3*Mvir)/((4*np.pi*100*rhocrit)))**(1/3))   #scale radius (kpc)
rho_s = (100/3)*((c**3)/(np.log(1+c)-(c/(1+c))))*rhocrit #characteristic density
h_gamma = 0

################################
######### Black Hole ###########
################################

def bh_v(r,M=2.7e9): #M in solar masses, r in kpc
  return np.sqrt(G*M/r)

################################
########### Bulge ##############
################################

#I'm not sure how many of these we need to be defined -- I kept everything that was called outside of another function.
#We can condense the number of functions once we know for certain if there are things we don't need again.

def b_gammafunc(x):
    f = ss.gammainc(2*n,x)*ss.gamma(2*n)-0.5*ss.gamma(2*n)
    return so.brentq(f,0,500000,rtol=0.000001,maxiter=100) #come within 1% of exact root within 100 iterations

def b_innerintegral(m):
    f = lambda x,m: ((np.exp(-np.power(x/r0, (1/n))))*(np.power(x/r0, ((1/n)-1))))/(np.sqrt((x**2)-(m**2))); #Inner function
    return si.quad(f, m, np.inf,args=(m,))[0]
b_innerintegralv = np.vectorize(b_innerintegral)

def b_I0(n=2.7):
    return L*(b_gammafunc**(2*n))/(((re**2)*2*np.pi*n)*ss.gamma(2*n))
def b_r0(n=2.7):
    return re/b_gammafunc**n

def b_vsquare(r,n=2.7):
    C = (4*G*q*ups*I0)/(r0*np.float(n))*(np.sqrt((np.sin(i)**2)+(1/(q**2))*(np.cos(i)**2)))
    e2 = 1-(q**2)
    h = lambda m,r: C*g(m)*(m**2)/(np.sqrt((r**2)-((m**2)*(e2)))) #integrate outer function
    y = si.quad(h, 0, r, args=(r,n))[0]
    return np.vectorize(y)
def b_v(r,n=2.7):
    return b_vsquare**(1/2)

################################
############ Halo ##############
################################

def h_rhat(r,z):               #r-hat from Casertano eq(9)
    return np.sqrt((r**2)+(z**2))

def h_rho(r,rho00=h_rho00,rc=h_rc): #Isothermal Density Profile
    return rho00*((1+((r/rc)**2))**(-1))
    
def h_vcasertano(r,z,rc=h_rc,rho00=h_rho00,gamma=h_gamma):                       #Velocity
    v0h = lambda r,rho00,rc,z: np.sqrt(h_rho(r,rho00,rc)*4*np.pi*G*(h_rhat(r,z)**2)) #eq 9 casertano
    return v0h(r,rho00,rc,z)*((r/rc)**gamma)                                     #eq 10 casertano

def h_vjimenez(r,rc=h_rc,rho00=h_rho00):
    return np.sqrt(4*np.pi*G*rho00*(rc**2)*(1-((rc/r)*np.arctan(r/rc))))

def h_v(r):
    rho = lambda r: rho_s/((r/rs)*((1+r/rs)**2))
    f = lambda R: 4*np.pi*rho(R)*(R**2)          #NFW Density Profile
    mdm = lambda r: si.quad(f, 0, r)[0]             #M(r)
    vdm2 = lambda r: (G*mdm(r))/r                #v^2: GM(r)/r
    vdm2v = np.vectorize(vdm2)
    return np.sqrt(vdm2v(r))

################################
############ Disk ##############
################################

def d_px(r,u,xi):       #Initial Function
    x = ((r**2)+(u**2)+(xi**2))/(2*r*u)
    return x-(np.sqrt((x(r,u,xi)**2)-1))

def d_rho0(r, R, h, d): #density piecewise function
    condlist = [r <= R, (r > R) & (r <= (R+d)), r > (R+d)]
    funclist = [lambda r: d_rho00*np.exp(-r/h), lambda r: d_rho00*np.exp(-R/h)*(1-((r-R)/d)), lambda r: 0]
    return np.piecewise(r, condlist, funclist)

def d_durho0(r, R, h, d): #partial derivative of rho(u,xi)
    condlist = [r <= R, (r > R) & (r <= (R+d)), r > (R+d)]
    funclist = [lambda r: -(1/h)*d_rho00*np.exp(-r/h), lambda r: -(1/d)*d_rho00*np.exp(-R/h), lambda r: 0]
    return np.piecewise(r, condlist, funclist)

#Disk Density Distribution
def d_rho_rz(r,z):
    return d_rho0(r, R, h, d)*(np.power((np.cosh(z/z0)), (-2)))
def d_drho_rz(r,z):
    return d_durho0(r, R, h, d)*(np.power(np.cosh(z/z0), -2))

#Disk Density Distribution (3D)
def d_rho_rz3D(r, z):
    return d_rho0(r, R, h, d)*(np.power(np.cosh(z/z0), (-2)))
def d_drho_rz3D(r, z):
    return d_durho0(r, R, h, d)*(np.power(np.cosh(z/z0), (-2)))

def d_K(r,u,xi): #Complete Elliptic Integral
    return ss.ellipk(px(r,u,xi)) - ss.ellipe(px(r,u,xi))

def d_innerfunc(z,r,u):  #Inner Function (3D)
    return u*d_drho_rz(u, z)*(2*K(r,u,z))/(np.pi*np.sqrt(r*u*px(r,u,z)))

def d_innerintegral(u,r): #Integrate Function
    return si.quad(f, 0, np.inf, args=(r,u,))[0]

def d_outerintegral(r): #Integrate Outer Function
    return si.nquad(intf, [[0.1, 125]], args=(r,),opts=[options,options])[0]

def d_Mintrho(r):
    rho_rz_r = lambda z,r: d_rho_rz(r,z)*r
    return si.quad(rho_rz_r, -125, 125, args=(r,))[0]
d_Mdblintrho = si.quad(d_Mintrho,0,125)[0]

def d_Fv(r): #multiplying by upsylon
    pref = epsdisk*(l0/Mdlbintrho)
    F = 4*np.pi*G*d_outerintegral(r)*pref
    return np.vectorize(F)   

def d_v(r): #velocity
    return np.sqrt(-r*Fv(r))