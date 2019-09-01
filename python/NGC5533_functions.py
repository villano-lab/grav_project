################################
########### Imports ############
################################

import numpy as np
import scipy.special as ss
import scipy.optimize as so
import scipy.integrate as si
try:
    import h5py as h5
    h5py = 1
except:
    h5py = 0
    print("Unable to load h5py. Datasets will not be able to be saved or loaded using NGC5533_functions.")

################################
########## Settings ############
################################    

options={'limit':100}     #Sets maximum subdivisions to 100 for integration instead of 50
x = np.linspace(0,10,100) #Default radius array

################################
########## Constants ###########
################################    

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
e2 = 1-(q**2)
i = 45*(np.pi/180)                #inclination angle
h = 8.9                           #radial scale-length (kpc)
z0 = 0.2*h                        #half-thickness (kpc)
R = 4*h                           #cut-off radius (kpc)
d = 0.2*h                         #cut-off length upper limits (kpc)
h_rc = 1.4                        #core radius (kpc)
h_rho00 = 0.31e9                  #central surface density (solar mass/kpc^3)
c = 1e-12                         #concentration parameter
Mvir = 1e11*((c/(11.7))**(-40/3)) #virial mass (in solar mass) solved from eq(5)
Mbh_def = 2.7e9                   #Black Hole mass

#---------Definitely Variable---------
n_c = 2.7

#---------Uncategorized-------------------
re = 2.6                                                 #1kpc
d_rho00 = 0.31e9                                         #prefactor that will cancel #central density
epsdisk = 5.0                                            #from Noordermeer's paper
rs = (1/c)*(((3*Mvir)/((4*np.pi*100*rhocrit)))**(1/3))   #scale radius (kpc)
rho_s = (100/3)*((c**3)/(np.log(1+c)-(c/(1+c))))*rhocrit #characteristic density
h_gamma = 0

################################
########### Saving #############
################################

def save(values,dataset,group,path='Inputs',file='x'+str(x[0])+'-'+str(x[len(x)-1])+'_'+str(len(x))):
    if h5py == 1:
        saved = h5.File(path+file)
        if group in ['disk', 'Disk', 'disc', 'Disc', 'd', 'D']:
            group = 'disk'
            print("Group name set to 'disk'.")
        if group in ['bh','Bh','BH','Black Hole','BlackHole','Blackhole,','blackhole','Black hole','black hole','Black Hole']:
            group = 'blackhole'
            print("Group name set to 'blackhole'.")
        if group in ['dm','DM','Dm','Dark Matter','Dark matter','dark matter','halo','h','H','Halo','darkmatter','Darkmatter','Dark Matter']:
            group = 'halo'
            print("Group name set to 'halo'.")
        if group in ['bulge','b','B','Bulge']:
            group = 'bulge'
            print("Group name set to 'bulge'.")
        if group in ['t','T','total','Total']:
            group = 'total'
            print("Group name set to 'total'.")
        try:
            grp = saved.create_group(group)
            grp.create_dataset(dataset,data=values)
        except KeyError:
            grp = saved[group]
            grp.create_dataset(dataset,data=values)
        saved.close()
        print("Saved.")
    #Placeholder; I will design this to store information at a later date.
    if h5py ==0:
        print("ERROR: h5py was not loaded.")
        return 1
    
def load(dataset,group,path='Inputs',file='x'+str(x[0])+'-'+str(x[len(x)-1])+'_'+str(len(x))):
    if h5py == 1:
        saved = h5.File(path+file)
        grp = saved[group]
        dset = grp[dataset]
        return dset[:]
        saved.close
    #Placeholder; I will design this to store information at a later date.
    if h5py ==0:
        print("ERROR: h5py was not loaded.")
        return 1

################################
######### Black Hole ###########
################################

def bh_v(r,M=Mbh_def,save=True): #M in solar masses, r in kpc
    return np.sqrt(G*M/r)
    if save is True:
        a = np.sqrt(G*M/x)
        save(a,'Mbh'+str(M),'blackhole')

################################
########### Bulge ##############
################################

#I'm not sure how many of these we need to be defined -- I kept everything that was called outside of another function.
#We can condense the number of functions once we know for certain if there are things we don't need again.

def b_gammafunc(x,n=n_c):
    return ss.gammainc(2*n,x)*ss.gamma(2*n)-0.5*ss.gamma(2*n)
b_root = so.brentq(b_gammafunc,0,500000,rtol=0.000001,maxiter=100) #come within 1% of exact root within 100 iterations

def b_I0(n=n_c):
    return L*(b_root**(2*n))/(re**2*2*np.pi*n*ss.gamma(2*n))
def b_r0(n=n_c):
    return re/np.power(b_root,n)

def b_innerintegral(m,n=n_c):
    f = lambda x,m,n: np.exp(-np.power(x/b_r0(n), (1/n)))*np.power(x/b_r0(n), 1/n-1)/(np.sqrt(x**2-m**2)) #Inner function
    return si.quad(f, m, np.inf,args=(m,n))[0]
b_innerintegralv = np.vectorize(b_innerintegral)

def b_vsquare(r,n=n_c):
    C = lambda n: (4*G*q*ups*b_I0(n))/(b_r0(n)*np.float(n))*(np.sqrt((np.sin(i)**2)+(1/(q**2))*(np.cos(i)**2)))
    h = lambda m,r,n: C(n)*b_innerintegral(m,n)*(m**2)/(np.sqrt((r**2)-((m**2)*(e2)))) #integrate outer function
    return si.quad(h, 0, r, args=(r,n))[0]
def b_vsquarev(r,n=n_c):
    return np.vectorize(b_vsquare, otypes=[np.float])
def b_v(r,n=n_c):
    return b_vsquare(r,n)**(1/2)
    #unvec = lambda r,n: b_vsquare(r,n)**(1/2)
    #return np.vectorize(unvec, otypes=[np.float])

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

def h_v(r,save=True):
    rho = lambda r: rho_s/((r/rs)*((1+r/rs)**2))
    f = lambda R: 4*np.pi*rho(R)*(R**2)          #NFW Density Profile
    mdm = lambda r: si.quad(f, 0, r)[0]          #M(r)
    vdm2 = lambda r: (G*mdm(r))/r                #v^2: GM(r)/r
    vdm2v = np.vectorize(vdm2)
    return np.sqrt(vdm2v(r))
    if save is True:
        savedata()

################################
############ Disk ##############
################################

def d_px(r,u,xi):       #Initial Function
    x = lambda r,u,xi: ((r**2)+(u**2)+(xi**2))/(2*r*u)
    return x(r,u,xi)-(np.sqrt((x(r,u,xi)**2)-1))

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
    return ss.ellipk(d_px(r,u,xi)) - ss.ellipe(d_px(r,u,xi))

def d_innerfunc(z,r,u):  #Inner Function (3D)
    return u*d_drho_rz(u, z)*(2*d_K(r,u,z))/(np.pi*np.sqrt(r*u*d_px(r,u,z)))

def d_innerintegral(u,r): #Integrate Function
    return si.quad(d_innerfunc, 0, np.inf, args=(r,u,))[0]

def d_outerintegral(r): #Integrate Outer Function
    return si.nquad(d_innerintegral, [[0.1, 125]], args=(r,),opts=[options,options])[0]

def d_Mintrho(r):
    rho_rz_r = lambda z,r: d_rho_rz(r,z)*r
    return si.quad(rho_rz_r, -125, 125, args=(r,))[0]
d_Mdblintrho = si.quad(d_Mintrho,0,125)[0]

def d_F(r): #multiplying by upsylon
    pref = epsdisk*(L0/d_Mdblintrho)
    return 4*np.pi*G*d_outerintegral(r)*pref
d_Fv = np.vectorize(d_F)

def d_v(r): #velocity
    return np.sqrt(-r*d_Fv(r))

################################
############ Total #############
################################
def v(r,M=Mbh_def,n=n_c): 
    return np.sqrt(h_v(r)**2+d_v(r)*d_v(r)+bh_v(r,M)**2+b_v(r,n)**2)