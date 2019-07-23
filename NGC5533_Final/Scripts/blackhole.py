#Imports
#import numpy as np

#Constants
#Mbh = 2.7e9                                                 #mass of the central black hole in (solar mass)
#G = 4.300e-6                                                #gravitational constant (kpc/solar mass*(km/s)^2)

r0id = "varray_"+str(r0[0])+"-"+str(r0[len(r0)-1])+"_"+str(len(r0))+".hdf5"
location = "Inputs/"+r0id
Mbhval = "Mbh"+str(n)

try:
    saved = h5.File(location,"w")
except OSError:
    saved = h5.File(location,"r")

try:
    grp = saved.create_group("blackhole")
except ValueError:
    grp = saved["blackhole"]

try:     
    #Open
    vb = grp[Mbhval]
    #Read
    vbr = vb[:]
except KeyError:
    #equation for orbital velocity
    def vbh(r):
        return np.sqrt((G*Mbh)/r)
    
    #Calculate
    vbhr = np.zeros(len(r0))
    for i,n in enumerate(r0):
        vbhr[i] = vbh(n)
    
    #Save
    grp.create_dataset(Mbhval,data=vbhr)