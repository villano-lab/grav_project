print("========================= Running Scripts/halo.py ==========================")
print()

####################################### Setup. May be done in parent code. ############################################
"""
#Imports
print("Importing libraries...")
import numpy as np
import h5py as h5

#Constants
print("Setting constants...")
G = 4.300e-6                                                #gravitational constant (kpc/solar mass*(km/s)^2)
rho00 = 0.31e9                                              #central surface density (solar mass/kpc^3)
rc = 1.4                                                    #core radius (kpc)

#"Variables"
print("Setting arrays...")
r0 = np.linspace(0.1, 10, 100)                              #radius over which plotting will be done

#Savedata name values
print("Setting name values...")
startval = str(r0[0])
endval = str(r0[len(r0)-1])
r0id = "varray_"+startval.replace('.','_')+"-"+endval.replace('.','_')+"_"+str(len(r0))+".hdf5"
location = "Inputs/"+r0id

#Load savedata if available, else create file
print("Setting up savedata file '"+r0id+"'...")
saved = h5.File(location,'a')
"""
######################################## Save vcdm for given rh00 #########################################################

r00val = "r00"+str(rho00)

try:
    grp = saved.create_group("halo")
    print("Created group 'halo'.")
except ValueError:
    grp = saved["halo"]
    print("Loaded group 'halo'.")

try:     
    #Open
    vcdmd = grp[r00val]
    #Read
    vcdmr = vcdmd[:]
    print("Loaded savedata.")
except KeyError:
    #using a velocity equation from the paper "Dark halo properties from rotation curves" by Jimenez et.al.
    print("Defining functions...")
    vcdm = lambda r: np.sqrt(4*np.pi*G*rho00*(rc**2)*(1-((rc/r)*np.arctan(r/rc))))       #eq 9 from Jimenez paper #rho00
    
    #Calculate
    print("Calculating...")
    vcdmr = np.zeros(len(r0))
    for i,n in enumerate(r0):
        vcdmr[i] = vcdm(n)
    
    #Save
    print("Saving dataset '"+r00val+"'...")
    grp.create_dataset(r00val,data=vcdmr)
    
#Close file. Comment out if opening file elsewhere.
#print("Closing file '"+r0id+"'...")
#saved.close()
    
print()
print("Finished running script for disk velocity.")
print(":)")
