print("========================= Running Scripts/blackhole.py ==========================")
print()

####################################### Setup. May be done in parent code. ############################################
"""
#Imports
print("Importing libraries...")
import numpy as np
import h5py as h5

#Constants
print("Setting constants...")
Mbh = 2.7e9                                                 #mass of the central black hole in (solar mass)
G = 4.300e-6                                                #gravitational constant (kpc/solar mass*(km/s)^2)

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
######################################## Save vbh for given Mbh ##########################################

Mbhval = "Mbh"+str(Mbh)

try:
    grp = saved.create_group("blackhole")
    print("Created group 'blackhole'")
except ValueError:
    grp = saved["blackhole"]
    print("Loaded group 'blackhole'")

try:     
    #Open
    vbh = grp[Mbhval]
    #Read
    vbhr = vbh[:]
    print("Loaded savedata.")
except KeyError:
    #equation for orbital velocity
    print("Defining functions...")
    def vbh(r):
        return np.sqrt((G*Mbh)/r)
    
    #Calculate
    print("Calculating...")
    vbhr = np.zeros(len(r0))
    for i,n in enumerate(r0):
        vbhr[i] = vbh(n)
    
    #Save
    print("Saving dataset '"+Mbhval+"'...")
    grp.create_dataset(Mbhval,data=vbhr)
    
#Close file. Comment out if opening file elsewhere.
#print("Closing file '"+r0id+"'...")
#saved.close()
    
print()
print("Finished running script for disk velocity.")
print(":)")