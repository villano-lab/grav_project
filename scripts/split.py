#Imports
import sys
sys.path.append('/home/kitty/notebooks/grav_project/python')
import NGC5533_functions as nf
import h5py as h5

#Setup; change as needed
location = '../fitting/'
splitting = h5.File(location+'Inputs.hdf5')
groups = ['blackhole','halo','bulge','disk',]#'total']
          
#Split. Further reorganization may be desired.
for grp in groups:         #For every group,
    group = splitting[grp]
    for dataset in group:    #For every dataset in that group,
        nf.savedata(group[dataset][0],group[dataset] [1],grp,dataset,path=location,file=grp+'.hdf5')
        #Save the data in a new file

#Close
splitting.close()