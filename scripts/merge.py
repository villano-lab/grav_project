import sys
#Change this based on your local machine's path to the repo
sys.path.append('/home/gerudo7/git/grav_project/python')
import NGC5533_functions as nf
import h5py as h5

#Change the below to whatever file you want to merge from:
fromfile = h5.File('../fitting/Static_Notebooks/Inputs.hdf5')
#Change the below to whatever file you want to merge to:
tofile = '../fitting/Inputs.hdf5'

#If you only want to keep certain groups or datasets, you'll need to change this
#Keep this in mind if you're merging because of corrupt old data.
for group in fromfile:             #for every group,
    for dataset in fromfile[group]:#for every dataset within that group,
        #save that data under the same group and dataset name.
        #print(dataset) #Troubleshooting; lots of output. ("verbose")
        dset = fromfile[group][dataset]
        #print(dset) #Troubleshooting; lots of output. ("verbose")
        xdata = dset[0]
        ydata = dset[1]
        nf.savedata(xdata,ydata,group[1:],dataset[1:],file=tofile)
fromfile.close()