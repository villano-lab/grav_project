import sys
sys.path.append('/home/gerudo7/git_repositories/grav_project/python/')
import NGC5533_functions as nf
import numpy as np

x = np.linspace(0.01,20,100)
rdat = np.linspace(0,125,100)

#For every pair of numbers in x,
for k in x:
    for j in x:
        #calculate and save the disk and bulge values.
        nf.d_v(rdat,k,j,load=True,path='../fitting/')
        nf.b_v(rdat,k,j,load=True,path='../fitting/')
        print(str(j)+':'+str(k)+' (maximum = 20:20)  ')
