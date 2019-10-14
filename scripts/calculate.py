import ../python/NGC5533_Functions as nf

x = np.linspace(0,20,100)
rdat = np.linspace(0,125,100)

#For every pair of numbers in x,
for k in x:
    for j in x:
        #calculate and save the disk and bulge values.
        nf.d_v(rdat,k,j,save=True,path='../fitting')
        nf.b_v(rdat,k,j,save=True,path='../fitting')
        print(str(j)+':'+str(k)+' (maximum = 20:20)', end='')