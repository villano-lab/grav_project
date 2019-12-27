# Traced curves from Noordermeer's paper: Rotation curves of flattened Sersic bulges, Figure 4.

################################
########### Imports ############
################################
import sys
sys.path.append('../../python/')

import dataPython as dp
import numpy as np
import scipy.interpolate as inter
import matplotlib.pyplot as plt

################################
######### Data files ###########
################################
data_total = dp.getXYdata('../data/final/nord-120kpc-total.txt')
data_bh = dp.getXYdata('../data/final/nord-120kpc-blackhole.txt')
data_bulge = dp.getXYdata('../data/final/nord-120kpc-bulge.txt')
data_disk = dp.getXYdata('../data/final/nord-120kpc-disk.txt')
data_halo = dp.getXYdata('../data/final/nord-120kpc-halo.txt')
data_gas = dp.getXYdata('../data/final/nord-120kpc-gas.txt')
data_greyb_bottom = dp.getXYdata('../data/final/nord-120kpc-bottomband.txt')
data_greyb_top = dp.getXYdata('../data/final/nord-120kpc-topband.txt')

################################
##### Measured data points #####
################################
data = dp.getXYdata_wXYerr('../data/100kpc_data.txt')
r_dat = np.asarray(data['xx'])
v_dat = np.asarray(data['yy'])
v_err0 = np.asarray(data['ex'])
v_err1 = np.asarray(data['ey'])

################################
####### Uncertainty band #######
################################
r_bottomband = np.asarray(data_greyb_bottom['xx'])
v_bottomband = np.asarray(data_greyb_bottom['yy'])
r_topband = np.asarray(data_greyb_top['xx'])
v_topband = np.asarray(data_greyb_top['yy'])

noord_greyb_bottom = inter.spline(r_bottomband,v_bottomband,rval,kind='smoothest')
noord_greyb_top = inter.spline(r_topband,v_topband,rval,kind='smoothest')

################################
######### Total curve ##########
################################
r_total = np.asarray(data_total['xx'])
v_total = np.asarray(data_total['yy'])

noord_total = inter.spline(r_total,v_total,rval,kind='smoothest')

################################
######### Black Hole ###########
################################
r_bh = np.asarray(data_bh['xx'])
v_bh = np.asarray(data_bh['yy'])

noord_blackhole = inter.spline(r_bh,v_bh,rval,kind='smoothest')

################################
############ Bulge #############
################################
r_bulge = np.asarray(data_bulge['xx'])
v_bulge = np.asarray(data_bulge['yy'])

noord_bulge = inter.spline(r_bulge,v_bulge,rval,kind='smoothest')

################################
############ Disk ##############
################################
r_disk = np.asarray(data_disk['xx'])
v_disk = np.asarray(data_disk['yy'])

noord_disk = inter.spline(r_disk,v_disk,rval,kind='smoothest')

################################
############ Halo ##############
################################
r_halo = np.asarray(data_halo['xx'])
v_halo = np.asarray(data_halo['yy'])

noord_halo = inter.spline(r_halo,v_halo,rval,kind='smoothest')

################################
############# Gas ##############
################################
r_gas = np.asarray(data_gas['xx'])
v_gas = np.asarray(data_gas['yy'])

noord_gas = inter.spline(r_gas,v_gas,rval,kind='smoothest')