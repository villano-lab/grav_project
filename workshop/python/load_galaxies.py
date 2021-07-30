#Data to fit to for each galaxy to be used in workshop

###############################
########## Imports ############
###############################
import dataPython         as dp
import numpy              as np
import scipy.interploate  as inter
import matplotlib.pylplot as plt

###############################
########### NGC5533 ###########
###############################

NGC5533 = {

    # Load data from files for Noordermeer's band and fitted curves
    # 'raw' in the sense that I haven't split everything up. Idk there's probably a better way to name these
    'raw_total'        : dp.getXYdata('data/NGC5533/noord-120kpc-total.txt'     );
    'raw_blackhole'    : dp.getXYdata('data/NGC5533/noord-120kpc-blackhole.txt' );
    'raw_bulge'        : dp.getXYdata('data/NGC5533/noord-120kpc-bulge.txt'     );
    'raw_disk'         : dp.getXYdata('data/NGC5533/noord-120kpc-disk.txt'      );
    'raw_halo'         : dp.getXYdata('data/NGC5533/noord-120kpc-halo.txt'      );
    'raw_gas'          : dp.getXYdata('data/NGC5533/noord-120kpc-gas.txt'       );
    'raw_band_btm' : dp.getXYdata('data/NGC5533/noord-120kpc-bottomband.txt');
    'raw_band_top' : dp.getXYdata('data/NGC5533/noord-120kpc-topband.txt'   );

    # Get data from 100kpc file
    'measured_data': dp.getXYdata('data/100kpc_data.txt');
        
}
    
#Organize 100kpc data
NGC5533['m_radii']      = np.asarray(NGC5533['measured_data']['xx'])
NGC5533['m_velocities'] = np.asarray(NGC5533['measured_data']['yy'])
NGC5533['m_r_errors']   = np.asarray(NGC5533['measured_data']['ex'])
NGC5533['m_v_errors']   = np.asarray(NGC5533['measured_data']['ey'])

#Organize band data
NGC5533['n_r_btmband']   = np.asarray(NGC5533['raw_band_bttm']['xx'])
NGC5533['n_v_btmband']   = np.asarray(NGC5533['raw_band_bttm']['yy'])
NGC5533['n_r_topband']   = np.asarray(NGC5533['raw_band_top']['xx'])
NGC5533['n_v_topband']   = np.asarray(NGC5533['raw_band_top']['yy'])
NGC5533['n_v_bandwidth'] = NGC5533['n_v_topband'] - NGC5533['n_v_btmband']
NGC5533['n_v_bandwidth'] = NGC5533['n_v_bandwidth'][0::28] #For weights, v_errors and band must line up.
NGC5533['n_v_bandwidth'] = NGC5533['n_v_bandwidth'][1:]
    
# Smoothing
NGC5533['n_tb'], NGC5533['n_cb'], NGC5533['n_kb'] = inter.splrep(NGC5533['n_r_btmband'],NGC5533['n_v_btmband'])
NGC5533['n_tt'], NGC5533['n_ct'], NGC5533['n_kt'] = inter.splrep(NGC5533['n_r_topband'],NGC5533['n_v_topband'])
NGC5533['n_band_btm'] = inter.BSpline(NGC5533['n_tb'], NGC5533['n_cb'], NGC5533['n_kb'])
NGC5533['n_band_top'] = inter.BSpline(NGC5533['n_tt'], NGC5533['n_ct'], NGC5533['n_kt'])

# Total Curve #######################
NGC5533['total'] = {
    'r' : np.asarray(NGC5533['raw_total']['xx']);
    'v' : np.asarray(NGC5533['raw_total']['yy']);
}
NGC5533['total']['t'], NGC5533['total']['c'], NGC5533['total']['k'] = inter.splrep(NGC5533['total']['r'], NGC5533['total']['v'])
NGC5533['total']['spline'] = inter.BSpline(NGC5533['total']['t'], NGC5533['total']['c'], NGC5533['total']['k'])

# Black Hole ########################
NGC5533['blackhole'] = {
    'r' : np.asarray(NGC5533['raw_blackhole']['xx']);
    'v' : np.asarray(NGC5533['raw_blackhole']['yy']);
}
NGC5533['blackhole']['t'], NGC5533['blackhole']['c'], NGC5533['blackhole']['k'] = inter.splrep(NGC5533['blackhole']['r'], NGC5533['blackhole']['v'])
NGC5533['blackhole']['spline'] = inter.BSpline(NGC5533['blackhole']['t'], NGC5533['blackhole']['c'], NGC5533['blackhole']['k'])

# Bulge #############################
NGC5533['bulge'] = {
    'r' : np.asarray(NGC5533['raw_bulge']['xx']);
    'v' : np.asarray(NGC5533['raw_bulge']['yy']);
}
NGC5533['bulge']['t'], NGC5533['bulge']['c'], NGC5533['bulge']['k'] = inter.splrep(NGC5533['bulge']['r'], NGC5533['bulge']['v'])
NGC5533['bulge']['spline'] = inter.BSpline(NGC5533['bulge']['t'], NGC5533['bulge']['c'], NGC5533['bulge']['k'])

# Disk ##############################
NGC5533['disk'] = {
    'r' : np.asarray(NGC5533['raw_disk']['xx']);
    'v' : np.asarray(NGC5533['raw_disk']['yy']);
}
NGC5533['disk']['t'], NGC5533['disk']['c'], NGC5533['disk']['k'] = inter.splrep(NGC5533['disk']['r'], NGC5533['disk']['v'])
NGC5533['disk']['spline'] = inter.BSpline(NGC5533['disk']['t'], NGC5533['disk']['c'], NGC5533['disk']['k'])

# Halo ##############################
NGC5533['halo'] = {
    'r' : np.asarray(NGC5533['raw_halo']['xx']);
    'v' : np.asarray(NGC5533['raw_halo']['yy']);
}
NGC5533['halo']['t'], NGC5533['halo']['c'], NGC5533['halo']['k'] = inter.splrep(NGC5533['halo']['r'], NGC5533['halo']['v'])
NGC5533['halo']['spline'] = inter.BSpline(NGC5533['halo']['t'], NGC5533['halo']['c'], NGC5533['halo']['k'])

# Gas ###############################
NGC5533['gas'] = {
    'r' : np.asarray(NGC5533['raw_gas']['xx']);
    'v' : np.asarray(NGC5533['raw_gas']['yy']);
}
NGC5533['gas']['t'], NGC5533['gas']['c'], NGC5533['gas']['k'] = inter.splrep(NGC5533['gas']['r'], NGC5533['gas']['v'])
NGC5533['gas']['spline'] = inter.BSpline(NGC5533['gas']['t'], NGC5533['gas']['c'], NGC5533['gas']['k'])

###############################
########### NGC0891 ###########
###############################

NGC0891 = {
    'measured_data' : dp.getXYdata_wXYerr('data/891/891_data');
}

#Organize measured data
NGC0891['m_radii']      = NGC0891['measured_data']['xx']
NGC0891['m_velocities'] = NGC0891['measured_data']['yy']
NGC0891['m_v_errors']   = NGC0891['measured_data']['ey']
