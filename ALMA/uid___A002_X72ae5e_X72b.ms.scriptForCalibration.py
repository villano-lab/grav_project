# ALMA Data Reduction Script

# Calibration

thesteps = []
step_title = {0: 'Import of the ASDM',
              1: 'Fix of SYSCAL table times',
              2: 'listobs',
              3: 'A priori flagging',
              4: 'Generation and time averaging of the WVR cal table',
              5: 'Generation of the Tsys cal table',
              6: 'Generation of the antenna position cal table',
              7: 'Application of the WVR, Tsys and antpos cal tables',
              8: 'Split out science SPWs and time average',
              9: 'Listobs, clear pointing table, and save original flags',
              10: 'Initial flagging',
              11: 'Putting a model for the flux calibrator(s)',
              12: 'Save flags before bandpass cal',
              13: 'Bandpass calibration',
              14: 'Save flags before gain cal',
              15: 'Gain calibration',
              16: 'Save flags before applycal',
              17: 'Restore flags before applycal',
              18: 'Application of the bandpass and gain cal tables'}

if 'applyonly' not in globals(): applyonly = False
try:
  print 'List of steps to be executed ...', mysteps
  thesteps = mysteps
except:
  print 'global variable mysteps not set.'
if (thesteps==[]):
  thesteps = range(0,len(step_title))
  print 'Executing all steps: ', thesteps

# The Python variable 'mysteps' will control which steps
# are executed when you start the script using
#   execfile('scriptForCalibration.py')
# e.g. setting
#   mysteps = [2,3,4]# before starting the script will make the script execute
# only steps 2, 3, and 4
# Setting mysteps = [] will make it execute all steps.

import re

import os

if applyonly != True: es = aU.stuffForScienceDataReduction() 


if re.search('^4.1', casadef.casa_version) == None:
 sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.1')


# CALIBRATE_AMPLI: Ganymede
# CALIBRATE_ATMOSPHERE: Ganymede,J1037-2934,J1058+0133,hz9
# CALIBRATE_BANDPASS: J1037-2934
# CALIBRATE_FLUX: Ganymede
# CALIBRATE_FOCUS: 
# CALIBRATE_PHASE: J1058+0133
# CALIBRATE_POINTING: J0750+1231,J1037-2934,J1058+0133
# OBSERVE_TARGET: hz9

# Using reference antenna = DV04

# Import of the ASDM
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  if os.path.exists('uid___A002_X72ae5e_X72b.ms') == False:
    importasdm('uid___A002_X72ae5e_X72b', asis='Antenna Station CalWVR Receiver Source')

# Fix of SYSCAL table times
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  from recipes.almahelpers import fixsyscaltimes
  fixsyscaltimes(vis = 'uid___A002_X72ae5e_X72b.ms')

print "# A priori calibration"

# listobs
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.listobs')
  listobs(vis = 'uid___A002_X72ae5e_X72b.ms',
    listfile = 'uid___A002_X72ae5e_X72b.ms.listobs')
  
  

# A priori flagging
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  flagdata(vis = 'uid___A002_X72ae5e_X72b.ms',
    mode = 'manual',
    spw = '1~24',
    autocorr = T,
    flagbackup = F)
  
  flagdata(vis = 'uid___A002_X72ae5e_X72b.ms',
    mode = 'manual',
    intent = '*POINTING*,*SIDEBAND_RATIO*,*ATMOSPHERE*',
    flagbackup = F)
  
  flagcmd(vis = 'uid___A002_X72ae5e_X72b.ms',
    inpmode = 'table',
    action = 'plot')
  
  flagcmd(vis = 'uid___A002_X72ae5e_X72b.ms',
    inpmode = 'table',
    action = 'apply')
  

# Generation and time averaging of the WVR cal table
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.wvr') 
  
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_X72ae5e_X72b.ms.wvrgcal')
  
  wvrgcal(vis = 'uid___A002_X72ae5e_X72b.ms',
    caltable = 'uid___A002_X72ae5e_X72b.ms.wvr',
    toffset = 0,
    tie = ['hz9,J1058+0133'],
    statsource = 'hz9')
  
  casalog.setlogfile(mylogfile)
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.wvr.smooth') 
  
  smoothcal(vis = 'uid___A002_X72ae5e_X72b.ms',
    tablein = 'uid___A002_X72ae5e_X72b.ms.wvr',
    caltable = 'uid___A002_X72ae5e_X72b.ms.wvr.smooth',
    smoothtype = 'mean',
    smoothtime = 2.016)
  
  
  if applyonly != True: aU.plotWVRSolutions(caltable='uid___A002_X72ae5e_X72b.ms.wvr.smooth', spw='17', antenna='DA44',
    yrange=[-180,180],subplot=22, interactive=False,
    figfile='uid___A002_X72ae5e_X72b.ms.wvr.smooth.plots/uid___A002_X72ae5e_X72b.ms.wvr.smooth') 
  
  #Note: If you see wraps in these plots, try changing yrange or unwrap=True 
  #Note: If all plots look strange, it may be a bad WVR on the reference antenna.
  #      To check, you can set antenna='' to show all baselines.
  
  
  flagdata(vis = 'uid___A002_X72ae5e_X72b.ms',
    mode = 'manual',
    spw = '0',
    autocorr = T,
    flagbackup = F)
  
  

# Generation of the Tsys cal table
mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.tsys') 
  gencal(vis = 'uid___A002_X72ae5e_X72b.ms',
    caltable = 'uid___A002_X72ae5e_X72b.ms.tsys',
    caltype = 'tsys')
  
  if applyonly != True: aU.plotbandpass(caltable='uid___A002_X72ae5e_X72b.ms.tsys', overlay='time', 
    xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,
    showatm=True,pwv='auto',chanrange='5~123',showfdm=True, 
    field='', figfile='uid___A002_X72ae5e_X72b.ms.tsys.plots.overlayTime/uid___A002_X72ae5e_X72b.ms.tsys') 
  
  
  if applyonly != True: es.checkCalTable('uid___A002_X72ae5e_X72b.ms.tsys', msName='uid___A002_X72ae5e_X72b.ms', interactive=False) 
  

# Generation of the antenna position cal table
mystep = 6
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Position for antenna DV19 is derived from baseline run made on 2013-08-17 08:36:51.
  
  # Position for antenna DV18 is derived from baseline run made on 2013-08-17 08:36:51.
  
  # Note: no baseline run found for antenna DA64.
  
  # Position for antenna DA48 is derived from baseline run made on 2013-03-08 03:29:55.
  
  # Position for antenna DA45 is derived from baseline run made on 2013-03-08 03:29:55.
  
  # Position for antenna DA44 is derived from baseline run made on 2013-03-08 03:29:55.
  
  # Position for antenna DA46 is derived from baseline run made on 2013-08-17 08:36:51.
  
  # Position for antenna DV15 is derived from baseline run made on 2013-08-17 08:36:51.
  
  # Position for antenna DA43 is derived from baseline run made on 2013-03-08 03:29:55.
  
  # Position for antenna DA42 is derived from baseline run made on 2013-03-08 03:29:55.
  
  # Note: no baseline run found for antenna DA62.
  
  # Position for antenna PM01 is derived from baseline run made on 2013-05-10 03:34:01.
  
  # Position for antenna PM04 is derived from baseline run made on 2013-05-10 03:34:01.
  
  # Position for antenna DV11 is derived from baseline run made on 2013-08-17 08:36:51.
  
  # Position for antenna DV20 is derived from baseline run made on 2013-08-17 08:36:51.
  
  # Position for antenna DA51 is derived from baseline run made on 2013-03-05 08:03:59.
  
  # Position for antenna DV25 is derived from baseline run made on 2013-08-17 08:36:51.
  
  # Position for antenna DV08 is derived from baseline run made on 2013-03-08 03:29:55.
  
  # Note: no baseline run found for antenna DV09.
  
  # Position for antenna DV07 is derived from baseline run made on 2013-08-17 08:36:51.
  
  # Position for antenna DV03 is derived from baseline run made on 2013-08-17 08:36:51.
  
  # Note: no baseline run found for antenna DV01.
  
  # Position for antenna DV17 is derived from baseline run made on 2013-08-17 08:36:51.
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.antpos') 
  gencal(vis = 'uid___A002_X72ae5e_X72b.ms',
    caltable = 'uid___A002_X72ae5e_X72b.ms.antpos',
    caltype = 'antpos',
    antenna = 'DV19,DV18,DA51,DA48,DV08,DA45,DA44,DA46,DV15,DV03,DA43,DA42,DV11,PM01,DV20,DV17,DV25,DV07,PM04',
    parameter = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
  #  parameter = [-0.000453608576208,0.00103624816984,-8.55620019138e-05,-0.000290210596096,0.000500840221623,0.000122694263593,0.000183368908713,-0.000450953651932,-0.000106764446923,-0.000279455019968,-0.000114494208283,-0.000351774029861,0.000681170379647,0.00049271482177,-0.000117141018864,-0.0008147232099,0.000443949533267,-4.03596724729e-05,0.000138477416171,-0.000600100629702,-0.000470491487109,-3.4172211787e-05,0.00029221096951,-6.70448300752e-05,-0.000176969249007,0.000250389621818,4.57407823213e-06,-0.000205369131134,0.000259081692311,0.000129727393741,-0.000183955207432,-0.000445721339222,6.14651886934e-05,-0.000324754071755,0.000268328647425,0.000213224130552,0.000203554357668,-0.000146355884478,-8.53907404398e-05,-0.000384319492324,0.000268677305476,9.47788654794e-05,0.000170314656522,-4.65931306241e-05,-8.45103422819e-05,0.000847928335973,-0.00116970769294,-0.000657849868993,6.39669597149e-05,0.000639797188342,-0.000138476956636,0.000463083386421,-0.000467116013169,-0.000287898816168,-0.000229949713758,0.000208240039609,4.43194657724e-05])
  

# Application of the WVR, Tsys and antpos cal tables
mystep = 7
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  
  from recipes.almahelpers import tsysspwmap
  tsysmap = tsysspwmap(vis = 'uid___A002_X72ae5e_X72b.ms', tsystable = 'uid___A002_X72ae5e_X72b.ms.tsys', tsysChanTol=1)
  
  
  
  applycal(vis = 'uid___A002_X72ae5e_X72b.ms',
    field = '0',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X72ae5e_X72b.ms.tsys', 'uid___A002_X72ae5e_X72b.ms.wvr.smooth', 'uid___A002_X72ae5e_X72b.ms.antpos'],
    gainfield = ['0', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  # Note: J0750+1231 didn't have any Tsys measurement, and I couldn't find any close measurement. But this is not a science target, so this is probably Ok.
  
  applycal(vis = 'uid___A002_X72ae5e_X72b.ms',
    field = '2',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X72ae5e_X72b.ms.tsys', 'uid___A002_X72ae5e_X72b.ms.wvr.smooth', 'uid___A002_X72ae5e_X72b.ms.antpos'],
    gainfield = ['2', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  applycal(vis = 'uid___A002_X72ae5e_X72b.ms',
    field = '3',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X72ae5e_X72b.ms.tsys', 'uid___A002_X72ae5e_X72b.ms.wvr.smooth', 'uid___A002_X72ae5e_X72b.ms.antpos'],
    gainfield = ['3', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  applycal(vis = 'uid___A002_X72ae5e_X72b.ms',
    field = '4',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X72ae5e_X72b.ms.tsys', 'uid___A002_X72ae5e_X72b.ms.wvr.smooth', 'uid___A002_X72ae5e_X72b.ms.antpos'],
    gainfield = ['4', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  if applyonly != True: es.getCalWeightStats('uid___A002_X72ae5e_X72b.ms') 
  

# Split out science SPWs and time average
mystep = 8
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split') 
  split(vis = 'uid___A002_X72ae5e_X72b.ms',
    outputvis = 'uid___A002_X72ae5e_X72b.ms.split',
    datacolumn = 'corrected',
    spw = '17,19,21,23',
    keepflags = T)
  
  

print "# Calibration"

# Listobs, clear pointing table, and save original flags
mystep = 9
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.listobs')
  listobs(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    listfile = 'uid___A002_X72ae5e_X72b.ms.split.listobs')
  
  tb.open('uid___A002_X72ae5e_X72b.ms.split/POINTING', nomodify = False)
  a = tb.rownumbers()
  tb.removerows(a)
  tb.close()
  
  if not os.path.exists('uid___A002_X72ae5e_X72b.ms.split.flagversions/Original.flags'):
    flagmanager(vis = 'uid___A002_X72ae5e_X72b.ms.split',
      mode = 'save',
      versionname = 'Original')
  
  

# Initial flagging
mystep = 10
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Flagging shadowed data
  flagdata(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    mode = 'shadow',
    flagbackup = F)
  
  # Flagging edge channels
  flagdata(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    mode = 'manual',
    spw = '0:0~7;120~127,1:0~7;120~127,2:0~7;120~127,3:0~7;120~127',
    flagbackup = F)

  # slightly low amplitudes at beginning of some scans
  flagdata(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    mode = 'quack',
    quackinterval = 4.0,
    quackmode = 'beg',
    flagbackup = F)

  # Flagging outlier antennas (optional)
  flagdata(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    mode = 'manual',
    antenna = 'DV07,DV19,DV25',
    flagbackup = F)
  
  

# Putting a model for the flux calibrator(s)
mystep = 11
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  delmod('uid___A002_X72ae5e_X72b.ms.split',field='2')

  setjy(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    field = '2', # Ganymede
    spw = '0,1,2,3',
    standard = 'Butler-JPL-Horizons 2012')
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.setjy.field*.png') 
  for i in ['2']:
    plotms(vis = 'uid___A002_X72ae5e_X72b.ms.split',
      xaxis = 'uvdist',
      yaxis = 'amp',
      ydatacolumn = 'model',
      field = i,
      spw = '0,1,2,3',
      avgchannel = '9999',
      coloraxis = 'spw',
      plotfile = 'uid___A002_X72ae5e_X72b.ms.split.setjy.field'+i+'.png')
  

# Save flags before bandpass cal
mystep = 12
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    mode = 'save',
    versionname = 'BeforeBandpassCalibration')
  
  

# Bandpass calibration
mystep = 13
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.ap_pre_bandpass') 
  
  gaincal(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    caltable = 'uid___A002_X72ae5e_X72b.ms.split.ap_pre_bandpass',
    field = '0', # J1037-2934
    spw = '0:51~76,1:51~76,2:51~76,3:51~76',
    solint = 'int',
    refant = 'DV04',
    calmode = 'p')
  
  if applyonly != True: es.checkCalTable('uid___A002_X72ae5e_X72b.ms.split.ap_pre_bandpass', msName='uid___A002_X72ae5e_X72b.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.bandpass') 
  bandpass(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    caltable = 'uid___A002_X72ae5e_X72b.ms.split.bandpass',
    field = '0', # J1037-2934
    solint = 'inf',
    combine = 'scan',
    refant = 'DV04',
    solnorm = T,
    bandtype = 'B',
    gaintable = 'uid___A002_X72ae5e_X72b.ms.split.ap_pre_bandpass')
  
  if applyonly != True: es.checkCalTable('uid___A002_X72ae5e_X72b.ms.split.bandpass', msName='uid___A002_X72ae5e_X72b.ms.split', interactive=False) 
  

# Save flags before gain cal
mystep = 14
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    mode = 'save',
    versionname = 'BeforeGainCalibration')
  
  

# Gain calibration
mystep = 15
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Note: the Solar system object used for flux calibration is highly resolved on some baselines.
  # Note: we will first determine the flux of the phase calibrator(s) on a subset of antennas.
  
  delmod('uid___A002_X72ae5e_X72b.ms.split',field='0')
  delmod('uid___A002_X72ae5e_X72b.ms.split',field='3')
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.phase_short_int') 
  gaincal(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    caltable = 'uid___A002_X72ae5e_X72b.ms.split.phase_short_int',
    field = '2', # Ganymede
    selectdata = T,
    antenna = 'DA44,DA62,DV08,DV11,DV20&',
    solint = 'int',
    refant = 'DV04',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_X72ae5e_X72b.ms.split.bandpass')
  
  gaincal(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    caltable = 'uid___A002_X72ae5e_X72b.ms.split.phase_short_int',
    field = '0,3', # J1037-2934,J1058+0133
    selectdata = T,
    solint = 'int',
    refant = 'DV04',
    gaintype = 'G',
    calmode = 'p',
    append = T,
    gaintable = 'uid___A002_X72ae5e_X72b.ms.split.bandpass')
  
  if applyonly != True: es.checkCalTable('uid___A002_X72ae5e_X72b.ms.split.phase_short_int', msName='uid___A002_X72ae5e_X72b.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.ampli_short_inf') 
  gaincal(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    caltable = 'uid___A002_X72ae5e_X72b.ms.split.ampli_short_inf',
    field = '0,2,3', # J1037-2934,Ganymede,J1058+0133
    selectdata = T,
    solint = 'inf',
    refant = 'DV04',
    gaintype = 'T',
    calmode = 'a',
    gaintable = ['uid___A002_X72ae5e_X72b.ms.split.bandpass', 'uid___A002_X72ae5e_X72b.ms.split.phase_short_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_X72ae5e_X72b.ms.split.ampli_short_inf', msName='uid___A002_X72ae5e_X72b.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.flux_short_inf') 
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.fluxscale') 
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_X72ae5e_X72b.ms.split.fluxscale')
  
  fluxscale(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    caltable = 'uid___A002_X72ae5e_X72b.ms.split.ampli_short_inf',
    fluxtable = 'uid___A002_X72ae5e_X72b.ms.split.flux_short_inf',
    reference = '2') # Ganymede
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: es.fluxscale2(caltable = 'uid___A002_X72ae5e_X72b.ms.split.ampli_short_inf', removeOutliers=True, msName='uid___A002_X72ae5e_X72b.ms', writeToFile=True, preavg=10000)
  
  f = open('uid___A002_X72ae5e_X72b.ms.split.fluxscale')
  fc = f.readlines()
  f.close()
  
  for phaseCalName in ['J1058+0133','J1037-2934']:
    for i in range(len(fc)):
      if fc[i].find('Flux density for '+phaseCalName) != -1 and re.search('in SpW=[0-9]+(?: \(ref SpW=[0-9]+\))? is: [0-9]+\.[0-9]+', fc[i]) != None:
        line = (re.search('in SpW=[0-9]+(?: \(ref SpW=[0-9]+\))? is: [0-9]+\.[0-9]+', fc[i])).group(0)
        spwId = (line.split('='))[1].split()[0]
        flux = float((line.split(':'))[1].split()[0])
        setjy(vis = 'uid___A002_X72ae5e_X72b.ms.split',
          field = phaseCalName.replace(';','*;').split(';')[0],
          spw = spwId,
          fluxdensity = [flux,0,0,0])
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.phase_int') 
  gaincal(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    caltable = 'uid___A002_X72ae5e_X72b.ms.split.phase_int',
    field = '0,2,3', # J1037-2934,Ganymede,J1058+0133
    solint = 'int',
    refant = 'DV04',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_X72ae5e_X72b.ms.split.bandpass')
  
  if applyonly != True: es.checkCalTable('uid___A002_X72ae5e_X72b.ms.split.phase_int', msName='uid___A002_X72ae5e_X72b.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.flux_inf') 
  gaincal(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    caltable = 'uid___A002_X72ae5e_X72b.ms.split.flux_inf',
    field = '0,2,3', # J1037-2934,Ganymede,J1058+0133
    solint = 'inf',
    refant = 'DV04',
    gaintype = 'T',
    calmode = 'a',
    gaintable = ['uid___A002_X72ae5e_X72b.ms.split.bandpass', 'uid___A002_X72ae5e_X72b.ms.split.phase_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_X72ae5e_X72b.ms.split.flux_inf', msName='uid___A002_X72ae5e_X72b.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.phase_inf') 
  gaincal(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    caltable = 'uid___A002_X72ae5e_X72b.ms.split.phase_inf',
    field = '0,2,3', # J1037-2934,Ganymede,J1058+0133
    solint = 'inf',
    refant = 'DV04',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_X72ae5e_X72b.ms.split.bandpass')
  
  if applyonly != True: es.checkCalTable('uid___A002_X72ae5e_X72b.ms.split.phase_inf', msName='uid___A002_X72ae5e_X72b.ms.split', interactive=False) 
  

# Save flags before applycal
mystep = 16
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    mode = 'save',
    versionname = 'BeforeApplycal')
  
  

# Restore flags before applycal
mystep = 17
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  if applyonly == True:
    flagmanager(vis = 'uid___A002_X72ae5e_X72b.ms.split',
      mode = 'restore',
      versionname = 'BeforeApplycal')
  
  

# Application of the bandpass and gain cal tables
mystep = 18
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i in ['0', '2']: # J1037-2934,Ganymede
    applycal(vis = 'uid___A002_X72ae5e_X72b.ms.split',
      field = i,
      gaintable = ['uid___A002_X72ae5e_X72b.ms.split.bandpass', 'uid___A002_X72ae5e_X72b.ms.split.phase_int', 'uid___A002_X72ae5e_X72b.ms.split.flux_inf'],
      gainfield = ['', i, i],
      interp = 'linear,linear',
      calwt = F,
      flagbackup = F)
  
  applycal(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    field = '3,4', # hz9
    gaintable = ['uid___A002_X72ae5e_X72b.ms.split.bandpass', 'uid___A002_X72ae5e_X72b.ms.split.phase_inf', 'uid___A002_X72ae5e_X72b.ms.split.flux_inf'],
    gainfield = ['', '3', '3'], # J1058+0133
    interp = 'linear,linear',
    calwt = F,
    flagbackup = F)
  
  os.system('rm -rf uid___A002_X72ae5e_X72b.ms.split.cal')
  split(vis = 'uid___A002_X72ae5e_X72b.ms.split',
    outputvis = 'uid___A002_X72ae5e_X72b.ms.split.cal',
    datacolumn = 'corrected')
