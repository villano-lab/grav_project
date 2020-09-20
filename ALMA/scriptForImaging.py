import re

if re.search('^4.1', casadef.casa_version) == None:
 sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.1')


print "# Running clean."

listobs(vis='uid___A002_X72ae5e_X72b.ms.split.cal', listfile='uid___A002_X72ae5e_X72b.ms.split.cal.listobs')

os.system('rm -rf J1037.*')
clean(vis = 'uid___A002_X72ae5e_X72b.ms.split.cal',
  imagename = 'J1037',
  field = '0', 
  spw = '0,1,2,3',
  mode = 'mfs',
  interactive = T,
  imsize = [256,256],
  cell = '0.1arcsec',
  phasecenter = 0,
  weighting = 'briggs',
  robust = 2.0)
# Beam is 0.67"x0.48"
# RMS = 0.2 mJy/beam

os.system('rm -rf J1058.*')
clean(vis = 'uid___A002_X72ae5e_X72b.ms.split.cal',
  imagename = 'J1058',
  field = '3', 
  spw = '0,1,2,3',
  mode = 'mfs',
  interactive = T,
  imsize = [256,256],
  cell = '0.1arcsec',
  phasecenter = 3,
  weighting = 'briggs',
  robust = 2.0)
# Beam is 0.69"x0.50"
# RMS = 0.3 mJy/beam

os.system('rm -rf hz9.*')
clean(vis = 'uid___A002_X72ae5e_X72b.ms.split.cal',
  imagename = 'hz9',
  field = '4', 
  spw = '0,1,2,3',
  mode = 'mfs',
  interactive = T,
  imsize = [540,540],
  cell = '0.075arcsec',
  phasecenter = 4,
  weighting = 'briggs',
  niter = 10000,
  robust = 2.0)
# Beam is 0.64"x0.52"
# RMS = 0.04 mJy/beam
# 38 minutes on-source
# median PWV = 1.5
# Number of good antennas = 24 
# theoretical rms 0.04  mJy/beam
# requested rms is 0.04 mJy/beam

exportfits(imagename = 'hz9.image', fitsimage = 'hz9.fits')

