# ALMA Data Reduction Script
# $Id: scriptForPI.py,v 1.7 2014/02/11 11:53:29 dpetry Exp $

# Calibration application

import os
import sys
import glob

applyonly = True

savingslevel=0
if globals().has_key("SPACESAVING"):
    print 'SPACESAVING =', SPACESAVING
    if (type(SPACESAVING)!=int or SPACESAVING<0):
        sys.exit('ERROR: SPACESAVING value \"'+str(SPACESAVING)+'\" not permitted, must be int>0.\n'
                 +'Valid values: 0 = no saving,\n'
                 + '              1 = delete *.ms.split,\n'
                 + '              2 = delete *.ms and *.ms.split,\n'
                 + '            >=3 = delete *.ms, *.ms.split, and if possible *.ms.split.cal')
        
    savingslevel = SPACESAVING

if (os.path.basename(os.getcwd()) != 'script'):
    sys.exit('ERROR: Please start this script in directory \"script\".')

scriptnames = glob.glob('uid*.ms.scriptForCalibration.py')

if len(scriptnames) == 0:
    sys.exit('ERROR: No calibration script found.')

try:
    os.chdir('../raw')
except:
    sys.exit('ERROR: directory \"raw\" not present.\n'
             '       Please download your raw data and unpack it to create and fill directory \"raw\".') 

# check available disk space
tmppipe = os.popen("df -P -m $PWD | awk '/[0-9]%/{print $(NF-2)}'")
avspace = int((tmppipe.readline()).rstrip('\n'))
tmppipe.close()
tmppipe = os.popen("du -sm ../../* | cut -f1")
packspace = int((tmppipe.readline()).rstrip('\n'))
tmppipe.close()

fcalpresent = 0
if os.path.exists('../script/scriptForFluxCalibration.py'):
    fcalpresent = 1

spaceneed = packspace*(11.+fcalpresent*3.)   
    
if (savingslevel==1):
    print 'Will delete intermediate MSs named *.ms.split to save disk space.'
    spaceneed = packspace*(7.+fcalpresent*3.)   
elif (savingslevel==2):
    print 'Will delete intermediate MSs named *.ms and *.ms.split to save disk space.'
    spaceneed = packspace*(3.+fcalpresent*3.)   
elif (savingslevel>=3):
    print 'Will delete all intermediate MSs to save disk space.'
    spaceneed = packspace*(3.+fcalpresent*3.)   


print 'Found ',avspace,' MB of available free disk space.'
print 'Expect to need ',spaceneed,' MB of free disk space.'
if(spaceneed>avspace):
    sys.exit('ERROR: not enough free disk space. Need at least '+str(spaceneed)+' MB.')

asdmnames = glob.glob('uid*.asdm.sdm')


if len(asdmnames) == 0:
    sys.exit('ERROR: No ASDM found in directory \"raw\".')

print 'Found the following ASDMs:', asdmnames

for i in range(len(asdmnames)):
    asdmnames[i] = asdmnames[i].replace('.asdm.sdm', '')


for i in range(len(scriptnames)):
    scriptnames[i] = scriptnames[i].replace('.ms.scriptForCalibration.py', '')

if sorted(asdmnames) != sorted(scriptnames):
    sys.exit('ERROR: Inconsistency between ASDMs and calibration scripts.')


if os.path.exists('../calibrated'):
    os.chdir('../calibrated')
    sys.exit('WARNING: will stop here since directory '+os.path.abspath(os.path.curdir)
             +' already exists.\nPlease delete it first and then try again.')
    
print 'Creating destination directory for calibrated data.'
os.mkdir('../calibrated')

os.chdir('../calibrated')


for asdmname in asdmnames:

    print 'Processing ASDM '+asdmname

    os.mkdir(asdmname+'.calibration')

    os.chdir(asdmname+'.calibration')

    if not os.path.exists('../../raw/'+asdmname+'.asdm.sdm'):
        sys.exit('ERROR: cannot find raw/'+asdmname+'.asdm.sdm')

    os.system('ln -sf ../../raw/'+asdmname+'.asdm.sdm '+asdmname)

    execfile('../../script/'+asdmname+'.ms.scriptForCalibration.py')

    if not os.path.exists(asdmname+'.ms.split.cal'):
        print 'ERROR: '+asdmname+'.ms.split.cal was not created.'
    else:
        print asdmname+'.ms.split.cal was produced successfully, moving it to \"calibrated\" directory.'
        os.system('mv '+asdmname+'.ms.split.cal ..')
        if (savingslevel>=2):
            print 'Deleting intermediate MS ', asdmname+'.ms'
            os.system('rm -rf '+asdmname+'.ms')
        if (savingslevel>=1):
            print 'Deleting intermediate MS ', asdmname+'.ms.split'
            os.system('rm -rf '+asdmname+'.ms.split')

    os.chdir('..')


if len(asdmnames) > 1:
    os.system('ln -sf ../calibrated')
    if fcalpresent>0:
        execfile('../script/scriptForFluxCalibration.py')
        if (savingslevel>=3):
            print 'Deleting intermediate MS ', asdmname+'.ms.split.cal'
            os.system('rm -rf '+asdmname+'.ms.split.cal')

print 'Done. Please find results in current directory.'
