import numpy as np
#import casacore.tables as pt
import matplotlib
import matplotlib.pyplot as plt
#from flagdata import flagdata # For some reason we get global variable 'flagdata' not defined
#from flagmanager import flagmanager # For some reason we get global variable 'flagmanager' not defined
import datetime
import os



def flag_outliers(sourcename, all_clipmax, blcuts=[0,200,1000,2000,4000,8000]):
    """
    Use CASA's flagdata function to flag the data per bl cut, where the visibilities
    have an amplitude above clipmax, defined in find_outliers(). 
    The flags are backed up every time 'flagdata' is called. 
    all_clipmax -- see outlier_per_baseline_crosscor()
    blcuts            -- see outlier_per_baseline_crosscor()
    """


    for i in range(len(blcuts)-1):
        uvrange = str(blcuts[i])+'~'+str(blcuts[i+1]) # default units meters
        print ("Flagging data in uvrange "+uvrange+". Flaggin all vis amp > %.4f"%all_clipmax[i])

        # First flag on RL and LR correlation only.
        flagdata(vis=sourcename, mode='clip', action='apply',display=''
            ,clipminmax=[0,all_clipmax[i]], uvrange=uvrange, clipoutside=True, field=''
            , spw='',antenna='',timerange='',correlation='XY,YX',scan=''
            ,intent='',array='',observation='',feed='',autocorr=False,flagbackup=True)
        # Then extend flags to the RR and LL correlations, thus polarization axis.
        flagdata(vis=sourcename, mode='extend', spw='', extendpols=True,
            action='apply', display='', field='',antenna='',timerange='',correlation='',
            scan='',intent='',array='', uvrange=uvrange, observation='',feed='',
            combinescans=False, growaround=False,flagneartime=False,flagnearfreq=False,
            flagbackup=True)


all_clipmax = np.load('all_clipmax.npy')
flag_outliers('BC_Averaged_III.ms', all_clipmax, blcuts=[0,200,1000,2000,4000,8000])
