import numpy as np
import casacore.tables as pt
import matplotlib
import matplotlib.pyplot as plt
#from flagdata import flagdata # For some reason we get global variable 'flagdata' not defined
#from flagmanager import flagmanager # For some reason we get global variable 'flagmanager' not defined
import datetime
import os

matplotlib.rcParams['agg.path.chunksize'] = 10000
plt.ioff()

def find_outliers(uvdist, amp, thresh=5.0):
    """
    Find outliers in the Amp vs UVdist plane that are > thresh*sigma away from mean
    This function averages all correlations and channels, so we give it a sliced
    version of the amp array when finding outliers in the cross correlation data.
    uvdist -- np.array: data column of uv distances (in meters)
    amp    -- np.array: amplitude of visibilities, shape: (?,6,n_corr) (?,channel,correlation)
    thresh -- float: sigma threshold
    Returns
    outlier_mask -- np.array: mask of size (?,6,4), for every channel (6) and correlation (4)
                 OR np.array: mask of size (?,6,2), for every channel (6) and RL,LR correlation (2)
    clipmax      -- float: mean(amp) + thresh*std(amp), for defining the clip in flagdata.
    """
    # Average over all axes
    meanamp = np.ma.mean(amp) 
    stdamp = np.ma.std(amp)
    clipmax = meanamp + thresh*stdamp 
    outlier_mask = (amp > clipmax) 

    return outlier_mask, clipmax


def plot_uv_distance(sourcename, blcuts=[0,200]):

    # Load columns from .MS into python
    with pt.table(sourcename) as table1:
        # Pre-select some columns with SQL like language, sorted by spw
        table1 = table1.query(query="", columns="DATA,UVW,FLAG,ANTENNA1,DATA_DESC_ID"
            , sortlist='DATA_DESC_ID')
        data = table1.getcol("DATA")
        flags = table1.getcol('FLAG')
        data = np.ma.masked_array(data, flags)
        amp = np.abs(data)
        uvw = table1.getcol("UVW")
        antenna1 = table1.getcol("ANTENNA1")
        spw = table1.getcol("DATA_DESC_ID")
        # Projected baseline separations.
        uvdist = np.sqrt(np.square(uvw[:,0]) + np.square(uvw[:,1]))
    print(amp.shape[1])
    for i in range(len(blcuts)-1):
        if not os.path.exists('UV_Images'):
           os.mkdir('UV_Images')
#        mask = (uvdist >= blcuts[i]) & (uvdist <= blcuts[i+1])
        for ch in range(amp.shape[1]):
            print('plotting channel: ', ch)
            plt.figure()
            plt.title('visibilities channel '+str(ch))
            plt.xlabel('uvdist [m]')
            plt.ylabel('amp')
            plt.plot(uvdist, amp[:,ch,1], 'k.', markersize=0.4, alpha = 0.5)
            plt.plot(uvdist, amp[:,ch,2], 'r.', markersize=0.4, alpha = 0.5)
            plt.savefig('UV_Images/uvdist_'+str(ch)+'.png')
            plt.close()

def outliers_per_chan(sourcename, blcuts=[0,200,1000,2000,4000,8000], thresh=5.0):
    # Load columns from .MS into python
    with pt.table(sourcename) as table1:
        # Pre-select some columns with SQL like language, sorted by spw
        table1 = table1.query(query="", columns="DATA,UVW,FLAG,ANTENNA1,DATA_DESC_ID"
            , sortlist='DATA_DESC_ID')
        data = table1.getcol("DATA")
        flags = table1.getcol('FLAG')
        data = np.ma.masked_array(data, flags)
        amp = np.abs(data)
        uvw = table1.getcol("UVW")
        antenna1 = table1.getcol("ANTENNA1")
        spw = table1.getcol("DATA_DESC_ID")
        # Projected baseline separations.
        uvdist = np.sqrt(np.square(uvw[:,0]) + np.square(uvw[:,1]))

    all_outlier_masks = [] # For each blcut, we get 16 (spw) outlier masks
    all_clipmax = [] # And for each blcut we get 16 (spw) clipmaxes

    for i in range(len(blcuts)-1):
        for ch in range(amp.shape[1]):
            # amplitudes and uv dist for this spw
            amp_spw = amp[:,ch,1]
            # Amp and uvdist for this spw and this baseline cut
            mask = (uvdist >= blcuts[i]) & (uvdist <= blcuts[i+1])

            outlier_mask, clipmax = find_outliers(uvdist[mask], amp_spw[mask], thresh)
            all_outlier_masks.append(outlier_mask)
            all_clipmax.append(clipmax)
            np.save('all_clipmax.npy', all_clipmax)

    return all_outlier_masks, all_clipmax, blcuts



outliers_per_chan('BC_Averaged_III.ms', blcuts=[0,200,1000,2000,4000,8000], thresh=5.0)


#plot_uv_distance('BC_Averaged_III.ms')
