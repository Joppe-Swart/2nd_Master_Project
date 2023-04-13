"""

Author: Joppe Swart
Created: March 2023
Last modified: March 2023
Description: This script helpf by making nice figures of Stokes I images

"""

import numpy as np
import astropy
from astropy.io import fits

def findrms(im,maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    mIn = np.ndarray.flatten(im)
    m=mIn[np.abs(mIn)>maskSup]
    rmsold=np.nanstd(m)
#    print('std = ', rmsold)
    diff=1e-2
    cut=5.
    bins=np.arange(np.nanmin(m),np.nanmax(m),(np.nanmax(m)-np.nanmin(m))/30.)
    med=np.nanmedian(m)
    mean = np.nanmean(m)
#    print('the medians are', med, mean)
    for i in range(10):
        ind=np.where(np.abs(m-med)<rmsold*cut)[0]
        rms=np.nanstd(m[ind])
#        print(i, 'std = ', rms)
        if np.abs((rms-rmsold)/rmsold)<diff: break
        rmsold=rms
    print(f"The rms of the noise is given by {rms}")
    return rms

def calculate_rms(channels):
    """
    calculate_rms: This function calculates the rms of all channels
    INPUTS:
        sourcename: Data block of interest. 
        channels: Number of channels for which we need the rms.
    OUTPUTS:
        Three npy files with the noise in the stokes channels.
    """

    all_chan = np.arange(channels)
    all_noiseQ = [] 
    all_noiseU = []
    all_noiseI = [] 

    teller = 0
    pb_hdu = fits.open('0055-I-pb_model.fits')[0]
    pb_data = pb_hdu.data[0]
    pbmask = np.full(pb_data.shape, np.nan)
    pbmask[np.where(pb_data>0.9)] = 1
    for i in range(channels):
        if teller in [20, 21, 23, 36,37,38,39,40,73,74,75]:
            teller +=1
            continue
        if teller == channels:
            break
        print(f"finding the noise for channel {i}")
        Q_image = fits.open(f"stokes_q_corr/{teller:04d}-Q-image-pb.smoothed.fits")[0].data[0]
        U_image = fits.open(f"stokes_u_corr/{teller:04d}-U-image-pb.smoothed.fits")[0].data[0]
        I_image = fits.open(f"stokes_i_corr/{teller:04d}-I-image-pb.smoothed.fits")[0].data[0]
        print(Q_image[0].shape, U_image[0].shape, I_image[0].shape)

        rmsQ = findrms(Q_image[0]*pbmask)
        rmsU = findrms(U_image[0]*pbmask)
        rmsI = findrms(I_image[0]*pbmask)

        all_noiseQ.append(rmsQ)
        all_noiseU.append(rmsU)
        all_noiseI.append(rmsI)
        teller += 1

    all_noiseQ = np.array([all_noiseQ])
    all_noiseU = np.array([all_noiseU])
    all_noiseI = np.array([all_noiseI])

    np.save(f'all_noiseQ_corr.npy',all_noiseQ)
    np.save(f'all_noiseU_corr.npy',all_noiseU)
    np.save(f'all_noiseI_corr.npy',all_noiseI)

calculate_rms(126)

calculate_rms(126)
