"""

Author: Joppe Swart
Created: January 2023
Last modified: January 2023
Description: This script checks channels and flaggs the bad ones

"""

import matplotlib.pylab as plt
import numpy as np
import astropy
from astropy.io import fits
from astropy import wcs
import os

import pyfits
import numpy as np

def RMS(image):
   data = pyfits.getdata(image,ignore_missing_end=True)
   rms = data.std()/data.mean()
   print(rms)
   return rms

def calculate_rms(sourcename, channels):
    """
    calculate_rms: This function calculates the rms of all channels
    INPUTS:
        sourcename: Data block of interest. 
        channels: Number of channels for which we need the rms.
    OUTPUTS:
        Three npy files with the noise in the stokes channels.
    """

    import numpy as np

    all_chan = np.arange(channels)
    all_noiseQ = np.empty(channels) 
    all_noiseU = np.empty(channels) 
    all_noiseI = np.empty(channels) 

    # Use hinges-fences with fence=1.0 so we cut the highest and lowest values
    # from the img
    box = '3496,3496,4696,4696' # Calculate central rms
    # algorithm = 'classic'
    algorithm = 'hinges-fences'

    teller = 0
    for i in range(channels):
        if teller in failedchannels:
            all_noiseQ[teller] = np.nan
            all_noiseU[teller] = np.nan
            all_noiseI[teller] = np.nan
            teller +=1
            continue
        if teller == channels:
            break
        rmsQ = imstat(f"stokes_q/{sourcename}-{teller:04d}-Q-image-pb.smoothed.fits",
                      algorithm=algorithm,fence=1.0,box=box)['rms'][0]
        rmsU = imstat(f"stokes_u/{sourcename}-{teller:04d}-U-image-pb.smoothed.fits",
                      algorithm=algorithm,fence=1.0,box=box)['rms'][0]
        rmsI = imstat(f"stokes_i/{sourcename}-{teller:04d}-I-image-pb.smoothed.fits",
                      algorithm=algorithm,fence=1.0,box=box)['rms'][0]

        if rmsQ > 2e-3 or rmsU > 2e-3: 
            # If RMS above 2 mJy then say it's a failed channel
            # Usually RMS is around 200 microJansky (in single channel)
            failedchannels = np.append(failedchannels,teller)
            teller += 1

        else:
            all_noiseQ[teller] = rmsQ
            all_noiseU[teller] = rmsU
            all_noiseI[teller] = rmsI
            teller += 1

    failedchannels = np.sort(failedchannels)
    np.save(f'{sourcename}_failedchannels.npy',failedchannels)
    print(f'All failed channnels are {failedchannels}')

    np.save(f'{sourcename}_all_noiseQ.npy',all_noiseQ)
    np.save(f'{sourcename}_all_noiseU.npy',all_noiseU)
    np.save(f'{sourcename}_all_noiseI.npy',all_noiseI)


