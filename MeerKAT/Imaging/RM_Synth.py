"""

Author: Joppe Swart
Created: January 2023
Last modified: January 2023
Description: This script makes the parameter and invariance files needed for RM synthesis.
             The script needs to be executed in CASA.

"""

import matplotlib.pylab as plt
import numpy as np
import astropy
from astropy.io import fits
from astropy import wcs
import os

def calculate_rms(channels):
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

    # Ignore the failed channels
    failedchannels = np.load(f'failedchannels.npy')
    
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
        rmsQ = imstat(f"stokes_q/{teller:04d}-Q-image-pb.smoothed.fits",
                      algorithm=algorithm,fence=1.0,box=box)['rms'][0]
        rmsU = imstat(f"stokes_u/{teller:04d}-U-image-pb.smoothed.fits",
                      algorithm=algorithm,fence=1.0,box=box)['rms'][0]
        rmsI = imstat(f"stokes_i/{teller:04d}-I-image-pb.smoothed.fits",
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
    np.save(f'failedchannels.npy',failedchannels)
    print(f'All failed channnels are {failedchannels}')

    np.save(f'all_noiseQ.npy',all_noiseQ)
    np.save(f'all_noiseU.npy',all_noiseU)
    np.save(f'all_noiseI.npy',all_noiseI)


def create_inverse_variance_weights(sourcename):
    """
    create_inverse_variance_weights: This function makes an inverse weight file for all channels
    INPUTS:
        sourcename: Data block of interest. 
    OUTPUTS:
        A text file containing the invariance weight for each channel.
    """

    import numpy as np

    all_noiseQ = np.load(f'all_noiseQ.npy')
    all_noiseU = np.load(f'all_noiseU.npy')
    # Remove failed channels if any
    nanmask = np.invert(np.isnan(all_noiseQ))
    
    all_noiseQ = all_noiseQ[nanmask]
    all_noiseU = all_noiseU[nanmask]

    averagerms = (all_noiseQ + all_noiseU)/2.

    weights = (1/averagerms)**2 # INVERSE VARIANCE WEIGHING

    # Normalize so that they sum to 1
    weights /= np.sum(weights)

    # save the weights to a file
    weightfile = open('inv_var_weights.txt','w')
    for i, weight in enumerate(weights):
        weightfile.write(str(weight))
        if i != len(weights)-1:
            weightfile.write('\n')
    weightfile.close()

def create_parameterfile(sourcename):
    """
    create_parameterfile: This function creates a parameter file for RM synthesis.
    INPUTS:
        sourcename: Data block of interest. 
    OUTPUTS:
        A text file containing the parameterds needed for RM synthesis.
    """

    with open(f'dphi5rmsynth.par','w') as file:
        # This is a comment
        file.write(r'% Parameter file for rmsynthesis python code')
        file.write('\n')
        file.write('\n')

        file.write(r'% ra and dec min and max of the subimage to process, given in pixels')
        file.write('\n')
        file.write(r'% a value of -1 means to use the bound of the image')
        file.write('\n')
        file.write(r'dec_min 3496')
        file.write('\n')
        file.write(r'dec_max 4696')
        file.write('\n')
        file.write(r'ra_min 3496')
        file.write('\n')
        file.write(r'ra_max 4696')
        file.write('\n')
        file.write('\n')

        file.write(r'% Define the phi axis, dphi in rad/m/m')
        file.write('\n')
        file.write(r'phi_min -1000')
        file.write('\n')
        file.write(r'nphi 400')
        file.write('\n')
        file.write(r'dphi 5')
        file.write('\n')
        file.write('\n')

        file.write(r'% Clean parameters. Gain is the loop gain, niter is the number of clean iterations')
        file.write('\n')
        file.write(r'do_clean False')
        file.write('\n')
        file.write(r'gain 0.05')
        file.write('\n')
        file.write(r'niter 50000')
        file.write('\n')
        file.write(r'cutoff 2e-5')
        file.write('\n')
        file.write('\n')

        file.write(r'% Weighting parameter. Give the name of the weight file (located in the input_dir). ')
        file.write('\n')
        file.write(r'do_weight inv_var_weights.txt')
        file.write('\n')
        file.write('\n')

        file.write(r'% Detection threshold on polarized intensity map')
        file.write('\n')
        file.write(r'threshold 5e-5')
        file.write('\n')
        file.write('\n')

        file.write(r'% output file')
        file.write('\n')
        file.write(f'outputfn dphi5/rmsynth')
        file.write('\n')
        file.write('\n')

        file.write(r'% directory where the input fits file can be found')
        file.write('\n')
        file.write(f'input_dir /net/rijn9/data2/swart/DATA/MeerKAT_DATA/Imaging_final/Images/')
        file.write('\n')
#calculate_rms(channels=126)
#create_inverse_variance_weights('BC')
create_parameterfile('BC')


