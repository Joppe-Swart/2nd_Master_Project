"""
This script performs RM synthesis on the QU cubes for data block II
Run in Casa
First we make a parameter file
Secondly we apply this file to the data
"""

import matplotlib.pylab as plt
import numpy as np
import astropy
from astropy.io import fits
from astropy import wcs
import os

def calculate_rms(sourcename, channels =25):
    """
    Function to calculate the rms per channel
    """
    all_chan = np.arange(channels)

    all_noiseQ = np.empty(channels) 
    all_noiseU = np.empty(channels) 
    all_noiseI = np.empty(channels) 

    # Ignore the failed channels
    failedchannels = np.load(f'{sourcename}_failedchannels.npy')
    
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


def create_inverse_variance_weights(sourcename):
    """
    Given the rms noise levels that were just calculated. Translate these
    to weights.
    """


    all_noiseQ = np.load(f'{sourcename}_all_noiseQ.npy')
    all_noiseU = np.load(f'{sourcename}_all_noiseU.npy')
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
    with open(f'{sourcename}_rmsynth.par','w') as file:
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
        file.write(r'phi_min -10000')
        file.write('\n')
        file.write(r'nphi 400')
        file.write('\n')
        file.write(r'dphi 50')
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
        file.write(f'outputfn {sourcename}_rmsynth')
        file.write('\n')
        file.write('\n')

        file.write(r'% directory where the input fits file can be found')
        file.write('\n')
        file.write(f'input_dir /net/rijn9/data2/swart/DATA/MeerKAT_DATA/Imaging/Block_II/')
        file.write('\n')


calculate_rms(sourcename = 'BC_II', channels =25)
create_inverse_variance_weights(sourcename = 'BC_II')
create_parameterfile(sourcename = 'BC_II')

#run = "python /net/bovenrijn/data1/digennaro/software/pyrmsynth/rmsynthesis.py -s"+"BC_II_rmsynth.par"
#os.system(run)





