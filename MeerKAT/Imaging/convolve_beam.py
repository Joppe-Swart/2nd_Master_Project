"""
Script to calculate primary beams for given frequencies
This script needs python 3.6 or higher.
"""

import matplotlib.pylab as plt
import numpy as np
from astropy.io import fits
from astropy import wcs
import os


# Find the poorest bmaj and bmin from all the images
def find_beam(max_beam = 10, channels = 25):
    """
    Find the poorest bmaj and bmin for each channel.
    This code needs to be run where the fits files are stored
    """
    teller = 0
    worstbmaj = 0
    worstbmin = 0
    allbmaj = []
    allbmin = []
    failedchannels = []

    for i in range(channels):
        print(f'finding beam for channel {i}')
        Q_image = f"BC_III_QU-{i:04d}-Q-image.fits"
        with fits.open(Q_image) as hdu:
            bmaj = hdu[0].header['BMAJ']*3600 # convert deg to arcsec
            bmin = hdu[0].header['BMIN']*3600 # convert deg to arcsec

            if bmaj < 0.8 or bmin < 0.8:
                # We don't get sub-arcsec resolution, so likely the image is empty
                failedchannels.append(teller)
            elif bmaj > max_beam or bmin > max_beam:
                failedchannels.append(teller)
            else:
                if bmaj > worstbmaj :
                    worstbmaj = bmaj
                if bmin > worstbmin:
                    worstbmin = bmin
            print(f"channel {i} has the following axis: bmaj = {bmaj} and bmin = {bmin}")

            allbmaj.append(bmaj)
            allbmin.append(bmin)
        teller += 1

    plt.plot(range(0,teller),allbmaj,label='Bmaj')
    plt.plot(range(0,teller),allbmin,label='Bmin')
    plt.xlabel('Channel')
    plt.ylabel('Beam size (arcsec)')
    plt.title('beam sizes')
    plt.legend()
    plt.savefig('all_beams.png')
    plt.close()

    print(f"The failled channels are {failedchannels}")

    return worstbmaj, worstbmin, failedchannels


# Convolve each model to the same primary beam
def convolve_beam(channels = 25, enlrgfac = 1.2):
    """
    Make all images the same resolution
    """
    bmaj, bmin, failedchannels = find_beam()
    print(f"All channels will be convolved to the same resolution: bmaj = {bmaj} and bmin = {bmin}.")

    bmaj *= enlrgfac
    ### CIRCULAR RESOLUTION
    bmin = bmaj 
    bpa = 0

    teller = 0
    for i in range(channels):
        if teller in failedchannels:
            teller +=1
            continue
        if teller == channels:
            break
        print(f"smoothing channel {teller}")
        Q_image = f"BC_III_QU-{teller:04d}-Q-image.fits"
        U_image = f"BC_III_QU-{teller:04d}-U-image.fits"
        I_image = f"BC_III_I-{teller:04d}-image.fits"

        imsmooth(imagename=Q_image,kernel='gauss'
                ,major='%.2f'%(bmaj)+'arcsec',minor='%.2f'%(bmin)+'arcsec'
                ,pa='%.2f'%(bpa)+'deg',targetres=True,overwrite=True,outfile=f"{Q_image.split('.fits')[0]}.smoothed")
        imsmooth(imagename=U_image,kernel='gauss'
                ,major='%.2f'%(bmaj)+'arcsec',minor='%.2f'%(bmin)+'arcsec'
                ,pa='%.2f'%(bpa)+'deg',targetres=True,overwrite=True,outfile=f"{U_image.split('.fits')[0]}.smoothed")
        imsmooth(imagename=I_image,kernel='gauss'
                ,major='%.2f'%(bmaj)+'arcsec',minor='%.2f'%(bmin)+'arcsec'
                ,pa='%.2f'%(bpa)+'deg',targetres=True,overwrite=True,outfile=f"{I_image.split('.fits')[0]}.smoothed")

        print(f"Convolving channel {teller} to fits")
        smoothQ = f"{Q_image.split('.fits')[0]}.smoothed"
        fitsQ = f"{smoothQ}.fits"
        exportfits(imagename=smoothQ,overwrite=True,fitsimage=fitsQ)
        smoothU = f"{U_image.split('.fits')[0]}.smoothed"
        fitsU = f"{smoothU}.fits"
        exportfits(imagename=smoothU,overwrite=True,fitsimage=fitsU)
        smoothI = f"{I_image.split('.fits')[0]}.smoothed"
        fitsI = f"{smoothI}.fits"
        exportfits(imagename=smoothI,overwrite=True,fitsimage=fitsI)
        teller +=1


convolve_beam(channels = 25, enlrgfac = 1.2)





