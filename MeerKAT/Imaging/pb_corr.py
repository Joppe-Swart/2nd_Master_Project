"""
Script to calculate primary beams for given frequencies
This script needs python 3.6 or higher.
"""

import matplotlib.pylab as plt
from katbeam import JimBeam
import numpy as np
from astropy.io import fits
from astropy import wcs
import os

size = 8192
scale = 1.1 #arcsecond
channels = 25

# Make the primary beam
def showbeam(beam,freqMHz=1000,pol='I',beamextent=10., i = 0):
      """
      Function to calculate the pb corrections
      """
      print(f"Calculating primary beam for Stokes {pol} at {freqMHz} MHz.")
      margin=np.linspace(-beamextent/2.,beamextent/2.,8192)
      x,y=np.meshgrid(margin,margin)
      if pol=='H':
          beampixels=beam.HH(x,y,freqMHz)
      elif pol=='V':
          beampixels=beam.VV(x,y,freqMHz)
      else:
          beampixels=beam.I(x,y,freqMHz)
          pol='I'
      Q_hdu = fits.open(f"BC_III_QU-{i:04d}-Q-image.fits")[0]
      Q_data = Q_hdu.data[0]
      hdu = fits.PrimaryHDU()
      hdu.header = Q_hdu.header
      hdu.data = np.reshape(beampixels, (-1,size, size))
      print(f'shape of the pb {hdu.data.shape}')
      hdu.writeto(f"BC_III_{freqMHz:.1f}MHz-{i:04d}-{pol}-pb_model.fits", overwrite=True)

# Find the poorest bmaj and bmin from all the images
def find_beam(max_beam = 10, channels = channels):
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

#    plt.plot(range(0,teller),allbmaj,label='Bmaj')
#    plt.plot(range(0,teller),allbmin,label='Bmin')
#    plt.xlabel('Channel')
#    plt.ylabel('Beam size (arcsec)')
#    plt.title('beam sizes')
#    plt.legend()
#    plt.savefig(figpath+'all_cube_beams.png')
#    plt.close()

    print(f"The failled channels are {failedchannels}")

    return worstbmaj, worstbmin, failedchannels

# Apply the Primary beam corrections
def pbcor(size = size, scale = 1.1, channels = channels):
    """
    Function that applies the primary beam to each channel
    """
    _,_,failedchannels = find_beam()
    print(f'Applying Primary beam corrections')
    teller = 0
    for i in range((channels)):
        if teller in failedchannels:
            print(f"channel {teller} failed.")
            teller +=1
            continue
        if teller == channels:
            break
        print(f'Calculating and applying primary beam corrections on channel {teller}')
        I_hdu = fits.open(f"BC_III_I-{teller:04d}-image.smoothed.fits")[0]
        Q_hdu = fits.open(f"BC_III_QU-{teller:04d}-Q-image.smoothed.fits")[0]
        U_hdu = fits.open(f"BC_III_QU-{teller:04d}-U-image.smoothed.fits")[0]

        I_header = I_hdu.header
        Q_header = Q_hdu.header
        U_header = U_hdu.header

        I_data = I_hdu.data[0]
        Q_data = Q_hdu.data[0]
        U_data = U_hdu.data[0]

        freq = Q_header['CRVAL3']/1e6 #convert to MHz
        if not os.path.exists(f"BC_III_{freq:.1f}MHz-{teller:04d}-I-pb_model.fits"):
            lbeam=JimBeam('MKAT-AA-L-JIM-2020')
            beam_extend = size*scale/3600. #convert to deg
            showbeam(lbeam,freq,'I',beam_extend, teller)
        pb_model_hdu = fits.open(f"BC_III_{freq:.1f}MHz-{teller:04d}-I-pb_model.fits")[0]
        pb_model_header = pb_model_hdu.header
        pb_model_data = pb_model_hdu.data
        print(f"shape of the pb that is used {pb_model_data.shape}")

        apply_beam_Q = Q_data/pb_model_data #data
        apply_beam_U = U_data/pb_model_data
        apply_beam_I = I_data/pb_model_data

        hdu_Q_pb = fits.PrimaryHDU()
        hdu_Q_pb.header = Q_header
        hdu_Q_pb.data = apply_beam_Q
        hdu_Q_pb.writeto(f"BC_III_QU-{teller:04d}-Q-image-pb.smoothed.fits", overwrite=True)

        hdu_U_pb = fits.PrimaryHDU()
        hdu_U_pb.header = U_header
        hdu_U_pb.data = apply_beam_U
        hdu_U_pb.writeto(f"BC_III_QU-{teller:04d}-U-image-pb.smoothed.fits", overwrite=True)

        hdu_I_pb = fits.PrimaryHDU()
        hdu_I_pb.header = I_header
        hdu_I_pb.data = apply_beam_I
        hdu_I_pb.writeto(f"BC_III_I-{teller:04d}-image-pb.smoothed.fits", overwrite=True)
        teller +=1

    print(f'Pb-corrections applied.')


pbcor(size = size, scale = scale, channels = channels)




