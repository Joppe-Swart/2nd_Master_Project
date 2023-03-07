"""

Author: Joppe Swart
Created: December 2022
Last modified: January 2023
Description: This script calculates the primary beam corrections for each channel.

"""

def showbeam(beam,freqMHz,pol,beamextent, i, size):
    """
    showbeam: This function makes the primary beam per channel using JimBeam
    INPUTS:
        beam: Standard beam package from JimBeam.
        sourcename: Name of the data block of interest.
        freqMHz: Frequency in MHz where we want the beam.
        pol: Polarization for which we want the beam.
        beamextent: Width of the image.
        i: Index for the file name.
    OUTPUTS:
        fits file containing the primary beam.
    """

    import numpy as np
    from astropy.io import fits

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
    Q_hdu = fits.open(f"{i:04d}-Q-image.fits")[0]
    Q_data = Q_hdu.data[0]
    hdu = fits.PrimaryHDU()
    hdu.header = Q_hdu.header
    hdu.data = np.reshape(beampixels, (-1,size, size))
    print(f'shape of the pb {hdu.data.shape}')
    hdu.writeto(f"{i:04d}-{pol}-pb_model.fits", overwrite=True)


def calculate_pb(size, scale, channels):
    """
    calculate_pb: Calculate the primary beam per channel
    INPUTS:
        sourcename: Name of the data block of interest.
        size: Pixel size of the image for which we want the pb. Assumes square image.
        scale: Scale of each pixel in arcseconds.
        channels: Number of channels for which we want the primary beam.
    """

    from katbeam import JimBeam
    import numpy as np
    from astropy.io import fits
    import os

    failedchannels = np.load(f'failedchannels.npy')
    print(f'Calculating primary beam corrections')
    teller = 0
    for i in range((channels)):
        if teller in failedchannels:
            print(f"channel {teller} failed.")
            teller +=1
            continue
        if teller == channels:
            break
        print(f'Calculating the primary beam corrections on channel {teller}')
        I_hdu = fits.open(f"{teller:04d}-I-image.smoothed.fits")[0]

        freq = I_hdu.header['CRVAL3']/1e6 #convert to MHz
        if not os.path.exists(f"{teller:04d}-I-pb_model.fits"):
            lbeam=JimBeam('MKAT-AA-L-JIM-2020')
            beam_extend = size*scale/3600. #convert to deg
            showbeam(lbeam,freq,'I',beam_extend, teller, size)

        teller +=1

    print(f'primary beam corrections calculated, apply the corrections with the apply script')

calculate_pb(size=8192, scale=1.1, channels=126)




