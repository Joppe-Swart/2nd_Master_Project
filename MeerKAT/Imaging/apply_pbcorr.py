"""

Author: Joppe Swart
Created: January 2023
Last modified: January 2023
Description: This script applies the primary beam corrections calculated previously.
             This script needs to be executed in CASA

"""

def apply_pb(channels):
    """
    apply_pb:Function to apply the primary beam corrections. The primary beam files
          should be in the same directory.
    INPUTS:
        sourcename: Name of the data block to apply pb corrections to.
        channels: Number of channels that need to be corrected.
    OUTPUTS:
        Three folders with stokes .fits images and all the stokes as .im files.
    """

    import numpy as np
    from astropy.io import fits
    import os

    failedchannels = np.load(f'failedchannels.npy')
    print(f'Applying Primary beam corrections')
    teller = 124
    for i in range((channels)):
        if teller in failedchannels:
            print(f"channel {teller} is a failed channel")
            teller +=1
            continue
        if teller == channels:
            break
        print(f'Calculating and applying primary beam corrections on channel {teller}')

        I_hdu = fits.open(f"{teller:04d}-I-image.smoothed.fits")[0]
        freq = I_hdu.header['CRVAL3']/1e6 #convert Hz to MHz

        if not os.path.exists('stokes_i'):
            os.mkdir('stokes_i')
        if not os.path.exists('stokes_q'):
            os.mkdir('stokes_q')
        if not os.path.exists('stokes_u'):
            os.mkdir('stokes_u')

        print("Applying Q")
        impbcor(imagename=f"{teller:04d}-Q-image.smoothed.fits", \
                pbimage=f"{teller:04d}-I-pb_model.fits", \
                outfile=f"{teller:04d}-Q-image-pb.smoothed.im", mode='divide', overwrite=True)

        print("Applying U")
        impbcor(imagename=f"{teller:04d}-U-image.smoothed.fits", \
                pbimage=f"{teller:04d}-I-pb_model.fits", \
                outfile=f"{teller:04d}-U-image-pb.smoothed.im", mode='divide', overwrite=True)

        print("Applying I")
        impbcor(imagename=f"{teller:04d}-I-image.smoothed.fits", \
                pbimage=f"{teller:04d}-I-pb_model.fits", \
                outfile=f"{teller:04d}-I-image-pb.smoothed.im", mode='divide', overwrite=True)

        print("Exporting Q to fits in the stokes_q directory")
        exportfits(imagename=f"{teller:04d}-Q-image-pb.smoothed.im", \
                   fitsimage=f"stokes_q/{teller:04d}-Q-image-pb.smoothed.fits")

        print("Exporting U to fits in the stokes_u directory")
        exportfits(imagename=f"{teller:04d}-U-image-pb.smoothed.im", \
                   fitsimage=f"stokes_u/{teller:04d}-U-image-pb.smoothed.fits")

        print("Exporting I to fits in the stokes_i directory")
        exportfits(imagename=f"{teller:04d}-I-image-pb.smoothed.im", \
                   fitsimage=f"stokes_i/{teller:04d}-I-image-pb.smoothed.fits")

        teller +=1

    print(f'Primary beam corrections applied.')

apply_pb(channels=126)



