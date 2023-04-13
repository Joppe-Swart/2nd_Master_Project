"""

Author: Joppe Swart
Created: March 2023
Last modified: March 2023
Description: This scroipt convolves a given imagem to a givne bmaj c.
             The script needs to be executed in CASA.

"""


def convolve(imagename, refimage):
    """
    convolve: this funciton convolves the image
    INPUTS:
        imagename: name of the image yhou want convolved
        refimage: image you want ot convolve to 
    OUTPUTS:
        convolved fits image
    """

    import numpy as np
    from astropy.io import fits

    ref_hdu = fits.open(refimage)[0]
    ref_header = ref_hdu.header
    bmaj = ref_header['bmaj']*3600 #convert to arcseconds
    bmin = ref_header['bmin']*3600
    bpa = ref_header['bpa']
    imsmooth(imagename=imagename,kernel='gauss'
                ,major='%.2f'%(bmaj)+'arcsec',minor='%.2f'%(bmin)+'arcsec'
                ,pa='%.2f'%(bpa)+'deg',targetres=True,overwrite=True,outfile=f"{imagename.split('.fits')[0]}.smoothed")

    exportfits(imagename=f"{imagename.split('.fits')[0]}.smoothed",overwrite=True,fitsimage=f"{imagename.split('.fits')[0]}.smoothed.fits")
    return

convolve('selfcal_bulletcluster_003-MFS-image.fits','polint_Ricean_corr.fits')


