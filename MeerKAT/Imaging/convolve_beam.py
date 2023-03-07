"""

Author: Joppe Swart
Created: December 2022
Last modified: January 2023
Description: This script convolves all channels to the worst bmaj, bmin and bpa.
             The script needs to be executed in CASA.

"""



def find_beam(channels, max_beam = 10):
    """
    find_beam: This function finds the beams of all channels
    INPUTS:
        channels: Number of channels for which the beam needs to be convolved.
        max_beam: The maximum beam size in arcsec. 
    OUTPUTS:
        worstbmaj: Largest major beam.
        worstbmin: Largest minor beam.
        worstbpa: Largest Position angle.
        failedchannels: The channels that contain a to large or to small beam.
    """

    import matplotlib.pylab as plt
    import numpy as np
    from astropy.io import fits

    teller = 0
    worstbmaj = 0
    worstbmin = 0
    worstbpa = 0
    allbmaj = []
    allbmin = []
    failedchannels = []
    allbpa = []

    for i in range(channels):
        print(f'finding beam for channel {i}')
        Q_image = f"{i:04d}-Q-image.fits"
        with fits.open(Q_image) as hdu:
            bmaj = hdu[0].header['BMAJ']*3600 # convert deg to arcsec
            bmin = hdu[0].header['BMIN']*3600 # convert deg to arcsec
            bpa = hdu[0].header['BPA']

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
                if bpa > worstbpa:
                    worstbpa = bpa
            print(f"channel {i} has the following axis: bmaj = {bmaj} and bmin = {bmin}")

            allbmaj.append(bmaj)
            allbmin.append(bmin)
            allbpa.append(bpa)
        teller += 1

    plt.plot(range(0,teller),allbmaj,'b.', label='Bmaj')
    plt.plot(range(0,teller),allbmin,'r.', label='Bmin')
    plt.xlabel('Channel')
    plt.ylabel('Beam size (arcsec)')
    plt.title('beam sizes')
    plt.legend()
    plt.savefig('all_beams.png')
    plt.close()
    print(f'All major beams are {allbmaj}')

    print(f"The failled channels are {failedchannels}")
    np.save(f'failedchannels.npy', failedchannels)
    return worstbmaj, worstbmin, worstbpa, failedchannels


def convolve_beam(channels, enlrgfac = 1.05):
    """
    convolve_beam: This function convolves all beams to the same (worst) resolution.
    INPUTS:
        channels: Number of channels we want to convolve.
        enlrgfac: enlarge factor to prevent CASA errors.
    OUTPUTS:
        A smoothed casa file and an smoothed fits file for each channel and polarisation.
    """

    import numpy as np
    from astropy.io import fits

    nchan = channels
    bmaj, bmin, bpa, failedchannels = find_beam(channels, max_beam = 20)
    bmaj = 40
    bmin = 40
    bpa = 0
    print(f"All channels will be convolved to the same resolution: bmaj = {bmaj} and bmin = {bmin} and bpa = {bpa}.")

#    bmaj *= enlrgfac
#    bmin *= enlrgfac

    teller = 0
    for i in range(channels):
        if teller in failedchannels:
            teller +=1
            continue
        if teller == channels:
            break
        print(f"smoothing channel {teller}")
        Q_image = f"{teller:04d}-Q-image.fits"
        U_image = f"{teller:04d}-U-image.fits"
        I_image = f"{teller:04d}-I-image.fits"

        if not os.path.exists('res_40arcsec_test'):
            os.mkdir('res_40arcsec_test')
        imsmooth(imagename=Q_image,kernel='gauss'
                ,major='%.2f'%(bmaj)+'arcsec',minor='%.2f'%(bmin)+'arcsec'
                ,pa='%.2f'%(bpa)+'deg',targetres=True,overwrite=True,outfile=f"res_40arcsec_test/{Q_image.split('.fits')[0]}.smoothed")
        imsmooth(imagename=U_image,kernel='gauss'
                ,major='%.2f'%(bmaj)+'arcsec',minor='%.2f'%(bmin)+'arcsec'
                ,pa='%.2f'%(bpa)+'deg',targetres=True,overwrite=True,outfile=f"res_40arcsec_test/{U_image.split('.fits')[0]}.smoothed")
        imsmooth(imagename=I_image,kernel='gauss'
                ,major='%.2f'%(bmaj)+'arcsec',minor='%.2f'%(bmin)+'arcsec'
                ,pa='%.2f'%(bpa)+'deg',targetres=True,overwrite=True,outfile=f"res_40arcsec_test/{I_image.split('.fits')[0]}.smoothed")

        print(f"Convolving channel {teller} to fits")
        smoothQ = f"res_40arcsec_test/{Q_image.split('.fits')[0]}.smoothed"
        fitsQ = f"{smoothQ}.fits"
        exportfits(imagename=smoothQ,overwrite=True,fitsimage=fitsQ)
        smoothU = f"res_40arcsec_test/{U_image.split('.fits')[0]}.smoothed"
        fitsU = f"{smoothU}.fits"
        exportfits(imagename=smoothU,overwrite=True,fitsimage=fitsU)
        smoothI = f"res_40arcsec_test/{I_image.split('.fits')[0]}.smoothed"
        fitsI = f"{smoothI}.fits"
        exportfits(imagename=smoothI,overwrite=True,fitsimage=fitsI)
        teller +=1

convolve_beam(channels=126)






