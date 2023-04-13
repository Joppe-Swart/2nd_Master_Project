"""

Author: Joppe Swart
Created: March 2023
Last modified: March 2023
Description: This script helpf by making nice figures of Stokes I images

"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
import astropy.units as u
import scipy.constants as c
from matplotlib import style
import scipy.stats as stats
import emcee
import corner
import astropy.units as u
from matplotlib import cm
from matplotlib.colors import ListedColormap
palette = cm.get_cmap("tab20", 20)
import aplpy
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

plt.style.use(['science'])


# Make the funcitions that find the noise in a given image
def findrms(im,maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    mIn = np.ndarray.flatten(im)
    m=mIn[np.abs(mIn)>maskSup]
    rmsold=np.nanstd(m)
    print('std = ', rmsold)
    diff=1e-2
    cut=5.
    bins=np.arange(np.nanmin(m),np.nanmax(m),(np.nanmax(m)-np.nanmin(m))/30.)
    med=np.nanmedian(m)
    mean = np.nanmean(m)
    print('the medians are', med, mean)
    for i in range(10):
        ind=np.where(np.abs(m-med)<rmsold*cut)[0]
        rms=np.nanstd(m[ind])
        print(i, 'std = ', rms)
        if np.abs((rms-rmsold)/rmsold)<diff: break
        rmsold=rms
    return rms

def rms_of_the_noise(im):#, lower = 3000, upper = 5000):
    im = np.ndarray.flatten(im[np.abs(im)>1e-7])#[lower:upper,lower:upper])
    std = np.nanstd(im)
#    print('the lenght is ', len(im))
    print('this is std', std)
    med=np.nanmedian(im)
    mean = np.nanmean(im)
    print(med, mean)
    index = np.where(np.abs(im) < 5*std)[0]
    N = np.sum(~np.isnan(im))
    print(N)
    rms = np.sqrt(1/N*np.nansum(im[np.where(np.abs(im-mean) < std)]**2))
    return rms


# Open all the fits images

# Pb corrections to find a suitable mask
pb_hdu = fits.open('0055-I-pb_model.fits')[0]
pb_data = pb_hdu.data
print(pb_data.shape)

# Import the fits file for Stokes I image noise
pol_I_hdu = fits.open('selfcal_bulletcluster_003-MFS-image.fits')[0]
I_header = pol_I_hdu.header
I_data = pol_I_hdu.data
print(I_data.shape)
I_freqmean = np.nanmean(I_data, axis=0)
I_wcs = WCS(I_header)

# Import the fits file for Stokes I primary beam corrected
pol_I_hdu_pb = fits.open('bullet_cluster_pb_corr.smoothed.fits')[0]
I_header_pb = pol_I_hdu_pb.header
I_data_pb = pol_I_hdu_pb.data
print(I_data_pb.shape)

# Import the polarised intensity map
polint_hdu = fits.open('rmsynth_polint.fits')[0]
polint_header = polint_hdu.header
polint_data = polint_hdu.data
np.reshape(polint_data, I_data_pb.shape)
print(polint_data.shape)

#Make the mask where we want to calculate the noise
pbmask = np.full(pb_data.shape, np.nan)
pbmask[np.where(pb_data>0.9)] = 1

#Correct for the ricean bias and make a corrected fits file
Ricean_bias = np.nanstd(np.ndarray.flatten(polint_data*pbmask))
print(Ricean_bias)
polint_data_corr = polint_data
print(polint_data_corr.shape)
polint_data_corr[np.where(polint_data<4*Ricean_bias)] = np.sqrt(np.abs(polint_data_corr[np.where(polint_data<4*Ricean_bias)]**2 - Ricean_bias**2))
polint_data_corr[np.where(polint_data>4*Ricean_bias)] = np.sqrt(np.abs(polint_data_corr[np.where(polint_data>4*Ricean_bias)]**2 - 2.3*Ricean_bias**2))
#polint_data_corr = np.sqrt(np.abs(polint_data_corr[np.where(polint_data_corr>4*Ricean_bias)]**2 - 2.3*Ricean_bias**2)) #polint_data - 8.8e-5 #
hdu_PI_corr = fits.PrimaryHDU()
hdu_PI_corr.header = polint_header
hdu_PI_corr.data = polint_data_corr
hdu_PI_corr.writeto(f"polint_Ricean_corr.fits", overwrite=True)


#Calculate the noise 
imagenoise = findrms(I_data*pbmask)#rms_of_the_noise(I_data[0]*pbmask)# 
print('the noise of stokes I is given by', imagenoise)

imagenoise_pb = findrms(I_data_pb*pbmask)#rms_of_the_noise(I_data_pb[0]*pbmask)#findrms(I_data_pb[0],lower = 0, upper = 8192)#
print('the noise of stokes I pb is given by', imagenoise_pb)

imagenoise_polint = findrms(polint_data_corr*pbmask)#findrms(polint_data,lower = 0, upper = 8192)#rms_of_the_noise(polint_data*pbmask)#
print('the image noise of the polarised intensity is given by', imagenoise_polint)

imagename = 'Bullet_Cluster_Polarisation_degree'
# Calculate the polint map and make a fits
print(polint_data_corr.shape)
hdu = fits.PrimaryHDU()
hdu.header = I_header
hdu.data = 100*polint_data_corr/I_data_pb
print(f'shape of the pb {hdu.data.shape}')
hdu.writeto(f"polarization_fraction_map.fits", overwrite=True)

mask=np.full(I_data.shape, np.nan)
print(mask.shape,'mask shape')
print(hdu.data.shape, 'polarisation_map shape')
mask[np.where((I_data>5*imagenoise) & (hdu.data>1) & (polint_data_corr>5*imagenoise_polint))] =1

hdu_masked_polfrac = fits.PrimaryHDU()
hdu_masked_polfrac.header = I_header
hdu_masked_polfrac.data = hdu.data*mask
print

hdu_masked_polfrac.writeto(f"polarization_fraction_map_masked.fits", overwrite=True)

hdu_estimat = fits.PrimaryHDU()
hdu_estimat.header = I_header
hdu_estimat.data = 100*5*imagenoise_polint/I_data_pb
hdu_estimat.writeto(f"upper_limitmap.fits", overwrite=True)


lev_factor = 3.
levs = np.sqrt([1.,16.,64, 256,1024]) # 3, 12, 48, 96
levelsr = np.ndarray.tolist(lev_factor*imagenoise_pb*levs)

sigma_levels=[1,3,5,10,20]

#plt.figure(figsize=(6,6))
#plt.subplot(projection=I_wcs, slices=('x', 'y', 0, 0))
#map1 = plt.imshow(hdu.data[0,0,:,:]*mask[0,0,:,:], cmap='rainbow',vmin=0, vmax=40, interpolation='nearest')
#plt.colorbar(map1, label = r'\%p')
#plt.xlabel('RA')
#plt.ylabel('DEC')
##plt.xlim(4100, 4400)
##plt.ylim(4350, 4650)
#plt.xlim(3496, 4696)
#plt.ylim(3496, 4696)
#plt.title(r'Polarisation fraction')
#plt.savefig(r'Polarisation_fraction.pdf', dpi=400)
#plt.contour(I_data_pb[0,0,:,:], levels=levelsr, colors='black', alpha=0.7, linewidths=0.35)
#plt.savefig(r'Polarisation_fraction_contour.pdf', dpi=400)
#plt.show()

upper_mask=np.full(I_data.shape, np.nan)
upper_mask[np.where((I_data>3*imagenoise) & (polint_data_corr>0))] =1

plt.figure(figsize=(6,6))
plt.subplot(projection=I_wcs, slices=('x', 'y', 0, 0))
map1 = plt.imshow(hdu_estimat.data[0,0,:,:]*upper_mask[0,0,:,:], cmap='rainbow',vmin=0, vmax=20, interpolation='nearest')
plt.colorbar(map1, label = r'\%p')
plt.xlabel('RA')
plt.ylabel('DEC')
#plt.xlim(4100, 4400)
#plt.ylim(4350, 4650)
plt.xlim(3496, 4696)
plt.ylim(3496, 4696)
plt.title(r'polarisation needed for 5sigma detection')
plt.savefig(r'Upper_Polarisation_fraction.pdf', dpi=400)
plt.contour(I_data_pb[0,0,:,:], levels=levelsr, colors='black', alpha=0.8, linewidths=0.35)
plt.savefig(r'Upper_Polarisation_fraction_contour.pdf', dpi=400)
plt.show()

#plt.figure()
#plt.subplot(projection=I_wcs, slices=('x', 'y', 0, 0))
#map1 = plt.imshow(1e6*I_data[0], cmap='bone',vmin=-24, vmax=64, interpolation='nearest')
#plt.colorbar(map1, label = r'$\mu$Jy/beam')
#plt.xlabel('RA')
#plt.ylabel('DEC')
#plt.xlim(3496, 4696)
#plt.ylim(3496, 4696)
## apertures.plot(color='white')
## annulus.plot(color='white')
##plt.title(r'Stokes I')
##plt.savefig(r'Stokes_I_zoom.pdf', dpi=400)
#plt.close()

#plt.figure()
#plt.subplot(projection=I_wcs, slices=('x', 'y', 0, 0))
#map1 = plt.imshow(1e6*I_data[0], cmap='bone',vmin=-24, vmax=64, interpolation='nearest')
#plt.colorbar(map1, label = r'$\mu$Jy/beam')
#plt.xlabel('RA')
#plt.ylabel('DEC')
##plt.xlim(3496, 4696)
##plt.ylim(3496, 4696)

## apertures.plot(color='white')
## annulus.plot(color='white')
##plt.title(r'Stokes I')
##plt.savefig(r'Stokes_I.pdf', dpi=400)
#plt.close()



