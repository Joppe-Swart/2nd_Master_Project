"""

Author: Joppe Swart
Created: April 2023
Last modified: April 2023
Description: Make the radio image with the contours of the x-ray image

"""

import numpy as np
import math
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
import scipy
import astropy.units as u
from matplotlib import cm
from matplotlib import colors
from matplotlib.colors import ListedColormap
palette = cm.get_cmap("tab20", 20)
import aplpy
from astropy.cosmology import FlatLambdaCDM
from mpl_toolkits.axes_grid1.anchored_artists import (AnchoredSizeBar)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

plt.style.use(['science'])

def flatten(f):
	"""
	Flatten a fits file so that it becomes a 2D image.
	Return new header and data.
	Taken from Jort Boxelaar
	"""

	naxis=f[0].header['NAXIS']
	if naxis<2:
		raise RadioError('Can\'t make map from this')
	if naxis == 2:
		return fits.PrimaryHDU(header=f[0].header,data=f[0].data)

	w  = wcs.WCS(f[0].header)
	wn = wcs.WCS(naxis=2)

	wn.wcs.crpix[0]=w.wcs.crpix[0]
	wn.wcs.crpix[1]=w.wcs.crpix[1]
	wn.wcs.cdelt=w.wcs.cdelt[0:2]
	wn.wcs.crval=w.wcs.crval[0:2]
	wn.wcs.ctype[0]=w.wcs.ctype[0]
	wn.wcs.ctype[1]=w.wcs.ctype[1]

	header = wn.to_header()
	header["NAXIS"]=2
	copy=('EQUINOX','EPOCH','BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
	for k in copy:
		r=f[0].header.get(k)
		if r is not None:
			header[k]=r

	slice=[]
	for i in range(naxis,0,-1):
		if i<=2:
			slice.append(np.s_[:],)
		else:
			slice.append(0)

	hdu = fits.PrimaryHDU(header=header,data=f[0].data[tuple(slice)])
	return hdu

# Make the funcitions that find the noise in a given image
def findrms(im,maskSup=1e-7):
	"""
	find the rms of an array, from Cycil Tasse/kMS
	"""
	mIn = np.ndarray.flatten(im)
	m=mIn[np.abs(mIn)>maskSup]
	rmsold=np.nanstd(m)
	#print('std = ', rmsold)
	diff=1e-2
	cut=3.
	bins=np.arange(np.nanmin(m),np.nanmax(m),(np.nanmax(m)-np.nanmin(m))/30.)
	med=np.nanmedian(m)
	mean = np.nanmean(m)
	#print('the medians are', med, mean)
	for i in range(15):
		ind=np.where(np.abs(m-med)<rmsold*cut)[0]
		rms=np.nanstd(m[ind])
		print(i, 'std = ', rms)
		if np.abs((rms-rmsold)/rmsold)<diff: break
		rmsold=rms
	return rms

save_loc = 'Results/Polarisation_maps'

pol_I_hdu = fits.open('DATA/selfcal_bulletcluster_003-MFS-image.fits')
pol_I_hdu = flatten(pol_I_hdu)
I_header = pol_I_hdu.header
I_data = pol_I_hdu.data
I_wcs = WCS(I_header)

xray_hdu = fits.open('DATA/bullet_0.5-2.0_flux.img')[0]
xray_header = xray_hdu.header
xray_data = xray_hdu.data
xray_wcs = WCS(xray_header)

#print(f"the shape of the xray data is {xray_data.shape}")
#imagenoise_xray = findrms(xray_data)
#print(f"the noise of the xray image is {imagenoise_xray}")

smooth_xray=scipy.ndimage.gaussian_filter(xray_data, 10/(2*math.sqrt(2*math.log(2))))
print(np.nanstd(smooth_xray))
imagenoise_xray = np.nanstd(smooth_xray)

smooth_xray1=scipy.ndimage.gaussian_filter(xray_data, 5/(2*math.sqrt(2*math.log(2))))
print(np.nanstd(smooth_xray1))
imagenoise_xray1 = np.nanstd(smooth_xray1)

lev_factor = 1
levs = np.sqrt([1, 4, 16, 64, 256, 1024, 4048]) #1, 5, 10, 20
#levs1 = np.sqrt([64, 256, 1024,4048])
levelsXray = np.ndarray.tolist(lev_factor*imagenoise_xray*levs)
#levelsXray1 = np.ndarray.tolist(lev_factor*imagenoise_xray1*levs1)

pixelscale = 1.1 #arcsec/pixels
bmaj = 6.4/pixelscale # in pixels
bmin = 6/pixelscale
bpa = 50.5 #in degrees
distance_mpc = cosmo.angular_diameter_distance(z=0.296).value #Mpc
#	print(distance_mpc)
physical_scale_mpc = 0.5 #Mpc
# Calculate the angular size of the scale bar in degrees
angular_size_deg = (physical_scale_mpc / distance_mpc) * (180.0 / np.pi)
# Convert the angular size to arcseconds
angular_size_arcsec = angular_size_deg * 3600.0
#	print(angular_size_arcsec/1.1)
scale_length_pixels = angular_size_arcsec/pixelscale

fig_width = 4.134 # half the width of A4 page in inches
fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
fig = plt.figure(figsize=(fig_width, fig_height))
ax = plt.subplot(projection=xray_wcs, slices=('x', 'y'))
map1 = ax.imshow(np.clip(xray_data, 1e-9, 1), cmap='inferno',norm=colors.LogNorm(vmin=1e-9, vmax=3e-7), interpolation='nearest')
cbar = fig.colorbar(map1, pad=0.01)
cbar.set_label(label = r'Photon counts [cm$^{-2}$ s$^{-1}$]', size=10, labelpad=-0.5)
cbar.ax.tick_params(labelsize=8)
ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
ax.tick_params(axis='both', which='major', labelsize=8)
#ax.set_xlim(700, 1100)
#ax.set_ylim(700, 1100)
ax.set_xlim(600, 1150)
ax.set_ylim(550, 1100)
#		ax.set_title(f'X-ray (${{\sigma_x =}}$ {(imagenoise_xray*1e6):.2f} ${{\mu}}$Jy)')
plt.savefig(f'{save_loc}/X-ray.pdf', dpi=200)
plt.close()


fig_width = 4.134 # half the width of A4 page in inches
fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
fig = plt.figure(figsize=(fig_width, fig_height))
ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
map1 = ax.imshow(1e6*np.clip(I_data, 1e-8, 1), cmap='viridis',norm=colors.LogNorm(vmin=3, vmax=1500), interpolation='nearest')
cbar = fig.colorbar(map1, pad=0.01)
cbar.set_label(label = r'Intensity [$\mu$Jy beam$^{-1}$]', size=10, labelpad=-0.5)
cbar.ax.tick_params(labelsize=8)
bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='white')
ax.add_artist(bar)
ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
ax.tick_params(axis='both', which='major', labelsize=8)
ax.set_xlim(3850, 4350)
ax.set_ylim(3800, 4300)
beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin, height=bmaj, angle=bpa, linewidth=0.5, facecolor=palette(3), edgecolor=palette(2))
ax.add_patch(beam_ellipse)

ax.contour(smooth_xray, levels=levelsXray, colors=['white'], alpha = 1, linewidths = 0.1, transform=ax.get_transform(xray_wcs), extent=[0, xray_data.shape[1], 0, xray_data.shape[0]])
#ax.contour(smooth_xray1, levels=levelsXray1, colors=['red'], alpha = 1, linewidths = 0.2, transform=ax.get_transform(xray_wcs), extent=[0, xray_data.shape[1], 0, xray_data.shape[0]])
#		ax.set_title(r'Stokes I')
plt.tight_layout()
plt.savefig(f'{save_loc}/StokesI_with_xray_contour.pdf', dpi=200)
plt.close()



