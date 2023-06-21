"""

Author: Joppe Swart
Created: March 2023
Last modified: March 2023
Description: This script helpf by making nice figures of Stokes I images

"""

import pyregion
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Ellipse
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
from matplotlib import style
from matplotlib import cm, colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1.anchored_artists import (AnchoredSizeBar)
palette = cm.get_cmap("tab20", 20)
import os
import requests
from time import sleep,time

from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

plt.style.use(['science'])

#def uncertainty_cmap():
#	colors = [(0.0, 'indigo'),
#		(1/6, 'blueviolet'),
#		(2/6, 'royalblue'),
#		(3/6, 'turquoise'),
#		(4/6, 'green'),
#		(5/6, 'yellow'),
#		(1.0, 'darkred')]


#	return LinearSegmentedColormap.from_list('cyclic_rainbow', colors)

def uncertainty_cmap():
	colors = [(0.0, 'indigo'),
		(1/3, 'magenta'),
		(2/3, 'yellow'),
		(1.0, 'red')]
	return LinearSegmentedColormap.from_list('cyclic_rainbow', colors)

def cyclic_rainbow_cmap():
	colors = [(0.0, 'red'),
		(1/8, 'blue'),
		(3/8, 'cyan'),
		(5/8, 'lightgreen'),
		(7/8, 'orange'),
		(1.0, 'red')]


	return LinearSegmentedColormap.from_list('cyclic_rainbow', colors)


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
	print('std = ', rmsold)
	diff=1e-1
	cut=3
	bins=np.arange(np.nanmin(m),np.nanmax(m),(np.nanmax(m)-np.nanmin(m))/30.)
	med=np.nanmedian(m)
	mean = np.nanmean(m)
	#print('the medians are', med, mean)
	for i in range(10):
		ind=np.where(np.abs(m-med)<rmsold*cut)[0]
		rms=np.nanstd(m[ind])
		print(i, 'std = ', rms)
		if np.abs((rms-rmsold)/rmsold)<diff: break
		rmsold=rms
	return rms

def Ricean_bias(save_loc='DATA/polarisation_maps'):
	if not os.path.exists(f'{save_loc}'):
		os.makedirs(f'{save_loc}')

	# Pb corrections to find a suitable mask for the noise calculation
	pb_hdu = fits.open('DATA/0055-I-pb_model.fits')[0]
	pb_data = pb_hdu.data
	pbmask = np.full(pb_data.shape, np.nan)
	pbmask[np.where(pb_data>0.9)] = 1

	# Import the polarised intensity map
	polint_hdu = fits.open('DATA/rmsynth/Results_full_field_polint.fits')[0]
	polint_header = polint_hdu.header
	polint_data = polint_hdu.data
	np.reshape(polint_data, pb_data.shape)

	#Correct for the ricean bias and make a corrected fits file
	Ricean_bias = np.sqrt(np.nanmean(polint_data[3835:3890,4275:4330])**2)
	print(f"the correction is given by {Ricean_bias}")
	polint_data_corr = polint_data
	polint_data_corr[np.where(polint_data<4*Ricean_bias)] = polint_data_corr[np.where(polint_data<4*Ricean_bias)] - Ricean_bias
	#polint_data_corr[np.where((polint_data>=1.5*Ricean_bias) & (polint_data<4))] = np.sqrt(polint_data_corr[np.where((polint_data>=1.5*Ricean_bias) & (polint_data<4))]**2 - (1.5*Ricean_bias)**2)
	polint_data_corr[np.where(polint_data>4*Ricean_bias)] = np.sqrt(polint_data_corr[np.where(polint_data>4*Ricean_bias)]**2 - 2.3*Ricean_bias**2)
	hdu_PI_corr = fits.PrimaryHDU()
	hdu_PI_corr.header = polint_header
	hdu_PI_corr.data = polint_data_corr
	hdu_PI_corr.writeto(f"{save_loc}/polint_Ricean_corr.fits", overwrite=True)


def make_polint_fits(save_loc='DATA/polarisation_maps'):
	if not os.path.exists(f'{save_loc}'):
		os.makedirs(f'{save_loc}')

	polint_data = fits.open(f"{save_loc}/polint_Ricean_corr.fits")[0].data
	I_hdu = fits.open(f"DATA/bullet_cluster_pb_corr.smoothed.fits")[0]
	I_header = I_hdu.header
	I_data = I_hdu.data

	I_ori_hdu = fits.open('DATA/selfcal_bulletcluster_003-MFS-image.fits')[0]
	I_ori_header = I_ori_hdu.header
	I_ori_data = I_ori_hdu.data

	polint_ori_data = fits.open('DATA/rmsynth/Results_full_field_polint.fits')[0].data
	#imagenoise_polint = np.sqrt(np.nanmean(polint_data[3835:3890,4275:4330])**2)

	# Pb corrections to find a suitable mask for the noise calculation
	pb_hdu = fits.open('DATA/0055-I-pb_model.fits')[0]
	pb_data = pb_hdu.data
	pbmask = np.full(pb_data.shape, np.nan)
	pbmask[np.where(pb_data>0.9)] = 1
	np.reshape(polint_data, pb_data.shape)

	# Calculate the noise
	imagenoise_I =  findrms(I_data*pbmask)
	imagenoise_polint = findrms(polint_data*pbmask)

	print(f"the noise of the polarised intensity map is {imagenoise_polint}")
	print(f"the noise of Stokes I is {imagenoise_I}")

	# Make the polarised maps 
	hdu = fits.PrimaryHDU()
	hdu.header = I_header
	hdu.data = 100*polint_data/I_data
	hdu.writeto(f"{save_loc}/polarization_fraction_map.fits", overwrite=True)

	#Make a mask to determine which pixels are relevant
	mask=np.full(I_data.shape, np.nan)
	mask[np.where((I_data>5*imagenoise_I) & (hdu.data>0.5) & (polint_data>5*imagenoise_polint))] =1

	hdu_masked_polfrac = fits.PrimaryHDU()
	hdu_masked_polfrac.header = I_header
	hdu_masked_polfrac.data = hdu.data*mask
	hdu_masked_polfrac.writeto(f"{save_loc}/polarization_fraction_map_masked.fits", overwrite=True)

	hdu_upper = fits.PrimaryHDU()
	hdu_upper.header = I_header
	hdu_upper.data = 100*5*imagenoise_polint/I_data
	hdu_upper.writeto(f"{save_loc}/upper_limitmap.fits", overwrite=True)

	# Make a mask to only use pixels with 3 sigma detection
	upper_mask=np.full(I_data.shape, np.nan)
	upper_mask[np.where((I_data>5*imagenoise_I))] =1

	hdu_upper_masked = fits.PrimaryHDU()
	hdu_upper_masked.header = I_header
	hdu_upper_masked.data = 100*5*imagenoise_polint/I_data*upper_mask
	hdu_upper_masked.writeto(f"{save_loc}/upper_limitmap_masked.fits", overwrite=True)

def make_plots(save_loc='Results/Polarisation_maps', I = True, polint=True, polfrac = True, Xray = True, optical = True, diffem = True):
	if not os.path.exists(f'{save_loc}'):
		os.makedirs(f'{save_loc}')

	# Define the beam of the unsmoothed field
	pixelscale = 1.1 #arcsec/pixels
	bmaj = 6.4/pixelscale # in pixels
	bmin = 6/pixelscale
	bpa = 50.5 #in degrees

	bmaj_convolved = 10.4/pixelscale
	bmin_convolved = 9.6/pixelscale
	bpa_convolved = 61.1
	
	distance_mpc = cosmo.angular_diameter_distance(z=0.296).value #Mpc
#	print(distance_mpc)
	physical_scale_mpc = 0.5 #Mpc
	# Calculate the angular size of the scale bar in degrees
	angular_size_deg = (physical_scale_mpc / distance_mpc) * (180.0 / np.pi)
	# Convert the angular size to arcseconds
	angular_size_arcsec = angular_size_deg * 3600.0
#	print(angular_size_arcsec/1.1)
	scale_length_pixels = angular_size_arcsec/pixelscale
	
	# Pb corrections to find a suitable mask for the noise calculation
	pb_data = fits.open('DATA/0055-I-pb_model.fits')[0].data
	pbmask = np.full(pb_data.shape, np.nan)
	pbmask[np.where(pb_data>0.9)] = 1

	I_ori_hdu = fits.open('DATA/selfcal_bulletcluster_003-MFS-image.fits')[0]
	I_ori_header = I_ori_hdu.header
	I_ori_data = I_ori_hdu.data
	imagenoise_I_ori = findrms(I_ori_data*pbmask)

	# Import the Stokes I radio image
	I_hdul = fits.open(f"DATA/bullet_cluster_pb_corr.smoothed.fits")
	I_hdu = I_hdul[0]
	I_header = I_hdu.header
	I_data = I_hdu.data
	I_wcs = WCS(flatten(I_hdul).header)

	# Import the polarised intensity map
	polint_corr_hdu = fits.open('DATA/polarisation_maps/polint_Ricean_corr.fits')[0]
	polint_data_corr = polint_corr_hdu.data

	polint_hdu = fits.open('DATA/rmsynth/Results_full_field_polint.fits')[0]
	polint_data = polint_hdu.data

	# Import the upper limit map
	upper_limitmap_hdu = fits.open('DATA/polarisation_maps/upper_limitmap_masked.fits')[0]
	upper_limitmap_data = upper_limitmap_hdu.data

	# Import the masked polarisation fraction map
	pol_frac_hdu_masked = fits.open('DATA/polarisation_maps/polarization_fraction_map_masked.fits')[0]
	pol_frac_data_masked = pol_frac_hdu_masked.data

	# Import the polarisation map
	pol_frac_hdu = fits.open('DATA/polarisation_maps/polarization_fraction_map.fits')[0]
	pol_frac_data = pol_frac_hdu.data


	# Import the X ray image
	xray_hdu = fits.open('DATA/bullet_0.5-2.0_flux.img')[0]
	xray_header = xray_hdu.header
	xray_data = xray_hdu.data
	xray_wcs = WCS(xray_header)

	# Import the optical image
	optical_hdu = fits.open('DATA/cutout_104.5833_-55.8250.fits')
	optical_header = optical_hdu[0].header
	optical_data = optical_hdu[0].data
	optical_wcs = WCS(optical_header)

	# Import the region
	region_bullet = "Regions/11_bullet_regions.reg"
	r_bullet = pyregion.open(region_bullet).as_imagecoord(I_header)
	patch_list_bullet, artist_list_bullet = r_bullet.get_mpl_patches_texts()
	
	region_diffem = "Regions/5_diffem_regions.reg"
	r_diffem = pyregion.open(region_diffem).as_imagecoord(I_header)
	patch_list_diffem, artist_list_diffem = r_diffem.get_mpl_patches_texts()
	#print(patch_list)
	#artist_list = ["1","2","3","4","5","6","7","8","9","10"]

	region_correction = "Correction/5brightsources.reg"
	r_correction = pyregion.open(region_correction).as_imagecoord(I_header)
	patch_list_correction, artist_list_correction = r_correction.get_mpl_patches_texts()

	#Calculate the noise
	imagenoise_I = findrms(I_data*pbmask)#findrms(I_data_pb[0],lower = 0, upper = 8192)#background_noise(I_data_pb*pbmask)#
	print('the noise of stokes I pb is given by', imagenoise_I)

	imagenoise_polint = findrms(polint_data_corr*pbmask)#imagenoise_polint =
	print('the image noise of the polarised intensity is given by', imagenoise_polint)

	imagenoise_xray = findrms(xray_data)
	print(f"the image noise of the xray image is {imagenoise_xray}")

	# Make the levels for the radio contours
	lev_factor = 5.
	lev_radio = np.sqrt([1, 4, 16, 64, 256, 1024, 4096, 16384, 65536]) #5, 25,100
	level_radio = np.ndarray.tolist(lev_factor*imagenoise_I*lev_radio)

	lev_xray = np.sqrt([25]) # 5
	levelsXray = np.ndarray.tolist(imagenoise_xray*lev_xray)

	if I == True:# Make the Stokes I images
		fig_width = 4.134 # half the width of A4 page in inches
		fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
		fig = plt.figure(figsize=(2*fig_width, 2*fig_height))
		ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
		map1 = ax.imshow(1e6*(np.clip(I_ori_data[0,0,:,:], 1e-8, 1)), cmap='viridis',norm=colors.LogNorm(vmin=3, vmax=300), interpolation='nearest')
		ax.set_xlim(3850, 4430)
		ax.set_ylim(3850, 4600)
		plt.savefig(f'{save_loc}/Bullet_StokesI_titel_page.pdf', dpi=200)
		bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='white')
		ax.add_artist(bar)
		cbar = fig.colorbar(map1, pad=0.01)
		cbar.set_label(label = r'Intensity [$\mu$Jy beam$^{-1}$]', size=10, labelpad=-1.5)
		cbar.ax.tick_params(labelsize=8)
		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
		ax.tick_params(axis='both', which='major', labelsize=8)
		beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
		beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
		beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin, height=bmaj, angle=bpa, linewidth=0.5, facecolor=palette(3), edgecolor=palette(2))
		ax.add_patch(beam_ellipse)
		plt.tight_layout()
#		ax.set_title(f'Stokes I (${{\sigma_I =}}$ {(imagenoise_I_ori*1e6):.2f} ${{\mu}}$Jy)')
		plt.savefig(f'{save_loc}/Bullet_StokesI_zoom.pdf', dpi=200)
		plt.arrow(4150,4400,100,70,color='white', head_width = 10, head_length = 5)
		plt.text(4085, 4375, 'Diffuse emission', fontsize=10, color='white')
		plt.arrow(3920,4200,40,-90,color='white', head_width = 10, head_length = 5)
		plt.text(3900, 4205, 'Relic', fontsize=10, color='white')
		plt.arrow(4200,4300,-40,-120,color='white', head_width = 10, head_length = 5)
		plt.text(4180,4305, 'Halo', fontsize=10, color='white')
		plt.arrow(4060,3900,0,85,color='white', head_width = 10, head_length = 5)
		plt.text(4005,3875, 'J06587-5558', fontsize=10, color='white')
		plt.savefig(f'{save_loc}/Bullet_StokesI_zoom_arrows.pdf', dpi=150)
		plt.close()
		
		fig = plt.figure(figsize=(fig_width, fig_height))
		ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
		map1 = ax.imshow(1e6*(np.clip(I_ori_data[0,0,:,:], 1e-8, 1)), cmap='viridis',norm=colors.LogNorm(vmin=3, vmax=300), interpolation='nearest')
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
		for i, patch in enumerate(patch_list_bullet):
			patch.set_linewidth(0.7)
#			patch.set_edgecolor(palette(4))
			ax.add_patch(patch)
			# calculate center of polygon using vertices
			if isinstance(patch, patches.Polygon):
				vertices = patch.get_path().vertices
				x = [p[0] for p in vertices]
				y = [p[1] for p in vertices]
				center = (np.sum(x)/len(vertices), np.sum(y)/len(vertices))
			elif isinstance(patch, patches.Circle):
				center = (patch.center[0], patch.center[1])
			elif isinstance(patch, patches.Ellipse):
				center = (patch.center[0], patch.center[1])
			else:
				raise TypeError(f"Unsupported patch type {type(patch)}")
			artist_list_bullet.append(ax.text(center[0], center[1], str(i+1), color='red')) # add number to artist list

		for t in artist_list_bullet:
			t.set_fontsize(10)
			ax.add_artist(t)
		plt.tight_layout()
		plt.savefig(f'{save_loc}/Bullet_StokesI_regions.pdf', dpi=200)
		plt.close()

		fig = plt.figure(figsize=(fig_width, fig_height))
		ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
		map1 = ax.imshow(1e6*(np.clip(I_ori_data[0,0,:,:], 1e-8, 1)), cmap='viridis',norm=colors.LogNorm(vmin=3, vmax=50), interpolation='nearest')
		cbar = fig.colorbar(map1, pad=0.01)
		cbar.set_label(label = r'Intensity [$\mu$Jy beam$^{-1}$]', size=10, labelpad=-0.5)
		cbar.ax.tick_params(labelsize=8)
		bar = AnchoredSizeBar(ax.transData, 0.5*scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='white')
		ax.add_artist(bar)
		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
		ax.tick_params(axis='both', which='major', labelsize=8)
		plt.ylim(4400, 4630)
		plt.xlim(4120, 4350)
		beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
		beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
		beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin, height=bmaj, angle=bpa, linewidth=0.5, facecolor=palette(3), edgecolor=palette(2))
		ax.add_patch(beam_ellipse)
		for i, patch in enumerate(patch_list_diffem):
			patch.set_linewidth(0.7)
#			patch.set_edgecolor(palette(4))
			ax.add_patch(patch)
			# calculate center of polygon using vertices
			if isinstance(patch, patches.Polygon):
				vertices = patch.get_path().vertices
				x = [p[0] for p in vertices]
				y = [p[1] for p in vertices]
				center = (np.sum(x)/len(vertices), np.sum(y)/len(vertices))
			elif isinstance(patch, patches.Circle):
				center = (patch.center[0], patch.center[1])
			elif isinstance(patch, patches.Ellipse):
				center = (patch.center[0], patch.center[1])
			else:
				raise TypeError(f"Unsupported patch type {type(patch)}")
			artist_list_diffem.append(ax.text(center[0], center[1], str(i+10), color='red')) # add number to artist list

		for t in artist_list_diffem:
			t.set_fontsize(10)
			ax.add_artist(t)
		plt.tight_layout()
		plt.savefig(f'{save_loc}/diffem_StokesI_regions.pdf', dpi=200)
		plt.close()

		fig_width = 2*4.134 # half the width of A4 page in inches
		fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
		fig = plt.figure(figsize=(fig_width, fig_height))
		ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
		map1 = ax.imshow(1e6*(np.clip(I_ori_data[0,0,:,:], 1e-8, 1)), cmap='viridis',norm=colors.LogNorm(vmin=3, vmax=300), interpolation='nearest')
		cbar = fig.colorbar(map1, pad=0.01)
		cbar.set_label(label = r'Intensity [$\mu$Jy beam$^{-1}$]', size=10)#, labelpad=-1.5)
		cbar.ax.tick_params(labelsize=8)
		bar = AnchoredSizeBar(ax.transData, 2*scale_length_pixels, '1 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='white')
		ax.add_artist(bar)
		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
		ax.tick_params(axis='both', which='major', labelsize=8)
		beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
		beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
		beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin, height=bmaj, angle=bpa, linewidth=0.5, facecolor=palette(3), edgecolor=palette(2))
		ax.add_patch(beam_ellipse)
		plt.xlim(3096, 5096)
		plt.ylim(3096, 5096)
#		plt.title(f'Stokes I (${{\sigma_I =}}$ {(imagenoise_I_ori*1e6):.2f} ${{\mu}}$Jy)')
		plt.savefig(f'{save_loc}/Bullet_StokesI_field.pdf', dpi=200)
		plt.close()

		fig = plt.figure(figsize=(fig_width, fig_height))
		ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
		map1 = ax.imshow(1e6*(np.clip(I_ori_data[0,0,:,:], 1e-8, 1)), cmap='viridis',norm=colors.LogNorm(vmin=3, vmax=10000), interpolation='nearest')
		bar = AnchoredSizeBar(ax.transData, 2*scale_length_pixels, '1 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='white')
		ax.add_artist(bar)
		cbar = fig.colorbar(map1, pad=0.01)
		cbar.set_label(label = r'Intensity [$\mu$Jy beam$^{-1}$]', size=10, labelpad=-0.5)
		cbar.ax.tick_params(labelsize=8)
		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
		ax.tick_params(axis='both', which='major', labelsize=8)
		beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
		beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
		beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin, height=bmaj, angle=bpa, linewidth=0.5, facecolor=palette(3), edgecolor=palette(2))
		ax.add_patch(beam_ellipse)
		plt.xlim(3096, 5096)
		plt.ylim(3096, 5096)
		for i, patch in enumerate(patch_list_correction):
			patch.set_linewidth(0.7)
#			patch.set_edgecolor(palette(4))
			ax.add_patch(patch)
			# calculate center of polygon using vertices
			if isinstance(patch, patches.Polygon):
				vertices = patch.get_path().vertices
				x = [p[0] for p in vertices]
				y = [p[1] for p in vertices]
				center = (np.sum(x)/len(vertices), np.sum(y)/len(vertices))
			elif isinstance(patch, patches.Circle):
				center = (patch.center[0], patch.center[1])
			elif isinstance(patch, patches.Ellipse):
				center = (patch.center[0], patch.center[1])
			else:
				raise TypeError(f"Unsupported patch type {type(patch)}")
			artist_list_correction.append(ax.text(center[0], center[1], str(i+1), color='red')) # add number to artist list

		for t in artist_list_correction:
			t.set_fontsize(12)
			ax.add_artist(t)
		plt.tight_layout()
		plt.savefig(f'{save_loc}/Bullet_StokesI_field_correction_regions.pdf', dpi=200)
		plt.close()


	if polint == True:# Make the Stokes I images
		fig_width = 4.134 # half the width of A4 page in inches
		fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
		fig = plt.figure(figsize=(2*fig_width, 2*fig_height))
		ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
		map1 = plt.imshow(1e6*np.clip(polint_data_corr, 1e-8, 1), cmap='viridis',norm=colors.LogNorm(vmin=1, vmax=50), interpolation='nearest')
		cbar = fig.colorbar(map1, pad=0.01)
		cbar.set_label(label = r'Polarised intensity [$\mu$Jy beam$^{-1}$]', size=10, labelpad=-0.5)
		cbar.ax.tick_params(labelsize=8)
		bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='white')
		ax.add_artist(bar)
		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
		ax.tick_params(axis='both', which='major', labelsize=8)
		ax.set_xlim(3850, 4430)
		ax.set_ylim(3850, 4600)
		beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
		beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
		beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5, facecolor=palette(3), edgecolor=palette(2))
		ax.add_patch(beam_ellipse)
#		plt.title(f'Polarised intensity (${{\sigma_p =}}$ {(imagenoise_polint*1e6):.2f} ${{\mu}}$Jy)')
		plt.tight_layout()
		plt.savefig(f'{save_loc}/Bullet_polarised_intensity.pdf', dpi=200)
		plt.close()

	if polfrac == True: # Make the polarised fraction maps
		fig_width = 4.134 # half the width of A4 page in inches
		fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
		fig = plt.figure(figsize=(fig_width, fig_height))
		ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
		map1 = plt.imshow(pol_frac_data_masked[0,0,:,:], cmap='rainbow',vmin=0, vmax=20, interpolation='nearest')
		cbar = fig.colorbar(map1, pad=0.01)
		cbar.ax.tick_params(labelsize=8)
		cbar.set_label(label = r'Polarisation fraction [\%]', size=10)#, labelpad=-0.5)
		cbar.ax.tick_params(labelsize=8)
		bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
		ax.add_artist(bar)
		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
		ax.tick_params(axis='both', which='major', labelsize=8)
		plt.xlim(3850, 4350)
		plt.ylim(3800, 4300)
#		plt.xlim(3496, 4696)
#		plt.ylim(3496, 4696)
		beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
		beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
		beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
		ax.add_patch(beam_ellipse)
#		plt.title(r'Polarisation fraction')
		plt.tight_layout()
		plt.savefig(f'{save_loc}/Polarisation_fraction.pdf', dpi=200)
		plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
		plt.savefig(f'{save_loc}/Polarisation_fraction_contour.pdf', dpi=200)
		plt.close()


		fig_width = 4.134 # half the width of A4 page in inches
		fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
		fig = plt.figure(figsize=(fig_width, fig_height))
		ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
		map1 = plt.imshow(upper_limitmap_data[0,0,:,:], cmap='rainbow',vmin=0, vmax=20, interpolation='nearest')
		cbar = fig.colorbar(map1, pad=0.01)
		cbar.ax.tick_params(labelsize=8)
		cbar.set_label(label = r'Polarisation fraction [\%]', size=10)#, labelpad=-0.5)
		cbar.ax.tick_params(labelsize=8)
		bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
		ax.add_artist(bar)
		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
		ax.tick_params(axis='both', which='major', labelsize=8)
		plt.xlim(3850, 4350)
		plt.ylim(3800, 4300)
#		plt.xlim(3496, 4696)
#		plt.ylim(3496, 4696)
		beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
		beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
		beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5, facecolor='black', edgecolor='black')
		ax.add_patch(beam_ellipse)
#		plt.title(f'polarisation needed for ${{3\sigma}}$ detection')
		plt.tight_layout()
		plt.savefig(f'{save_loc}/Upper_Polarisation_fraction.pdf', dpi=200)
		plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
		plt.savefig(f'{save_loc}/Upper_Polarisation_fraction_contour.pdf', dpi=200)
		plt.close()

#	if Xray ==True:
#		fig_width = 4.134 # half the width of A4 page in inches
#		fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
#		fig = plt.figure(figsize=(fig_width, fig_height))
#		ax = plt.subplot(projection=xray_wcs, slices=('x', 'y'))
#		map1 = ax.imshow(xray_data, cmap='inferno',vmin=0, vmax=1e-7, interpolation='nearest')
#		cbar = fig.colorbar(map1, pad=0.01)
#		cbar.set_label(label = r'Intensity [Arbitray units]', size=10, labelpad=-1.5)
#		cbar.ax.tick_params(labelsize=8)
#		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
#		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
#		ax.tick_params(axis='both', which='major', labelsize=8)
#		ax.set_xlim(700, 1100)
#		ax.set_ylim(700, 1100)
##		ax.set_title(f'X-ray (${{\sigma_x =}}$ {(imagenoise_xray*1e6):.2f} ${{\mu}}$Jy)')
#		plt.savefig(f'{save_loc}/X-ray.pdf', dpi=200)
#		plt.close()


#		fig_width = 4.134 # half the width of A4 page in inches
#		fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
#		fig = plt.figure(figsize=(fig_width, fig_height))
#		ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
#		map1 = ax.imshow(1e6*I_ori_data[0,0,:,:], cmap='bone',vmin=-24, vmax=64, interpolation='nearest')
#		cbar = fig.colorbar(map1, pad=0.01)
#		cbar.set_label(label = r'$\mu$Jy/beam', size=10, labelpad=-1.5)
#		cbar.ax.tick_params(labelsize=8)
#		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
#		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
#		ax.tick_params(axis='both', which='major', labelsize=8)
#		ax.set_xlim(3850, 4350)
#		ax.set_ylim(3800, 4300)
#		beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
#		beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
#		beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin, height=bmaj, angle=bpa, linewidth=0.5, facecolor=palette(3), edgecolor=palette(2))
#		ax.add_patch(beam_ellipse)
#		ax.contour(xray_data, levels=levelsXray, colors=['r'], alpha = 1, linewidths = 0.2, transform=ax.get_transform(xray_wcs), extent=[0, xray_data.shape[1], 0, xray_data.shape[0]])
##		ax.set_title(r'Stokes I')
#		plt.tight_layout()
#		plt.savefig(f'{save_loc}/StokesI_with_xray_contour.pdf', dpi=200)
#		plt.close()

	if optical == True:
		lev_fac_op = 3.
		levs_opt = np.sqrt([1, 4, 16, 64, 256, 1024, 4096]) #5, 15, 60, 120
		level_optical = np.ndarray.tolist(lev_fac_op*imagenoise_I*levs_opt)
		
		fig_width = 4.134 # half the width of A4 page in inches
		fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
		fig = plt.figure(figsize=(fig_width, fig_height))
		ax = plt.subplot(projection=optical_wcs, slices=('x', 'y', 0))
		map1 = ax.imshow(optical_data[0,:,:], cmap='hot',vmin=-0.01, vmax=0.18, interpolation='nearest')
		cbar = fig.colorbar(map1, pad=0.01)
		cbar.set_label(label = r'Fluxscale [Arbitrary units]', size=10)#, labelpad=-0.5)
		cbar.ax.tick_params(labelsize=8)
		bar = AnchoredSizeBar(ax.transData, 0.5*1.1/0.262*scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='white')
		ax.add_artist(bar)
		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
		ax.tick_params(axis='both', which='major', labelsize=8)
		ax.set_xlim(500, 1500)
		ax.set_ylim(500, 1500)
		beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
		beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
		beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved*1.1/0.262, height=bmaj_convolved*1.1/0.262, angle=bpa_convolved, linewidth=0.5, facecolor = palette(5), edgecolor=palette(4))
		ax.add_patch(beam_ellipse)
#		ax.set_title(r'dss2 red')
		ax.contour(I_data[0,0,:,:], levels=level_optical, colors=['white'], alpha = 0.8, linewidths = 0.2, transform=ax.get_transform(I_wcs))
		plt.tight_layout()
		plt.savefig(f'{save_loc}/optical_radio_contour.pdf', dpi=200)
		plt.close()

	if diffem == True:
		mask_diffem=np.full(I_data.shape, np.nan)
		mask_diffem[np.where((I_data>3*imagenoise_I) & (pol_frac_data>0.5) & (polint_data>7*imagenoise_polint))] =1

		lev_factor = 3.
		levs = np.sqrt([1, 4, 16, 64, 256, 1024, 4096]) #3, 6, 9
		level_diffem = np.ndarray.tolist(lev_factor*imagenoise_I*levs)


		fig_width = 4.134 # half the width of A4 page in inches
		fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
		fig = plt.figure(figsize=(fig_width, fig_height))
		ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
		map1 = plt.imshow(pol_frac_data[0,0,:,:]*mask_diffem[0,0,:,:], cmap='rainbow',vmin=0, vmax=25, interpolation='nearest')
		cbar = fig.colorbar(map1, pad=0.01)
		cbar.set_label(label = r'Polarisation fraction [\%]', size=10, labelpad=-0.5)
		cbar.ax.tick_params(labelsize=8)
		bar = AnchoredSizeBar(ax.transData, 0.5*scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
		ax.add_artist(bar)
		ax.set_xlabel('RA (J2000)', fontsize=10, labelpad=0.5)
		ax.set_ylabel('DEC (J2000)',fontsize=10, labelpad=-1)
		ax.tick_params(axis='both', which='major', labelsize=8)
		plt.ylim(4400, 4630)
		plt.xlim(4120, 4350)
		beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
		beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
		beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
		ax.add_patch(beam_ellipse)
#		plt.title(r'Polarisation fraction')
		plt.tight_layout()
		plt.savefig(f'{save_loc}/Polarisation_fraction_diffem.pdf', dpi=200)
		plt.contour(I_data[0,0,:,:], levels=level_diffem, colors='black', alpha=0.4, linewidths=0.07)
		plt.savefig(f'{save_loc}/Polarisation_fraction_diffem_contour.pdf', dpi=200)
		plt.close()

#		hdu_masked_polfrac_difem = fits.PrimaryHDU()
#		hdu_masked_polfrac_difem.header = I_header
#		hdu_masked_polfrac_difem.data = pol_frac_data*mask_diffem
#		hdu_masked_polfrac_difem.writeto(f"DATA/polarisation_maps/diffem_polarization_fraction_map_masked.fits", overwrite=True)

def plotRMTF(rmtfile, saveloc):
	rmtf =	 np.loadtxt(rmtfile,delimiter=' ')
	fig_width = 4.134 # half the width of A4 page in inches
	fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
	fig = plt.figure(figsize=(fig_width, fig_height))
	plt.plot(rmtf[:,0],rmtf[:,1], linewidth=0.5, color=palette(0), label='Real',ls='dashed')
	plt.plot(rmtf[:,0],rmtf[:,2], linewidth=0.5,color=palette(2),label='Imaginary',ls='dashed')
	plt.plot(rmtf[:,0],np.sqrt(rmtf[:,1]**2+rmtf[:,2]**2), linewidth=0.5,color=palette(4), label='Amplitude')
	plt.legend(fontsize=10,loc='best', frameon=True, shadow= True)
	plt.xticks(fontsize=8)
	plt.yticks(fontsize=8)
	plt.xlabel(r'$\Phi$ [rad m$^{-2}$]',fontsize=10, labelpad=0.5)
	plt.ylabel(r'$|$F$|$ [Jy beam$^{-1}$ RMSF$^{-1}$]',fontsize=10, labelpad=0.5)
	plt.tight_layout()
	plt.savefig(f'{saveloc}/rmtf.pdf', dpi=150)
	plt.show()
	plt.close()

#Ricean_bias()
#make_polint_fits()
make_plots()
#plotRMTF('/net/rijn9/data2/swart/DATA/MeerKAT_DATA/Analysis/DATA/rmsynth/Results_full_field_rmsf.txt','Results/Polarisation_maps')





