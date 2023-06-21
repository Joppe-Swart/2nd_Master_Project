"""

Author: Joppe Swart
Created: May 2023
Last modified: May 2023
Description: This script helpf by making nice figures of the fitted images

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
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
palette = cm.get_cmap("tab20", 20)
import os
import requests
from time import sleep,time
from mpl_toolkits.axes_grid1.anchored_artists import (AnchoredSizeBar)
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

plt.style.use(['science'])

def cyclic_rainbow_cmap():
	colors = [(0.0, 'red'),
		(1/8, 'blue'),
		(3/8, 'cyan'),
		(5/8, 'lightgreen'),
		(7/8, 'orange'),
		(1.0, 'red')]


	return LinearSegmentedColormap.from_list('cyclic_rainbow', colors)

def uncertainty_cmap():
	colors = [(0.0, 'indigo'),
		(1/2, 'magenta'),
		(1.0, 'yellow')]
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

def plot_bullet(channels, x_min = 3950, x_max = 3955, y_min = 4220, y_max = 4225, loc=f'Results/Images_Fitted_Pixels/all_pixels_final'):
	"""
	plot: Function that loops over all pixels and fits the p0, chi, rm and sigma rm
	"""
	I_hdul = fits.open(f"DATA/bullet_cluster_pb_corr.smoothed.fits")
	I_hdu = I_hdul[0]
	I_header = I_hdu.header
	I_data = I_hdu.data
	I_wcs = WCS(flatten(I_hdul).header)
	
	pb_data = fits.open('DATA/0055-I-pb_model.fits')[0].data
	pbmask = np.full(pb_data.shape, np.nan)
	pbmask[np.where(pb_data>0.9)] = 1
	

	imagenoise_I = findrms(I_data*pbmask)

	lev_factor = 5.
	lev_radio = np.sqrt([1, 4, 16, 64, 256, 1024, 4096]) #5, 25,100
	level_radio = np.ndarray.tolist(lev_factor*imagenoise_I*lev_radio)


	datap = fits.open(f'{loc}/fitted_p0.fits')[0].data
	datachi = fits.open(f'{loc}/fitted_chi0.fits')[0].data
	datarm = fits.open(f'{loc}/fitted_rm.fits')[0].data
	datasigmarm = fits.open(f'{loc}/fitted_sigmarm.fits')[0].data
	dataalpha = fits.open(f'{loc}/fitted_a.fits')[0].data
	dataalpha_uncrt = fits.open(f'{loc}/fitted_a_uncrt.fits')[0].data
	datap_uncrt = fits.open(f'{loc}/fitted_p0_uncrt.fits')[0].data
	datachi_uncrt = fits.open(f'{loc}/fitted_chi0_uncrt.fits')[0].data
	datarm_uncrt = fits.open(f'{loc}/fitted_rm_uncrt.fits')[0].data
	datasigmarm_uncrt = fits.open(f'{loc}/fitted_sigmarm_uncrt.fits')[0].data

	loc=loc+"/figures"
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

	fig_width = 4.134 # half the width of A4 page in inches
	fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(100*datap[0,0,:,:], cmap='rainbow',vmin=0, vmax=60, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'$p_0$ [\%]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted Polarisation fraction')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_p0.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_p0_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datarm[0,0,:,:], cmap='rainbow',vmin=-70, vmax=20, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'RM [rad m$^{-2}$]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted RM')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_RM.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_RM_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datachi[0,0,:,:], cmap=cyclic_rainbow_cmap(),vmin=0, vmax=3.14, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'$\chi_0$ [rad]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\chi$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_chi.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_chi_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datasigmarm[0,0,:,:]**0.5, cmap='rainbow',vmin=0, vmax=30, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'$\sigma_{\mathrm{RM}}$ [rad m$^{-2}$]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\sigma_{RM}$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_sigmarm.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_sigmarm_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(dataalpha[0,0,:,:], cmap='rainbow',vmin=-2, vmax=0, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'$\alpha$', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\sigma_{RM}$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_a.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_a_contour.pdf', dpi=150)
	plt.close()

# Now the uncertainty maps
	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(dataalpha_uncrt[0,0,:,:], cmap=uncertainty_cmap(),vmin=0, vmax=2, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'Uncertainty $\alpha$', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\sigma_{RM}$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_a_uncrt.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_a_uncrt_contour.pdf', dpi=150)
	plt.close()

	fig_width = 4.134 # half the width of A4 page in inches
	fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(100*datap_uncrt[0,0,:,:], cmap=uncertainty_cmap(),vmin=0, vmax=20, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'Uncertainty $p_0$ [\%]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted Polarisation fraction')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_p0_uncrt.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_p0_uncrt_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datarm_uncrt[0,0,:,:], cmap=uncertainty_cmap(),vmin=0, vmax=15, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'Uncertainty RM [rad m$^{-2}$]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted RM')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_RM_uncrt.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_RM_uncrt_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datachi_uncrt[0,0,:,:], cmap=uncertainty_cmap(),vmin=0, vmax=1, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'Uncertainty $\chi_0$ [rad]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\chi$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_chi_uncrt.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_chi_uncrt_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datasigmarm_uncrt[0,0,:,:]**0.5, cmap=uncertainty_cmap(),vmin=0, vmax=15, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'Uncertainty $\sigma_{\mathrm{RM}}$ [rad m$^{-2}$]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.5 Mpc', 2,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\sigma_{RM}$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_sigmarm_uncrt.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_sigmarm_uncrt_contour.pdf', dpi=150)
	plt.close()



def plot_diffem(channels, x_min = 3950, x_max = 3955, y_min = 4220, y_max = 4225, loc=f'Results/Images_Fitted_Pixels/diffem_without_chi'):
	"""
	plot: Function that loops over all pixels and fits the p0, chi, rm and sigma rm
	"""
	I_hdul = fits.open(f"DATA/bullet_cluster_pb_corr.smoothed.fits")
	I_hdu = I_hdul[0]
	I_header = I_hdu.header
	I_data = I_hdu.data
	I_wcs = WCS(flatten(I_hdul).header)
	
	pb_data = fits.open('DATA/0055-I-pb_model.fits')[0].data
	pbmask = np.full(pb_data.shape, np.nan)
	pbmask[np.where(pb_data>0.9)] = 1
	

	imagenoise_I = findrms(I_data*pbmask)

	lev_factor = 3.
	lev_radio = np.sqrt([1, 4, 16, 64, 256, 1024, 4096]) #5, 25,100
	level_radio = np.ndarray.tolist(lev_factor*imagenoise_I*lev_radio)


	datap = fits.open(f'{loc}/fitted_p0.fits')[0].data
	datachi = fits.open(f'{loc}/fitted_chi0.fits')[0].data
	datarm = fits.open(f'{loc}/fitted_rm.fits')[0].data
	datasigmarm = fits.open(f'{loc}/fitted_sigmarm.fits')[0].data
	dataalpha = fits.open(f'{loc}/fitted_a.fits')[0].data
	dataalpha_uncrt = fits.open(f'{loc}/fitted_a_uncrt.fits')[0].data
	datap_uncrt = fits.open(f'{loc}/fitted_p0_uncrt.fits')[0].data
	datachi_uncrt = fits.open(f'{loc}/fitted_chi0_uncrt.fits')[0].data
	datarm_uncrt = fits.open(f'{loc}/fitted_rm_uncrt.fits')[0].data
	datasigmarm_uncrt = fits.open(f'{loc}/fitted_sigmarm_uncrt.fits')[0].data

	loc=loc+"/figures"
	pixelscale = 1.1 #arcsec/pixels
	bmaj = 6.4/pixelscale # in pixels
	bmin = 6/pixelscale
	bpa = 50.5 #in degrees

	bmaj_convolved = 10.4/pixelscale
	bmin_convolved = 9.6/pixelscale
	bpa_convolved = 61.1
	
	distance_mpc = cosmo.angular_diameter_distance(z=0.296).value #Mpc
#	print(distance_mpc)
	physical_scale_mpc = 0.25 #Mpc
	# Calculate the angular size of the scale bar in degrees
	angular_size_deg = (physical_scale_mpc / distance_mpc) * (180.0 / np.pi)
	# Convert the angular size to arcseconds
	angular_size_arcsec = angular_size_deg * 3600.0
#	print(angular_size_arcsec/1.1)
	scale_length_pixels = angular_size_arcsec/pixelscale

	fig_width = 4.134 # half the width of A4 page in inches
	fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(100*datap[0,0,:,:], cmap='rainbow',vmin=20, vmax=85, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'$p_0$ [\%]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted Polarisation fraction')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_p0.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_p0_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datarm[0,0,:,:], cmap='rainbow',vmin=-50, vmax=20, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'RM [rad m$^{-2}$]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted RM')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_RM.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_RM_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datachi[0,0,:,:], cmap=cyclic_rainbow_cmap(),vmin=0, vmax=3.14, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'$\chi_0$ [rad]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\chi$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_chi.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_chi_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datasigmarm[0,0,:,:]**0.5, cmap='rainbow',vmin=0, vmax=30, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'$\sigma_{\mathrm{RM}}$ [rad m$^{-2}$]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\sigma_{RM}$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_sigmarm.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_sigmarm_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(dataalpha[0,0,:,:], cmap='rainbow',vmin=-2, vmax=0, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'$\alpha$', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\sigma_{RM}$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_a.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_a_contour.pdf', dpi=150)
	plt.close()

# Now the uncertainty maps
	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(dataalpha_uncrt[0,0,:,:], cmap=uncertainty_cmap(),vmin=0, vmax=2, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'Uncertainty $\alpha$', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\sigma_{RM}$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_a_uncrt.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_a_uncrt_contour.pdf', dpi=150)
	plt.close()

	fig_width = 4.134 # half the width of A4 page in inches
	fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(100*datap_uncrt[0,0,:,:], cmap=uncertainty_cmap(),vmin=0, vmax=30, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'Uncertainty $p_0$ [\%]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted Polarisation fraction')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_p0_uncrt.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_p0_uncrt_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datarm_uncrt[0,0,:,:], cmap=uncertainty_cmap(),vmin=0, vmax=20, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'Uncertainty RM [rad m$^{-2}$]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted RM')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_RM_uncrt.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_RM_uncrt_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datachi_uncrt[0,0,:,:], cmap=uncertainty_cmap(),vmin=0, vmax=2, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'Uncertainty $\chi_0$ [rad]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\chi$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_chi_uncrt.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_chi_uncrt_contour.pdf', dpi=150)
	plt.close()

	fig = plt.figure(figsize=(fig_width, fig_height))
	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
	map1 = plt.imshow(datasigmarm_uncrt[0,0,:,:]**0.5, cmap=uncertainty_cmap(),vmin=0, vmax=20, interpolation='nearest')
	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
	cbar.set_label(label = r'Uncertainty $\sigma_{\mathrm{RM}}$ [rad m$^{-2}$]', size=9)#, labelpad=-0.5)
	cbar.ax.tick_params(labelsize=8)
	bar = AnchoredSizeBar(ax.transData, scale_length_pixels, '0.25 Mpc', 1,pad=0.5,sep=5, borderpad=0.5, frameon=False,size_vertical=0.5, color='black')
	ax.add_artist(bar)
	ax.set_xlabel('RA (J2000)', fontsize=9, labelpad=0.5)
	ax.set_ylabel('DEC (J2000)',fontsize=9, labelpad=-0.5)
	ax.tick_params(axis='both', which='major', labelsize=8)
	plt.ylim(y_min, y_max)
	plt.xlim(x_min, x_max)
	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,  facecolor='black', edgecolor='black')
	ax.add_patch(beam_ellipse)
	#plt.title(r'Fitted $\sigma_{RM}$')
	plt.tight_layout()
	plt.savefig(f'{loc}/fitted_sigmarm_uncrt.pdf', dpi=150)
	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.4, linewidths=0.07)
	plt.savefig(f'{loc}/fitted_sigmarm_uncrt_contour.pdf', dpi=150)
	plt.close()

plot_diffem(115, x_min = 4120, x_max = 4350, y_min = 4400, y_max = 4630, loc =f'Results/Images_Fitted_Pixels/diff_em_without_chi')
plot_bullet(115, x_min = 3850, x_max = 4350, y_min = 3800, y_max = 4300, loc =f'Results/Images_Fitted_Pixels/all_pixels_final')
