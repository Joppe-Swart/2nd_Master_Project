"""

Author: Joppe Swart
Created: April 2023
Last modified: April 2023
Description: This script calculate the correction factor for each channel from the brightest source.

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
import scipy.constants as c
from matplotlib import style
import scipy.stats as stats
import emcee
import corner
import astropy.units as u
from matplotlib import cm
from matplotlib.colors import ListedColormap
import os
palette = cm.get_cmap("tab20", 20)

plt.style.use(['science'])

def make_fits_cube(channels):
	teller = 0
	teller_rms = 0
	Q_hdu = fits.HDUList()
	U_hdu = fits.HDUList()
	I_hdu = fits.HDUList()
	
	for i in range(channels):
		if teller in [20, 21, 23, 36,37,38,39,40,73,74,75]:
			teller +=1
			continue
		if teller == channels:
			break
		print(f'appending channel {i}')
		Q_image_hdu = fits.open(f"DATA/stokes_q_corr/{teller:04d}-Q-image-pb.smoothed.fits")[0]
		U_image_hdu = fits.open(f"DATA/stokes_u_corr/{teller:04d}-U-image-pb.smoothed.fits")[0]
		I_image_hdu = fits.open(f"DATA/stokes_i_corr/{teller:04d}-I-image-pb.smoothed.fits")[0]


		Q_hdu.append(Q_image_hdu)
		U_hdu.append(U_image_hdu)
		I_hdu.append(I_image_hdu)

		teller += 1

	Q_hdu.writeto('DATA/Q_cube.fits', overwrite=True)
	U_hdu.writeto('DATA/U_cube.fits', overwrite=True)
	I_hdu.writeto('DATA/I_cube.fits', overwrite=True)

make_fits_cube(126)




