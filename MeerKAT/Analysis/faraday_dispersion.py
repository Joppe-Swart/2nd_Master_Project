"""

Author: Joppe Swart
Created: June 2023
Last modified: June 2023
Description: Script to make the faraday dispersion functions

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
palette = cm.get_cmap("tab20", 20)
import os
import requests
from time import sleep,time

from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
plt.style.use(['science'])	

for i in range(1, 13):
	print(f"plotting region {i}")
	data = f'/net/rijn9/data2/swart/DATA/MeerKAT_DATA/Analysis/Results/Images_Regions/region/Final_region_analysis/faraday_dispersion/faraday_dispersion_region{i}.dat'
	phi, intensity = np.loadtxt(data, unpack=True)
	fig_width = 4.134 # half the width of A4 page in inches
	fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
	fig = plt.figure(figsize=(fig_width, fig_height))
	plt.grid(color='lightgrey', linestyle=':', linewidth=0.5)
	plt.xlabel(r'RM [rad m$^{-2}$]', fontsize=10)
	plt.ylabel(r'Average polarised intensity [$\mu$Jy beam$^{-1}$]', fontsize=10)
	plt.plot(phi, 10**6*intensity, 'k--', linewidth=0.5)
	plt.tight_layout()
	plt.savefig(f'/net/rijn9/data2/swart/DATA/MeerKAT_DATA/Analysis/Results/Images_Regions/region/Final_region_analysis/faraday_dispersion/faraday_dispersion_region{i}.pdf', dpi=200)
	plt.close()

