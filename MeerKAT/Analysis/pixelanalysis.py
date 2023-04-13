"""

Author: Joppe Swart
Created: March 2023
Last modified: March 2023
Description: This script helpf by making nice figures of Stokes I images

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
#palette = cm.get_cmap("Vega20", 20)

plt.style.use(['science'])

# Make a catalog of pixels with a certain value
def make_cat(xmin_cut, xmax_cut, ymin_cut, ymax_cut, pol_per):
    # Now we need to find pixels with a cerain value from polint map
    polfrac_hdu = fits.open('polarization_fraction_map_masked.fits')[0]
    polfrac_header = polfrac_hdu.header
    polfrac_data = polfrac_hdu.data
#    print(polfrac_data.shape)
#    mask=np.zeros(polint_data.shape, dtype=bool)
#    mask[min_cut:max_cut, min_cut:max_cut] =True
#    print(np.std(polint_data*mask))
#    std = np.std(polint_data*mask)
#   plt.figure()
    # plt.subplot(projection=I_wcs, slices=('x', 'y', 0, 0))
#    map1 = plt.imshow(1e6*polint_data*mask, cmap='bone',vmin=-24, vmax=64, interpolation='nearest')
#    plt.colorbar(map1, label = r'$\mu$Jy/beam')
#    plt.xlabel('RA')
#    plt.ylabel('DEC')
#    #plt.xlim(3496, 4696)
#    #plt.ylim(3496, 4696)

#    # apertures.plot(color='white')
#    # annulus.plot(color='white')
#    plt.title(r'Stokes I')
#    #plt.savefig(r'Stokes_I.pdf')
#    plt.show()
#    rows, cols = polint_data.shape
#    print(rows, cols)
    teller = 0
    pol = []
    for i in range(xmin_cut,xmax_cut):
        for j in range(ymin_cut,ymax_cut):
            if polfrac_data[i,j] >= pol_per:
#                print(f'pixel [{i},{j}] contains a >20% polarization detection of {polfrac_data[i,j]}%')
                teller +=1
                pol.append([polfrac_data[i,j], int(i), int(j)])
    pol = np.array(pol)
#    print(pol.shape)
    np.save(f'pol_20_pixels.npy', pol)
#    print(teller)

    return pol

#make_cat()

#maak een functie die in de linpol ding zoekt
def find_pixels(channel):
    pol_I_hdu = fits.open(f'stokes_i/{channel:04d}-I-image.pb.smoothed.fits')[0]
    pol_Q_hdu = fits.open(f'stokes_q/{channel:04d}-Q-image.pb.smoothed.fits')[0]
    pol_U_hdu = fits.open(f'stokes_u/{channel:04d}-U-image.pb.smoothed.fits')[0]

    I_header = pol_I_hdu.header
    Q_header = pol_Q_hdu.header
    U_header = pol_U_hdu.header

    I_data = pol_I_hdu.data[0]
    Q_data = pol_Q_hdu.data[0]
    U_data = pol_U_hdu.data[0]

    I_wcs = WCS(I_header)
    Q_wcs = WCS(Q_header)
    U_wcs = WCS(U_header)

    if not os.path.exists('linpol'):
        os.mkdir('linpol')
    lin_pol = np.sqrt(Q_data**2+U_data**2)
#    lin_pol_freqmean = np.nanmean(lin_pol, axis=0)
    linpol_hdu = fits.PrimaryHDU(lin_pol)
    linpol_hdu.writeto(r'linpol/{channel:04d}-linpol-image.fits', overwrite=True)

    return

def find_high_pol_pix(xmin_cut, xmax_cut, ymin_cut, ymax_cut, pol_per):
    pol = pd.DataFrame(data=make_cat(xmin_cut, xmax_cut, ymin_cut, ymax_cut, pol_per))
    row = pol.loc[[pol[0].idxmax()],:]
#    print((row))
    return row

#find_high_pol_pix(3920, 3970,4170, 4240 ,20)

def make_lin_pol(channel, x, y):
    pol_I_hdu = fits.open(f'stokes_i/{channel:04d}-I-image-pb.smoothed.fits')[0]
    pol_Q_hdu = fits.open(f'stokes_q/{channel:04d}-Q-image-pb.smoothed.fits')[0]
    pol_U_hdu = fits.open(f'stokes_u/{channel:04d}-U-image-pb.smoothed.fits')[0]

    I_header = pol_I_hdu.header
    Q_header = pol_Q_hdu.header
    U_header = pol_U_hdu.header

    I_data = pol_I_hdu.data[0]
    Q_data = pol_Q_hdu.data[0]
    U_data = pol_U_hdu.data[0]

    I_wcs = WCS(I_header)
    Q_wcs = WCS(Q_header)
    U_wcs = WCS(U_header)

    freq = I_header['CRVAL3']
    wave2 = (c.c/freq)**2


    if not os.path.exists('linpol'):
        os.mkdir('linpol')
    if os.path.exists(f'linpol/{channel:04d}-linpol-image.fits'):
        linpol_hdu = fits.open(f'linpol/{channel:04d}-linpol-image.fits')[0]
        lin_pol = linpol_hdu.data
#        print(lin_pol.shape)
    if not os.path.exists(f'linpol/{channel:04d}-linpol-image.fits'):
        print(f'Making the linear polarization fits for channel {channel}')
        lin_pol = np.sqrt(Q_data[0]**2+U_data[0]**2)# - 8.8e-6
        linpol_hdu = fits.PrimaryHDU(lin_pol)
        linpol_hdu.writeto(f'linpol/{channel:04d}-linpol-image.fits', overwrite=True)

#    if not os.path.exists('polfrac'):
#        os.mkdir('polfrac')
#    print(f'Making the polarization fraction fits for channel {channel}')
#    np.clip(I_data, 1e-8,1)
#    pol_frac = lin_pol/I_data[0]
#    pol_frac_hdu = fits.PrimaryHDU(pol_frac)
#    pol_frac_hdu.header = I_header
#    pol_frac_hdu.writeto(f'polfrac/{channel:04d}-polfrac.fits', overwrite=True)

#    if not os.path.exists('polfrac_image'):
#        os.mkdir('polfrac_image')

#    print(f'Makin the polarization fraction image for channel {channel}')
#    mask=np.full(pol_frac.shape, np.nan)
#    mask[np.where((pol_frac >0.0)&(pol_frac <0.05))] =1

#    plt.figure()
##    plt.subplot(projection=I_wcs, slices=('x', 'y', 0, 0))
#    map1 = plt.imshow(pol_frac*mask, cmap='rainbow',vmin=0, vmax=0.05, interpolation='nearest')
#    plt.colorbar(map1, label = r'\%p')
#    plt.xlabel('RA')
#    plt.ylabel('DEC')
#    plt.xlim(3496, 4696)
#    plt.ylim(3496, 4696)
#    plt.title(f'Polarization fraction channel {channel}')
#    plt.savefig(f'polfrac_image/{channel:04d}_polfrac.pdf', dpi=400)
#    plt.close()


    Q = Q_data[0][x, y]
    U = U_data[0][x, y]
    I = I_data[0][x, y]
    linpol = lin_pol[x,y]
    pol_frac = np.sqrt(Q**2+U**2)/I
    chi = 0.5*np.arctan2(U_data[0][x, y], Q_data[0][x, y])

    return freq, wave2, I, Q, U, linpol, pol_frac, chi

def make_plots(channels, x, y):
    freq_list = []
    wave2_list = []
    I_list = []
    Q_list = []
    U_list = []
    linpol_list = []
    pol_frac_list = []
    chi_list = []

    teller = 0
    for i in range(channels):
        if teller in [20, 21, 23, 36,37,38,39,40,73,74,75]:
            teller +=1
            continue
        if teller == channels:
            break
        print(f'Retrieving data for channel {teller}')
        freq, wave2, I, Q, U, linpol, pol_frac, chi = make_lin_pol(teller, x, y)
#        print(np.sqrt(Q**2 + U**2)/I, pol_frac)
        freq_list.append(freq)
        wave2_list.append(wave2)
        I_list.append(I)
        Q_list.append(Q)
        U_list.append(U)
        linpol_list.append(linpol)
        pol_frac_list.append(pol_frac)
        chi_list.append(chi)
        teller +=1

    if not os.path.exists('Analysis_Images'):
        os.mkdir('Analysis_Images')

    if not os.path.exists(f'Analysis_Images/Pixel{x}_{y}'):
        os.mkdir(f'Analysis_Images/Pixel{x}_{y}')

    plt.figure()
#    plt.plot(wave2_list, I_list, 'k-.', alpha=0.5)
    plt.plot(wave2_list, I_list, 'k.', alpha=0.9)
    plt.ylabel('Flux [Jy]')
    plt.xlabel(r'$\lambda^2$ [m$^2$]]')
    plt.title(f'Stokes I [{x}, {y}]')
    plt.savefig(fname=f'Analysis_Images/Pixel{x}_{y}/flux_stokesI_{x}_{y}.pdf', dpi = 150)
    plt.close()

    plt.figure()
#    plt.plot(wave2_list, Q_list, 'k-.', alpha=0.5)
    plt.plot(wave2_list, Q_list, 'k.', alpha=0.9)
    plt.ylabel('Flux [Jy]')
    plt.xlabel(r'$\lambda^2$ [m$^2$]]')
    plt.title(f'Stokes Q [{x}, {y}]')
    plt.savefig(fname=f'Analysis_Images/Pixel{x}_{y}/flux_stokesQ_{x}_{y}.pdf', dpi = 150)
    plt.close()

    plt.figure()
#    plt.plot(wave2_list, U_list, 'k-.', alpha=0.5)
    plt.plot(wave2_list, U_list, 'k.', alpha=0.9)
    plt.ylabel('Flux [Jy]')
    plt.xlabel(r'$\lambda^2$ [m$^2$]]')
    plt.title(f'Stokes U [{x}, {y}]')
    plt.savefig(fname=f'Analysis_Images/Pixel{x}_{y}/flux_stokesU_{x}_{y}.pdf', dpi = 150)
    plt.close()

    plt.figure()
#    plt.plot(wave2_list, np.array(pol_frac_list)*100, 'k-.', alpha=0.5)
    plt.plot(wave2_list, np.array(pol_frac_list)*100, 'k.', alpha=0.9)
    plt.ylabel(r'p [$\%$]')
    plt.xlabel(r'$\lambda^2$ [m$^2$]]')
#    plt.ylim(0, 15)
    plt.title(f'Polarization degree [{x}, {y}]')
    plt.savefig(fname=f'Analysis_Images/Pixel{x}_{y}/pol_degree_{x}_{y}.pdf', dpi = 150)
    plt.close()

    plt.figure()
    plt.plot(freq_list, np.array(pol_frac_list)*100, 'r.', alpha=0.9)
#    plt.plot(wave2_list, np.array(pol_frac_list)*100, 'k.', alpha=0.9)
    plt.ylabel(r'p [$\%$]')
    plt.xlabel(r'$\lambda^2$ [m$^2$]]')
#    plt.ylim(0, 15)
    plt.title(f'Polarization degree [{x}, {y}]')
    plt.savefig(fname=f'Analysis_Images/Pixel{x}_{y}/pol_degree_freq_{x}_{y}.pdf', dpi = 150)
    plt.close()

    plt.figure()
#    plt.plot(wave2_list, chi_list, 'k-.', alpha=0.5)
    plt.plot(wave2_list, chi_list, 'k.', alpha=0.9)
    plt.ylabel(r'$\chi$ [deg]')
    plt.xlabel(r'$\lambda^2$ [m$^2$]]')
    plt.title(f'Angle of polarization [{x}, {y}]')
    plt.savefig(fname=f'Analysis_Images/Pixel{x}_{y}/pol_angle_{x}_{y}.pdf', dpi = 150)
    plt.close()

    plt.figure()
    plt.plot(np.array([freq_list])*1e-9, np.array([I_list])*1e3, 'r.', alpha=0.9, label='Stokes I')
    plt.plot(np.array([freq_list])*1e-9, np.array([linpol_list])*1e3, 'b.', alpha=0.9, label='Linear Polarized')
    plt.ylabel(r'Flux [mJy]')
    plt.xlabel(r'Frequency [GHz]')
    plt.title(f'Fluxes for pixel [{x}, {y}]')
    plt.yscale('log')
    plt.savefig(fname=f'Analysis_Images/Pixel{x}_{y}/fluxes.pdf', dpi = 150)
    plt.close()

# Think of the pixels in reverse
#make_plots(126, 4083, 3992)
#make_plots(126, 4063, 4007)
#make_plots(126, 4075, 4000)

def make_plot(xmin_cut, xmax_cut, ymin_cut, ymax_cut, pol_per):
    pol_pix = find_high_pol_pix(xmin_cut, xmax_cut, ymin_cut, ymax_cut, pol_per)
#    print(pol_pix)
    x = int(pol_pix[1])
    y = int(pol_pix[2])
    make_plots(126, x, y)

    return print(f'images of pixel [{x}, {y}] are done')

#make_plot(3920, 3970,4170, 4240 ,20)
make_plots(126, 3955, 4222)











