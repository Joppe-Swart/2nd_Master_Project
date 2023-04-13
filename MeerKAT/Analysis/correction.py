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

def stokes_I_mcmc(flux, flux_err, freq, wave2, save_loc, I_0_guess = 0.05, a_guess = 0.1, plot = True):

	freq_mcmc = np.array(freq)/1e9
	freq = np.array(freq)/1e9


	freq = freq # Lijst met frequenties in GHz
	y_obs = flux # data punten 
	dy = flux_err # error in de data
#	print(flux_err)

	print(f'defining the parameters for mcmc')
	def lnL_StokesI(theta, freq, y, yerr):
		I_0, a = theta
		model_stokes_I = I_0*freq**(a )#+ b*np.log10(freq))#*u.uJy
		inv_sigma2 = 1.0/(np.power(yerr,2))

		return -0.5*(np.sum((y-model_stokes_I)**2*inv_sigma2))

	def lnprior_StokesI(theta):
		I_0, a = theta
		if 0 < I_0 < 10 and -3 < a < 3:#and -5 < b < 5:
			return 0.0
		return -np.inf

	def lnprob_StokesI(theta, freq, y, yerr):
		lp_StokesI = lnprior_StokesI(theta)
		if not np.isfinite(lp_StokesI):
			return -np.inf
		return lp_StokesI + lnL_StokesI(theta, freq, y, yerr)

	print(f'Start the actual mcmc')
	ndim_StokesI, nwalkers = 2, 150
	theta_StokesI_guess = np.array([I_0_guess, a_guess])
	pos_StokesI = [theta_StokesI_guess + 1e-4*np.random.randn(ndim_StokesI) \
		   for i in range(nwalkers)]
	sampler_StokesI = emcee.EnsembleSampler(nwalkers, ndim_StokesI, 
				lnprob_StokesI, args=(freq, y_obs, dy))
	tmp = sampler_StokesI.run_mcmc(pos_StokesI, 650)

	if plot == True:
		if not os.path.exists(f'{save_loc}'):
			os.makedirs(f'{save_loc}')

		fig, axes = plt.subplots(ncols=1, nrows=2)
		fig.set_size_inches(12,12)
		axes[0].plot(sampler_StokesI.chain[:, :, 0].transpose(), color='black', alpha=0.3)
		axes[0].set_ylabel(r'$I_0$')
		axes[0].axvline(150, ls='dashed', color='red')
		axes[1].plot(sampler_StokesI.chain[:, :, 1].transpose(), color='black', alpha=0.3)
		axes[1].set_ylabel(r'$a$')
		axes[1].axvline(150, ls='dashed', color='red')
		fig.savefig(f'{save_loc}/chain_StokesI.pdf')
		plt.close()

	samples_StokesI = sampler_StokesI.chain[:, 150:, :].reshape((-1, 2))

	if plot == True:
		fig = corner.corner(samples_StokesI, labels=[r"$I_0$", r"$a$", r"b"], \
		quantiles=[0.16, 0.50, 0.84], show_titles=True)#truths=[S_pl_guess, a_pl_guess],)
		fig.savefig(f'{save_loc}/corner_StokesI.pdf')
		plt.close()

	median_I_0_StokesI = np.percentile(samples_StokesI[:, 0], 50.0)
	median_a_StokesI = np.percentile(samples_StokesI[:, 1], 50.0)

	p16_I_0_StokesI = np.percentile(samples_StokesI[:, 0], 16)
	p16_a_StokesI = np.percentile(samples_StokesI[:, 1], 16)

	p84_I_0_StokesI = np.percentile(samples_StokesI[:, 0], 84)
	p84_a_StokesI = np.percentile(samples_StokesI[:, 1], 84)

	sigma_I_0_StokesI = 0.5*(p84_I_0_StokesI-p16_I_0_StokesI)
	sigma_a_StokesI = 0.5*(p84_a_StokesI-p16_a_StokesI)


	MCMC_StokesI=np.empty(shape=[len(samples_StokesI[:,0]), 0])
	for i in range(len(freq_mcmc)):
		MCMC_StokesI = np.append(MCMC_StokesI, (samples_StokesI[:,0]*freq_mcmc[i]**\
		(samples_StokesI[:,1])).reshape(len(samples_StokesI[:,0]), 1), axis = 1)# + 

	pl_16 = []
	pl_84 = []
	for i in range(len(freq_mcmc)):
		pl_16.append(np.percentile(np.sort(MCMC_StokesI[:,i]),16))
		pl_84.append(np.percentile(np.sort(MCMC_StokesI[:,i]),84))

	StokesI_16 = np.array(pl_16)#*u.mJy
	StokesI_84 = np.array(pl_84)#*u.mJy


	if plot == True:
		plt.figure()
		plt.errorbar(freq, flux,yerr= flux_err, color = 'black', capsize= 3, 
		 capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5, label ='Observations')
		# plt.plot(freq_m/1E6, I_table[:,0], 'k-.', alpha=0.5)
		plt.plot(freq_mcmc, median_I_0_StokesI*(freq_mcmc)**(median_a_StokesI), \
			 color = palette(0), alpha = 1, label='MCMC fit')
		plt.fill_between(freq_mcmc,StokesI_16, StokesI_84, facecolor = palette(1), \
				 edgecolor = palette(1), alpha = 0.3)
		plt.ylabel('Flux [Jy]')
		plt.xlabel('Frequency [GHz]')
		plt.title('Stokes I')
		plt.yscale('log')
		plt.xscale('log')
		plt.xticks([0.9,1.2,1.5],[0.9,1.2,1.5])
		plt.legend(loc='best', frameon=True, shadow= True)
		plt.savefig(f'{save_loc}/flux_stokesI_MCMC.pdf')
		plt.close()

	np.save('correctionparams.npy', np.array([median_I_0_StokesI, median_a_StokesI]))

	return np.array([median_I_0_StokesI, median_a_StokesI])

def find_flux_region(channel, regionfile, region, rmsI):

	import Integrate_flux as intf

	I_image = f'stokes_i/{channel:04d}-I-image-pb.smoothed.fits'

	fluxI, NbeamsI = intf.integratedflux(I_image, regionfile, region, hdul=None)
	uncertaintyI = intf.uncertainty_flux(I_image, fluxI, NbeamsI, region, regionfile_rms = None, rms = rmsI, delta_cal=0.0, hdul=None)

	pol_I_hdu = fits.open(f'stokes_i/{channel:04d}-I-image-pb.smoothed.fits')[0]
	I_header = pol_I_hdu.header

	freq = I_header['CRVAL3']
	wave2 = (c.c/freq)**2


	I = fluxI
	I_err = uncertaintyI

	return freq, wave2, I, I_err

def make_list_region(channels, regionfile, region):
	freq_list = []
	wave2_list = []
	I_list = []
	I_list_err = []

	rmsI = np.load('all_noiseI.npy')[0]

	teller = 0
	teller_rms = 0
	for i in range(channels):
		if teller in [20, 21, 23, 36,37,38,39,40,73,74,75]:
			teller +=1
			continue
		if teller == channels:
			break
		print(f'Retrieving data for channel {teller}')
		freq, wave2, I, I_err = find_flux_region(teller, regionfile, region, rmsI[teller_rms])
	#print(np.sqrt(Q**2 + U**2)/I, pol_frac)
		freq_list.append(freq)
		wave2_list.append(wave2)
		I_list.append(I)
		I_list_err.append(I_err)
		teller +=1
		teller_rms +=1
	np.save('freq_list.npy', freq_list)
	np.save('brightest_stokes_I.npy', I_list)
	return freq_list, wave2_list, I_list, I_list_err

def run_mcmc_region(channels, regionfile, region, save_loc, plot = True, I_0_guess = 0.05, a_guess = 0.1):
	freq_list, wave2_list, I_list, I_err = make_list_region(channels, regionfile, region)

	save_loc = save_loc
	stokes_I_mcmc(np.array(I_list), np.array(I_err), np.array(freq_list), np.array(wave2_list), save_loc, I_0_guess = I_0_guess, a_guess = a_guess, plot = plot)

def find_correction():
	run_mcmc_region(126, 'brightestsource1.reg', 0, save_loc = f"Analysis_Images/correction", plot = True, I_0_guess = 0.05, a_guess = 0.1)
	stokes_I = np.load('brightest_stokes_I.npy')
	freq = np.load('freq_list.npy')
	fit_params = np.load('correctionparams.npy')

	fit = fit_params[0]*(freq/1e9)**(fit_params[1])
#	print(fit)
#	print(stokes_I)
	correction = fit/stokes_I
	np.save('correction_factors.npy', correction)
#	print(correction)
	return correction

def apply_correction(channels):
	correction = np.load('correction_factors.npy')
	teller = 0
	teller_corr = 0

	if not os.path.exists('stokes_i_corr'):
		os.mkdir('stokes_i_corr')
	if not os.path.exists('stokes_q_corr'):
		os.mkdir('stokes_q_corr')
	if not os.path.exists('stokes_u_corr'):
		os.mkdir('stokes_u_corr')
		
	for i in range(channels):
		if teller in [20, 21, 23, 36,37,38,39,40,73,74,75]:
			teller +=1
			continue
		if teller == channels:
			break
		print(f"applying correction on channel {i}")
		Q_image_hdu = fits.open(f"stokes_q/{teller:04d}-Q-image-pb.smoothed.fits")[0]
		U_image_hdu = fits.open(f"stokes_u/{teller:04d}-U-image-pb.smoothed.fits")[0]
		I_image_hdu = fits.open(f"stokes_i/{teller:04d}-I-image-pb.smoothed.fits")[0]


		Q_header = Q_image_hdu.header
		Q_data = Q_image_hdu.data
		Q_hdu_corr = fits.PrimaryHDU()
		Q_hdu_corr.header = Q_header
		Q_hdu_corr.data = Q_data*correction[teller_corr]
		Q_hdu_corr.writeto(f"stokes_q_corr/{teller:04d}-Q-image-pb.smoothed.fits", overwrite=True)
		
		U_header = U_image_hdu.header
		U_data = U_image_hdu.data
		U_hdu_corr = fits.PrimaryHDU()
		U_hdu_corr.header = U_header
		U_hdu_corr.data = U_data*correction[teller_corr]
		U_hdu_corr.writeto(f"stokes_u_corr/{teller:04d}-U-image-pb.smoothed.fits", overwrite=True)
		
		I_header = I_image_hdu.header
		I_data = I_image_hdu.data
		I_hdu_corr = fits.PrimaryHDU()
		I_hdu_corr.header = I_header
		I_hdu_corr.data = I_data*correction[teller_corr]
		I_hdu_corr.writeto(f"stokes_i_corr/{teller:04d}-I-image-pb.smoothed.fits", overwrite=True)

		teller_corr +=1
		teller += 1
	return
apply_correction(126)











