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
from matplotlib import gridspec

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

	#np.save(f'Correction/fitparamsregion{i}.npy', np.array([median_I_0_StokesI, median_a_StokesI]))

	return np.array([median_I_0_StokesI, median_a_StokesI])

def find_flux_region(channel, regionfile, region, rmsI,corr):

	import Integrate_flux as intf
	if corr==True:
		I_image = f'DATA/stokes_i_corr/{channel:04d}-I-image-pb.smoothed.fits'
	if corr==False:
		I_image = f'DATA/stokes_i/{channel:04d}-I-image-pb.smoothed.fits'

	fluxI, NbeamsI = intf.integratedflux(I_image, regionfile, region, hdul=None)
	uncertaintyI = intf.uncertainty_flux(I_image, fluxI, NbeamsI, region, regionfile_rms = None, rms = rmsI, delta_cal=0.0, hdul=None)

	pol_I_hdu = fits.open(f'DATA/stokes_i/{channel:04d}-I-image-pb.smoothed.fits')[0]
	I_header = pol_I_hdu.header

	freq = I_header['CRVAL3']
	wave2 = (c.c/freq)**2


	I = fluxI
	I_err = uncertaintyI

	return freq, wave2, I, I_err

def make_list_region(channels, regionfile, region, corr):
	freq_list = []
	wave2_list = []
	I_list = []
	I_list_err = []
	if corr ==True:
		rmsI = np.load('DATA/all_noiseI_corr.npy')[0]
	if corr ==False:
		rmsI = np.load('DATA/all_noiseI.npy')[0]

	teller = 0
	teller_rms = 0
	for i in range(channels):
		if teller in [20, 21, 23, 36,37,38,39,40,73,74,75]:
			teller +=1
			continue
		if teller == channels:
			break
		print(f'Retrieving data for channel {teller}')
		freq, wave2, I, I_err = find_flux_region(teller, regionfile, region, rmsI[teller_rms], corr)
	#print(np.sqrt(Q**2 + U**2)/I, pol_frac)
		freq_list.append(freq)
		wave2_list.append(wave2)
		I_list.append(I)
		I_list_err.append(I_err)
		teller +=1
		teller_rms +=1
	#np.save('Correction/freq_list.npy', freq_list)
	if corr ==True:
		np.save(f'Correction/Stokes_I_corr_region{region}.npy', I_list)
	if corr ==False:
		np.save(f'Correction/Stokes_I_region{region}.npy', I_list)
#	np.save('Correction/Uncorr_rms_I.npy', rmsI)
	return np.array(freq_list), np.array(wave2_list), np.array(I_list), np.array(I_list_err)

def run_mcmc_region(channels, regionfile, region, save_loc, plot = True, I_0_guess = 0.05, a_guess = 0.1):
	freq_list, wave2_list, I_list, I_err = make_list_region(channels, regionfile, region)

	save_loc = save_loc
	stokes_I_mcmc(np.array(I_list), np.array(I_err), np.array(freq_list), np.array(wave2_list), save_loc, I_0_guess = I_0_guess, a_guess = a_guess, plot = plot)

def find_correction():
	#run_mcmc_region(126, 'Regions/brightestsource.reg', 0, save_loc = f"Correction", plot = True, I_0_guess = 0.05, a_guess = 0.1)
	stokes_I = np.load(f'Correction/Stokes_I_uncorr.npy')
	freq = np.load(f'Correction/freq_list.npy')
	fit_params = np.load(f'Correction/correctionparams.npy')

	fit = fit_params[0]*(freq/1e9)**(fit_params[1])
#	print(fit)
#	print(stokes_I)
	correction = fit/stokes_I
	np.save(f'Correction/correction_factors.npy', correction)
#	print(correction)
	return correction

def apply_correction(channels):
	correction = np.load('Correction/correction_factors.npy')
	teller = 0
	teller_corr = 0

	if not os.path.exists('DATA/stokes_i_corr'):
		os.mkdir('DATA/stokes_i_corr')
	if not os.path.exists('DATA/stokes_q_corr'):
		os.mkdir('DATA/stokes_q_corr')
	if not os.path.exists('DATA/stokes_u_corr'):
		os.mkdir('DATA/stokes_u_corr')
		
	for i in range(channels):
		if teller in [20, 21, 23, 36,37,38,39,40,73,74,75]:
			teller +=1
			continue
		if teller == channels:
			break
		print(f"applying correction on channel {i}")
		Q_image_hdu = fits.open(f"DATA/stokes_q/{teller:04d}-Q-image-pb.smoothed.fits")[0]
		U_image_hdu = fits.open(f"DATA/stokes_u/{teller:04d}-U-image-pb.smoothed.fits")[0]
		I_image_hdu = fits.open(f"DATA/stokes_i/{teller:04d}-I-image-pb.smoothed.fits")[0]


		Q_header = Q_image_hdu.header
		Q_data = Q_image_hdu.data
		Q_hdu_corr = fits.PrimaryHDU()
		Q_hdu_corr.header = Q_header
		Q_hdu_corr.data = Q_data*correction[teller_corr]
		Q_hdu_corr.writeto(f"DATA/stokes_q_corr/{teller:04d}-Q-image-pb.smoothed.fits", overwrite=True)
		
		U_header = U_image_hdu.header
		U_data = U_image_hdu.data
		U_hdu_corr = fits.PrimaryHDU()
		U_hdu_corr.header = U_header
		U_hdu_corr.data = U_data*correction[teller_corr]
		U_hdu_corr.writeto(f"DATA/stokes_u_corr/{teller:04d}-U-image-pb.smoothed.fits", overwrite=True)
		
		I_header = I_image_hdu.header
		I_data = I_image_hdu.data
		I_hdu_corr = fits.PrimaryHDU()
		I_hdu_corr.header = I_header
		I_hdu_corr.data = I_data*correction[teller_corr]
		I_hdu_corr.writeto(f"DATA/stokes_i_corr/{teller:04d}-I-image-pb.smoothed.fits", overwrite=True)

		teller_corr +=1
		teller += 1
	return
#find_correction()
#apply_correction(126)

def plot_corr(I, I_err, I_corr, I_err_corr, freq, wave2,fitparams, fitparams_corr, region,regionfile='Correction/5brightsources.reg'):

	def freq_to_lambda(freq):
		wave_sq = (3.e8/freq)**2
		return(wave_sq)

	freqs = np.array([1.0,1.2,1.5])*1.e9
	freqticks = np.array([freq_to_lambda(freq) for freq in freqs])
	stokes_I = np.array(I)
	stokes_I_err = np.array(I_err)
	stokes_I_corr = np.array(I_corr)
	stokes_I_err_corr = np.array(I_err_corr)
	freq = freq
	wave2 = wave2

	freq_mcmc = np.linspace(0.89, 1.81, 150)
	wave_mcmc = (c.c/(freq_mcmc*1e9))**2
	freq = np.array(freq)/1e9
	print(wave2.shape, stokes_I.shape)
	
	fit_params = fitparams
	fit_params_corr = fitparams_corr
	
	fig = plt.figure(figsize=(10,14))
	gs = gridspec.GridSpec(4,1,height_ratios = [4,1,4,1], hspace=0)
	ax1 = fig.add_subplot(gs[0])
	ax2 = fig.add_subplot(gs[2], sharex = ax1)
	ax3 = ax1.twiny()

	ax1.grid(color='lightgrey', linestyle=':', linewidth=0.7)
	ax2.grid(color='lightgrey', linestyle=':', linewidth=0.7)

	ax1.set_ylabel(r'$I$ [mJy beam$^{-1}$]', fontsize=28)
	ax2.set_ylabel(r'$I_{corr}$ [mJy beam$^{-1}$]', fontsize=28) 
	ax3.set_xlabel(r'$\nu$ [GHz]', fontsize=28)

	ax1.tick_params(axis='y',labelsize=20)
	ax2.tick_params(axis='y',labelsize=20)
	ax3.tick_params(axis='x',labelsize=20) 

	ax1.yaxis.set_ticks_position('both')
	ax2.yaxis.set_ticks_position('both')
	ax2.xaxis.set_ticks_position('both')

	ax1.set_xlim([0.0275, 0.105])
	ax3.set_xlim([0.0275, 0.105])
##	ax1.set_ylim([-0.5*np.pi-0.5, 0.5*np.pi+0.5])
#	ax2.set_ylim([-1, 100*median_p0 +1])
	#ax2.set_ylim([-0.05,0.65])
	#ax2.set_yticks([0.,0.1,0.2,0.3,0.4,0.5,0.6])

	ax3.set_xticks(freqticks)
	ax3.set_xticklabels('{:.1f}'.format(freq) for freq in freqs/1.e9)

	xticklabels1 = ax1.get_xticklabels()
	xticklabels2 = ax2.get_xticklabels()
	plt.setp(xticklabels1, visible=False)
	plt.setp(xticklabels2, visible=False)
	#plt.suptitle(sourcetitles[row], fontsize=18)



	ax1.errorbar(wave2, stokes_I*1e3, yerr=stokes_I_err*1e3, color = 'black', capsize= 3, capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5)
#	ax1.fill_between(lambda_mcmc,pol_ang_lower, pol_ang_upper, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
	ax2.errorbar(wave2, stokes_I_corr*1e3, yerr=stokes_I_err_corr*1e3, color = 'black', capsize= 3, capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5)
#	ax2.fill_between(lambda_mcmc,100*pol_frac_lower, 100*pol_frac_upper, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)

	ax1.plot(wave_mcmc, 1e3*fit_params[0]*(freq_mcmc)**fit_params[1], '-', color = palette(0),linewidth=2., zorder=len(wave2)+1)
	ax2.plot(wave_mcmc, 1e3*fit_params_corr[0]*(freq_mcmc)**fit_params_corr[1], '-', color = palette(0),linewidth=2., zorder=len(wave2)+1)

	# print res
	ax1res = fig.add_subplot(gs[1],sharex = ax1)
	ax2res = fig.add_subplot(gs[3],sharex = ax1)

	ax1res.set_xlim([0.0275, 0.105])
	ax2res.set_xlim([0.0275, 0.105])
	#ax1res.set_ylim([-1, 1])
	#ax2res.set_ylim([-1.5, 1.5])

	ax1res.grid(color='lightgrey', linestyle=':', linewidth=0.7)
	ax2res.grid(color='lightgrey', linestyle=':', linewidth=0.7)

	ax1res.set_ylabel('res', fontsize=24)
	ax2res.set_ylabel('res', fontsize=24)
	ax2res.set_xlabel(r'$\lambda^2$ [m$^2$]', fontsize=28)

	ax1res.tick_params(axis='y',labelsize=16)
	ax2res.tick_params(axis='y',labelsize=16)
	ax2res.tick_params(axis='x',labelsize=20)

	ax1res.xaxis.set_ticks_position('both')
	ax2res.xaxis.set_ticks_position('both')
	ax1res.yaxis.set_ticks_position('both')
	ax2res.yaxis.set_ticks_position('both')

	xticklabels1 = ax1res.get_xticklabels()
	xticklabels2 = ax2res.get_xticklabels()
	plt.setp(xticklabels1, visible=False)


	ax1res.errorbar(wave2, (1e3*stokes_I - 1e3*fit_params[0]*(freq)**fit_params[1]), yerr=1e3*stokes_I_err,fmt='.',markersize=4,color='k',alpha=0.75)
	ax1res.plot(wave_mcmc,[0]*len(wave_mcmc),'--',color=palette(0),linewidth=1., zorder=len(wave2)+1, alpha=0.8)
	ax2res.errorbar(wave2, (1e3*stokes_I_corr - 1e3*fit_params[0]*(freq)**fit_params[1]), yerr=1e3*stokes_I_err,fmt='.',markersize=4,color='k',alpha=0.75)
	ax2res.plot(wave_mcmc,[0]*len(wave_mcmc),'--',color=palette(0),linewidth=1., zorder=len(wave2)+1, alpha=0.8)

	plt.tight_layout()
	plt.savefig(f'Correction/manual_bandpass_region{i}.pdf', dpi = 150)
	plt.close()
	

for i in range(5):
	freq, wave2, I_list, I_err = make_list_region(126, regionfile='Correction/5brightsources.reg', region=i, corr=False)
	_, _, I_list_corr, I_err_corr = make_list_region(126, regionfile='Correction/5brightsources.reg', region=i, corr=True)
	fitparams = stokes_I_mcmc(I_list, I_err, freq, wave2, save_loc='Correction', I_0_guess = 0.05, a_guess = 0.1, plot = False)
	fitparams_corr = stokes_I_mcmc(I_list_corr, I_err_corr, freq, wave2, save_loc='Correction', I_0_guess = 0.05, a_guess = 0.1, plot = False)
	plot_corr(I_list, I_err, I_list_corr, I_err_corr, freq, wave2,fitparams,fitparams_corr, region=i, regionfile='Correction/5brightsources.reg')








