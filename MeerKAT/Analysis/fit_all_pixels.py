"""

Author: Joppe Swart
Created: March 2023
Last modified: March 2023
Description: This script fits all polarized pixels

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


def stokes_qu_mcmc(flux, flux_err, fluxU, fluxU_err, fluxQ, fluxQ_err, freq, \
		   wave2, save_loc, mcmc_I = True, mcmc_QU = True, I_0_guess = 0.05, a_guess = 0.1, \
		   p0_guess = 0.3, chi0_guess = 0.2, rm_guess = -25,sigma_rm_geuss = 10, plot = True):
	if mcmc_I == True:
	
		pol_ang = 0.5*np.arctan2(fluxU, fluxQ) # One of the polar coordinates
		linear_pol = np.sqrt(fluxQ**2+fluxU**2)
		deg_pol = linear_pol/flux

#		print(fluxQ.shape, fluxQ_err.shape, linear_pol.shape)
		pol_ang_err = np.sqrt(((((fluxU**4)/(fluxQ**6))*fluxQ_err**2)+(((fluxU**2)/(fluxQ**4))*fluxU_err**2))/((1+(fluxU/fluxQ)**2)**2))#fluxQ*fluxU_err-fluxU*fluxQ_err)/(linear_pol**2)
		linear_pol_err = np.sqrt(fluxU**2*fluxU_err**2+fluxQ**2*fluxQ_err**2)/linear_pol
		deg_pol_err = np.sqrt((linear_pol_err/flux)**2 + (linear_pol**2/flux**4)*flux_err**2)

		freq_mcmc = np.array(freq)/1e9
		freq = np.array(freq)/1e9


		freq = freq # Lijst met frequenties in GHz
		y_obs = flux # data punten 
		dy = flux_err # error in de data
#		print(flux_err)

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
			(samples_StokesI[:,1])).reshape(len(samples_StokesI[:,0]), 1), axis = 1)

		pl_16 = []
		pl_84 = []
		for i in range(len(freq_mcmc)):
			pl_16.append(np.percentile(np.sort(MCMC_StokesI[:,i]),16))
			pl_84.append(np.percentile(np.sort(MCMC_StokesI[:,i]),84))

		StokesI_16 = np.array(pl_16)#*u.mJy
		StokesI_84 = np.array(pl_84)#*u.mJy

		def red_chi_StokesI(I_0, a):
			chi_StokesI = 0
			model = I_0*freq**(a)
			for i in range(len(freq)):
				chi_StokesI += ((y_obs[i]-model[i])/dy[i])**2
				red_chi_StokesI = chi_StokesI/(len(freq)-2)
			return red_chi_StokesI

		if plot == True:
			axis = [0.9, 1.2, 1.5]
			plt.figure()
			plt.errorbar(freq, flux,yerr= flux_err, color = 'black', capsize= 3, 
			 capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5, label ='Observations')
			plt.plot(freq_mcmc, median_I_0_StokesI*(freq_mcmc)**(median_a_StokesI), \
				 color = palette(0), alpha = 1, label='MCMC fit')
			plt.fill_between(freq_mcmc,StokesI_16, StokesI_84, facecolor = palette(1), \
					 edgecolor = palette(1), alpha = 0.3)
			plt.ylabel('Flux [Jy]')
			plt.xlabel('Frequency [GHz]')
			plt.xticks(axis, axis)
			plt.title('Stokes I')
			plt.yscale('log')
			plt.xscale('log')
			plt.legend(loc='best', frameon=True, shadow= True)
			plt.savefig(f'{save_loc}/flux_stokesI_MCMC.pdf')
			plt.close()

#		np.save('stokes_i_params.npy', np.array([median_I_0_StokesI, median_a_StokesI]))



	if mcmc_QU == True:
		Stokes_I = median_I_0_StokesI*freq**(median_a_StokesI)# + median_b_StokesI*np.log10(freq))
		Stokes_I_mcmc = median_I_0_StokesI*freq**(median_a_StokesI)# + median_b_StokesI*np.log10(freq_mcmc))

		# start with QU
		wave2 = wave2

		lambda_mcmc = wave2
		x_lambda = wave2 # Lijst met frequenties in GHz
		y_U = fluxU # data punten 
		dy_U = fluxU_err# error in de data
		y_Q = fluxQ # data punten 
		dy_Q = fluxQ_err # error in de data


		def lnLQU(thetaQU, x_lambda, y1, yerr1, y2, yerr2):
			p0, chi0, rm, sigma_rm = thetaQU
			model1 = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x_lambda**2)) * \
				 np.sin(2*(chi0 + rm*x_lambda)) #Stokes U
			model2 = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x_lambda**2)) * \
				 np.cos(2*(chi0 + rm*x_lambda)) #Stokes Q

			return -0.5*(np.sum(((y1-model1)/yerr1)**2 + ((y2-model2)/yerr2)**2))

		def lnpriorQU(thetaQU):
			p0, chi0, rm, sigma_rm = thetaQU
			if 0.01 < p0 < 1 and -100 < rm < 100 and  sigma_rm > 0 and 0 < chi0 < np.pi:
				return 0.0
			return -np.inf

		def lnprobQU(thetaQU, x_lambda, y1, yerr1, y2, yerr2):
			lp = lnpriorQU(thetaQU)
			if not np.isfinite(lp):
				return -np.inf
			return lp + lnLQU(thetaQU, x_lambda, y1, yerr1, y2, yerr2)

		ndim, nwalkers = 4, 500
		theta_guess = np.array([p0_guess, chi0_guess, rm_guess, sigma_rm_geuss])
		pos = [theta_guess + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
		samplerQU = emcee.EnsembleSampler(nwalkers, ndim, lnprobQU, args=(x_lambda, y_U, dy_U, y_Q, dy_Q))
		tmpQU = samplerQU.run_mcmc(pos, 700)

		if plot == True:
			if not os.path.exists(f'{save_loc}'):
				os.makedirs(f'{save_loc}')

			fig, axes = plt.subplots(ncols=1, nrows=4)
			fig.set_size_inches(12,12)
			axes[0].plot(samplerQU.chain[:, :, 0].transpose(), color='black', alpha=0.3)
			axes[0].set_ylabel(r'$p_0$')
			axes[0].axvline(200, ls='dashed', color='red')
			axes[1].plot(samplerQU.chain[:, :, 1].transpose(), color='black', alpha=0.3)
			axes[1].set_ylabel(r'$\chi_0$')
			axes[1].axvline(200, ls='dashed', color='red')
			axes[2].plot(samplerQU.chain[:, :, 2].transpose(), color='black', alpha=0.3)
			axes[2].set_ylabel(r'$RM$')
			axes[2].axvline(200, ls='dashed', color='red')
			axes[3].plot(samplerQU.chain[:, :, 3].transpose(), color='black', alpha=0.3)
			axes[3].set_ylabel(r'$\sigma_{RM}$')
			axes[3].axvline(200, ls='dashed', color='red')
			fig.savefig(f'{save_loc}/chain_StokesQU.pdf')
			plt.close()

		samplesQU = samplerQU.chain[:, 400:, :].reshape((-1, 4))


		if plot == True:
			fig = corner.corner(samplesQU, labels=[r"$p_0$", r"$\chi_0$", r"RM", r"$\sigma_{RM}$"], 
			quantiles=[0.16, 0.50, 0.84], show_titles=True)
			fig.savefig(f'{save_loc}/corner_StokesQU.pdf')
			plt.close()

		median_p0 = np.percentile(samplesQU[:, 0], 50.0)
		median_chi0 = np.percentile(samplesQU[:, 1], 50.0)
		median_RM = np.percentile(samplesQU[:, 2], 50.0)
		median_sigmaRM = np.percentile(samplesQU[:, 3], 50.0)

		p16_p0 = np.percentile(samplesQU[:, 0], 16)
		p16_chi0 = np.percentile(samplesQU[:, 1], 16)
		p16_RM = np.percentile(samplesQU[:, 2], 16)
		p16_sigmaRM = np.percentile(samplesQU[:, 3], 16)

		p84_p0 = np.percentile(samplesQU[:, 0], 84)
		p84_chi0 = np.percentile(samplesQU[:, 1], 84)
		p84_RM = np.percentile(samplesQU[:, 2], 84)
		p84_sigmaRM = np.percentile(samplesQU[:, 3], 84)

		sigma_p0 = 0.5*(p84_p0-p16_p0)
		sigma_chi0 = 0.5*(p84_chi0-p16_chi0)
		sigma_RM = 0.5*(p84_RM-p16_RM)
		sigma_sigmaRM = 0.5*(p84_sigmaRM-p16_sigmaRM)

		MCMC_U=np.empty(shape=[len(samplesQU[:,0]), 0])
		MCMC_Q=np.empty(shape=[len(samplesQU[:,0]), 0])
		MCMC_pol_ang=np.empty(shape=[len(samplesQU[:,0]), 0])
		for i in range(len(lambda_mcmc)):
			MCMC_U = np.append(MCMC_U, (Stokes_I_mcmc[i] * samplesQU[:,0] * np.exp(-2*(samplesQU[:,3])*(lambda_mcmc[i]**2)) * 
						np.sin(2*(samplesQU[:,1] + samplesQU[:,2]*lambda_mcmc[i]))).reshape(len(samplesQU[:,0]), 1), 
						axis = 1)
			MCMC_Q = np.append(MCMC_Q, (Stokes_I_mcmc[i] * samplesQU[:,0] * np.exp(-2*(samplesQU[:,3])*(lambda_mcmc[i]**2)) * 
						np.cos(2*(samplesQU[:,1] + samplesQU[:,2]*lambda_mcmc[i]))).reshape(len(samplesQU[:,0]), 1), 
						axis = 1)
			MCMC_pol_ang = np.append(MCMC_pol_ang, (0.5*np.arctan2(Stokes_I_mcmc[i] * samplesQU[:,0] * np.exp(-2*(samplesQU[:,3])*
						(lambda_mcmc[i]**2)) * np.sin(2*(samplesQU[:,1] + samplesQU[:,2]*lambda_mcmc[i])),
						Stokes_I_mcmc[i] * samplesQU[:,0] * np.exp(-2*(samplesQU[:,3])*(lambda_mcmc[i]**2)) * 
						np.cos(2*(samplesQU[:,1] + samplesQU[:,2]*lambda_mcmc[i])))).reshape(len(samplesQU[:,0]), 1), axis = 1)

		U_16 = []
		U_84 = []
		Q_16 = []
		Q_84 = []
		pol_ang_16 = []
		pol_ang_84 = []

		for i in range(len(lambda_mcmc)):
			U_16.append(np.percentile(np.sort(MCMC_U[:,i]),16))
			U_84.append(np.percentile(np.sort(MCMC_U[:,i]),84))
			Q_16.append(np.percentile(np.sort(MCMC_Q[:,i]),16))
			Q_84.append(np.percentile(np.sort(MCMC_Q[:,i]),84))
			pol_ang_16.append(np.percentile(np.sort(MCMC_pol_ang[:,i]),16))
			pol_ang_84.append(np.percentile(np.sort(MCMC_pol_ang[:,i]),84))

		StokesQ_16 = np.array(Q_16)#*u.mJy
		StokesQ_84 = np.array(Q_84)#*u.mJy
		StokesU_16 = np.array(U_16)#*u.mJy
		StokesU_84 = np.array(U_84)#*u.mJy
		pol_ang_16_StoeksUQ = np.array(pol_ang_16)#*u.mJy
		pol_ang_84_StokesUQ = np.array(pol_ang_84)#*u.mJy


		def red_chi_U(p0, chi0, RM, sigma_RM):
			chi = 0
			model = Stokes_I * p0 * np.exp(-2*(sigma_RM)*(x_lambda**2)) * np.sin(2*(chi0 + RM*x_lambda))
			for i in range(len(x_lambda)):
				chi += ((y_U[i]-model[i])/dy_U[i])**2
				red_chi = chi/(len(x_lambda)-4)
			return red_chi

		def red_chi_Q(p0, chi0, RM, sigma_RM):
			chi = 0
			model = Stokes_I * p0 * np.exp(-2*(sigma_RM)*(x_lambda**2)) * np.cos(2*(chi0 + RM*x_lambda))
			for i in range(len(x_lambda)):
				chi += ((y_Q[i]-model[i])/dy_Q[i])**2
				red_chi = chi/(len(x_lambda)-4)
			return red_chi

		def Q(p0, chi0, rm , sigma_rm, Stokes_I, x):
			Q = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x**2)) * np.cos(2*(chi0 + rm*x))
			return Q

		def U(p0, chi0, rm , sigma_rm, Stokes_I, x):
			U = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x**2)) * np.sin(2*(chi0 + rm*x))
			return U

		if plot == True:
			axis = [0.9, 1.2, 1.5]
			plt.figure()
			plt.errorbar(wave2, fluxU,yerr= fluxU_err, color = 'black', capsize= 3, capthick=0.5, 
			 fmt='.', markersize= 7, elinewidth = 0.5, label ='Observations')
			plt.plot(lambda_mcmc, Stokes_I*median_p0*np.exp(-2*median_sigmaRM*lambda_mcmc**2)*np.sin(2*(median_chi0+
			 median_RM*lambda_mcmc)), color = palette(0), alpha = 1, label='MCMC fit')
			plt.fill_between(lambda_mcmc,StokesU_16, StokesU_84, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
			plt.ylabel('Flux [Jy]')
			plt.xlabel(r' $\lambda^2 \, [m^2]$')
			plt.title('Stokes U')
			plt.legend(loc='best', frameon=True, shadow= True)
			plt.savefig(f'{save_loc}/flux_stokesU_MCMC.pdf', dpi = 150)
			plt.close()

			plt.figure()
			plt.errorbar(wave2, fluxQ,yerr= fluxQ_err, color = 'black', capsize= 3, capthick=0.5, 
			 fmt='.', markersize= 7, elinewidth = 0.5, label ='Observations')
			plt.plot(lambda_mcmc, Stokes_I*median_p0*np.exp(-2*median_sigmaRM*lambda_mcmc**2)*np.cos(2*(median_chi0+
			 median_RM*lambda_mcmc)), color = palette(0), alpha = 1, label='MCMC fit')
			plt.fill_between(lambda_mcmc,StokesQ_16, StokesQ_84, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
			plt.ylabel('Flux [Jy]')
			plt.xlabel(r' $\lambda^2 \, [m^2]$')
			plt.title('Stokes Q')
			plt.legend(loc='best', frameon=True, shadow= True)
			plt.savefig(f'{save_loc}/flux_stokesQ_MCMC.pdf', dpi = 150)
			plt.close()

			plt.figure()
			plt.errorbar(wave2, pol_ang ,yerr= np.abs(pol_ang_err), color = 'black', capsize= 3, capthick=0.5, \
			 fmt='.', markersize= 7, elinewidth = 0.5, label ='Observations')
			plt.plot(lambda_mcmc, 0.5*np.arctan2(U(median_p0, median_chi0, median_RM , median_sigmaRM, Stokes_I, lambda_mcmc),\
			 Q(median_p0, median_chi0, median_RM , median_sigmaRM, Stokes_I, lambda_mcmc)),\
			 color = palette(0), alpha = 1, label='MCMC fit')
			plt.fill_between(lambda_mcmc,pol_ang_16, pol_ang_84, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
			plt.ylabel(r'Polarization Angle $(\chi)$ [rad]')
			plt.ylim(-np.pi, np.pi)
			plt.xlabel(r' $\lambda^2 \, [m^2]$')
			plt.title('Polarization angle from individual mcmc procedure')
			plt.legend(loc='best', frameon=True, shadow= True)
			plt.savefig(f'{save_loc}/pol_ang_MCMC.pdf', dpi = 150)
			plt.close()


			plt.figure()
			plt.errorbar(freq, 100*deg_pol,yerr= 100*deg_pol_err, color = 'black', \
					capsize= 3, capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5)
			plt.plot(freq_mcmc, 100*np.sqrt(U(median_p0, median_chi0, median_RM , median_sigmaRM, Stokes_I, lambda_mcmc)**2 +
			Q(median_p0, median_chi0, median_RM , median_sigmaRM, Stokes_I, lambda_mcmc)**2)/
			Stokes_I_mcmc , color = palette(0), alpha = 1, label='MCMC fit')
			plt.ylabel('p [%]')
			plt.xlabel('Frequency [GHz]')
			plt.title('Degree of polarization')
			plt.legend(loc = 'best', frameon=True, shadow= True)
			plt.savefig(f'{save_loc}/degree_of_pol.pdf', dpi = 150)
			plt.close()



def fit_pixels(channels):
	"""
	fit_pixels: funciton that fits all the pixels in a given refgion
	"""
	
	#first find which pixels we want to fit
	pol_hdu = fits.open('polarization_fraction_map_masked.fits')[0]
#	print(pol_hdu.data.shape)
	mask = np.full(pol_hdu.data.shape, np.nan)
#	print(mask.shape)
	mask[0,0,3950:3955,4220:4225] = 1
	data = pol_hdu.data*mask
#	print(data.shape)

	Q_cube = fits.open('Q_cube.fits')
	U_cube = fits.open('U_cube.fits')
	I_cube = fits.open('I_cube.fits')


	j = np.where(data[~np.isnan(data)])
	print(np.array(j).size)
#	print(data[~np.isnan(data)][0])
#	print(f'j = {j}')
	fluxQ = []
	fluxU = []
	fluxI = []
	freq = []

	for i in range(channels):
		freq.append(I_cube[i].header['CRVAL3'])
		
		Q_data = Q_cube[i].data*mask
		Q = Q_data[~np.isnan(Q_data)]
		fluxQ.append(Q[j])
		
		U_data = U_cube[i].data*mask
		U = U_data[~np.isnan(U_data)]
		fluxU.append(U[j])
		
		I_data = I_cube[i].data*mask
		I = I_data[~np.isnan(I_data)]
		fluxI.append(I[j])
	
	Q_cube.close()
	U_cube.close()
	I_cube.close()

	fluxI = np.array(fluxI)
	fluxU = np.array(fluxU)
	fluxQ = np.array(fluxQ)

	with fits.open('rmsynth_phi.fits') as rm:
		rm_data = rm[0].data*mask
		RM = rm_data[~np.isnan(rm_data)]
	print(f'the rm is {RM[j]}')
	freq = np.array(freq)
	wave2 = (c.c/freq)**2

	I_err = np.load('all_noiseI_corr.npy')[0]
	Q_err = np.load('all_noiseQ_corr.npy')[0]
	U_err = np.load('all_noiseU_corr.npy')[0]
	
	best_params = np.ones([])
	for i in range(np.array(j).size):
		save_loc = f'fit_pixel/test/pixel{i}'
		print(i)
		stokes_I = fluxI[:,i]
		stokes_Q = fluxQ[:,i]
		stokes_U = fluxU[:,i]
#		print(stokes_I.shape)
		stokes_qu_mcmc(stokes_I, I_err, stokes_U, U_err, stokes_Q, Q_err, freq, \
		   wave2, save_loc, mcmc_I = True, mcmc_QU = True, I_0_guess = np.mean(stokes_I), a_guess = 0.1, \
		   p0_guess = (data[~np.isnan(data)][i])/100, chi0_guess = np.abs(np.mean(0.5*np.arctan2(stokes_U,stokes_Q))), rm_guess = RM[i], sigma_rm_geuss = 1, plot = True)

fit_pixels(115)












