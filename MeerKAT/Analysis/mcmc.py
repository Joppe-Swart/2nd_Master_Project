"""

Author: Joppe Swart
Created: March 2023
Last modified: March 2023
Description: This script help by making nice figures of Stokes I images                   

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


def make_lin_pol(channel, x, y):
	"""
	make_lin_pol: function to find valaues per channel
	"""
	pol_I_hdu = fits.open(f'stokes_i_corr/{channel:04d}-I-image-pb.smoothed.fits')[0]
	pol_Q_hdu = fits.open(f'stokes_q_corr/{channel:04d}-Q-image-pb.smoothed.fits')[0]
	pol_U_hdu = fits.open(f'stokes_u_corr/{channel:04d}-U-image-pb.smoothed.fits')[0]

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


	if not os.path.exists('linpol_corr'):
		os.makedirss('linpol_corr')
	if os.path.exists(f'linpol_corr/{channel:04d}-linpol-image.fits'):
		linpol_hdu = fits.open(f'linpol_corr/{channel:04d}-linpol-image.fits')[0]
		lin_pol = linpol_hdu.data
		#print(lin_pol.shape)
	if not os.path.exists(f'linpol_corr/{channel:04d}-linpol-image.fits'):
		print(f'Making the linear polarization fits for channel {channel}')
		lin_pol = np.sqrt(Q_data[0]**2+U_data[0]**2)# - 8.8e-6
		linpol_hdu = fits.PrimaryHDU(lin_pol)
		linpol_hdu.writeto(f'linpol_corr/{channel:04d}-linpol-image.fits', 
				   overwrite=True)

	Q = Q_data[0][x, y]
	U = U_data[0][x, y]
	I = I_data[0][x, y]
	linpol = lin_pol[x,y]
	pol_frac = np.sqrt(Q**2+U**2)/I
	chi = 0.5*np.arctan2(U_data[0][x, y], Q_data[0][x, y])

	return freq, wave2, I, Q, U, linpol, pol_frac, chi

def make_lists(channels, x, y):
	freq_list = []
	wave2_list = []
	I_list = []
	Q_list = []
	U_list = []
	I_list_err = []
	Q_list_err = []
	U_list_err = []
	linpol_list = []
	pol_frac_list = []
	chi_list = []

	phi_hdu = fits.open(f'rmsynth_phi.fits')[0]
	phi_data = phi_hdu.data
	print(phi_data.shape)
	phi = phi_data[x, y]

	teller = 0
	for i in range(channels):
		if teller in [20, 21, 23, 36,37,38,39,40,73,74,75]:
			teller +=1
			continue
		if teller == channels:
			break
		print(f'Retrieving data for channel {teller}')
		freq, wave2, I, Q, U, linpol, pol_frac, chi = make_lin_pol(teller, x, y)
	#print(np.sqrt(Q**2 + U**2)/I, pol_frac)
		freq_list.append(freq)
		wave2_list.append(wave2)
		I_list.append(I)
		Q_list.append(Q)
		U_list.append(U)
		I_list_err.append(0.03*I)
		Q_list_err.append(0.03*Q)
		U_list_err.append(0.03*U)
		linpol_list.append(linpol)
		pol_frac_list.append(pol_frac)
		chi_list.append(chi)
		teller +=1
	return freq_list, wave2_list, I_list, I_list_err, Q_list, Q_list_err, U_list, \
		U_list_err, linpol_list, pol_frac_list, chi_list, phi


def stokes_qu_mcmc(flux, flux_err, fluxU, fluxU_err, fluxQ, fluxQ_err, pol_ang, deg_pol, freq, \
		   wave2, save_loc, mcmc_I = True, mcmc_QU = True, I_0_guess = 0.05, a_guess = 0.1, \
		   p0_guess = 0.3, chi0_guess = 0.2, rm_guess = -25,sigma_rm_geuss = 10, plot = True):
	if mcmc_I == True:
	
		pol_ang = 0.5*np.arctan2(fluxU, fluxQ) # One of the polar coordinates
		linear_pol = np.sqrt(fluxQ**2+fluxU**2)
		deg_pol = linear_pol/flux

		print(fluxQ.shape, fluxQ_err.shape, linear_pol.shape)
		pol_ang_err = np.sqrt(((((fluxU**4)/(fluxQ**6))*fluxQ_err**2)+(((fluxU**2)/(fluxQ**4))*fluxU_err**2))/((1+(fluxU/fluxQ)**2)**2))#fluxQ*fluxU_err-fluxU*fluxQ_err)/(linear_pol**2)
		linear_pol_err = np.sqrt(fluxU**2*fluxU_err**2+fluxQ**2*fluxQ_err**2)/linear_pol
		deg_pol_err = np.sqrt((linear_pol_err/flux)**2 + (linear_pol**2/flux**4)*flux_err**2)

		freq_mcmc = np.array(freq)/1e9
		freq = np.array(freq)/1e9


		freq = freq # Lijst met frequenties in GHz
		y_obs = flux # data punten 
		dy = flux_err # error in de data
		print(flux_err)

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
			plt.figure()
			plt.xticks([0.9,1.2,1.5],[0.9,1.2,1.5])
			plt.errorbar(freq, flux,yerr= flux_err, color = 'black', capsize= 3, 
			 capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5, label ='Observations')
			plt.plot(freq_mcmc, median_I_0_StokesI*(freq_mcmc)**(median_a_StokesI), \
				 color = palette(0), alpha = 1, label='MCMC fit')
			plt.fill_between(freq_mcmc,StokesI_16, StokesI_84, facecolor = palette(1), \
					 edgecolor = palette(1), alpha = 0.3)
			plt.ylabel('Flux [Jy]')
			plt.xlabel('Frequency [GHz]')
			plt.title('Stokes I')
			plt.yscale('log')
			plt.xscale('log')
			plt.legend(loc='best', frameon=True, shadow= True)
			plt.savefig(f'{save_loc}/flux_stokesI_MCMC.pdf')
			plt.close()

		np.save('stokes_i_params.npy', np.array([median_I_0_StokesI, median_a_StokesI]))



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
			if  0 < p0 < 1 and -40 < rm < 40 and  sigma_rm > 0 and -np.pi < chi0 < np.pi:
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

def run_stokes_mcmc_pixel(channels, x, y, mcmc_I = True, mcmc_QU = False):
	freq_list, wave2_list, I_list, I_err, Q_list, Q_err, U_list, U_err, \
	linpol_list, pol_frac_list, chi_list, phi = make_lists(channels, x, y)
	I_err = np.load('all_noiseI_corr.npy')[0]
	Q_err = np.load('all_noiseQ_corr.npy')[0]
	U_err = np.load('all_noiseU_corr.npy')[0]
	print(f"the shape of the error is {U_err.shape}")
	save_loc = f"Analysis_Images/Pixel{x}_{y}/"
	stokes_qu_mcmc(np.array(I_list), I_err, np.array(U_list), U_err, np.array(Q_list), \
	Q_err, np.array(chi_list), np.array(pol_frac_list), np.array(freq_list), \
	np.array(wave2_list), save_loc, mcmc_I = mcmc_I, mcmc_QU = mcmc_QU)

#run_stokes_mcmc_pixel(126, 3955, 4222, mcmc_I = True, mcmc_QU = False)

def find_flux_region(channel, regionfile, region, rmsI, rmsQ, rmsU):

	import Integrate_flux as intf

	I_image = f'stokes_i_corr/{channel:04d}-I-image-pb.smoothed.fits'
	Q_image = f'stokes_q_corr/{channel:04d}-Q-image-pb.smoothed.fits'
	U_image = f'stokes_u_corr/{channel:04d}-U-image-pb.smoothed.fits'

	fluxI, NbeamsI = intf.integratedflux(I_image, regionfile, region, hdul=None)
	uncertaintyI = intf.uncertainty_flux(I_image, fluxI, NbeamsI, region, regionfile_rms = None, rms = rmsI, delta_cal=0.0, hdul=None)
	fluxQ, NbeamsQ = intf.integratedflux(Q_image, regionfile, region, hdul=None)
	uncertaintyQ = intf.uncertainty_flux(Q_image, fluxQ, NbeamsQ, region, regionfile_rms = None, rms = rmsQ, delta_cal=0.0, hdul=None)
	fluxU, NbeamsU = intf.integratedflux(U_image, regionfile, region, hdul=None)
	uncertaintyU = intf.uncertainty_flux(U_image, fluxU, NbeamsU, region, regionfile_rms = None, rms = rmsU, delta_cal=0.0, hdul=None)

	pol_I_hdu = fits.open(f'stokes_i_corr/{channel:04d}-I-image-pb.smoothed.fits')[0]
	I_header = pol_I_hdu.header

	freq = I_header['CRVAL3']
	wave2 = (c.c/freq)**2

	Q = fluxQ
	U = fluxU
	I = fluxI
	Q_err = uncertaintyQ
	U_err = uncertaintyU
	I_err = uncertaintyI
	linpol = np.sqrt(Q**2+U**2)
	pol_frac = np.sqrt(Q**2+U**2)/I
	chi = 0.5*np.arctan2(U, Q)

	return freq, wave2, I, Q, U, I_err, Q_err, U_err, linpol, pol_frac, chi

def make_list_region(channels, regionfile, region):
	freq_list = []
	wave2_list = []
	I_list = []
	Q_list = []
	U_list = []
	I_list_err = []
	Q_list_err = []
	U_list_err = []
	linpol_list = []
	pol_frac_list = []#def make_fits_cube(channels):
#	

	chi_list = []

	rmsI = np.load('all_noiseI_corr.npy')[0]
	rmsQ = np.load('all_noiseQ_corr.npy')[0]
	rmsU = np.load('all_noiseU_corr.npy')[0]

	teller = 0
	teller_rms = 0
	for i in range(channels):
		if teller in [20, 21, 23, 36,37,38,39,40,73,74,75]:
			teller +=1
			continue
		if teller == channels:
			break
		print(f'Retrieving data for channel {teller}')
		print(rmsI[teller_rms], rmsQ[teller_rms], rmsU[teller_rms])
		freq, wave2, I, Q, U, I_err, Q_err, U_err, linpol, pol_frac, chi = \
		find_flux_region(teller, regionfile, region, rmsI[teller_rms], rmsQ[teller_rms], rmsU[teller_rms])
	#print(np.sqrt(Q**2 + U**2)/I, pol_frac)
		freq_list.append(freq)
		wave2_list.append(wave2)
		I_list.append(I)
		Q_list.append(Q)
		U_list.append(U)
		I_list_err.append(I_err)
		Q_list_err.append(Q_err)
		U_list_err.append(U_err)
		linpol_list.append(linpol)
		pol_frac_list.append(pol_frac)
		chi_list.append(chi)
		teller +=1
		teller_rms +=1
	return freq_list, wave2_list, I_list, I_list_err, Q_list, Q_list_err, \
		U_list, U_list_err, linpol_list, pol_frac_list, chi_list

def run_mcmc_region(channels, regionfile, region, save_loc, mcmc_I = True, mcmc_QU = False, plot = True, I_0_guess = 0.05, \
			a_guess = 0.1, p0_guess = 0.3, chi0_guess = 0.2, rm_guess = -15,sigma_rm_geuss = 10):
	freq_list, wave2_list, I_list, I_err, Q_list, Q_err, U_list, U_err, linpol_list, \
	pol_frac_list, chi_list = make_list_region(channels, regionfile, region)

	save_loc = save_loc
	stokes_qu_mcmc(np.array(I_list), np.array(I_err), np.array(U_list), np.array(U_err), \
	np.array(Q_list), np.array(Q_err), np.array(chi_list), np.array(pol_frac_list), \
	np.array(freq_list), np.array(wave2_list), save_loc, mcmc_I = mcmc_I, mcmc_QU = mcmc_QU, \
		 plot = plot, I_0_guess = I_0_guess, a_guess = a_guess, p0_guess = p0_guess, \
		 chi0_guess = chi0_guess, rm_guess = rm_guess, sigma_rm_geuss = sigma_rm_geuss)


#for i in range(4):
#	run_mcmc_region(126, '5regions.reg', i, save_loc = f"Analysis_Images/region/test{i}", mcmc_I = True, mcmc_QU = True, plot = True)

#run_mcmc_region(126, '5regions.reg', 3, save_loc = f"Analysis_Images/region/test{3}", mcmc_I = True, \
#		mcmc_QU = True, plot = True, I_0_guess = 0.05, a_guess = 0.1, p0_guess = 0.2, \
#		chi0_guess = 2, rm_guess = -15,sigma_rm_geuss = 25)

#run_mcmc_region(126, 'brightestsource1.reg', 0, save_loc = f"Analysis_Images/region/brightsource1", mcmc_I = True, \
#		mcmc_QU = False, plot = True, I_0_guess = 0.05, a_guess = 0.1, p0_guess = 0.2, \
#		chi0_guess = 2, rm_guess = -15,sigma_rm_geuss = 25)

run_mcmc_region(126, 'posisble_shock.reg', 0, save_loc = f"Analysis_Images/region/shock", mcmc_I = True, mcmc_QU = True, plot = True, I_0_guess = 0.05, \
			a_guess = 0.1, p0_guess = 0.3, chi0_guess = 2, rm_guess = -15,sigma_rm_geuss = 10)










