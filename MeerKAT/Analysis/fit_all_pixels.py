"""

Author: Joppe Swart
Created: March 2023
Last modified: March 2023
Description: This script fits all polarized pixels

"""

import matplotlib.patches as patches
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
import scipy.constants as c
import scipy.optimize
from matplotlib import style
import matplotlib.ticker as ticker
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

	w= wcs.WCS(f[0].header)
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

def circular_mean_std(angles):
	"""
	Calculate the mean and std of angles distributed between 0 and pi. (e.g. pol angles)
	"""
	# Convert angles to vectors on the unit circle complex plane
	vectors = np.exp(1j*angles)
	# But now vectors with angle 0 and pi are completely opposite.
	# In polarisation angle, they are the same. So mirror those on the left of the x-axis (i.e do +np.pi)
	vectors[np.real(vectors)<0] = vectors[np.real(vectors)<0] * -1
	mean_vector = np.mean(vectors)
	mean_angle = np.angle(mean_vector) % np.pi
	# Calculate standard deviation, as the deviation in the spread of the angles
	deviations = np.angle(vectors) - mean_angle
	circular_std = np.sqrt(np.mean(np.square((deviations))))
	return mean_angle, circular_std


def stokes_qu_mcmc(flux, flux_err, fluxU, fluxU_err, fluxQ, fluxQ_err, freq, \
		 wave2, save_loc, mcmc_I = True, mcmc_QU = True, I_0_guess = 0.05, a_guess = 0.1, \
		 p0_guess = 0.3, chi0_guess = 0.2, rm_guess = -25,sigma_rm_geuss = 10, plot = True):
	"""
	stokes_qu_mcmc: This function runs mcmc on a given list of flux values. 
	INPUTS:	flux: list of stokes I fluxes
		flux_err: erro on those stokes I fluxes
		fluxU: Stokes U fluxes
		fluxU_err: Stokes U error fluxes
		fluxQ: Stokes Q fluxes
		fluxQ_err: stokes Q errors
		freq: List of frequencies
		wave2: list of wavelength squared
		save_loc: direction where to save the output files
		mcmc_I: If True, run mcmc on stokes I
		mcmc_QU: If true run mcmc on stokes QU
		I_0_guess: initial guess I_0
		a_guess: initial guess a
		p0_guess: initial guess p_0
		chi0_guess: initial guess chi_0
		rm_guess: initial guess rm
		sigma_rm_geuss0: initial guess sigma rm
		plot: If true, make the plots
	OUTPUTS: best fit parameters and a bunch of images
	"""
#	def Q(p0, chi0, rm , sigma_rm, Stokes_I, x):
#		Q = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x**2)) * np.cos(2*(chi0 + rm*x))
#		return Q

#	def U(p0, chi0, rm , sigma_rm, Stokes_I, x):
#		U = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x**2)) * np.sin(2*(chi0 + rm*x))
#		return U

#	def I(I_0, a, x):
#		I = I_0*x**(a)
#		return I

	if mcmc_I == True:
		
		pol_ang = 0.5*np.arctan2(fluxU, fluxQ) # One of the polar coordinates
		linear_pol = np.sqrt(fluxQ**2+fluxU**2)
		deg_pol = linear_pol/flux



#		print(fluxQ.shape, fluxQ_err.shape, linear_pol.shape)
#		pol_ang_err = np.sqrt(((((fluxU**4)/(fluxQ**6))*fluxQ_err**2)+(((fluxU**2)/(fluxQ**4))*fluxU_err**2))/((1+(fluxU/fluxQ)**2)**2))#fluxQ*fluxU_err-fluxU*fluxQ_err)/(linear_pol**2)
#		linear_pol_err = np.sqrt(fluxU*fluxU_err+fluxQ*fluxQ_err)/linear_pol
		

		pol_ang_err = 0.5*np.sqrt((fluxQ*fluxU_err)**2+(fluxU*fluxQ_err)**2)/(linear_pol**2)
		linear_pol_err = (fluxU*fluxU_err+fluxQ*fluxQ_err)/linear_pol
#		deg_pol_err = np.sqrt((linear_pol_err/flux)**2 + (linear_pol**2/flux**4)*flux_err**2)
		deg_pol_err =1/(flux**2*deg_pol)*np.sqrt((deg_pol**2*flux*flux_err)**2+(fluxQ*fluxQ_err)**2+(fluxU*fluxU_err)**2)#np.sqrt(linear_pol_err**2 + flux_err**2)

		freq_mcmc = np.linspace(0.89, 1.81, 150)
		wave_mcmc = (c.c/(freq_mcmc*1e9))**2
		freq = np.array(freq)/1e9

		func = lambda p, i: p[0] * i**(p[1])
		guess= [1.0,-1.0] #norm, a

		errfunc = lambda p, i, j, err: (j - func(p, i)) / err
		out = scipy.optimize.leastsq(errfunc, guess, args=(freq, flux, flux_err), full_output=True)

		coeff = out[0]
		I_0_guess = coeff[0]
		a_guess = coeff[1]


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
			if I_0 >= 0. and -3 < a < 3:#and -5 < b < 5:
				return 0.0
			return -np.inf

		def lnprob_StokesI(theta, freq, y, yerr):
			lp_StokesI = lnprior_StokesI(theta)
			if not np.isfinite(lp_StokesI):
				return -np.inf
			return lp_StokesI + lnL_StokesI(theta, freq, y, yerr)

		print(f'Start the actual mcmc')
		ndim_StokesI, nwalkers = 2, 100
		theta_StokesI_guess = np.array([I_0_guess, a_guess])
		pos_StokesI = [theta_StokesI_guess + 1e-4*np.random.randn(ndim_StokesI) \
				for i in range(nwalkers)]
		sampler_StokesI = emcee.EnsembleSampler(nwalkers, ndim_StokesI, 
					lnprob_StokesI, args=(freq, y_obs, dy))
		tmp = sampler_StokesI.run_mcmc(pos_StokesI, 300)

		if plot == True:
			if not os.path.exists(f'{save_loc}'):
				os.makedirs(f'{save_loc}')

			fig, axes = plt.subplots(ncols=1, nrows=2)
			fig.set_size_inches(12,12)
			axes[0].plot(sampler_StokesI.chain[:, :, 0].transpose(), color='black', alpha=0.3)
			axes[0].set_ylabel(r'$I_0$')
			axes[0].axvline(100, ls='dashed', color='red')
			axes[1].plot(sampler_StokesI.chain[:, :, 1].transpose(), color='black', alpha=0.3)
			axes[1].set_ylabel(r'$a$')
			axes[1].axvline(100, ls='dashed', color='red')
			fig.savefig(f'{save_loc}/chain_StokesI.pdf',dpi = 50)
			plt.close()

		samples_StokesI = sampler_StokesI.chain[:, 100:, :].reshape((-1, 2))

		if plot == True:
			fig = corner.corner(samples_StokesI, labels=[r"$I_0$", r"$a$", r"b"], \
			quantiles=[0.16, 0.50, 0.84], show_titles=True)#truths=[S_pl_guess, a_pl_guess],)
			fig.savefig(f'{save_loc}/corner_StokesI.pdf', dpi=100)
			plt.close()

		median_I_0_StokesI = np.percentile(samples_StokesI[:, 0], 50.0)
		median_a_StokesI = np.percentile(samples_StokesI[:, 1], 50.0)


		p16_I_0_StokesI = np.percentile(samples_StokesI[:, 0], 16)
		p16_a_StokesI = np.percentile(samples_StokesI[:, 1], 16)


		p84_I_0_StokesI = np.percentile(samples_StokesI[:, 0], 84)
		p84_a_StokesI = np.percentile(samples_StokesI[:, 1], 84)


		sigma_I_0_StokesI = 0.5*(p84_I_0_StokesI-p16_I_0_StokesI)
		sigma_a_StokesI = 0.5*(p84_a_StokesI-p16_a_StokesI)

		if plot == True:
			MCMC_StokesI=np.empty(shape=[len(samples_StokesI[:,0]), 0])
			for i in range(len(freq_mcmc)):
				MCMC_StokesI = np.append(MCMC_StokesI, (samples_StokesI[:,0]*freq_mcmc[i]**\
				(samples_StokesI[:,1])).reshape(len(samples_StokesI[:,0]), 1), axis = 1)

			pl_16 = []
			pl_84 = []
			pl_50 = []
			for i in range(len(freq_mcmc)):
				pl_16.append(np.percentile(np.sort(MCMC_StokesI[:,i]),16))
				pl_84.append(np.percentile(np.sort(MCMC_StokesI[:,i]),84))
				pl_50.append(np.percentile(np.sort(MCMC_StokesI[:,i]),50))
				

			StokesI_16 = np.array(pl_16)#*u.mJy
			StokesI_84 = np.array(pl_84)#*u.mJy
			StokesI_50 = np.array(pl_50)

#		def red_chi_StokesI(I_0, a):
#			chi_StokesI = 0
#			model = I_0*freq**(a)
#			for i in range(len(freq)):
#				chi_StokesI += ((y_obs[i]-model[i])/dy[i])**2
#				red_chi_StokesI = chi_StokesI/(len(freq)-2)
#			return red_chi_StokesI

		if plot == True:
			I_med = StokesI_50 
			I_lower = StokesI_16 
			I_upper = StokesI_84
			fig_width = 5.78851 # half the width of A4 page in inches
			fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
			plt.figure()
			plt.errorbar(freq, 1e3*flux,yerr= 1e3*flux_err, color = 'black', capsize= 3, 
			 capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5, label ='Observations Stokes I')
			plt.plot(freq_mcmc, 1e3*I_med, \
				 color = palette(0), alpha = 1, label='MCMC fit Stokes I')
			plt.fill_between(freq_mcmc,1e3*I_lower, 1e3*I_upper, facecolor = palette(1), \
					 edgecolor = palette(1), alpha = 0.3)
			plt.ylabel('Flux [Jy]',fontsize=12)
			plt.xlabel('Frequency [GHz]', fontsize=12)
#			plt.yscale('log')
#			plt.xscale('log')
			plt.xticks(fontsize=12)
			plt.yticks(fontsize=12)
#			plt.title('Stokes I')
			plt.xlim(0.9, 1.8)
			plt.xticks([0.9,1,1.1,1.2,1.3,1.4, 1.5, 1.6, 1.7, 1.8], ['0.9', ' ', '1.1', ' ', '1.3', ' ', '1.5', ' ', '1.7', ' '], fontsize = 12)
			plt.legend(loc='best', frameon=True, shadow= True)
#			plt.xticks([0.9, 1.2, 1.6], [0.9, 1.2, 1.6], fontsize = 12)
			plt.savefig(f'{save_loc}/flux_stokesI_MCMC.pdf', dpi = 50)
			plt.close()

#		np.save('stokes_i_params.npy', np.array([median_I_0_StokesI, median_a_StokesI]))



	if mcmc_QU == True:

		Stokes_I = median_I_0_StokesI*freq**(median_a_StokesI)# + median_b_StokesI*np.log10(freq))
		Stokes_I_mcmc = median_I_0_StokesI*freq_mcmc**(median_a_StokesI)# + median_b_StokesI*np.log10(freq_mcmc))

		func = lambda p, x: 2.*(p[1] + p[2]*x) 

		guess = [0.6,1,rm_guess,1] # p0, chi0, RM, sigma_RM

		# Using Eq 6 in Gabri paper
		errfunc = lambda p, x1, y1, err1, x2, y2, err2, stokes_i: \
			abs( (y1 - (stokes_i* (p[0]*np.exp(-2*(p[3])*(x1**2))) * np.cos(func(p,x1))) )/err1 ) + \
			abs( (y2 - (stokes_i * (p[0]*np.exp(-2*(p[3])*(x1**2))) * np.sin(func(p,x2))) )/err2 )

		out = scipy.optimize.leastsq(errfunc, guess, args=(wave2,fluxQ,fluxQ_err,wave2,fluxU,fluxU_err,Stokes_I))
		coeff = out[0]

		# boundaries on p0
		if coeff[0] < 0:
			coeff[0] = -1*coeff[0]
			coeff[1] =np.pi/2 - coeff[1]
		elif coeff[0] > 1:
			coeff[0] = 1

		# boundaries on chi0
		if coeff[1] >= np.pi or coeff[1] < 0:
			coeff[1] = coeff[1] % np.pi

		p0_guess = coeff[0]
		chi0_guess = coeff[1]
		rm_guess = coeff[2]
		sigma_rm_geuss = abs(coeff[3])

		# start with QU
		wave2 = wave2

		lambda_mcmc = wave_mcmc
		x_lambda = wave2 # Lijst met frequenties in GHz
		y_U = fluxU # data punten 
		dy_U = fluxU_err# error in de data
		y_Q = fluxQ # data punten 
		dy_Q = fluxQ_err # error in de data


		def lnLQU(thetaQU, x_lambda, y1, yerr1, y2, yerr2):
			p0, chi0, rm, sigma_rm = thetaQU
			model1 = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x_lambda**2)) * np.sin(2*(chi0 + rm*x_lambda)) #Stokes U
			model2 = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x_lambda**2)) * np.cos(2*(chi0 + rm*x_lambda)) #Stokes Q

			return -0.5*(np.sum(((y1-model1)/yerr1)**2 + ((y2-model2)/yerr2)**2))

		def lnpriorQU(thetaQU):
			p0, chi0, rm, sigma_rm = thetaQU
			if 0.0 <= p0 <= 1 and -1000 <= rm <= 1000 and sigma_rm >= 0: # and 0 < chi0 < np.pi:
				return 0.0
			return -np.inf

		def lnprobQU(thetaQU, x_lambda, y1, yerr1, y2, yerr2):
			lp = lnpriorQU(thetaQU)
			if not np.isfinite(lp):
				return -np.inf
			return lp + lnLQU(thetaQU, x_lambda, y1, yerr1, y2, yerr2)

		ndim, nwalkers = 4, 100
		theta_guess = np.array([p0_guess, chi0_guess, rm_guess, sigma_rm_geuss])
		pos = [theta_guess + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
		samplerQU = emcee.EnsembleSampler(nwalkers, ndim, lnprobQU, args=(x_lambda, y_U, dy_U, y_Q, dy_Q))
		tmpQU = samplerQU.run_mcmc(pos, 500)

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
			fig.savefig(f'{save_loc}/chain_StokesQU.pdf', dpi = 50)
			plt.close()

		samplesQU = samplerQU.chain[:, 200:, :].reshape((-1, 4))
		samplesQU[:,1] = samplesQU[:,1] % np.pi # Fold back the chi0 to [0,pi]
#		print(samplesQU[:,1])
		mean_chi0 = stats.circmean(samplesQU[:,1], high=np.pi, low = 0)
		std_chi0 = stats.circstd(samplesQU[:,1], high=np.pi, low = 0)#mean_chi0, std_chi0 = circular_mean_std(samplesQU[:,1])
#		print(f"the mean chi is {mean_chi0}")


		if plot == True:
			fig = corner.corner(samplesQU, labels=[r"$p_0$", r"$\chi_0$ [rad]", r"RM [rad m$^{-2}$]", r"$\sigma^2_{RM}$[rad$^2$ m$^{-4}$]"], 
			quantiles=[0.16, 0.50, 0.84], show_titles=False, label_kwargs={"fontsize": 20},labelpad=-0.1, title_kwargs={"fontsize": 14}, title_fmt='.3f')
			fig.savefig(f'{save_loc}/corner_StokesQU.pdf', dpi = 150)
			plt.close()

		median_p0 = np.percentile(samplesQU[:, 0], 50.0)
		median_chi0 = mean_chi0 #np.percentile(samplesQU[:, 1], 50.0)
		median_RM = np.percentile(samplesQU[:, 2], 50.0)
		median_sigmaRM = np.percentile(samplesQU[:, 3], 50.0)

		p16_p0 = np.percentile(samplesQU[:, 0], 16)
		p16_chi0 = mean_chi0 - std_chi0 #np.percentile(samplesQU[:, 1], 16)
		p16_RM = np.percentile(samplesQU[:, 2], 16)
		p16_sigmaRM = np.percentile(samplesQU[:, 3], 16)

		p84_p0 = np.percentile(samplesQU[:, 0], 84)
		p84_chi0 = mean_chi0 + std_chi0 #np.percentile(samplesQU[:, 1], 84)
		p84_RM = np.percentile(samplesQU[:, 2], 84)
		p84_sigmaRM = np.percentile(samplesQU[:, 3], 84)

		sigma_p0 = 0.5*(p84_p0-p16_p0)
		sigma_chi0 = 0.5*(p84_chi0-p16_chi0)
		sigma_RM = 0.5*(p84_RM-p16_RM)
		sigma_sigmaRM = 0.5*(p84_sigmaRM-p16_sigmaRM)
		
		if plot == True:
			def Q(p0, chi0, rm , sigma_rm, Stokes_I, x):
				Q = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x**2)) * np.cos(2*(chi0 + rm*x))
				return Q

			def U(p0, chi0, rm , sigma_rm, Stokes_I, x):
				U = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x**2)) * np.sin(2*(chi0 + rm*x))
				return U

			def fit_pol_frac(p0, chi0, rm , sigma_rm, Stokes_I, x):
				pol_frac = np.sqrt(U(p0, chi0, rm , sigma_rm, Stokes_I, x)**2+Q(p0, chi0, rm , sigma_rm, Stokes_I, x)**2)/Stokes_I
				return pol_frac

			def fit_pol_ang(p0, chi0, rm , sigma_rm, Stokes_I, x):
				pol_ang = 0.5*np.arctan2(U(p0, chi0, rm , sigma_rm, Stokes_I, x),\
					Q(p0, chi0, rm , sigma_rm, Stokes_I, x))
				return pol_ang

			MCMC_U=np.empty(shape=[len(samplesQU[:,0]), 0])
			MCMC_Q=np.empty(shape=[len(samplesQU[:,0]), 0])
			MCMC_pol_ang=np.empty(shape=[len(samplesQU[:,0]), 0])
			MCMC_pol_frac=np.empty(shape=[len(samplesQU[:,0]), 0])
			for i in range(len(lambda_mcmc)):
				MCMC_U = np.append(MCMC_U, (U(samplesQU[:,0],samplesQU[:,1],samplesQU[:,2],samplesQU[:,3],Stokes_I_mcmc[i], lambda_mcmc[i])).reshape(len(samplesQU[:,0]), 1), 
							axis = 1)
				MCMC_Q = np.append(MCMC_Q, (Q(samplesQU[:,0],samplesQU[:,1],samplesQU[:,2],samplesQU[:,3],Stokes_I_mcmc[i], lambda_mcmc[i])).
						reshape(len(samplesQU[:,0]), 1), axis = 1)
				MCMC_pol_ang = np.append(MCMC_pol_ang, (fit_pol_ang(samplesQU[:,0],samplesQU[:,1],samplesQU[:,2],samplesQU[:,3],Stokes_I_mcmc[i], lambda_mcmc[i])).
						reshape(len(samplesQU[:,0]), 1), axis = 1)
				MCMC_pol_frac = np.append(MCMC_pol_frac, (fit_pol_frac(samplesQU[:,0],samplesQU[:,1],samplesQU[:,2],samplesQU[:,3],Stokes_I_mcmc[i], lambda_mcmc[i])).
						reshape(len(samplesQU[:,0]), 1), axis = 1)
#			print('the shape is',MCMC_U.shape)
			U_16 = []
			U_84 = []
			U_50 = []
			Q_16 = []
			Q_84 = []
			Q_50 = []
			pol_ang_16 = []
			pol_ang_84 = []
			pol_ang_50 = []
			pol_frac_16 = []
			pol_frac_84 = []
			pol_frac_50 = []

			for i in range(len(lambda_mcmc)):
				U_16.append(np.percentile(np.sort(MCMC_U[:,i]),16))
				U_84.append(np.percentile(np.sort(MCMC_U[:,i]),84))
				U_50.append(np.percentile(np.sort(MCMC_U[:,i]),50))
				Q_16.append(np.percentile(np.sort(MCMC_Q[:,i]),16))
				Q_84.append(np.percentile(np.sort(MCMC_Q[:,i]),84))
				Q_50.append(np.percentile(np.sort(MCMC_Q[:,i]),50))
				pol_ang_16.append(np.percentile(np.sort(MCMC_pol_ang[:,i]),16))
				pol_ang_84.append(np.percentile(np.sort(MCMC_pol_ang[:,i]),84))
				pol_ang_50.append(np.percentile(np.sort(MCMC_pol_ang[:,i]),50))
				pol_frac_16.append(np.percentile(np.sort(MCMC_pol_frac[:,i]),16))
				pol_frac_84.append(np.percentile(np.sort(MCMC_pol_frac[:,i]),84))
				pol_frac_50.append(np.percentile(np.sort(MCMC_pol_frac[:,i]),50))

			StokesQ_16 = np.array(Q_16)#*u.mJy
			StokesQ_84 = np.array(Q_84)#*u.mJy
			StokesQ_50 = np.array(Q_50)
			StokesU_16 = np.array(U_16)#*u.mJy
			StokesU_84 = np.array(U_84)#*u.mJy
			StokesU_50 = np.array(U_50)
			pol_ang_16 = np.array(pol_ang_16)#*u.mJy
			pol_ang_84 = np.array(pol_ang_84)#*u.mJy
			pol_ang_50 = np.array(pol_ang_50)
			pol_frac_16 = np.array(pol_frac_16)#*u.mJy
			pol_frac_84 = np.array(pol_frac_84)#*u.mJy
			pol_frac_50 = np.array(pol_frac_50)


#		def red_chi_U(p0, chi0, RM, sigma_RM):
#			chi = 0
#			model = Stokes_I * p0 * np.exp(-2*(sigma_RM)*(x_lambda**2)) * np.sin(2*(chi0 + RM*x_lambda))
#			for i in range(len(x_lambda)):
#				chi += ((y_U[i]-model[i])/dy_U[i])**2
#				red_chi = chi/(len(x_lambda)-4)
#			return red_chi

#		def red_chi_Q(p0, chi0, RM, sigma_RM):
#			chi = 0
#			model = Stokes_I * p0 * np.exp(-2*(sigma_RM)*(x_lambda**2)) * np.cos(2*(chi0 + RM*x_lambda))
#			for i in range(len(x_lambda)):
#				chi += ((y_Q[i]-model[i])/dy_Q[i])**2
#				red_chi = chi/(len(x_lambda)-4)
#			return red_chi
#		cov = np.cov(samplesQU, rowvar=False)
#		delta = 1.96 * np.sqrt(np.diag(cov))# 95% confidence interval

		def I(I_0, a, x):
			I = I_0*x**(a)
			return I

		def Q(p0, chi0, rm , sigma_rm, Stokes_I, x):
			Q = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x**2)) * np.cos(2*(chi0 + rm*x))
			return Q

		def U(p0, chi0, rm , sigma_rm, Stokes_I, x):
			U = Stokes_I * p0 * np.exp(-2*(sigma_rm)*(x**2)) * np.sin(2*(chi0 + rm*x))
			return U

		def fit_pol_frac(p0, chi0, rm , sigma_rm, Stokes_I, x):
			pol_frac = np.sqrt(U(p0, chi0, rm , sigma_rm, Stokes_I, x)**2+Q(p0, chi0, rm , sigma_rm, Stokes_I, x)**2)/Stokes_I
			return pol_frac

		def fit_pol_ang(p0, chi0, rm , sigma_rm, Stokes_I, x):
			pol_ang = 0.5*np.arctan2(U(median_p0, median_chi0, median_RM , median_sigmaRM, Stokes_I, x),\
				Q(median_p0, median_chi0, median_RM , median_sigmaRM, Stokes_I, x))
			return pol_ang

#		def find_upper(samples, percentile):
	
				
		if plot == True:
			Q_med = StokesQ_50 #Q(median_p0, median_chi0, median_RM, median_sigmaRM, Stokes_I_mcmc, lambda_mcmc)
			Q_lower = StokesQ_16 #Q(median_p0+delta[0], median_chi0+delta[1], median_RM+delta[2], median_sigmaRM+delta[3], Stokes_I_mcmc, lambda_mcmc)
			Q_upper = StokesQ_84 #Q(median_p0-delta[0], median_chi0-delta[1], median_RM-delta[2], median_sigmaRM-delta[3], Stokes_I_mcmc, lambda_mcmc)

			U_med = StokesU_50 #U(median_p0, median_chi0, median_RM, median_sigmaRM, Stokes_I_mcmc, lambda_mcmc)
			U_lower = StokesU_16 #U(median_p0+delta[0], median_chi0+delta[1], median_RM+delta[2], median_sigmaRM+delta[3], Stokes_I_mcmc, lambda_mcmc)
			U_upper = StokesU_84 #U(median_p0-delta[0], median_chi0-delta[1], median_RM-delta[2], median_sigmaRM-delta[3], Stokes_I_mcmc, lambda_mcmc)

			pol_frac_med = pol_frac_50 #fit_pol_frac(median_p0, median_chi0, median_RM, median_sigmaRM, Stokes_I_mcmc, lambda_mcmc)
			pol_frac_lower = pol_frac_16 #fit_pol_frac(median_p0+delta[0], median_chi0+delta[1], median_RM+delta[2], median_sigmaRM+delta[3], Stokes_I_mcmc, lambda_mcmc)
			pol_frac_upper = pol_frac_84 #fit_pol_frac(median_p0-delta[0], median_chi0-delta[1], median_RM-delta[2], median_sigmaRM-delta[3], Stokes_I_mcmc, lambda_mcmc)

			pol_ang_med = pol_ang_50 #fit_pol_ang(median_p0, median_chi0, median_RM, median_sigmaRM, Stokes_I_mcmc, lambda_mcmc)
			pol_ang_lower = pol_ang_16 #fit_pol_ang(median_p0+delta[0], median_chi0+delta[1], median_RM+delta[2], median_sigmaRM+delta[3], Stokes_I_mcmc, lambda_mcmc)
			pol_ang_upper = pol_ang_84 #fit_pol_ang(median_p0-delta[0], median_chi0-delta[1], median_RM-delta[2], median_sigmaRM-delta[3], Stokes_I_mcmc, lambda_mcmc)

#			freq_axis = [0.9, 1.2, 1.5]
#			plt.figure(figsize=(12, 8))
#			plt.subplot(2,2,1)
#			plt.errorbar(wave2, fluxU,yerr= fluxU_err, color = 'black', capsize= 3, capthick=0.5, 
#			 fmt='.', markersize= 7, elinewidth = 0.5, label ='Observations Stokes U')
#			plt.plot(lambda_mcmc, U_med, color = palette(0), alpha = 1, label='MCMC fit Stokes U')
#			plt.fill_between(lambda_mcmc,U_lower,U_upper, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
#			plt.ylabel('Flux [Jy]')
#			plt.xlabel(r' $\lambda^2 \, [m^2]$')
##			plt.title('Stokes U')
#			plt.xlim(0.0275, 0.105)
#			plt.legend(loc='best', frameon=True, shadow= True)
##			plt.savefig(f'{save_loc}/flux_stokesU_MCMC.pdf', dpi = 50)
##			plt.close()

#			plt.subplot(2,2,2)
#			plt.errorbar(wave2, fluxQ,yerr= fluxQ_err, color = 'black', capsize= 3, capthick=0.5, 
#			 fmt='.', markersize= 7, elinewidth = 0.5, label ='Observations Stokes Q')
#			plt.plot(lambda_mcmc, Q_med, color = palette(0), alpha = 1, label='MCMC fit Stokes Q')
#			plt.fill_between(lambda_mcmc,Q_lower, Q_upper, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
#			plt.ylabel('Flux [Jy]')
#			plt.xlabel(r' $\lambda^2 \, [m^2]$')
##			plt.title('Stokes Q')
#			plt.xlim(0.0275, 0.105)
#			plt.legend(loc='best', frameon=True, shadow= True)
##			plt.savefig(f'{save_loc}/flux_stokesQ_MCMC.pdf', dpi = 50)
##			plt.close()

#			plt.subplot(2,2,3)
#			plt.errorbar(wave2, pol_ang ,yerr= pol_ang_err, color = 'black', capsize= 3, capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5, label =r'Observations $\chi$')
#			plt.plot(lambda_mcmc, pol_ang_med, color = palette(0), alpha = 1, label=r'MCMC fit $\chi$')
#			plt.fill_between(lambda_mcmc,pol_ang_lower, pol_ang_upper, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
#			plt.ylabel(r'Polarization Angle $(\chi)$ [rad]')
##			plt.ylim(-np.pi, np.pi)
#			plt.xlabel(r' $\lambda^2 \, [m^2]$')
##			plt.title('Polarization angle from individual mcmc procedure')
#			plt.xlim(0.0275, 0.105)
#			plt.legend(loc='best', frameon=True, shadow= True)
##			plt.savefig(f'{save_loc}/pol_ang_MCMC.pdf', dpi = 50)
##			plt.close()


#			plt.subplot(2,2,4)
#			plt.errorbar(wave2, 100*deg_pol,yerr= 100*deg_pol_err, color = 'black', \
#					capsize= 3, capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5, label=r'Polarisation degree')
#			plt.plot(lambda_mcmc, 100*pol_frac_med, color = palette(0), alpha = 1, label='MCMC fit polarisation degree')
#			plt.fill_between(lambda_mcmc,100*pol_frac_lower, 100*pol_frac_upper, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
#			plt.ylabel(r'p [$\%$]')
#			plt.xlabel(r' $\lambda^2 \, [m^2]$')
##			plt.title('Degree of polarization')
#			plt.xlim(0.0275, 0.105)
#			plt.legend(loc = 'best', frameon=True, shadow= True)
#			plt.savefig(f'{save_loc}/fit_results.pdf', dpi = 50)
#			plt.close()
			
			def models(x, norm, a, p0, X0, RM, sigma, obs_i, err_i, obs_q, err_q, obs_u, err_u):
				freq = (3.e8/np.sqrt(x))*1.e-9
				mod_i = norm * freq**(a)

				if len(x) == len(obs_i):
					dof = len(x)-2
					chisqred_i = (sum( ((obs_i - mod_i)/err_i)**2 ))/dof
				else:
					chisqred_i = None


				mod_q = p0*np.exp(-2*(sigma)*(x**2))*mod_i *np.cos(2.* (X0 + RM*x))
				mod_u = p0*np.exp(-2*(sigma)*(x**2))*mod_i *np.sin(2.* (X0 + RM*x))

				if len(x) == len(obs_q):
					dof = 2*(len(x)-4)
					chisqred_qu = (sum( ((obs_q - mod_q)/err_q)**2 + ((obs_u - mod_u)/err_u)**2 ))/dof
				else:
					chisqred_qu = None


				mod_X = 0.5*np.arctan2(mod_u,mod_q)
				mod_p = np.sqrt(mod_q**2 + mod_u**2) / mod_i

				return(mod_i, mod_q, mod_u, chisqred_i, chisqred_qu, mod_X, mod_p)

			def freq_to_lambda(freq):
				wave_sq = (3.e8/freq)**2
				return(wave_sq)

			freqs = np.array([1.0,1.2,1.5])*1.e9
			freqticks = np.array([freq_to_lambda(freq) for freq in freqs])
			
			fig = plt.figure(figsize=(10,14))
			gs = gridspec.GridSpec(6,1,height_ratios = [4,1,4,1,4,1], hspace=0.)
			ax1 = fig.add_subplot(gs[0])
			ax2 = fig.add_subplot(gs[2],sharex = ax1)
			ax3 = fig.add_subplot(gs[4],sharex = ax1)
			ax4 = ax1.twiny() 
			plt.tight_layout(pad=6)

			ax1.grid(color='lightgrey', linestyle=':', linewidth=0.7)
			ax2.grid(color='lightgrey', linestyle=':', linewidth=0.7)
			ax3.grid(color='lightgrey', linestyle=':', linewidth=0.7)

			ax1.set_ylabel(r'$I$ [$\mu$Jy]', fontsize=28)
			ax2.set_ylabel(r'$Q$ [$\mu$Jy]', fontsize=28)
			ax3.set_ylabel(r'$U$ [$\mu$Jy]', fontsize=28)
			ax4.set_xlabel(r'$\nu$ [GHz]', fontsize=28)

			ax1.tick_params(axis='y',labelsize=20)
			ax2.tick_params(axis='y',labelsize=20)
			ax3.tick_params(axis='y',labelsize=20)
			ax4.tick_params(axis='x',labelsize=20)

			ax2.xaxis.set_ticks_position('both')
			ax3.xaxis.set_ticks_position('both')
			ax1.yaxis.set_ticks_position('both')
			ax2.yaxis.set_ticks_position('both')
			ax3.yaxis.set_ticks_position('both')
			#ax2.yaxis.set_ticks([-0.2,-0.1,0,0.1,0.2])


			ax1.set_xlim([0.0275, 0.105])
			ax4.set_xlim([0.0275, 0.105])

			ax4.set_xticks(freqticks)
			ax4.set_xticklabels('{:.1f}'.format(freq) for freq in freqs/1.e9)

			xticklabels1 = ax1.get_xticklabels()
			xticklabels2 = ax2.get_xticklabels()
			xticklabels3 = ax3.get_xticklabels()
			plt.setp(xticklabels1, visible=True)
			plt.setp(xticklabels2, visible=False)
			plt.setp(xticklabels3, visible=False)
			#plt.suptitle(sourcetitles[row], fontsize=18)

			ax1.errorbar(x_lambda, 1e6*flux, yerr=1e6*flux_err, color = 'black', capsize= 3, capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5)
			ax1.fill_between(lambda_mcmc,1e6*I_lower, 1e6*I_upper, interpolate=True, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
			ax2.errorbar(x_lambda, 1e6*fluxQ, yerr=1e6*fluxQ_err, color = 'black', capsize= 3, capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5)
			ax2.fill_between(lambda_mcmc,1e6*Q_lower, 1e6*Q_upper, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
			ax3.errorbar(x_lambda, 1e6*fluxU, yerr=1e6*fluxU_err, color = 'black', capsize= 3, capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5)
			ax3.fill_between(lambda_mcmc,1e6*U_lower, 1e6*U_upper, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)


			# print overlay models
			x = np.linspace(0.0275, 0.105, 200)
			fit_i = models(x, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[0]
			fit_q = models(x, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[1]
			fit_u = models(x, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[2]

			chisqred_i= models(wave2, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[3]
			chisqred_qu = models(wave2, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[4]

			leg1, = ax1.plot(x, 1e6*fit_i,'-',color=palette(0),linewidth=2.,zorder=len(flux)+1)
			leg2, = ax2.plot(x, 1e6*fit_q,'-',color=palette(0),linewidth=2.,zorder=len(fluxQ)+1)
			ax3.plot(x, 1e6*fit_u,'-',color=palette(0),linewidth=2.,zorder=len(fluxU)+1)


			## print res
			ax1res = fig.add_subplot(gs[1],sharex = ax1)
			ax2res = fig.add_subplot(gs[3],sharex = ax1)
			ax3res = fig.add_subplot(gs[5],sharex = ax1)

			ax1res.set_xlim([0.0275, 0.105]) 
			ax2res.set_xlim([0.0275, 0.105]) 
			ax3res.set_xlim([0.0275, 0.105]) 
			#axres.set_ylim([-0.3, 0.3])
			#ax2res.set_ylim([-0.3, 0.3])
			#ax3res.set_ylim([-0.15, 0.15])


			ax1res.grid(color='lightgrey', linestyle=':', linewidth=0.7)
			ax2res.grid(color='lightgrey', linestyle=':', linewidth=0.7)
			ax3res.grid(color='lightgrey', linestyle=':', linewidth=0.7)

			ax1res.set_ylabel('res', fontsize=24)
			ax2res.set_ylabel('res', fontsize=24)
			ax3res.set_ylabel('res', fontsize=24)
			ax3res.set_xlabel(r'$\lambda^2$ [m$^2$]', fontsize=28)

			ax1res.tick_params(axis='y',labelsize=16)
			ax2res.tick_params(axis='y',labelsize=16)
			ax3res.tick_params(axis='y',labelsize=16)
			ax3res.tick_params(axis='x',labelsize=20)

			ax1res.xaxis.set_ticks_position('both')
			ax2res.xaxis.set_ticks_position('both')
			ax3res.xaxis.set_ticks_position('both')
			ax1res.yaxis.set_ticks_position('both')
			ax2res.yaxis.set_ticks_position('both')
			ax3res.yaxis.set_ticks_position('both')

			xticklabels1 = ax1res.get_xticklabels()
			xticklabels2 = ax2res.get_xticklabels()
			xticklabels3 = ax3res.get_xticklabels()
			plt.setp(xticklabels1, visible=False)
			plt.setp(xticklabels2, visible=False)

			model_i = models(wave2, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[0]
			model_q = models(wave2, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[1]
			model_u = models(wave2, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[2]


			ax1res.errorbar(wave2, 1e6*(flux-model_i), yerr=1e6*flux,fmt='.',markersize=4,color='k',alpha=0.75)
			ax1res.plot(wave_mcmc,[0]*len(wave_mcmc),'--',color=palette(0),linewidth=1., zorder=len(flux)+1, alpha=0.8)
			ax2res.errorbar(wave2, 1e6*(fluxQ-model_q), yerr=1e6*fluxQ_err,fmt='.',markersize=4,color='k',alpha=0.75)
			ax2res.plot(wave_mcmc,[0]*len(wave_mcmc),'--',color=palette(0),linewidth=1., zorder=len(fluxQ)+1, alpha=0.8)
			ax3res.errorbar(wave2, 1e6*(fluxU-model_u), yerr=1e6*fluxU,fmt='.',markersize=4,color='k',alpha=0.75)
			ax3res.plot(wave_mcmc,[0]*len(wave_mcmc),'--',color=palette(0),linewidth=1., zorder=len(fluxU)+1, alpha=0.8)


			plt.tight_layout()
			plt.savefig(f'{save_loc}/fit_resultsIQU.pdf', dpi = 50)
			plt.close()



			fig = plt.figure(figsize=(10,14))
			gs = gridspec.GridSpec(4,1,height_ratios = [4,1,4,1], hspace=0)
			ax1 = fig.add_subplot(gs[0])
			ax2 = fig.add_subplot(gs[2], sharex = ax1)
			ax3 = ax1.twiny()

			ax1.grid(color='lightgrey', linestyle=':', linewidth=0.7)
			ax2.grid(color='lightgrey', linestyle=':', linewidth=0.7)

			ax1.set_ylabel(r'$\chi = \frac{1}{2} \arctan \left ( \frac{U}{Q} \right )$ [rad]', fontsize=28)
			ax2.set_ylabel(r'$p = \frac{\sqrt{Q^2 + U^2}}{I}$', fontsize=28) 
			ax3.set_xlabel(r'$\nu$ [GHz]', fontsize=28)

			ax1.tick_params(axis='y',labelsize=20)
			ax2.tick_params(axis='y',labelsize=20)
			ax3.tick_params(axis='x',labelsize=20) 

			ax1.yaxis.set_ticks_position('both')
			ax2.yaxis.set_ticks_position('both')
			ax2.xaxis.set_ticks_position('both')

			ax1.set_xlim([0.0275, 0.105])
			ax3.set_xlim([0.0275, 0.105])
			ax1.set_ylim([-0.5*np.pi-0.5, 0.5*np.pi+0.5])
#			ax2.set_ylim([-1, 100*median_p0 +1])
			#ax2.set_ylim([-0.05,0.65])
			#ax2.set_yticks([0.,0.1,0.2,0.3,0.4,0.5,0.6])

			ax3.set_xticks(freqticks)
			ax3.set_xticklabels('{:.1f}'.format(freq) for freq in freqs/1.e9)

			xticklabels1 = ax1.get_xticklabels()
			xticklabels2 = ax2.get_xticklabels()
			plt.setp(xticklabels1, visible=False)
			plt.setp(xticklabels2, visible=False)
			#plt.suptitle(sourcetitles[row], fontsize=18)



			ax1.errorbar(wave2, pol_ang, yerr=pol_ang_err, color = 'black', capsize= 3, capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5)
			ax1.fill_between(lambda_mcmc,pol_ang_lower, pol_ang_upper, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)
			ax2.errorbar(x_lambda, 100*deg_pol, yerr=100*deg_pol_err, color = 'black', capsize= 3, capthick=0.5, fmt='.', markersize= 7, elinewidth = 0.5)
			ax2.fill_between(lambda_mcmc,100*pol_frac_lower, 100*pol_frac_upper, facecolor = palette(1), edgecolor = palette(1), alpha = 0.3)

			fit_X = models(x, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[5]
			fit_p = models(x, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[6]

			ax1.plot(x, fit_X, '-', color = palette(0),linewidth=2., zorder=len(pol_ang_med)+1)
			ax2.plot(x, 100*fit_p, '-', color = palette(0),linewidth=2.,zorder=len(pol_frac_med)+1)

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

			model_X = models(wave2, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[5]
			model_p = models(wave2, median_I_0_StokesI, median_a_StokesI, median_p0, median_chi0, median_RM, median_sigmaRM, flux, flux_err, fluxQ, fluxQ_err, fluxU, fluxU_err)[6]

			ax1res.errorbar(wave2, (pol_ang - model_X), yerr=pol_ang_err,fmt='.',markersize=4,color='k',alpha=0.75)
			ax1res.plot(wave_mcmc,[0]*len(wave_mcmc),'--',color=palette(0),linewidth=1., zorder=len(pol_ang)+1, alpha=0.8)
			ax2res.errorbar(wave2, 100*(deg_pol - model_p),yerr=100*deg_pol_err,fmt='.',markersize=4,color='k',alpha=0.75)
			ax2res.plot(wave_mcmc,[0]*len(wave_mcmc),'--',color=palette(0),linewidth=1., zorder=len(deg_pol)+1, alpha=0.8)

			plt.tight_layout()
			plt.savefig(f'{save_loc}/fit_resultspolfrac_ang.pdf', dpi = 50)
			plt.close()


	print(f"the medium that stokesQU returns = {median_chi0}")
	return median_I_0_StokesI, median_a_StokesI, median_RM, median_chi0, median_p0, median_sigmaRM, sigma_I_0_StokesI, sigma_a_StokesI, sigma_p0, sigma_chi0, sigma_RM, sigma_sigmaRM



#def fit_pixels(channels):
#	"""
#	fit_pixels: funciton that fits all the pixels in a given refgion
#	"""
#	
#	#first find which pixels we want to fit
#	pol_hdu = fits.open('DATA/polarization_fraction_map_masked.fits')[0]
##	print(pol_hdu.data.shape)
#	mask = np.full(pol_hdu.data.shape, np.nan)
##	print(mask.shape)
#	mask[0,0,3950:3955,4220:4225] = 1
#	data = pol_hdu.data*mask
##	print(data.shape)

#	Q_cube = fits.open('DATA/Q_cube.fits')
#	U_cube = fits.open('DATA/U_cube.fits')
#	I_cube = fits.open('DATA/I_cube.fits')


#	j = np.where(data[~np.isnan(data)])
#	print(np.array(j).size)
##	print(data[~np.isnan(data)][0])
##	print(f'j = {j}')
#	fluxQ = []
#	fluxU = []
#	fluxI = []
#	freq = []

#	for i in range(channels):
#		freq.append(I_cube[i].header['CRVAL3'])
#		
#		Q_data = Q_cube[i].data*mask
#		Q = Q_data[~np.isnan(Q_data)]
#		fluxQ.append(Q[j])
#		
#		U_data = U_cube[i].data*mask
#		U = U_data[~np.isnan(U_data)]
#		fluxU.append(U[j])
#		
#		I_data = I_cube[i].data*mask
#		I = I_data[~np.isnan(I_data)]
#		fluxI.append(I[j])
#	
#	Q_cube.close()
#	U_cube.close()
#	I_cube.close()

#	fluxI = np.array(fluxI)
#	fluxU = np.array(fluxU)
#	fluxQ = np.array(fluxQ)

#	with fits.open('DATA/rmsynth_DATA/rmsynth_phi.fits') as rm:
#		rm_data = rm[0].data*mask
#		RM = rm_data[~np.isnan(rm_data)]
#	print(f'the rm is {RM[j]}')
#	freq = np.array(freq)
#	wave2 = (c.c/freq)**2

#	I_err = np.load('DATA/all_noiseI_corr.npy')[0]
#	Q_err = np.load('DATA/all_noiseQ_corr.npy')[0]
#	U_err = np.load('DATA/all_noiseU_corr.npy')[0]
#	
#	best_params = np.ones([])
#	for i in range(np.array(j).size):
#		save_loc = f'Resultsfit_pixel/test/pixel{i}'
#		print(i)
#		stokes_I = fluxI[:,i]
#		stokes_Q = fluxQ[:,i]
#		stokes_U = fluxU[:,i]
##		print(stokes_I.shape)
#		stokes_qu_mcmc(stokes_I, I_err, stokes_U, U_err, stokes_Q, Q_err, freq, \
#		 wave2, save_loc, mcmc_I = True, mcmc_QU = True, I_0_guess = np.mean(stokes_I), a_guess = 0.1, \
#		 p0_guess = (data[~np.isnan(data)][i])/100, chi0_guess = np.abs(np.mean(0.5*np.arctan2(stokes_U,stokes_Q))), rm_guess = RM[i], sigma_rm_geuss = 1, plot = True)

#fit_pixels(115)

def loop_pixels(channels, x_min = 3950, x_max = 3955, y_min = 4220, y_max = 4225, loc =f'Results/Images_Fitted_Pixels/all_pixels_final/', plot = True):
	"""
	loop_pixels: Function that loops over all pixels and fits the p0, chi, rm and sigma rm
	"""
	# Read in the input fits file
#	input_file = 'DATA/polarisation_maps/polarization_fraction_map_masked.fits'
	input_file = 'DATA/polarisation_maps/diffem_polarization_fraction_map_masked.fits'
	hdulist = fits.open(input_file)
	data = hdulist[0].data
	
	inputmask = np.full(data.shape, np.nan)
	inputmask[0,0,x_min:x_max, y_min:y_max] = 1

	data = data*inputmask
	
	dataI = data.copy()
	data_a = data.copy()
	dataI_uncrt = data.copy()
	data_a_uncrt = data.copy()
	datap = data.copy()
	datachi = data.copy()
	datarm = data.copy()
	datasigmarm = data.copy()
	datap_uncrt = data.copy()
	datachi_uncrt = data.copy()
	datarm_uncrt = data.copy()
	datasigmarm_uncrt = data.copy()

	x_min = x_min
	x_max = x_max
	y_min = y_min
	y_max = y_max

	Q_cube = fits.open('DATA/Q_cube.fits')
	U_cube = fits.open('DATA/U_cube.fits')
	I_cube = fits.open('DATA/I_cube.fits')

	I_err = np.load('DATA/all_noiseI_corr.npy')[0]
	Q_err = np.load('DATA/all_noiseQ_corr.npy')[0]
	U_err = np.load('DATA/all_noiseU_corr.npy')[0]

	rm_synth = fits.open('DATA/rmsynth/Results_full_field_phi.fits')
	rm_data = rm_synth[0].data
#	print(rm_data.shape)



	# Loop over all pixels in the image
	pixels = 0
	for i in range(x_min, x_max):
		for j in range(y_min, y_max):
			if not np.isnan(data[0,0,i,j]):
				pixels +=1
				fluxQ = []
				fluxU = []
				fluxI = []
				freq = []
				for k in range(channels):
					freq.append(I_cube[k].header['CRVAL3'])

					Q_data = Q_cube[k].data[0,0,i,j]
					fluxQ.append(Q_data)

					U_data = U_cube[k].data[0,0,i,j]
					fluxU.append(U_data)

					I_data = I_cube[k].data[0,0,i,j]
					fluxI.append(I_data)
				print(f"these pixels are polarization{pixels}")
				freq = np.array(freq)
				wave2 = (c.c/freq)**2
				fluxI = np.array(fluxI)
				fluxU = np.array(fluxU)
				fluxQ = np.array(fluxQ)

				save_loc = f'{loc}/pixel{pixels}'
				median_I_0_StokesI, median_a_StokesI, median_RM, median_chi0, median_p0, \
				median_sigmaRM, sigma_I_0_StokesI, sigma_a_StokesI, sigma_p0, sigma_chi0, \
				sigma_RM, sigma_sigmaRM = stokes_qu_mcmc(fluxI, I_err, fluxU, U_err, fluxQ, \
				Q_err, freq, wave2, save_loc, mcmc_I=True, mcmc_QU=True, plot=plot)#,I_0_guess=np.mean(fluxI), a_guess=-0.1,\
#				p0_guess=np.clip(data[0,0,i,j]/100, 0.00001, 0.999999), \
#				chi0_guess=np.clip(np.abs(np.mean(0.5*np.arctan2(fluxU,fluxQ))), 0.00001, np.pi-0.00001),\
#				rm_guess=np.clip(rm_data[i,j], -99.9999, 99.9999), sigma_rm_geuss=1, plot=plot)
				print(f"The median that gets append = {median_chi0}")

				#Write the mcmc loop here
				datap[0,0,i,j] = median_p0
				datachi[0,0,i,j] = median_chi0
				datarm[0,0,i,j] = median_RM
				datasigmarm[0,0,i,j] = median_sigmaRM
				dataI[0,0,i,j] = median_I_0_StokesI
				data_a[0,0,i,j] = median_a_StokesI
				dataI_uncrt[0,0,i,j] = sigma_I_0_StokesI
				data_a_uncrt[0,0,i,j] = sigma_a_StokesI
				datap_uncrt[0,0,i,j] = sigma_p0
				datachi_uncrt[0,0,i,j] = sigma_chi0
				datarm_uncrt[0,0,i,j] = sigma_RM
				datasigmarm_uncrt[0,0,i,j] = sigma_sigmaRM

	Q_cube.close()
	U_cube.close()
	I_cube.close()
	rm_synth.close()



	# Write the new pixel values to a new fits file with the same dimensions as the input file
	output_fileI_uncrt = f'{loc}/fitted_I0_uncrt.fits'
	hduI_uncrt = fits.PrimaryHDU(dataI_uncrt)
	hduI_uncrt.header = hdulist[0].header
	hduI_uncrt.writeto(output_fileI_uncrt, overwrite=True)

	output_file_a_uncrt = f'{loc}/fitted_a_uncrt.fits'
	hdu_a_uncrt = fits.PrimaryHDU(data_a_uncrt)
	hdu_a_uncrt.header = hdulist[0].header
	hdu_a_uncrt.writeto(output_file_a_uncrt, overwrite=True)

	output_fileI = f'{loc}/fitted_I0.fits'
	hduI = fits.PrimaryHDU(dataI)
	hduI.header = hdulist[0].header
	hduI.writeto(output_fileI, overwrite=True)

	output_file_a = f'{loc}/fitted_a.fits'
	hdu_a = fits.PrimaryHDU(data_a)
	hdu_a.header = hdulist[0].header
	hdu_a.writeto(output_file_a, overwrite=True)

	output_filep = f'{loc}/fitted_p0.fits'
	hdup = fits.PrimaryHDU(datap)
	hdup.header = hdulist[0].header
	hdup.writeto(output_filep, overwrite=True)

	output_filechi = f'{loc}/fitted_chi0.fits'
	hduchi = fits.PrimaryHDU(datachi)
	hduchi.header = hdulist[0].header
	hduchi.writeto(output_filechi, overwrite=True)

	output_filerm = f'{loc}/fitted_rm.fits'
	hdurm = fits.PrimaryHDU(datarm)
	hdurm.header = hdulist[0].header
	hdurm.writeto(output_filerm, overwrite=True)

	output_file_sigmarm = f'{loc}/fitted_sigmarm.fits'
	hdusigmarm = fits.PrimaryHDU(datasigmarm)
	hdusigmarm.header = hdulist[0].header
	hdusigmarm.writeto(output_file_sigmarm, overwrite=True)

	output_filep_uncrt = f'{loc}/fitted_p0_uncrt.fits'
	hdup_uncrt = fits.PrimaryHDU(datap_uncrt)
	hdup_uncrt.header = hdulist[0].header
	hdup_uncrt.writeto(output_filep_uncrt, overwrite=True)

	output_filechi_uncrt = f'{loc}/fitted_chi0_uncrt.fits'
	hduchi_uncrt = fits.PrimaryHDU(datachi_uncrt)
	hduchi_uncrt.header = hdulist[0].header
	hduchi_uncrt.writeto(output_filechi_uncrt, overwrite=True)

	output_filerm_uncrt = f'{loc}/fitted_rm_uncrt.fits'
	hdurm_uncrt = fits.PrimaryHDU(datarm_uncrt)
	hdurm_uncrt.header = hdulist[0].header
	hdurm_uncrt.writeto(output_filerm_uncrt, overwrite=True)

	output_file_sigmarm_uncrt = f'{loc}/fitted_sigmarm_uncrt.fits'
	hdusigmarm_uncrt = fits.PrimaryHDU(datasigmarm_uncrt)
	hdusigmarm_uncrt.header = hdulist[0].header
	hdusigmarm_uncrt.writeto(output_file_sigmarm_uncrt, overwrite=True)

#	I_hdul = fits.open(f"DATA/bullet_cluster_pb_corr.smoothed.fits")
#	I_hdu = I_hdul[0]
#	I_header = I_hdu.header
#	I_data = I_hdu.data
#	I_wcs = WCS(flatten(I_hdul).header)
#	
#	pb_data = fits.open('DATA/0055-I-pb_model.fits')[0].data
#	pbmask = np.full(pb_data.shape, np.nan)
#	pbmask[np.where(pb_data>0.9)] = 1

#	imagenoise_I = findrms(I_data*pbmask)

#	lev_factor = 3.
#	lev_radio = np.sqrt([1, 4, 16, 64, 256, 1024, 4096]) #5, 25,100
#	level_radio = np.ndarray.tolist(lev_factor*imagenoise_I*lev_radio)


#	pixelscale = 1.1 #arcsec/pixels
#	bmaj = 6.4/pixelscale # in pixels
#	bmin = 6/pixelscale
#	bpa = 50.5 #in degrees

#	bmaj_convolved = 10.4/pixelscale
#	bmin_convolved = 9.6/pixelscale
#	bpa_convolved = 61.1


#	fig_width = 4.134 # half the width of A4 page in inches
#	fig_height = fig_width * 0.75 # height = width * 0.75 to maintain aspect ratio
#	fig = plt.figure(figsize=(fig_width, fig_height))
#	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
#	map1 = plt.imshow(100*datap[0,0,:,:], cmap='rainbow',vmin=15, vmax=90, interpolation='nearest')
#	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
#	cbar.ax.tick_params(labelsize=9)
#	ax.set_xlabel('RA (J2000)', fontsize=12, labelpad=0.5)
#	ax.set_ylabel('DEC (J2000)',fontsize=12, labelpad=-0.5)
#	ax.tick_params(axis='both', which='major', labelsize=9)
#	plt.ylim(y_min, y_max)
#	plt.xlim(x_min, x_max)
#	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
#	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
#	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,facecolor='black', edgecolor='black')
#	ax.add_patch(beam_ellipse)
#	#plt.title(r'Fitted Polarisation fraction')
#	plt.tight_layout()
#	plt.savefig(f'{loc}/fitted_p0.pdf', dpi=150)
#	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.5, linewidths=0.1)
#	plt.savefig(f'{loc}/fitted_p0_contour.pdf', dpi=150)
#	plt.close()

#	fig = plt.figure(figsize=(fig_width, fig_height))
#	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
#	map1 = plt.imshow(datarm[0,0,:,:], cmap='rainbow',vmin=-50, vmax=50, interpolation='nearest')
#	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
#	cbar.ax.tick_params(labelsize=9)
#	ax.set_xlabel('RA (J2000)', fontsize=12, labelpad=0.5)
#	ax.set_ylabel('DEC (J2000)',fontsize=12, labelpad=-0.5)
#	ax.tick_params(axis='both', which='major', labelsize=9)
#	plt.ylim(y_min, y_max)
#	plt.xlim(x_min, x_max)
#	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
#	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
#	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,facecolor='black', edgecolor='black')
#	ax.add_patch(beam_ellipse)
#	#plt.title(r'Fitted RM')
#	plt.tight_layout()
#	plt.savefig(f'{loc}/fitted_RM.pdf', dpi=150)
#	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.5, linewidths=0.1)
#	plt.savefig(f'{loc}/fitted_RM_contour.pdf', dpi=150)
#	plt.close()

#	fig = plt.figure(figsize=(fig_width, fig_height))
#	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
#	map1 = plt.imshow(datachi[0,0,:,:], cmap='rainbow',vmin=0, vmax=3.14, interpolation='nearest')
#	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
#	cbar.ax.tick_params(labelsize=9)
#	ax.set_xlabel('RA (J2000)', fontsize=12, labelpad=0.5)
#	ax.set_ylabel('DEC (J2000)',fontsize=12, labelpad=-0.5)
#	ax.tick_params(axis='both', which='major', labelsize=9)
#	plt.ylim(y_min, y_max)
#	plt.xlim(x_min, x_max)
#	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
#	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
#	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,facecolor='black', edgecolor='black')
#	ax.add_patch(beam_ellipse)
#	#plt.title(r'Fitted $\chi$')
#	plt.tight_layout()
#	plt.savefig(f'{loc}/fitted_chi.pdf', dpi=150)
#	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.5, linewidths=0.1)
#	plt.savefig(f'{loc}/fitted_chi_contour.pdf', dpi=150)
#	plt.close()

#	fig = plt.figure(figsize=(fig_width, fig_height))
#	ax = plt.subplot(projection=I_wcs, slices=('x', 'y'))
#	map1 = plt.imshow(datasigmarm[0,0,:,:], cmap='rainbow',vmin=0, vmax=350, interpolation='nearest')
#	cbar = plt.colorbar(map1, pad=0.01)#label = r'\%p')
#	cbar.ax.tick_params(labelsize=9)
#	ax.set_xlabel('RA (J2000)', fontsize=12, labelpad=0.5)
#	ax.set_ylabel('DEC (J2000)',fontsize=12, labelpad=-0.5)
#	ax.tick_params(axis='both', which='major', labelsize=9)
#	plt.ylim(y_min, y_max)
#	plt.xlim(x_min, x_max)
#	beam_xpos = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1]-ax.get_xlim()[0])
#	beam_ypos = ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1]-ax.get_ylim()[0])
#	beam_ellipse = plt.matplotlib.patches.Ellipse((beam_xpos, beam_ypos), width=bmin_convolved, height=bmaj_convolved, angle=bpa_convolved, linewidth=0.5,facecolor='black', edgecolor='black')
#	ax.add_patch(beam_ellipse)
#	#plt.title(r'Fitted $\sigma_{RM}$')
#	plt.tight_layout()
#	plt.savefig(f'{loc}/fitted_sigmarm.pdf', dpi=150)
#	plt.contour(I_data[0,0,:,:], levels=level_radio, colors='black', alpha=0.5, linewidths=0.1)
#	plt.savefig(f'{loc}/fitted_sigmarm_contour.pdf', dpi=150)
#	plt.close()
#	# Close the input fits file
	hdulist.close()


#loop_pixels(115)
loop_pixels(115, x_min = 4400, x_max = 4630, y_min = 4120, y_max = 4350, loc =f'Results/Images_Fitted_Pixels/diff_em_without_chi', plot=True)
#loop_pixels(115, x_min = 3850, x_max = 4350, y_min = 3800, y_max = 4300, loc =f'Results/Images_Fitted_Pixels/all_pixels_final', plot=False)

#loop_pixels(115, x_min = 3960, x_max = 3985, y_min = 4210, y_max = 4235, loc =f'Results/Images_Fitted_Pixels/test_without_chi', plot=True)
#



