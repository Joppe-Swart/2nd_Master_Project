"""

Author: Joppe Swart
Created: May 2023
Last modified: May 2023
Description:This script makes the brightest fit image

"""





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
			ax2.set_ylim([-1, 100*median_p0 +1])
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

			ax1.plot(lambda_mcmc, pol_ang_med, '-', color = palette(0),linewidth=2., zorder=len(pol_ang_med)+1)
			ax2.plot(lambda_mcmc, 100*pol_frac_med, '-', color = palette(0),linewidth=2.,zorder=len(pol_frac_med)+1)

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


			ax1res.errorbar(wave2, (pol_ang - fit_pol_ang(median_p0, median_RM, median_chi0, median_sigmaRM, Stokes_I, freq)), yerr=pol_ang_err,fmt='.',markersize=4,color='k',alpha=0.75)
			ax1res.plot(wave_mcmc,[0]*len(wave_mcmc),'--',color=palette(0),linewidth=1., zorder=len(pol_ang)+1, alpha=0.8)
			ax2res.errorbar(wave2, 100*(deg_pol - fit_pol_frac(median_p0, median_RM, median_chi0, median_sigmaRM, Stokes_I, freq)),yerr=100*deg_pol_err,fmt='.',markersize=4,color='k',alpha=0.75)
			ax2res.plot(wave_mcmc,[0]*len(wave_mcmc),'--',color=palette(0),linewidth=1., zorder=len(deg_pol)+1, alpha=0.8)

			plt.tight_layout()
			plt.savefig(f'{save_loc}/fit_resultspolfrac_ang.pdf', dpi = 150)
			plt.close()
