"""
Script to manually selfcal the data, the final images are made with facet_selfcal
"""


#Start Selfcal 

flagdata(vis='BC_Target_Block_I.ms', mode='manual', flagbackup=True, spw='0:712~732')

#first clean
wsclean = "wsclean -reorder -name Bullet_Cluster_selfcal1 -size 8192 8192 -scale 1.1arcsec -data-column DATA -weight briggs -0.5 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -join-channels -gap-channel-division -channels-out 6 -pol I -multiscale -parallel-deconvolution 1024 -use-wgridder -mem 50 BC_Target_Block_I.ms BC_Target_Block_II.ms BC_Target_Block_III.ms"

singularity = "singularity shell -B /tmp,/dev/shm,/disks/paradata,/data1,/data2,/net/rijn,/net/rijn9/,/software /net/achterrijn/data1/sweijen/software/containers/lofar_sksp_rijnX.sif --noprofile --norc"

subprocess.call(singularity+" "+wsclean, shell = True)


# First round of selfcal
gaincal(vis='BC_Target_Block_I.ms',caltable='BC_Target_Block_I.G0', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_II.ms',caltable='BC_Target_Block_II.G0', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_III.ms',caltable='BC_Target_Block_III.G0', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

#gaincal(vis='Bullet_Cluster_concat.ms/',caltable='Bullet_Cluster_concat.G1', field='',\
#	spw='', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=True, \
#	gaintable=['Bullet_Cluster_concat.G0'], interp=['linear'], parang=True)

#os.makedirs('Calibration_Files/selfclean/G1/Amp', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G0/93-277/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G0/504-738/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G0/892-932/Phase', exist_ok = True)


plotms(vis='BC_Target_Block_I.G0',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G0/93-277/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_II.G0',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G0/504-738/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_III.G0',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G0/892-932/Phase/phase', expformat='png', exprange='all', overwrite=True)


applycal(vis='BC_Target_Block_I.ms', field='',\
	gaintable=['BC_Target_Block_I.G0'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_II.ms', field='',\
	gaintable=['BC_Target_Block_II.G0'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_III.ms', field='',\
	gaintable=['BC_Target_Block_III.G0'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)


wsclean = "wsclean -reorder -name Bullet_Cluster_selfcal2 -size 8192 8192 -scale 1.1arcsec -data-column CORRECTED_DATA -weight briggs -0.5 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -join-channels -gap-channel-division -channels-out 6 -pol I -multiscale -parallel-deconvolution 1024 -use-wgridder -mem 50 BC_Target_Block_I.ms BC_Target_Block_II.ms BC_Target_Block_III.ms"

singularity = "singularity shell -B /tmp,/dev/shm,/disks/paradata,/data1,/data2,/net/rijn,/net/rijn9/,/software /net/achterrijn/data1/sweijen/software/containers/lofar_sksp_rijnX.sif --noprofile --norc"

subprocess.call(singularity+" "+wsclean, shell = True)


# Second round of Selfcal 
gaincal(vis='BC_Target_Block_I.ms',caltable='BC_Target_Block_I.G1', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_II.ms',caltable='BC_Target_Block_II.G1', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_III.ms',caltable='BC_Target_Block_III.G1', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

#gaincal(vis='Bullet_Cluster_concat.ms/',caltable='Bullet_Cluster_concat.G1', field='',\
#	spw='', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=True, \
#	gaintable=['Bullet_Cluster_concat.G0'], interp=['linear'], parang=True)

#os.makedirs('Calibration_Files/selfclean/G1/Amp', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G1/93-277/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G1/504-738/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G1/892-932/Phase', exist_ok = True)


plotms(vis='BC_Target_Block_I.G1',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G1/93-277/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_II.G1',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G1/504-738/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_III.G1',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G1/892-932/Phase/phase', expformat='png', exprange='all', overwrite=True)


applycal(vis='BC_Target_Block_I.ms', field='',\
	gaintable=['BC_Target_Block_I.G1'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_II.ms', field='',\
	gaintable=['BC_Target_Block_II.G1'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_III.ms', field='',\
	gaintable=['BC_Target_Block_III.G1'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)

wsclean = "wsclean -reorder -name Bullet_Cluster_selfcal3 -size 8192 8192 -scale 1.1arcsec -data-column CORRECTED_DATA -weight briggs -0.5 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -join-channels -gap-channel-division -channels-out 6 -pol I -multiscale -parallel-deconvolution 1024 -use-wgridder -mem 50 BC_Target_Block_I.ms BC_Target_Block_II.ms BC_Target_Block_III.ms"

singularity = "singularity shell -B /tmp,/dev/shm,/disks/paradata,/data1,/data2,/net/rijn,/net/rijn9/,/software /net/achterrijn/data1/sweijen/software/containers/lofar_sksp_rijnX.sif --noprofile --norc"

subprocess.call(singularity+" "+wsclean, shell = True)


# Third round of selfcal 
gaincal(vis='BC_Target_Block_I.ms',caltable='BC_Target_Block_I.G2', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_II.ms',caltable='BC_Target_Block_II.G2', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_III.ms',caltable='BC_Target_Block_III.G2', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

#gaincal(vis='Bullet_Cluster_concat.ms/',caltable='Bullet_Cluster_concat.G2', field='',\
#	spw='', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=True, \
#	gaintable=['Bullet_Cluster_concat.G0'], interp=['linear'], parang=True)

#os.makedirs('Calibration_Files/selfclean/G2/Amp', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G2/93-277/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G2/504-738/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G2/892-932/Phase', exist_ok = True)


plotms(vis='BC_Target_Block_I.G2',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G2/93-277/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_II.G2',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G2/504-738/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_III.G2',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G2/892-932/Phase/phase', expformat='png', exprange='all', overwrite=True)


applycal(vis='BC_Target_Block_I.ms', field='',\
	gaintable=['BC_Target_Block_I.G2'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_II.ms', field='',\
	gaintable=['BC_Target_Block_II.G2'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_III.ms', field='',\
	gaintable=['BC_Target_Block_III.G2'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)


wsclean = "wsclean -reorder -name Bullet_Cluster_selfcal4 -size 8192 8192 -scale 1.1arcsec -data-column CORRECTED_DATA -weight briggs -0.5 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -join-channels -gap-channel-division -channels-out 6 -pol I -multiscale -parallel-deconvolution 1024 -use-wgridder -mem 50 BC_Target_Block_I.ms BC_Target_Block_II.ms BC_Target_Block_III.ms"

singularity = "singularity shell -B /tmp,/dev/shm,/disks/paradata,/data1,/data2,/net/rijn,/net/rijn9/,/software /net/achterrijn/data1/sweijen/software/containers/lofar_sksp_rijnX.sif --noprofile --norc"

subprocess.call(singularity+" "+wsclean, shell = True)

# fourth round of selfcal 
gaincal(vis='BC_Target_Block_I.ms',caltable='BC_Target_Block_I.G3', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False,\
	gaintable=['BC_Target_Block_I.G2'], interp=['linear'],parang=True)

gaincal(vis='BC_Target_Block_II.ms',caltable='BC_Target_Block_II.G3', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, \
	gaintable=['BC_Target_Block_II.G2'], interp=['linear'], parang=True)

gaincal(vis='BC_Target_Block_III.ms',caltable='BC_Target_Block_III.G3', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False,\
	gaintable=['BC_Target_Block_III.G2'], interp=['linear'],parang=True)

#gaincal(vis='Bullet_Cluster_concat.ms/',caltable='Bullet_Cluster_concat.G3', field='',\
#	spw='', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=True, \
#	gaintable=['Bullet_Cluster_concat.G0'], interp=['linear'], parang=True)

#os.makedirs('Calibration_Files/selfclean/G3/Amp', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G3/93-277/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G3/504-738/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G3/892-932/Phase', exist_ok = True)


plotms(vis='BC_Target_Block_I.G3',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G3/93-277/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_II.G3',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G3/504-738/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_III.G3',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G3/892-932/Phase/phase', expformat='png', exprange='all', overwrite=True)


applycal(vis='BC_Target_Block_I.ms', field='',\
	gaintable=['BC_Target_Block_I.G3'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_II.ms', field='',\
	gaintable=['BC_Target_Block_II.G3'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_III.ms', field='',\
	gaintable=['BC_Target_Block_III.G3'],\
                  gainfield=[''],
                  interp=['linear'], calwt=[False], parang=True)


wsclean = "wsclean -reorder -name Bullet_Cluster_selfcal5 -size 8192 8192 -scale 1.1arcsec -data-column CORRECTED_DATA -weight briggs -0.5 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -join-channels -gap-channel-division -channels-out 6 -pol I -multiscale -parallel-deconvolution 1024 -use-wgridder -mem 50 BC_Target_Block_I.ms BC_Target_Block_II.ms BC_Target_Block_III.ms"

singularity = "singularity shell -B /tmp,/dev/shm,/disks/paradata,/data1,/data2,/net/rijn,/net/rijn9/,/software /net/achterrijn/data1/sweijen/software/containers/lofar_sksp_rijnX.sif --noprofile --norc"

subprocess.call(singularity+" "+wsclean, shell = True)

# fifth round of selfcal with amplitude
gaincal(vis='BC_Target_Block_I.ms',caltable='BC_Target_Block_I_p.G4', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_II.ms',caltable='BC_Target_Block_II_p.G4', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_III.ms',caltable='BC_Target_Block_III_p.G4', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_I.ms',caltable='BC_Target_Block_I_ap.G4', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=False,\
	gaintable=['BC_Target_Block_I_p.G4'], interp=['linear'],parang=True)

gaincal(vis='BC_Target_Block_II.ms',caltable='BC_Target_Block_II_ap.G4', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=False, \
	gaintable=['BC_Target_Block_II_p.G4'], interp=['linear'], parang=True)

gaincal(vis='BC_Target_Block_III.ms',caltable='BC_Target_Block_III_ap.G4', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=False,\
	gaintable=['BC_Target_Block_III_p.G4'], interp=['linear'],parang=True)


os.makedirs('Calibration_Files/selfclean/G4/93-277/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G4/504-738/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G4/892-932/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G4/93-277/Amp', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G4/504-738/Amp', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G4//Amp', exist_ok = True)

plotms(vis='BC_Target_Block_I_ap.G4',xaxis='time',yaxis='amp',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,0,0.3],\
	showgui=False, plotfile='Calibration_Files/selfclean/G4/93-277/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_II_ap.G4',xaxis='time',yaxis='amp',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,0,0.3],\
	showgui=False, plotfile='Calibration_Files/selfclean/G4/504-738/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_III_ap.G4',xaxis='time',yaxis='amp',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,0,0.3],\
	showgui=False, plotfile='Calibration_Files/selfclean/G4/892-932/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_I_p.G4',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G4/93-277/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_II_p.G4',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G4/504-738/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_III_p.G4',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G4/892-932/Phase/phase', expformat='png', exprange='all', overwrite=True)


applycal(vis='BC_Target_Block_I.ms', field='',\
	gaintable=['BC_Target_Block_I_p.G4',
                  'BC_Target_Block_I_ap.G4'],\
                  gainfield=['',''],
                  interp=['linear','linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_II.ms', field='',\
	gaintable=['BC_Target_Block_II_p.G4',
                  'BC_Target_Block_II_ap.G4'],\
                  gainfield=['',''],
                  interp=['linear','linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_III.ms', field='',\
	gaintable=['BC_Target_Block_III_p.G4',
                  'BC_Target_Block_III_ap.G4'],\
                  gainfield=['',''],
                  interp=['linear','linear'], calwt=[False], parang=True)


wsclean = "wsclean -reorder -name Bullet_Cluster_selfcal6 -size 8192 8192 -scale 1.1arcsec -data-column CORRECTED_DATA -weight briggs -0.5 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -join-channels -gap-channel-division -channels-out 6 -pol I -multiscale -parallel-deconvolution 1024 -use-wgridder -mem 50 BC_Target_Block_I.ms BC_Target_Block_II.ms BC_Target_Block_III.ms"

singularity = "singularity shell -B /tmp,/dev/shm,/disks/paradata,/data1,/data2,/net/rijn,/net/rijn9/,/software /net/achterrijn/data1/sweijen/software/containers/lofar_sksp_rijnX.sif --noprofile --norc"

subprocess.call(singularity+" "+wsclean, shell = True)

# sixth round of selfcal with amplitude
gaincal(vis='BC_Target_Block_I.ms',caltable='BC_Target_Block_I_p.G5', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_II.ms',caltable='BC_Target_Block_II_p.G5', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_III.ms',caltable='BC_Target_Block_III_p.G5', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='p',solnorm=False, parang=True)

gaincal(vis='BC_Target_Block_I.ms',caltable='BC_Target_Block_I_ap.G5', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=False,\
	gaintable=['BC_Target_Block_I_p.G5'], interp=['linear'],parang=True)

gaincal(vis='BC_Target_Block_II.ms',caltable='BC_Target_Block_II_ap.G5', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=False, \
	gaintable=['BC_Target_Block_II_p.G5'], interp=['linear'], parang=True)

gaincal(vis='BC_Target_Block_III.ms',caltable='BC_Target_Block_III_ap.G5', field='',\
	spw='', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=False,\
	gaintable=['BC_Target_Block_III_p.G5'], interp=['linear'],parang=True)


os.makedirs('Calibration_Files/selfclean/G5/93-277/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G5/504-738/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G5/892-932/Phase', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G5/93-277/Amp', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G5/504-738/Amp', exist_ok = True)
os.makedirs('Calibration_Files/selfclean/G5//Amp', exist_ok = True)

plotms(vis='BC_Target_Block_I_ap.G5',xaxis='time',yaxis='amp',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,0,0.3],\
	showgui=False, plotfile='Calibration_Files/selfclean/G5/93-277/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_II_ap.G5',xaxis='time',yaxis='amp',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,0,0.3],\
	showgui=False, plotfile='Calibration_Files/selfclean/G5/504-738/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_III_ap.G5',xaxis='time',yaxis='amp',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,0,0.3],\
	showgui=False, plotfile='Calibration_Files/selfclean/G5/892-932/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_I_p.G5',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G5/93-277/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_II_p.G5',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G5/504-738/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='BC_Target_Block_III_p.G5',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/selfclean/G5/892-932/Phase/phase', expformat='png', exprange='all', overwrite=True)


applycal(vis='BC_Target_Block_I.ms', field='',\
	gaintable=['BC_Target_Block_I_p.G5',
                  'BC_Target_Block_I_ap.G5'],\
                  gainfield=['',''],
                  interp=['linear','linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_II.ms', field='',\
	gaintable=['BC_Target_Block_II_p.G5',
                  'BC_Target_Block_II_ap.G5'],\
                  gainfield=['',''],
                  interp=['linear','linear'], calwt=[False], parang=True)

applycal(vis='BC_Target_Block_III.ms', field='',\
	gaintable=['BC_Target_Block_III_p.G5',
                  'BC_Target_Block_III_ap.G5'],\
                  gainfield=['',''],
                  interp=['linear','linear'], calwt=[False], parang=True)


wsclean = "wsclean -reorder -name Bullet_Cluster_selfcal7 -size 8192 8192 -scale 1.1arcsec -data-column CORRECTED_DATA -weight briggs -0.5 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -join-channels -gap-channel-division -channels-out 6 -pol I -multiscale -parallel-deconvolution 1024 -use-wgridder -mem 50 BC_Target_Block_I.ms BC_Target_Block_II.ms BC_Target_Block_III.ms"

singularity = "singularity shell -B /tmp,/dev/shm,/disks/paradata,/data1,/data2,/net/rijn,/net/rijn9/,/software /net/achterrijn/data1/sweijen/software/containers/lofar_sksp_rijnX.sif --noprofile --norc"

subprocess.call(singularity+" "+wsclean, shell = True)
