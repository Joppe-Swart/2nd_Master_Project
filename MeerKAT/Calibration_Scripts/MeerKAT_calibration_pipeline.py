
"""
Author: Joppe Swart
Description: Script to calibrate the data
Last modified: November 2022
"""

import os

hanningsmooth(vis='1529816457_sdp_l0.ms', outputvis='Bullet_Cluster_HS.ms')

measuramentset = 'Bullet_Cluster_HS.ms'


# Make a directory for the figures etc
os.makedirs('Calibration_Files', exist_ok = True)

# Make a listobs
listobs(vis='Bullet_Cluster_HS.ms', listfile='Calibration_Files/Bullet_listobs.txt', overwrite=True)

# Quack flagging
flagdata(vis='Bullet_Cluster_HS.ms/', mode='quack', quackinterval=10.0, quackmode='beg')

# Manual flagging of the known rfi bands in the MeerKat L band

# GPS
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:167~335; 3230~3325; 1564~1660; 2320~2378; 1320~1416')

# GLONASS
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:3359~3445; 1684~1718; 1493', uvrange='>1km')

# GALILEO
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:3279; 1260~1300; 1517; 1444; 1860; 1770~1960', uvrange='>1km')

# IRIDIUM and INMARSAT
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:3474~3521;3043~3177', uvrange='>1km')

# For all the data larger blocks of channels need to be flagged
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, timerange='10:57:00~11:02:00')
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, timerange='07:14:00~07:18:00')
#Block I
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:0~335')
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:655~683')
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:1009~1023')
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:712~722')
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:900~1000')

#Block II
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:1172~1985')
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, timerange='08:30:00~08:33:00', spw='0:2248~2296')
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, timerange='07:14:00~07:24:00', spw='0:1985~2033')
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, timerange='14:46:00~14:57:00', spw='0:1985~2033')
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, timerange='14:21:00~14:40:00', spw='0:2961~2985')
#Block III
flagdata(vis='Bullet_Cluster_HS.ms/', mode='manual', flagbackup=True, spw='0:2990~3526')

# Setjy on the flux calibrators
execfile('setjy_0408-658.py')

execfile('setjy_3C286.py')

# Initial phase calibration
gaincal(vis='Bullet_Cluster_HS.ms/', caltable='Bullet_Cluster.G0', field='0408-658', \
	refant='m014', spw='0:2500~2510', calmode='p', solint='10s', minsnr=5, parang=True)

os.makedirs('Calibration_Files/Phase_cal/G0', exist_ok=True)

plotms(vis='Bullet_Cluster.G0',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G0/phase', expformat='png', exprange='all', overwrite=True)

# Delay calibration
gaincal(vis='Bullet_Cluster_HS.ms/',caltable='Bullet_Cluster.K0', field='0408-658',\
	refant='m014',spw='0:5~3717',gaintype='K',  solint='inf',combine='scan',minsnr=5, \
   	gaintable=['Bullet_Cluster.G0'])

os.makedirs('Calibration_Files/Delay', exist_ok = True)

plotms(vis='Bullet_Cluster.K0/',xaxis='antenna1',yaxis='delay',coloraxis='baseline', \
	showgui=False, plotfile='Calibration_Files/Delay/delay', expformat='png', exprange='all', overwrite=True)

# Bandpass calibration
bandpass(vis='Bullet_Cluster_HS.ms/',caltable='Bullet_Cluster.B0', field='0408-658',\
	spw='',refant='m014',combine='scan',solint='inf',bandtype='B',\
	gaintable=['Bullet_Cluster.G0','Bullet_Cluster.K0'], parang=True)


os.makedirs('Calibration_Files/Bandpass/B0/Amp', exist_ok = True)
os.makedirs('Calibration_Files/Bandpass/B0/Phase', exist_ok = True)


plotms(vis='Bullet_Cluster.B0/',field='0408-658', xaxis='chan',yaxis='amp',coloraxis='corr', \
	plotrange=[-1,-1,0,3],iteraxis='antenna',\
	showgui=False, plotfile='Calibration_Files/Bandpass/B0/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.B0/',field='0408-658', xaxis='chan',yaxis='phase',coloraxis='corr',\
	plotrange=[-1,-1,-180,180], iteraxis='antenna',\
	showgui=False, plotfile='Calibration_Files/Bandpass/B0/Phase/phase', expformat='png', exprange='all', overwrite=True)

# Phase calibration 
gaincal(vis='Bullet_Cluster_HS.ms/',caltable='Bullet_Cluster.G1', field='0,1,3,4',\
	spw='0:5~3720', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=False, \
	gaintable=['Bullet_Cluster.K0/','Bullet_Cluster.B0/'], interp=['linear','nearest'], parang=True)

os.makedirs('Calibration_Files/Phase_cal/G1/Amp', exist_ok = True)
os.makedirs('Calibration_Files/Phase_cal/G1/Phase', exist_ok = True)
os.makedirs('Calibration_Files/Phase_cal/G1/Diff', exist_ok = True)

plotms(vis='Bullet_Cluster.G1',xaxis='time',yaxis='amp',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G1/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.G1',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G1/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.G1',xaxis='time',yaxis='phase',\
	iteraxis='corr',coloraxis='baseline',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G1/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.G1',xaxis='time',yaxis='amp', iteraxis='corr',coloraxis='baseline',\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G1/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.G1', xaxis='time', yaxis='phase', correlation='/', coloraxis='baseline', plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G1/Diff/diff', expformat='png', exprange='all', overwrite=True)

 Kross hand delays
gaincal(vis='Bullet_Cluster_HS.ms', caltable='Bullet_Cluster.Kcross',\
        field='3C286', spw='0:5~3720',\
        gaintype='KCROSS', solint='inf', combine='scan', refant='m014',\
        gaintable=['Bullet_Cluster.K0',\
                   'Bullet_Cluster.B0',\
                   'Bullet_Cluster.G1'],\
        gainfield=['','',''], parang=True)

os.makedirs('Calibration_Files/Pol_Cal/Delay', exist_ok = True)

plotms(vis='Bullet_Cluster.Kcross',xaxis='antenna1',yaxis='delay',coloraxis='baseline', \
	showgui=False, plotfile='Calibration_Files/Pol_Cal/Delay/delay', expformat='png', exprange='all', overwrite=True)

# Solve the leakage terms
polcal(vis='Bullet_Cluster_HS.ms',caltable='Bullet_Cluster.D0823-500',\
       field='0823-500',spw='0:5~3720', refant='m014',poltype='Df+QU',solint='inf',combine='scan',\
       gaintable=['Bullet_Cluster.K0',\
                  'Bullet_Cluster.B0',\
                  'Bullet_Cluster.G1',\
                  'Bullet_Cluster.Kcross'],\
                  gainfield=['','','',''])

polcal(vis='Bullet_Cluster_HS.ms',caltable='Bullet_Cluster.D0647-475',\
       field='0647-475',spw='0:5~3720', refant='m014',poltype='Df+QU',solint='inf',combine='scan',\
       gaintable=['Bullet_Cluster.K0',\
                  'Bullet_Cluster.B0',\
                  'Bullet_Cluster.G1',\
                  'Bullet_Cluster.Kcross'],\
                  gainfield=['','','',''])

polcal(vis='Bullet_Cluster_HS.ms',caltable='Bullet_Cluster.D0408-658',\
       field='0408-658',spw='0:5~3720', refant='m014',poltype='Df',solint='inf',combine='scan',\
       gaintable=['Bullet_Cluster.K0',\
                  'Bullet_Cluster.B0',\
                  'Bullet_Cluster.G1',\
                  'Bullet_Cluster.Kcross'],\
                  gainfield=['','','',''])

os.makedirs('Calibration_Files/Pol_Cal/D0823-500/Amp', exist_ok = True)
os.makedirs('Calibration_Files/Pol_Cal/D0823-500/Phase', exist_ok = True)
os.makedirs('Calibration_Files/Pol_Cal/D0823-500/Antenna', exist_ok = True)

os.makedirs('Calibration_Files/Pol_Cal/D0647-475/Amp', exist_ok = True)
os.makedirs('Calibration_Files/Pol_Cal/D0647-475/Phase', exist_ok = True)
os.makedirs('Calibration_Files/Pol_Cal/D0647-475/Antenna', exist_ok = True)

os.makedirs('Calibration_Files/Pol_Cal/D0408-658/Amp', exist_ok = True)
os.makedirs('Calibration_Files/Pol_Cal/D0408-658/Phase', exist_ok = True)
os.makedirs('Calibration_Files/Pol_Cal/D0408-658/Antenna', exist_ok = True)


plotms(vis='Bullet_Cluster.D0823-500',xaxis='chan',yaxis='amp', iteraxis='antenna',coloraxis='corr',plotrange=[-1,-1,0,0.05],\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D0823-500/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.D0823-500',xaxis='chan',yaxis='phase', iteraxis='antenna',coloraxis='corr',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D0823-500/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.D0823-500',xaxis='antenna1',yaxis='amp',\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D0823-500/Antenna/antenna', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.D0647-475',xaxis='chan',yaxis='amp', iteraxis='antenna',coloraxis='corr',plotrange=[-1,-1,0,0.05],\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D0647-475/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.D0647-475',xaxis='chan',yaxis='phase', iteraxis='antenna',coloraxis='corr',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D0647-475/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.D0647-475',xaxis='antenna1',yaxis='amp',\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D0647-475/Antenna/antenna', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.D0408-658',xaxis='chan',yaxis='amp', iteraxis='antenna',coloraxis='corr',plotrange=[-1,-1,0,0.05],\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D0408-658/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.D0408-658',xaxis='chan',yaxis='phase', iteraxis='antenna',coloraxis='corr',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D0408-658/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.D0408-658',xaxis='antenna1',yaxis='amp',\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D0408-658/Antenna/antenna', expformat='png', exprange='all', overwrite=True)

# Solve for the pol angle
polcal(vis='Bullet_Cluster_HS.ms',caltable='Bullet_Cluster.X1',\
       field='3C286',combine='scan', refant='m014',poltype='Xf',solint='inf',\
       gaintable=['Bullet_Cluster.K0',\
                  'Bullet_Cluster.B0',\
                  'Bullet_Cluster.G1',\
                  'Bullet_Cluster.Kcross',\
                  'Bullet_Cluster.D0408-658'],\
       gainfield=['','','','',''])

os.makedirs('Calibration_Files/Pol_Cal/X1', exist_ok = True)

plotms(vis='Bullet_Cluster.X1',xaxis='chan',yaxis='phase', plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/X1/phase', expformat='png', exprange='all', overwrite=True)

# Another phase cal
# Phase calibration


gaincal(vis='Bullet_Cluster_HS.ms/',caltable='Bullet_Cluster.G2', field='0,1,3,4',\
	spw='0:5~3720', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=False, \
	gaintable=['Bullet_Cluster.K0',\
                  'Bullet_Cluster.B0',\
                  'Bullet_Cluster.Kcross',\
                  'Bullet_Cluster.D0408-658',\
                  'Bullet_Cluster.X1'],\
                  interp=['linear','nearest','','',''], parang=True)

os.makedirs('Calibration_Files/Phase_cal/G2/Amp', exist_ok = True)
os.makedirs('Calibration_Files/Phase_cal/G2/Phase', exist_ok = True)
os.makedirs('Calibration_Files/Phase_cal/G2/Diff', exist_ok = True)

plotms(vis='Bullet_Cluster.G2',xaxis='time',yaxis='amp',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,0,0.3],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G2/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.G2',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G2/Phase/phase', expformat='png', exprange='all', overwrite=True)


plotms(vis='Bullet_Cluster.G2',xaxis='time',yaxis='phase',\
	iteraxis='corr',coloraxis='baseline',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G2/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.G2',xaxis='time',yaxis='amp', iteraxis='corr',coloraxis='baseline',\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G2/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster.G2', xaxis='time', yaxis='phase', correlation='/', coloraxis='baseline', plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G2/Diff/diff', expformat='png', exprange='all', overwrite=True)

smoothcal(vis='Bullet_Cluster_HS.ms', tablein='Bullet_Cluster.G2/', caltable='BC_Smoothcal_60S.G2', smoothtime=60) 
smoothcal(vis='Bullet_Cluster_HS.ms', tablein='Bullet_Cluster.G2/', caltable='BC_Smoothcal_600S.G2', smoothtime=600)
smoothcal(vis='Bullet_Cluster_HS.ms', tablein='Bullet_Cluster.G2/', caltable='BC_Smoothcal_3600S.G2', smoothtime=3600)

flagmanager(vis='Bullet_Cluster_HS.ms/', mode='restore', versionname='applycal_1')

smoothcal(vis='Bullet_Cluster_HS.ms', tablein='Bullet_Cluster.G2/', caltable='BC_Smoothcal_8100S.G2', smoothtime=8100)

#Set the flux scale for all calibrators
myscale = fluxscale(vis='Bullet_Cluster_HS.ms/',\
                    caltable='BC_Smoothcal_8100S.G2', \
                    fluxtable='Bullet_Cluster.fluxscale1', \
                    reference=['0408-658'],\
                    transfer=['0647-475,0823-500'],\
                    incremental=False)

os.makedirs('Calibration_Files/Fluxscale', exist_ok = True)

plotms(vis='Bullet_Cluster.fluxscale1',xaxis='time',yaxis='amp',\
       correlation='X',coloraxis='baseline',\
       showgui=False, plotfile='Calibration_Files/Fluxscale/X', expformat='png', exprange='all', overwrite=True)
plotms(vis='Bullet_Cluster.fluxscale1',xaxis='time',yaxis='amp',\
       correlation='Y',coloraxis='baseline',\
       showgui=False, plotfile='Calibration_Files/Fluxscale/Y', expformat='png', exprange='all', overwrite=True)

## For apply cal use parang = True

applycal(vis='Bullet_Cluster_HS.ms/', field='3C286',\
	gaintable=['Bullet_Cluster.fluxscale1',\
                  'Bullet_Cluster.K0',\
                  'Bullet_Cluster.B0',\
                  'Bullet_Cluster.Kcross',\
                  'Bullet_Cluster.D0408-658',\
                  'Bullet_Cluster.X1'],\
                  gainfield=['','','','','',''],
                  interp=['linear','nearest','nearest','nearest','nearest','nearest'], calwt=[False], parang=True)

applycal(vis='Bullet_Cluster_HS.ms/', field='0408-658',\
	gaintable=['Bullet_Cluster.fluxscale1',\
                  'Bullet_Cluster.K0',\
                  'Bullet_Cluster.B0',\
                  'Bullet_Cluster.Kcross',\
                  'Bullet_Cluster.D0408-658',\
                  'Bullet_Cluster.X1'],\
                  gainfield=['','','','','',''],
                  interp=['linear','nearest','nearest','nearest','nearest','nearest'], calwt=[False], parang=True)

applycal(vis='Bullet_Cluster_HS.ms/', field='0647-475',\
	gaintable=['Bullet_Cluster.fluxscale1',\
                  'Bullet_Cluster.K0',\
                  'Bullet_Cluster.B0',\
                  'Bullet_Cluster.Kcross',\
                  'Bullet_Cluster.D0408-658',\
                  'Bullet_Cluster.X1'],\
                  gainfield=['','','','','',''],
                  interp=['linear','nearest','nearest','nearest','nearest','nearest'], calwt=[False], parang=True)

applycal(vis='Bullet_Cluster_HS.ms/', field='0823-500',\
	gaintable=['Bullet_Cluster.fluxscale1',\
                  'Bullet_Cluster.K0',\
                  'Bullet_Cluster.B0',\
                  'Bullet_Cluster.Kcross',\
                  'Bullet_Cluster.D0408-658',\
                  'Bullet_Cluster.X1'],\
                  gainfield=['','','','','',''],
                  interp=['linear','nearest','nearest','nearest','nearest','nearest'], calwt=[False], parang=True)
flagmanager(vis='Bullet_Cluster_HS.ms/', mode='restore', versionname='applycal_2')

applycal(vis='Bullet_Cluster_HS.ms/', field='',\
	gaintable=['Bullet_Cluster.fluxscale1',\
                  'Bullet_Cluster.K0',\
                  'Bullet_Cluster.B0',\
                  'Bullet_Cluster.Kcross',\
                  'Bullet_Cluster.D0408-658',\
                  'Bullet_Cluster.X1'],\
                  gainfield=['','','','','',''],
                  interp=['linear','nearest','nearest','nearest','nearest','nearest'], calwt=[False], parang=True)


os.makedirs('Calibration_Files/Applycal/3C286', exist_ok = True)
os.makedirs('Calibration_Files/Applycal/0408-658', exist_ok = True)
os.makedirs('Calibration_Files/Applycal/0647-475', exist_ok = True)
os.makedirs('Calibration_Files/Applycal/0823-500', exist_ok = True)
os.makedirs('Calibration_Files/Applycal/bullet', exist_ok = True)

plotms(vis='Bullet_Cluster_HS.ms',field='3C286',correlation='XX, YY',
       timerange='',antenna='',avgtime='60',
       xaxis='channel',yaxis='amp',ydatacolumn='corrected',
       coloraxis='corr',
       showgui=False, plotfile='Calibration_Files/Applycal/3C286/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster_HS.ms',field='3C286',correlation='XX, YY',
       timerange='',antenna='',avgtime='60',
       xaxis='channel',yaxis='phase',ydatacolumn='corrected',
       plotrange=[-1,-1,-180,180],coloraxis='corr',
       showgui=False, plotfile='Calibration_Files/Applycal/3C286/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster_HS.ms',field='0408-658',correlation='XX, YY',
       timerange='',antenna='',avgtime='60',
       xaxis='channel',yaxis='amp',ydatacolumn='corrected',
       coloraxis='corr',
       showgui=False, plotfile='Calibration_Files/Applycal/0408-658/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster_HS.ms',field='0408-658',correlation='XX, YY',
       timerange='',antenna='',avgtime='60',
       xaxis='channel',yaxis='phase',ydatacolumn='corrected',
       plotrange=[-1,-1,-180,180],coloraxis='corr',
       showgui=False, plotfile='Calibration_Files/Applycal/0408-658/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster_HS.ms',field='0647-475',correlation='XX, YY',
       timerange='',antenna='',avgtime='60',
       xaxis='channel',yaxis='amp',ydatacolumn='corrected',
       coloraxis='corr',
       showgui=False, plotfile='Calibration_Files/Applycal/0647-475/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster_HS.ms',field='0647-475',correlation='XX, YY',
       timerange='',antenna='',avgtime='60',
       xaxis='channel',yaxis='phase',ydatacolumn='corrected',
       plotrange=[-1,-1,-180,180],coloraxis='corr',
       showgui=False, plotfile='Calibration_Files/Applycal/0647-475/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster_HS.ms',field='0823-500',correlation='XX, YY',
       timerange='',antenna='',avgtime='60',
       xaxis='channel',yaxis='amp',ydatacolumn='corrected',
       coloraxis='corr',
       showgui=False, plotfile='Calibration_Files/Applycal/0823-500/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster_HS.ms',field='0823-500',correlation='XX, YY',
       timerange='',antenna='',avgtime='60',
       xaxis='channel',yaxis='phase',ydatacolumn='corrected',
       plotrange=[-1,-1,-180,180],coloraxis='corr',
       showgui=False, plotfile='Calibration_Files/Applycal/0823-500/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster_HS.ms',field='bullet',correlation='XX, YY',
       timerange='',antenna='',avgtime='60',
       xaxis='channel',yaxis='amp',ydatacolumn='corrected',
       coloraxis='corr',
       showgui=False, plotfile='Calibration_Files/Applycal/bullet/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_Cluster_HS.ms',field='bullet',correlation='XX, YY',
       timerange='',antenna='',avgtime='60',
       xaxis='channel',yaxis='phase',ydatacolumn='corrected',
       plotrange=[-1,-1,-180,180],coloraxis='corr',
       showgui=False, plotfile='Calibration_Files/Applycal/bullet/phase', expformat='png', exprange='all', overwrite=True)

split(vis='Bullet_Cluster_HS.ms/', outputvis='BC_Test_Block_I.ms', spw='0:336~1171', field = 'bullet', width = 4)
split(vis='Bullet_Cluster_HS.ms/', outputvis='BC_Test_Block_II.ms', spw='0:1986~2989', field = 'bullet', width = 4)
split(vis='Bullet_Cluster_HS.ms/', outputvis='BC_Test_Block_III.ms', spw='0:3516~3721', field = 'bullet', width = 4)






