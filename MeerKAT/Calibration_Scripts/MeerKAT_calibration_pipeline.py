"""
Author: Joppe Swart
Description: Script to calibrate the data
Last modified: November 2022
"""

import os

#measuramentset = 'Bullet_CLuster_Hanningsmooth.ms'


# Make a directory for the figures etc
os.makedirs('Calibration_Files', exist_ok = True)

# Make a listobs
listobs(vis='Bullet_CLuster_Hanningsmooth.ms', listfile='Calibration_Files/Bullet_listobs.txt', overwrite=True) 

## Manual flagging of the known rfi bands in the MeerKat L band
## GSM Up and Down link
#flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:0~119;167~335')

## GPS
#flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:167~335; 3230~3325; 1564~1660; 2320~2378; 1320~1416')

## GLONASS
#flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:3359~3445; 1684~1718; 1493')

## GALILEO
#flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:3279; 1270; 1517; 1444; 1860')

## IRIDIUM and INMARSAT
#flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:3474~3521;3043~3177')

# Setjy on the flux calibrator
#execfile('setjy_manual.py')

#execfile('setjy_polarisation_coefficients.py')

# Initial phase calibration
gaincal(vis='Bullet_CLuster_Hanningsmooth.ms/', caltable='Bullet_CLuster_Hanningsmooth.G0/', field='0408-658', \
	refant='m014', spw='0:2000~2003', calmode='p', solint='10s', minsnr=5)

os.makedirs('Calibration_Files/Phase_cal/G0', exist_ok=True)

plotms(vis='Bullet_CLuster_Hanningsmooth.G0',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G0/phase', expformat='png', exprange='all', overwrite=True)

# Delay calibration
gaincal(vis='Bullet_CLuster_Hanningsmooth.ms/',caltable='Bullet_CLuster_Hanningsmooth.K0', field='0408-658',\
	refant='m014',spw='0:5~3717',gaintype='K',  solint='inf',combine='scan',minsnr=5, \
   	gaintable=['Bullet_CLuster_Hanningsmooth.G0'])

os.makedirs('Calibration_Files/Delay', exist_ok = True)

plotms(vis='Bullet_CLuster_Hanningsmooth.K0/',xaxis='antenna1',yaxis='delay',coloraxis='baseline', \
	showgui=False, plotfile='Calibration_Files/Delay/delay', expformat='png', exprange='all', overwrite=True)

# Bandpass calibration
bandpass(vis='Bullet_CLuster_Hanningsmooth.ms/',caltable='Bullet_CLuster_Hanningsmooth.B0', field='0408-658',\
	spw='',refant='m014',combine='scan',solint='inf',bandtype='B',\
	gaintable=['Bullet_CLuster_Hanningsmooth.G0','Bullet_CLuster_Hanningsmooth.K0'])


os.makedirs('Calibration_Files/Bandpass/B0/Amp', exist_ok = True)
os.makedirs('Calibration_Files/Bandpass/B0/Phase', exist_ok = True)


plotms(vis='Bullet_CLuster_Hanningsmooth.B0/',field='0408-658', xaxis='chan',yaxis='amp',coloraxis='corr', \
	plotrange=[-1,-1,0,3],iteraxis='antenna',\
	showgui=False, plotfile='Calibration_Files/Bandpass/B0/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_CLuster_Hanningsmooth.B0/',field='0408-658', xaxis='chan',yaxis='phase',coloraxis='corr',\
	plotrange=[-1,-1,-180,180], iteraxis='antenna',\
	showgui=False, plotfile='Calibration_Files/Bandpass/B0/Phase/phase', expformat='png', exprange='all', overwrite=True)

# Phase calibration 
gaincal(vis='Bullet_CLuster_Hanningsmooth.ms/',caltable='Bullet_CLuster_Hanningsmooth.G1', field='0,1,3,4',\
	spw='0:5~3720', solint='10s',refant='m014',gaintype='G',calmode='ap',solnorm=False, \
	gaintable=['Bullet_CLuster_Hanningsmooth.K0/','Bullet_CLuster_Hanningsmooth.B0/'], interp=['linear','nearest'])

os.makedirs('Calibration_Files/Phase_cal/G1/Amp', exist_ok = True)
os.makedirs('Calibration_Files/Phase_cal/G1/Phase', exist_ok = True)
os.makedirs('Calibration_Files/Phase_cal/G1/Diff', exist_ok = True)

plotms(vis='Bullet_CLuster_Hanningsmooth.G1',xaxis='time',yaxis='phase',\
	iteraxis='corr',coloraxis='baseline',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G1/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_CLuster_Hanningsmooth.G1',xaxis='time',yaxis='amp', iteraxis='corr',coloraxis='baseline',\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G1/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_CLuster_Hanningsmooth.G1', xaxis='time', yaxis='phase', correlation='/', coloraxis='baseline', plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G1/Diff', expformat='png', exprange='all', overwrite=True)

# Kross hand delays
gaincal(vis='Bullet_CLuster_Hanningsmooth.ms', caltable='Bullet_CLuster_Hanningsmooth.Kcross',\
        field='3C286', spw='0:5~3720',\
        gaintype='KCROSS', solint='inf', combine='scan', refant='m014',\
        gaintable=['Bullet_CLuster_Hanningsmooth.K0',\
                   'Bullet_CLuster_Hanningsmooth.B0',\
                   'Bullet_CLuster_Hanningsmooth.G1'],\
        gainfield=['','',''], parang=True)

os.makedirs('Calibration_Files/Pol_Cal/Delay', exist_ok = True)

plotms(vis='Bullet_CLuster_Hanningsmooth.Kcross',xaxis='antenna',yaxis='delay',\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D1/Delay/delay', expformat='png', exprange='all', overwrite=True)

# Solve the leakage terms
polcal(vis='Bullet_CLuster_Hanningsmooth.ms',caltable='Bullet_CLuster_Hanningsmooth.D1',\
       field='0823-500',spw='0:5~3720', refant='m014',poltype='Df+QU',solint='inf',combine='scan',\
       gaintable=['Bullet_CLuster_Hanningsmooth.K0',\
                  'Bullet_CLuster_Hanningsmooth.B0',\
                  'Bullet_CLuster_Hanningsmooth.G1',\
                  'Bullet_CLuster_Hanningsmooth.Kcross'],\
                  gainfield=['','','',''])

os.makedirs('Calibration_Files/Pol_Cal/D1/Amp', exist_ok = True)
os.makedirs('Calibration_Files/Pol_Cal/D1/Phase', exist_ok = True)


plotms(vis='Bullet_CLuster_Hanningsmooth.D1',xaxis='chan',yaxis='amp', iteraxis='antenna',coloraxis='corr', field='0647-475',\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D1/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_CLuster_Hanningsmooth.D1',xaxis='chan',yaxis='phase', iteraxis='antenna',coloraxis='corr',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D1/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_CLuster_Hanningsmooth.D1',xaxis='antenna1',yaxis='amp',coloraxis='corr', field='0647-475',\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/D1/antenna', expformat='png', exprange='all', overwrite=True)

# Solve for the pol angle
polcal(vis='Bullet_CLuster_Hanningsmooth.ms',caltable='Bullet_CLuster_Hanningsmooth.X1',\
       field='3C286',combine='scan', refant='m014',poltype='Xf',solint='inf',\
       gaintable=['Bullet_CLuster_Hanningsmooth.K0',\
                  'Bullet_CLuster_Hanningsmooth.B0',\
                  'Bullet_CLuster_Hanningsmooth.G1',\
                  'Bullet_CLuster_Hanningsmooth.Kcross',\
                  'Bullet_CLuster_Hanningsmooth.D1'],\
       gainfield=['','','','',''])

os.makedirs('Calibration_Files/Pol_Cal/X1', exist_ok = True)

plotms(vis='Bullet_CLuster_Hanningsmooth.X1',xaxis='chan',yaxis='phase',\
	showgui=False, plotfile='Calibration_Files/Pol_Cal/X1/phase', expformat='png', exprange='all', overwrite=True)

# Another phase cal
# Phase calibration 
gaincal(vis='Bullet_CLuster_Hanningsmooth.ms/',caltable='Bullet_CLuster_Hanningsmooth.G2', field='0,1,3,4',\
	spw='0:5~3720', solint='60s',refant='m014',gaintype='G',calmode='ap',solnorm=False, \
	gaintable=['Bullet_CLuster_Hanningsmooth.K0',\
                  'Bullet_CLuster_Hanningsmooth.B0',\
                  'Bullet_CLuster_Hanningsmooth.Kcross',\
                  'Bullet_CLuster_Hanningsmooth.D1',\
                  'Bullet_CLuster_Hanningsmooth.X1'],\
                  interp=['linear','nearest','','',''])

os.makedirs('Calibration_Files/Phase_cal/G2/Amp', exist_ok = True)
os.makedirs('Calibration_Files/Phase_cal/G2/Phase', exist_ok = True)
os.makedirs('Calibration_Files/Phase_cal/G2/Diff', exist_ok = True)

plotms(vis='Bullet_CLuster_Hanningsmooth.G2',xaxis='time',yaxis='phase',\
	iteraxis='corr',coloraxis='baseline',plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G2/Phase/phase', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_CLuster_Hanningsmooth.G2',xaxis='time',yaxis='amp', iteraxis='corr',coloraxis='baseline',\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G2/Amp/amp', expformat='png', exprange='all', overwrite=True)

plotms(vis='Bullet_CLuster_Hanningsmooth.G2', xaxis='time', yaxis='phase', correlation='/', coloraxis='baseline', plotrange=[-1,-1,-180,180],\
	showgui=False, plotfile='Calibration_Files/Phase_cal/G2/Diff', expformat='png', exprange='all', overwrite=True)


# For apply cal use parang = True
