"""
Casa commands for calibrating the data
Author: Joppe Swart

Tasks performed in casa	


"""
# Listobs
listobs(vis='1529816457_sdp_l0.ms/', listfile='Bullet_listobs.txt') 

# Hanning Smoothing
hanningsmooth(vis='1529816457_sdp_l0.ms', outputvis='Bullet_CLuster_Hanningsmooth.ms') 

#Some initial cleaning
tclean(vis='Bullet_CLuster_Hanningsmooth.ms/', imagename='Stokes_I_try',field='2',imsize=[512], \
	cell=['1.4arcsec'],phasecenter='J2000 06:58:32.7 -55.57.19.0',stokes='I', \
	specmode='mfs', gridder='standard',niter=10000,interactive=True,threshold='15uJy',\
	noise='1uJy',robust=0.5,cycleniter=100,parallel=True)
	
tclean(vis='Bullet_CLuster_Hanningsmooth.ms/', imagename='Stokes_U_try',field='2',imsize=[512], \
	cell=['1.4arcsec'],phasecenter='J2000 06:58:32.7 -55.57.19.0',stokes='I', \
	specmode='mfs', gridder='standard',niter=10000,interactive=True,threshold='100uJy',\
	noise='1uJy',robust=0.5,cycleniter=100,parallel=True)
	
	
# Plotting the amp vs uv dist
plotms(vis='Bullet_CLuster_Hanningsmooth.ms/', selectdata=True, correlation='XX, YY', averagedata=True, \
	avgchannel='3723', xaxis='uvdist', yaxis='amp', field='2')
	
# Plot the Antennas
plotants(vis='Bullet_CLuster_Hanningsmooth.ms/', figfile='Bullet_CLuster_antenna_pos.png') 


# Quack flagging
flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='quack', quackinterval=10.0, quackmode='beg')

# Manual flagging of teh known rfi bands in the MeerKat L band
# GSM Up and Down link
flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:0~119;167~335')

#GPS
flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:167~335; 3230~3325; 1564~1660; 2320~2378; 1320~1416')

#GLONASS
flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:3359~3445; 1684~1718; 1493')

#GALILEO
flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:3279; 1270; 1517; 1444; 1860')

#IRIDIUM and INMARSAT
flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:3474~3521;3043~3177')

# Setjy on the flux calibrator
execfile('setjy_manual.py')

execfile('setjy_polarisation_coefficients.py')
setjy(vis='Bullet_CLuster_Hanningsmooth.ms',field='3C286',standard='manual',\
      fluxdensity=[15.77,0,0,0], spix=[-0.48757786, -0.1840742, -0.03037325, 0], reffreq='1280MHz',\
      polindex=[0.09598155, 0.02465338, -0.08372889, 0.19818536, 0], \
      polangle=[25.48946014, 8.48508625, -11.05647654, 1.3602341, 0],\
      usescratch=False,scalebychan=True,spw='')
      
#setjy(vis='Bullet_CLuster_Hanningsmooth.ms/',field='3C286',model='3C286_L.im',usescratch=False,scalebychan=True,spw='')

# Initial Phase Calibration
gaincal(vis='Bullet_CLuster_Hanningsmooth.ms/', caltable='Bullet_CLuster_Hanningsmooth.G0all/', field='0,1,3,4', \
	refant='m014', spw='0:1200~1203', calmode='p', solint='int', minsnr=5)
	
gaincal(vis='Bullet_CLuster_Hanningsmooth.ms/', caltable='Bullet_CLuster_Hanningsmooth.G0/', field='0408-658', \
	refant='m014', spw='0:1000~1003', calmode='p', solint='int', minsnr=5)
	
plotms(vis='Bullet_CLuster_Hanningsmooth.G0',xaxis='time',yaxis='phase',coloraxis='corr',iteraxis='antenna',plotrange=[-1,-1,-180,180])
	
#Flagging some bad data
flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, scan='8, 15', \
	antenna='m001, m004, m007, m008, m009, m010, m018, m039, m050, m057, m059, m062, m063') 
	
flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, scan='57', \
	antenna='m022, m023,m024, m025,m041, m042,m059') 
	
#Delay Calibration
gaincal(vis='Bullet_CLuster_Hanningsmooth.ms/',caltable='Bullet_CLuster_Hanningsmooth.K0', field='0408-658',\
	refant='m014',spw='0:5~3717',gaintype='K',  solint='inf',combine='scan',minsnr=5, \
   	gaintable=['Bullet_CLuster_Hanningsmooth.G0'])
   	
plotms(vis='Bullet_CLuster_Hanningsmooth.K0/',xaxis='antenna1',yaxis='delay',coloraxis='baseline')

	
# Bandpass calibration only on the bandpass calibrators
bandpass(vis='Bullet_CLuster_Hanningsmooth.ms/',caltable='Bullet_CLuster_Hanningsmooth.B0', field='0408-658',\
	spw='',refant='m014',combine='scan',solint='inf',bandtype='B',\
	gaintable=['Bullet_CLuster_Hanningsmooth.G0','Bullet_CLuster_Hanningsmooth.K0']) 
	
plotms(vis='Bullet_CLuster_Hanningsmooth.B0/',field='0408-658', xaxis='chan',yaxis='amp',coloraxis='corr', \
	plotrange=[-1,-1,0,3],iteraxis='antenna',gridrows=2,gridcols=2)
	
plotms(vis='Bullet_CLuster_Hanningsmooth.B0/',field='0408-658', xaxis='chan',yaxis='phase',coloraxis='corr',\
	plotrange=[-1,-1,-180,180], iteraxis='antenna',gridrows=2,gridcols=2)
	
# Some alternative bad band flagging

flagdata(vis='Bullet_CLuster_Hanningsmooth.ms/', mode='manual', flagbackup=True, spw='0:660~680')   

flagdata(vis='Bullet_CLuster_Hanningsmooth.ms', mode='manual', flagbackup=True, field='0408-658', \
	spw='0:660~680; 925~975; 1400~1945; 2950~3530; 1000~1050')                                         
	
# Gain Calibration (amplitude and phase)
gaincal(vis='Bullet_CLuster_Hanningsmooth.ms/',caltable='Bullet_CLuster_Hanningsmooth.G1', field='0,1,3,4',\
	spw='0:5~3720', solint='inf',refant='m014',gaintype='G',calmode='ap',solnorm=False, \
	gaintable=['Bullet_CLuster_Hanningsmooth.K0/','Bullet_CLuster_Hanningsmooth.B0/'], interp=['linear','nearest'])
	
plotms(vis='Bullet_CLuster_Hanningsmooth.G1',xaxis='time',yaxis='phase'\
       ,iteraxis='corr',coloraxis='baseline',plotrange=[-1,-1,-180,180])
       
plotms(vis='Bullet_CLuster_Hanningsmooth.G1',xaxis='time',yaxis='amp', iteraxis='corr',coloraxis='baseline')

plotms(vis='Bullet_CLuster_Hanningsmooth.G1', xaxis='time', yaxis='phase', correlation='/', coloraxis='baseline', plotrange=[-1,-1,-180,180])


# Polarization calibration (set the polarization with 3C286)
#This does not work
#setjy(vis='Bullet_CLuster_Hanningsmooth.ms',field='3C286',standard='Perley-Butler 2017',
#      model='3C286_L.im', polindex = [9.60, 1.93e-3, -5.11e-6, 9.45e-9, 0], polangle = [25.49, 6.63e-3, 6.75e-6, 6.49e-10, 0], \
#      rotmeas = 1.50, usescratch=False,scalebychan=True,spw='')

#alpha = log(15.7726/18.639)/log(1.27913/0.889997) = -0.46
#i0=15.77 # Stokes I value at 1.28 GHz
#c0= 0.086 ~ 0.099 # Fractional polarization=8.6~9.9%
#d0=33*pi/180 # polarization angle of 33 degrees converted to radians
#setjy(vis='Bullet_CLuster_Hanningsmooth.ms', field='3C286', standard='manual',\
#      spw='0', fluxdensity=[i0,0,0,0], spix=[alpha,0], reffreq='1.28GHz',\
#      polindex=[c0,0], polangle=[d0,0], scalebychan=True, usescratch=False)

## Check if it works      
#i0=15.77 # Stokes I value for spw 0 ch 0
#p0=0.086~0.099*i0 # Fractional polarization=11.2%
#q0=p0*cos(66*pi/180) # Stokes Q for spw 0 for pang = 33 deg (Q+iU phase = 66 deg)
#u0=p0*sin(66*pi/180) # Stokes U for spw 0 for pang = 33 deg (Q+iU phase = 66 deg)


	
# Solving Cross hand delays

gaincal(vis='Bullet_CLuster_Hanningsmooth.ms', caltable='Bullet_CLuster_Hanningsmooth.Kcross',\
        field='3C286', spw='0:5~3720',\
        gaintype='KCROSS', solint='inf', combine='scan', refant='m014',\
        gaintable=['Bullet_CLuster_Hanningsmooth.K0',\
                   'Bullet_CLuster_Hanningsmooth.B0',\
                   'Bullet_CLuster_Hanningsmooth.G1'],
        gainfield=['','',''], parang=True)
        

plotms(vis='Bullet_CLuster_Hanningsmooth.Kcross',xaxis='antenna1',yaxis='delay',coloraxis='corr')

# Solving Leakage terms 

polcal(vis='Bullet_CLuster_Hanningsmooth.ms',caltable='Bullet_CLuster_Hanningsmooth.D1',\
       field='0647-475',spw='0:5~3720', refant='m014',poltype='Df+QU',solint='inf',combine='scan',\
       gaintable=['Bullet_CLuster_Hanningsmooth.K0',\
                  'Bullet_CLuster_Hanningsmooth.B0',\
                  'Bullet_CLuster_Hanningsmooth.G1',\
                  'Bullet_CLuster_Hanningsmooth.Kcross'],\
                  gainfield=['','','0647-475',''])
                  
plotms(vis='Bullet_CLuster_Hanningsmooth.D1',xaxis='chan',yaxis='amp', iteraxis='antenna',coloraxis='corr', field='0647-475')

plotms(vis='Bullet_CLuster_Hanningsmooth.D1',xaxis='chan',yaxis='phase', iteraxis='antenna',coloraxis='corr',plotrange=[-1,-1,-180,180])

plotms(vis='Bullet_CLuster_Hanningsmooth.D1',xaxis='antenna1',yaxis='amp',coloraxis='corr', field='0647-475')

# use this for a source with unkown polarization

polcal(vis='Bullet_CLuster_Hanningsmooth.ms',caltable='Bullet_CLuster_Hanningsmooth.D2',\
       field='polarized source',spw='0:5~3720', refant='m014',poltype='Df+QU',solint='inf',combine='scan',\
       gaintable=['Bullet_CLuster_Hanningsmooth.K0',\
                  'Bullet_CLuster_Hanningsmooth.B0',\
                  'Bullet_CLuster_Hanningsmooth.G1',\
                  'Bullet_CLuster_Hanningsmooth.Kcross'],\
       gainfield=['','','polarized source',''])

# Solving the polarization angle

polcal(vis='Bullet_CLuster_Hanningsmooth.ms',caltable='Bullet_CLuster_Hanningsmooth.X1',\
       field='3C286',combine='scan', refant='m014',poltype='Xf',solint='inf',\
       gaintable=['Bullet_CLuster_Hanningsmooth.K0',\
                  'Bullet_CLuster_Hanningsmooth.B0',\
                  'Bullet_CLuster_Hanningsmooth.G1',\
                  'Bullet_CLuster_Hanningsmooth.Kcross',\
                  'Bullet_CLuster_Hanningsmooth.D1'],\
       gainfield=['','','','',''])
       
plotms(vis='Bullet_CLuster_Hanningsmooth.X1',xaxis='chan',yaxis='phase')

# Scaling the amplitude gains (setting the flux of all calibrators

nog een keer gaincal ie. G2 maken 

myscale = fluxscale(vis='Bullet_CLuster_Hanningsmooth.ms',\
                    caltable='Bullet_CLuster_Hanningsmooth.G1',\ 
                    fluxtable='Bullet_CLuster_Hanningsmooth.fluxscale1',\ 
                    reference=['06...'],\
                    transfer=['gain calibrators'],\
                    incremental=False)

plotms(vis='Bullet_CLuster_Hanningsmooth.fluxscale1',xaxis='time',yaxis='amp', correlation='X',coloraxis='baseline')
plotms(vis='Bullet_CLuster_Hanningsmooth.fluxscale1',xaxis='time',yaxis='amp', correlation='Y',coloraxis='baseline')


# Apply the calibration 

applycal(vis='Bullet_CLuster_Hanningsmooth.ms',\
         field='3C286',\
         gaintable=['Bullet_CLuster_Hanningsmooth.fluxscale1',\
                    'Bullet_CLuster_Hanningsmooth.K0',\
                    'Bullet_CLuster_Hanningsmooth.B0',\
                    'Bullet_CLuster_Hanningsmooth.Kcross', \
                    'Bullet_CLuster_Hanningsmooth.D2',\
                    'Bullet_CLuster_Hanningsmooth.X1'],\
         gainfield=['3C286','','','','',''], \
         interp=['nearest','','','','',''],\
         calwt=[False], parang=True)



   	
   	


