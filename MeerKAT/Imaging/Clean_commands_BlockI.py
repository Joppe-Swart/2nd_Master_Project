"""
here we collect the commands we use.


"""
# average the data (100 channels left)
DP3 msin=BC_Block_I_Calibrated.ms msout=BC_Averaged_I.ms steps=[av] av.type=averager av.freqstep=2 av.timestep=2 msin.datacolumn=DATA msout.storagemanager=dysco

# split the data into 4
DP3 msin=BC_Averaged_I.ms filter.nchan=[26,26,26,26] filter.starchan=[0,26,52,78] msout.name=[BC_Averaged_I_1.ms, BC_Averaged_I_2.ms, BC_Averaged_I_3.ms, BC_Averaged_I_4.ms] steps=[filter, split] split.replaceparms=[filter.nchan, filter.startchan, msout.name] msin.datacolumn=DATA msout.storagemanager=dysco

DP3 msin=BC_Averaged_I.ms steps=[split] split.type=split split. split.steps=[out] msin.datacolumn=DATA msout.storagemanager=dysco

# Wsclean command for QU imaging 
wsclean -reorder -name BC_I -pol QU -size 8192 8192 -scale 1.1arcsec -weight briggs 0 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -parallel-deconvolution 1024 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -data-column DATA -multiscale -multiscale-max-scales 4 -join-polarizations -join-channels -squared-channel-joining -channels-out 25 BC_Averaged_I.ms/

#wsclean for I
wsclean -reorder -name BC_I -pol I -size 8192 8192 -scale 1.1arcsec -weight briggs 0 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -parallel-deconvolution 1024 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -data-column DATA -multiscale -multiscale-max-scales 4 -join-channels -squared-channel-joining -channels-out 25 BC_Averaged_I.ms/
