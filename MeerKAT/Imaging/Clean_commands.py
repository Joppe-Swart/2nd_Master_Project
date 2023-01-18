"""
here we collect the commands we use.


"""

# Wsclean command for QU imaging 
wsclean -reorder -name BC_III -pol QU -size 8192 8192 -scale 1.1arcsec -weight briggs 0 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -parallel-deconvolution 1024 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -data-column CORRECTED_DATA -join-polarizations -join-channels -squared-channel-joining -channels-out 25 BC_Block_III_Calibrated.ms.avg
