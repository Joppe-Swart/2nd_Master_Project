"""
This script performs some more manual flagging on the
solutions from teh facet selfcal
"""

from facetselfcal import flaglowamps, flaghighamps, applycal
import losoto
import losoto.lib_operations
from losoto import h5parm

#flag the lower solutions
#flaglowamps('merged_selfcalcyle003_BC_Test_Block_I.ms.avg.h5', lowampval=0.60, flagging=True, setweightsphases=True)
#flaghighamps('merged_selfcalcyle003_BC_Test_Block_I.ms.avg.h5', highampval=1.50, flagging=True, setweightsphases=True)

#flaglowamps('merged_selfcalcyle003_BC_Test_Block_II.ms.avg.h5', lowampval=0.60, flagging=True, setweightsphases=True)
#flaghighamps('merged_selfcalcyle003_BC_Test_Block_II.ms.avg.h5', highampval=1.50, flagging=True, setweightsphases=True)

#flaglowamps('merged_selfcalcyle003_BC_Test_Block_III.ms.avg.h5', lowampval=0.60, flagging=True, setweightsphases=True)
#flaghighamps('merged_selfcalcyle003_BC_Test_Block_III.ms.avg.h5', highampval=1.50, flagging=True, setweightsphases=True)

#losoto merged_selfcalcyle003_BC_Test_Block_I.ms.avg.h5 losoto_flag_apgrid.parset

applycal('BC_Test_Block_I.ms.avg', 'merged_selfcalcyle003_BC_Test_Block_I.ms.avg.h5', \
         msincol='DATA', msoutcol='CORRECTED_DATA', msout='.', dysco=True)

applycal('BC_Test_Block_II.ms.avg', 'merged_selfcalcyle003_BC_Test_Block_II.ms.avg.h5', \
         msincol='DATA', msoutcol='CORRECTED_DATA', msout='.', dysco=True)

applycal('BC_Test_Block_III.ms.avg', 'merged_selfcalcyle003_BC_Test_Block_III.ms.avg.h5', \
         msincol='DATA', msoutcol='CORRECTED_DATA', msout='.', dysco=True)

#wsclean -reorder -name BC_Selfcal -size 8192 8192 -scale 1.1arcsec -data-column CORRECTED_DATA -weight briggs -0.5 -weighting-rank-filter 3 -auto-mask 3.0 -auto-threshold 1 -niter 50000 -clean-border 2 -mgain 0.8 -fit-beam -join-channels -gap-channel-division -channels-out 12 -pol I -multiscale -parallel-deconvolution 1024 -use-wgridder -mem 50 BC_Test_Block_I.ms.avg BC_Test_Block_II.ms.avg BC_Test_Block_III.ms.avg


