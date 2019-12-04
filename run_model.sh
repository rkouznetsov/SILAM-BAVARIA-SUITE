#!/bin/sh
set -u
set -e
cd $scriptdir
START_TIME=`date -u -d $fcdate +"%Y %m %d 00 00 0."`
export START_TIME

# mpirun -n $ny $silam_binary $scriptdir/MoscowICON.ctrl 2>&1 | tee $outputdir/${fcdate}/outlog_`date -u +"%Y%m%d%H%M%S"`.log
$silam_binary $scriptdir/MoscowICON.ctrl 2>&1 | tee $outputdir/${fcdate}/outlog_`date -u +"%Y%m%d%H%M%S"`.log



