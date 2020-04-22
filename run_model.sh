#!/bin/bash
set -u
set -e
set -o pipefail
cd $scriptdir
START_TIME=`date -u -d $fcdate +"%Y %m %d 00 00 0."`
export START_TIME

# mpirun -n $ny $silam_binary $scriptdir/MoscowICON.ctrl 2>&1 | tee $outputdir/${fcdate}/outlog_`date -u +"%Y%m%d%H%M%S"`.log
mkdir -p $OUTPUTDIR/${fcdate}
/usr/bin/time -v $silam_binary $scriptdir/BavarianPollen.ctrl 2>&1 | tee $OUTPUTDIR/${fcdate}/outlog_`date -u +"%Y%m%d%H%M%S"`.log
OMP_NUM_THREADS=1
#gdb -ex 'break globals::set_error'  -ex 'break globals::ooops' -ex r --args ${silam_binary}_debug $scriptdir/BavarianPollenDebug.ctrl



