#!/bin/bash
set -u
set -e
set -o pipefail
cd $scriptdir
START_TIME=`date -u -d $fcdate +"%Y %m %d 00 00 0."`
export START_TIME

# mpirun -n $ny $silam_binary $scriptdir/MoscowICON.ctrl 2>&1 | tee $outputdir/${fcdate}/outlog_`date -u +"%Y%m%d%H%M%S"`.log

runFailed=false

mkdir -p $OUTPUTDIR/${fcdate}
outdump=$OUTPUTDIR/${fcdate}/outlog_`date -u +"%Y%m%d%H%M%S"`.log

/usr/bin/time -v $silam_binary $scriptdir/BavarianPollen.ctrl 2>&1 | tee $outdump || runFailed=true
if  $runFailed; then
   msg="Failure on exit of model run Bavaria at Haze"
   echo ====================================================================
   echo
   echo 
   echo $msg
   echo sent to $MAILTO
   echo
   echo
   echo ====================================================================

   echo "$msg " | mailx -s "Model run failed  for Bavaria at Haze" -a $outdump $MAILTO
   exit 254

fi

###echo "Run for Bavaria at Haze OK" | mailx -s "Model run ok, please reset initialisation" $MAILTO




#export OMP_NUM_THREADS=1
#gdb -ex 'break globals::set_error'  -ex 'break globals::ooops' -ex r --args ${silam_binary}_debug $scriptdir/BavarianPollenDebug.ctrl



