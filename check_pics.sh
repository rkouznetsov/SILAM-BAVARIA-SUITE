#!/bin/bash
set -u
set -e
set -o pipefail
cd $scriptdir


testfile=output/webloads/pollen/$fcdate/ragweed_srf_048.png

for i in `seq 720`; do
  sleep 20
  echo `date`, sleeping 
  [ -f $testfile ] && exit 0
done  

if [ ! -f $testfile ] ; then
   msg="Moded  run Bavaria at Haze has not complete within 2 hours, missing file after 2 hours $testfile"
   echo ====================================================================
   echo
   echo 
   echo $msg
   echo sent to $MAILTO
   echo
   echo
   echo ====================================================================

   echo "$msg " | mailx -s "Model run incomplete for Bavaria at Haze" $MAILTO
   exit 0 # Do not break the suite

fi




