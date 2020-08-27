#!/bin/bash
set -u
set -e
set -o pipefail
cd $scriptdir


testfile=output/webloads/pollen/$fcdate/ragweed_srf_048.png

sleep 7200

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




