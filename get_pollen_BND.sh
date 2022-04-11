#!/bin/bash

# Script to download silam european pollen runs from FMI SILAM thredds setup
# (c)opyleft Roux aka Rostislav DOT Kouznetsov AT fmi.fi 
#Updated  for v5_5 Jun 2020

# Enjoy!

email="test-download-bavaria"
#email=test

set -e

if [ -z "$1" ]; then # get the argument (does not make much sense since only yesterdays data are available)
     basedate=yesterday
else 
     basedate=$1
fi
basedate=`date -u -d "$basedate" +%Y%m%d`

set -u 
ncks=ncks

BND_PATH=SILAM-bnd

targetdir=$BND_PATH/`date -u -d "$basedate" +%Y%m%d00`
mkdir -p $targetdir
cd $targetdir

echo `date` Getting  boundaries to  $targetdir



#echo $SHELL
#exit 1

#thredds/ncss/silam_europe_pollen_v5_7_1-TopSecret/runs/silam_europe_pollen_v5_7_1-TopSecret_RUN_2021-06-07T00:00:00Z
runpref=silam_europe_pollen_v5_8_RUN_
urlbase="http://silam.fmi.fi/thredds/ncss/silam_europe_pollen_v5_8/runs/$runpref"

pollens="POLLEN_ALDER_m22 POLLEN_BIRCH_m22 POLLEN_GRASS_m32 POLLEN_MUGWORT_m18 POLLEN_MUGW1_m18 POLLEN_MUGW2_m18 POLLEN_MUGW3_m18 POLLEN_MUGW4_m18 POLLEN_MUGW5_m18 POLLEN_OLIVE_m28 POLLEN_RAGWEED_m18"

pollenvars="Poll_Rdy2fly heatsum poll_left poll_tot_m2 pollen_corr cnc"
varlist="daymean_temp2m,temp_2m_acc"

for pol in $pollens; do
  for var in $pollenvars; do
    newvar=${var}_${pol}
    echo $newvar | grep 'heatsum_POLLEN_GRASS_m32\|heatsum_POLLEN_MUGW\|poll_tot_m2_POLLEN_MUGW\|pollen_corr_POLLEN_MUGW' >/dev/null && continue

    varlist="$varlist,${newvar}"

  done
done

echo $varlist
#exit 0





# German domain
bbox="spatial=bb&north=56&west=4&east=16&south=46"



maxjobs=4

# make dates
run=`date -u -d $basedate +"%FT00:00:00Z"`



for try  in `seq 0 10`; do
   missfiles=""
   for hr in `seq 24 96` ; do
        time=`date -u -d"$basedate + $hr hours" +"%FT%H:00:00Z"`
        outf=`date -u -d"$basedate + $hr hours" +"SILAM4DE${run}_%Y%m%d%H.nc"`
        [ -f $outf ] && continue
         missfiles="$missfiles $outf"


        URL="$urlbase$run?var=$varlist&$bbox&temporal=range&time_start=$time&time_end=$time&accept=netcdf&email=$email"
#        echo wget \"$URL\"
#        exit
        #
        # For some reason thredds does not supply _CoordinateModelRunDate anymore
        attcmd="-a _CoordinateModelRunDate,global,c,c,$run ${outf}.tmp -a history,global,d,,, -a history_of_appended_files,global,d,,,"
        (wget -q $URL -O ${outf}.tmp && ncatted -h $attcmd && $ncks -h --mk_rec_dmn time ${outf}.tmp $outf && rm ${outf}.tmp && echo  $outf done!) &
        while [ `jobs | wc -l` -ge $maxjobs ]; do sleep 1; done

   done
   wait
   [ -z "$missfiles" ] && break
done


if [ -z "$missfiles" ]; then
  echo "`date` Finished okay after $try attempts"
  exit 0
else
  echo `date` Failed!
  exit 255
fi


