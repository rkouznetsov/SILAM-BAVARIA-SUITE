#!/bin/bash

# Script to download silam european pollen runs from FMI SILAM thredds setup
# (c)opyleft Roux aka Rostislav DOT Kouznetsov AT fmi.fi 
#Updated  for v5_5 Jun 2020

# Enjoy!

email="ROBERT.GEBAUER@gmx.net"
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
#cd `/bin/pwd` ## Containers have hard time with home directories

echo `date` Getting  boundaries to  $targetdir



#echo $SHELL
#exit 1

runpref=silam_europe_pollen_v5_8_1_RUN_
urlbase="http://thredds.silam.fmi.fi/thredds/ncss/grid/silam_europe_pollen_v5_8_1/runs/$runpref"; ncver="3"
#urlbase="http://thredds.silam.fmi.fi/thredds/ncss/grid/silam_europe_pollen_v5_8/runs/$runpref"; ncver="3"
#urlbase="https://silam.fmi.fi/thredds/ncss/silam_europe_pollen_v5_8_1/runs/$runpref"; ncver=""
#urlbase="https://silam.fmi.fi/thredds/ncss/silam_europe_pollen_v5_8/runs/$runpref"; ncver=""

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



maxjobs=16

## Fri 23 Sep 2022 02:15:55 PM EEST, thredds 5.5 after fix 
# https://github.com/Unidata/tds/issues/286
# maxjobs, Faied files  at first round, total time
#  1           0   3:37 
#  2           1   2:28
#  4           4   1:18
#  8           10  0:49
#
#

## After fix https://github.com/Unidata/tds/issues/293
# maxjobs, Faied files  at first round, total time
#  4           0   1:01
#  8           0  0:39
# 16           0  0:34

ncks --version
ncatted --version


# make dates
run=`date -u -d $basedate +"%FT00:00:00Z"`



for try  in `seq 0 10`; do
   missfiles=""
   for hr in `seq 24 96` ; do
        time=`date -u -d"$basedate + $hr hours" +"%FT%H:00:00Z"`
        outf=`date -u -d"$basedate + $hr hours" +"SILAM4DE${run}_%Y%m%d%H.nc"`
        [ -f $outf ] && continue
        [ -f ${outf}.tmp ] && rm ${outf}.tmp 
         missfiles="$missfiles $outf"


        URL="$urlbase$run?var=$varlist&$bbox&temporal=range&time_start=$time&time_end=$time&accept=netcdf${ncver}&email=$email"
#        echo wget \"$URL\"
#        exit
        #
        # For some reason thredds does not supply _CoordinateModelRunDate anymore
        echo Launching ${outf} 
        attcmd="-a _CoordinateModelRunDate,global,c,c,$run -a history,global,d,,, -a history_of_appended_files,global,d,,,"
        (wget  $URL -O ${outf}.tmp && ncatted -h $attcmd ${outf}.tmp && $ncks -4 -L5 --cnk_plc g2d -h --mk_rec_dmn time ${outf}.tmp $outf && rm ${outf}.tmp && echo  $outf done!) &
        if [ $maxjobs -gt 1 ]; then 
          while [ `jobs | wc -l` -ge $maxjobs ]; do sleep 1; done
        else
          wait
        fi

   done
   wait
   [ -z "$missfiles" ] && break
   echo "Missing files after try $try:"
   echo $missfiles
done


if [ -z "$missfiles" ]; then
  echo "`date` Finished okay after $try attempts"
  exit 0
else
  echo `date` Failed!
  exit 255
fi


