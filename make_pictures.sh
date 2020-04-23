#!/bin/bash

source ./environment

set -u
set -e

set -o pipefail

umask 0002


# environment: version, family, publish


gradsscriptdir=$scriptdir/grads


rsync="/usr/bin/rsync -v"
rotate_dates=$scriptdir/rotate_dates.sh


outputdir=${OUTPUTDIR}

# The new transfers through shared drive:
#fmi_data_path=/fmi/data/silam.fmi.fi/AQ/$suite/$family/

picture_dir=$outputdir/webloads/$fctype/$fcdate${outsuff}

analysis_time=${fcdate}_${fchour}

[ -d $picture_dir ] || mkdir -p $picture_dir

cd $gradsscriptdir
pwd

   mpdset="cosmo_d2"
   PUTLOGOCMD=$PUTLOGOCMDPOLLEN
   #pollen run parallelised by species: reason: long variables, different scripts for taxons
   d=`date -u -d"$maxhours hours $fcdate" +"%Y%m%d%H"`
   ctl=$outputdir/${fcdate}${outsuff}/ALL_SRCS_${fctype}_${d}.nc4.ctl
   ls $ctl || exit 10

   for a in "b"; do #dummy loop
     # cincentrations made with .ctl file  to access proper names
     grads_script="${gradsscriptdir}/draw_pollen.gs"
     for taxon in aphids alder birch grass mugwort olive ragweed; do
        echo "\"run $grads_script $ctl $analysis_time $taxon $picture_dir/${taxon}_srf $mpdset\"" 
        #grads -bpc "run $grads_script $ctl $analysis_time $taxon $picture_dir/${taxon}_srf $mpdset"
     done

#    # Pollen indices are all in hourly files, no need for .ctl
#    ioff=0
#    grads_script=$gradsscriptdir/main_process_pollenindex.gs
#    for h in `seq 1 $maxhours`; do ## Do not plot 0-th hour
#      d=`date -u -d"$h hours $fcdate" +"%Y%m%d%H"`
#      pidx_binary=$outputdir/$fcdate${outsuff}/ALL_SRCS_${fctype}_PI_${d}.nc
#      echo "\"run $grads_script   $analysis_time ${pidx_binary}  $picture_dir $ioff $mpdset\"" 
#      ioff=`expr $ioff + 1`
#    done 

   done | time xargs -t -l -P 6 grads -bpc

exit 0

# put logo if corresponding command is provided
#cd $scriptdir/delme || exit 234
if [  -z ${PUTLOGOCMD+x}  ]; then
   echo PUTLOGOCMD is not set
else
   echo Putting logos
   ls ${picture_dir}/*.png |grep -v AQI_ |grep -v POLI_| time aprun -n1 -d28  xargs -n 1 -I XXX -P 28 ${PUTLOGOCMD} XXX XXX  
   if [  -n ${PUTLOGOCMDINDEX+x}  ]; then
      #separate logo for AQI
      ls ${picture_dir}/*[AO][QL]I_???.png | time aprun -n1 -d28  xargs -n 1 -I XXX -P 28 ${PUTLOGOCMDINDEX} XXX XXX  
   fi
   #ls ${picture_dir}/*.png | time xargs -n 1 -I XXX -P 20 ${PUTLOGOCMD} XXX XXX  
fi
#echo waiting..
#wait



echo pics Done!
[ -n "$outsuff" ] && exit 0 #No publish for modified runs


if $publish; then
    keepdays=7
    pushd $picture_dir/..
    for d in `find . -type d -name '20??????'|sort|head -n -$keepdays`; do
       echo removing $d
       rm -rf $d
    done
    popd

    echo Rotatedirs
    $rotate_dates $outputdir/webloads/$fctype 3

    echo Syncing $outputdir/webloads/$fctype to $fmi_data_path
    mkdir -p $fmi_data_path
    $rsync -a --delete  $outputdir/webloads/$fctype $fmi_data_path
fi
exit 0

