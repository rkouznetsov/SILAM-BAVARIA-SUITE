#!/bin/bash
set -u -e

set -o pipefail

umask 0002

nproc=`nproc`
#nproc="1 echo"  # this disables all xargs calls

# environment: version, family, publish


gradsscriptdir=$scriptdir/grads


rsync="/usr/bin/rsync -v"

outputdir=${OUTPUTDIR}

#outputdir=${OUTPUT_DIR}

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

   done | time xargs -t -l -P $nproc grads -bpc


# put logo if corresponding command is provided
#cd $scriptdir/delme || exit 234
if [  -z ${PUTLOGOCMD+x}  ]; then
   echo PUTLOGOCMD is not set
else
   echo Putting logos
   ls ${picture_dir}/*.png |grep -v AQI_ |grep -v POLI_| xargs  -I XXX -P $nproc ${PUTLOGOCMD} XXX XXX  
#   if [  -n ${PUTLOGOCMDINDEX+x}  ]; then
#      #separate logo for AQI
#      ls ${picture_dir}/*[AO][QL]I_???.png | xargs  -I XXX -P $nproc ${PUTLOGOCMDINDEX} XXX XXX  
#   fi
   echo compressing pics
   ls ${picture_dir}/*.png  | xargs  -I XXX -P $nproc convert XXX PNG8:XXX
fi
#echo waiting..
#wait



echo Done with logos!
[ -n "$outsuff" ] && exit 0 #No publish for modified runs


if $publish; then
    keepdays=7
    pushd $picture_dir/..
    for d in `find . -type d -name '20??????'|sort|head -n -$keepdays`; do
       echo removing $d
       rm -rf $d
    done
    ii=0
    for d in `find . -type d -name '20??????'|sort -r`; do
       linkname=`printf  %03d $ii`
       rm -f $linkname
       echo ln -s  $d $linkname
       ln -sf  $d $linkname
       ii=`expr $ii + 1`  
    done

    #deploy animation if not yet...
    rsync -av $scriptdir/www/*.html .
    if [ !  -d Napit ]; then
     tar -xvf  $scriptdir/www/Napit.tar
    fi

    popd
    fmi_data_path=eslogin:/fmi/data/silam.fmi.fi/partners/Bavaria
    echo Syncing $outputdir/webloads to $fmi_data_path
#    mkdir -p $fmi_data_path
    $rsync -a --delete  $outputdir/webloads/$fctype/* $fmi_data_path/
fi
exit 0

