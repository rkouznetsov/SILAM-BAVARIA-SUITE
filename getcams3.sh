#/bin/sh

#. environment
#exec 5> debug_output.txt
#BASH_XTRACEFD="5"
set -u
set -e
#set -x 

fcdate=`date -u -d "4 hours ago" +%Y%m%d%H`
#fcdate=2019120512

#adjust it to  00Z or 12Z fcst
fch=`expr \(  $fcdate % 100 \) / 12 \* 12` || true  # Prevent from failure if 0 returned
fch=`printf %02d $fch`
fcdate=`echo $fcdate |cut -b 1-8`$fch

echo  fcdate=$fcdate


urlpref=https://opendata.dwd.de/weather/nwp/icon-eu/grib/${fch}


outdir=$ICON_PATH/$fcdate
mkdir -p $outdir


##working directory for temporary files
scratchdir=/tmp/getICON${fcdate}
rm -rf $scratchdir
mkdir -p $scratchdir
cd $scratchdir



ncurls=60

var3d="p qv t tke u v w"
#clc onlu at pressure levels!


staticvar="fr_lake fr_land hsurf plcov soiltyp"
staticvarnn="soiltyp"
staticnotuse="elat elon clat clon"

static3d="hhl"



notuse="
fi #Geopotential at pressure level
relhum  # Rh @pressure level
"

#varsfc="alb_rad alhfl_s asob_s asob_t aswdifd_s aswdifu_s aswdir_s cape_con cape_ml clch clcl clcm clct clct_mod cldepth fr_ice h_snow hbas_con htop_con htop_dc hzerocl pmsl ps qv_s rain_con rain_gsp relhum_2m rho_snow snow_con snow_gsp t_2m t_g t_snow td_2m tmax_2m tmin_2m tot_prec u_10m v_10m vmax_10m w_snow ww z0 "
varsfc="alb_rad alhfl_s asob_s asob_t aswdifd_s aswdifu_s aswdir_s cape_con cape_ml clch clcl clcm clct clct_mod cldepth h_snow hbas_con htop_con htop_dc hzerocl pmsl ps qv_s rain_con rain_gsp relhum_2m rho_snow snow_con snow_gsp t_2m t_g t_snow td_2m tmax_2m tmin_2m tot_prec u_10m v_10m vmax_10m w_snow ww z0 "


varsoil=""
for soillev in 0 1 3 9 27 81 729; do ###Layers in soil
  varsoil="$varsoil ${soillev}_w_so"
done
for soillev in 0 2 6 18 54 162 486 1458; do  ##levels in soil
  varsoil="$varsoil ${soillev}_t_so"
done


pref_static=icon-eu_europe_regular-lat-lon_time-invariant_
pref_3d=icon-eu_europe_regular-lat-lon_model-level_
pref_3dh=icon-eu_europe_regular-lat-lon_hybrid-level_
pref_sfc=icon-eu_europe_regular-lat-lon_single-level_
pref_soil=icon-eu_europe_regular-lat-lon_soil-level_
cutgrib="cdo sellonlatbox,32.2,43,50.2,62" 

#export https_proxy=http://wwwproxy.fmi.fi:8080
##
##Get Static
filepref=$pref_static
staticfile=$outdir/${filepref}${fcdate}.grib2.bz2
if [ ! -f $staticfile ]; then
  for itry in `seq 1 10`; do
    for v in $staticvar; do
      V=`echo $v | tr  '[:lower:]' '[:upper:]'`
      destfile=${filepref}${fcdate}_${V}.grib2.bz2
      [ -s $destfile  ] && continue
      echo "curl -s -f --connect-timeout 5 --max-time 10 --retry 5 --retry-delay 2 --retry-max-time 70 $urlpref/$v/$destfile -o ${destfile}.tmp && mv ${destfile}.tmp ${destfile} && echo ${destfile} Done, try $itry || echo ${destfile} Failed, try $itry.."
      #curl --verbose   -f $urlpref/$v/$destfile -o $destfile
    done |xargs -P $ncurls  -I XXX sh -c "XXX"
    base=`basename $staticfile .bz2`
    nmsg=`ls ${filepref}${fcdate}_*.grib2.bz2| wc -w` 
    if [ $nmsg ==  5 ]; then
      break
    else
      echo  "nfiles =  $nmsg Retry $itry"

    fi
    #nmsg=`grib_ls $base.tmp |grep "total messages" |awk '{print $1;}'`
  done
  if [ $itry -ge 9 ]; then
    echo "Failed after 9 tries"
    exit 255
  fi
  lbzip2 -dc ${filepref}${fcdate}_*.grib2.bz2 > $base.tmp
  rm ${filepref}${fcdate}_*.grib2.bz2
  
  $cutgrib $base.tmp $base
  lbzip2 $base
  rm $base.tmp 
  mv $base.bz2 $staticfile
fi


steplist="`seq -f %03.0f 1 78` `seq -f%03.0f 81 3 120`"
#steplist="`seq -f %03.0f 1 12`"
#steplist="001 002 003 "

for step in $steplist; do
  file3dh=$outdir/${pref_3dh}${fcdate}+${step}.grib2.bz2 #final file in hybrid
  if [ ! -f $file3dh ]; then
    for itry in `seq 1 10`; do
      for v in $var3d; do
        V=`echo $v | tr  '[:lower:]' '[:upper:]'`
        for lev in `seq 1 60`; do
          destfile=${pref_3d}${fcdate}_${step}_${lev}_${V}.grib2.bz2
          [ -s $destfile  ] && continue
          echo "curl -s --max-time 40 -f $urlpref/$v/$destfile -o ${destfile}.tmp && mv ${destfile}.tmp ${destfile} && echo ${destfile} Done, try $itry || echo ${destfile} Failed, try $itry.."
        done
      done | xargs -P $ncurls  -I XXX sh -c "XXX" || echo Some curls for $file3dh failed
      
      file3d=${pref_3d}${fcdate}+${step}.grib2
      base=`basename $file3dh .bz2`
      nmsg=`ls ${pref_3d}${fcdate}_${step}_*.grib2.bz2|wc -w`
      if [ $nmsg ==  420 ]; then
        break
      else
        echo  "nmsg =  $nmsg Retry $itry"
      fi
      # nmsg=`grib_ls $file3d.tmp |grep "total messages" |awk '{print $1;}'`
    done
    if [ $itry -ge 9 ]; then
      echo "Failed after 9 tries"
      exit 255
    fi
    cat ${pref_3d}${fcdate}_${step}_*.grib2.bz2 | lbzip2 -dvc > $file3d.tmp 
    rm ${pref_3d}${fcdate}_*.grib2.bz2
    $cutgrib $file3d.tmp  $file3d
    $scriptdir/tools/ICONEU2hybrid $file3d $base
    lbzip2 $base
    rm $file3d.tmp  $file3d
    mv $base.bz2 $file3dh
  fi

  filepref=$pref_sfc
  filesfc=$outdir/${filepref}${fcdate}+${step}.grib2.bz2
  if [ ! -f $filesfc ]; then
    for itry in `seq 1 10`; do
      for v in $varsfc; do
        V=`echo $v | tr  '[:lower:]' '[:upper:]'`
          destfile=${filepref}${fcdate}_${step}_${V}.grib2.bz2
          [ -f $destfile ] || echo "curl -s -f --connect-timeout 5 --max-time 10 --retry 5 --retry-delay 2 --retry-max-time 70 $urlpref/$v/$destfile -o $destfile.tmp && mv ${destfile}.tmp ${destfile}"
      done | xargs -P $ncurls -t -I XXX sh -c "XXX" || echo Some curls for $filesfc failed
      base=`basename $filesfc .bz2`
      nmsg=`ls ${filepref}${fcdate}_${step}_*.grib2.bz2|wc -w`
      if [ $nmsg ==  42 ]; then
        break
      else
        echo  "nmsg_sfc =  $nmsg Retry $itry"
      fi
      # nmsg=`grib_ls $file3d.tmp |grep "total messages" |awk '{print $1;}'`
    done
    if [ $itry -ge 9 ]; then
      echo "Failed after 9 tries"
      exit 255
    fi
    cat  ${filepref}${fcdate}_*.grib2.bz2 | lbzip2 -dc > $base.tmp &&  rm ${filepref}${fcdate}_*.grib2.bz2
    $cutgrib $base.tmp $base
    lbzip2 $base
    rm $base.tmp
    mv $base.bz2 $filesfc

  fi

  filepref=$pref_soil
  filesoil=$outdir/${filepref}${fcdate}+${step}.grib2.bz2
  if [ ! -f $filesoil ]; then
    for itry in `seq 1 10`; do
      for v in $varsoil; do
        V=`echo $v | tr  '[:lower:]' '[:upper:]'`
          destfile=${filepref}${fcdate}_${step}_${V}.grib2.bz2
          vtmp=`echo $v| sed -e 's/.*_w_so/w_so/' -e 's/.*_t_so/t_so/'`  ##Cut soil level from varname
          [ -f $destfile ] || echo "curl -s -f --connect-timeout 5 --max-time 10 --retry 5 --retry-delay 2 --retry-max-time 70 $urlpref/$vtmp/$destfile -o $destfile.tmp && mv ${destfile}.tmp ${destfile}"
      done | xargs -P $ncurls -t -I XXX sh -c "XXX" || echo Some curls for $filesoil failed

      base=`basename $filesoil .bz2`
      nmsg=`ls ${filepref}${fcdate}_${step}_*.grib2.bz2|wc -w`
      if [ $nmsg ==  15 ]; then
        break
      else
        echo  "nmsg_soil =  $nmsg Retry $itry"
      fi
      # nmsg=`grib_ls $file3d.tmp |grep "total messages" |awk '{print $1;}'`
    done
    if [ $itry -ge 9 ]; then
      echo "Failed after 9 tries"
      exit 255
    fi
   
    cat ${filepref}${fcdate}_*.grib2.bz2 | lbzip2 -dc> $base.tmp &&  rm ${filepref}${fcdate}_*.grib2.bz2
    $cutgrib $base.tmp $base
    lbzip2 $base
    rm $base.tmp
    mv $base.bz2 $filesoil
  fi
done


exit 0
#

