#/bin/sh

#. environment
#exec 5> debug_output.txt
#BASH_XTRACEFD="5"
set -u
set -e
#set -x 

fcdate=`date -u -d "3 hours ago" +%Y%m%d%H`

fch=`echo $fcdate |cut -b 9-10`
if [ $fch -lt 2 ]; then
  msg="$fch is less than 2 hours sonce the forecast lead time."
  msg1="The DWD directory might be updating now. Please try latar."
    echo $msg
    echo $msg1

    echo "$msg " | mailx -s "getCOSMO.sh failed  for Bavaria at Haze"  $MAILTO
fi

# -k to ugnore the failure of ssl sertificte
curl="curl -s -f --connect-timeout 5 --max-time 10 --retry 5 --retry-delay 2 --retry-max-time 70"
#curl="curl -k -s -f --connect-timeout 5 --max-time 10 --retry 5 --retry-delay 2 --retry-max-time 70"

#fcdate=2019120512

#03 H is the longest forecast
fch=03
fcdate=`echo $fcdate |cut -b 1-8`${fch}


echo  fcdate=$fcdate
#https://opendata.dwd.de/weather/nwp/icon-d2/grib/03/alb_rad/icon-d2_germany_regular-lat-lon_single-level_2021021103_000_2d_alb_rad.grib2.bz2

urlpref=https://opendata.dwd.de/weather/nwp/icon-d2/grib/${fch}


outdir=$ICON_PATH/$fcdate
mkdir -p $outdir


##working directory for temporary files
scratchdir=/tmp/getICON${fcdate}
rm -rf $scratchdir
mkdir -p $scratchdir
cd $scratchdir



ncurls=60

var3d="clc p qv t tke u v w"


staticvar="fr_lake fr_land fr_ice hsurf soiltyp plcov lai rootdp depth_lk"
staticvarnn="soiltyp"
staticnotuse="elat elon clat clon"

static3d="hhl"



notuse="
fi #Geopotential at pressure level
relhum  # Rh @pressure level
"

#varsfc="alb_rad alhfl_s asob_s asob_t aswdifd_s aswdifu_s aswdir_s cape_con cape_ml clch clcl clcm clct clct_mod cldepth fr_ice h_snow hbas_con htop_con htop_dc hzerocl pmsl ps qv_s rain_con rain_gsp relhum_2m rho_snow snow_con snow_gsp t_2m t_g t_snow td_2m tmax_2m tmin_2m tot_prec u_10m v_10m vmax_10m w_snow ww z0 "
varsfc="alb_rad alhfl_s apab_s ashfl_s asob_s asob_t aswdifd_s aswdifu_s aswdir_s athb_s athb_t aumfl_s avmfl_s c_t_lk cape_ml ceiling cin_ml  clch clcl clcm clct clct_mod cldepth dbz_850 dbz_cmax echotop freshsnw grau_gsp h_ice h_ml_lk h_snow hbas_sc htop_sc htop_dc hzerocl lpi lpi_max mh pmsl ps qv_s rain_con rain_gsp relhum_2m rho_snow runoff_g runoff_s sdi_2 snow_con snow_gsp snowc snowlmt synmsg_bt_cl_ir10.8 synmsg_bt_cl_wv6.2 t_2m t_bot_lk t_g t_ice t_mnw_lk t_snow t_wml_lk tcond10_mx tcond_max td_2m tmax_2m tmin_2m tot_prec tqc tqc_dia tqg tqi tqi_dia tqr tqs tqv tqv_dia twater u_10m v_10m vmax_10m w_i w_snow ww z0 "
#varsfc="alb_rad alhfl_s apab_s ashfl_s asob_s asob_t aswdifd_s aswdifu_s aswdir_s athb_s athb_t aumfl_s avmfl_s c_t_lk cape_ml ceiling cin_ml clc clch clcl clcm clct clct_mod cldepth dbz_850 dbz_cmax echotop freshsnw grau_gsp h_ice h_ml_lk h_snow hbas_sc htop_sc htop_dc hzerocl lpi lpi_max mh pmsl ps qv_s rain_con rain_gsp relhum_2m rho_snow runoff_g runoff_s sdi_2 snow_con snow_gsp snowc snowlmt synmsg_bt_cl_ir10.8 synmsg_bt_cl_wv6.2 t_2m t_bot_lk t_g t_ice t_mnw_lk t_snow t_wml_lk tch tcm tcond10_mx tcond_max td_2m tmax_2m tmin_2m tot_prec tqc tqc_dia tqg tqi tqi_dia tqr tqs tqv tqv_dia twater u_10m v_10m vmax_10m w_i w_snow ww z0 "

#varsfc="alb_rad alhfl_s ashfl_s asob_s asob_t aswdifd_s aswdifu_s aswdir_s cape_ml cin_ml  clch clcl clcm clct clct_mod cldepth dbz_850 dbz_cmax  grau_gsp h_snow hzerocl  mh  pmsl prs_gsp ps  qv_s relhum_2m rho_snow rootdp runoff_g runoff_s sdi_1 sdi_2 snow_gsp snowlmt t_2m t_g t_s t_snow  tch tcm td_2m tot_prec u_10m v_10m w_snow  ww z0"

varsoil=""
for soillev in 0 1 3 9 27 81 243 729; do ###Layers in soil
  varsoil="$varsoil ${soillev}_w_so ${soillev}_w_so_ice ${soillev}_smi"
done
for soillev in 0 2 6 18 54 162 486 1458; do  ##levels in soil
  varsoil="$varsoil ${soillev}_t_so"
done


pref_static=icon-d2_germany_regular-lat-lon_time-invariant_
pref_3d=icon-d2_germany_regular-lat-lon_model-level_
pref_3dh=icon-d2_germany_regular-lat-lon_hybrid-level_
pref_sfc=icon-d2_germany_regular-lat-lon_single-level_
pref_soil=icon-d2_germany_regular-lat-lon_soil-level_
#cutgrib="cdo sellonlatbox,32.2,43,50.2,62" 

#export https_proxy=http://wwwproxy.fmi.fi:8080
##
##Get Static
filepref=$pref_static
staticfile=$outdir/${filepref}${fcdate}.grib2
if [ ! -f $staticfile ]; then
  for itry in `seq 1 10`; do
    for v in $staticvar; do
      V=$v #`echo $v | tr  '[:lower:]' '[:upper:]'`
      destfile=${filepref}${fcdate}_000_0_${V}.grib2.bz2
      [ -s $destfile  ] && continue

      echo "$curl $urlpref/$v/$destfile -o ${destfile}.tmp && mv ${destfile}.tmp ${destfile} && echo ${destfile} Done, try $itry || echo ${destfile} Failed, try $itry.."
      >&2 echo "curl --verbose   -f $urlpref/$v/$destfile -o $destfile"
    done |xargs -P $ncurls  -I XXX sh -c "XXX"
    base=`basename $staticfile`
    nmsg=`ls ${filepref}${fcdate}_*.grib2.bz2| wc -w` 
#https://opendata.dwd.de/weather/nwp/icon-d2/grib/03/plcov/icon-d2_germany_regular-lat-lon_time-invariant_2021021103_000_0_plcov.grib2.bz2
#https://opendata.dwd.de/weather/nwp/icon-d2/grib/03/plcov/icon-d2_germany_regular-lat-lon_time-invariant_2021021103_plcov.grib2.bz2
    if [ $nmsg ==  9 ]; then
      break
    else
      echo  "nfiles =  $nmsg Retry $itry"

    fi
    #nmsg=`grib_ls $base.tmp |grep "total messages" |awk '{print $1;}'`
  done
  if [ $itry -ge 9 ]; then
    msg="Failed static after 9 tries"
    echo $msg
    echo "$msg " | mailx -s "getCOSMO.sh failed  for Bavaria at Haze"  $MAILTO
    exit 255
  fi
  lbzip2 -dc ${filepref}${fcdate}_*.grib2.bz2 > $base.tmp
  rm ${filepref}${fcdate}_*.grib2.bz2
  
  grib_set -s packingType=grid_ccsds $base.tmp $base
  rm $base.tmp 
  mv $base $staticfile
fi


steplist="`seq -f %03.0f 1 45`"
#steplist="`seq -f %03.0f 1 12`"
#steplist="001 002 003 "

for step in $steplist; do
  file3dh=$outdir/${pref_3dh}${fcdate}+${step}.grib2 #final file in hybrid
  file3d=$outdir/${pref_3d}${fcdate}+${step}.grib2
  if [ ! -f $file3dh ]; then 
    echo Getting stuff to `pwd`
    for itry in `seq 1 10`; do
      for v in $var3d; do
        V=$v #`echo $v | tr  '[:lower:]' '[:upper:]'`
#        https://opendata.dwd.de/weather/nwp/icon-d2/grib/03/u/icon-d2_germany_regular-lat-lon_model-level_2021021103_001_44_U.grib2.bz2

#        https://opendata.dwd.de/weather/nwp/icon-d2/grib/03/u/icon-d2_germany_regular-lat-lon_model-level_2021021103_010_3_u.grib2.bz2
        maxlev=65
        [ "$V" == "w" ] && maxlev=66
        for lev in `seq 1 $maxlev`; do
          destfile=${pref_3d}${fcdate}_${step}_${lev}_${V}.grib2.bz2
          [ -s $destfile  ] && continue
          
          patchfile=$scriptdir/`basename ${destfile}`
          [ -f $patchfile ] && cp $patchfile ${destfile} && continue ### FIXME A hook to pick up a patch from $scriptdir
# patch for missing file           
#          if [ $itry -gt 3 ]; then 
#            if [ $destfile == "icon-d2_germany_regular-lat-lon_model-level_2021102003_015_63_tke.grib2.bz2" ]; then
#              moldfile=icon-d2_germany_regular-lat-lon_model-level_2021102003_015_62_tke.grib2.bz2
#              destbase=`basename  $destfile .bz2` 
#              echo "lbzip2 -dc $moldfile > delme && grib_set -s scaledValueOfFirstFixedSurface=63  delme $destbase && lbzip2 $destbase  && rm delme && echo Recovered $destfile"
#              continue
#            fi
#          fi 
          echo "$curl $urlpref/$v/$destfile -o ${destfile}.tmp && mv ${destfile}.tmp ${destfile} && echo ${destfile} Done, try $itry || echo ${destfile} Failed, try $itry.."
        done
      done | xargs -P $ncurls -I XXX  sh -c "XXX" || echo Some curls for $file3dh failed
      
      base=`basename $file3dh`
      nmsg=`ls ${pref_3d}${fcdate}_${step}_*.grib2.bz2|wc -w`
      if [ $nmsg ==  521 ]; then
        break
      else
        echo  "nmsg =  $nmsg Retry $itry"
      fi
      # nmsg=`grib_ls $file3d.tmp |grep "total messages" |awk '{print $1;}'`
    done
    if [ $itry -ge 9 ]; then
      msg="Failed 3D after 9 tries"
      echo $msg
      echo "$msg " | mailx -s "getCOSMO.sh failed  for Bavaria at Haze"  $MAILTO
      exit 255
    fi
    cat ${pref_3d}${fcdate}_${step}_*.grib2.bz2 | lbzip2 -dvc > $base
    grib_set -w edition=2 -s packingType=grid_ccsds  $base $file3d.tmp
    mv $file3d.tmp $file3d
    rm ${pref_3d}${fcdate}_*.grib2.bz2
    $scriptdir/tools/ICON-d2-2hybrid $file3d $base
    mv $base $file3dh
  fi

  filepref=$pref_sfc
  filesfc=$outdir/${filepref}${fcdate}+${step}.grib2
  if [ ! -f $filesfc ]; then
    for itry in `seq 1 10`; do
      for v in $varsfc; do
        V=$v #`echo $v | tr  '[:lower:]' '[:upper:]'`
          destfile=${filepref}${fcdate}_${step}_2d_${V}.grib2.bz2
          [ -f $destfile ] || echo "$curl $urlpref/$v/$destfile -o $destfile.tmp && mv ${destfile}.tmp ${destfile}"
      done | xargs -P $ncurls -t -I XXX sh -c "XXX" || echo Some curls for $filesfc failed
      base=`basename $filesfc`
      nmsg=`ls ${filepref}${fcdate}_${step}_*.grib2.bz2|wc -w`
      if [ $nmsg ==  84 ]; then
        break
      else
        echo  "nmsg_sfc =  $nmsg Retry $itry"
      fi
      # nmsg=`grib_ls $file3d.tmp |grep "total messages" |awk '{print $1;}'`
    done
    if [ $itry -ge 9 ]; then
      msg="Failed sfc after 9 tries"
      echo $msg
      echo "$msg " | mailx -s "getCOSMO.sh failed  for Bavaria at Haze"  $MAILTO
      exit 255
    fi
    cat  ${filepref}${fcdate}_*.grib2.bz2 | lbzip2 -dvc > $base.tmp 
    rm ${filepref}${fcdate}_*.grib2.bz2
    # get rid of 15-min steps
    stepmin=`expr ${step} \* 60`

    grib_set -w stepRange=${stepmin} -S -s packingType=grid_ccsds $base.tmp $base
    grib_set -w stepRange=0-${stepmin} -S -s packingType=grid_ccsds $base.tmp $base.cum
    cat $base.cum >> $base
    rm $base.tmp $base.cum
    mv $base $filesfc

  fi

  if [ $step -gt 27 ]; then
    echo Soil-level disabled for step $step
    continue
  fi

  filepref=$pref_soil
  filesoil=$outdir/${filepref}${fcdate}+${step}.grib2
  if [ ! -f $filesoil ]; then
    for itry in `seq 1 10`; do
      for v in $varsoil; do
        V=$v #`echo $v | tr  '[:lower:]' '[:upper:]'`
          destfile=${filepref}${fcdate}_${step}_${V}.grib2.bz2
          vtmp=`echo $v| sed -e 's/.*_w_so/w_so/' -e 's/.*_t_so/t_so/'`  ##Cut soil level from varname
          [ -f $destfile ] || echo "$curl $urlpref/$vtmp/$destfile -o $destfile.tmp && mv ${destfile}.tmp ${destfile}"
      done | xargs -P $ncurls -t -I XXX sh -c "XXX" || echo Some curls for $filesoil failed

      base=`basename $filesoil .bz2`
      nmsg=`ls ${filepref}${fcdate}_${step}_*.grib2.bz2|wc -w`
      if [ $nmsg ==  24 ]; then
        break
      else
        echo  "nmsg_soil =  $nmsg Retry $itry"
      fi
      # nmsg=`grib_ls $file3d.tmp |grep "total messages" |awk '{print $1;}'`
    done
    if [ $itry -ge 9 ]; then
      msg="Failed soil after 9 tries"
      echo $msg
      echo "$msg " | mailx -s "getCOSMO.sh failed  for Bavaria at Haze"  $MAILTO
      exit 255
    fi
   
    cat ${filepref}${fcdate}_*.grib2.bz2 | lbzip2 -dc> $base.tmp &&  rm ${filepref}${fcdate}_*.grib2.bz2
    grib_set -s packingType=grid_ccsds $base.tmp $base
    rm $base.tmp
    mv $base $filesoil
  fi
done
echo Cleanup..
cd $scriptdir
rm -rvf $scratchdir
echo Done! 
#https://opendata.dwd.de/weather/nwp/icon-d2/grib/03/z0/
#icon-d2_germany_regular-lat-lon_single-level_2021021103_001_2d_z0.grib2.bz2
#icon-d2_germany_regular-lat-lon_single-level_2021021103_045_2d_z0.grib2.bz2
exit 0
#
icon-d2_germany_regular-lat-lon_soil-level_2021021103_001_0_t_so.grib2.bz2
icon-d2_germany_regular-lat-lon_soil-level_2021021103_001_2d_0_t_so.grib2.bz2
