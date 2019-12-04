#!/bin/sh


set -e
set -u
export OMP_NUM_THREADS=20


#ncap2="aprun -d20 -n1 -cc none ncap2"
#ncap2=ncap2
#
# Makes total PM and other measured things
#
datadir=$outputdir/$fcdate/
cd $datadir


inoutfiles=""
for f in $datadir/oper_????????.nc4; do
	outfile=`echo $f |sed -e s/oper_/oper_pm_/`
	inoutfiles="$inoutfiles $f $outfile"
done


scriptfile=ncalscript-$$ 


  dustFine="dust_m_30*1e9f dust_m1_5*1e9f"
  dustCoarse="dust_m6_0*1e9f"
  #pmGFAS="PM_GFAS_m_50*1e9f"
  pmGFAS="PM_FRP_m_17*1e9f"
 
  
  # Add stuff from CB4 
  script=""
  script="$script cnc_NO2_gas=cnc_NO2_gas*46e6f; cnc_NO2_gas@units=\"ug/m3\";"
  script="$script cnc_NO_gas=cnc_NO_gas*30e6f;   cnc_NO_gas@units=\"ug/m3\";"
  script="$script cnc_CO_gas=cnc_CO_gas*28e6f;   cnc_CO_gas@units=\"ug/m3\";"
  script="$script cnc_O3_gas=cnc_O3_gas*48e6f;   cnc_O3_gas@units=\"ug/m3\";"
  script="$script cnc_SO2_gas=cnc_SO2_gas*64e6f; cnc_SO2_gas@units=\"ug/m3\";"
  #
  #script="$script air_dens=air_dens;"
  script="$script BLH=BLH;  Kz_1m = Kz_1m; MO_len_inv = MO_len_inv; air_dens=air_dens;"

#  script="$script cnc_NO2_gas=cnc_NO2_gas*46e6f; cnc_NO2_gas@units=\"ug/m3\";"

  bvbspecies="BVB0_m_50*1e9f BVB1e0_m_50*1e9f BVB1e1_m_50*1e9f BVB1e2_m_50*1e9f BVB1e3_m_50*1e9f"

  pm25_species="SO4_m_20*96e6f SO4_m_70*96e6f NH415SO4_m_20*117e6f NH415SO4_m_70*117e6f NH4NO3_m_70*80e6f sslt_m_05*1e9f sslt_m_50*1e9f $dustFine $pmGFAS AVB0_m_50*1e9f AVB1e0_m_50*1e9f AVB1e1_m_50*1e9f AVB1e2_m_50*1e9f AVB1e3_m_50*1e9f ${bvbspecies} EC_m_50*1e9f mineral_m_50*1e9f"
  pm10_species="NO3_c_m3_0*62e6f sslt_m3_0*1e9f $dustCoarse PM_m6_0*1e9f"
 
  #
  #Make PM
  var="cnc_"
  units=ug/m3
  long_name="Concentration in air"

   #Accumulator and rp vars if needed
   script="$script if (exists(rp)) {rp=rp;  *tmp[\$time, \$height, \$rlat, \$rlon] = 0.f;} else {*tmp[\$time, \$height, \$lat, \$lon] = 0.f;}"
   script="$script tmp@units=\"$units\";"
   
   components=""
   for spec in $pm25_species; do 
        script="$script tmp += ${var}${spec};"
        components="$components ${var}${spec}"
   done
   if [ ! -z "$pm25_species" ]; then
     script="$script ${var}PM2_5 = tmp; ${var}PM2_5@long_name=\"${long_name} ${var}PM2_5\"; ${var}PM2_5@components=\"$components\"; where(${var}PM2_5 < 0.5) ${var}PM2_5=0.5f; ${var}PM2_5@substance_name=\"PM2_5\"; ${var}PM2_5@silam_amount_unit=\"kg\";  ${var}PM2_5@mode_distribution_type=\"NO_MODE\";"
   fi

   for spec in $pm10_species; do 
      script="$script tmp += ${var}${spec};"
      components="$components ${var}${spec}"
   done

   if [ ! -z "$pm10_species" ]; then
     script="$script ${var}PM10 = tmp; ${var}PM10@long_name=\"${long_name} ${var}PM10\"; ${var}PM10@components=\"$components\"; where(${var}PM10 < 0.5) ${var}PM10=0.5f;  ${var}PM10@substance_name=\"PM10\"; ${var}PM10@silam_amount_unit=\"kg\";  ${var}PM10@mode_distribution_type=\"NO_MODE\";"
   fi


  echo $script |sed -e 's/;/;\n/g'> $scriptfile
  echo -n $inoutfiles  |  xargs -n 2 -P 8 -t  ncap2 -4 -L5 -O -v -S $scriptfile
  echo $scriptfile
  echo Done!


#  rm  $scriptfile

exit 0
  #
  echo Compressing logs
  #
  for logfile in q*log; do 
    echo "lbzip2 -n 14 -c $logfile > $outdir/`basename $logfile`.bz2"
  done  |  xargs -I XXX -P 4 -t sh -c "XXX"

  last_dumpfile=`ls $indir/*log |grep -v "run_ALL_SRCS" |tail -n 1`
  head -n 100000 $last_dumpfile > $outdir/`basename $last_dumpfile .log`.head.log
  tail -n 100000 $last_dumpfile > $outdir/`basename $last_dumpfile .log`.tail.log





