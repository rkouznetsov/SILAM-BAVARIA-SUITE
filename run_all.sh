#!/bin/bash




cd /home/yury/SILAM/Apps
echo Starting run_all `date`

. ./environment

set -e
set -u 

d=`date -u +%H`
if [ $d -lt 12 ]; then
  bash ./getcams3.sh 
  echo Only meteo complete at `date`
  
  exit 0
else
  bash ./getcams3.sh
  echo Meteo complete at `date`
fi


./run_model.sh

echo run complete at `date`
./make_total_pm_VBS.sh

echo PM done at `date`

./make_pictures.sh
echo Pictures done at `date`

./movetohtml.sh

./extract_timeseries.sh

