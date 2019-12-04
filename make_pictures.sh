#!/bin/sh

set -u
set -e 

infile='/home/yury/Oper_out/20191203/oper_pm_20191204.nc4.ctl'
outDir='/home/yury/Oper_out/pics/20191203/'

for sp in O3 NO2 CO NO PM10 PM2_5; do
  echo "run plotonevar.gs $infile $outDir $sp"
done | xargs -P8 -IXXX -t sh -c "grads -bpc \"XXX\""

