#!/bin/sh

set -u
set -e 


infile=`ls  /home/yury/Oper_out/${fcdate}/oper_pm_*.nc4.ctl |tail -n 1`
outDir="/home/yury/Oper_out/pics/${fcdate}/"

for sp in O3 NO2 CO NO PM10 PM2_5 SO2; do
  echo "run plotonevar.gs $infile $outDir $sp"
done | xargs -P8 -IXXX -t sh -c "grads -blc \"XXX\""

