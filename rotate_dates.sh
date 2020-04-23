#!/bin/bash
path_local=$1
days_back=$2

# This scripts creates symlinks in a given directory as follows:
# - assume today is 15.6.2001, days_back = 2
# 000 -> 20010615
# 001 -> 20010614
# 002 -> 20010613
# Existing links are overwritten.

if [[ -z $2 ]]; then
    echo 'Usage: rotate_dates.sh directory'
    exit 1
fi


#if [ $fcdate == `date -u -d "3 hours" +"%Y%m%d"` ]; then 
   pushd $path_local
   for i in `seq 0 $days_back`; do
       now=`date -u -d" -$i days $fcdate" +%Y%m%d`
       i3=`printf '%03i' $i`
       rm -f ${i3}
       echo ln -s $now ${i3}
       ln -s $now ${i3}
   done
   popd
#else
#  echo The fcdate $fcdate seems not to be the latest. Not touching latest symlinks.
#fi

