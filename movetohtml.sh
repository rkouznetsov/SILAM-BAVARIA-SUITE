#!/bin/bash

set -u
set -e 

rm -f /var/www/html/silam/00?

for d in 2 1 0; do
 today=`date -u -d "-$d day $fcdate"  '+%Y%m%d'`;
 rsync -av /home/yury/Oper_out/pics/${today} /var/www/html/silam/
 ln -s -f ${today} /var/www/html/silam/00${d}
done
echo clean up 
find  /var/www/html/silam/ -type d -name '????????'| sort | head -n -3| xargs -rt rm -r  # -n 1 rm -R
exit 
