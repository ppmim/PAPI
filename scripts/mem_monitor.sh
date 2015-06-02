#!/bin/bash

# Memory monitoring

p1=`ps -e | fgrep geirs_disp | awk '{print $1}'`
if [[  "$p1" ]]; then
   geirs=`pmap -d $p1 |grep mapped`
   echo `date --iso-8601=seconds` "GEIRS : " $geirs >> /tmp/mem_monitor.txt
fi

p2=`ps -e | fgrep java | awk '{print $1}'`
if [ "$p2" ]; then
   java=`pmap -d $p2 |grep mapped`
   echo `date --iso-8601=seconds` "JAVA  : " $java >> /tmp/mem_monitor.txt
fi

# More memory monitoring

echo " " >> /tmp/mem_ps.txt
echo `date --iso-8601=seconds` "GEIRS :"   >> /tmp/mem_ps.txt
ps uax|grep geirs | grep -v grep | awk '{print $5, $6, $11}' >> /tmp/mem_ps.txt
echo "---" >> /tmp/mem_ps.txt
echo `date --iso-8601=seconds` "JAVA_OT :"  >> /tmp/mem_ps.txt
ps uax|grep java | grep -v grep | awk '{print $5, $6, $11}' >> /tmp/mem_ps.txt

memory=`free |grep Mem`
cache=`free |grep cache:`
echo `date --iso-8601=seconds` "MEMORY  : " $memory >> /tmp/mem_monitor.txt
echo `date --iso-8601=seconds` "CACHE  : " $cache >> /tmp/mem_monitor.txt
