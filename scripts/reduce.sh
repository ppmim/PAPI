#!/bin/bash
# Matilde's project : O2k data reduction 


PAPI=/home/panic/DEVELOP/PIPELINE/PANIC/trunk/papi.py
DIRS_JAN="120103 120104 120105 120106 120107 120108 120109 120110 120111 120112 120113 120114 120115 120118 120119 120127 120128 120129 120130 120131"
DIRS_FEB="120218 120215 120213 120212 120211 120210 120209 120208 120201"

for dir in $DIRS_JAN
do
    if [ ! -d /data2/out/${dir} ]
    then
        mkdir -p /data2/out/${dir}
    fi
    ${PAPI} -c /home/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/papi_suse11.cfg -s /data/O2K/Jan.2012/${dir} -g filter -d /data2/out/${dir} -R lemon
done


for dir in $DIRS_FEB
do
    if [ ! -d /data2/out/${dir} ]
    then
        mkdir -p /data2/out/${dir}
    fi
    ${PAPI} -c /home/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/papi_suse11.cfg -s /data/O2K/Feb.2012/${dir} -g filter -d /data2/out/${dir} -R lemon
done

exit
