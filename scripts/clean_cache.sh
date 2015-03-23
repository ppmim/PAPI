#!/bin/bash
# next command is used to drop memory cache and try to recover GEIRS from 
# its 'Buffer allocation error"
# It has to be executed as root user
sync
echo 3 > /proc/sys/vm/drop_caches 
