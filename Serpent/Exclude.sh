#!/bin/sh
# This script produces the list of (un)suitable Farux nodes for a specific problem (n CPUs, mem Gb of memory)
# These actions are necessary to overcome the problems of large heterogeneity of this cluster
# Excluded : (Total node memory [Mb]) - truncate((Number of CPUs)/(n CPUs))*(mem Gb)*(1024 Mb/Gb) < 0
# Included : (Total node memory [Mb]) - truncate((Number of CPUs)/(n CPUs))*(mem Gb)*(1024 Mb/Gb) > 0
# Usage : ./Exclude.sh
n=5 # CPUs
mem=30 # Gb of memory required
sinfo -Nel | tail -n +3 | awk '{if ($7-int($5/3)*20*1024<0) print $1}' | sort -n > Exclude.txt
sinfo -Nel | tail -n +3 | awk '{if ($7-int($5/3)*20*1024>0) print $1}' | sort -n > Include.txt
