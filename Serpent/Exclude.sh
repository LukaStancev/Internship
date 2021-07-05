#!/bin/sh
# This script produces the list of (un)suitable Farux nodes for a specific problem that requires "n" CPUs and "mem" Gb of memory
# These actions are necessary to overcome the problems of large heterogeneity of this cluster
# Excluded : (Total node memory [Mb]) - truncate((Number of CPUs)/(n CPUs))*(mem Gb)*(1024 Mb/Gb) < 0
# Included : (Total node memory [Mb]) - truncate((Number of CPUs)/(n CPUs))*(mem Gb)*(1024 Mb/Gb) > 0
# Usage : ./Exclude.sh
n=3 # CPUs
mem=20 # Gb of memory required
sinfo -Nel | tail -n +3 | awk -v n="$n" -v mem="$mem" '{if ($7-int($5/n)*mem*1024<0) print $1}' | sort -n > Exclude.txt
sinfo -Nel | tail -n +3 | awk -v n="$n" -v mem="$mem" '{if ($7-int($5/n)*mem*1024>0) print $1}' | sort -n > Include.txt
