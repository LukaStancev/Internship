#!/bin/bash
#
# author : A. Hebert
# use    : rdonjon <file.x2m> [<compiler>] [-w]
# note   : <file.x2m> must be located on directory ./data/
#          If <file.access> exists, it is executed.
#          <compiler> is the compiler used to generate the executable
#          -w     to execute in console
#          -quiet quiet execution for regression testing
#
if [ $# = 0 ]
   then
   echo "usage: rdonjon file [absoft]" 1>&2
   exit 1
fi

System=`uname -s`
Sysx="`echo $System | cut -b -6`"
if [ $Sysx = "CYGWIN" ]; then
   MACH=`uname -o`
elif [ $Sysx = "AIX" ]; then
   MACH=`uname -s`
else
   MACH=`uname -sm | sed 's/[ ]/_/'`
fi

xxx=0
term=0
quiet=0
log=0

for param in $*
do
   case $param in
      -w) echo Execute in terminal
          term=1
          ;;
      -quiet)
          quiet=1
          ;;
      *) if [ $log = 0 ]; then
            xxx=`basename $param .x2m`
            log=1
            data=$param
            typ='custom'
         else
            re='^[0-9]+$'
            if ! [[ $param =~ $re ]] ; then
               typ=$param
            else
               intg=$param
            fi
         fi ;;
   esac
done
Code=Donjon
if [ $quiet = 0 ]; then
    echo 'execute' $xxx 'with' $Code 'on' $(hostname) 'system (' $MACH ') with' $typ 'compiler'
fi

if [ -d ../"$MACH" ]; then
  if [ $quiet = 0 ]; then
    echo 'use the existing directory' $MACH
  fi
else
  echo 'creation of directory' $MACH
  mkdir ../"$MACH"
fi
CodeDir=$(dirname "$(readlink -f "$0")")/Version5/Donjon
DataDir=$PWD

if [ $Sysx = "AIX" ]; then
  Tmpdir=/usr/tmp
elif [ $Sysx = "SunOS" ]; then
  Tmpdir=/var/tmp
else
  Tmpdir=/tmp
fi
inum=1
while [ -d $Tmpdir/rundir$inum ]
  do
  inum=`expr $inum + 1 `
done
Rundir=$Tmpdir/rundir$inum$intg
mkdir $Rundir
if [ $quiet = 0 ]; then
  echo "RunDirectory:" $Rundir
fi
cd $Rundir

if [ $typ = 'custom' ]; then
  cp "$CodeDir"/bin/"$MACH"/$Code ./code
else
  cp "$CodeDir"/bin/"$MACH"'_'$typ/$Code ./code
fi
cp "$DataDir"/$data ./mydata

if [ -d "$DataDir"/`echo $xxx`_proc ]; then
  cp "$DataDir"/`echo $xxx`_proc/*.c2m .
fi
if [ -f "$DataDir"/$xxx.access ]; then
  if [ $quiet = 0 ]; then
    "$DataDir"/$xxx.access "$DataDir"/..
  else
    "$DataDir"/$xxx.access "$DataDir"/.. > /dev/null
  fi
fi
if [ -f "$DataDir"/assertS.c2m ]; then
  cp "$DataDir"/assertS.c2m .
fi
if [ -f "$DataDir"/assertV.c2m ]; then
  cp "$DataDir"/assertV.c2m .
fi
before=$(date +%s)
if [ $term = 0 ]; then
  ./code <mydata >$xxx.result
else
  ./code <mydata
fi
if [ $quiet = 0 ]; then
  time=$(( $(date +%s) - before))
  printf 'End of execution. Total execution time: %dh %dmin %ds\n' \
    $(($time/3600)) $(($time%3600/60)) $(($time%60))
fi
if [ -f "$DataDir"/$xxx.save ]; then
  if [ $quiet = 0 ]; then
    "$DataDir"/$xxx.save "$DataDir"/..
  else
    "$DataDir"/$xxx.save "$DataDir"/.. > /dev/null
  fi
fi
mv $xxx.result "$DataDir"/../"$MACH"
if [ $quiet = 0 ]; then
  echo 'the listing is located on ../'$MACH
fi

cd "$DataDir"/../"$MACH"
if [ $quiet = 0 ]; then
   tail -45 $xxx.result
else
  RED='\033[0;31m'
  GREEN='\033[0;32m'
  NC='\033[0m' # No Color
  if tail $xxx.result | grep -q "normal end" ; then
    printf "${GREEN}[OK]${NC}\n"
  else
    printf "${RED}[FAILED]${NC}\n"
  fi
fi
chmod -R 777 $Rundir
/bin/rm -r -f $Rundir
cd ..
