#!/bin/bash
echo "Submitting flatten.sh"

numevents=${1}
ntupledir=${2}
generator=${3}

flatdir=$PWD

memory=8000
diskspace=8000
masterqueue="1nd"

if [ ! "$CMSSW_BASE" ]; then
  echo "-------> error: define cms environment."
  exit 1
fi
cmsbase=$CMSSW_BASE

mkdir -p logs/

chmod 744 flatten.sh #permissions

if [ -n "$4" ]
  then
    f=${4}
    echo "File -> $f"
  	bsub -oo logs/"$f"_log.txt -q ${masterqueue} -C 0  -R "rusage[mem=${memory}:pool=${diskspace}]" "export PRODHOME=`pwd`; flatten.sh" $f $numevents $ntupledir $generator $flatdir $cmsbase 
  else
    for f in confs/$generator/*; do
  	echo "File -> $f"
  	bsub -oo logs/"$f"_log.txt -q ${masterqueue} -C 0  -R "rusage[mem=${memory}:pool=${diskspace}]" "export PRODHOME=`pwd`; flatten.sh" $f $numevents $ntupledir $generator $flatdir $cmsbase 
done
fi

rm *.so *.d *.pcm

echo
echo "Jobs submitted"
echo
