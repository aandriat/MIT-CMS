#!/bin/bash
echo "Submitting flatten.sh"

numevents=0
generator="Madgraph5"
options="Photos"
cme="13TeV"

ntupledir="/eos/user/a/aandriat/wz/$cme/ntuples/$generator/$options"

flatdir=$PWD

confsdir="$flatdir/confs/$cme/$generator/$options"



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

if [ -n "$1" ]
  then
    f=${1}
    echo "File -> $f"
    bsub -oo logs/"$f"_log.txt -q ${masterqueue} -C 0  -R "rusage[mem=${memory}:pool=${diskspace}]" "export PRODHOME=`pwd`; flatten.sh" $f $numevents $ntupledir $confsdir $flatdir $cmsbase 
  else
    cd $confsdir
    for f in *; do
     cd $flatdir
  	 echo "File -> $f"
  	 bsub -oo logs/"$f"_log.txt -q ${masterqueue} -C 0  -R "rusage[mem=${memory}:pool=${diskspace}]" "export PRODHOME=`pwd`; flatten.sh" $f $numevents $ntupledir $confsdir $flatdir $cmsbase 
    done
fi

cd $flatdir
rm *.so *.d *.pcm

echo
echo "Jobs submitted"
echo
