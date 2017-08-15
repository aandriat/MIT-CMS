#!/bin/bash

filename=${1}
process=${2}

if [ -n "$3" ]
  then
    workdir=${3}
  else
    workdir=$PWD    
fi

cd $workdir

numJobs=1000
generator="Powheg"
options="nominal"
cme="13TeV"
bacon_base="/afs/cern.ch/work/a/aandriat/public/wz_analysis/CMSSW_7_6_3_patch2/src"
grid_base="/afs/cern.ch/work/a/aandriat/public/wz_analysis/CMSSW_7_1_28/src"

sampledir="$grid_base/tarballs/$cme/$generator"
gridgendir="$grid_base/grid-gen"
outputdir="/eos/user/a/aandriat/wz/$cme/samples/$generator/$options"

memory=12000
diskspace=12000
masterqueue="2nd"

confdir="$gridgendir/config/$cme/$generator/$options"

# for f in *; do
#   echo "File -> $f"
# done

echo
echo "CHMOD"
echo
chmod 744 run_gridlhe.sh

mkdir -p logs/$cme/$generator/$options/$process

echo
echo "Submitting run_gridlhe.sh"
echo
val=0
while((val<numJobs)); do
    echo $val
    bsub -oo logs/$cme/$generator/$options/$process/"$process"_log_"$val".txt -q ${masterqueue} -C 0  -R "rusage[mem=${memory}:pool=${diskspace}]" -J $process$val "export PRODHOME=`pwd`; run_gridlhe.sh" $filename $process $val $sampledir $gridgendir $outputdir $grid_base $bacon_base $confdir
    val=$(($val+1))
done

echo 
echo "Submitted"
echo 



