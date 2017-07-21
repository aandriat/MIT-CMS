#!/bin/bash

# if [ ! "$CMSSW_BASE" ]; then
#   echo "-------> error: define cms environment."
#   exit 1
# fi

grid_base="/afs/cern.ch/work/a/aandriat/public/wz_analysis/CMSSW_7_1_28/src/"
bacon_base="/afs/cern.ch/work/a/aandriat/public/wz_analysis/CMSSW_7_6_3_patch2/src"

sampledir="/afs/cern.ch/work/a/aandriat/public/wz_analysis/CMSSW_7_1_28/src/tarballs"
gridgendir="/afs/cern.ch/work/a/aandriat/public/wz_analysis/CMSSW_7_1_28/src/grid-gen"
outputdir="/eos/user/a/aandriat/wz/13TeV/samples/Powheg/"

filename=${1}
process=${2}
numJobs=${3}

memory=4000
diskspace=4000
masterqueue="1nh"

# for f in *; do
#   echo "File -> $f"
# done

echo
echo "CHMOD"
echo
chmod 744 run_gridlhe.sh

mkdir -p logs/$process

echo
echo "Submitting run_gridlhe.sh"
echo
val=0
while((val<numJobs)); do
    echo $val
    bsub -oo logs/$process/"$process"_log_"$val".txt -q ${masterqueue} -C 0  -R "rusage[mem=${memory}:pool=${diskspace}]" "export PRODHOME=`pwd`; run_gridlhe.sh" $filename $process $val $sampledir $gridgendir $outputdir $grid_base $bacon_base
    val=$(($val+1))
done

echo 
echo "Submitted"
echo 
exit 0



