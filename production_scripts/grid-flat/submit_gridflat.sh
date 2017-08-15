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

cme="13TeV"
generator="Powheg"
options="nominal"
pdfset="NNPDF31"
grid_base="/afs/cern.ch/work/a/aandriat/public/wz_analysis/CMSSW_7_1_28/src"
bacon_base="/afs/cern.ch/work/a/aandriat/public/wz_analysis/CMSSW_7_6_3_patch2/src"


sampledir="$grid_base/tarballs/$cme/$generator/$pdfset"
outputdir="/eos/user/a/aandriat/wz"

memory=4000
diskspace=4000
masterqueue="8nh"


# for f in *; do
#   echo "File -> $f"
# done

echo
echo "CHMOD"
echo
chmod 744 run_gridflat.sh

mkdir -p logs/$cme/$generator/$options/$pdfset/$process

echo
echo "Submitting run_gridflat.sh"
echo
val=0
while((val<numJobs)); do
    echo $val
    bsub -oo logs/$cme/$generator/$options/$pdfset/$process/"$process"_log_"$val".txt -q ${masterqueue} -C 0  -R "rusage[mem=${memory}:pool=${diskspace}]" -J $process$val "export PRODHOME=`pwd`; run_gridflat.sh" $filename $process $cme $generator $options $pdfset $val $sampledir $outputdir $grid_base $bacon_base
    val=$(($val+1))
done

echo 
echo "Submitted"
echo
exit 0



