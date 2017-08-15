#!/bin/bash
filename=${1}
process=${2}
val=${3}
sampledir=${4}
gridgendir=${5}
outputdir=${6}
grid_base=${7}
bacon_base=${8}
confdir=${9}

echo
workDir=$PWD
echo $workDir
echo

function back {
	cd $workDir
}

echo 
echo $filename
echo

export EOS_MGM_URL=root://eosuser.cern.ch

mkdir -p temp
eos mkdir -p $outputdir/$process/
cp $confdir/* temp
cd temp

echo
echo "Substituting actual sample path in config file"
echo
sed -i 's:TARBALL:'$sampledir/$filename':g' cmsRun_cfg.py

cd $grid_base
eval `scramv1 runtime -sh`
export PYTHONHOME=`scram tool info python | grep PYTHON_BASE | sed 's/PYTHON_BASE=//'`
back

cd temp
echo
echo "Starting cmsRun"
echo
cmsRun cmsRun_cfg.py
echo
echo "Finished cmsRun"
echo

for f in *; do
  echo "File -> $f"
done

echo
echo "Changing reference file in makingGenBacon_MC.py"
echo
sed -i 's:GENFILE:'lhegensim.root':g' makingGenBacon_MC.py

cd $bacon_base
eval `scramv1 runtime -sh`
export PYTHONHOME=`scram tool info python | grep PYTHON_BASE | sed 's/PYTHON_BASE=//'`
back

cd temp
echo
echo "Starting cmsRun"
echo
cmsRun makingGenBacon_MC.py
echo
echo "Finished cmsRun"
echo

for f in *; do
  echo "File -> $f"
done

echo
echo "Saving $filename to $outputdir"
echo


xrdcopy -f Output.root root://eosuser.cern.ch/$outputdir/$process/"$process"_bacon_"$val".root

echo
echo "Done"
echo

back
rm -rf temp/

cd $gridgendir
echo
echo "$PWD"
echo
echo "Exit"
echo
