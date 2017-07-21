#!/bin/bash
echo "Starting flatten.sh"

f=${1}
numevents=${2}
ntupledir=${3}
generator=${4}
flatdir=${5}
cmsbase=${6}

workDir=$PWD
function back {
	cd $workDir
}

cd $cmsbase
eval `scramv1 runtime -sh`
back
if [ ! "$CMSSW_BASE" ]; then
  echo "-------> error: failed to define cms environment."
  exit 1
fi

export EOS_MGM_URL=root://eosuser.cern.ch
eos mkdir -p $ntupledir/$generator
mkdir -p $flatdir/logs/

mkdir temp
cp -r $flatdir/../* temp/
cd temp/FlattenBacon

for a in *; do
  echo "File -> $a"
done

echo "$PWD"
echo "Flatten Bacon"
root -l flatten_gen.C+\(\"confs/$generator/"$f"\",\".\",$numevents\) -q |& tee $flatdir/logs/"$f"_root_log.txt
echo "Done Flattening"

for b in *; do
  echo "File -> $b"
done

echo "Moving file"
xrdcopy -f *.root root://eosuser.cern.ch/$ntupledir/$generator/
echo "File moved to /$ntupledir/$generator/"
back

rm -rf temp/

echo
echo "Done"
echo 