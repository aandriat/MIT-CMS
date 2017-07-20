#!/bin/bash
echo "Starting flatten.sh"

process=${1}
generator="Madgraph5"
topdir=$PWD
ntupledir="/data/t3home000/aandriat/13TeV/ntuples"

if [ ! "$CMSSW_BASE" ]; then
  echo "-------> error: define cms environment."
  exit 1
fi

function back {
	cd $topdir
}

echo "Make ntuple dir"
mkdir -p $ntupledir/$generator

echo "Clear process log"
mkdir -p $topdir/logs/
rm $topdir/logs/"$process"_log.txt

echo "Flatten Bacon"
root -l flatten_gen.C+\(\"flatten_bacon_"$process".conf\",\"$ntupledir/$generator\",0\) -q |& tee -a $topdir/logs/"$process"_log.txt

rm *.so *.d *.pcm
