#!/bin/bash
filename=${1}
gensampledir="/afs/cern.ch/work/a/aandriat/public/powheg_testing/CMSSW_7_1_26/src/gensamples"
outputdir="/afs/cern.ch/work/a/aandriat/public/TestTeV"
workdir=`pwd`

rm -rf temp_process
mkdir -p temp_process

cp makingGenBacon_MC.py temp_process/
cd temp_process/

echo
echo "Changing reference file in makingGenBacon_MC.py"
echo
sed -i 's:GENFILE:'$gensampledir/$filename.root':g' makingGenBacon_MC.py

echo
echo "Starting cmsRun"
echo
cmsRun makingGenBacon_MC.py


echo
echo "Saving $filename to $outputdir"
echo
mv Output.root $outputdir/"$filename"_bacon.root

echo
echo "Done"
echo