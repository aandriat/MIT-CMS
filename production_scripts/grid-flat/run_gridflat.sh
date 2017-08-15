#!/bin/bash
filename=${1}
process=${2}
cme=${3}
generator=${4}
options=${5}
pdfset=${6}
val=${7}
sampledir=${8}
outputdir=${9}
grid_base=${10}
bacon_base=${11}

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
mkdir -p temp/confs
mkdir -p $grid_base/grid-flat/logs/$cme/$generator/$options/$pdfset/$process
eos mkdir -p $outputdir/$cme/samples/$generator/$options/$pdfset/$process/
eos mkdir -p $outputdir/$cme/ntuples/$generator/$options/$pdfset/$process/
cp $grid_base/grid-flat/config/$cme/$generator/$options/* temp/confs/
cp -r $bacon_base/MIT-CMS/test_13TeV/* temp/
cp $grid_base/grid-flat/flatten_bacon_parallel.C temp/FlattenBacon/

cd temp/confs

echo
echo "Substituting actual sample path in config file"
echo
sed -i 's:TARBALL:'$sampledir/$filename':g' cmsRun_cfg.py

cd $grid_base
eval `scramv1 runtime -sh`
export PYTHONHOME=`scram tool info python | grep PYTHON_BASE | sed 's/PYTHON_BASE=//'`
back

cd temp/confs
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

cd temp/confs
echo
echo "Starting cmsRun"
echo
cmsRun makingGenBacon_MC.py
echo
echo "Finished cmsRun"
echo

# echo
# echo "Saving $filename to $outputdir"
# echo
# xrdcopy -f Output.root $outputdir/$cme/samples/$generator/$options/$process/"$process"_bacon_"$val".root

echo "$PWD"
for f in *; do
  echo "File -> $f"
done
back

mv temp/confs/Output.root temp/FlattenBacon/
cd temp/FlattenBacon/
for f in *; do
  echo "File -> $f"
done
echo "Flatten Bacon"
root -l flatten_bacon_parallel.C+\(\"Output.root\",\"$process\",\".\",0\) -q |& tee $grid_base/grid-flat/logs/$cme/$generator/$options/$pdfset/$process/"$process"_flatten_"$val".txt
echo "Done Flattening"

echo "Moving file"
xrdcopy -f flat_output.root $outputdir/$cme/ntuples/$generator/$options/$pdfset/$process/"$process"_gen_"$val".root
echo "File $process_gen_$val.root moved to $outputdir/$cme/ntuples/$generator/$pdfset/$options"

back
rm -rf temp/

echo
echo "$PWD"
echo
echo "Exit"
echo
exit 0