#!/bin/bash
filename=${1}
process=${2}
val=${3}
sampledir="/afs/cern.ch/work/a/aandriat/public/powheg_testing/CMSSW_7_1_26/src/genproductions/bin/Powheg"
gridgendir="/afs/cern.ch/work/a/aandriat/public/powheg_testing/CMSSW_7_1_26/src/grid-gen/"
gensampledir="/afs/cern.ch/work/a/aandriat/public/powheg_testing/CMSSW_7_1_26/src/gensamples/"

echo 
echo $filename
echo

cd $gridgendir

export PYTHONHOME=`scram tool info python | grep PYTHON_BASE | sed 's/PYTHON_BASE=//'`

rm -rf $process
mkdir -p $process
mkdir -p $gensampledir
cp cmsRun_cfg.py $process/
cd $process

echo
echo "Substituting actual sample path in config file"
echo
sed -i 's:/cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc481/13TeV/madgraph/V5_2.4.2/AnomalousCouplingsSMP/QCD/v1/WMhadZlepJJ_QCD_LO_SM_mjj100_pTj10_tarball.tar.xz:'$sampledir/$filename':g' cmsRun_cfg.py

echo
echo "Starting cmsRun"
echo
cmsRun cmsRun_cfg.py
echo
echo "Finished cmsRun"
echo

echo
echo "Copying output LHE to gensamples"
echo
mv SMP-RunIISummer15wmLHEGS-00128.root $gensampledir/"$process"_genLHE_"$val".root
echo "Done"

cd $gridgendir
echo
echo "$PWD"
echo
echo "Exit"
echo
