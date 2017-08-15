#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc481
if [ -r CMSSW_7_1_28/src ] ; then 
 echo release CMSSW_7_1_28 already exists
else
scram p CMSSW CMSSW_7_1_28
fi
cd CMSSW_7_1_28/src
eval `scram runtime -sh`

export X509_USER_PROXY=$HOME/private/personal/voms_proxy.cert
curl -s --insecure https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/SMP-RunIISummer15wmLHEGS-00128 --retry 2 --create-dirs -o Configuration/GenProduction/python/SMP-RunIISummer15wmLHEGS-00128-fragment.py 
[ -s Configuration/GenProduction/python/SMP-RunIISummer15wmLHEGS-00128-fragment.py ] || exit $?;

scram b
cd ../../
cmsDriver.py Configuration/GenProduction/python/SMP-RunIISummer15wmLHEGS-00128-fragment.py --fileout file:lhegensim.root --mc --eventcontent RAWSIM,LHE --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring,Configuration/GenProduction/randomizeSeeds.randomizeSeeds --datatier GEN-SIM,LHE --conditions MCRUN2_71_V1::All --beamspot Realistic50ns13TeVCollision --step LHE,GEN --magField 38T_PostLS1 --python_filename cmsRun_cfg.py --no_exec -n 1000 || exit $? ; 

echo "Made SMP Configuration"

# echo "nothing" ;cmsDriver.py Configuration/GenProduction/python/SMP-RunIISummer15wmLHEGS-00128-fragment.py --fileout file:SMP-RunIISummer15wmLHEGS-00128.root --mc --eventcontent DQM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --datatier DQM --conditions MCRUN2_71_V1::All --beamspot Realistic50ns13TeVCollision --step LHE,GEN,VALIDATION:genvalid_all --magField 38T_PostLS1  --fileout file:SMP-RunIISummer15wmLHEGS-00128_genvalid.root --mc -n 1000 --python_filename SMP-RunIISummer15wmLHEGS-00128_genvalid.py --dump_python --no_exec || exit $? ;

# cmsDriver.py step2 --filein file:SMP-RunIISummer15wmLHEGS-00128_genvalid.root --conditions MCRUN2_71_V1::All --mc -s HARVESTING:genHarvesting --harvesting AtJobEnd --python_filename SMP-RunIISummer15wmLHEGS-00128_genvalid_harvesting.py --no_exec || exit $? ; 

# mv SMP-RunIISummer15wmLHEGS-00128_1_cfg.py ../CMSSW_7_1_26/src/grid-gen/cmsRun_cfg.py
