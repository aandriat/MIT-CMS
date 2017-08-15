#!/bin/bash

filename="Z_slc6_amd64_gcc481_CMSSW_7_1_28_my_Zmm.tgz"
process="zmm"

memory=1000
diskspace=1000
masterqueue="1nh"

workdir=$PWD

chmod 744 submit_gridlhe.sh

echo
echo "Submitting submit_gridlhe.sh"
echo

bsub -oo logs/metasubmit_"$filename".txt -q ${masterqueue} -C 0  -R "rusage[mem=${memory}:pool=${diskspace}]" "export PRODHOME=`pwd`; submit_gridlhe.sh" $filename $process $workdir

echo 
echo "Submitted"
echo 