#!/bin/bash
memory=${1}
diskspace=${2}
masterqueue=${3}
filename=${4}
process=${5}

if [ ! "$CMSSW_BASE" ]; then
  echo "-------> error: define cms environment."
  exit 1
fi

echo
echo "CHMOD"
echo
chmod 744 run_gridlhe.sh

echo
echo "Submitting run_gridlhe.sh"
echo
val=0
numJobs=1
while((val<numJobs)); do
    echo $val
    bsub -oo "$process"_log_"$val".txt -q ${masterqueue} -C 0  -R "rusage[mem=${memory}:pool=${diskspace}]" "export PRODHOME=`pwd`; run_gridlhe.sh" $filename $process $val
    val=$(($val+1))
done

echo 
echo "Submitted"
echo 
exit 0



