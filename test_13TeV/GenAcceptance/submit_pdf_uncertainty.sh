#!/bin/bash
echo "Submitting pdf_uncertainty.sh"

numevents=0
generator="Powheg"
options="nominal"
cme="13TeV"
pdfset="NNPDF31"
ntupledir="/eos/user/a/aandriat/wz/$cme/ntuples/$generator/$options/$pdfset"
acceptdir=$PWD


if [ ! "$CMSSW_BASE" ]; then
  echo "-------> error: define cms environment."
  exit 1
fi
cmsbase=$CMSSW_BASE

memory=30000
diskspace=30000
masterqueue="2nw"


mkdir -p logs/

chmod 744 pdf_uncertainty.sh #permissions

# echo "Calculate pdf acceptance uncertainty"
# root -l pdf_uncertainty_acceptance.C+\(\"$ntupledir\",$numevents\) -q |& tee -a process_log.txt

bsub -oo logs/"$generator"_"$options"_"$cme"_"$pdfset"_submit.txt -q ${masterqueue} -C 0  -R "rusage[mem=${memory}:pool=${diskspace}]" "export PRODHOME=`pwd`; pdf_uncertainty.sh" $numevents $generator $options $cme $pdfset $ntupledir $acceptdir $cmsbase

echo "Pdf_uncertainty_acceptacne submitted"
echo
echo

rm *.so *.d *.pcm

echo
echo "Job submitted"
echo
