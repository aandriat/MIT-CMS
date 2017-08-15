#!/bin/bash
echo "Starting pdf_uncertainty.sh"

numevents=${1}
generator=${2}
options=${3}
cme=${4}
pdfset=${5}
ntupledir=${6}
acceptdir=${7}
cmsbase=${8}


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

mkdir -p $acceptdir/acceptances/

mkdir temp
cp -r $acceptdir/pdf_uncertainty_acceptance.C temp/
cd temp/

echo "Calculate pdf acceptance uncertainty"
root -l pdf_uncertainty_acceptance.C+\(\"$ntupledir\",$numevents\) -q |& tee "$generator"_"$options"_"$cme"_"$pdfset"_log.txt

rm $acceptdir/acceptances/"$generator"_"$options"_"$cme"_"$pdfset"_acceptance.txt
rm $acceptdir/acceptances/"$generator"_"$options"_"$cme"_"$pdfset"_acceptance.root

mv pdf_acceptances.txt $acceptdir/acceptances/"$generator"_"$options"_"$cme"_"$pdfset"_acceptance.txt
mv pdf_acceptances.root $acceptdir/acceptances/"$generator"_"$options"_"$cme"_"$pdfset"_acceptance.root

echo "Gen_Acceptance Complete"
echo
echo

rm *.so *.d *.pcm