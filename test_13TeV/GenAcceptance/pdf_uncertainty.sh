#!/bin/bash
echo "Starting pdf_uncertainty.sh"

generator=${1}
numevents=${2}

echo "Make ntuple dir"
ntupledir="/eos/user/a/aandriat/wz/13TeV/ntuples/$generator/"
export EOS_MGM_URL=root://eosuser.cern.ch
eos mkdir -p $ntupledir

echo "Clear process log"
rm process_log.txt

echo "Calculate pdf acceptance uncertainty"
root -l pdf_uncertainty_acceptance.C+\(\"$ntupledir\",$numevents\) -q |& tee -a process_log.txt

echo "Gen_Acceptance Complete, log in process_log.txt"
echo
echo

rm *.so *.d *.pcm