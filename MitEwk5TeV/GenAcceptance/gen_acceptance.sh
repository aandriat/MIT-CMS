#!/bin/bash
echo "Starting gen_acceptance.sh"
echo "Clear ntuples"
rm -rf /data/t3home000/aandriat/5TeV/ntuples
mkdir -p /data/t3home000/aandriat/5TeV/ntuples

echo "Clear acceptances"
rm pdf_acceptances.root
rm pdf_acceptances.txt

echo "Clear process log"
rm process_log.txt

echo "Flatten Bacon"
root -l flatten_gen.C+\(\"flatten_bacon.conf\",\"/data/t3home000/aandriat/5TeV/ntuples\",0\) -q |& tee -a process_log.txt

echo "Calculate pdf acceptance uncertainty"
root -l pdf_uncertainty_acceptance.C+\(\"/data/t3home000/aandriat/5TeV/ntuples\",0\) -q |& tee -a process_log.txt

echo "Gen_Acceptance Complete, log in process_log.txt"
echo
echo

rm *.so *.d *.pcm
