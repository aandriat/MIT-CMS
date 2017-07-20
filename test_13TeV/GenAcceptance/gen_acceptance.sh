#!/bin/bash
echo "Starting gen_acceptance.sh"

echo "Make ntuple dir"
mkdir -p /data/t3home000/aandriat/13TeV/ntuples/Madgraph5

echo "Clear process log"
rm process_log.txt

echo "Flatten Bacon"
root -l flatten_gen.C+\(\"flatten_bacon.conf\",\"/data/t3home000/aandriat/13TeV/ntuples/Madgraph5\",0\) -q |& tee -a process_log.txt

echo "Calculate pdf acceptance uncertainty"
root -l pdf_uncertainty_acceptance.C+\(\"/data/t3home000/aandriat/13TeV/ntuples/Madgraph5\",0\) -q |& tee -a process_log.txt

echo "Gen_Acceptance Complete, log in process_log.txt"
echo
echo

rm *.so *.d *.pcm
