#!/bin/bash
echo "Starting flatten.sh"

echo "Make ntuple dir"
mkdir -p /data/t3home000/aandriat/13TeV/ntuples/Madgraph5

echo "Clear process log"
rm process_log.txt

echo "Flatten Bacon"
root -l flatten_gen.C+\(\"flatten_bacon.conf\",\"/data/t3home000/aandriat/13TeV/ntuples/Madgraph5\",0\) -q |& tee -a process_log.txt

rm *.so *.d *.pcm
