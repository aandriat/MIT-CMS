#!/bin/bash

rm -rf /data/t3home000/aandriat/5TeV/ntuples
mkdir -p /data/t3home000/aandriat/5TeV/ntuples

rm acceptances.txt
rm pdf_acceptances.txt

root -l flatten_gen.C+\(\"acceptance.conf\",\"/data/t3home000/aandriat/5TeV/ntuples\",0\) -q
root -l acceptance.C+\(\"acceptance.conf\",\"/data/t3home000/aandriat/5TeV/ntuples\",0\) -q
root -l pdf_uncertainty_acceptance.C+\(\"acceptance.conf\",\"/data/t3home000/aandriat/5TeV/ntuples\",0\) -q


rm *.so *.d *.pcm
