#!/bin/bash
echo "Clear acceptances"
rm pdf_acceptances.root
rm pdf_acceptances.txt

echo "Clear process log"
rm process_log.txt

echo "Calculate pdf acceptance uncertainty"
root -l pdf_uncertainty_acceptance.C+\(\"/data/t3home000/aandriat/13TeV/ntuples/powheg/\",0\) -q |& tee -a process_log.txt


rm *.so *.d *.pcm
