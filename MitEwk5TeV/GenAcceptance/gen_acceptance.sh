#!/bin/bash

rm -rf ntuples/
mkdir -p ntuples/

rm acceptances.txt
rm pdf_acceptances.txt

root -l flatten_gen.C+\(\"acceptance.conf\",\".\",0\) -q
root -l acceptance.C+\(\"acceptance.conf\",\"ntuples\",0\) -q
root -l pdf_uncertainty_acceptance.C+\(\"acceptance.conf\",\"ntuples\",0\) -q


rm *.so *.d *.pcm
