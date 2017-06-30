#!/bin/bash

rm pdf_acceptances.root
rm pdf_acceptances.txt

root -l pdf_uncertainty_acceptance.C+\(\"/data/t3home000/aandriat/5TeV/ntuples\"\) -q


rm *.so *.d *.pcm
