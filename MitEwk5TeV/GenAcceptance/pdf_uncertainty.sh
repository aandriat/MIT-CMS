#!/bin/bash

rm pdf_acceptances.txt

root -l pdf_uncertainty_acceptance.C+\(\"acceptance.conf\",\"/data/t3home000/aandriat/5TeV/ntuples\",0\) -q

rm *.so *.d *.pcm
