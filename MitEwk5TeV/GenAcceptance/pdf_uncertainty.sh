#!/bin/bash

rm pdf_acceptances.txt

root -l pdf_uncertainty_acceptance.C+\(\"acceptance.conf\",\"ntuples\",0\) -q

rm *.so *.d *.pcm
