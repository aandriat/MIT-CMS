#!/bin/bash

rm acceptances.txt

root -l acceptance.C+\(\"/data/t3home000/aandriat/13TeV/ntuples\",0\) -q

rm *.so *.d *.pcm
