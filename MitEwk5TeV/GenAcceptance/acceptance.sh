#!/bin/bash

rm acceptances.txt

root -l acceptance.C+\(\"acceptance.conf\",\"/data/t3home000/aandriat/5TeV/ntuples\",0\) -q

rm *.so *.d *.pcm