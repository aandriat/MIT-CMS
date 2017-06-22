#!/bin/bash

rm -rf ntuples/
mkdir -p ntuples/

rm acceptances.txt

root -l flatten_gen.C+\(\"acceptance.conf\",\".\",1000\) -q
root -l acceptance.C+\(\"acceptance.conf\",\"ntuples\",0\) -q

rm *.so *.d *.pcm
