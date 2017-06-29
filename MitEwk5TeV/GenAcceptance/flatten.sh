#!/bin/bash

rm -rf /data/t3home000/aandriat/5TeV/ntuples
mkdir -p /data/t3home000/aandriat/5TeV/ntuples

root -l flatten_gen.C+\(\"acceptance.conf\",\"/data/t3home000/aandriat/5TeV/ntuples\",0\) -q

rm *.so *.d *.pcm
