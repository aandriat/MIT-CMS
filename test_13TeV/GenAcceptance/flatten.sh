#!/bin/bash

rm -rf /data/t3home000/aandriat/13TeV/ntuples
mkdir -p /data/t3home000/aandriat/13TeV/ntuples

root -l flatten_gen.C+\(\"flatten_bacon.conf\",\"/data/t3home000/aandriat/13TeV/ntuples\",0\) -q

rm *.so *.d *.pcm
