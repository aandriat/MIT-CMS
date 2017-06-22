#!/bin/bash

rm -rf ntuples/
mkdir -p ntuples/

root -l flatten_gen.C+\(\"acceptance.conf\",\".\",1000\) -q

rm *.so *.d *.pcm
