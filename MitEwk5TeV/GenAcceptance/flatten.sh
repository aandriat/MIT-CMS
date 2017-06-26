#!/bin/bash

rm -rf ntuples/
mkdir -p ntuples/

root -l flatten_gen.C+\(\"acceptance.conf\",\".\",0\) -q

rm *.so *.d *.pcm
