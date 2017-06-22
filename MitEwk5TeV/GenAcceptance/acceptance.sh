#!/bin/bash

rm acceptances.txt

root -l acceptance.C+\(\"acceptance.conf\",\"ntuples\",0\) -q

rm *.so *.d *.pcm