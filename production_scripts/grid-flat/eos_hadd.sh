#!/bin/bash

export EOS_MGM_URL=root://eosuser.cern.ch
cd /eos/user/a/aandriat/wz/13TeV/ntuples/Powheg/nominal/NNPDF31

kill -9 $(jobs -p)
rm *.root
echo "Deleting failed files"
find . -type f -name "*.root" -size -400000c
find . -type f -name "*.root" -size -400000c -delete
echo "Done cleaning"

jobs
echo "Adding paralell samples"
hadd zmm_gen.root zmm/*.root &> zmm_compression.txt &
hadd zee_gen.root zee/*.root &> zee_compression.txt &
hadd wpe_gen.root wpe/*.root &> wpe_compression.txt &
hadd wpm_gen.root wpm/*.root &> wpm_compression.txt &
hadd wme_gen.root wme/*.root &> wme_compression.txt &
hadd wmm_gen.root wmm/*.root &> wmm_compression.txt &
echo "Waiting to finish"
jobs
wait
echo "Finished adding"

echo "Adding wm and wme"
hadd wm_gen.root wmm_gen.root wpm_gen.root &> wm_compression.txt &
hadd we_gen.root wme_gen.root wpe_gen.root &> we_compression.txt &
echo "Waiting to finish"
jobs
wait
echo "Finished adding"

echo 
echo "Done compressing"
echo 
exit 0
