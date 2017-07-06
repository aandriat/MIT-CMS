#!/bin/bash
echo "Starting ratios.sh"

echo "Clear acceptances"
rm pdf_acceptance_ratios.root
rm pdf_acceptance_ratios.txt

echo "pdf_acceptance_ratios.C"
root -l pdf_acceptance_ratios.C+ -q

rm *.so *.d *.pcm
