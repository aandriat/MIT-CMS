Alexander Andriatis 06-23-17

GenAcceptance calculates the generator-level event acceptance from MC Bacon samples.

flatten_gen.C takes Bacon MC samples and stores the 4-vectors of gen-level particles from select final states as individual flat-ntuple files. Important to note it is dependent on helper functions toolbox::flavor() and toolbox::fillGen() located in ../Utils/MyTools.hh

acceptance.C reads flat-ntuple files and applies kinematic and geometric cuts to caclulate gen-level detector acceptance. Stores fiducial events in flat-ntuples

pdf_uncertainty.C uses fiducial flat n-tuples to calculate acceptance error from PDF uncertainty through replica pdf weights stored in lheweights

acceptance.conf defines the samples and final states being considered

ConfParse.hh parses the conf file

gen_acceptance.sh runs the full script


