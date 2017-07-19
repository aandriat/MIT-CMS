#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include "TLorentzVector.h"         // 4-vector class
#include "TH1D.h"
#include "TRandom.h"

#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
#include "../Utils/LeptonCorr.hh"   // muon scale and resolution corrections

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions
#endif

void pdf_uncertainty_acceptance(  const TString outputDir=".", const Int_t n_events=0   // ntuple directory
            ){
  
  gBenchmark->Start("pdf_uncertainty_acceptance");
    // Create output directory
  gSystem->mkdir(outputDir,kTRUE);

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
    TString string_err = "_err";

    std::map<TString, Int_t> sample_list = { // List of signals to calculate acceptances for: 0 to skip, 1 if signal process with process_gen.root file, 2 if dervied quantity
      { "wpe" , 1},
      { "wme" , 1},
      { "we" , 1},
      { "wpewme", 2},
      { "zee", 1},
      { "wpezee", 2},
      { "wmezee", 2},
      { "wezee", 2},
      { "wpm", 1},
      { "wmm", 1},
      { "wm", 1},
      { "wpmwmm", 2},
      { "zmm", 1},
      { "wpmzmm", 2},
      { "wmmzmm", 2},
      { "wmzmm", 2},
    };
    std::map<TString, Int_t>::iterator setparameters; // Iterator for looping over sample_list

    std::map<TString, std::pair<TString, TString>> ratio_list = { // Map of derived quantities (ratios) and the signal processes they come from
      { "wpewme", std::pair<TString, TString> ("wpe","wme")},
      { "wpezee", std::pair<TString, TString> ("wpe","zee")},
      { "wmezee", std::pair<TString, TString> ("wme","zee")},
      { "wezee", std::pair<TString, TString> ("we","zee")},
      { "wpmwmm", std::pair<TString, TString> ("wpm","wmm")},
      { "wpmzmm", std::pair<TString, TString> ("wpm","zmm")},
      { "wmmzmm", std::pair<TString, TString> ("wmm","zmm")},
      { "wmzmm", std::pair<TString, TString> ("wm","zmm")},
    };
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // Declare input ntuple variables
  Double_t weightGen;  
  TLorentzVector *glep1=0, *glep2=0, *gvec=0;
  Int_t glepq1, glepq2, fiducial;
  std::vector<float> *lheweight = new std::vector<float>();

  // Declare input file
  TString infilename;
  TFile *infile=0;
  TTree *intree=0;

  // Declare output file
  TString outfilename;
  TFile *outFile=0;
  TTree *outTree=0;

  // Useful strings
  TString parname;
  TString parname_err;
  TString type, vartype;

  // Map of parameters which stores the process name and acceptance as a double, also makes entry for error on acceptance
  std::map<TString, Double_t> *parameters = new std::map<TString, Double_t>;
  std::map<TString, Double_t>::iterator setbranch;
  for ( setparameters = sample_list.begin(); setparameters != sample_list.end(); setparameters++){
    if (setparameters->second==0) continue;
      parname = setparameters->first;
      parname_err = parname + string_err;
      parameters->insert( std::pair<TString, Double_t> (parname, 0.0));  
      parameters->insert( std::pair<TString, Double_t> (parname_err, 0.0));
  }

  // Declares file to store calculated acceptances to
  outfilename = "pdf_acceptances.root";
  outFile = new TFile(outfilename,"RECREATE"); 
  outTree = new TTree("Events","Events");
  for (setbranch = parameters->begin(); setbranch != parameters->end(); setbranch++){
    parname = setbranch->first;
    type = "/D";
    vartype = parname + type;
    outTree->Branch(parname, &(*parameters)[parname],vartype);
    cout << "Made branch " << parname << endl;
  }

  Int_t max_events=0;
  max_events = n_events;
  Double_t lhe=0.0, totalWeightGen=0.0, ntotal=0.0, fiducialWeightGen=0.0, acceptance=0.0, accerr=0.0;

  // Calculates acceptances for each LHE weight and associated error
  cout << "Calculating Acceptances " << endl;
  for (Int_t i=0; i<112; i++){
    cout << "PDF Number " << i << endl;
    for ( setparameters = sample_list.begin(); setparameters != sample_list.end(); setparameters++){ // Quantities from process_gen.root files
      max_events = n_events;
      lhe=0.0;
      totalWeightGen=0.0;
      fiducialWeightGen=0.0;
      acceptance=0.0;
      accerr=0.0;
      ntotal=0.0;

      parname = setparameters->first;
      parname_err = parname + string_err;
      if (setparameters->second == 1){
        cout << parname << endl;
        infilename = outputDir + TString("/") + parname + TString("_gen.root");
          infile = TFile::Open(infilename);         assert(infile);
          intree = (TTree*)infile->Get("Events"); assert(intree);
          intree->SetBranchAddress("glepq1",   &glepq1);                  // lepton1 charge 
          intree->SetBranchAddress("glepq1",   &glepq1);                  // lepton2 charge
          intree->SetBranchAddress("glep1",   &glep1);                    // lepton1 4-vector
          intree->SetBranchAddress("glep2",   &glep2);                    // lepton2 4-vector
          intree->SetBranchAddress("gvec",   &gvec);                      // boson 4-vector
          intree->SetBranchAddress("weightGen",   &weightGen);            // event weights
          intree->SetBranchAddress("lheweight",   &lheweight);            // lheweights
          intree->SetBranchAddress("fiducial",    &fiducial);             // if in fiducial region =1

        if (n_events==0 || n_events > intree->GetEntries()){
          max_events = intree->GetEntries();
        }

        for(Int_t ientry=0; ientry<max_events; ientry++) {
          intree->GetEntry(ientry);
          lhe = (*lheweight)[i];
          if (ientry==0){
            cout << "LHE Weight " << lhe << endl;
          }
          if (i==0 && lhe==0) {
            lhe=1;
          }
          totalWeightGen+=(weightGen*lhe);
          ntotal += weightGen/abs(weightGen);
          if (fiducial==1){
            fiducialWeightGen += (weightGen*lhe);
          }
        }
        acceptance = fiducialWeightGen/totalWeightGen;
        accerr = sqrt(acceptance*(1.0-acceptance)/ntotal);

        (*parameters)[parname] = acceptance;
        parname_err = setparameters->first + string_err;
        (*parameters)[parname_err] = accerr;
      }
    }

    for ( setparameters = sample_list.begin(); setparameters != sample_list.end(); setparameters++){ // Derived quantities
      if (setparameters->second == 2){
        TString ratio_name, num_name, den_name, ratio_name_err, num_name_err, den_name_err;

        ratio_name = setparameters->first;
        num_name = ratio_list[ratio_name].first;
        den_name = ratio_list[ratio_name].second;

        ratio_name_err = ratio_name + string_err;
        num_name_err = num_name + string_err;
        den_name_err = den_name + string_err;

        (*parameters)[ratio_name] = (*parameters)[num_name]/(*parameters)[den_name];
        (*parameters)[ratio_name_err] = (*parameters)[ratio_name]*sqrt(TMath::Power((*parameters)[num_name_err]/(*parameters)[num_name],2)+TMath::Power((*parameters)[den_name_err]/(*parameters)[den_name],2));
      }
    }
    outTree->Fill();
  }
  outFile->Write();
  outFile->Close(); 

  // Reads the file of parameter acceptances and errors that was just created
  infilename = "pdf_acceptances.root";
  infile = TFile::Open(infilename);         assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);
  for (setbranch = parameters->begin(); setbranch != parameters->end(); setbranch++){
    intree->SetBranchAddress(setbranch->first,  &(*parameters)[setbranch->first]);
  }

  cout << "Calculating Errors " << endl;
  ofstream txtfile;
  txtfile.open("pdf_acceptances.txt", ios::out);
  txtfile << "Process " << "Nominal_Acceptance " << "Nominal_Acceptance_Error " << "Mean_Acceptance " << "PDF_Uncertainty " << "Alphas_Acceptance " << "Alphas_Acceptance_Deviation " << endl;

  for ( setparameters = sample_list.begin(); setparameters != sample_list.end(); setparameters++){
    parname = setparameters->first;
    parname_err = parname + string_err;

    if (setparameters->second==0) continue;

    cout << parname << endl;

    Double_t n_pdfs = 0.0, nomacc = 0.0, nomacc_err = 0.0, mean = 0.0, err = 0.0, mean2 = 0.0, alphas_acc=0.0,  alphas_acc_err = 0.0;

      intree->GetEntry(0); // The nominal acceptance
      nomacc = (*parameters)[parname];
      nomacc_err = (*parameters)[parname_err];

    n_pdfs = 100;
    for(Int_t ientry=10; ientry<110; ientry++) { // entry 9 is the nominal value again, then 100 replicas
      intree->GetEntry(ientry);      
      mean += (*parameters)[parname]/n_pdfs;
      mean2 += TMath::Power((*parameters)[parname],2)/n_pdfs;
    }
    err = sqrt((n_pdfs/(n_pdfs-1))*(mean2 - TMath::Power(mean,2)));

    n_pdfs = 2; //actually 2
    for(Int_t ientry=110; ientry<112; ientry++) { // Alpha s up and down uncertainty
      intree->GetEntry(ientry);      
      alphas_acc += (*parameters)[parname]/n_pdfs;
    }
    alphas_acc_err = abs(nomacc - alphas_acc);

    txtfile << parname << " " << nomacc << " " << nomacc_err << " " << mean << " " << err << " " <<  alphas_acc << " " << alphas_acc_err << endl;
  }
  txtfile.close();

  cout << endl;
  cout << "  <> Output saved in " << "pdf_acceptances.txt" << endl;    
  cout << endl;  
      
  gBenchmark->Show("pdf_uncertainty_acceptance");

}