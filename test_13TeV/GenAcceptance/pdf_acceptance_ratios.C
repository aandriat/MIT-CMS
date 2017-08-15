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

void pdf_acceptance_ratios(){
  
  gBenchmark->Start("pdf_acceptance_ratios");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  TString string_err = "_err";

  std::map<TString, Int_t> sample_list = { // Processes for which to calculate acceptances: 0 to skip, 1 if gen process with a process_gen.root file, 2 if derived quantity
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
  std::map<TString, Int_t>::iterator setparameters; // Iterator over sample_list

  // Declare 13 TeV input file
  TString infilename_13;
  TFile *infile_13=0;
  TTree *intree_13=0;

  // Delcare 5 TeV input file
  TString infilename_5;
  TFile *infile_5=0;
  TTree *intree_5=0;

  // Delcare output file of acceptance ratios
  TString outfilename;
  TFile *outFile=0;
  TTree *outTree=0;

  // Useful strings
  TString parname;
  TString parname_err;
  TString type, vartype;


  // Makes map of parameters with names and blank doubles to store acceptance values and errors
  std::map<TString, Double_t> *parameters = new std::map<TString, Double_t>;
  std::map<TString, Double_t>::iterator setbranch;
  for ( setparameters = sample_list.begin(); setparameters != sample_list.end(); setparameters++){
    if (setparameters->second==0) continue;
      parname = setparameters->first;
      parname_err = parname + string_err;
      parameters->insert( std::pair<TString, Double_t> (parname, 0.0));  
      parameters->insert( std::pair<TString, Double_t> (parname_err, 0.0));
  }

  // Declare ratio output file with branch for each parameter and error
  outfilename = "pdf_acceptance_ratios.root";
  outFile = new TFile(outfilename,"RECREATE"); 
  outTree = new TTree("Events","Events");
  for (setbranch = parameters->begin(); setbranch != parameters->end(); setbranch++){
    parname = setbranch->first;
    type = "/D";
    vartype = parname + type;
    outTree->Branch(parname, &(*parameters)[parname],vartype);
  }

  // Declares the 13 TeV input file with acceptances and acceptance errors, should already exist from previous runs of pdf_uncertainty_acceptance.C
  std::map<TString, Double_t> *parameters_13 = new std::map<TString, Double_t>;
  for ( setparameters = sample_list.begin(); setparameters != sample_list.end(); setparameters++){
    if (setparameters->second==0) continue;
      parname = setparameters->first;
      parname_err = parname + string_err;
      parameters_13->insert( std::pair<TString, Double_t> (parname, 0.0));  
      parameters_13->insert( std::pair<TString, Double_t> (parname_err, 0.0));
  }
  infilename_13 = "/afs/cern.ch/work/a/aandriat/public/wz_analysis/CMSSW_7_6_3_patch2/src/MIT-CMS/test_13TeV/GenAcceptance/acceptances/Powheg_nominal_13TeV_acceptance.root";
  infile_13 = TFile::Open(infilename_13);         assert(infile_13);
  intree_13 = (TTree*)infile_13->Get("Events"); assert(intree_13);
  for (setbranch = parameters_13->begin(); setbranch != parameters_13->end(); setbranch++){
    intree_13->SetBranchAddress(setbranch->first,  &(*parameters_13)[setbranch->first]);
  }

  // Declares the 5 TeV input file with acceptances and acceptance errors, should already exist from previous runs of pdf_uncertainty_acceptance.C
  std::map<TString, Double_t> *parameters_5 = new std::map<TString, Double_t>;
  for ( setparameters = sample_list.begin(); setparameters != sample_list.end(); setparameters++){
    if (setparameters->second==0) continue;
      parname = setparameters->first;
      parname_err = parname + string_err;
      parameters_5->insert( std::pair<TString, Double_t> (parname, 0.0));  
      parameters_5->insert( std::pair<TString, Double_t> (parname_err, 0.0));
  }
  infilename_5 = "/afs/cern.ch/work/a/aandriat/public/wz_analysis/CMSSW_7_6_3_patch2/src/MIT-CMS/test_13TeV/GenAcceptance/acceptances/Powheg_nominal_5TeV_acceptance.root";
  infile_5 = TFile::Open(infilename_5);         assert(infile_5);
  intree_5 = (TTree*)infile_5->Get("Events"); assert(intree_5);
  for (setbranch = parameters_5->begin(); setbranch != parameters_5->end(); setbranch++){
    intree_5->SetBranchAddress(setbranch->first,  &(*parameters_5)[setbranch->first]);
  }

  // Makes sure each file has same number of replica PDFs, sets number of PDFs
  Double_t n_pdfs_13 = intree_13->GetEntries();
  Double_t n_pdfs_5 = intree_5->GetEntries();
  if (n_pdfs_13!=n_pdfs_5) cout << "Numbers of PDFs across energies do not match!!!" << endl;
  Double_t n_pdfs = intree_5->GetEntries();
  cout << "Number of pdfs: " << n_pdfs << endl;

  // Calculates ratio of acceptance of 5 TeV / 13 TeV, acceptance error given by sum in quadrature of individual errors.
  cout << "Calculating Acceptances " << endl;
  for (Int_t i=0; i<n_pdfs; i++){
    cout << "PDF Number " << i << endl;
    for ( setparameters = sample_list.begin(); setparameters != sample_list.end(); setparameters++){
      (*parameters)[parname] = (*parameters_5)[parname]/(*parameters_13)[parname];
      (*parameters)[parname_err] = (*parameters)[parname]*sqrt(TMath::Power((*parameters_5)[parname_err]/(*parameters_5)[parname],2)+TMath::Power((*parameters_13)[parname_err]/(*parameters_13)[parname],2));
    }
    outTree->Fill(); // Writes acceptance ratio to tree
  }
  outFile->Write();
  outFile->Close(); 

  // Reads the acceptance ratio tree that was just made
  infilename = "pdf_acceptance_ratios.root";
  infile = TFile::Open(infilename);         assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);
  for (setbranch = parameters->begin(); setbranch != parameters->end(); setbranch++){
    intree->SetBranchAddress(setbranch->first,  &(*parameters)[setbranch->first]);
  }

  // Calculates various quantities of interest
  cout << "Calculating Errors " << endl;
  ofstream txtfile;
  txtfile.open("pdf_acceptance_ratios.txt", ios::app);
  txtfile << "Process " << "Nominal_Acceptance " << "Nominal_Acceptance_Error " << "Mean_Acceptance " << "PDF_Uncertainty " << "Alphas_Acceptance " << "Alphas_Acceptance_Deviation " << endl;

  for ( setparameters = sample_list.begin(); setparameters != sample_list.end(); setparameters++){
    parname = setparameters->first;
    parname_err = parname + string_err;

    cout << parname << endl;

    Double_t n_pdfs = 0.0, nomacc = 0.0, nomacc_err = 0.0, mean = 0.0, err = 0.0, mean2 = 0.0, alphas_acc=0.0,  alphas_acc_err = 0.0;

      intree->GetEntry(0);
      nomacc = (*parameters)[parname];
      nomacc_err = (*parameters)[parname_err];

    n_pdfs = 100;
    for(Int_t ientry=9; ientry<109; ientry++) {
      intree->GetEntry(ientry);      
      mean += (*parameters)[parname]/n_pdfs;
      mean2 += TMath::Power((*parameters)[parname],2)/n_pdfs;
    }
    err = sqrt((n_pdfs/(n_pdfs-1))*(mean2 - TMath::Power(mean,2)));

    n_pdfs = 1; //actually 2
    for(Int_t ientry=109; ientry<110; ientry++) { // Actually 111 but Alphas_up is broken
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
      
  gBenchmark->Show("pdf_acceptance_ratios");
  
}