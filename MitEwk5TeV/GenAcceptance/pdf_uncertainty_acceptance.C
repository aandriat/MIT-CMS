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

void pdf_uncertainty_acceptance(const TString conf="acceptance.conf", // input file
                  const TString outputDir="ntuples",   // ntuple directory
                  const Int_t n_events=0 //Number of events, 0 for all
            ){
  
  gBenchmark->Start("pdf_uncertainty_acceptance");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  confParse(conf, snamev, samplev); // parse .conf file

  // Declare ntuple variables
  Double_t weightGen, totalWeightGen;  
  TLorentzVector *glep1=0, *glep2=0, *gvec=0;
  Int_t glepq1, glepq2;
  std::vector<float> *lheweight = new std::vector<float>();

  // Print the samples and files that will be read from
  for(Int_t isam=0; isam<samplev.size(); isam++) {
    CSample* samp = samplev[isam];
    cout << "Sample name: " << endl;
    cout << snamev[isam] << endl;
    const Int_t nfiles = samp->fnamev.size();
    for(Int_t ifile=0; ifile<nfiles; ifile++) {
      cout <<"File name: " << endl;  
      cout << samp->fnamev[ifile] << endl; 
    }
  }

  // loop over samples
  for(Int_t isam=0; isam<samplev.size(); isam++) {

    Double_t mean_acc, acc_error; // Output variables of calculated values

    CSample* samp = samplev[isam]; // Read sample

    // Read from input ntuple of signal events from flattened Bacon
    TString infilename = outputDir + TString("/") + snamev[isam] + TString("_selected.root");
    TFile *infile = TFile::Open(infilename);         assert(infile);
    TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
    intree->SetBranchAddress("glepq1",   &glepq1);                  // lepton1 charge 
    intree->SetBranchAddress("glepq1",   &glepq1);                  // lepton2 charge
    intree->SetBranchAddress("glep1",   &glep1);                    // lepton1 4-vector
    intree->SetBranchAddress("glep2",   &glep2);                    // lepton2 4-vector
    intree->SetBranchAddress("gvec",   &gvec);                      // boson 4-vector
    intree->SetBranchAddress("weightGen",   &weightGen);            // event weights
    intree->SetBranchAddress("lheweight",   &lheweight);            // lheweights
    intree->SetBranchAddress("totalWeightGen",   &totalWeightGen);  // totalWeight
    
    cout << "Processing " << snamev[isam] << endl;

    Double_t inclusive_weight=0;
    intree->GetEntry(0);
    inclusive_weight = totalWeightGen;

    Int_t max_events = n_events;
    if (n_events==0){
      max_events = intree->GetEntries();
    }

    // Calculates mean acceptance
    Double_t numpdfs = 100;
    Double_t sum_acc=0;
    for (Int_t i=9; i<109; i++){
      Double_t acceptance=0, fiducial_weight=0, pdftotalWeightGen=0, temp_lheweight=0;
      for(Int_t ientry=0; ientry<max_events; ientry++) {
        if(ientry%1000000==0) cout << "Processing PDF " << i << ", event " << ientry << ". " << (double)ientry/(double)max_events*100 << " percent done with this file." << endl;
        Double_t pdfweight=0;
        intree->GetEntry(ientry);
        temp_lheweight = (*lheweight)[i];
        pdfweight = weightGen;
        pdfweight *= temp_lheweight;
        fiducial_weight += pdfweight;
      }
      pdftotalWeightGen = inclusive_weight*temp_lheweight;
      acceptance = fiducial_weight/pdftotalWeightGen;
      sum_acc += acceptance;
    }
    mean_acc = (1/numpdfs) * sum_acc; //mean acceptance

    // Calculates error on the acceptance
    Double_t sum2_acc_mean=0;
    for (Int_t i=9; i<109; i++){
      Double_t acceptance=0, fiducial_weight=0, pdftotalWeightGen=0, temp_lheweight=0;
      for(Int_t ientry=0; ientry<max_events; ientry++) {
        if(ientry%1000000==0) cout << "Processing PDF " << i << ", event " << ientry << ". " << (double)ientry/(double)max_events*100 << " percent done with this file." << endl;
        Double_t pdfweight=0;
        intree->GetEntry(ientry);
        temp_lheweight = (*lheweight)[i];
        pdfweight = weightGen;
        pdfweight *= temp_lheweight;
        fiducial_weight += pdfweight;
      }
      pdftotalWeightGen = inclusive_weight*(*lheweight)[i];
      acceptance = fiducial_weight/pdftotalWeightGen;
      acceptance -= mean_acc;
      acceptance *= acceptance;
      sum2_acc_mean+=acceptance;
    }
    acc_error = sqrt(1/(numpdfs-1) * sum2_acc_mean); // Error on the acceptance
    // Prints pdf acceptance
    cout << "Mean Acceptance = " << mean_acc << endl;
    cout << "PDF Acceptance Error = " << acc_error << endl;

    // Saves pdf acceptance as a text file
    ofstream txtfile;
    txtfile.open("pdf_acceptances.txt", ios::app);
    txtfile << "*" << endl;
    txtfile << "For " << snamev[isam] << endl;
    txtfile << "Mean Acceptance = " << mean_acc << endl;
    txtfile << "PDF Acceptance Error = " << acc_error << endl;
    txtfile.close();
  }

  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("pdf_uncertainty_acceptance");

}