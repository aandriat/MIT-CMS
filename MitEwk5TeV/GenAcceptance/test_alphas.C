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

void test_alphas(  const TString outputDir="."   // ntuple directory
            ){
  
  gBenchmark->Start("test_alphas");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // Declare input ntuple variables
  Double_t weightGen;  
  TLorentzVector *glep1=0, *glep2=0, *gvec=0;
  Int_t glepq1, glepq2;
  std::vector<float> *lheweight = new std::vector<float>();

  // Declare output ntuple variables
  Double_t wpe, wme, we, zee, wpm, wmm, wm, zmm;
  Double_t wpewme, wpezee, wmezee, wezee;
  Double_t wpmwmm, wpmzmm, wmmzmm, wmzmm; 

  // Set up output tree of acceptances
  TString outfilename = "pdf_acceptances.root";
  TFile *outFile = new TFile(outfilename,"RECREATE"); 
  TTree *outTree = new TTree("Events","Events");
  outTree->Branch("wpe",        &wpe,         "wpe/D");
  outTree->Branch("wme",        &wme,         "wme/D");
  outTree->Branch("we",         &we,           "we/D");
  outTree->Branch("zee",        &zee,         "zee/D");
  outTree->Branch("wpm",        &wpm,         "wpm/D");
  outTree->Branch("wmm",        &wmm,         "wmm/D");
  outTree->Branch("wm",         &wm,           "wm/D");
  outTree->Branch("zmm",        &zmm,         "zmm/D");
  outTree->Branch("wpewme",        &wpewme,       "wpewme/D");
  outTree->Branch("wpezee",        &wpezee,       "wpezee/D");
  outTree->Branch("wmezee",        &wmezee,       "wmezee/D");
  outTree->Branch("wezee",         &wezee,        "wezee/D");
  outTree->Branch("wpmwmm",        &wpmwmm,       "wpmwmm/D");
  outTree->Branch("wpmzmm",        &wpmzmm,       "wpmzmm/D");
  outTree->Branch("wmmzmm",        &wmmzmm,       "wmmzmm/D");
  outTree->Branch("wmzmm",         &wmzmm,        "wmzmm/D");

  TFile *infile = 0;
  TTree *intree = 0;
  TString infilename;
  Double_t fiducialWeightGen;
  Double_t totalWeightGen;

  for (Int_t i=0; i<110; i++){

    cout << "LHE Weight " << i << endl;

    cout << "wpe" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("wpe") + TString("_gen.root");
      infile = TFile::Open(infilename);         assert(infile);
      intree = (TTree*)infile->Get("Events"); assert(intree);
      intree->SetBranchAddress("glepq1",   &glepq1);                  // lepton1 charge 
      intree->SetBranchAddress("glepq1",   &glepq1);                  // lepton2 charge
      intree->SetBranchAddress("glep1",   &glep1);                    // lepton1 4-vector
      intree->SetBranchAddress("glep2",   &glep2);                    // lepton2 4-vector
      intree->SetBranchAddress("gvec",   &gvec);                      // boson 4-vector
      intree->SetBranchAddress("weightGen",   &weightGen);            // event weights
      intree->SetBranchAddress("lheweight",   &lheweight);            // lheweights
      
      // Declare variables used to store weight information
      fiducialWeightGen=0;
      totalWeightGen=0;

      intree->GetEntry(0);
      cout << "Alpha S weight " << (*lheweight)[i] << endl; 

      // // Calculate total number of (weighted) events in sample
      // for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      //   intree->GetEntry(ientry);
      //   totalWeightGen+=(weightGen*(*lheweight)[i]);
      //   if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
      //     fiducialWeightGen += (weightGen*(*lheweight)[i]);
      //   }
      // }
      // wpe = fiducialWeightGen/totalWeightGen;      
      // delete intree;
      // delete infile;
  }

  cout << endl;
  cout << "  <> Output saved in " << "pdf_acceptances.txt" << endl;    
  cout << endl;  
      
  gBenchmark->Show("test_alphas");

}


