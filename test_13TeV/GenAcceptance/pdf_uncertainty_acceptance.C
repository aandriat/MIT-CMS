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

void pdf_uncertainty_acceptance(  const TString outputDir="."   // ntuple directory
            ){
  
  gBenchmark->Start("pdf_uncertainty_acceptance");

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

  for (Int_t i=8; i<108; i++){

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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
        if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      wpe = fiducialWeightGen/totalWeightGen;
      cout << wpe << endl;      
      delete intree;
      delete infile;

    cout << "wme" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("wme") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      wme = fiducialWeightGen/totalWeightGen;
      cout << wme << endl;      
      delete intree;
      delete infile;

    cout << "we" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("we") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      we = fiducialWeightGen/totalWeightGen;
      cout << we << endl;      
      delete intree;
      delete infile;

    cout << "zee" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("zee") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (glep2->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 ||(TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5)) && (TMath::Abs(glep2->Eta()) < 1.4442 || (TMath::Abs(glep2->Eta()) > 1.566 && TMath::Abs(glep2->Eta()) < 2.5)) && (gvec->M() > 60 && gvec->M() < 120)){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      zee = fiducialWeightGen/totalWeightGen;
      cout << zee << endl;      
      delete intree;
      delete infile;

    cout << "wpm" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("wpm") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      wpm = fiducialWeightGen/totalWeightGen;
      cout << wpm << endl;      
      delete intree;
      delete infile;

    cout << "wmm" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("wmm") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      wmm = fiducialWeightGen/totalWeightGen;
      cout << wmm << endl;      
      delete intree;
      delete infile;

    cout << "wm" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("wm") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      wm = fiducialWeightGen/totalWeightGen;
      cout << wm << endl;      
      delete intree;
      delete infile;

    cout << "zmm" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("zmm") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (glep2->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4) && (TMath::Abs(glep2->Eta()) < 2.4) && (gvec->M() > 60 && gvec->M() < 120)){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      zmm = fiducialWeightGen/totalWeightGen;
      cout << zmm << endl;      
      delete intree;
      delete infile;

    cout << "wpewme" << endl;
      wpewme = wpe/wme;
    cout << wpewme << endl;

    cout << "wpezee" << endl;
      wpezee = wpe/zee;
    cout << wpezee << endl;

    cout << "wmezee" << endl;
      wmezee = wme/zee;
    cout << wmezee << endl;

    cout << "wezee" << endl;
      wezee = we/zee;
    cout << wezee << endl;

    cout << "wpmwmm" << endl;
      wpmwmm = wpm/wmm;
    cout << wpmwmm << endl;

    cout << "wpmzmm" << endl;
      wpmzmm = wpm/zmm;
    cout << wpmzmm << endl;

    cout << "wmmzmm" << endl;
      wmmzmm = wmm/zmm;
    cout << wmmzmm << endl;

    cout << "wmzmm" << endl;
      wmzmm = wm/zmm;
    cout << wmzmm << endl;

    outTree->Fill();
  }
  outFile->Write();
  outFile->Close(); 

  // Read from input ntuple of signal events from flattened Bacon
  infilename = "pdf_acceptances.root";
  infile = TFile::Open(infilename);         assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);
  intree->SetBranchAddress("wpe",        &wpe );
  intree->SetBranchAddress("wme",        &wme );
  intree->SetBranchAddress("we",         &we  );
  intree->SetBranchAddress("zee",        &zee );
  intree->SetBranchAddress("wpm",        &wpm );
  intree->SetBranchAddress("wmm",        &wmm );
  intree->SetBranchAddress("wm",         &wm  );
  intree->SetBranchAddress("zmm",        &zmm );
  intree->SetBranchAddress("wpewme",        &wpewme   );
  intree->SetBranchAddress("wpezee",        &wpezee   );
  intree->SetBranchAddress("wmezee",        &wmezee   );
  intree->SetBranchAddress("wezee",         &wezee   );
  intree->SetBranchAddress("wpmwmm",        &wpmwmm   );
  intree->SetBranchAddress("wpmzmm",        &wpmzmm   );
  intree->SetBranchAddress("wmmzmm",        &wmmzmm   );
  intree->SetBranchAddress("wmzmm",         &wmzmm   );

  Double_t n_pdfs = intree->GetEntries();
  cout << "Number of pdfs" << endl;
  cout << n_pdfs << endl;

  Double_t wpe_mean=0.0, wpe_mean2=0.0, wpe_err=0.0;
  Double_t wme_mean=0.0, wme_mean2=0.0, wme_err=0.0;
  Double_t we_mean=0.0, we_mean2=0.0, we_err=0.0;
  Double_t zee_mean=0.0, zee_mean2=0.0, zee_err=0.0;
  Double_t wpm_mean=0.0, wpm_mean2=0.0, wpm_err=0.0;
  Double_t wmm_mean=0.0, wmm_mean2=0.0, wmm_err=0.0;
  Double_t wm_mean=0.0, wm_mean2=0.0, wm_err=0.0;
  Double_t zmm_mean=0.0, zmm_mean2=0.0, zmm_err=0.0;
  Double_t wpewme_mean=0.0, wpewme_mean2=0.0, wpewme_err=0.0;
  Double_t wpezee_mean=0.0, wpezee_mean2=0.0, wpezee_err=0.0;
  Double_t wmezee_mean=0.0, wmezee_mean2=0.0, wmezee_err=0.0;
  Double_t wezee_mean=0.0, wezee_mean2=0.0, wezee_err=0.0;
  Double_t wpmwmm_mean=0.0, wpmwmm_mean2=0.0, wpmwmm_err=0.0;
  Double_t wpmzmm_mean=0.0, wpmzmm_mean2=0.0, wpmzmm_err=0.0;
  Double_t wmmzmm_mean=0.0, wmmzmm_mean2=0.0, wmmzmm_err=0.0;
  Double_t wmzmm_mean=0.0, wmzmm_mean2=0.0, wmzmm_err=0.0;


  for(Int_t ientry=0; ientry<n_pdfs; ientry++) {
    intree->GetEntry(ientry);
    wpe_mean += wpe/n_pdfs;
    wme_mean += wme/n_pdfs;
    we_mean += we/n_pdfs;
    zee_mean += zee/n_pdfs;
    wpm_mean += wpm/n_pdfs;
    wmm_mean += wmm/n_pdfs;
    wm_mean += wm/n_pdfs;
    zmm_mean += zmm/n_pdfs;
    wpewme_mean += wpewme/n_pdfs;
    wpezee_mean += wpezee/n_pdfs;
    wmezee_mean += wmezee/n_pdfs;
    wezee_mean += wezee/n_pdfs;
    wpmwmm_mean += wpmwmm/n_pdfs;
    wpmzmm_mean += wpmzmm/n_pdfs;
    wmmzmm_mean += wmmzmm/n_pdfs;
    wmzmm_mean += wmzmm/n_pdfs;

    wpe_mean2 += TMath::Power(wpe,2)/n_pdfs;
    wme_mean2 += TMath::Power(wme,2)/n_pdfs;
    we_mean2 += TMath::Power(we,2)/n_pdfs;
    zee_mean2 += TMath::Power(zee,2)/n_pdfs;
    wpm_mean2 += TMath::Power(wpm,2)/n_pdfs;
    wmm_mean2 += TMath::Power(wmm,2)/n_pdfs;
    wm_mean2 += TMath::Power(wm,2)/n_pdfs;
    zmm_mean2 += TMath::Power(zmm,2)/n_pdfs;
    wpewme_mean2 += TMath::Power(wpewme,2)/n_pdfs;
    wpezee_mean2 += TMath::Power(wpezee,2)/n_pdfs;
    wmezee_mean2 += TMath::Power(wmezee,2)/n_pdfs;
    wezee_mean2 += TMath::Power(wezee,2)/n_pdfs;
    wpmwmm_mean2 += TMath::Power(wpmwmm,2)/n_pdfs;
    wpmzmm_mean2 += TMath::Power(wpmzmm,2)/n_pdfs;
    wmmzmm_mean2 += TMath::Power(wmmzmm,2)/n_pdfs;
    wmzmm_mean2 += TMath::Power(wmzmm,2)/n_pdfs;
  }

  wpe_err = sqrt((n_pdfs/(n_pdfs-1))*(wpe_mean2 - TMath::Power(wpe_mean,2)));
  wme_err = sqrt((n_pdfs/(n_pdfs-1))*(wme_mean2 - TMath::Power(wme_mean,2)));
  we_err = sqrt((n_pdfs/(n_pdfs-1))*(we_mean2 - TMath::Power(we_mean,2)));
  zee_err = sqrt((n_pdfs/(n_pdfs-1))*(zee_mean2 - TMath::Power(zee_mean,2)));
  wpm_err = sqrt((n_pdfs/(n_pdfs-1))*(wpm_mean2 - TMath::Power(wpm_mean,2)));
  wmm_err = sqrt((n_pdfs/(n_pdfs-1))*(wmm_mean2 - TMath::Power(wmm_mean,2)));
  wm_err = sqrt((n_pdfs/(n_pdfs-1))*(wm_mean2 - TMath::Power(wm_mean,2)));
  zmm_err = sqrt((n_pdfs/(n_pdfs-1))*(zmm_mean2 - TMath::Power(zmm_mean,2)));
  wpewme_err = sqrt((n_pdfs/(n_pdfs-1))*(wpewme_mean2 - TMath::Power(wpewme_mean,2)));
  wpezee_err = sqrt((n_pdfs/(n_pdfs-1))*(wpezee_mean2 - TMath::Power(wpezee_mean,2)));
  wmezee_err = sqrt((n_pdfs/(n_pdfs-1))*(wmezee_mean2 - TMath::Power(wmezee_mean,2)));
  wezee_err = sqrt((n_pdfs/(n_pdfs-1))*(wezee_mean2 - TMath::Power(wezee_mean,2)));
  wpmwmm_err = sqrt((n_pdfs/(n_pdfs-1))*(wpmwmm_mean2 - TMath::Power(wpmwmm_mean,2)));
  wpmzmm_err = sqrt((n_pdfs/(n_pdfs-1))*(wpmzmm_mean2 - TMath::Power(wpmzmm_mean,2)));
  wmmzmm_err = sqrt((n_pdfs/(n_pdfs-1))*(wmmzmm_mean2 - TMath::Power(wmmzmm_mean,2)));
  wmzmm_err = sqrt((n_pdfs/(n_pdfs-1))*(wmzmm_mean2 - TMath::Power(wmzmm_mean,2)));


    cout << "Process " << "Acceptance " << "Acceptance_Error " << endl;
      cout << "wpe" << " " << wpe_mean2 << " " << TMath::Power(wpe_mean,2) << " " << wpe_err << endl;
      cout << "wme" << " " << wme_mean2 << " " << TMath::Power(wme_mean,2) << " " << wme_err << endl;
      cout << "we" << " " << we_mean2 << " " << TMath::Power(we_mean,2) << " " << we_err << endl;
      cout << "zee" << " " << zee_mean2 << " " << TMath::Power(zee_mean,2) << " " << zee_err << endl;
      cout << "wpm" << " " << wpm_mean2 << " " << TMath::Power(wpm_mean,2) << " " << wpm_err << endl;
      cout << "wmm" << " " << wmm_mean2 << " " << TMath::Power(wmm_mean,2) << " " << wmm_err << endl;
      cout << "wm" << " " << wm_mean2 << " " << TMath::Power(wm_mean,2) << " " << wm_err << endl;
      cout << "zmm" << " " << zmm_mean2 << " " << TMath::Power(zmm_mean,2) << " " << zmm_err << endl;
      cout << "wpewme" << " " << wpewme_mean2 << " " << TMath::Power(wpewme_mean,2) << " " << wpewme_err << endl;
      cout << "wpezee" << " " << wpezee_mean2 << " " << TMath::Power(wpezee_mean,2) << " " << wpezee_err << endl;
      cout << "wmezee" << " " << wmezee_mean2 << " " << TMath::Power(wmezee_mean,2) << " " << wmezee_err << endl;
      cout << "wezee" << " " << wezee_mean2 << " " << TMath::Power(wezee_mean,2) << " " << wezee_err << endl;
      cout << "wpmwmm" << " " << wpmwmm_mean2 << " " << TMath::Power(wpmwmm_mean,2) << " " << wpmwmm_err << endl;
      cout << "wpmzmm" << " " << wpmzmm_mean2 << " " << TMath::Power(wpmzmm_mean,2) << " " << wpmzmm_err << endl;
      cout << "wmmzmm" << " " << wmmzmm_mean2 << " " << TMath::Power(wmmzmm_mean,2) << " " << wmmzmm_err << endl;
      cout << "wmzmm" << " " << wmzmm_mean2 << " " << TMath::Power(wmzmm_mean,2) << " " << wmzmm_err << endl;

  ofstream txtfile;
  txtfile.open("pdf_acceptances.txt", ios::app);
    txtfile << "Process " << "Acceptance " << "Acceptance_Error " << endl;
      txtfile << "wpe" << " " << wpe_mean << " " << wpe_err << endl;
      txtfile << "wme" << " " << wme_mean << " " << wme_err << endl;
      txtfile << "we" << " " << we_mean << " " << we_err << endl;
      txtfile << "zee" << " " << zee_mean << " " << zee_err << endl;
      txtfile << "wpm" << " " << wpm_mean << " " << wpm_err << endl;
      txtfile << "wmm" << " " << wmm_mean << " " << wmm_err << endl;
      txtfile << "wm" << " " << wm_mean << " " << wm_err << endl;
      txtfile << "zmm" << " " << zmm_mean << " " << zmm_err << endl;
      txtfile << "wpewme" << " " << wpewme_mean << " " << wpewme_err << endl;
      txtfile << "wpezee" << " " << wpezee_mean << " " << wpezee_err << endl;
      txtfile << "wmezee" << " " << wmezee_mean << " " << wmezee_err << endl;
      txtfile << "wezee" << " " << wezee_mean << " " << wezee_err << endl;
      txtfile << "wpmwmm" << " " << wpmwmm_mean << " " << wpmwmm_err << endl;
      txtfile << "wpmzmm" << " " << wpmzmm_mean << " " << wpmzmm_err << endl;
      txtfile << "wmmzmm" << " " << wmmzmm_mean << " " << wmmzmm_err << endl;
      txtfile << "wmzmm" << " " << wmzmm_mean << " " << wmzmm_err << endl;
  txtfile.close();

  // Alpha S Uncertainty Contribution

for (Int_t i=108; i<109; i++){

    cout << "Alpha_s " << i << endl;

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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
        if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      wpe = fiducialWeightGen/totalWeightGen;      
      delete intree;
      delete infile;

    cout << "wme" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("wme") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      wme = fiducialWeightGen/totalWeightGen;      
      delete intree;
      delete infile;

    cout << "we" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("we") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      we = fiducialWeightGen/totalWeightGen;      
      delete intree;
      delete infile;

    cout << "zee" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("zee") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (glep2->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 ||(TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5)) && (TMath::Abs(glep2->Eta()) < 1.4442 || (TMath::Abs(glep2->Eta()) > 1.566 && TMath::Abs(glep2->Eta()) < 2.5)) && (gvec->M() > 60 && gvec->M() < 120)){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      zee = fiducialWeightGen/totalWeightGen;      
      delete intree;
      delete infile;

    cout << "wpm" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("wpm") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      wpm = fiducialWeightGen/totalWeightGen;      
      delete intree;
      delete infile;

    cout << "wmm" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("wmm") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      wmm = fiducialWeightGen/totalWeightGen;      
      delete intree;
      delete infile;

    cout << "wm" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("wm") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      wm = fiducialWeightGen/totalWeightGen;      
      delete intree;
      delete infile;

    cout << "zmm" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("zmm") + TString("_gen.root");
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

      // Calculate total number of (weighted) events in sample
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
        totalWeightGen+=(weightGen*(*lheweight)[i]);
          if ((glep1->Pt() > 25) && (glep2->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4) && (TMath::Abs(glep2->Eta()) < 2.4) && (gvec->M() > 60 && gvec->M() < 120)){
          fiducialWeightGen += (weightGen*(*lheweight)[i]);
        }
      }
      zmm = fiducialWeightGen/totalWeightGen;      
      delete intree;
      delete infile;

    cout << "wpewme" << endl;
      wpewme = wpe/wme;

    cout << "wpezee" << endl;
      wpezee = wpe/zee;

    cout << "wmezee" << endl;
      wmezee = wme/zee;

    cout << "wezee" << endl;
      wezee = we/zee;

    cout << "wpmwmm" << endl;
      wpmwmm = wpm/wmm;

    cout << "wpmzmm" << endl;
      wpmzmm = wpm/zmm;

    cout << "wmmzmm" << endl;
      wmmzmm = wmm/zmm;

    cout << "wmzmm" << endl;
      wmzmm = wm/zmm;

    ofstream txtfile2;
      txtfile2.open("alphas_acceptance.txt", ios::app);
      txtfile2 << "Process " << "Acceptance " << "Acceptance_Error " << endl;
      txtfile2 << "wpe" << " " << wpe << " " << fabs(wpe - wpe_mean) << endl;
      txtfile2 << "wme" << " " << wme << " " << fabs(wme - wme_mean) << endl;
      txtfile2 << "we" << " " << we << " " << fabs(we - we_mean) << endl;
      txtfile2 << "zee" << " " << zee << " " << fabs(zee - zee_mean) << endl;
      txtfile2 << "wpm" << " " << wpm << " " << fabs(wpm - wpm_mean) << endl;
      txtfile2 << "wmm" << " " << wmm << " " << fabs(wmm - wmm_mean) << endl;
      txtfile2 << "wm" << " " << wm << " " << fabs(wm - wm_mean) << endl;
      txtfile2 << "zmm" << " " << zmm << " " << fabs(zmm - zmm_mean) << endl;
      txtfile2 << "wpewme" << " " << wpewme << " " << fabs(wpewme - wpewme_mean) << endl;
      txtfile2 << "wpezee" << " " << wpezee << " " << fabs(wpezee - wpezee_mean) << endl;
      txtfile2 << "wmezee" << " " << wmezee << " " << fabs(wmezee - wmezee_mean) << endl;
      txtfile2 << "wezee" << " " << wezee << " " << fabs(wezee - wezee_mean) << endl;
      txtfile2 << "wpmwmm" << " " << wpmwmm << " " << fabs(wpmwmm - wpmwmm_mean) << endl;
      txtfile2 << "wpmzmm" << " " << wpmzmm << " " << fabs(wpmzmm - wpmzmm_mean) << endl;
      txtfile2 << "wmmzmm" << " " << wmmzmm << " " << fabs(wmmzmm - wmmzmm_mean) << endl;
      txtfile2 << "wmzmm" << " " << wmzmm << " " << fabs(wmzmm - wmzmm_mean) << endl;
    txtfile2.close();
  }


// End Alpha S Contribution





  cout << endl;
  cout << "  <> Output saved in " << "pdf_acceptances.txt" << endl;    
  cout << endl;  
      
  gBenchmark->Show("pdf_uncertainty_acceptance");

}


