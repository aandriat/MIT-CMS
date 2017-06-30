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

void pdf_uncertainty_acceptance(  const TString outputDir=".",   // ntuple directory
            ){
  
  gBenchmark->Start("pdf_uncertainty_acceptance");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // Declare input ntuple variables
  Double_t weightGen, totalWeightGen;  
  TLorentzVector *glep1=0, *glep2=0, *gvec=0;
  Int_t glepq1, glepq2;
  std::vector<float> *lheweight = new std::vector<float>();

  // Declare output ntuple variables
  Double_t wpe, wme, w, zee, wpm, wmm, w, zmm;
  Double_t wpe-wme, wpe-zee, wme-zee, we-zee;
  Double_t wpm-wmm, wpm-zmm, wmm-zmm, wm-zmm; 

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
  outTree->Branch("wpe-wme",        &wpe-wme,       "wpe-wme/D");
  outTree->Branch("wpe-zee",        &wpe-zee,       "wpe-zee/D");
  outTree->Branch("wme-zee",        &wme-zee,       "wme-zee/D");
  outTree->Branch("we-zee",         &we-zee,        "we-zee/D");
  outTree->Branch("wpm-wmm",        &wpm-wmm,       "wpm-wmm/D");
  outTree->Branch("wmp-zmm",        &wmp-zmm,       "wmp-zmm/D");
  outTree->Branch("wmm-zmm",        &wmm-zmm,       "wmm-zmm/D");
  outTree->Branch("wm-zmm",         &wm-zmm,        "wm-zmm/D");

  TFile *infile = 0;
  TTree *intree = 0;
  TString infile;
  Double_t fiducialWeightGen;
  Double_t totalWeightGen;

  for (Int_t i=9; i<109; i++){

    cout << "LHE Weight " << i << endl;

    cout << "wpe" << endl;
      // Read from input ntuple of signal events from flattened Bacon
      infilename = outputDir + TString("/") + TString("wpe") + TString("_gen.root");
      *infile = TFile::Open(infilename);         assert(infile);
      *intree = (TTree*)infile->Get("Events"); assert(intree);
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
      *infile = TFile::Open(infilename);         assert(infile);
      *intree = (TTree*)infile->Get("Events"); assert(intree);
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
      *infile = TFile::Open(infilename);         assert(infile);
      *intree = (TTree*)infile->Get("Events"); assert(intree);
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
      *infile = TFile::Open(infilename);         assert(infile);
      *intree = (TTree*)infile->Get("Events"); assert(intree);
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
      *infile = TFile::Open(infilename);         assert(infile);
      *intree = (TTree*)infile->Get("Events"); assert(intree);
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
      *infile = TFile::Open(infilename);         assert(infile);
      *intree = (TTree*)infile->Get("Events"); assert(intree);
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
      *infile = TFile::Open(infilename);         assert(infile);
      *intree = (TTree*)infile->Get("Events"); assert(intree);
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
      *infile = TFile::Open(infilename);         assert(infile);
      *intree = (TTree*)infile->Get("Events"); assert(intree);
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

    cout << "wpe-wme" << endl;
      wpe-wme = wpe/wme;

    cout << "wpe-zee" << endl;
      wpe-zee = wpe/zee;

    cout << "wme-zee" << endl;
      wme-zee = wme/zee;

    cout << "we-zee" << endl;
      we-zee = we/zee;

    cout << "wpm-wmm" << endl;
      wpm-wmm = wpm/wmm;

    cout << "wpm-zmm" << endl;
      wpm-zmm = wpm/zmm;

    cout << "wmm-zmm" << endl;
      wmm-zmm = wmm/zmm;

    cout << "wm-zmm" << endl;
      wm-zmm = wm/zmm;

    outTree->Fill();
  }
  outFile->Write();
  outFile->Close(); 

  // Read from input ntuple of signal events from flattened Bacon
  infilename = "pdf_acceptances.root";
  *infile = TFile::Open(infilename);         assert(infile);
  *intree = (TTree*)infile->Get("Events"); assert(intree);
  intree->SetBranchAddress("wpe",        &wpe );
  intree->SetBranchAddress("wme",        &wme );
  intree->SetBranchAddress("we",         &we  );
  intree->SetBranchAddress("zee",        &zee );
  intree->SetBranchAddress("wpm",        &wpm );
  intree->SetBranchAddress("wmm",        &wmm );
  intree->SetBranchAddress("wm",         &wm  );
  intree->SetBranchAddress("zmm",        &zmm );
  intree->SetBranchAddress("wpe-wme",        &wpe-wme   );
  intree->SetBranchAddress("wpe-zee",        &wpe-zee   );
  intree->SetBranchAddress("wme-zee",        &wme-zee   );
  intree->SetBranchAddress("we-zee",         &we-zee   );
  intree->SetBranchAddress("wpm-wmm",        &wpm-wmm   );
  intree->SetBranchAddress("wpm-zmm",        &wpm-zmm   );
  intree->SetBranchAddress("wmm-zmm",        &wmm-zmm   );
  intree->SetBranchAddress("wm-zmm",         &wm-zmm   );

  Int_t n_pdfs = intree->GetEntries();

  Double_t wpe_mean, wpe_mean2, wpe_err;
  Double_t wme_mean, wme_mean2, wme_err;
  Double_t we_mean, we_mean2, we_err;
  Double_t zee_mean, zee_mean2, zee_err;
  Double_t wpm_mean, wpm_mean2, wpm_err;
  Double_t wmm_mean, wmm_mean2, wmm_err;
  Double_t wm_mean, wm_mean2, wm_err;
  Double_t zmm_mean, zmm_mean2, zmm_err;
  Double_t wpe-wme_mean, wpe-wme_mean2, wpe-wme_err;
  Double_t wpe-zee_mean, wpe-zee_mean2, wpe-zee_err;
  Double_t wme-zee_mean, wme-zee_mean2, wme-zee_err;
  Double_t we-zee_mean, we-zee_mean2, we-zee_err;
  Double_t wpm-wmm_mean, wpm-wmm_mean2, wpm-wmm_err;
  Double_t wpm-zmm_mean, wpm-zmm_mean2, wpm-zmm_err;
  Double_t wmm-zmm_mean, wmm-zmm_mean2, wmm-zmm_err;
  Double_t wm-zmm_mean, wm-zmm_mean2, wm-zmm_err;


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
    wpe-wme_mean += wpe-wme/n_pdfs;
    wpe-zee_mean += wpe-zee/n_pdfs;
    wme-zee_mean += wme-zee/n_pdfs;
    we-zee_mean += we-zee/n_pdfs;
    wpm-wmm_mean += wpm-wmm/n_pdfs;
    wpm-zmm_mean += wpm-zmm/n_pdfs;
    wmm-zmm_mean += wmm-zmm/n_pdfs;
    wm-zmm_mean += wm-zmm/n_pdfs;

    wpe_mean2 += power(wpe,2)/n_pdfs;
    wme_mean2 += power(wme,2)/n_pdfs;
    we_mean2 += power(we,2)/n_pdfs;
    zee_mean2 += power(zee,2)/n_pdfs;
    wpm_mean2 += power(wpm,2)/n_pdfs;
    wmm_mean2 += power(wmm,2)/n_pdfs;
    wm_mean2 += power(wm,2)/n_pdfs;
    zmm_mean2 += power(zmm,2)/n_pdfs;
    wpe-wme_mean2 += power(wpe-wme,2)/n_pdfs;
    wpe-zee_mean2 += power(wpe-zee,2)/n_pdfs;
    wme-zee_mean2 += power(wme-zee,2)/n_pdfs;
    we-zee_mean2 += power(we-zee,2)/n_pdfs;
    wpm-wmm_mean2 += power(wpm-wmm,2)/n_pdfs;
    wpm-zmm_mean2 += power(wpm-zmm,2)/n_pdfs;
    wmm-zmm_mean2 += power(wmm-zmm,2)/n_pdfs;
    wm-zmm_mean2 += power(wm-zmm,2)/n_pdfs;
  }

  wpe_err = sqrt((n_pdfs/(n_pdfs-1))*(wpe_mean2-power(wpe_mean,2)));
  wme_err = sqrt((n_pdfs/(n_pdfs-1))*(wme_mean2-power(wme_mean,2)));
  we_err = sqrt((n_pdfs/(n_pdfs-1))*(we_mean2-power(we_mean,2)));
  zee_err = sqrt((n_pdfs/(n_pdfs-1))*(zee_mean2-power(zee_mean,2)));
  wpm_err = sqrt((n_pdfs/(n_pdfs-1))*(wpm_mean2-power(wpm_mean,2)));
  wmm_err = sqrt((n_pdfs/(n_pdfs-1))*(wmm_mean2-power(wmm_mean,2)));
  wm_err = sqrt((n_pdfs/(n_pdfs-1))*(wm_mean2-power(wm_mean,2)));
  zmm_err = sqrt((n_pdfs/(n_pdfs-1))*(zmm_mean2-power(zmm_mean,2)));
  wpe-wme_err = sqrt((n_pdfs/(n_pdfs-1))*(wpe-wme_mean2-power(wpe-wme_mean,2)));
  wpe-zee_err = sqrt((n_pdfs/(n_pdfs-1))*(wpe-zee_mean2-power(wpe-zee_mean,2)));
  wme-zee_err = sqrt((n_pdfs/(n_pdfs-1))*(wme-zee_mean2-power(wme-zee_mean,2)));
  we-zee_err = sqrt((n_pdfs/(n_pdfs-1))*(we-zee_mean2-power(we-zee_mean,2)));
  wpm-wmm_err = sqrt((n_pdfs/(n_pdfs-1))*(wpm-wmm_mean2-power(wpm-wmm_mean,2)));
  wpm-zmm_err = sqrt((n_pdfs/(n_pdfs-1))*(wpm-zmm_mean2-power(wpm-zmm_mean,2)));
  wmm-zmm_err = sqrt((n_pdfs/(n_pdfs-1))*(wmm-zmm_mean2-power(wmm-zmm_mean,2)));
  wm-zmm_err = sqrt((n_pdfs/(n_pdfs-1))*(wm-zmm_mean2-power(wm-zmm_mean,2)));

  ofstream txtfile;
  txtfile.open("pdf_acceptances.txt", ios::app);
      txtfile << "*" << endl;
      txtfile << "For " << "wpe"  << endl;
      txtfile << "Acceptance = " << wpe_mean << endl;
      txtfile << "Acceptance Error = " << wpe_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wme"  << endl;
      txtfile << "Acceptance = " << wme_mean << endl;
      txtfile << "Acceptance Error = " << wme_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "we"  << endl;
      txtfile << "Acceptance = " << we_mean << endl;
      txtfile << "Acceptance Error = " << we_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "zee"  << endl;
      txtfile << "Acceptance = " << zee_mean << endl;
      txtfile << "Acceptance Error = " << zee_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wpm"  << endl;
      txtfile << "Acceptance = " << wpm_mean << endl;
      txtfile << "Acceptance Error = " << wpm_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wmm"  << endl;
      txtfile << "Acceptance = " << wmm_mean << endl;
      txtfile << "Acceptance Error = " << wmm_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wm"  << endl;
      txtfile << "Acceptance = " << wm_mean << endl;
      txtfile << "Acceptance Error = " << wm_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "zmm"  << endl;
      txtfile << "Acceptance = " << zmm_mean << endl;
      txtfile << "Acceptance Error = " << zmm_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wpe-wme"  << endl;
      txtfile << "Acceptance = " << wpe-wme_mean << endl;
      txtfile << "Acceptance Error = " << wpe-wme_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wpe-zee"  << endl;
      txtfile << "Acceptance = " << wpe-zee_mean << endl;
      txtfile << "Acceptance Error = " << wpe-zee_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wme-zee"  << endl;
      txtfile << "Acceptance = " << wme-zee_mean << endl;
      txtfile << "Acceptance Error = " << wme-zee_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "we-zee"  << endl;
      txtfile << "Acceptance = " << we-zee_mean << endl;
      txtfile << "Acceptance Error = " << we-zee_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wpm-wmm"  << endl;
      txtfile << "Acceptance = " << wpm-wmm_mean << endl;
      txtfile << "Acceptance Error = " << wpm-wmm_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wpm-zmm"  << endl;
      txtfile << "Acceptance = " << wpm-zmm_mean << endl;
      txtfile << "Acceptance Error = " << wpm-zmm_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wmm-zmm"  << endl;
      txtfile << "Acceptance = " << wmm-zmm_mean << endl;
      txtfile << "Acceptance Error = " << wmm-zmm_err << endl;

      txtfile << "*" << endl;
      txtfile << "For " << "wm-zmm"  << endl;
      txtfile << "Acceptance = " << wm-zmm_mean << endl;
      txtfile << "Acceptance Error = " << wm-zmm_err << endl;
  txtfile.close();

  cout << endl;
  cout << "  <> Output saved in " << "pdf_acceptances.txt" << endl;    
  cout << endl;  
      
  gBenchmark->Show("pdf_uncertainty_acceptance");

}


