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

void acceptance(const TString conf="acceptance.conf", // input file
                  const TString outputDir=".",   // ntuple directory
                  const Int_t n_events=0 //Number of events, 0 for all
            ){
  
  gBenchmark->Start("acceptance");

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

    CSample* samp = samplev[isam]; // Read sample

    // Set up output ntuple of events in fiducial region
    TString outfilename = outputDir + TString("/") + snamev[isam] + TString("_selected.root");
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");
    outTree->Branch("glepq1",     &glepq1,     "glepq1/I");         // lepton1 charge
    outTree->Branch("glepq2",     &glepq2,     "glepq2/I");         // lepton2 charge
    outTree->Branch("glep1",       "TLorentzVector",  &glep1);      // lepton1 4-vector
    outTree->Branch("glep2",       "TLorentzVector",  &glep2);      // lepton2 4-vector
    outTree->Branch("gvec",       "TLorentzVector",  &gvec);        // boson 4-vector
    outTree->Branch("weightGen",   &weightGen,   "weightGen/D");    // event weight (MC)
    outTree->Branch("totalWeightGen",   &totalWeightGen,   "totalWeightGen/D");    // total event weight (MC) (same for all events)
    outTree->Branch("lheweight",  &lheweight);                      // lheweights vector (same for all events)

    // Read from input ntuple of signal events from flattened Bacon
    TString infilename = outputDir + TString("/") + snamev[isam] + TString("_gen.root");
    TFile *infile = TFile::Open(infilename);         assert(infile);
    TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
    intree->SetBranchAddress("glepq1",   &glepq1);                  // lepton1 charge 
    intree->SetBranchAddress("glepq1",   &glepq1);                  // lepton2 charge
    intree->SetBranchAddress("glep1",   &glep1);                    // lepton1 4-vector
    intree->SetBranchAddress("glep2",   &glep2);                    // lepton2 4-vector
    intree->SetBranchAddress("gvec",   &gvec);                      // boson 4-vector
    intree->SetBranchAddress("weightGen",   &weightGen);            // event weights
    intree->SetBranchAddress("lheweight",   &lheweight);            // lheweights
    

    cout << "Processing " << snamev[isam] << endl;

    // Declare variables used to store weight information
    Double_t fiducialWeightGen=0;
    Double_t ntotal=0;
    totalWeightGen=0;

    // Set number of events to loop through
    Int_t max_events = n_events;
    if (n_events==0){
      max_events = intree->GetEntries();
    }

    // Calculate total number of (weighted) events in sample
    for(Int_t ientry=0; ientry<max_events; ientry++) {
      intree->GetEntry(ientry);
        totalWeightGen+=weightGen;
        ntotal += weightGen/abs(weightGen);
    }

    Bool_t pass_selection; // flag to accept only fiducial events

    // Loop through events
    for(Int_t ientry=0; ientry<max_events; ientry++) {

      pass_selection=kFALSE; // reset flag

      intree->GetEntry(ientry); // Store variables from input tree
      if(ientry%10000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;

      if (snamev[isam]=="zmm"){
        if ((glep1->Pt() > 25) && (glep2->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4) && (TMath::Abs(glep2->Eta()) < 2.4) && (gvec->M() > 60 && gvec->M() < 120)){
          pass_selection=kTRUE;
        }
      }
      else if (snamev[isam]=="zee"){
        if ((glep1->Pt() > 25) && (glep2->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 ||(TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5)) && (TMath::Abs(glep2->Eta()) < 1.4442 || (TMath::Abs(glep2->Eta()) > 1.566 && TMath::Abs(glep2->Eta()) < 2.5)) && (gvec->M() > 60 && gvec->M() < 120)){
          pass_selection=kTRUE;
        }
      }
      else if (snamev[isam]=="wpm"){
        if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
          pass_selection=kTRUE;
        }
      }
      else if (snamev[isam]=="wpe"){
        if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
          pass_selection=kTRUE;
        }
      }
      else if (snamev[isam]=="wmm"){
        if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
          pass_selection=kTRUE;
        }
      }
      else if (snamev[isam]=="wme"){
        if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
          pass_selection=kTRUE;
        }
      }
      else{
        cout << "Unsupported Channel" << endl;
      }

      // If the acceptance criteria are met, fills output tree with selected events
      if (pass_selection==kTRUE){
        fiducialWeightGen += weightGen;
        outTree->Fill();
      }
    }
    // Saves output trees to root file in ntuples
    outFile->Write();
    outFile->Close(); 

    // Calculates the gen-level acceptance
    Double_t acceptance = fiducialWeightGen/totalWeightGen;
    Double_t accerr = sqrt(acceptance*(1.-acceptance)/ntotal);

    // Prints gen-level acceptance
    cout << "Fiducial Weight = " << fiducialWeightGen << endl;
    cout << "Inclusive Weight = " << totalWeightGen << endl;
    cout << "Acceptance = " << acceptance << endl;
    cout << "Acceptance Error = " << accerr << endl;

    // Saves gen-level acceptance as a text file
    ofstream txtfile;
    txtfile.open("acceptances.txt", ios::app);
    txtfile << "*" << endl;
    txtfile << "For " << snamev[isam] << endl;
    txtfile << "Fiducial Weight = " << fiducialWeightGen << endl;
    txtfile << "Inclusive Weight = " << totalWeightGen << endl;
    txtfile << "Acceptance = " << acceptance << endl;
    txtfile << "Acceptance Error = " << accerr << endl;
    txtfile.close();
  }

  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("acceptance");

}