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
                  const TString outputDir="ntuples",   // ntuple directory
                  const Int_t n_events=0 //Number of events, 0 for all
            ){
  
  gBenchmark->Start("flatten_gen");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

// (TString input="/home/aandriat/5TeV_Samples/DYJetsToLL_TuneCUETP8M1_5020GeV-amcatnloFXFX-pythia8.root",
//         TString output="ntuples/all_Zmm_gen.root") {


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  //
  // parse .conf file
  //
  confParse(conf, snamev, samplev);
  //
  // Declare ntuple variables
  //
  Double_t weightGen;
  TLorentzVector *glep1=0, *glep2=0, *gvec=0;
  //TLorentzVector vecsum;
  Int_t glepq1, glepq2;
  std::vector<float> *lheweight = new std::vector<float>();

  TFile *infile=0;
  TTree *intree=0;

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
    CSample* samp = samplev[isam];
    // Set up output ntuple
    TString infilename = outputDir + TString("/") + snamev[isam] + TString("_gen.root");

    TFile *infile = TFile::Open(infilename);         assert(infile);
    TTree *intree = (TTree*)infile->Get("Events"); assert(intree);

    intree->SetBranchAddress("glepq1",   &glepq1);     // lepton1 charge 
    intree->SetBranchAddress("glepq1",   &glepq1);     // lepton2 charge
    intree->SetBranchAddress("glep1",   &glep1);     // lepton1 4-vector
    intree->SetBranchAddress("glep2",   &glep2);     // lepton2 4-vector
    intree->SetBranchAddress("gvec",   &gvec);     // boson 4-vector
    intree->SetBranchAddress("weightGen",   &weightGen);     // event weights
    //  intree->SetBranchAddress("lheweight",   &lheweight);     // lheweights
    
    cout << "Processing " << snamev[isam] << endl;

      Double_t totalWeightGen=0;
      Double_t fiducialWeightGen=0;
      Double_t ntotal=0;
      for(Int_t ientry=0; ientry<intree->GetEntries(); ientry++) {
        intree->GetEntry(ientry);
          totalWeightGen+=weightGen;
          ntotal += weightGen/abs(weightGen);
      }

      //
      // loop over events
      //
      Int_t max_events = n_events;
      if (n_events==0){
        max_events = intree->GetEntries();
      }
        if (snamev[isam]=="zmm"){
            for(Int_t ientry=0; ientry<max_events; ientry++) {
              intree->GetEntry(ientry);   
              if(ientry%10000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;

              if (glep1->Pt() < 25) continue;
              if (glep2->Pt() < 25) continue;
              if (TMath::Abs(glep1->Eta()) > 2.4) continue;
              if (TMath::Abs(glep2->Eta()) > 2.4) continue;
              //vecsum = *glep1 + *glep2;
              if (gvec->M() < 60 || gvec->M() > 120) continue;
              //if (vecsum.M() < 60 || vecsum.M() > 120) continue;
              fiducialWeightGen += weightGen;
            }
        }
        else if (snamev[isam]=="zee"){
            for(Int_t ientry=0; ientry<max_events; ientry++) {
              intree->GetEntry(ientry);   
              if(ientry%10000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;

              if (glep1->Pt() < 25) continue;
              if (glep2->Pt() < 25) continue;
              if ((TMath::Abs(glep1->Eta()) > 1.4442 and TMath::Abs(glep1->Eta()) < 1.566) or TMath::Abs(glep1->Eta()) > 2.5) continue;
              if ((TMath::Abs(glep2->Eta()) > 1.4442 and TMath::Abs(glep2->Eta()) < 1.566) or TMath::Abs(glep2->Eta()) > 2.5) continue;

              //vecsum = *glep1 + *glep2;
              if (gvec->M() < 60 || gvec->M() > 120) continue;
              //if (vecsum.M() < 60 || vecsum.M() > 120) continue;

              fiducialWeightGen += weightGen;
            }
        }
        else if (snamev[isam]=="wpm"){
            for(Int_t ientry=0; ientry<max_events; ientry++) {
              intree->GetEntry(ientry);   
              if(ientry%10000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;

              if (glepq1!=1) continue;
              if (glep1->Pt() < 25) continue;
              if (TMath::Abs(glep1->Eta()) > 2.4) continue;

              fiducialWeightGen += weightGen;
            }
        }
        else if (snamev[isam]=="wpe"){
            for(Int_t ientry=0; ientry<max_events; ientry++) {
              intree->GetEntry(ientry);   
              if(ientry%10000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;

              if (glepq1!=1) continue;
              if (glep1->Pt() < 25) continue;
              if ((TMath::Abs(glep1->Eta()) > 1.4442 and TMath::Abs(glep1->Eta()) < 1.566) or TMath::Abs(glep1->Eta()) > 2.5) continue;

              fiducialWeightGen += weightGen;
            }
        }
        else if (snamev[isam]=="wmm"){
            for(Int_t ientry=0; ientry<max_events; ientry++) {
              intree->GetEntry(ientry);   
              if(ientry%10000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;

              if (glepq1!=-1) continue;
              if (glep1->Pt() < 25) continue;
              if (TMath::Abs(glep1->Eta()) > 2.4) continue;

              fiducialWeightGen += weightGen;
            }
        }
        else if (snamev[isam]=="wme"){
            for(Int_t ientry=0; ientry<max_events; ientry++) {
              intree->GetEntry(ientry);   
              if(ientry%10000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;

              if (glepq1!=-1) continue;
              if (glep1->Pt() < 25) continue;
              if ((TMath::Abs(glep1->Eta()) > 1.4442 and TMath::Abs(glep1->Eta()) < 1.566) or TMath::Abs(glep1->Eta()) > 2.5) continue;

              fiducialWeightGen += weightGen;
            }
        }
        else{
          cout << "Unsupported Channel" << endl;
        }

      Double_t acceptance = fiducialWeightGen/totalWeightGen;
      Double_t accerr = sqrt(acceptance*(1.+acceptance)/ntotal);

      cout << "Fiducial Weight = " << fiducialWeightGen << endl;
      cout << "Inclusive Weight = " << totalWeightGen << endl;
      cout << "Acceptance = " << acceptance << endl;
      cout << "Acceptance Error = " << accerr << endl;


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
}