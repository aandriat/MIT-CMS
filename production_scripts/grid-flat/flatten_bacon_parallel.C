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

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "../Utils/MyTools.hh"      // various helper functions
#endif

//=== MAIN MACRO ================================================================================================= 

void flatten_bacon_parallel(const TString infilename="bxx.root", // input file followed by _bacon.root
                  const TString process="bxx", // process type
                  const TString outputDir=".",   // ntuple directory
                  const Int_t n_events=0 //Number of events, 0 for all
            ) {
  gBenchmark->Start("flatten_bacon_parallel");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // // Create output directory
  // gSystem->mkdir(outputDir,kTRUE);

  // Declare output ntuple variables
  Double_t weightGen;
  TLorentzVector *glep1=0, *glep2=0, *gvec=0;
  Int_t glepq1, glepq2, fiducial;
  std::vector<float> *lheweight = new std::vector<float>();

  // Data structures to store info from Bacon TTrees
  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr      = new TClonesArray("baconhep::TGenParticle");
  
  // Data structure to hold Bacon TFile
  TFile *infile=0;
  TTree *eventTree=0;

    // Set up output ntuple
    TString outfilename = "flat_output.root";
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");
    outTree->Branch("glepq1",     &glepq1,     "glepq1/I");      // lepton1 charge
    outTree->Branch("glepq2",     &glepq2,     "glepq2/I");      // lepton2 charge
    outTree->Branch("glep1",       "TLorentzVector",  &glep1);     // lepton1 4-vector
    outTree->Branch("glep2",       "TLorentzVector",  &glep2);     // lepton2 4-vector
    outTree->Branch("gvec",       "TLorentzVector",  &gvec);     // boson 4-vector
    outTree->Branch("weightGen",   &weightGen,   "weightGen/D");    // event weight (MC)
    outTree->Branch("lheweight",  &lheweight); // lhe weights vector (same for all events)
    outTree->Branch("fiducial", &fiducial, "fiducial/I"); // If event is in fiducial region of detector, store as 1

      // Read input file and get the TTrees
      cout << "Processing " << infilename << endl;
      infile = TFile::Open(infilename); 
      assert(infile);
      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);  

      eventTree->SetBranchAddress("GenEvtInfo", &gen); TBranch *genBr = eventTree->GetBranch("GenEvtInfo");
      eventTree->SetBranchAddress("GenParticle",&genPartArr); TBranch *genPartBr = eventTree->GetBranch("GenParticle");

      // Define number of events to loop over (n_events=0 = all events)
      Int_t max_events = n_events;
      if (n_events==0){
        max_events = eventTree->GetEntries();
      }

      Double_t totalWeightGen=0;
      for(Int_t ientry=0; ientry<max_events; ientry++) {
        genBr->GetEntry(ientry);
          totalWeightGen+=gen->weight;
      }

      //
      // loop over events
      //
      for(Int_t ientry=0; ientry<max_events; ientry++) {
        genBr->GetEntry(ientry);
        genPartArr->Clear(); genPartBr->GetEntry(ientry);
        if(ientry%100000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;

        if (genPartArr->GetEntries()==0){
          cout << "gen_particles not found" << endl;
        }

        glep1->SetPtEtaPhiM(0,0,0,0);
        glep2->SetPtEtaPhiM(0,0,0,0);
        gvec->SetPtEtaPhiM(0,0,0,0);
        glepq1=-99;
        glepq2=-99;
        fiducial=0;
        Int_t BOSON_ID=0;
        Int_t LEPTON_ID=0;
        Bool_t isSignal=kFALSE;
        Bool_t isWrongFlavor=kFALSE;
        Bool_t mass_Zcut=kTRUE;

          if (process=="zmm"){
            LEPTON_ID = 13;
            BOSON_ID = 23;                        
            isSignal = process.Contains("zmm",TString::kIgnoreCase);// Assume signal sample is given name "zmm" - flag to store GEN Z kinematics
            isWrongFlavor = (process.CompareTo("zxx",TString::kIgnoreCase)==0);// flag to reject Z->mm events when selecting at wrong-flavor background events
            if (isWrongFlavor && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;// veto wrong channel
            else if (isSignal && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
            toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2, &glepq1, &glepq2,0);
            if (gvec->M() < 60 || gvec->M() > 120){
              mass_Zcut=kFALSE;
            }
            // cout << "Event checked" << endl;
            if ((glep1->Pt() > 25) && (glep2->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4) && (TMath::Abs(glep2->Eta()) < 2.4) && (gvec->M() > 60 && gvec->M() < 120)){
              fiducial=1;
              // cout << "fiducial" << endl;
            }


          }
          else if (process=="zee"){
            LEPTON_ID = 11;
            BOSON_ID = 23;                        
            isSignal = process.Contains("zee",TString::kIgnoreCase);// Assume signal sample is given name "zmm" - flag to store GEN Z kinematics
            isWrongFlavor = (process.CompareTo("zxx",TString::kIgnoreCase)==0);// flag to reject Z->mm events when selecting at wrong-flavor background events

            // cout << "Event in zee" << endl;

            if (isWrongFlavor && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;// veto wrong channel
            else if (isSignal && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
            toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2, &glepq1, &glepq2,0);
            if (gvec->M() < 60 || gvec->M() > 120){
              mass_Zcut=kFALSE;
            }
            // cout << "Event checked" << endl;
            if ((glep1->Pt() > 25) && (glep2->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 ||(TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5)) && (TMath::Abs(glep2->Eta()) < 1.4442 || (TMath::Abs(glep2->Eta()) > 1.566 && TMath::Abs(glep2->Eta()) < 2.5)) && (gvec->M() > 60 && gvec->M() < 120)){
              fiducial=1;
              // cout << "fiducial" << endl;
            }
          } 
          else if (process=="wpm"){
            LEPTON_ID = -13;
            BOSON_ID = 24;                        
            isSignal = process.Contains("wpm",TString::kIgnoreCase);// Assume signal sample is given name "zmm" - flag to store GEN Z kinematics
            isWrongFlavor = (process.CompareTo("wpx",TString::kIgnoreCase)==0);// flag to reject Z->mm events when selecting at wrong-flavor background events
            if (isWrongFlavor && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;// veto wrong channel
            else if (isSignal && toolbox::flavor(genPartArr, BOSON_ID)!=LEPTON_ID) continue;
            toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2, &glepq1, &glepq2,0);
            // cout << "Event checked" << endl;
            if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
              fiducial=1;
              // cout << "fiducial" << endl;
            }
          }
          else if (process=="wpe"){
            LEPTON_ID = -11;
            BOSON_ID = 24;                        
            isSignal = process.Contains("wpe",TString::kIgnoreCase);// Assume signal sample is given name "zmm" - flag to store GEN Z kinematics
            isWrongFlavor = (process.CompareTo("wpx",TString::kIgnoreCase)==0);// flag to reject Z->mm events when selecting at wrong-flavor background events
            if (isWrongFlavor && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;// veto wrong channel
            else if (isSignal && toolbox::flavor(genPartArr, BOSON_ID)!=LEPTON_ID) continue;
            toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2, &glepq1, &glepq2,0);
            // cout << "Event checked" << endl;
            if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
              fiducial=1;
              // cout << "fiducial" << endl;
            }
          }
          else if (process=="wmm"){
            LEPTON_ID = 13;
            BOSON_ID = -24;                        
            isSignal = process.Contains("wmm",TString::kIgnoreCase);// Assume signal sample is given name "zmm" - flag to store GEN Z kinematics
            isWrongFlavor = (process.CompareTo("wmx",TString::kIgnoreCase)==0);// flag to reject Z->mm events when selecting at wrong-flavor background events
            if (isWrongFlavor && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;// veto wrong channel
            else if (isSignal && toolbox::flavor(genPartArr, BOSON_ID)!=LEPTON_ID) continue;
            toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2, &glepq1, &glepq2,0);
            // cout << "Event checked" << endl;
            if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
              fiducial=1;
              // cout << "fiducial" << endl;
            }
          }
          else if (process=="wme"){
            LEPTON_ID = 11;
            BOSON_ID = -24;                        
            isSignal = process.Contains("wme",TString::kIgnoreCase);// Assume signal sample is given name "zmm" - flag to store GEN Z kinematics
            isWrongFlavor = (process.CompareTo("wmx",TString::kIgnoreCase)==0);// flag to reject Z->mm events when selecting at wrong-flavor background events
            if (isWrongFlavor && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;// veto wrong channel
            else if (isSignal && toolbox::flavor(genPartArr, BOSON_ID)!=LEPTON_ID) continue;
            toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2, &glepq1, &glepq2,0);
            // cout << "Event checked" << endl;
            if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
              fiducial=1;
              // cout << "fiducial" << endl;
            }
          }
          else if (process=="wm"){
            LEPTON_ID = 13;
            BOSON_ID = 24;                        
            isSignal = process.Contains("wm",TString::kIgnoreCase);// Assume signal sample is given name "zmm" - flag to store GEN Z kinematics
            isWrongFlavor = (process.CompareTo("wx",TString::kIgnoreCase)==0);// flag to reject Z->mm events when selecting at wrong-flavor background events
            if (isWrongFlavor && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;// veto wrong channel
            else if (isSignal && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
            toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2, &glepq1, &glepq2,1);
            // cout << "Selected event!" << endl;
            // if (glepq2!=-99){
            //   cout << "Selected neutrino!" << endl;
            //   cout << "Neutrino pt is " << glep2->Pt() << endl;
            // }
            // cout << "Event checked" << endl;
            if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 2.4)){
              fiducial=1;
              // cout << "fiducial" << endl;
            }

          }
          else if (process=="we"){
            LEPTON_ID = 11;
            BOSON_ID = 24;                        
            isSignal = process.Contains("we",TString::kIgnoreCase);// Assume signal sample is given name "zmm" - flag to store GEN Z kinematics
            isWrongFlavor = (process.CompareTo("wx",TString::kIgnoreCase)==0);// flag to reject Z->mm events when selecting at wrong-flavor background events
            if (isWrongFlavor && fabs(toolbox::flavor(genPartArr, BOSON_ID))==LEPTON_ID) continue;// veto wrong channel
            else if (isSignal && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
            toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2, &glepq1, &glepq2,1);
            // cout << "Event checked" << endl;
            // cout << "Selected event!" << endl;
            // if (glepq2!=-99){
            //   cout << "Selected neutrino!" << endl;
            //   cout << "Neutrino pt is " << glep2->Pt() << endl;
            // }
            if ((glep1->Pt() > 25) && (TMath::Abs(glep1->Eta()) < 1.4442 || (TMath::Abs(glep1->Eta()) > 1.566 && TMath::Abs(glep1->Eta()) < 2.5))){
              fiducial=1;
              // cout << "fiducial" << endl;
            }            
          }
          else{
            cout << "Unsupported Channel" << endl;
          }
        
        if (mass_Zcut==kFALSE) continue; // definition of Z boson
        weightGen = gen->weight;
        //cout << "weightGen " << weightGen << endl;

        // Stores lheweights in vector
        // cout << "LHE weights " << endl;
        if (ientry==0){
          cout << "LHE Weight size " << gen->lheweight.size() << endl;
        }

        lheweight->clear();
        for (Int_t j = 0; j< (Int_t) gen->lheweight.size(); j++)
         {
           lheweight->push_back(gen->lheweight[j]);
           // cout << gen->lheweight[j] << endl; 
         }

        // Saves events from correct signal to flattened ntuple
        outTree->Fill();
      }
      cout << "Total GenWeight for sample " << infilename << " is " << totalWeightGen << endl;
      outFile->Write();
      outFile->Close();
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //=============================================================================================================

  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("flatten_bacon_parallel"); 
}
