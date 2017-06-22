#if !defined(__CINT__) //|| defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <TRandom3.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <sstream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <TChain.h>
#include <TH1.h>
#include "TLorentzVector.h"

#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
Double_t totalWeightGen=0;
Double_t MASS_Z = 91.2;
using namespace std;

#endif



void makeFlat(TString input="/data/t3home000/sabrandt/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
        TString output="all_zmumu.root",
        Int_t vid=13) {
  
  TChain chain("Events");
  chain.Add(input);
  
  // Data structures to store info from TTrees
  baconhep::TGenEventInfo *info = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr      = new TClonesArray("baconhep::TGenParticle");
  
  chain.SetBranchAddress("GenEvtInfo",  &info);        TBranch *infoBr     = chain.GetBranch("GenEvtInfo");
  chain.SetBranchAddress("GenParticle", &genPartArr);  TBranch *partBr     = chain.GetBranch("GenParticle");

  Double_t dimuon_mass, weightGen;

  TLorentzVector *muon1=0, *muon2=0;
  TLorentzVector dimuon;

  std::vector<float> *lheweight = new std::vector<float>();

  TFile *ofile = new TFile(output, "recreate");
  TTree *otree = new TTree("Events", "Events");

  otree->Branch("dimuon_mass",     &dimuon_mass,     "dimuon_mass/D");      // dilepton mass
  otree->Branch("muon1",       "TLorentzVector",  &muon1);     // lepton2 4-vector
  otree->Branch("muon2",       "TLorentzVector",  &muon2);     // lepton2 4-vector
//  otree->Branch("lheweight",  &lheweight);
  otree->Branch("weightGen",   &weightGen,   "weightGen/D");    // event weight (MC)

  Int_t parentnum1, parentnum2;

  //for (Int_t ie=0; ie<chain.GetEntries(); ie++) {
  for (Int_t ie=0; ie<1000000; ie++) {
     if(ie%1000==0) cout << "Processing event " << ie << ". " << (double)ie/(double)chain.GetEntries()*100 << " percent done with this file." << endl;

    infoBr->GetEntry(ie);
    genPartArr->Clear(); partBr->GetEntry(ie);
    if (genPartArr->GetEntries()==0) continue;
    muon1->SetPtEtaPhiM(0,0,0,0);
    muon2->SetPtEtaPhiM(0,0,0,0);
    dimuon.SetPtEtaPhiM(0,0,0,0);
    parentnum1=0;
    parentnum2=0;
    TLorentzVector part1, part2, disum;
    Int_t parent1, parent2;
    Double_t sum_mass, mass_diff;
    Double_t min_diff = MASS_Z;

    Bool_t isz_mumu=kFALSE, has_Z=kFALSE, mufromz=kFALSE, mu2fromz=kFALSE;
    for (Int_t i=0; i<genPartArr->GetEntries(); i++) {
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
      if (genloop->pdgId!=23 || genloop->status!=22) continue;
//      cout << "Event number: "<< ie << ", " << "Particle number: " << i << ", " << "PID: " << genloop->pdgId << ", " << "Parent: " << genloop->parent << ", " << "Status: " << genloop->status << ", " << "PT: " << genloop->pt << ", " << "Eta: " << genloop->eta << endl;
      has_Z = kTRUE;
    }
    for (Int_t i=0; i<genPartArr->GetEntries(); i++) {
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
      if (genloop->parent < 0) continue;
      const baconhep::TGenParticle* genloop2 = (baconhep::TGenParticle*) ((*genPartArr)[genloop->parent]);
      if (genloop->pdgId==vid && genloop->status==23 && genloop2->pdgId==23){
        mufromz = kTRUE;
      }
    }

    for (Int_t i=0; i<genPartArr->GetEntries(); i++) {
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
      if (genloop->parent < 0) continue;
      const baconhep::TGenParticle* genloop2 = (baconhep::TGenParticle*) ((*genPartArr)[genloop->parent]);
      if (genloop->pdgId==-vid && genloop->status==23 && genloop2->pdgId==23){
        mu2fromz = kTRUE;
      }
    }
    if (has_Z==kTRUE && mufromz==kTRUE && mu2fromz==kTRUE){
        isz_mumu = kTRUE;
    }
    if(isz_mumu == kFALSE) continue;

    for (Int_t i=0; i<genPartArr->GetEntries(); i++) {
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
//      cout << "Event number: "<< ie << ", " << "Particle number: " << i << ", " << "PID: " << genloop->pdgId << ", " << "Parent: " << genloop->parent << ", " << "Status: " << genloop->status << ", " << "PT: " << genloop->pt << ", " << "Eta: " << genloop->eta << endl;
      if (genloop->pdgId!=vid || genloop->status!=1) continue;
        part1.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
//        cout << "Event number: "<< ie << ", " << "Particle number: " << i << ", " << "PID: " << genloop->pdgId << ", " << "Parent: " << genloop->parent << ", " << "Status: " << genloop->status << ", " << "PT: " << genloop->pt << ", " << "Eta: " << genloop->eta << endl;
        parent1=genloop->parent;
      for (Int_t j=0; j<genPartArr->GetEntries(); j++) {
        const baconhep::TGenParticle* genloop2 = (baconhep::TGenParticle*) ((*genPartArr)[j]);
        if (genloop2->pdgId!=-vid || genloop2->status!=1) continue;
        
        part2.SetPtEtaPhiM(genloop2->pt, genloop2->eta, genloop2->phi, genloop2->mass);
//      cout << "Event number: "<< ie << ", " << "Particle number: " << j << ", " << "PID: " << genloop2->pdgId << ", " << "Parent: " << genloop2->parent << ", " << "Status: " << genloop2->status << ", " << "PT: " << genloop2->pt << ", " << "Eta: " << genloop2->eta << endl;
        parent2=genloop2->parent;
        
        disum = part1 + part2;
        sum_mass = disum.M();
        mass_diff = TMath::Abs(MASS_Z - sum_mass); // Selects the muon pair with mass closest to Z mass
            if(mass_diff < min_diff){
                *muon1=part1;
                *muon2=part2;
                parentnum1=parent1;
                parentnum2=parent2;
                min_diff = mass_diff;
            }
      }
    }

    if(parentnum1==0 || parentnum2==0) continue;

    dimuon = *muon1 + *muon2;
    dimuon_mass = dimuon.M();
//    cout << dimuon_mass << endl;

    weightGen = info->weight;
//    cout << weightGen << endl;

    // lheweight->clear();
    // for (int j = 0; j<109; j++)
    //  {
    //    lheweight->push_back(info->lheweight[j]);
    //  }

    otree->Fill();
    totalWeightGen+=info->weight;
  }

  cout << "Total Gen Weight = " << totalWeightGen << endl;
  ofile->Write();
  ofile->Close();
}