#ifndef MYTOOLS_HH
#define MYTOOLS_HH

#include <iostream>
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
namespace toolbox 
{
  Double_t deltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2);
  
  Double_t deltaPhi(const Double_t phi1, const Double_t phi2);
  
  Int_t roundToInt(const Double_t x);
  
  Int_t flavor(TClonesArray *genPartArr, Int_t vid);
  
  void fillGen(TClonesArray *genPartArr, Int_t vid, TLorentzVector* &vec, TLorentzVector* &lep1, TLorentzVector* &lep2, Int_t* lep1q, Int_t* lep2q, Int_t absM);

}

//------------------------------------------------------------------------------------------------------------------------
Double_t toolbox::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);
    
  Double_t deta = eta1-eta2;
  
  return sqrt(dphi*dphi + deta*deta);
}

//------------------------------------------------------------------------------------------------------------------------
Double_t toolbox::deltaPhi(Double_t phi1, Double_t phi2) 
{
  // Compute dPhi between two given angles. Results is in [0,pi].
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);

  return dphi;
}

//------------------------------------------------------------------------------------------------------------------------
Int_t toolbox::roundToInt(Double_t x)
{
  if(x>0)
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)floor(x) : (Int_t)ceil(x);
  else
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)ceil(x) : (Int_t)floor(x);
}

//------------------------------------------------------------------------------------------------------------------------
Int_t toolbox::flavor(TClonesArray *genPartArr, Int_t vid) {

  for (Int_t i=0; i<genPartArr->GetEntries(); i++) {
    const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
    Int_t pdgId=fabs(genloop->pdgId);
    Int_t parentPdgId=fabs(dynamic_cast<baconhep::TGenParticle*>(genPartArr->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId);
    if ( (pdgId==11||pdgId==13||pdgId==15) && (parentPdgId==vid || genloop->status==23) ) return genloop->pdgId;
  }
  return 0;
}

//------------------------------------------------------------------------------------------------------------------------
void toolbox::fillGen(TClonesArray *genPartArr, Int_t vid, TLorentzVector* &vec, TLorentzVector* &lep1, TLorentzVector* &lep2, Int_t* lep1q, Int_t* lep2q, Int_t absM) 
{
  Int_t iv=-1, iv1=-1, iv2=-1;
  TLorentzVector *lepPos=0, *lepNeg=0;
  TLorentzVector *preLepPos=0, *preLepNeg=0;
  for (Int_t i=0; i<genPartArr->GetEntries(); i++) { //loops through all particles in event
    const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]); //gets the ith particle
    if(fabs(genloop->pdgId)==22) continue; //ignores photons
    //cout << i << " " << genloop->pdgId << " " << genloop->status << " " << genloop->parent << " " << genloop->pt << " " << genloop->mass << std::endl;
    if (genloop->status==23 && (fabs(genloop->pdgId)==15 || fabs(genloop->pdgId)==13 || fabs(genloop->pdgId)==11)) { //outgoing lepton
      if (genloop->pdgId<0 && lepPos==0) {//initial fill of gen_lepton
        //cout << "pre lepPos filled" << endl;
      	lepPos=new TLorentzVector(0,0,0,0);
      	lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
      	preLepPos=new TLorentzVector(0,0,0,0);
      	preLepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
      	iv1=i;
        //cout << "lepPos filled" << endl;
      }
      else if (genloop->pdgId>0 && lepNeg==0) {//intial fill of gen_lepton
        //cout << "pre lepNeg filled" << endl;
      	lepNeg=new TLorentzVector(0,0,0,0);
      	lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
      	preLepNeg=new TLorentzVector(0,0,0,0);
      	preLepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
      	iv2=i;
        //cout << "lepNeg filled" << endl;
      }
    }
    else if ((absM==0 && genloop->pdgId==vid && (genloop->status==3||genloop->status==22)) || (absM==1 && fabs(genloop->pdgId)==fabs(vid) && (genloop->status==3||genloop->status==62))) {
      //cout << "pre Boson filled" << endl;
      vec=new TLorentzVector(0,0,0,0); //finds the vector boson
      vec->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
      iv=i;
      //cout << "Boson filled" << endl;
    }
    else if (iv!=-1 && genloop->parent==iv) { //if the boson has propogated into another boson
      if ((absM==0 && genloop->pdgId==vid)||(absM==1 && fabs(genloop->pdgId)==fabs(vid))) { //if it's still a boson re-write the boson
        //cout << "pre Boson refilled" << endl;
        vec->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
        iv=i;
        //cout << "Boson refilled" << endl;
      }
      else if (fabs(genloop->pdgId)==15 || fabs(genloop->pdgId)==13 || fabs(genloop->pdgId)==11) { //if it decayed into leptons save those leptons instead
      	if (genloop->pdgId<0 && lepPos==0) {
          //cout << "pre lepPos refilled" << endl;
          lepPos=new TLorentzVector(0,0,0,0);
          lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
          preLepPos=new TLorentzVector(0,0,0,0);
          preLepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
          iv1=i;
          //cout << "lepPos refilled" << endl;
      	}
      	else if (genloop->pdgId>0 && lepNeg==0) {
          //cout << "pre lepNeg refilled" << endl;
          lepNeg=new TLorentzVector(0,0,0,0);
          lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
      	  preLepNeg=new TLorentzVector(0,0,0,0);
      	  preLepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
          iv2=i;
          //cout << "lepNeg refilled" << endl;
      	}
      }
    }
    else if (iv1!=-1 && genloop->parent==iv1) {// if the old lepton propagated and is still that lepton re-write the lepton
      //cout << "pre lepPos refilled 2" << endl;
      lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass); 
      iv1=i;
      //cout << "lepPos refilled 2" << endl;
    }
    else if (iv2!=-1 && genloop->parent==iv2) {
      //cout << "pre lepNeg refilled 2" << endl;
      lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
      iv2=i;
      //cout << "lepNeg refilled 2" << endl;
    }
  }//end particle loop

  if (iv==-1 && preLepNeg && preLepPos) { //if the boson wasn't found make a pseudo_boson out of the two immediate leptons
    TLorentzVector temp = *preLepNeg + *preLepPos;
    vec->SetPtEtaPhiM(temp.Pt(), temp.Eta(), temp.Phi(), temp.M());
    //vec->SetPtEtaPhiM(0,0,0,temp.M());
    //cout << "No Boson " << preLepNeg->Pt() << " " << preLepPos->Pt() << " " << temp.M() << " " << temp.Pt() << " " << temp.Eta() << " " << temp.Phi() << " " <<  vec->M() << std::endl;
  }

  if (lepPos && lepNeg) { //if two leptons are selected for with opposite charges
    if (lepPos->Pt()>lepNeg->Pt()) { //puts the highest pt one first
      //cout << "pre leptons ordered +-" << endl;
      lep1->SetPtEtaPhiM(lepPos->Pt(), lepPos->Eta(), lepPos->Phi(), lepPos->M());
      lep2->SetPtEtaPhiM(lepNeg->Pt(), lepNeg->Eta(), lepNeg->Phi(), lepNeg->M());
      *lep1q=1;
      *lep2q=-1;
      //cout << "leptons ordered +-" << endl;
    }
    else {
      //cout << "pre leptons ordered -+" << endl;
      lep1->SetPtEtaPhiM(lepNeg->Pt(), lepNeg->Eta(), lepNeg->Phi(), lepNeg->M());
      lep2->SetPtEtaPhiM(lepPos->Pt(), lepPos->Eta(), lepPos->Phi(), lepPos->M());
      *lep1q=-1;
      *lep2q=1;
      //cout << "leptons ordered -+" << endl;
    }
  }
  else if (lepPos) //in case only one lepton is selected
    {
      //cout << "pre only positive lepton" << endl;
      lep1->SetPtEtaPhiM(lepPos->Pt(), lepPos->Eta(), lepPos->Phi(), lepPos->M());
      *lep1q=1;
      //cout << "only positive lepton" << endl;
    }
  else if (lepNeg) //in case only one lepton is selected
    {
      //cout << "pre only negative lepton" << endl;
      lep1->SetPtEtaPhiM(lepNeg->Pt(), lepNeg->Eta(), lepNeg->Phi(), lepNeg->M());
      *lep1q=-1;
      //cout << "only negative lepton" << endl;
    }
  //std:://cout << "Vector boson " << vec->Pt() << " " << vec->Eta() << std::endl;

  //clears intermediate variables, leaves only lep1, lep2, and Vector
  delete preLepNeg;
  delete preLepPos;
  delete lepNeg;
  delete lepPos;
  
  preLepNeg=0; preLepPos=0;
  lepNeg=0; lepPos=0;
  //cout << "Done gen fill" << endl;

  }
#endif
