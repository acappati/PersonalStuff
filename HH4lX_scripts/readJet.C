// **************************************************
// 
// example of how to pair jets and compute inv mass
//
// usage:    specify input file
//           run with: 
//                    root -l -b -q readJet.C++
//
// **************************************************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"



void readJet(){

  Short_t nCleanedJets            = 0;
  Short_t nCleanedJetsPt30        = 0;
  Short_t nCleanedJetsPt30BTagged = 0;
  vector<Float_t> *JetPt     = 0;
  vector<Float_t> *JetEta    = 0;
  vector<Float_t> *JetPhi    = 0;
  vector<Float_t> *JetMass   = 0;


  TFile* inputFile = new TFile("ZZXAnalysis.root");
  TTree* inputTree = (TTree*)inputFile->Get("ZZTree/candTree");


  inputTree->SetBranchAddress("nCleanedJets",            &nCleanedJets           );
  inputTree->SetBranchAddress("nCleanedJetsPt30",        &nCleanedJetsPt30       );
  inputTree->SetBranchAddress("nCleanedJetsPt30BTagged", &nCleanedJetsPt30BTagged);
  inputTree->SetBranchAddress("JetPt",   &JetPt);
  inputTree->SetBranchAddress("JetEta",  &JetEta);
  inputTree->SetBranchAddress("JetPhi",  &JetPhi);
  inputTree->SetBranchAddress("JetMass", &JetMass);


  
  // --- loop over tree entries
  Long64_t entries = inputTree->GetEntries();
  for (Long64_t z=0; z<entries; ++z){

    inputTree->GetEntry(z);


    // define a TLorenzVector for each Jet in the event
    // --> per each event we have a vector of TLorenzVectors
    vector<TLorentzVector> JetLorentzVec;
    
    for(UInt_t i=0; i<JetPt->size(); i++){
    
      TLorentzVector temp;
      temp.SetPtEtaPhiM(JetPt->at(i), JetEta->at(i), JetPhi->at(i), JetMass->at(i));
      JetLorentzVec.push_back(temp);  
      
    }
    

    // example of how to get Pt, Eta, Phi, Mass of this TLorentzVector
    for(UInt_t i=0; i<JetPt->size(); i++){
      cout<<JetLorentzVec.at(i).Pt() <<endl;
      cout<<JetLorentzVec.at(i).Eta()<<endl;
      cout<<JetLorentzVec.at(i).Phi()<<endl;
      cout<<JetLorentzVec.at(i).M()<<endl;
    }


    // example of how to pair two TLorentzVectors
    if(JetLorentzVec.size() > 1){
      cout<<JetLorentzVec.at(0).M()<<" "<<JetLorentzVec.at(1).M()<<endl;
      float jetPairMass = (JetLorentzVec.at(0)+JetLorentzVec.at(1)).M(); // invariant mass of the jet pair
      cout<<jetPairMass<<endl;     
    }

  }


}
