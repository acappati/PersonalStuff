// **************************************************
// 
// example of how to pair taus and compute inv mass
//
// usage:    specify input file
//           run with: 
//                    root 
//                    .L yields4lX.C++
//                    yields4lX()
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
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"



void yields4lX(){

  Float_t lumi = 140.;

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Float_t overallEventWeight;
  Float_t PUWeight;
  Float_t genHEPMCweight;
  Float_t dataMCWeight;
  Float_t trigEffWeight;

  Short_t ZZsel;
  Float_t genBR;
  Float_t xsec;
  vector<Float_t> *LepEta = 0;
  Short_t nTaus = 0;
  vector<Float_t> *tauPt     = 0;
  vector<Float_t> *tauEta    = 0;
  vector<Float_t> *tauPhi    = 0;
  vector<Float_t> *tauMass   = 0;
  vector<Float_t> *photonEta = 0;
  vector<Float_t> *JetEta    = 0;
  vector<float>* JetIsBtagged;




  //TFile* inputFile = new TFile("mc_4ltautau/ZZXAnalysis.root");
  //TFile* inputFile = new TFile("mc_4lgammagamma/ZZXAnalysis.root");
  TFile* inputFile = new TFile("mc_4lbb/ZZXAnalysis.root");
  
  TH1F* hCounters = (TH1F*)inputFile->Get("ZZTree/Counters");
  Float_t gen_sumWeights = (Float_t)hCounters->GetBinContent(40);
  Float_t partialSampleWeight = lumi * 1000. / gen_sumWeights ;
  Float_t gen_sumPUWeight = (Float_t)hCounters->GetBinContent(42);
  Float_t gen_sumGenMCweight = (Float_t)hCounters->GetBinContent(41);


  TTree* inputTree = (TTree*)inputFile->Get("ZZTree/candTree");
  inputTree->SetBranchAddress("RunNumber",   &nRun);
  inputTree->SetBranchAddress("EventNumber", &nEvent);
  inputTree->SetBranchAddress("LumiNumber",  &nLumi);
  inputTree->SetBranchAddress("ZZsel", &ZZsel);
  inputTree->SetBranchAddress("overallEventWeight", &overallEventWeight);
  inputTree->SetBranchAddress("PUWeight", &PUWeight);
  inputTree->SetBranchAddress("genHEPMCweight", &genHEPMCweight);
  inputTree->SetBranchAddress("dataMCWeight", &dataMCWeight);
  inputTree->SetBranchAddress("trigEffWeight", &trigEffWeight);
  inputTree->SetBranchAddress("xsec", &xsec);  
  inputTree->SetBranchAddress("genBR", &genBR);  
  inputTree->SetBranchAddress("LepEta", &LepEta);
  inputTree->SetBranchAddress("tauEta", &tauEta);
  inputTree->SetBranchAddress("photonEta", &photonEta);
  inputTree->SetBranchAddress("JetEta", &JetEta);
  inputTree->SetBranchAddress("JetIsBtagged", &JetIsBtagged);

  
  Float_t yield=0.;

  
  // --- loop over tree entries
  Long64_t entries = inputTree->GetEntries();
  for (Long64_t z=0; z<entries; ++z){

    inputTree->GetEntry(z);

    //if(z>100) break;  // debug


    if(LepEta->size()!=4){
        cout<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", stored "<<LepEta->size()<<" leptons instead of 4"<<endl;
        continue;
    }


    if( !(ZZsel>=90) ) continue;


    // // tau sel
    // if(tauEta->size()<2) continue;

    // // photon sel
    // if(photonEta->size()<2) continue;

    // jet sel
    if(JetEta->size()<2) continue;
    if(!JetIsBtagged->at(0) && !JetIsBtagged->at(1)) continue; // at least 1bjet
    cout<<JetIsBtagged->at(0)<<" "<<JetIsBtagged->at(1)<<endl;





    Float_t eventWeight = partialSampleWeight * xsec * genBR * overallEventWeight;
 
    yield += eventWeight;

    
 
    cout<<"-------------------------------------"<<endl;
    cout<<"eventWeight:               "<<eventWeight<<endl;
    cout<<"partialSampleWeight:       "<<partialSampleWeight<<endl;
    cout<<"gen_sumWeights(bin40):     "<<gen_sumWeights<<endl;
    cout<<"gen_sumPUWeight(bin42):    "<<gen_sumPUWeight<<endl;
    cout<<"gen_sumGenMCweight(bin41): "<<gen_sumGenMCweight<<endl;
    cout<<"overallEventWeight:        "<<overallEventWeight<<endl;
    cout<<"prod:                      "<<PUWeight*genHEPMCweight*dataMCWeight*trigEffWeight<<endl;
    cout<<"        "<<endl;
    cout<<"PUWeight:                  "<<PUWeight<<endl;
    cout<<"genHEPMCweight:            "<<genHEPMCweight<<endl;
    cout<<"dataMCWeight:              "<<dataMCWeight<<endl;
    cout<<"trigEffWeight:             "<<trigEffWeight<<endl;
    cout<<"       "<<endl;
    cout<<" BR: "<<genBR<<endl;
    

  }

  
  cout<<" ---"<<endl;
  cout<<" yield: "<<yield<<endl;
  cout<<" ---"<<endl;

  }


  


