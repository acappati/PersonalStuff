// *****************************************
//
// this macro reads a root file of data and 
// selects Z->4mu events
//
// to run:
//     root -l -b -q takeEvent.C++ 
//
// *****************************************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

#include "TString.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TAxis.h"



using namespace std;


void takeEvent(){


  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;
  Short_t ZZsel;
  Float_t ZZMass;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  vector<Float_t> *LepLepId = 0;

  TFile* inputFile = new TFile("/data3/Higgs/170907_Data2016/AllData/ZZ4lAnalysis2016.root");
  TTree* inputTree = (TTree*)inputFile->Get("ZZTree/candTree");

  inputTree->SetBranchAddress("RunNumber", &nRun);
  inputTree->SetBranchAddress("EventNumber", &nEvent);
  inputTree->SetBranchAddress("LumiNumber", &nLumi);
  inputTree->SetBranchAddress("ZZsel", &ZZsel);     
  inputTree->SetBranchAddress("ZZMass", &ZZMass);
  inputTree->SetBranchAddress("LepPt", &LepPt);
  inputTree->SetBranchAddress("LepEta", &LepEta);
  inputTree->SetBranchAddress("LepPhi", &LepPhi);
  inputTree->SetBranchAddress("LepLepId", &LepLepId);


  // --- loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; z++){

    inputTree->GetEntry(z);


    if(ZZsel < 0.) continue; // skip events that do not pass the trigger

    if(LepLepId->at(0) == -13 && LepLepId->at(1) == 13 && LepLepId->at(2) == -13 && LepLepId->at(3) == 13){   //ZZ->4mu events


      if(fabs(LepEta->at(0)) >= 0.8 || fabs(LepEta->at(1)) >= 0.8 || fabs(LepEta->at(2)) >= 0.8 || fabs(LepEta->at(3)) >= 0.8) continue;


      // --- Z->4mu events
      if(ZZMass >= 90.3 && ZZMass < 92.){

        cout<<"----------------------------"<<endl;
        cout<<"Z->4mu"<<endl;
        cout<<"RunNumber:LumiNumber:EventNumber = "<<nRun<<":"<<nLumi<<":"<<nEvent<<endl;
        cout<<"Mass   = "<<ZZMass<<endl;
        cout<<"LepID  = "<<LepLepId->at(0)<<setw(12)<<LepLepId->at(1)<<setw(12)<<LepLepId->at(2)<<setw(12)<<LepLepId->at(3)<<endl;
        cout<<"LepPt  = "<<LepPt->at(0)<<setw(12)<<LepPt->at(1)<<setw(12)<<LepPt->at(2)<<setw(12)<<LepPt->at(3)<<endl;
        cout<<"LepEta = "<<LepEta->at(0)<<setw(12)<<LepEta->at(1)<<setw(12)<<LepEta->at(2)<<setw(12)<<LepEta->at(3)<<endl;
        cout<<"LepPhi = "<<LepPhi->at(0)<<setw(12)<<LepPhi->at(1)<<setw(12)<<LepPhi->at(2)<<setw(12)<<LepPhi->at(3)<<endl;     
        cout<<"----------------------------"<<endl;
      }


      // --- H->4mu events
      if(ZZMass >= (125.09-1.1) && ZZMass <= (125.09+1.1)){

        cout<<"----------------------------"<<endl;
        cout<<"H->4mu"<<endl;
        cout<<"RunNumber:LumiNumber:EventNumber = "<<nRun<<":"<<nLumi<<":"<<nEvent<<endl;
        cout<<"Mass   = "<<ZZMass<<endl;
        cout<<"LepID  = "<<LepLepId->at(0)<<setw(12)<<LepLepId->at(1)<<setw(12)<<LepLepId->at(2)<<setw(12)<<LepLepId->at(3)<<endl;
        cout<<"LepPt  = "<<LepPt->at(0)<<setw(12)<<LepPt->at(1)<<setw(12)<<LepPt->at(2)<<setw(12)<<LepPt->at(3)<<endl;
        cout<<"LepEta = "<<LepEta->at(0)<<setw(12)<<LepEta->at(1)<<setw(12)<<LepEta->at(2)<<setw(12)<<LepEta->at(3)<<endl;
        cout<<"LepPhi = "<<LepPhi->at(0)<<setw(12)<<LepPhi->at(1)<<setw(12)<<LepPhi->at(2)<<setw(12)<<LepPhi->at(3)<<endl;     
        cout<<"----------------------------"<<endl;
      }
      

    }
    

  }


}
