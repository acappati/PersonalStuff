// **********************************
// 
// usage:
//     - specify parameters (input/output directory, luminosity) at the end of this file
//     - run with:
//                 root -l -b -q ComputeLeptonScaleSyst_ControlStudyOnMC.C++
//
// **********************************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TFrame.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "RooAbsCollection.h"
#include "RooNumConvPdf.h"
#include "TLorentzVector.h"

#include "Plotter/CMS_lumi.C" // CMS official label definition
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"  // contains DCB definition


using namespace std;
using namespace RooFit ;

#define REDO2lHISTOS 1
#define REDOTHE2lFIT 0
#define COMPUTE2lSCALE 0
#define COMPARE2lDATAMCFIT 0

#define WRITEEXTRATEXTONPLOTS 1 // draw Preliminary on Plots


// *** global definitions
Float_t Zmass_nominalPDG  = 91.1876; // GeV/c^2 (http://pdg.lbl.gov/2017/listings/rpp2017-list-z-boson.pdf)
Float_t Zwidth_nominalPDG = 2.4952;  // GeV/c^2 (http://pdg.lbl.gov/2017/listings/rpp2017-list-z-boson.pdf)
Float_t mass_ele_nominalPDG = 0.0005; //GeV/c^2   
Float_t mass_mu_nominalPDG  = 0.1057; //GeV/c^2


const int nDatasets = 1;
string datasets[nDatasets] = {
    "DYJetsToLL_M50"
  };


// categories of lepton eta
const int nCatEta = 6;
enum CategEta { eleEta1st = 0, eleEta2nd = 1, eleEta3rd = 2, muEta1st = 3, muEta2nd = 4, muEta3rd = 5};
string sCategEta[nCatEta] = {
  "ele_Eta_0-0.8",
  "ele_Eta_0.8-1.5",
  "ele_Eta_1.5-2.5",
  "mu_Eta_0-0.9",
  "mu_Eta_0.9-1.4",
  "mu_Eta_1.4-2.4"
};

// categories of lepton pT 
const int nCatpT = 5;
enum CategpT {pTmin20 = 0, pT2030 = 1, pT3040 = 2, pT4050 = 3, pT50100 = 4};
string sCategpT[nCatpT] = {
  "pT_min-20",
  "pT_20-30",
  "pT_30-40",
  "pT_40-50",
  "pT_50-100"
};

// processes
enum Process {DY=0};

//final state
const int nFinalStates = 3;
enum FinalState {fs4e = 0, fs4mu = 1, fs2e2mu = 2};
string sFinalState[nFinalStates] = {
  "fs_4e",
  "fs_4mu",
  "fs_2e2mu"
};





// *** compute invariant mass function
float compute2lInvMass(float pTLep1, float etaLep1, float phiLep1, int IdLep1, float pTLep2, float etaLep2, float phiLep2)
{

  float mass_lep = 0.;
  // condition on 1st lepton is valid also for the 2nd since Z->ll
  if(int(fabs(IdLep1) == 11)) {mass_lep = mass_ele_nominalPDG;}
  if(int(fabs(IdLep1) == 13)) {mass_lep = mass_mu_nominalPDG;}

  TLorentzVector lep1;
  TLorentzVector lep2;
  lep1.SetPtEtaPhiM(pTLep1 , etaLep1, phiLep1, mass_lep);
  lep2.SetPtEtaPhiM(pTLep2 , etaLep2, phiLep2, mass_lep);

  return (lep1 + lep2).M();

} // end of compute invariant mass function



// *** read file and do histograms (with cut on leading lepton and categories based on subleading lepton pT and eta)
void do2lHistograms_cutsOnLeadingLep(string inputPathMC_DY, string inputPathData, float lumi)
{

  TH1::SetDefaultSumw2(true);

  
  TFile* inputFile[nDatasets];
  TTree* inputTree[nDatasets];
  TH1F* hCounters[nDatasets];
  Double_t gen_sumWeights[nDatasets];
  Float_t partialSampleWeight[nDatasets];

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;

  Short_t Zsel;
  Float_t ZMass;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  vector<Float_t> *LepLepId = 0;
  Float_t xsec;
  Float_t overallEventWeight;

  
  
  // define histos 
  TH1F* hist[nDatasets][nCatEta][nCatpT];
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){
        hist[dat][catEta][catPt] = new TH1F(Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hist_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 120, 60.,120.);
        hist[dat][catEta][catPt]->Sumw2(true);
      }
    }
  }
  

  int currentProcess = DY;
  int currentCategEta = -1;
  int currentCategPt = -1;


  // loop over datasets (only 1 dataset by now: DY MC)
  for(int d=0; d<nDatasets; d++){


    string inputFileName = string(Form("%s%s/ZZ4lAnalysis.root",inputPathMC_DY.c_str(),datasets[d].c_str()));
    inputFile[d] = TFile::Open(inputFileName.c_str());


    hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");    
    gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(40);
    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d];


    inputTree[d] = (TTree*)inputFile[d]->Get("ZTree/candTree");
    inputTree[d]->SetBranchAddress("RunNumber", &nRun);
    inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
    inputTree[d]->SetBranchAddress("LumiNumber", &nLumi);
    inputTree[d]->SetBranchAddress("Zsel", &Zsel);      // WARNING: works with the Z inclusive tree 
    inputTree[d]->SetBranchAddress("ZMass", &ZMass);    // WARNING: works with the Z inclusive tree
    inputTree[d]->SetBranchAddress("LepPt", &LepPt);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);
    inputTree[d]->SetBranchAddress("LepPhi", &LepPhi);
    inputTree[d]->SetBranchAddress("LepLepId", &LepLepId);
    inputTree[d]->SetBranchAddress("xsec", &xsec); 
    inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);

    
    // process tree 
    Long64_t entries = inputTree[d]->GetEntries();
    cout<<"Processing dataset "<<datasets[d]<<" ("<<entries<<" entries) ..."<<endl;
    
    for (Long64_t z=0; z<entries; ++z){
      
      inputTree[d]->GetEntry(z);
      
      if( Zsel < 0. ) continue;  // skip events that do not pass the trigger 
      if(ZMass < 40.) continue;  // applying Z1 mass request 


      // *********************************************************************************************
      // *********************************************************************************************
      // intro known distorsion on the leptons pT  (choose which distorsion introduce)
      
      // no distorsion
      float LepPt0 = LepPt->at(0);
      float LepPt1 = LepPt->at(1);

      // uniform up distorsion of 0.001
      if(false){
        LepPt0 = LepPt->at(0) * (1.+ 0.001);
        LepPt1 = LepPt->at(1) * (1.+ 0.001);
      }
      
      // pT dependent distorsion (stair)
      if(false){
        if(LepPt->at(0) < 20. ) { LepPt0 = LepPt->at(0) * (1.- 0.002);}
        else if(LepPt->at(0) >= 20. && LepPt->at(0) < 30.   ) {LepPt0 = LepPt->at(0) * (1.- 0.001);}
        else if(LepPt->at(0) >= 30. && LepPt->at(0) < 40.   ) {LepPt0 = LepPt->at(0) * (1.+ 0.000);}
        else if(LepPt->at(0) >= 40. && LepPt->at(0) < 50.   ) {LepPt0 = LepPt->at(0) * (1.+ 0.001);}
        else if(LepPt->at(0) >= 50. && LepPt->at(0) <= 100. ) {LepPt0 = LepPt->at(0) * (1.+ 0.002);}
        else {LepPt0 = LepPt->at(0);}

        if(LepPt->at(1) < 20. ) { LepPt1 = LepPt->at(1) * (1.- 0.002);}
        else if(LepPt->at(1) >= 20. && LepPt->at(1) < 30.   ) {LepPt1 = LepPt->at(1) * (1.- 0.001);}
        else if(LepPt->at(1) >= 30. && LepPt->at(1) < 40.   ) {LepPt1 = LepPt->at(1) * (1.+ 0.000);}
        else if(LepPt->at(1) >= 40. && LepPt->at(1) < 50.   ) {LepPt1 = LepPt->at(1) * (1.+ 0.001);}
        else if(LepPt->at(1) >= 50. && LepPt->at(1) <= 100. ) {LepPt1 = LepPt->at(1) * (1.+ 0.002);}
        else {LepPt1 = LepPt->at(1);}
      }
      // **********************************************************************************************
      // **********************************************************************************************


      // CUTS PROPOSED TO REDUCE BKG 
      // (leading lep pT > 20 GeV, ZMass > 40 GeV , categ made according to subleading lep pT and Eta)
      // choose leading lepton to cut on lepton pT (>20 GeV) as in ZZ
      // and subleading lepton to divide events in pT categories 
      float LepPtLeading = 0.;    
      float LepEtaLeading; 
      float LepPhiLeading;
      int LepIDLeading;
      float LepPtSubLeading; 
      float LepEtaSubLeading;
      float LepPhiSubLeading;
      int LepIDSubLeading; 
      
      if(LepPt0 >= LepPt1){
        
        LepPtLeading  = LepPt0;
        LepEtaLeading = LepEta->at(0);
        LepPhiLeading = LepPhi->at(0);
        LepIDLeading  = LepLepId->at(0);
        LepPtSubLeading  = LepPt1;	
        LepEtaSubLeading = LepEta->at(1);
        LepPhiSubLeading = LepPhi->at(1);	
        LepIDSubLeading  = LepLepId->at(1);

      } else if (LepPt0 < LepPt1){
        
        LepPtLeading  = LepPt1;
        LepEtaLeading = LepEta->at(1);
        LepPhiLeading = LepPhi->at(1);
        LepIDLeading  = LepLepId->at(1);
        LepPtSubLeading  = LepPt0;	
        LepEtaSubLeading = LepEta->at(0);
        LepPhiSubLeading = LepPhi->at(0);	
        LepIDSubLeading  = LepLepId->at(0);
      }
      
      if(LepPtLeading < 20.) continue;  
      
      

      // define event weight 
      Double_t eventWeight = partialSampleWeight[d] * xsec * overallEventWeight ;

      
      // Eta categories 
      //Z->ee histos 
      if(int(fabs(LepIDSubLeading)) == 11 ){
      
        if(fabs(LepEtaSubLeading) >= 0. && fabs(LepEtaSubLeading) < 0.8 ) currentCategEta = eleEta1st;
        else if(fabs(LepEtaSubLeading) >= 0.8 && fabs(LepEtaSubLeading) < 1.5 ) currentCategEta = eleEta2nd;
        else if(fabs(LepEtaSubLeading) >= 1.5 && fabs(LepEtaSubLeading) <= 2.5 ) currentCategEta = eleEta3rd;
        else cerr<<"error: wrong eta!"<<endl;
      }

      //Z->mumu histos 
      if(int(fabs(LepIDSubLeading)) == 13 ){
      
        if(fabs(LepEtaSubLeading) >= 0. && fabs(LepEtaSubLeading) < 0.9 ) currentCategEta = muEta1st;
        else if(fabs(LepEtaSubLeading) >= 0.9 && fabs(LepEtaSubLeading) < 1.4 ) currentCategEta = muEta2nd;
        else if(fabs(LepEtaSubLeading) >= 1.4 && fabs(LepEtaSubLeading) <= 2.4 ) currentCategEta = muEta3rd;
        else cerr<<"error: wrong eta!"<<endl;  
      }

      
      // pT categories 
      if(LepPtSubLeading < 20. ) currentCategPt = pTmin20;
      else if(LepPtSubLeading >= 20. && LepPtSubLeading < 30. ) currentCategPt = pT2030;
      else if(LepPtSubLeading >= 30. && LepPtSubLeading < 40. ) currentCategPt = pT3040;
      else if(LepPtSubLeading >= 40. && LepPtSubLeading < 50. ) currentCategPt = pT4050;
      else if(LepPtSubLeading >= 50. && LepPtSubLeading <= 100. ) currentCategPt = pT50100;
      else continue;

      
      
      if(currentCategEta < 0 || currentCategPt < 0) continue;

      
      
      // new mass: pTLep1, etaLep1, phiLep1, IdLep1, pTLep2, etaLep2, phiLep2
      float ZMass_new = compute2lInvMass(LepPtLeading, LepEtaLeading, LepPhiLeading, LepIDLeading, LepPtSubLeading, LepEtaSubLeading, LepPhiSubLeading);
      cout<<ZMass<<" "<<ZMass_new<<endl;

      // fill histos
      hist[currentProcess][currentCategEta][currentCategPt]->Fill(ZMass_new,eventWeight);

      
    } //end loop over tree entries 

    
  }//end loop over datasets


  // write histos in a file 
  TFile* fOutHistos = new TFile("file_DataMCHistos_varStair.root","recreate");
  fOutHistos->cd();
  for(int d=0; d<nDatasets; d++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){
        hist[d][catEta][catPt]->Write(hist[d][catEta][catPt]->GetName());
        delete hist[d][catEta][catPt];
      }
    }
  }
  fOutHistos->Close();
  delete fOutHistos;


}// end doHistograms function 




// perform the fit with DCBxBW (double fit) 
void doThe2lFit_DCBfit(string outputPathFitResultsPlots, string lumiText) 
{

  // read file with histos
  TFile* fInHistos = TFile::Open("file_DataMCHistos_varStair.root");

  // define input histos
  TH1F* inputhist[nDatasets][nCatEta][nCatpT];

  // define histos to store fit results 
  TH1F* hfitResults_poleBW[nDatasets][nCatEta][nCatpT];
  TH1F* hfitResults_widthBW[nDatasets][nCatEta][nCatpT];

  TH1F* hfitResults_meanDCB[nDatasets][nCatEta][nCatpT];
  TH1F* hfitResults_sigmaDCB[nDatasets][nCatEta][nCatpT];
  TH1F* hfitResults_a1DCB[nDatasets][nCatEta][nCatpT];
  TH1F* hfitResults_n1DCB[nDatasets][nCatEta][nCatpT];
  TH1F* hfitResults_a2DCB[nDatasets][nCatEta][nCatpT];
  TH1F* hfitResults_n2DCB[nDatasets][nCatEta][nCatpT];     

  
  // loop over datasets, eta cat and pt cat 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        // take histos from file
        inputhist[dat][catEta][catPt] = (TH1F*)fInHistos->Get(Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        //cout<<inputhist[dat][catEta][catPt]->GetName()<<endl; //debug

        // define histos to store fit results
        hfitResults_poleBW[dat][catEta][catPt] = new TH1F(Form("hfitResults_poleBW_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hfitResults_poleBW_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1);
        hfitResults_widthBW[dat][catEta][catPt] = new TH1F(Form("hfitResults_widthBW_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hfitResults_widthBW_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1);
	hfitResults_meanDCB[dat][catEta][catPt] = new TH1F(Form("hfitResults_meanDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hfitResults_meanDCB_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1);
	hfitResults_sigmaDCB[dat][catEta][catPt] = new TH1F(Form("hfitResults_sigmaDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hfitResults_sigmaDCB_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1); 
	hfitResults_a1DCB[dat][catEta][catPt] = new TH1F(Form("hfitResults_a1DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hfitResults_a1DCB_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1);
	hfitResults_n1DCB[dat][catEta][catPt] = new TH1F(Form("hfitResults_n1DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hfitResults_n1DCB_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1);
	hfitResults_a2DCB[dat][catEta][catPt] = new TH1F(Form("hfitResults_a2DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hfitResults_a2DCB_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1);
	hfitResults_n2DCB[dat][catEta][catPt] = new TH1F(Form("hfitResults_n2DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hfitResults_n2DCB_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1);
        
        

        // define roofit variable for the fit
        RooRealVar mll = RooRealVar("mll","m_{l^{+}l^{-}}",60,120);

        // set fit range
        mll.setRange("range80100gev",80,100); 

        // make roofit datasets from root histos
        RooDataHist dh("dh","dh",mll,Import(*inputhist[dat][catEta][catPt]));

        
        // *** FIT ***  
        // convolution of a BW function and a DoubleCrystalBall resolution function 
        // ---Breit-Wigner function (physical mass distribution)
        RooRealVar pole_BW("pole_BW","pole_BW",Zmass_nominalPDG);
        RooRealVar width_BW("width_BW","width_BW",Zwidth_nominalPDG);

        RooBreitWigner BW_pdf("BW_pdf","Breit-Wigner function",mll,pole_BW,width_BW);

        pole_BW.setConstant(1);    //fix to physical value (setConstant take bool as argument, default is true)
        width_BW.setConstant(1);  //fix to physical value

        // ---Double Crystal Ball function (experimental resolution)
        RooRealVar mean_DCB("mean_DCB","mean_DCB",0.,-2.,2.); 
        RooRealVar sigma_DCB("sigma_DCB","sigma_DCB",1.,0.0001,5.);
        RooRealVar a1_DCB("a1_DCB","a1_DCB",1.,0.,10);
        RooRealVar n1_DCB("n1_DCB","n1_DCB",2.,0.,10);
        RooRealVar a2_DCB("a2_DCB","a2_DCB",1.,0.,10);
        RooRealVar n2_DCB("n2_DCB","n2_DCB",2.,0.,10);
            
        RooDoubleCB DCBall_pdf("DCBall_pdf","Double Crystal ball function",mll,mean_DCB,sigma_DCB,a1_DCB,n1_DCB,a2_DCB,n2_DCB);

        //---Convolution function
        RooNumConvPdf tot_pdf("tot_pdf","Total PDF",mll,DCBall_pdf,BW_pdf);
                


        // do the 1st fit 
        //mean_DCB.setConstant(1);
                
        tot_pdf.chi2FitTo(dh, Range("range80100gev"));
        

        double fitres1[8] = {pole_BW.getVal(), width_BW.getVal(), mean_DCB.getVal(), sigma_DCB.getVal(), a1_DCB.getVal(), n1_DCB.getVal(), a2_DCB.getVal(), n2_DCB.getVal()};
        double fitres1_err[8] = {pole_BW.getError(), width_BW.getError(), mean_DCB.getError(), sigma_DCB.getError(), a1_DCB.getError(), n1_DCB.getError(), a2_DCB.getError(), n2_DCB.getError()};

        // do the 2nd fit 
        //mean_DCB.setConstant(0);
        //a1_DCB.setConstant(1);
        //n1_DCB.setConstant(1);
        //a2_DCB.setConstant(1);
        //n2_DCB.setConstant(1);
                
        //RooNumConvPdf tot_pdf("tot_pdf","Total PDF",mll,DCBall_pdf,BW_pdf);
        tot_pdf.chi2FitTo(dh, Range("range80100gev"));
        

        double fitres2[8] = {pole_BW.getVal(), width_BW.getVal(), mean_DCB.getVal(), sigma_DCB.getVal(), a1_DCB.getVal(), n1_DCB.getVal(), a2_DCB.getVal(), n2_DCB.getVal()};
        double fitres2_err[8] = {pole_BW.getError(), width_BW.getError(), mean_DCB.getError(), sigma_DCB.getError(), a1_DCB.getError(), n1_DCB.getError(), a2_DCB.getError(), n2_DCB.getError()};


                
        // plot data on the frame
        RooPlot* frame = mll.frame();
        frame->SetName(inputhist[dat][catEta][catPt]->GetName());  // name is the name which appears in the root file
        frame->SetTitle(inputhist[dat][catEta][catPt]->GetName()); // title is the title on the canvas  
        dh.plotOn(frame,DataError(RooAbsData::SumW2)); 
        tot_pdf.plotOn(frame, NormRange("range80100gev"), LineColor(kBlue)); //NormRange needed to normalize pdf to data in the fitting range 
        


        // plot on the canvas and save plots
        TCanvas* c = new TCanvas(inputhist[dat][catEta][catPt]->GetName(),inputhist[dat][catEta][catPt]->GetName()); 
        c->cd();
        frame->Draw();

        // draw fit results on canvas
        TPaveText* pv = new TPaveText(0.64,0.52,0.95,0.87,"brNDC");
        pv->AddText(Form("BW pole: %.3f #pm %.3f", pole_BW.getVal(), pole_BW.getError()));
        pv->AddText(Form("BW width: %.3f #pm %.3f", width_BW.getVal(), width_BW.getError()));
        pv->AddText(Form("DCB mean: %.3f #pm %.3f", mean_DCB.getVal(), mean_DCB.getError()));
        pv->AddText(Form("DCB sigma: %.3f #pm %.3f", sigma_DCB.getVal(), sigma_DCB.getError()));
        pv->AddText(Form("DCB a1: %.3f #pm %.3f", a1_DCB.getVal(), a1_DCB.getError()));
        pv->AddText(Form("DCB n1: %.3f #pm %.3f", n1_DCB.getVal(), n1_DCB.getError()));
        pv->AddText(Form("DCB a2: %.3f #pm %.3f", a2_DCB.getVal(), a2_DCB.getError()));
        pv->AddText(Form("DCB n2: %.3f #pm %.3f", n2_DCB.getVal(), n2_DCB.getError()));
        pv->SetFillColor(kWhite);
        pv->SetBorderSize(1);
        pv->SetTextFont(42);
        pv->SetTextSize(0.037);
        pv->SetTextAlign(12); // text left aligned 
        pv->Draw();

        // draw fit results on canvas
        TPaveText* pv1 = new TPaveText(0.10,0.52,0.39,0.87,"brNDC");
        pv1->AddText(Form("BW pole: %.3f #pm %.3f", fitres1[0], fitres1_err[0]));
        pv1->AddText(Form("BW width: %.3f #pm %.3f", fitres1[1], fitres1_err[1]));
        pv1->AddText(Form("DCB mean: %.3f #pm %.3f", fitres1[2], fitres1_err[2]));
        pv1->AddText(Form("DCB sigma: %.3f #pm %.3f", fitres1[3], fitres1_err[3]));
        pv1->AddText(Form("DCB a1: %.3f #pm %.3f", fitres1[4], fitres1_err[4]));
        pv1->AddText(Form("DCB n1: %.3f #pm %.3f", fitres1[5], fitres1_err[5]));
        pv1->AddText(Form("DCB a2: %.3f #pm %.3f", fitres1[6], fitres1_err[6]));
        pv1->AddText(Form("DCB n2: %.3f #pm %.3f", fitres1[7], fitres1_err[7]));
        pv1->SetFillColor(kWhite);
        pv1->SetBorderSize(1);
        pv1->SetTextFont(42);
        pv1->SetTextSize(0.037);
        pv1->SetTextAlign(12); // text left aligned 
        pv1->Draw();

        // print official CMS label and lumi 
        writeExtraText = WRITEEXTRATEXTONPLOTS;
        extraText  = "Preliminary";
        lumi_sqrtS = lumiText + " (13 TeV)";
        cmsTextSize = 0.42;
        lumiTextSize = 0.35;
        extraOverCmsTextSize = 0.72;
        relPosX = 0.12;
        CMS_lumi(c,0,0);


        c->Update();

        c->SaveAs((outputPathFitResultsPlots + "/" + Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()) + ".pdf").c_str());
        c->SaveAs((outputPathFitResultsPlots + "/" + Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()) + ".png").c_str());

        

        // save fit values in histos 
        hfitResults_poleBW[dat][catEta][catPt]->Fill(0.5,pole_BW.getVal());
        hfitResults_poleBW[dat][catEta][catPt]->SetBinError(1,pole_BW.getError());
        hfitResults_widthBW[dat][catEta][catPt]->Fill(0.5,width_BW.getVal());  
        hfitResults_widthBW[dat][catEta][catPt]->SetBinError(1,width_BW.getError());
	                                                  
	hfitResults_meanDCB[dat][catEta][catPt]->Fill(0.5,mean_DCB.getVal());  
        hfitResults_meanDCB[dat][catEta][catPt]->SetBinError(1,mean_DCB.getError());
	hfitResults_sigmaDCB[dat][catEta][catPt]->Fill(0.5,sigma_DCB.getVal()); 
        hfitResults_sigmaDCB[dat][catEta][catPt]->SetBinError(1,sigma_DCB.getError());
	hfitResults_a1DCB[dat][catEta][catPt]->Fill(0.5,a1_DCB.getVal());
        hfitResults_a1DCB[dat][catEta][catPt]->SetBinError(1,a1_DCB.getError());
	hfitResults_n1DCB[dat][catEta][catPt]->Fill(0.5,n1_DCB.getVal());
        hfitResults_n1DCB[dat][catEta][catPt]->SetBinError(1,n1_DCB.getError());
	hfitResults_a2DCB[dat][catEta][catPt]->Fill(0.5,a2_DCB.getVal());
        hfitResults_a2DCB[dat][catEta][catPt]->SetBinError(1,a2_DCB.getError());
	hfitResults_n2DCB[dat][catEta][catPt]->Fill(0.5,n2_DCB.getVal());    
        hfitResults_n2DCB[dat][catEta][catPt]->SetBinError(1,n2_DCB.getError());    


        // // write fit plot frame in a file 
        // fOutFitResults_plotFrame->cd();
        // frame->Write();

      }
    }
  }

  
  // write fit results histos in a file 
  TFile* fOutFitResults = new TFile("file_FitResults_varStair.root","recreate");
  fOutFitResults->cd();
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){
        hfitResults_poleBW[dat][catEta][catPt]->Write(hfitResults_poleBW[dat][catEta][catPt]->GetName());
        hfitResults_widthBW[dat][catEta][catPt]->Write(hfitResults_widthBW[dat][catEta][catPt]->GetName());
        hfitResults_meanDCB[dat][catEta][catPt]->Write(hfitResults_meanDCB[dat][catEta][catPt]->GetName());
        hfitResults_sigmaDCB[dat][catEta][catPt]->Write(hfitResults_sigmaDCB[dat][catEta][catPt]->GetName());
        hfitResults_a1DCB[dat][catEta][catPt]->Write(hfitResults_a1DCB[dat][catEta][catPt]->GetName());
        hfitResults_n1DCB[dat][catEta][catPt]->Write(hfitResults_n1DCB[dat][catEta][catPt]->GetName());
        hfitResults_a2DCB[dat][catEta][catPt]->Write(hfitResults_a2DCB[dat][catEta][catPt]->GetName());
        hfitResults_n2DCB[dat][catEta][catPt]->Write(hfitResults_n2DCB[dat][catEta][catPt]->GetName());
        delete hfitResults_poleBW[dat][catEta][catPt];
        delete hfitResults_widthBW[dat][catEta][catPt];
        delete hfitResults_meanDCB[dat][catEta][catPt];
        delete hfitResults_sigmaDCB[dat][catEta][catPt];
        delete hfitResults_a1DCB[dat][catEta][catPt];
        delete hfitResults_n1DCB[dat][catEta][catPt];
        delete hfitResults_a2DCB[dat][catEta][catPt];
        delete hfitResults_n2DCB[dat][catEta][catPt];
      }
    }
  }
  fOutFitResults->Close();
  delete fOutFitResults;
  
  // fOutFitResults_plotFrame->Close();
  // delete fOutFitResults_plotFrame;


}// end doTheFit DCB with prefit function 




// perform the fit 
void doThe2lFit_gaussFit(string outputPathFitResultsPlots, string lumiText) 
{

  // read file with histos
  TFile* fInHistos = TFile::Open("file_DataMCHistos_varUp_frompT0_varpT001_CUT.root");

  // define input histos
  TH1F* inputhist[nDatasets][nCatEta][nCatpT];

  // define histos to store fit results 
  TH1F* hfitResults_meanDCB[nDatasets][nCatEta][nCatpT];
  TH1F* hfitResults_sigmaDCB[nDatasets][nCatEta][nCatpT];
  
  
  // loop over datasets, eta cat and pt cat 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        // take histos from file
        inputhist[dat][catEta][catPt] = (TH1F*)fInHistos->Get(Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        //cout<<inputhist[dat][catEta][catPt]->GetName()<<endl; //debug

        // define histos to store fit results
        hfitResults_meanDCB[dat][catEta][catPt] = new TH1F(Form("hfitResults_meanDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hfitResults_meanDCB_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1);
	hfitResults_sigmaDCB[dat][catEta][catPt] = new TH1F(Form("hfitResults_sigmaDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("hfitResults_sigmaDCB_%s_%s_%s", datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1); 
    

        // define roofit variable for the fit
        RooRealVar mll = RooRealVar("mll","m_{l^{+}l^{-}}",60,120);

        // make roofit datasets from root histos
        RooDataHist dh("dh","dh",mll,Import(*inputhist[dat][catEta][catPt]));
        

        // **** PREFIT ****
        // set range for the prefit
        mll.setRange("range80100gev",80,100); 

        // define prefit function: gaussian
        RooRealVar mean_prefit_DCB("mean_prefit_DCB","mean_DCB",Zmass_nominalPDG,80.,100.);
        RooRealVar sigma_prefit_DCB("sigma_prefit_DCB","sigma_DCB",Zwidth_nominalPDG,0.0001,5.);
        
        RooGaussian gauss_prefit_pdf("gauss_prefit_pdf","gaussian function",mll,mean_prefit_DCB,sigma_prefit_DCB);

        // do the fit 
        gauss_prefit_pdf.chi2FitTo(dh,Range("range80100gev"));



        // **** FIT ****
        // define fit range 
        float fitRange_min = mean_prefit_DCB.getVal() - sigma_prefit_DCB.getVal();
        float fitRange_max = mean_prefit_DCB.getVal() + sigma_prefit_DCB.getVal();

        // set range for the fit
	mll.setRange("range_fit", fitRange_min, fitRange_max); 
        
        // define fit function: gaussian
        RooRealVar mean_DCB("mean_DCB","mean_DCB",Zmass_nominalPDG, fitRange_min, fitRange_max);
        RooRealVar sigma_DCB("sigma_DCB","sigma_DCB",Zwidth_nominalPDG, 0.0001, 5.);
        
        RooGaussian gauss_pdf("gauss_pdf","gaussian function",mll,mean_DCB,sigma_DCB);

        // do the fit 
        gauss_pdf.chi2FitTo(dh,Range("range_fit"));
                


        // plot data on the frame
        RooPlot* frame = mll.frame();
        frame->SetName(inputhist[dat][catEta][catPt]->GetName());  // name is the name which appears in the root file
        frame->SetTitle(inputhist[dat][catEta][catPt]->GetName()); // title is the title on the canvas  

        dh.plotOn(frame,DataError(RooAbsData::SumW2));
        gauss_prefit_pdf.plotOn(frame,LineColor(kBlue)); 
        gauss_pdf.plotOn(frame,LineColor(kRed)); 



        // plot on the canvas and save plots
        TCanvas* c = new TCanvas(inputhist[dat][catEta][catPt]->GetName(),inputhist[dat][catEta][catPt]->GetName()); 
        c->cd();
        frame->Draw();

        // draw fit results on canvas
        TPaveText* pv = new TPaveText(0.64,0.65,0.95,0.87,"brNDC");
        pv->AddText(Form("mean prefit: %.3f #pm %.3f", mean_prefit_DCB.getVal(), mean_prefit_DCB.getError())); ((TText*)pv->GetListOfLines()->Last())->SetTextColor(kBlue);
        pv->AddText(Form("sigma prefit: %.3f #pm %.3f", sigma_prefit_DCB.getVal(), sigma_prefit_DCB.getError())); ((TText*)pv->GetListOfLines()->Last())->SetTextColor(kBlue);
        pv->AddText(Form("mean fit: %.3f #pm %.3f", mean_DCB.getVal(), mean_DCB.getError())); ((TText*)pv->GetListOfLines()->Last())->SetTextColor(kRed);
        pv->AddText(Form("sigma fit: %.3f #pm %.3f", sigma_DCB.getVal(), sigma_DCB.getError())); ((TText*)pv->GetListOfLines()->Last())->SetTextColor(kRed);
        pv->SetFillColor(kWhite);
        pv->SetBorderSize(1);
        pv->SetTextFont(42);
        pv->SetTextSize(0.037);
        pv->SetTextAlign(12); // text left aligned 
        pv->Draw();

        // print official CMS label and lumi 
        writeExtraText = WRITEEXTRATEXTONPLOTS;
        extraText  = "Preliminary";
        lumi_sqrtS = lumiText + " (13 TeV)";
        cmsTextSize = 0.42;
        lumiTextSize = 0.35;
        extraOverCmsTextSize = 0.72;
        relPosX = 0.12;
        CMS_lumi(c,0,0);


        c->Update();

        c->SaveAs((outputPathFitResultsPlots + "/" + Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()) + ".pdf").c_str());
        c->SaveAs((outputPathFitResultsPlots + "/" + Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()) + ".png").c_str());


        // save fit values in histos 
      	hfitResults_meanDCB[dat][catEta][catPt]->Fill(0.5,mean_DCB.getVal());  
        hfitResults_meanDCB[dat][catEta][catPt]->SetBinError(1,mean_DCB.getError());
	hfitResults_sigmaDCB[dat][catEta][catPt]->Fill(0.5,sigma_DCB.getVal()); 
        hfitResults_sigmaDCB[dat][catEta][catPt]->SetBinError(1,sigma_DCB.getError());
	
        // // write fit plot frame in a file 
        // fOutFitResults_plotFrame->cd();
        // frame->Write();

      }
    }
  }

  
  // write fit results histos in a file 
  TFile* fOutFitResults = new TFile("file_FitResults_varUp_frompT0_varpT001_CUT_chi2fit_gauss.root","recreate");
  fOutFitResults->cd();
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){
        hfitResults_meanDCB[dat][catEta][catPt]->Write(hfitResults_meanDCB[dat][catEta][catPt]->GetName());
        hfitResults_sigmaDCB[dat][catEta][catPt]->Write(hfitResults_sigmaDCB[dat][catEta][catPt]->GetName());
        delete hfitResults_meanDCB[dat][catEta][catPt];
        delete hfitResults_sigmaDCB[dat][catEta][catPt];
      }
    }
  }
  fOutFitResults->Close();
  delete fOutFitResults;
  
  // fOutFitResults_plotFrame->Close();
  // delete fOutFitResults_plotFrame;


}// end doTheFit_gaussFit function 



// take fit results and compute dilepton scale 
void computeDileptonScale(string outputPathDileptonScalePlots, string lumiText)
{
  
 //***********************************************
 // *** read file with fit results nominal
  TFile* fInFitResults1 = TFile::Open("file_FitResults_nominal.root");

  // define input histos 
  TH1F* hinput_meanFitResults1[nDatasets][nCatEta][nCatpT]; 

  // read histos with fit results from file 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        hinput_meanFitResults1[dat][catEta][catPt] = (TH1F*)fInFitResults1->Get(Form("hfitResults_meanDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        //cout<<hinput_meanFitResults1[dat][catEta][catPt]->GetName()<<endl; // debug
        
      }
    }
  }

  //***********************************************
  // *** read file with fit results with variations
  TFile* fInFitResults2 = TFile::Open("file_FitResults_varStair.root");

  // define input histos 
  TH1F* hinput_meanFitResults2[nDatasets][nCatEta][nCatpT];  

  // read histos with fit results from file 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        hinput_meanFitResults2[dat][catEta][catPt] = (TH1F*)fInFitResults2->Get(Form("hfitResults_meanDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        //cout<<hinput_meanFitResults2[dat][catEta][catPt]->GetName()<<endl; // debug
        
      }
    }
  }



  // dilepton scale vector and histos
  Float_t vec_dileptonScale[nCatEta][nCatpT]; 
  Float_t vec_dileptonScale_err[nCatEta][nCatpT]; 
  TH1F* h_dileptonScale[nCatEta][nCatpT]; 

  // compute dilepton scale 
  for(int catEta=0; catEta<nCatEta; catEta++){
    for(int catPt=0; catPt<nCatpT; catPt++){

      vec_dileptonScale[catEta][catPt] = (hinput_meanFitResults2[0][catEta][catPt]->GetBinContent(1) - hinput_meanFitResults1[0][catEta][catPt]->GetBinContent(1) )/ Zmass_nominalPDG;  // MC_mean2 - MC_mean1 / MassZ_PDG

      vec_dileptonScale_err[catEta][catPt] = TMath::Power((TMath::Power(hinput_meanFitResults2[0][catEta][catPt]->GetBinError(1),2) + TMath::Power(hinput_meanFitResults1[0][catEta][catPt]->GetBinError(1),2)), 0.5) / Zmass_nominalPDG;  


      h_dileptonScale[catEta][catPt] = new TH1F(Form("h_dileptonScale_%s_%s",sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("h_dileptonScale_%s_%s",sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1);
      //h_dileptonScale[catEta][catPt]->Sumw2(true);        

      h_dileptonScale[catEta][catPt]->Fill(0.5,vec_dileptonScale[catEta][catPt]);
      h_dileptonScale[catEta][catPt]->SetBinError(1,vec_dileptonScale_err[catEta][catPt]);

      //cout<<vec_dileptonScale[catEta][catPt]<<" "<<vec_dileptonScale_err[catEta][catPt]<<endl;
      cout<<sCategEta[catEta]<<" "<<sCategpT[catPt]<<" "<<hinput_meanFitResults2[0][catEta][catPt]->GetBinContent(1)<<" "<<hinput_meanFitResults1[0][catEta][catPt]->GetBinContent(1)<<endl;
      cout<<sCategEta[catEta]<<" "<<sCategpT[catPt]<<" "<<hinput_meanFitResults2[0][catEta][catPt]->GetBinError(1)<<" "<<hinput_meanFitResults1[0][catEta][catPt]->GetBinError(1)<<" "<<Zmass_nominalPDG<<endl;
    }
  }
  
  // --- dilepton scale plots 

  // define x axis values: mean value of each pT bin 
  Float_t vec_ptValues_ele[nCatpT] = {13.5, 25., 35., 45., 75.}; // first bin for ele: 7-20 GeV
  Float_t vec_ptValues_mu[nCatpT] = {12.5, 25., 35., 45., 75.};  // first bin for mu: 5-20 GeV
  // define x axis uncertainties
  Float_t vec_ptValues_ele_err[nCatpT] = {6.5, 5., 5., 5., 25.};
  Float_t vec_ptValues_mu_err[nCatpT] = {7.5, 5., 5., 5., 25.};

  // define y axis values: dilepton scale for each Eta bin 
  Float_t vec_dileptonScale_eleCatEta0[nCatpT];
  Float_t vec_dileptonScale_eleCatEta1[nCatpT];
  Float_t vec_dileptonScale_eleCatEta2[nCatpT];
  Float_t vec_dileptonScale_muCatEta3[nCatpT];
  Float_t vec_dileptonScale_muCatEta4[nCatpT];
  Float_t vec_dileptonScale_muCatEta5[nCatpT];
  // define y axis uncertainties
  Float_t vec_dileptonScale_eleCatEta0_err[nCatpT];
  Float_t vec_dileptonScale_eleCatEta1_err[nCatpT];
  Float_t vec_dileptonScale_eleCatEta2_err[nCatpT];
  Float_t vec_dileptonScale_muCatEta3_err[nCatpT];
  Float_t vec_dileptonScale_muCatEta4_err[nCatpT];
  Float_t vec_dileptonScale_muCatEta5_err[nCatpT];

  for(int catPt=0; catPt<nCatpT; catPt++){

    // y values
    vec_dileptonScale_eleCatEta0[catPt] = vec_dileptonScale[0][catPt];
    vec_dileptonScale_eleCatEta1[catPt] = vec_dileptonScale[1][catPt];
    vec_dileptonScale_eleCatEta2[catPt] = vec_dileptonScale[2][catPt];
    vec_dileptonScale_muCatEta3[catPt]  = vec_dileptonScale[3][catPt];
    vec_dileptonScale_muCatEta4[catPt]  = vec_dileptonScale[4][catPt];
    vec_dileptonScale_muCatEta5[catPt]  = vec_dileptonScale[5][catPt];

    // y uncertainties
    vec_dileptonScale_eleCatEta0_err[catPt] = vec_dileptonScale_err[0][catPt];
    vec_dileptonScale_eleCatEta1_err[catPt] = vec_dileptonScale_err[1][catPt];
    vec_dileptonScale_eleCatEta2_err[catPt] = vec_dileptonScale_err[2][catPt];
    vec_dileptonScale_muCatEta3_err[catPt]  = vec_dileptonScale_err[3][catPt];
    vec_dileptonScale_muCatEta4_err[catPt]  = vec_dileptonScale_err[4][catPt];
    vec_dileptonScale_muCatEta5_err[catPt]  = vec_dileptonScale_err[5][catPt];

  }

    
  // electron plots
  TGraphErrors* graph_eleCatEta0 = new TGraphErrors(nCatpT,vec_ptValues_ele,vec_dileptonScale_eleCatEta0,vec_ptValues_ele_err,vec_dileptonScale_eleCatEta0_err); 
  graph_eleCatEta0->SetMarkerColor(kBlue);
  graph_eleCatEta0->SetMarkerStyle(25);
  graph_eleCatEta0->SetMarkerSize(0.7);
  TGraphErrors* graph_eleCatEta1 = new TGraphErrors(nCatpT,vec_ptValues_ele,vec_dileptonScale_eleCatEta1,vec_ptValues_ele_err,vec_dileptonScale_eleCatEta1_err);
  graph_eleCatEta1->SetMarkerColor(kBlack);
  graph_eleCatEta1->SetMarkerStyle(25);
  graph_eleCatEta1->SetMarkerSize(0.7);
  TGraphErrors* graph_eleCatEta2 = new TGraphErrors(nCatpT,vec_ptValues_ele,vec_dileptonScale_eleCatEta2,vec_ptValues_ele_err,vec_dileptonScale_eleCatEta2_err);
  graph_eleCatEta2->SetMarkerColor(kRed);
  graph_eleCatEta2->SetMarkerStyle(25);
  graph_eleCatEta2->SetMarkerSize(0.7);

  TMultiGraph* mg_ele = new TMultiGraph();
  mg_ele->Add(graph_eleCatEta0);
  mg_ele->Add(graph_eleCatEta1);
  mg_ele->Add(graph_eleCatEta2);

  //mg_ele->SetMinimum(-0.02);
  //mg_ele->SetMaximum(0.02);

  TCanvas* can_ele = new TCanvas("can_ele","can_ele",600,600);
  mg_ele->Draw("ap"); 
  mg_ele->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
  mg_ele->GetXaxis()->SetTitleFont(43);
  mg_ele->GetXaxis()->SetTitleSize(18);
  mg_ele->GetXaxis()->SetTitleOffset(1.3);
  mg_ele->GetXaxis()->SetLabelFont(43);
  mg_ele->GetXaxis()->SetLabelSize(13);
  mg_ele->GetXaxis()->SetLimits(0.,100.); //set x axis limits (same as SetRangeUser but for TGraphs)
  mg_ele->GetYaxis()->SetTitle("(m^{peak}_{var} - m^{peak}_{nominal})/m_{PDG}");
  mg_ele->GetYaxis()->SetTitleFont(43);
  mg_ele->GetYaxis()->SetTitleSize(18);
  mg_ele->GetYaxis()->SetTitleOffset(1.5);
  mg_ele->GetYaxis()->SetLabelFont(43);
  mg_ele->GetYaxis()->SetLabelSize(13);
  gPad->SetTickx(); // set ticks also on upper x axis
  gPad->SetTicky(); // set ticks also on right y axis

  TLegend* legend_ele = new TLegend(0.58,0.15,0.8,0.35);
  legend_ele->AddEntry(graph_eleCatEta0,"Z,|#eta| 0.0-0.8","lp");
  legend_ele->AddEntry(graph_eleCatEta1,"Z,|#eta| 0.8-1.5","lp");
  legend_ele->AddEntry(graph_eleCatEta2,"Z,|#eta| 1.5-2.5","lp");
  legend_ele->SetTextFont(43);
  legend_ele->SetTextSize(14);
  legend_ele->SetLineColor(kWhite);
  legend_ele->Draw();

  // 0 level line
  TLine* line0_ele = new TLine(0,0,100,0);
  line0_ele->SetLineColor(kBlack);
  line0_ele->SetLineStyle(2); // dashed line
  line0_ele->Draw();

  // print official CMS label and lumi 
  writeExtraText = WRITEEXTRATEXTONPLOTS;
  extraText  = "Preliminary";
  lumi_sqrtS = lumiText + " (13 TeV)";
  cmsTextSize = 0.42;
  lumiTextSize = 0.35;
  extraOverCmsTextSize = 0.72;
  relPosX = 0.12;
  CMS_lumi(can_ele,0,0);


  can_ele->Update();

  can_ele->SaveAs((outputPathDileptonScalePlots + "/DileptonScale_electrons.png").c_str());
  can_ele->SaveAs((outputPathDileptonScalePlots + "/DileptonScale_electrons.pdf").c_str());

   

  // muon plots
  TGraphErrors* graph_muCatEta3 = new TGraphErrors(nCatpT,vec_ptValues_mu,vec_dileptonScale_muCatEta3,vec_ptValues_mu_err,vec_dileptonScale_muCatEta3_err); 
  graph_muCatEta3->SetMarkerColor(kBlue);
  graph_muCatEta3->SetMarkerStyle(25);
  graph_muCatEta3->SetMarkerSize(0.7);
  TGraphErrors* graph_muCatEta4 = new TGraphErrors(nCatpT,vec_ptValues_mu,vec_dileptonScale_muCatEta4,vec_ptValues_mu_err,vec_dileptonScale_muCatEta4_err);
  graph_muCatEta4->SetMarkerColor(kBlack);
  graph_muCatEta4->SetMarkerStyle(25);
  graph_muCatEta4->SetMarkerSize(0.7);
  TGraphErrors* graph_muCatEta5 = new TGraphErrors(nCatpT,vec_ptValues_mu,vec_dileptonScale_muCatEta5,vec_ptValues_mu_err,vec_dileptonScale_muCatEta5_err);
  graph_muCatEta5->SetMarkerColor(kRed);
  graph_muCatEta5->SetMarkerStyle(25);
  graph_muCatEta5->SetMarkerSize(0.7);

  TMultiGraph* mg_mu = new TMultiGraph();
  mg_mu->Add(graph_muCatEta3);
  mg_mu->Add(graph_muCatEta4);
  mg_mu->Add(graph_muCatEta5);

  //mg_mu->SetMinimum(-0.02);
  //mg_mu->SetMaximum(0.02);
 
  TCanvas* can_mu = new TCanvas("can_mu","can_mu",600,600);
  mg_mu->Draw("ap");
  mg_mu->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
  mg_mu->GetXaxis()->SetTitleFont(43);
  mg_mu->GetXaxis()->SetTitleSize(18);
  mg_mu->GetXaxis()->SetTitleOffset(1.3);
  mg_mu->GetXaxis()->SetLabelFont(43);
  mg_mu->GetXaxis()->SetLabelSize(13);
  mg_mu->GetXaxis()->SetLimits(0.,100.); //set x axis limits (same as SetRangeUser but for TGraphs)
  mg_mu->GetYaxis()->SetTitle("(m^{peak}_{var} - m^{peak}_{nominal})/m_{PDG}");
  mg_mu->GetYaxis()->SetTitleFont(43);
  mg_mu->GetYaxis()->SetTitleSize(18);
  mg_mu->GetYaxis()->SetTitleOffset(1.5);
  mg_mu->GetYaxis()->SetLabelFont(43);
  mg_mu->GetYaxis()->SetLabelSize(13);
  gPad->SetTickx(); // set ticks also on upper x axis
  gPad->SetTicky(); // set ticks also on right y axis

  TLegend* legend_mu = new TLegend(0.58,0.15,0.8,0.35);
  legend_mu->AddEntry(graph_eleCatEta0,"Z,|#eta| 0.0-0.9","lp");
  legend_mu->AddEntry(graph_eleCatEta1,"Z,|#eta| 0.9-1.4","lp");
  legend_mu->AddEntry(graph_eleCatEta2,"Z,|#eta| 1.4-2.4","lp");
  legend_mu->SetTextFont(43);
  legend_mu->SetTextSize(14);
  legend_mu->SetLineColor(kWhite);
  legend_mu->Draw();

  // 0 level line
  TLine* line0_mu = new TLine(0,0,100,0);
  line0_mu->SetLineColor(kBlack);
  line0_mu->SetLineStyle(2); // dashed line
  line0_mu->Draw();

  // print official CMS label and lumi 
  writeExtraText = WRITEEXTRATEXTONPLOTS;
  extraText  = "Preliminary";
  lumi_sqrtS = lumiText + " (13 TeV)";
  cmsTextSize = 0.42;
  lumiTextSize = 0.35;
  extraOverCmsTextSize = 0.72;
  relPosX = 0.12;
  CMS_lumi(can_mu,0,0);
  

  can_mu->Update();

  can_mu->SaveAs((outputPathDileptonScalePlots + "/DileptonScale_muons.png").c_str());
  can_mu->SaveAs((outputPathDileptonScalePlots + "/DileptonScale_muons.pdf").c_str());


}// end computeDileptonScale function



// comparison between Data and MC 2l fit
void compareDataMCfitPlots(string outputPathCompare2lDataMcFit, string lumiText)
{

  //***********************************************
  // *** read file with fit results nominal
  TFile* fIn_nominal = TFile::Open("file_DataMCHistos_nominal.root");

  // define input histos 
  TH1F* hinput_nominal[nDatasets][nCatEta][nCatpT];   

  // read histos with fit results from file 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        hinput_nominal[dat][catEta][catPt] = (TH1F*)fIn_nominal->Get(Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        //cout<<hinput_nominal[dat][catEta][catPt]->GetName()<<endl; // debug

      }
    }
  }

  //***********************************************
  // *** read file with fit results with variations
  TFile* fIn_var = TFile::Open("file_DataMCHistos_varStair.root");

  // define input histos 
  TH1F* hinput_var[nDatasets][nCatEta][nCatpT];   

  // read histos with fit results from file 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        hinput_var[dat][catEta][catPt] = (TH1F*)fIn_var->Get(Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        //cout<<hinput_var[dat][catEta][catPt]->GetName()<<endl; // debug

      }
    }
  }

  
  //plot
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        TCanvas *can = new TCanvas(hinput_var[dat][catEta][catPt]->GetName(),hinput_var[dat][catEta][catPt]->GetName());
        can->cd();

        hinput_nominal[dat][catEta][catPt]->SetLineColor(kBlue);
        hinput_nominal[dat][catEta][catPt]->Draw("hist");

        hinput_var[dat][catEta][catPt]->SetLineColor(kRed);
        hinput_var[dat][catEta][catPt]->Draw("samel");

        gStyle->SetOptStat(0);

        
        TLegend* legend = new TLegend(0.74,0.68,0.94,0.87);
        legend->AddEntry(hinput_nominal[dat][catEta][catPt],"nominal","f");
        legend->AddEntry(hinput_var[dat][catEta][catPt],"var","l");
        legend->SetTextFont(43);
        legend->SetTextSize(14);
        legend->SetLineColor(kBlack);
        legend->SetFillColor(kWhite);
        legend->Draw();

        
        // print official CMS label and lumi 
        writeExtraText = WRITEEXTRATEXTONPLOTS;
        extraText  = "Preliminary";
        lumi_sqrtS = lumiText + " (13 TeV)";
        cmsTextSize = 0.42;
        lumiTextSize = 0.35;
        extraOverCmsTextSize = 0.72;
        relPosX = 0.12;
        CMS_lumi(can,0,0);


        can->Update();
        can->SaveAs((outputPathCompare2lDataMcFit + "/" + Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()) + ".pdf").c_str());
        can->SaveAs((outputPathCompare2lDataMcFit + "/" + Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()) + ".png").c_str());
        

      }
    }
  }
  


}// end compareDataMCfitPlots function









// *** main function
void ComputeLeptonScaleSyst_ControlStudyOnMC()
{
 
  string inputPathMC_DY = "/data3/Higgs/180416/MC_main/";
  string inputPathData = "/data3/Higgs/180416/";
  string inputPathMC_ggH = "/data3/Higgs/180416/MC_main/";

  string outputPathFitResultsPlots = "plotsSysts_FitResults_varStair";
  string outputPathDileptonScalePlots = "plotsSysts_DileptonScale_varStair";
  string outputPathCompare2lDataMcFit = "plotsSysts_CompareDataMC2lFit_varStair";
  
  

  float lumi = 41.30; //fb-1
  string lumiText = "41.30 fb^{-1}";


  // create output directories
  if(REDOTHE2lFIT) gSystem->Exec(("mkdir -p "+outputPathFitResultsPlots).c_str());  //dir for fit results plots

  if(COMPUTE2lSCALE) gSystem->Exec(("mkdir -p "+outputPathDileptonScalePlots).c_str()); //dir for dilepton scale plots

  if(COMPARE2lDATAMCFIT) gSystem->Exec(("mkdir -p "+outputPathCompare2lDataMcFit).c_str()); //dir for comparison between Data and MC 2l fit

 


  // execute functions 
  if(REDO2lHISTOS) do2lHistograms_cutsOnLeadingLep(inputPathMC_DY, inputPathData, lumi); //read MC file and do histos separating events in categories
                                                                                         //cut leading lep pT>20 GeV; categories based on subleading lepton pT and eta

  if(REDOTHE2lFIT) doThe2lFit_DCBfit(outputPathFitResultsPlots, lumiText); // do fit with DCBxBW 

  //if(REDOTHE2lFIT) doThe2lFit_gaussFit(outputPathFitResultsPlots, lumiText);  // do gaussian fit

  if(COMPUTE2lSCALE) computeDileptonScale(outputPathDileptonScalePlots, lumiText);

  if(COMPARE2lDATAMCFIT) compareDataMCfitPlots(outputPathCompare2lDataMcFit, lumiText);

  

  
}



