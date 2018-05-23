// **********************************
// 
// usage:
//     - specify parameters (input/output directory, luminosity) at the end of this file
//     - specify the sample type (nominal, uniform variation up of 001, stair variation) at the beginning of the file 
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
using namespace RooFit;


#define WRITEEXTRATEXTONPLOTS 1 // draw Preliminary on Plots


// ********************************************
// *** bool for choosing sample type
// *** choose their value to change sample type
bool sampletype_nominal = true;
bool sampletype_varUnif001 = false;
bool sampletype_varStair = false;

// ********************************************


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




// *** read file and do histograms (events separated into categories based on pT and Eta of 1 of the 2 leptons, 
//     determined randomly, and integrating over the other) as explained in AN2016_442 (Section 9)
void do2lHistograms_AN(string inputPathMC_DY, string inputPathData, float lumi, string sampletype_name)
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
      

      // *********************************************************************************************
      // *********************************************************************************************
      // intro known distorsion on the leptons pT  (choose which distorsion introduce)
      
      float LepPt0; 
      float LepPt1; 
      
      if(sampletype_varUnif001){        // uniform up distorsion of 0.001
        LepPt0 = LepPt->at(0) * (1.+ 0.001);
        LepPt1 = LepPt->at(1) * (1.+ 0.001);
      }
      else if(sampletype_varStair){     // pT dependent distorsion (stair)
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
      else{     // no distorsion
        LepPt0 = LepPt->at(0);
	LepPt1 = LepPt->at(1);
      }
      // **********************************************************************************************
      // **********************************************************************************************



      // define event weight 
      Double_t eventWeight = partialSampleWeight[d] * xsec * overallEventWeight ;

      
      // Eta categories 
      //Z->ee histos 
      if(int(fabs(LepLepId->at(0))) == 11 ){
      
        if(fabs(LepEta->at(0)) >= 0. && fabs(LepEta->at(0)) < 0.8 ) currentCategEta = eleEta1st;
        else if(fabs(LepEta->at(0)) >= 0.8 && fabs(LepEta->at(0)) < 1.5 ) currentCategEta = eleEta2nd;
        else if(fabs(LepEta->at(0)) >= 1.5 && fabs(LepEta->at(0)) <= 2.5 ) currentCategEta = eleEta3rd;
        else cerr<<"error: wrong eta!"<<endl;
      }

      //Z->mumu histos 
      if(int(fabs(LepLepId->at(0))) == 13 ){
      
        if(fabs(LepEta->at(0)) >= 0. && fabs(LepEta->at(0)) < 0.9 ) currentCategEta = muEta1st;
        else if(fabs(LepEta->at(0)) >= 0.9 && fabs(LepEta->at(0)) < 1.4 ) currentCategEta = muEta2nd;
        else if(fabs(LepEta->at(0)) >= 1.4 && fabs(LepEta->at(0)) <= 2.4 ) currentCategEta = muEta3rd;
        else cerr<<"error: wrong eta!"<<endl;  
      }

      
      // pT categories 
      if(LepPt0 < 20. ) currentCategPt = pTmin20;
      else if(LepPt0 >= 20. && LepPt0 < 30. ) currentCategPt = pT2030;
      else if(LepPt0 >= 30. && LepPt0 < 40. ) currentCategPt = pT3040;
      else if(LepPt0 >= 40. && LepPt0 < 50. ) currentCategPt = pT4050;
      else if(LepPt0 >= 50. && LepPt0 <= 100. ) currentCategPt = pT50100;
      else continue;

      
      
      if(currentCategEta < 0 || currentCategPt < 0) continue;

      
      
      // new mass: pTLep1, etaLep1, phiLep1, IdLep1, pTLep2, etaLep2, phiLep2
      float ZMass_new = compute2lInvMass(LepPt0, LepEta->at(0), LepPhi->at(0), LepLepId->at(0), LepPt1, LepEta->at(1), LepPhi->at(1));
      cout<<ZMass<<" "<<ZMass_new<<endl;

      // fill histos
      hist[currentProcess][currentCategEta][currentCategPt]->Fill(ZMass_new,eventWeight);

      
    } //end loop over tree entries 

    
  }//end loop over datasets


  // write histos in a file 
  TFile* fOutHistos = new TFile(string(Form("file_DataMCHistos_%s.root",sampletype_name.c_str())).c_str(),"recreate");
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


} // end doHistograms_AN function 




// perform the fit with DCBxBW (double fit) 
void doThe2lFit_DCBfit(string outputPathFitResultsPlots, string lumiText, string sampletype_name) 
{

  // define output file for fit result plots 
  TFile* fOutFitResults_plotFrame = new TFile(string(Form("file_FitResultsPlots_%s.root",sampletype_name.c_str())).c_str(),"recreate");

  // read file with histos
  TFile* fInHistos = TFile::Open(string(Form("file_DataMCHistos_%s.root",sampletype_name.c_str())).c_str());

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
        RooRealVar sigma_DCB("sigma_DCB","sigma_DCB",1.,0.0001,10.);
        RooRealVar a1_DCB("a1_DCB","a1_DCB",1.,0.,50);
        RooRealVar n1_DCB("n1_DCB","n1_DCB",2.,0.,50);
        RooRealVar a2_DCB("a2_DCB","a2_DCB",1.,0.,50);
        RooRealVar n2_DCB("n2_DCB","n2_DCB",2.,0.,50);
            
        RooDoubleCB DCBall_pdf("DCBall_pdf","Double Crystal ball function",mll,mean_DCB,sigma_DCB,a1_DCB,n1_DCB,a2_DCB,n2_DCB);

        //---Convolution function 
        RooNumConvPdf tot_pdf("tot_pdf","Total PDF",mll,DCBall_pdf,BW_pdf);
                

        

        // do the 1st fit                 
        tot_pdf.chi2FitTo(dh, Range("range80100gev"));
        
        double fitres1[8] = {pole_BW.getVal(), width_BW.getVal(), mean_DCB.getVal(), sigma_DCB.getVal(), a1_DCB.getVal(), n1_DCB.getVal(), a2_DCB.getVal(), n2_DCB.getVal()};
        double fitres1_err[8] = {pole_BW.getError(), width_BW.getError(), mean_DCB.getError(), sigma_DCB.getError(), a1_DCB.getError(), n1_DCB.getError(), a2_DCB.getError(), n2_DCB.getError()};

        // do the 2nd fit 
        tot_pdf.chi2FitTo(dh, Range("range80100gev"));
       
        double fitres2[8] = {pole_BW.getVal(), width_BW.getVal(), mean_DCB.getVal(), sigma_DCB.getVal(), a1_DCB.getVal(), n1_DCB.getVal(), a2_DCB.getVal(), n2_DCB.getVal()};
        double fitres2_err[8] = {pole_BW.getError(), width_BW.getError(), mean_DCB.getError(), sigma_DCB.getError(), a1_DCB.getError(), n1_DCB.getError(), a2_DCB.getError(), n2_DCB.getError()};

               

                
        // plot data on the frame
        RooPlot* frame = mll.frame();
        frame->SetName(inputhist[dat][catEta][catPt]->GetName());  // name is the name which appears in the root file
        frame->SetTitle(""); // title is the title on the canvas  
        dh.plotOn(frame,DataError(RooAbsData::SumW2), MarkerStyle(sampletype_nominal?kOpenCircle:kFullCircle)); //open circle marker for nominal distrib 
        tot_pdf.plotOn(frame, NormRange("range80100gev"), LineColor(sampletype_nominal?kBlue:kRed)); //NormRange to normalize pdf to data in the fitting range 
                                                                                                     //blue line for nominal distrib
                                                                                                     
        

        // plot on the canvas and save plots
        TCanvas* c = new TCanvas(inputhist[dat][catEta][catPt]->GetName(),inputhist[dat][catEta][catPt]->GetName()); 
        c->cd();
        frame->Draw();

        // draw 1st fit results on canvas
        TPaveText* pv1 = new TPaveText(0.10,0.51,0.39,0.88,"brNDC");
        pv1->AddText("1st fit: ");
        pv1->AddText(Form("BW pole: %.3f #pm %.3f",   fitres1[0], fitres1_err[0]));
        pv1->AddText(Form("BW width: %.3f #pm %.3f",  fitres1[1], fitres1_err[1]));
        pv1->AddText(Form("DCB mean: %.3f #pm %.3f",  fitres1[2], fitres1_err[2]));
        pv1->AddText(Form("DCB sigma: %.3f #pm %.3f", fitres1[3], fitres1_err[3]));
        pv1->AddText(Form("DCB a1: %.3f #pm %.3f",    fitres1[4], fitres1_err[4]));
        pv1->AddText(Form("DCB n1: %.3f #pm %.3f",    fitres1[5], fitres1_err[5]));
        pv1->AddText(Form("DCB a2: %.3f #pm %.3f",    fitres1[6], fitres1_err[6]));
        pv1->AddText(Form("DCB n2: %.3f #pm %.3f",    fitres1[7], fitres1_err[7]));
        pv1->SetFillColor(kWhite);
        pv1->SetBorderSize(1);
        pv1->SetTextFont(42);
        pv1->SetTextSize(0.037);
        pv1->SetTextAlign(12); // text left aligned 
        pv1->Draw();

        // draw 2nd fit results on canvas
        TPaveText* pv2 = new TPaveText(0.66,0.51,0.95,0.88,"brNDC");
        pv2->AddText("2nd fit: ");
        pv2->AddText(Form("BW pole: %.3f #pm %.3f",   fitres2[0], fitres2_err[0]));
        pv2->AddText(Form("BW width: %.3f #pm %.3f",  fitres2[1], fitres2_err[1]));
        pv2->AddText(Form("DCB mean: %.3f #pm %.3f",  fitres2[2], fitres2_err[2]));
        pv2->AddText(Form("DCB sigma: %.3f #pm %.3f", fitres2[3], fitres2_err[3]));
        pv2->AddText(Form("DCB a1: %.3f #pm %.3f",    fitres2[4], fitres2_err[4]));
        pv2->AddText(Form("DCB n1: %.3f #pm %.3f",    fitres2[5], fitres2_err[5]));
        pv2->AddText(Form("DCB a2: %.3f #pm %.3f",    fitres2[6], fitres2_err[6]));
        pv2->AddText(Form("DCB n2: %.3f #pm %.3f",    fitres2[7], fitres2_err[7]));
        pv2->SetFillColor(kWhite);
        pv2->SetBorderSize(1);
        pv2->SetTextFont(42);
        pv2->SetTextSize(0.037);
        pv2->SetTextAlign(12); // text left aligned 
        pv2->Draw();

        
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


        // write fit plot frame in a file 
        fOutFitResults_plotFrame->cd();
        frame->Write();

      }
    }
  }

  
  // write fit results histos in a file 
  TFile* fOutFitResults = new TFile(string(Form("file_FitResults_%s.root",sampletype_name.c_str())).c_str(),"recreate");
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
  
  fOutFitResults_plotFrame->Close();
  delete fOutFitResults_plotFrame;


}// end doTheFit DCB with prefit function 







// take fit results and compute dilepton scale 
void computeDileptonScale(string outputPathDileptonScalePlots, string lumiText, string sampletype_name)
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
  TFile* fInFitResults2 = TFile::Open(string(Form("file_FitResults_%s.root",sampletype_name.c_str())).c_str());

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
void compareDataMCfitPlots(string outputPathCompare2lDataMcFit, string lumiText, string sampletype_name)
{

  //***********************************************
  // *** read file with rooplots
  TFile* fIn_nominal = TFile::Open("file_FitResultsPlots_nominal.root");

  // define input rooplots 
  RooPlot* inputFrame_nominal[nDatasets][nCatEta][nCatpT];   

  // read rooplots from file 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        inputFrame_nominal[dat][catEta][catPt] = (RooPlot*)fIn_nominal->Get(Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));

      }
    }
  }

  //***********************************************
  // *** read file with rooplots
  TFile* fIn_var = TFile::Open(string(Form("file_FitResultsPlots_%s.root",sampletype_name.c_str())).c_str());

  // define input rooplots 
  RooPlot* inputFrame_var[nDatasets][nCatEta][nCatpT];   

  // read rooplots from file 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        inputFrame_var[dat][catEta][catPt] = (RooPlot*)fIn_var->Get(Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
       
      }
    }
  }


  //***********************************************
  // *** read file with fit results
  TFile* fInFitRes_nominal = TFile::Open("file_FitResults_nominal.root");

  // define input histos 
  TH1F* hin_nom_poleBW[nDatasets][nCatEta][nCatpT];
  TH1F* hin_nom_widthBW[nDatasets][nCatEta][nCatpT];  
  TH1F* hin_nom_meanDCB[nDatasets][nCatEta][nCatpT];
  TH1F* hin_nom_sigmaDCB[nDatasets][nCatEta][nCatpT]; 
  TH1F* hin_nom_a1DCB[nDatasets][nCatEta][nCatpT]; 
  TH1F* hin_nom_n1DCB[nDatasets][nCatEta][nCatpT];  
  TH1F* hin_nom_a2DCB[nDatasets][nCatEta][nCatpT]; 
  TH1F* hin_nom_n2DCB[nDatasets][nCatEta][nCatpT]; 
  
  // read histos with fit results from file 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        hin_nom_poleBW[dat][catEta][catPt] = (TH1F*)fInFitRes_nominal->Get(Form("hfitResults_poleBW_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_nom_widthBW[dat][catEta][catPt] = (TH1F*)fInFitRes_nominal->Get(Form("hfitResults_widthBW_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_nom_meanDCB[dat][catEta][catPt] = (TH1F*)fInFitRes_nominal->Get(Form("hfitResults_meanDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_nom_sigmaDCB[dat][catEta][catPt] = (TH1F*)fInFitRes_nominal->Get(Form("hfitResults_sigmaDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_nom_a1DCB[dat][catEta][catPt] = (TH1F*)fInFitRes_nominal->Get(Form("hfitResults_a1DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_nom_n1DCB[dat][catEta][catPt] = (TH1F*)fInFitRes_nominal->Get(Form("hfitResults_n1DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_nom_a2DCB[dat][catEta][catPt] = (TH1F*)fInFitRes_nominal->Get(Form("hfitResults_a2DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_nom_n2DCB[dat][catEta][catPt] = (TH1F*)fInFitRes_nominal->Get(Form("hfitResults_n2DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        
      }
    }
  }

  //***********************************************
  // *** read file with fit results
  TFile* fInFitRes_var = TFile::Open(string(Form("file_FitResults_%s.root",sampletype_name.c_str())).c_str());

  // define input histos 
  TH1F* hin_var_poleBW[nDatasets][nCatEta][nCatpT];
  TH1F* hin_var_widthBW[nDatasets][nCatEta][nCatpT];  
  TH1F* hin_var_meanDCB[nDatasets][nCatEta][nCatpT];
  TH1F* hin_var_sigmaDCB[nDatasets][nCatEta][nCatpT]; 
  TH1F* hin_var_a1DCB[nDatasets][nCatEta][nCatpT]; 
  TH1F* hin_var_n1DCB[nDatasets][nCatEta][nCatpT];  
  TH1F* hin_var_a2DCB[nDatasets][nCatEta][nCatpT]; 
  TH1F* hin_var_n2DCB[nDatasets][nCatEta][nCatpT]; 
  
  // read histos with fit results from file 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        hin_var_poleBW[dat][catEta][catPt] = (TH1F*)fInFitRes_var->Get(Form("hfitResults_poleBW_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_var_widthBW[dat][catEta][catPt] = (TH1F*)fInFitRes_var->Get(Form("hfitResults_widthBW_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_var_meanDCB[dat][catEta][catPt] = (TH1F*)fInFitRes_var->Get(Form("hfitResults_meanDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_var_sigmaDCB[dat][catEta][catPt] = (TH1F*)fInFitRes_var->Get(Form("hfitResults_sigmaDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_var_a1DCB[dat][catEta][catPt] = (TH1F*)fInFitRes_var->Get(Form("hfitResults_a1DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_var_n1DCB[dat][catEta][catPt] = (TH1F*)fInFitRes_var->Get(Form("hfitResults_n1DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_var_a2DCB[dat][catEta][catPt] = (TH1F*)fInFitRes_var->Get(Form("hfitResults_a2DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_var_n2DCB[dat][catEta][catPt] = (TH1F*)fInFitRes_var->Get(Form("hfitResults_n2DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        
      }
    }
  }

  
  // *************************************  
  //plot
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        TCanvas *can = new TCanvas("can","can");
        can->cd();

        inputFrame_nominal[dat][catEta][catPt]->Draw();
        inputFrame_var[dat][catEta][catPt]->Draw("same");

        
        // nominal fit results
        TPaveText* pv1 = new TPaveText(0.10,0.51,0.39,0.88,"brNDC");
        pv1->AddText("nominal ");
        pv1->AddText(Form("BW pole: %.3f #pm %.3f",   hin_nom_poleBW[1][catEta][catPt]->GetBinContent(1),   hin_nom_poleBW[1][catEta][catPt]->GetBinError(1)));
        pv1->AddText(Form("BW width: %.3f #pm %.3f",  hin_nom_widthBW[1][catEta][catPt]->GetBinContent(1),  hin_nom_widthBW[1][catEta][catPt]->GetBinError(1)));
        pv1->AddText(Form("DCB mean: %.3f #pm %.3f",  hin_nom_meanDCB[1][catEta][catPt]->GetBinContent(1),  hin_nom_meanDCB[1][catEta][catPt]->GetBinError(1)));
        pv1->AddText(Form("DCB sigma: %.3f #pm %.3f", hin_nom_sigmaDCB[1][catEta][catPt]->GetBinContent(1), hin_nom_sigmaDCB[1][catEta][catPt]->GetBinError(1)));
        pv1->AddText(Form("DCB a1: %.3f #pm %.3f",    hin_nom_a1DCB[1][catEta][catPt]->GetBinContent(1),    hin_nom_a1DCB[1][catEta][catPt]->GetBinError(1)));
        pv1->AddText(Form("DCB n1: %.3f #pm %.3f",    hin_nom_n1DCB[1][catEta][catPt]->GetBinContent(1),    hin_nom_n1DCB[1][catEta][catPt]->GetBinError(1)));
        pv1->AddText(Form("DCB a2: %.3f #pm %.3f",    hin_nom_a2DCB[1][catEta][catPt]->GetBinContent(1),    hin_nom_a2DCB[1][catEta][catPt]->GetBinError(1)));
        pv1->AddText(Form("DCB n2: %.3f #pm %.3f",    hin_nom_n2DCB[1][catEta][catPt]->GetBinContent(1),    hin_nom_n2DCB[1][catEta][catPt]->GetBinError(1)));
        pv1->SetFillColor(kWhite);
        pv1->SetBorderSize(1);
        pv1->SetTextColor(kBlue);
        pv1->SetTextFont(42);
        pv1->SetTextSize(0.037);
        pv1->SetTextAlign(12); // text left aligned 
        pv1->Draw();

        // data fit results
        TPaveText* pv2 = new TPaveText(0.66,0.51,0.95,0.88,"brNDC");
        pv2->AddText(sampletype_name.c_str());
        pv2->AddText(Form("BW pole: %.3f #pm %.3f",   hin_var_poleBW[0][catEta][catPt]->GetBinContent(1),   hin_var_poleBW[0][catEta][catPt]->GetBinError(1)));
        pv2->AddText(Form("BW width: %.3f #pm %.3f",  hin_var_widthBW[0][catEta][catPt]->GetBinContent(1),  hin_var_widthBW[0][catEta][catPt]->GetBinError(1)));
        pv2->AddText(Form("DCB mean: %.3f #pm %.3f",  hin_var_meanDCB[0][catEta][catPt]->GetBinContent(1),  hin_var_meanDCB[0][catEta][catPt]->GetBinError(1)));
        pv2->AddText(Form("DCB sigma: %.3f #pm %.3f", hin_var_sigmaDCB[0][catEta][catPt]->GetBinContent(1), hin_var_sigmaDCB[0][catEta][catPt]->GetBinError(1)));
        pv2->AddText(Form("DCB a1: %.3f #pm %.3f",    hin_var_a1DCB[0][catEta][catPt]->GetBinContent(1),    hin_var_a1DCB[0][catEta][catPt]->GetBinError(1)));
        pv2->AddText(Form("DCB n1: %.3f #pm %.3f",    hin_var_n1DCB[0][catEta][catPt]->GetBinContent(1),    hin_var_n1DCB[0][catEta][catPt]->GetBinError(1)));
        pv2->AddText(Form("DCB a2: %.3f #pm %.3f",    hin_var_a2DCB[0][catEta][catPt]->GetBinContent(1),    hin_var_a2DCB[0][catEta][catPt]->GetBinError(1)));
        pv2->AddText(Form("DCB n2: %.3f #pm %.3f",    hin_var_n2DCB[0][catEta][catPt]->GetBinContent(1),    hin_var_n2DCB[0][catEta][catPt]->GetBinError(1)));
        pv2->SetFillColor(kWhite);
        pv2->SetBorderSize(1);
        pv2->SetTextColor(kRed);
        pv2->SetTextFont(42);
        pv2->SetTextSize(0.037);
        pv2->SetTextAlign(12); // text left aligned 
        pv2->Draw();


        
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

  // *** bool for executing functions
  bool REDO2lHISTOS;
  bool REDOTHE2lFIT;
  bool COMPUTE2lSCALE;
  bool COMPARE2lDATAMCFIT;

  string sampletype_name;

  if(sampletype_varUnif001){

    REDO2lHISTOS = 1;
    REDOTHE2lFIT = 1;
    COMPUTE2lSCALE = 1;
    COMPARE2lDATAMCFIT = 1;
    sampletype_name = "varUnif001";

  } else if(sampletype_varStair){

    REDO2lHISTOS = 1;
    REDOTHE2lFIT = 1;
    COMPUTE2lSCALE = 1;
    COMPARE2lDATAMCFIT = 1;
    sampletype_name = "varStair";

  } else{

    REDO2lHISTOS = 1;
    REDOTHE2lFIT = 1;
    COMPUTE2lSCALE = 0;
    COMPARE2lDATAMCFIT = 0;
    sampletype_name = "nominal";

  }
  // ***

 
  // *** choose data/mc input path 
  string inputPathMC_DY = "/data3/Higgs/180416/MC_main/";
  string inputPathData = "/data3/Higgs/180416/";
  string inputPathMC_ggH = "/data3/Higgs/180416/MC_main/";

  string outputPathFitResultsPlots = string(Form("plotsSysts_FitResults_%s",sampletype_name.c_str())).c_str();
  string outputPathDileptonScalePlots = string(Form("plotsSysts_DileptonScale_%s",sampletype_name.c_str())).c_str();
  string outputPathCompare2lDataMcFit = string(Form("plotsSysts_CompareDataMC2lFit_%s",sampletype_name.c_str())).c_str();
  
  
  // *** define luminosity
  float lumi = 41.30; //fb-1
  string lumiText = "41.30 fb^{-1}";


  // *** create output directories
  if(REDOTHE2lFIT) gSystem->Exec(("mkdir -p "+outputPathFitResultsPlots).c_str());  //dir for fit results plots

  if(COMPUTE2lSCALE) gSystem->Exec(("mkdir -p "+outputPathDileptonScalePlots).c_str()); //dir for dilepton scale plots

  if(COMPARE2lDATAMCFIT) gSystem->Exec(("mkdir -p "+outputPathCompare2lDataMcFit).c_str()); //dir for comparison between Data and MC 2l fit

 


  // *** execute functions 

  
  if(REDO2lHISTOS) do2lHistograms_AN(inputPathMC_DY, inputPathData, lumi, sampletype_name); //read MC file and do histograms (events separated into categories 
                                                                           //based on pT and Eta of 1 of the 2 leptons, determined randomly, and integrating
                                                                           //over the other) as explained in AN2016_442 (Section 9)
                                            
  if(REDOTHE2lFIT) doThe2lFit_DCBfit(outputPathFitResultsPlots, lumiText, sampletype_name); // do fit with DCBxBW 

  if(COMPUTE2lSCALE) computeDileptonScale(outputPathDileptonScalePlots, lumiText, sampletype_name);

  if(COMPARE2lDATAMCFIT) compareDataMCfitPlots(outputPathCompare2lDataMcFit, lumiText, sampletype_name);

  

  
}



