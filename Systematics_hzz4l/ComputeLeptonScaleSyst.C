// **********************************
// 
// usage:
//     - specify parameters (input/output directory, luminosity) at the end of this file
//     - run with:
//                 root -l -b -q ComputeLeptonScaleSyst.C++
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

#include "Plotter/CMS_lumi.C" // CMS official label definition
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"  // contains DCB definition


using namespace std;
using namespace RooFit ;

#define REDO2lHISTOS 1
#define REDOTHE2lFIT 1
#define COMPARE2lDATAMCFIT 1

#define WRITEEXTRATEXTONPLOTS 1 // draw Preliminary on Plots


// *** global definitions
Float_t Zmass_nominalPDG  = 91.1876; // GeV/c^2 (http://pdg.lbl.gov/2017/listings/rpp2017-list-z-boson.pdf)
Float_t Zwidth_nominalPDG = 2.4952;  // GeV/c^2 (http://pdg.lbl.gov/2017/listings/rpp2017-list-z-boson.pdf)


const int nDatasets = 2;
string datasets[nDatasets] = {
    "AllData",
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
  "mu_Eta_1.4-2.5"
};

// categories of lepton pT 
const int nCatpT = 6;
enum CategpT {pTmin20 = 0, pT2030 = 1, pT3040 = 2, pT4050 = 3, pT5060 = 4, pT60100 = 5};
string sCategpT[nCatpT] = {
  "pT_min-20",
  "pT_20-30",
  "pT_30-40",
  "pT_40_50",
  "pT_50-60",
  "pT_60-100"
};

// processes
enum Process {Data=0, DY=1};



// *** read file and do histograms
void do2lHistograms(string inputPathMC, string inputPathData, float lumi)
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
  

  int currentProcess = -1;
  int currentCategEta = -1;
  int currentCategPt = -1;


  // loop over datasets 
  for(int d=0; d<nDatasets; d++){

    // assign process to dataset
    if(datasets[d]=="AllData") currentProcess = Data;
    if(datasets[d]=="DYJetsToLL_M50") currentProcess = DY;


    string inputFileName = string(Form("%s%s/ZZ4lAnalysis.root",(currentProcess==Data?inputPathData:inputPathMC).c_str(),datasets[d].c_str()));
    inputFile[d] = TFile::Open(inputFileName.c_str());


    hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");    
    gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(40);
    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d];


    inputTree[d] = (TTree*)inputFile[d]->Get("ZTree/candTree");
    inputTree[d]->SetBranchAddress("RunNumber", &nRun);
    inputTree[d]->SetBranchAddress("EventNumber", &nEvent);
    inputTree[d]->SetBranchAddress("LumiNumber", &nLumi);
    inputTree[d]->SetBranchAddress("Zsel", &Zsel);      // WARNING: this script works with the Z inclusive tree 
    inputTree[d]->SetBranchAddress("ZMass", &ZMass);    // WARNING: this script works with the Z inclusive tree
    inputTree[d]->SetBranchAddress("LepPt", &LepPt);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);
    inputTree[d]->SetBranchAddress("LepLepId", &LepLepId);
    if(currentProcess!=Data){
      inputTree[d]->SetBranchAddress("xsec", &xsec); 
      inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    } 

    
    // process tree 
    Long64_t entries = inputTree[d]->GetEntries();
    cout<<"Processing dataset "<<datasets[d]<<" ("<<entries<<" entries) ..."<<endl;

    for (Long64_t z=0; z<entries; ++z){

      inputTree[d]->GetEntry(z);

      if( Zsel < 0. ) continue;  // skip events that do not pass the trigger
      if(ZMass < 40.) continue;  // applying Z1 mass request 


      // choose leading lepton to cut on lepton pT (>20 GeV) as in ZZ
      // and subleading lepton to divide events in pT categories 
      Float_t LepPtLeading = std::max(LepPt->at(0),LepPt->at(1));
      if(LepPtLeading < 20.) continue;  
      Float_t LepPtsubLeading = std::min(LepPt->at(0),LepPt->at(1));
      

      // define event weight 
      Double_t eventWeight;
      if(currentProcess==Data) eventWeight = 1.;
      else eventWeight = partialSampleWeight[d] * xsec * overallEventWeight ;

      
      // Eta categories 
      //Z->ee histos 
      if(int(fabs(LepLepId->at(0))) == 11 ){
      
        if(fabs(LepEta->at(0)) >= 0. && fabs(LepEta->at(0)) < 0.8 ) currentCategEta = eleEta1st;
        else if(fabs(LepEta->at(0)) >= 0.8 && fabs(LepEta->at(0)) < 1.5 ) currentCategEta = eleEta2nd;
        else if(fabs(LepEta->at(0)) >= 1.5 && fabs(LepEta->at(0)) <= 2.5 ) currentCategEta = eleEta3rd;
        else cerr<<"error: wrong eta!"<<endl;
      }

      //Z->ee histos 
      if(int(fabs(LepLepId->at(0))) == 13 ){
      
        if(fabs(LepEta->at(0)) >= 0. && fabs(LepEta->at(0)) < 0.9 ) currentCategEta = muEta1st;
        else if(fabs(LepEta->at(0)) >= 0.9 && fabs(LepEta->at(0)) < 1.4 ) currentCategEta = muEta2nd;
        else if(fabs(LepEta->at(0)) >= 1.4 && fabs(LepEta->at(0)) <= 2.5 ) currentCategEta = muEta3rd;
        else cerr<<"error: wrong eta!"<<endl;
      }

      
      // pT categories 
      if(LepPtsubLeading < 20. ) currentCategPt = pTmin20;
      else if(LepPtsubLeading >= 20. && LepPtsubLeading < 30. ) currentCategPt = pT2030;
      else if(LepPtsubLeading >= 30. && LepPtsubLeading < 40. ) currentCategPt = pT3040;
      else if(LepPtsubLeading >= 40. && LepPtsubLeading < 50. ) currentCategPt = pT4050;
      else if(LepPtsubLeading >= 50. && LepPtsubLeading < 60. ) currentCategPt = pT5060;
      else if(LepPtsubLeading >= 60. && LepPtsubLeading <= 100. ) currentCategPt = pT60100;
      else continue;

      
      
      if(currentCategEta < 0 || currentCategPt < 0) continue;
      
      // fill histos
      hist[currentProcess][currentCategEta][currentCategPt]->Fill(ZMass,eventWeight); 


    } //end loop over tree entries 


  }//end loop over datasets


  // write histos in a file 
  TFile* fOutHistos = new TFile("file_DataMCHistos.root","recreate");
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




// perform the fit 
void doThe2lFit(string outputPathFitResultsPlots, string lumiText)
{

  // define output file for fit result plots 
  TFile* fOutFitResults_plotFrame = new TFile("file_FitResultsPlots.root","recreate");

  // read file with histos
  TFile* fInHistos = TFile::Open("file_DataMCHistos.root");

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
        
        // hfitResults_poleBW[dat][catEta][catPt]->Sumw2(true);
        // hfitResults_widthBW[dat][catEta][catPt]->Sumw2(true);
        // hfitResults_meanDCB[dat][catEta][catPt]->Sumw2(true);
        // hfitResults_sigmaDCB[dat][catEta][catPt]->Sumw2(true);
        // hfitResults_a1DCB[dat][catEta][catPt]->Sumw2(true);
        // hfitResults_n1DCB[dat][catEta][catPt]->Sumw2(true);
        // hfitResults_a2DCB[dat][catEta][catPt]->Sumw2(true);
        // hfitResults_n2DCB[dat][catEta][catPt]->Sumw2(true);


        // define roofit variable for the fit
        RooRealVar mll = RooRealVar("mll","m_{l^{+}l^{-}}",60,120);

        // set fit range
        mll.setRange("range80100gev",80,100); 

        // make roofit datasets from root histos
        RooDataHist dh("dh","dh",mll,Import(*inputhist[dat][catEta][catPt]));

        // define fit function: 
        // convolution of a BW function and a DoubleCrystalBall resolution function 
        // ---Breit-Wigner function (physical mass distribution)
        RooRealVar pole_BW("pole_BW","pole_BW",Zmass_nominalPDG);
        RooRealVar width_BW("width_BW","width_BW",Zwidth_nominalPDG);

        pole_BW.setConstant(Zmass_nominalPDG);    //fix to physical value
        width_BW.setConstant(Zwidth_nominalPDG);  //fix to physical value

        RooBreitWigner BW_pdf("BW_pdf","Breit-Wigner function",mll,pole_BW,width_BW);

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
        
        
        // do the fit 
        tot_pdf.fitTo(dh,Range("range80100gev"));
                
        // plot data on the frame
        RooPlot* frame = mll.frame();
        frame->SetName(inputhist[dat][catEta][catPt]->GetName());  // name is the name which appears in the root file
        frame->SetTitle(inputhist[dat][catEta][catPt]->GetName()); // title is the title on the canvas  
        dh.plotOn(frame,DataError(RooAbsData::SumW2)); 
        tot_pdf.plotOn(frame,LineColor(dat==0?kBlue:kRed)); //blue line for data, red for MC 


        // plot on the canvas and save plots
        TCanvas* c = new TCanvas(inputhist[dat][catEta][catPt]->GetName(),inputhist[dat][catEta][catPt]->GetName()); 
        c->cd();
        frame->Draw();

        // draw fit results on canvas
        TPaveText* pv = new TPaveText(0.69,0.52,0.95,0.87,"brNDC");
        pv->AddText(("BW pole: " + std::to_string(pole_BW.getVal())).c_str());
        pv->AddText(("BW width: " + std::to_string(width_BW.getVal())).c_str());
        pv->AddText(("DCB mean: " + std::to_string(mean_DCB.getVal())).c_str());
        pv->AddText(("DCB sigma: " + std::to_string(sigma_DCB.getVal())).c_str());
        pv->AddText(("DCB a1: " + std::to_string(a1_DCB.getVal())).c_str());
        pv->AddText(("DCB n1: " + std::to_string(n1_DCB.getVal())).c_str());
        pv->AddText(("DCB a2: " + std::to_string(a2_DCB.getVal())).c_str());
        pv->AddText(("DCB n2: " + std::to_string(n2_DCB.getVal())).c_str());
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
  TFile* fOutFitResults = new TFile("file_FitResults.root","recreate");
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


}// end doTheFit function 



// take fit results and compute dilepton scale 
void computeDileptonScale(string outputPathDileptonScalePlots, string lumiText)
{
  
 // read file with fit results
  TFile* fInFitResults = TFile::Open("file_FitResults.root");

  // define input histos 
  TH1F* hinput_meanFitResults[nDatasets][nCatEta][nCatpT];

  // read histos with fit results from file 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        hinput_meanFitResults[dat][catEta][catPt] = (TH1F*)fInFitResults->Get(Form("hfitResults_meanDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        //cout<<hinput_meanFitResults[dat][catEta][catPt]->GetName()<<endl; // debug

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

      vec_dileptonScale[catEta][catPt] = (hinput_meanFitResults[0][catEta][catPt]->GetBinContent(1) - hinput_meanFitResults[1][catEta][catPt]->GetBinContent(1) )/ Zmass_nominalPDG;  // data_mean - MC_mean /ZmassPDG
      vec_dileptonScale_err[catEta][catPt] = TMath::Power((TMath::Power(hinput_meanFitResults[0][catEta][catPt]->GetBinError(1),2) + TMath::Power(hinput_meanFitResults[1][catEta][catPt]->GetBinError(1),2)), 0.5) / Zmass_nominalPDG;  


      h_dileptonScale[catEta][catPt] = new TH1F(Form("h_dileptonScale_%s_%s",sCategEta[catEta].c_str(),sCategpT[catPt].c_str()),Form("h_dileptonScale_%s_%s",sCategEta[catEta].c_str(),sCategpT[catPt].c_str()), 1,0,1);
      //h_dileptonScale[catEta][catPt]->Sumw2(true);        

      h_dileptonScale[catEta][catPt]->Fill(0.5,vec_dileptonScale[catEta][catPt]);
      h_dileptonScale[catEta][catPt]->SetBinError(1,vec_dileptonScale_err[catEta][catPt]);

    }
  }
  
  // --- dilepton scale plots 

  // define x axis values: mean value of each pT bin 
  Float_t vec_ptValues_ele[nCatpT] = {13.5, 25., 35., 45., 55., 80.}; // first bin for ele: 7-20 GeV
  Float_t vec_ptValues_mu[nCatpT] = {12.5, 25., 35., 45., 55., 80.};  // first bin for mu: 5-20 GeV
  // define x axis uncertainties
  Float_t vec_ptValues_ele_err[nCatpT] = {6.5, 5., 5., 5., 5., 20.};
  Float_t vec_ptValues_mu_err[nCatpT] = {7.5, 5., 5., 5., 5., 20.};

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

  mg_ele->SetMinimum(-0.02);
  mg_ele->SetMaximum(0.02);

  TCanvas* can_ele = new TCanvas("can_ele","can_ele",600,600);
  mg_ele->Draw("ap"); 
  mg_ele->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
  mg_ele->GetXaxis()->SetTitleFont(43);
  mg_ele->GetXaxis()->SetTitleSize(18);
  mg_ele->GetXaxis()->SetTitleOffset(1.3);
  mg_ele->GetXaxis()->SetLabelFont(43);
  mg_ele->GetXaxis()->SetLabelSize(13);
  mg_ele->GetXaxis()->SetLimits(0.,100.); //set x axis limits (same as SetRangeUser but for TGraphs)
  mg_ele->GetYaxis()->SetTitle("(m^{peak}_{data} - m^{peak}_{MC})/m_{PDG}");
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

  mg_mu->SetMinimum(-0.02);
  mg_mu->SetMaximum(0.02);
 
  TCanvas* can_mu = new TCanvas("can_mu","can_mu",600,600);
  mg_mu->Draw("ap");
  mg_mu->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
  mg_mu->GetXaxis()->SetTitleFont(43);
  mg_mu->GetXaxis()->SetTitleSize(18);
  mg_mu->GetXaxis()->SetTitleOffset(1.3);
  mg_mu->GetXaxis()->SetLabelFont(43);
  mg_mu->GetXaxis()->SetLabelSize(13);
  mg_mu->GetXaxis()->SetLimits(0.,100.); //set x axis limits (same as SetRangeUser but for TGraphs)
  mg_mu->GetYaxis()->SetTitle("(m^{peak}_{data} - m^{peak}_{MC})/m_{PDG}");
  mg_mu->GetYaxis()->SetTitleFont(43);
  mg_mu->GetYaxis()->SetTitleSize(18);
  mg_mu->GetYaxis()->SetTitleOffset(1.5);
  mg_mu->GetYaxis()->SetLabelFont(43);
  mg_mu->GetYaxis()->SetLabelSize(13);
  gPad->SetTickx(); // set ticks also on upper x axis
  gPad->SetTicky(); // set ticks also on right y axis

  TLegend* legend_mu = new TLegend(0.58,0.15,0.8,0.35);
  legend_mu->AddEntry(graph_eleCatEta0,"Z,|#eta| 0.0-0.9","lp");
  legend_mu->AddEntry(graph_eleCatEta1,"Z,|#eta| 0.9-1.5","lp");
  legend_mu->AddEntry(graph_eleCatEta2,"Z,|#eta| 1.5-2.5","lp");
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



  // --- save dilepton scale into a file 
  TFile* fOutDileptonScale = new TFile("file_DileptonScale.root","recreate");
  fOutDileptonScale->cd();
  for(int catEta=0; catEta<nCatEta; catEta++){
    for(int catPt=0; catPt<nCatpT; catPt++){

      h_dileptonScale[catEta][catPt]->Write(h_dileptonScale[catEta][catPt]->GetName());
      delete h_dileptonScale[catEta][catPt];

    }
  }
  fOutDileptonScale->Close();
  delete fOutDileptonScale;


}// end computeDileptonScale function



// comparison between Data and MC 2l fit
void compareDataMCfitPlots(string outputPathCompare2lDataMcFit, string lumiText)
{



}// end compareDataMCfitPlots function




// *** main function
void ComputeLeptonScaleSyst()
{
 
  string inputPathMC = "/data3/Higgs/180416/MC_main/";
  string inputPathData = "/data3/Higgs/180416/";

  string outputPathFitResultsPlots = "plotsSysts_FitResults";
  string outputPathDileptonScalePlots = "plotsSysts_DileptonScale";
  string outputPathCompare2lDataMcFit = "plotsSysts_CompareDataMC2lFit";

  float lumi = 41.30; //fb-1
  string lumiText = "41.30 fb^{-1}";


  // create output directories
  if(REDOTHE2lFIT) gSystem->Exec(("mkdir -p "+outputPathFitResultsPlots).c_str());  //dir for fit results plots

  gSystem->Exec(("mkdir -p "+outputPathDileptonScalePlots).c_str()); //dir for dilepton scale plots

  if(COMPARE2lDATAMCFIT) gSystem->Exec(("mkdir -p "+outputPathCompare2lDataMcFit).c_str()); //dir for comparison between Data and MC 2l fit


  // execute functions 
  if(REDO2lHISTOS) do2lHistograms(inputPathMC, inputPathData, lumi);

  if(REDOTHE2lFIT) doThe2lFit(outputPathFitResultsPlots, lumiText);

  computeDileptonScale(outputPathDileptonScalePlots, lumiText);

  if(COMPARE2lDATAMCFIT) compareDataMCfitPlots(outputPathCompare2lDataMcFit, lumiText);

  

}



