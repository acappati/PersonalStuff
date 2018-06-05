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
#include "TLorentzVector.h"

#include "Plotter/CMS_lumi.C" // CMS official label definition
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"  // contains DCB definition


using namespace std;
using namespace RooFit ;

#define REDO2lHISTOS 1
#define REDOTHE2lFIT 1
#define COMPUTE2lSCALE 1
#define COMPARE2lDATAMCFIT 1
#define REDO4lHISTOS 1
#define COMPUTE4lSCALE 1

#define WRITEEXTRATEXTONPLOTS 1 // draw Preliminary on Plots


// *** global definitions
Float_t Zmass_nominalPDG  = 91.1876; // GeV/c^2 (http://pdg.lbl.gov/2017/listings/rpp2017-list-z-boson.pdf)
Float_t Zwidth_nominalPDG = 2.4952;  // GeV/c^2 (http://pdg.lbl.gov/2017/listings/rpp2017-list-z-boson.pdf)
Float_t mass_ele_nominalPDG = 0.0005; //GeV/c^2   
Float_t mass_mu_nominalPDG  = 0.1057; //GeV/c^2


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
enum Process {Data=0, DY=1};

//final state
const int nFinalStates = 3;
enum FinalState {fs4e = 0, fs4mu = 1, fs2e2mu = 2};
string sFinalState[nFinalStates] = {
  "fs_4e",
  "fs_4mu",
  "fs_2e2mu"
};


//variations of the 4l distribution (per category variation)
const int nVariations4lDistr = 2;
enum variations4lDistr {distrNominal = 0, distrVar = 1};
string sVariations4lDistr[nVariations4lDistr] = {
  "distr_nominal",
  "distr_var"
};

//variations of the 4l distribution (max variation)  
const int nVariations4lDistr_maxVar = 3;
enum variations4lDistr_maxVar {nominal = 0, upVar = 1, dnVar = 2};
string sVariations4lDistr_maxVar[nVariations4lDistr_maxVar] = {
  "nominal",
  "upVar",
  "dnVar"
};



//************************************************

// *** read file and do histograms
void do2lHistograms(string inputPathMC_DY, string inputPathData, float lumi)
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


    string inputFileName = string(Form("%s%s/ZZ4lAnalysis.root",(currentProcess==Data?inputPathData:inputPathMC_DY).c_str(),datasets[d].c_str()));
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

      //Z->mumu histos 
      if(int(fabs(LepLepId->at(0))) == 13 ){
      
        if(fabs(LepEta->at(0)) >= 0. && fabs(LepEta->at(0)) < 0.9 ) currentCategEta = muEta1st;
        else if(fabs(LepEta->at(0)) >= 0.9 && fabs(LepEta->at(0)) < 1.4 ) currentCategEta = muEta2nd;
        else if(fabs(LepEta->at(0)) >= 1.4 && fabs(LepEta->at(0)) <= 2.4 ) currentCategEta = muEta3rd;
        else cerr<<"error: wrong eta!"<<endl;
      }

      
      // pT categories 
      if(LepPt->at(0) < 20. ) currentCategPt = pTmin20;
      else if(LepPt->at(0) >= 20. && LepPt->at(0) < 30. ) currentCategPt = pT2030;
      else if(LepPt->at(0) >= 30. && LepPt->at(0) < 40. ) currentCategPt = pT3040;
      else if(LepPt->at(0) >= 40. && LepPt->at(0) < 50. ) currentCategPt = pT4050;
      else if(LepPt->at(0) >= 50. && LepPt->at(0) <= 100. ) currentCategPt = pT50100;
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

        RooBreitWigner BW_pdf("BW_pdf","Breit-Wigner function",mll,pole_BW,width_BW);

        pole_BW.setConstant(1);   //fix to physical value (setConstant take bool as argument, default is true)
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
        
        
        // *** do the 1st fit 
        tot_pdf.chi2FitTo(dh,Range("range80100gev")); //a binned ML fit (done with fitTo) always assume Poisson error interpretation of data data 
                                                      //with non-unit wieghts can only be correctly fitted with a Chi2 fit (chi2FitTo) 

        double fitres1[8] = {pole_BW.getVal(), width_BW.getVal(), mean_DCB.getVal(), sigma_DCB.getVal(), a1_DCB.getVal(), n1_DCB.getVal(), a2_DCB.getVal(), n2_DCB.getVal()};
        double fitres1_err[8] = {pole_BW.getError(), width_BW.getError(), mean_DCB.getError(), sigma_DCB.getError(), a1_DCB.getError(), n1_DCB.getError(), a2_DCB.getError(), n2_DCB.getError()};


        // *** do the 2nd fit
        tot_pdf.chi2FitTo(dh,Range("range80100gev")); //a binned ML fit (done with fitTo) always assume Poisson error interpretation of data data 
                                                      //with non-unit wieghts can only be correctly fitted with a Chi2 fit (chi2FitTo) 

        double fitres2[8] = {pole_BW.getVal(), width_BW.getVal(), mean_DCB.getVal(), sigma_DCB.getVal(), a1_DCB.getVal(), n1_DCB.getVal(), a2_DCB.getVal(), n2_DCB.getVal()};
        double fitres2_err[8] = {pole_BW.getError(), width_BW.getError(), mean_DCB.getError(), sigma_DCB.getError(), a1_DCB.getError(), n1_DCB.getError(), a2_DCB.getError(), n2_DCB.getError()};


                
        // *** plot data on the frame
        RooPlot* frame = mll.frame();
        frame->SetName(inputhist[dat][catEta][catPt]->GetName());  // name is the name which appears in the root file
        frame->SetTitle(""); // title is the title on the canvas  
        dh.plotOn(frame,DataError(RooAbsData::SumW2), MarkerStyle(dat==0?kOpenCircle:kFullCircle)); 
        tot_pdf.plotOn(frame, NormRange("range80100gev"), LineColor(dat==0?kBlue:kRed)); //blue line for data, red for MC 
                                                                                         //NormRange needed to normalize pdf to data in the fitting range


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
        hfitResults_poleBW[dat][catEta][catPt]->Fill(0.5,fitres2[0]);
        hfitResults_poleBW[dat][catEta][catPt]->SetBinError(1,fitres2_err[0]);
        hfitResults_widthBW[dat][catEta][catPt]->Fill(0.5,fitres2[1]);  
        hfitResults_widthBW[dat][catEta][catPt]->SetBinError(1,fitres2_err[1]);
	                                                  
	hfitResults_meanDCB[dat][catEta][catPt]->Fill(0.5,fitres2[2]);  
        hfitResults_meanDCB[dat][catEta][catPt]->SetBinError(1,fitres2_err[2]);
	hfitResults_sigmaDCB[dat][catEta][catPt]->Fill(0.5,fitres2[3]); 
        hfitResults_sigmaDCB[dat][catEta][catPt]->SetBinError(1,fitres2_err[3]);
	hfitResults_a1DCB[dat][catEta][catPt]->Fill(0.5,fitres2[4]);
        hfitResults_a1DCB[dat][catEta][catPt]->SetBinError(1,fitres2_err[4]);
	hfitResults_n1DCB[dat][catEta][catPt]->Fill(0.5,fitres2[5]);
        hfitResults_n1DCB[dat][catEta][catPt]->SetBinError(1,fitres2_err[5]);
	hfitResults_a2DCB[dat][catEta][catPt]->Fill(0.5,fitres2[6]);
        hfitResults_a2DCB[dat][catEta][catPt]->SetBinError(1,fitres2_err[6]);
	hfitResults_n2DCB[dat][catEta][catPt]->Fill(0.5,fitres2[7]);    
        hfitResults_n2DCB[dat][catEta][catPt]->SetBinError(1,fitres2_err[7]);    


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

  //mg_ele->SetMinimum(-0.01);
  //mg_ele->SetMaximum(0.01);

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

  //mg_mu->SetMinimum(-0.005);
  //mg_mu->SetMaximum(0.005);
 
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

  // read file with rooplots 
  TFile* fIn_frames = TFile::Open("file_FitResultsPlots.root");

  // define input rooplots
  RooPlot* inframe[nDatasets][nCatEta][nCatpT];  

  // read rooplot from file
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        inframe[dat][catEta][catPt] = (RooPlot*)fIn_frames->Get(Form("hist_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
       
      }
    }
  }

  
  // read file with fit results
  TFile* fInFitResults = TFile::Open("file_FitResults.root");

  // define input histos 
  TH1F* hin_poleBW[nDatasets][nCatEta][nCatpT];  
  TH1F* hin_widthBW[nDatasets][nCatEta][nCatpT];  
  TH1F* hin_meanDCB[nDatasets][nCatEta][nCatpT];
  TH1F* hin_sigmaDCB[nDatasets][nCatEta][nCatpT]; 
  TH1F* hin_a1DCB[nDatasets][nCatEta][nCatpT]; 
  TH1F* hin_n1DCB[nDatasets][nCatEta][nCatpT];  
  TH1F* hin_a2DCB[nDatasets][nCatEta][nCatpT]; 
  TH1F* hin_n2DCB[nDatasets][nCatEta][nCatpT]; 
  
  // read histos with fit results from file 
  for(int dat=0; dat<nDatasets; dat++){
    for(int catEta=0; catEta<nCatEta; catEta++){
      for(int catPt=0; catPt<nCatpT; catPt++){

        hin_poleBW[dat][catEta][catPt] = (TH1F*)fInFitResults->Get(Form("hfitResults_poleBW_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_widthBW[dat][catEta][catPt] = (TH1F*)fInFitResults->Get(Form("hfitResults_widthBW_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_meanDCB[dat][catEta][catPt] = (TH1F*)fInFitResults->Get(Form("hfitResults_meanDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_sigmaDCB[dat][catEta][catPt] = (TH1F*)fInFitResults->Get(Form("hfitResults_sigmaDCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_a1DCB[dat][catEta][catPt] = (TH1F*)fInFitResults->Get(Form("hfitResults_a1DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_n1DCB[dat][catEta][catPt] = (TH1F*)fInFitResults->Get(Form("hfitResults_n1DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_a2DCB[dat][catEta][catPt] = (TH1F*)fInFitResults->Get(Form("hfitResults_a2DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        hin_n2DCB[dat][catEta][catPt] = (TH1F*)fInFitResults->Get(Form("hfitResults_n2DCB_%s_%s_%s",datasets[dat].c_str(),sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
        
      }
    }
  }
  

  // plot
  for(int catEta=0; catEta<nCatEta; catEta++){
    for(int catPt=0; catPt<nCatpT; catPt++){

      TCanvas *can = new TCanvas("can","can");
      can->cd();

      //fits
      inframe[1][catEta][catPt]->Draw(); // MC DY
      inframe[0][catEta][catPt]->Draw("same"); // data

      // MC fit results
      TPaveText* pv1 = new TPaveText(0.10,0.51,0.39,0.88,"brNDC");
      pv1->AddText("MC DY: ");
      pv1->AddText(Form("BW pole: %.3f #pm %.3f",   hin_poleBW[1][catEta][catPt]->GetBinContent(1),   hin_poleBW[1][catEta][catPt]->GetBinError(1)));
      pv1->AddText(Form("BW width: %.3f #pm %.3f",  hin_widthBW[1][catEta][catPt]->GetBinContent(1),  hin_widthBW[1][catEta][catPt]->GetBinError(1)));
      pv1->AddText(Form("DCB mean: %.3f #pm %.3f",  hin_meanDCB[1][catEta][catPt]->GetBinContent(1),  hin_meanDCB[1][catEta][catPt]->GetBinError(1)));
      pv1->AddText(Form("DCB sigma: %.3f #pm %.3f", hin_sigmaDCB[1][catEta][catPt]->GetBinContent(1), hin_sigmaDCB[1][catEta][catPt]->GetBinError(1)));
      pv1->AddText(Form("DCB a1: %.3f #pm %.3f",    hin_a1DCB[1][catEta][catPt]->GetBinContent(1),    hin_a1DCB[1][catEta][catPt]->GetBinError(1)));
      pv1->AddText(Form("DCB n1: %.3f #pm %.3f",    hin_n1DCB[1][catEta][catPt]->GetBinContent(1),    hin_n1DCB[1][catEta][catPt]->GetBinError(1)));
      pv1->AddText(Form("DCB a2: %.3f #pm %.3f",    hin_a2DCB[1][catEta][catPt]->GetBinContent(1),    hin_a2DCB[1][catEta][catPt]->GetBinError(1)));
      pv1->AddText(Form("DCB n2: %.3f #pm %.3f",    hin_n2DCB[1][catEta][catPt]->GetBinContent(1),    hin_n2DCB[1][catEta][catPt]->GetBinError(1)));
      pv1->SetFillColor(kWhite);
      pv1->SetBorderSize(1);
      pv1->SetTextColor(kRed);
      pv1->SetTextFont(42);
      pv1->SetTextSize(0.037);
      pv1->SetTextAlign(12); // text left aligned 
      pv1->Draw();

      // data fit results
      TPaveText* pv2 = new TPaveText(0.66,0.51,0.95,0.88,"brNDC");
      pv2->AddText("Data: ");
      pv2->AddText(Form("BW pole: %.3f #pm %.3f",   hin_poleBW[0][catEta][catPt]->GetBinContent(1),   hin_poleBW[0][catEta][catPt]->GetBinError(1)));
      pv2->AddText(Form("BW width: %.3f #pm %.3f",  hin_widthBW[0][catEta][catPt]->GetBinContent(1),  hin_widthBW[0][catEta][catPt]->GetBinError(1)));
      pv2->AddText(Form("DCB mean: %.3f #pm %.3f",  hin_meanDCB[0][catEta][catPt]->GetBinContent(1),  hin_meanDCB[0][catEta][catPt]->GetBinError(1)));
      pv2->AddText(Form("DCB sigma: %.3f #pm %.3f", hin_sigmaDCB[0][catEta][catPt]->GetBinContent(1), hin_sigmaDCB[0][catEta][catPt]->GetBinError(1)));
      pv2->AddText(Form("DCB a1: %.3f #pm %.3f",    hin_a1DCB[0][catEta][catPt]->GetBinContent(1),    hin_a1DCB[0][catEta][catPt]->GetBinError(1)));
      pv2->AddText(Form("DCB n1: %.3f #pm %.3f",    hin_n1DCB[0][catEta][catPt]->GetBinContent(1),    hin_n1DCB[0][catEta][catPt]->GetBinError(1)));
      pv2->AddText(Form("DCB a2: %.3f #pm %.3f",    hin_a2DCB[0][catEta][catPt]->GetBinContent(1),    hin_a2DCB[0][catEta][catPt]->GetBinError(1)));
      pv2->AddText(Form("DCB n2: %.3f #pm %.3f",    hin_n2DCB[0][catEta][catPt]->GetBinContent(1),    hin_n2DCB[0][catEta][catPt]->GetBinError(1)));
      pv2->SetFillColor(kWhite);
      pv2->SetBorderSize(1);
      pv2->SetTextColor(kBlue);
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
      can->SaveAs((outputPathCompare2lDataMcFit + "/" + Form("hist_%s_%s",sCategEta[catEta].c_str(),sCategpT[catPt].c_str()) + ".pdf").c_str());
      can->SaveAs((outputPathCompare2lDataMcFit + "/" + Form("hist_%s_%s",sCategEta[catEta].c_str(),sCategpT[catPt].c_str()) + ".png").c_str());

    }
  } 


  
}// end compareDataMCfitPlots function





// *** compute invariant mass function
float compute4lInvMass(vector<Float_t> *LepPt, vector<Float_t> *LepEta, vector<Float_t> *LepPhi, float* vec_lepMass, float* vec_scale2l)
{
  
  TLorentzVector lep0;
  TLorentzVector lep1;
  TLorentzVector lep2;
  TLorentzVector lep3;

  lep0.SetPtEtaPhiM( LepPt->at(0) * (1. + vec_scale2l[0]), LepEta->at(0), LepPhi->at(0), vec_lepMass[0] );
  lep1.SetPtEtaPhiM( LepPt->at(1) * (1. + vec_scale2l[1]), LepEta->at(1), LepPhi->at(1), vec_lepMass[1] );
  lep2.SetPtEtaPhiM( LepPt->at(2) * (1. + vec_scale2l[2]), LepEta->at(2), LepPhi->at(2), vec_lepMass[2] );
  lep3.SetPtEtaPhiM( LepPt->at(3) * (1. + vec_scale2l[3]), LepEta->at(3), LepPhi->at(3), vec_lepMass[3] );

  return (lep0 + lep1 + lep2 + lep3).M();

} // end of compute invariant mass function





// do 4l ggH histos  
void do4lHistograms_perCatVariation(string inputPathMC_ggH, float lumi)
{

  TH1::SetDefaultSumw2(true);

  TFile* inputFile;
  TTree* inputTree;
  TH1F* hCounters;
  Double_t gen_sumWeights;
  Float_t partialSampleWeight;

  Int_t nRun;
  Long64_t nEvent;
  Int_t nLumi;

  Short_t ZZsel;
  Float_t ZZMass;
  Short_t Z1Flav;
  Short_t Z2Flav;
  vector<Float_t> *LepPt = 0;
  vector<Float_t> *LepEta = 0;
  vector<Float_t> *LepPhi = 0;
  vector<Float_t> *LepLepId = 0;
  Float_t xsec;
  Float_t overallEventWeight;

  
  // define 4l histos 
  TH1F* h4lDistrib[nVariations4lDistr][nFinalStates]; 
  for(int var=0; var<nVariations4lDistr; var++){
    for(int fs=0; fs<nFinalStates; fs++){

      h4lDistrib[var][fs] = new TH1F(Form("h4lDistrib_%s_%s",sVariations4lDistr[var].c_str(),sFinalState[fs].c_str()),Form("h4lDistrib_%s_%s",sVariations4lDistr[var].c_str(),sFinalState[fs].c_str()),70,105.,140.);
      h4lDistrib[var][fs]->Sumw2(true);
    }
  }


  // read file with dilepton scale
  TFile* fIn_dileptonScale = TFile::Open("file_DileptonScale.root");

  TH1F* hIn_2lscale[nCatEta][nCatpT];
  
  for(int catEta=0; catEta<nCatEta; catEta++){
    for(int catPt=0; catPt<nCatpT; catPt++){

      hIn_2lscale[catEta][catPt] = (TH1F*)fIn_dileptonScale->Get(Form("h_dileptonScale_%s_%s",sCategEta[catEta].c_str(),sCategpT[catPt].c_str()));
    }
  }

 
  int currentFinalState = -1;
  int currentCategEta = -1;
  int currentCategPt = -1;


  string dataset = "ggH125";
  string inputFileName = string(Form("%s%s/ZZ4lAnalysis.root",inputPathMC_ggH.c_str(),dataset.c_str()));
  inputFile = TFile::Open(inputFileName.c_str());

  hCounters = (TH1F*)inputFile->Get("ZZTree/Counters");
  gen_sumWeights = (Long64_t)hCounters->GetBinContent(40);
  partialSampleWeight = lumi * 1000 / gen_sumWeights;


  inputTree = (TTree*)inputFile->Get("ZZTree/candTree");
  inputTree->SetBranchAddress("RunNumber", &nRun);
  inputTree->SetBranchAddress("EventNumber", &nEvent);
  inputTree->SetBranchAddress("LumiNumber", &nLumi);
  inputTree->SetBranchAddress("ZZsel", &ZZsel);      
  inputTree->SetBranchAddress("ZZMass", &ZZMass);   
  inputTree->SetBranchAddress("Z1Flav", &Z1Flav);
  inputTree->SetBranchAddress("Z2Flav", &Z2Flav);
  inputTree->SetBranchAddress("LepPt", &LepPt);
  inputTree->SetBranchAddress("LepEta", &LepEta);
  inputTree->SetBranchAddress("LepPhi", &LepPhi);
  inputTree->SetBranchAddress("LepLepId", &LepLepId);
  inputTree->SetBranchAddress("xsec", &xsec); 
  inputTree->SetBranchAddress("overallEventWeight", &overallEventWeight);
  
  
  //process tree 
  Long64_t entries = inputTree->GetEntries();
  cout<<"Processing dataset "<<dataset<<" ("<<entries<<" entries) ..."<<endl;
  
  for(Long64_t z=0; z<entries; ++z){

    inputTree->GetEntry(z);


    if(LepEta->size()!=4){
        cout<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", stored "<<LepEta->size()<<" leptons instead of 4"<<endl;
        continue;
    }


    if( !(ZZsel>=90) ) continue;


    Double_t eventWeight = partialSampleWeight * xsec * overallEventWeight ;


    //*** find final state 
    if(Z1Flav==-121){
	if(Z2Flav==-121)
	  currentFinalState = fs4e;
	else if(Z2Flav==-169)
	  currentFinalState = fs2e2mu;
	else
	  cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
      }else if(Z1Flav==-169){
	if(Z2Flav==-121)
	  currentFinalState = fs2e2mu;
	else if(Z2Flav==-169)
	  currentFinalState = fs4mu;
	else
	  cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z2Flav="<<Z2Flav<<endl;
      }else{
	cerr<<"error in event "<<nRun<<":"<<nLumi<<":"<<nEvent<<", Z1Flav="<<Z1Flav<<endl;
    }
     
    
    if(currentFinalState < 0) continue;

    //*** classify 4 leptons and assign variations
    float vec_scale2l[4] = {0.,0.,0.,0.};
    float vec_lepMass[4];

    for(int l=0; l<4; l++){

      // eta categories and lep mass 
      if(int(fabs(LepLepId->at(l))) == 11 ){

        vec_lepMass[l] = mass_ele_nominalPDG;

        if(fabs(LepEta->at(l)) >= 0. && fabs(LepEta->at(l)) < 0.8 ) currentCategEta = eleEta1st;
        else if(fabs(LepEta->at(l)) >= 0.8 && fabs(LepEta->at(l)) < 1.5 ) currentCategEta = eleEta2nd;
        else if(fabs(LepEta->at(l)) >= 1.5 && fabs(LepEta->at(l)) <= 2.5 ) currentCategEta = eleEta3rd;
        else cerr<<"error: wrong eta!"<<endl;
      }

      if(int(fabs(LepLepId->at(l))) == 13 ){

        vec_lepMass[l] = mass_mu_nominalPDG;
      
        if(fabs(LepEta->at(l)) >= 0. && fabs(LepEta->at(l)) < 0.9 ) currentCategEta = muEta1st;
        else if(fabs(LepEta->at(l)) >= 0.9 && fabs(LepEta->at(l)) < 1.4 ) currentCategEta = muEta2nd;
        else if(fabs(LepEta->at(l)) >= 1.4 && fabs(LepEta->at(l)) <= 2.4 ) currentCategEta = muEta3rd;
        else cerr<<"error: wrong eta!"<<endl;  
      } 

      // pT categories 
      if(LepPt->at(l) < 20. ) currentCategPt = pTmin20;
      else if(LepPt->at(l) >= 20. && LepPt->at(l) < 30. ) currentCategPt = pT2030;
      else if(LepPt->at(l) >= 30. && LepPt->at(l) < 40. ) currentCategPt = pT3040;
      else if(LepPt->at(l) >= 40. && LepPt->at(l) < 50. ) currentCategPt = pT4050;
      else if(LepPt->at(l) >= 50. && LepPt->at(l) <= 100. ) currentCategPt = pT50100;


      if(currentCategEta < 0 || currentCategPt < 0) continue;
           
      vec_scale2l[l] = hIn_2lscale[currentCategEta][currentCategPt]->GetBinContent(1);
    }

       
    //define ZZmass with lep pT scale variations  
    float ZZMass_var = compute4lInvMass(LepPt, LepEta, LepPhi, vec_lepMass, vec_scale2l); //apply to lepton pT the 2lscale 
    
    cout<<ZZMass<<" "<<ZZMass_var<<endl;
      

    //*** fill 4l histos 
    h4lDistrib[distrNominal][currentFinalState]->Fill(ZZMass,eventWeight);
    h4lDistrib[distrVar][currentFinalState]->Fill(ZZMass_var,eventWeight);

  } //end loop over tree entries



  // write histos in a file 
  TFile* fOut4lhist = new TFile("file_MC4lHistos_perCatVariation.root","recreate");
  fOut4lhist->cd();
  for(int var=0; var<nVariations4lDistr; var++){
    for(int fs=0; fs<nFinalStates; fs++){
       
      h4lDistrib[var][fs]->Write(h4lDistrib[var][fs]->GetName());
      delete h4lDistrib[var][fs];
    }
  }
  fOut4lhist->Close();
  delete fOut4lhist;



}// end do4lHistograms function




void compute4lScale_perCatVariation(string outputPath4lScaleFitPlots_perCat, string lumiText)
{ 

  // read file with 4l histos 
  TFile* fin4lhist = TFile::Open("file_MC4lHistos_perCatVariation.root");

  // input histos
  TH1F* hin4l[nVariations4lDistr][nFinalStates]; 

  // vec for store fit results and computing 4l scale
  float vec_4lfitRes[nVariations4lDistr][nFinalStates];

  // vec with fit results 
  double fitres[nVariations4lDistr][nFinalStates][6];
  double fitres_err[nVariations4lDistr][nFinalStates][6];

  // m4l rooplots 
  RooPlot* frame4l[nVariations4lDistr][nFinalStates];
  

  for(int var=0; var<nVariations4lDistr; var++){
    for(int fs=0; fs<nFinalStates; fs++){

      hin4l[var][fs] = (TH1F*)fin4lhist->Get(Form("h4lDistrib_%s_%s",sVariations4lDistr[var].c_str(),sFinalState[fs].c_str()));

      
      //*** FIT ***
      RooRealVar m4l = RooRealVar("m4l","m4l",105,140);
      
      RooDataHist dhm4l("dhm4l","dhm4l",m4l,Import(*hin4l[var][fs]));      

      RooRealVar mean_4lDCB("mean_4lDCB","mean_4lDCB",125.,120.,130.);
      RooRealVar sigma_4lDCB("sigma_4lDCB","sigma_4lDCB",1.6,0.001,30.);
      RooRealVar a1_4lDCB("a1_4lDCB","a1_4lDCB",1.46,0.5,50.);   
      RooRealVar n1_4lDCB("n1_4lDCB","n1_4lDCB",1.92,0.,50.);
      RooRealVar a2_4lDCB("a2_4lDCB","a2_4lDCB",1.46,0.,50.);
      RooRealVar n2_4lDCB("n2_4lDCB","n2_4lDCB",20.,0.,50.);

      RooDoubleCB DCB4l_pdf("DCB4l_pdf","Double Crystal ball function",m4l,mean_4lDCB,sigma_4lDCB,a1_4lDCB,n1_4lDCB,a2_4lDCB,n2_4lDCB);

      m4l.setRange("range115130gev",115,130);
   
      // do the fit
      DCB4l_pdf.chi2FitTo(dhm4l, Range("range115130gev"));

      // store fir res 
      fitres[var][fs][0] = mean_4lDCB.getVal();
      fitres[var][fs][1] = sigma_4lDCB.getVal();
      fitres[var][fs][2] = a1_4lDCB.getVal();
      fitres[var][fs][3] = n1_4lDCB.getVal();
      fitres[var][fs][4] = a2_4lDCB.getVal();
      fitres[var][fs][5] = n2_4lDCB.getVal();

      fitres_err[var][fs][0] = mean_4lDCB.getError();
      fitres_err[var][fs][1] = sigma_4lDCB.getError();
      fitres_err[var][fs][2] = a1_4lDCB.getError();
      fitres_err[var][fs][3] = n1_4lDCB.getError();
      fitres_err[var][fs][4] = a2_4lDCB.getError();
      fitres_err[var][fs][5] = n2_4lDCB.getError();
      

      // plot on frame 
      frame4l[var][fs] = m4l.frame();
      frame4l[var][fs]->SetName(hin4l[var][fs]->GetName());
      frame4l[var][fs]->SetTitle("");
      dhm4l.plotOn(frame4l[var][fs],DataError(RooAbsData::SumW2), MarkerStyle(var==distrNominal?kOpenCircle:kFullCircle));
      DCB4l_pdf.plotOn(frame4l[var][fs], NormRange("range115130gev"), LineColor(var==distrNominal?kBlue:kRed));

      
      TCanvas* c = new TCanvas(hin4l[var][fs]->GetName(),hin4l[var][fs]->GetName()); 
      c->cd();
      frame4l[var][fs]->Draw();

      // draw fit results on canvas
      TPaveText* pv1 = new TPaveText(0.10,0.55,0.41,0.88,"brNDC");
      pv1->AddText(Form("DCB mean: %.3f #pm %.3f",  fitres[var][fs][0], fitres_err[var][fs][0]));
      pv1->AddText(Form("DCB sigma: %.3f #pm %.3f", fitres[var][fs][1], fitres_err[var][fs][1]));
      pv1->AddText(Form("DCB a1: %.3f #pm %.3f",    fitres[var][fs][2], fitres_err[var][fs][2]));
      pv1->AddText(Form("DCB n1: %.3f #pm %.3f",    fitres[var][fs][3], fitres_err[var][fs][3]));
      pv1->AddText(Form("DCB a2: %.3f #pm %.3f",    fitres[var][fs][4], fitres_err[var][fs][4]));
      pv1->AddText(Form("DCB n2: %.3f #pm %.3f",    fitres[var][fs][5], fitres_err[var][fs][5]));
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

      c->SaveAs((outputPath4lScaleFitPlots_perCat + "/" + hin4l[var][fs]->GetName() + ".pdf").c_str());
      c->SaveAs((outputPath4lScaleFitPlots_perCat + "/" + hin4l[var][fs]->GetName() + ".png").c_str());

      
      // store 4l fit res 
      vec_4lfitRes[var][fs] = mean_4lDCB.getVal(); 

    }
  }


  //***********************************************************
  // compute 4l scale 
  float vec_4lScale[nFinalStates];
  for(int fs=0; fs<nFinalStates; fs++){

    vec_4lScale[fs] = ( vec_4lfitRes[distrVar][fs] - vec_4lfitRes[distrNominal][fs] ) / vec_4lfitRes[distrNominal][fs];

    cout<< sFinalState[fs] <<": "<<vec_4lScale[fs] <<endl;
  }
  //***********************************************************


  //**********************
  // compare 4l fit plots 
  for(int fs=0; fs<nFinalStates; fs++){

    TCanvas *canv = new TCanvas("canv","canv");
    canv->cd();
  
    frame4l[distrNominal][fs]->Draw();
    frame4l[distrVar][fs]->Draw("same");

    // draw fit results on canvas (nominal)
    TPaveText* pv1 = new TPaveText(0.10,0.51,0.41,0.88,"brNDC");
    pv1->AddText("Nominal distribution");
    pv1->AddText(Form("DCB mean: %.3f #pm %.3f",  fitres[distrNominal][fs][0], fitres_err[distrNominal][fs][0]));
    pv1->AddText(Form("DCB sigma: %.3f #pm %.3f", fitres[distrNominal][fs][1], fitres_err[distrNominal][fs][1]));
    pv1->AddText(Form("DCB a1: %.3f #pm %.3f",    fitres[distrNominal][fs][2], fitres_err[distrNominal][fs][2]));
    pv1->AddText(Form("DCB n1: %.3f #pm %.3f",    fitres[distrNominal][fs][3], fitres_err[distrNominal][fs][3]));
    pv1->AddText(Form("DCB a2: %.3f #pm %.3f",    fitres[distrNominal][fs][4], fitres_err[distrNominal][fs][4]));
    pv1->AddText(Form("DCB n2: %.3f #pm %.3f",    fitres[distrNominal][fs][5], fitres_err[distrNominal][fs][5]));
    pv1->SetFillColor(kWhite);
    pv1->SetBorderSize(1);
    pv1->SetTextColor(kBlue);
    pv1->SetTextFont(42);
    pv1->SetTextSize(0.037);
    pv1->SetTextAlign(12); // text left aligned 
    pv1->Draw();

    // draw fit results on canvas (var)
    TPaveText* pv2 = new TPaveText(0.66,0.51,0.97,0.88,"brNDC");
    pv2->AddText("var distribution");
    pv2->AddText(Form("DCB mean: %.3f #pm %.3f",  fitres[distrVar][fs][0], fitres_err[distrVar][fs][0]));
    pv2->AddText(Form("DCB sigma: %.3f #pm %.3f", fitres[distrVar][fs][1], fitres_err[distrVar][fs][1]));
    pv2->AddText(Form("DCB a1: %.3f #pm %.3f",    fitres[distrVar][fs][2], fitres_err[distrVar][fs][2]));
    pv2->AddText(Form("DCB n1: %.3f #pm %.3f",    fitres[distrVar][fs][3], fitres_err[distrVar][fs][3]));
    pv2->AddText(Form("DCB a2: %.3f #pm %.3f",    fitres[distrVar][fs][4], fitres_err[distrVar][fs][4]));
    pv2->AddText(Form("DCB n2: %.3f #pm %.3f",    fitres[distrVar][fs][5], fitres_err[distrVar][fs][5]));
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
    CMS_lumi(canv,0,0);

    // write 4lepton scale on plots
    TPaveText* pv3 = new TPaveText(0.75,0.3,0.95,0.4,"brNDC");
    pv3->AddText((sFinalState[fs] + " scale:").c_str());
    pv3->AddText(Form("%.6f",vec_4lScale[fs]));
    pv3->SetFillColor(kWhite);
    pv3->SetBorderSize(1);
    pv3->SetTextFont(42);
    pv3->SetTextSize(0.037);
    pv3->SetTextAlign(12); // text left aligned 
    pv3->Draw();

    canv->Update();

    canv->SaveAs((outputPath4lScaleFitPlots_perCat + "/scale_" + sFinalState[fs] + ".pdf").c_str());
    canv->SaveAs((outputPath4lScaleFitPlots_perCat + "/scale_" + sFinalState[fs] + ".png").c_str());


  } //end for on final states 
  //**********************

} // end compute 4l scale function




// *** main function
void ComputeLeptonScaleSyst()
{
 
  string inputPathData = "/data3/Higgs/180416/";
  string inputPathMC_DY = "/data3/Higgs/180416/MC_main/";
  string inputPathMC_ggH = "/data3/Higgs/180416/MC_main/";

  string outputPathFitResultsPlots = "plotsSysts_FitResults";
  string outputPathDileptonScalePlots = "plotsSysts_DileptonScale";
  string outputPathCompare2lDataMcFit = "plotsSysts_CompareDataMC2lFit";
  string outputPath4lScaleFitPlots_perCat = "plotsSysts_4leptonScaleFits_perCatVariation";
  string outputPath4lScaleFitPlots_maxVar = "plotsSysts_4leptonScaleFits_maxVariation";
  

  float lumi = 41.30; //fb-1
  string lumiText = "41.30 fb^{-1}";


  // create output directories
  if(REDOTHE2lFIT) gSystem->Exec(("mkdir -p "+outputPathFitResultsPlots).c_str());  //dir for fit results plots

  if(COMPUTE2lSCALE) gSystem->Exec(("mkdir -p "+outputPathDileptonScalePlots).c_str()); //dir for dilepton scale plots

  if(COMPARE2lDATAMCFIT) gSystem->Exec(("mkdir -p "+outputPathCompare2lDataMcFit).c_str()); //dir for comparison between Data and MC 2l fit

  if(COMPUTE4lSCALE) gSystem->Exec(("mkdir -p "+outputPath4lScaleFitPlots_perCat).c_str()); //dir for 4lepton scale fit plots per categ
  if(COMPUTE4lSCALE) gSystem->Exec(("mkdir -p "+outputPath4lScaleFitPlots_maxVar).c_str()); //dir for 4lepton scale fit plots max var 


  // execute functions 
  if(REDO2lHISTOS) do2lHistograms(inputPathMC_DY, inputPathData, lumi);

  if(REDOTHE2lFIT) doThe2lFit(outputPathFitResultsPlots, lumiText);

  if(COMPUTE2lSCALE) computeDileptonScale(outputPathDileptonScalePlots, lumiText);

  if(COMPARE2lDATAMCFIT) compareDataMCfitPlots(outputPathCompare2lDataMcFit, lumiText);

  // compute 4l scale assigning lepton pT variation according to lep pT and eta 
  if(REDO4lHISTOS) do4lHistograms_perCatVariation(inputPathMC_ggH, lumi);  
  if(COMPUTE4lSCALE) compute4lScale_perCatVariation(outputPath4lScaleFitPlots_perCat, lumiText); 

  // compute 4l scale assigning as lepton pT variation the maximum value obtained per category  
  // if(REDO4lHISTOS) do4lHistograms_maxVariation(inputPathMC_ggH, lumi); 
  // if(COMPUTE4lSCALE) compute4lScale_maxVariation(outputPath4lScaleFitPlots_maxVar, lumiText); 

  
}



