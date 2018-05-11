// **********************************
// 
// usage:
//     - specify parameters (input/output directory, luminosity) at the end of this file
//     - run with:
//                 root -l -b -q ComputeLeptonScaleSyst_TestOnFit.C++
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
#include "RooExtendPdf.h"
#include "RooBreitWigner.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "RooAbsCollection.h"
#include "RooNumConvPdf.h"
#include "TLorentzVector.h"

#include "Plotter/CMS_lumi.C" // CMS official label definition (must be in ZZAnalysis/AnalysisStep/test to access this)
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"  // contains DCB definition


using namespace std;
using namespace RooFit ;


#define WRITEEXTRATEXTONPLOTS 1 // draw Preliminary on Plots


// *** global definitions
Float_t Zmass_nominalPDG  = 91.1876; // GeV/c^2 (http://pdg.lbl.gov/2017/listings/rpp2017-list-z-boson.pdf)
Float_t Zwidth_nominalPDG = 2.4952;  // GeV/c^2 (http://pdg.lbl.gov/2017/listings/rpp2017-list-z-boson.pdf)



// perform the fit 
void doThe2lFit_roofit(string outputPathFitResultsPlots, string lumiText) 
{

  // read file with histos
  TFile* fInHistos = TFile::Open("file_DataMCHistos_varUp_frompT0_varpT001_CUT.root");

  // define input histos
  TH1F* inputhist = (TH1F*)fInHistos->Get("hist_DYJetsToLL_M50_ele_Eta_0.8-1.5_pT_40-50");
   
       
  // define roofit variable for the fit
  RooRealVar mll = RooRealVar("mll","m_{l^{+}l^{-}}",60,120);

  // make roofit datasets from root histos
  RooDataHist dh("dh","dh",mll,Import(*inputhist));
  

  // **** PREFIT ****
  // set range for the prefit
  mll.setRange("range80100gev",80,100); 

  // define prefit function: gaussian
  RooRealVar mean_prefit_DCB("mean_prefit_DCB","mean_DCB",Zmass_nominalPDG,80.,100.);
  RooRealVar sigma_prefit_DCB("sigma_prefit_DCB","sigma_DCB",Zwidth_nominalPDG,0.0001,5.);
  RooGaussian gauss_prefit_pdf("gauss_prefit_pdf","gaussian function",mll,mean_prefit_DCB,sigma_prefit_DCB);

  // do the fit 
  gauss_prefit_pdf.chi2FitTo(dh,Range("range80100gev")); //a binned ML fit (done with fitTo) always assume Poisson error interpretation of data
                                                         //data with non-unit wieghts can only be correctly fitted with a Chi2 fit (chi2FitTo)
  

  // **** FIT ****
  // set range for the fit 
  // define fit range 
  float fitRange_min = mean_prefit_DCB.getVal() - sigma_prefit_DCB.getVal();
  float fitRange_max = mean_prefit_DCB.getVal() + sigma_prefit_DCB.getVal();
  
  // set range for the fit
  mll.setRange("range_fit", fitRange_min, fitRange_max); 
        
  // define fit function: gaussian
  RooRealVar mean_DCB("mean_DCB","mean_DCB",Zmass_nominalPDG, fitRange_min, fitRange_max);
  RooRealVar sigma_DCB("sigma_DCB","sigma_DCB",Zwidth_nominalPDG,0.0001,5.);
  RooGaussian gauss_pdf("gauss_pdf","gaussian function",mll,mean_DCB,sigma_DCB);

  // do the fit 
  gauss_pdf.chi2FitTo(dh,Range("range_fit")); //a binned ML fit (done with fitTo) always assume Poisson error interpretation of data
                                              //data with non-unit wieghts can only be correctly fitted with a Chi2 fit (chi2FitTo)
                
  // plot data on the frame
  RooPlot* frame = mll.frame();
  dh.plotOn(frame,DataError(RooAbsData::SumW2));
  gauss_prefit_pdf.plotOn(frame,LineColor(kBlue)); 
  gauss_pdf.plotOn(frame,LineColor(kRed)); 


  // plot on the canvas and save plots
  TCanvas* c = new TCanvas(); 
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

  c->SaveAs((outputPathFitResultsPlots + "/hist_DYJetsToLL_M50_ele_Eta_0.8-1.5_pT_40-50_RooFit.png").c_str());

 
}// end doTheFit_roofit function 





//****************
// perform the fit 
void doThe2lFit_root(string outputPathFitResultsPlots, string lumiText)
{

  // read file with histos
  TFile* fInHistos = TFile::Open("file_DataMCHistos_varUp_frompT0_varpT001_CUT.root");

  // define input histos
  TH1F* inputhist = (TH1F*)fInHistos->Get("hist_DYJetsToLL_M50_ele_Eta_0.8-1.5_pT_40-50");

          
  // **** PREFIT ****
  // define fit function
  TF1* g1 = new TF1("g1","gaus", 80., 100.);
  g1->SetParameter(1,Zmass_nominalPDG);
  g1->SetParameter(2,Zwidth_nominalPDG);
  g1->SetLineColor(kBlue);
  
  // do the fit 
  inputhist->Fit("g1", "R");


  // **** FIT ****
  // set range for the fit 
  // define fit range 
  float fitRange_min = g1->GetParameter(1) - g1->GetParameter(2); //param 0 = norm, param 1 = mean, param 2 = sigma 
  float fitRange_max = g1->GetParameter(1) + g1->GetParameter(2);
   
  // define fit function 
  TF1* g2 = new TF1("g2","gaus", fitRange_min, fitRange_max);
  g2->SetParameter(1,Zmass_nominalPDG);
  g2->SetParameter(2,Zwidth_nominalPDG);
  g2->SetLineColor(kRed);
                
  // do the fit 
  inputhist->Fit("g2", "R"); 
 

  // plot on the canvas and save plots
  TCanvas* c = new TCanvas(); 
  c->cd();
  inputhist->SetMarkerStyle(20);
  inputhist->Draw("pe");
  g1->Draw("samel");
  g2->Draw("samel");


  gStyle->SetOptStat(0);
  

  // draw fit results on canvas
  TPaveText* pv = new TPaveText(0.64,0.65,0.95,0.87,"brNDC");
  pv->AddText(Form("mean prefit: %.3f #pm %.3f", g1->GetParameter(1), g1->GetParError(1))); ((TText*)pv->GetListOfLines()->Last())->SetTextColor(kBlue);
  pv->AddText(Form("sigma prefit: %.3f #pm %.3f", g1->GetParameter(2), g1->GetParError(2))); ((TText*)pv->GetListOfLines()->Last())->SetTextColor(kBlue);
  pv->AddText(Form("mean fit: %.3f #pm %.3f", g2->GetParameter(1), g2->GetParError(1))); ((TText*)pv->GetListOfLines()->Last())->SetTextColor(kRed);
  pv->AddText(Form("sigma fit: %.3f #pm %.3f", g2->GetParameter(2), g2->GetParError(2))); ((TText*)pv->GetListOfLines()->Last())->SetTextColor(kRed);
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
  
  c->SaveAs((outputPathFitResultsPlots + "/hist_DYJetsToLL_M50_ele_Eta_0.8-1.5_pT_40-50_Root.png").c_str());


  
}// end doTheFit_root function 






// *** main function
void ComputeLeptonScaleSyst_TestOnFit()
{
 
  string outputPathFitResultsPlots = "plotsSysts_FitResults_test";
    
  float lumi = 41.30; //fb-1
  string lumiText = "41.30 fb^{-1}";


  // create output directories
  gSystem->Exec(("mkdir -p "+outputPathFitResultsPlots).c_str());  //dir for fit results plots
  


  // execute functions 
  doThe2lFit_roofit(outputPathFitResultsPlots, lumiText);

  doThe2lFit_root(outputPathFitResultsPlots, lumiText);
  

  
}



