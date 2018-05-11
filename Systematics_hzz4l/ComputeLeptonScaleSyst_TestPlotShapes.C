//*************************
//
// usage:
//     root -l -b -q ComputeLeptonScaleSyst_TestPlotShapes.C++
//
//*************************

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

Float_t Zmass_nominalPDG  = 91.1876; // GeV/c^2 (http://pdg.lbl.gov/2017/listings/rpp2017-list-z-boson.pdf)
Float_t Zwidth_nominalPDG = 2.4952;  // GeV/c^2 (http://pdg.lbl.gov/2017/listings/rpp2017-list-z-boson.pdf)


void ComputeLeptonScaleSyst_TestPlotShapes()
{

  // define roofit variable for the fit
  RooRealVar mll = RooRealVar("mll","m_{l^{+}l^{-}}",60,120);

  // set fit range
  //mll.setRange("range80100gev",80,100); 

  
  // define fit function: 
  // convolution of a BW function and a DoubleCrystalBall resolution function 
  // ---Breit-Wigner function (physical mass distribution)
  RooRealVar pole_BW("pole_BW","pole_BW",Zmass_nominalPDG);
  RooRealVar width_BW("width_BW","width_BW",Zwidth_nominalPDG);

  pole_BW.setConstant(Zmass_nominalPDG);    //fix to physical value
  width_BW.setConstant(Zwidth_nominalPDG);  //fix to physical value

  RooBreitWigner BW_pdf("BW_pdf","Breit-Wigner function",mll,pole_BW,width_BW);

  // ---Double Crystal Ball function (experimental resolution)
  RooRealVar mean_DCB("mean_DCB","mean_DCB",0.);//,-2.,2.);
  RooRealVar sigma_DCB("sigma_DCB","sigma_DCB",1.);//,0.0001,5.);
  RooRealVar a1_DCB("a1_DCB","a1_DCB",1.);//,0.,10);
  RooRealVar n1_DCB("n1_DCB","n1_DCB",2.);//,0.,10);
  RooRealVar a2_DCB("a2_DCB","a2_DCB",1.);//,0.,10);
  RooRealVar n2_DCB("n2_DCB","n2_DCB",2.);//,0.,10);
  
  RooDoubleCB DCBall_pdf("DCBall_pdf","Double Crystal ball function",mll,mean_DCB,sigma_DCB,a1_DCB,n1_DCB,a2_DCB,n2_DCB);

  //---Convolution function
  RooNumConvPdf tot_pdf("tot_pdf","Total PDF",mll,DCBall_pdf,BW_pdf);


  RooDataSet *data = tot_pdf.generate(RooArgSet(mll),500);

  tot_pdf.fitTo(*data);


  // plot data on the frame
  RooPlot* frame = mll.frame();
  //frame->SetName(inputhist[dat][catEta][catPt]->GetName());  // name is the name which appears in the root file
  //frame->SetTitle(inputhist[dat][catEta][catPt]->GetName()); // title is the title on the canvas  
  //dh.plotOn(frame,DataError(RooAbsData::SumW2)); 
  data->plotOn(frame);
  tot_pdf.plotOn(frame,LineColor(kRed)); 
  

  // plot on the canvas and save plots
  TCanvas* c = new TCanvas(); 
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

  

  c->Update();

  gSystem->Exec("mkdir -p plotSysts_shapeControl");

  c->SaveAs("plotSysts_shapeControl/shape_BWandDCB_mean0sigma1.png");
  

}
