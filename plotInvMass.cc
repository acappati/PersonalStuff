// **************************
//
// purpose: plot 4l invariant mass (very basic version) 
//
// usage:
//        root
//        .L plotInvMass.cc++
//        plotInvMass()
//
// **************************

#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include <vector>
#include <iostream>


using namespace std;


//Function(s) prototype(s)
TGraphAsymmErrors* tgMaker(const TH1F * h);
void readFileMakeHisto(TString infile, TH1F *hist, Float_t L, Int_t k);



//Macro
void plotInvMass()
{
  Double_t L=35.9; //[{fb}^{-1}]
  Int_t nBin=225;//207
  Double_t xMin=70.; //70
  Double_t xMax=970.;//898
    
  //***************************
  //*********  qqZZ ***********
  //***************************
  TH1F *hmass0=new TH1F("hmass0","ZZMass",nBin,xMin,xMax);
  readFileMakeHisto("/data3/Higgs/170222/ZZTo4l/ZZ4lAnalysis.root",hmass0,L,0);
  //****************************
  
  
  //****************************
  //*********  ggH125  *********
  //****************************
  TH1F *hmass1=new TH1F("hmass1","ZZMass",nBin,xMin,xMax);
  readFileMakeHisto("/data3/Higgs/170222/ggH125/ZZ4lAnalysis.root",hmass1,L,1);
  //****************************
  

  //****************************
  //*********  Data  ***********
  //****************************
  TH1F *hmass2=new TH1F("hmass2","ZZMass",nBin,xMin,xMax);
  //readFileMakeHisto("/data3/Higgs/160803_complete_newEbE/AllData/ZZ4lAnalysis.root",hmass2,L,2); //L=12.9 [{fb}^{-1}]
  readFileMakeHisto("/data3/Higgs/170222/AllData/ZZ4lAnalysis.root",hmass2,L,2);
  //readFileMakeHisto("/data3/Higgs/170712_Data2017/AllData/ZZ4lAnalysis.root",hmass2,L,2);  
  //****************************


  //****************************
  //*********  ggZZ  ***********
  //****************************
  //hist from 3 to 8
  
  //*******  2e2mu  ************
  TH1F *hmass3=new TH1F("hmass3","ZZMass",nBin,xMin,xMax);
  readFileMakeHisto("/data3/Higgs/170222/ggTo2e2mu_Contin_MCFM701/ZZ4lAnalysis.root",hmass3,L,3);
 
  //*******  2e2tau  ************
  TH1F *hmass4=new TH1F("hmass4","ZZMass",nBin,xMin,xMax);
  readFileMakeHisto("/data3/Higgs/170222/ggTo2e2tau_Contin_MCFM701/ZZ4lAnalysis.root",hmass4,L,3);
  
  //*******  2mu2tau  ************
  TH1F *hmass5=new TH1F("hmass5","ZZMass",nBin,xMin,xMax);
  readFileMakeHisto("/data3/Higgs/170222/ggTo2mu2tau_Contin_MCFM701/ZZ4lAnalysis.root",hmass5,L,3);
  
  //*******  4e  ************
  TH1F *hmass6=new TH1F("hmass6","ZZMass",nBin,xMin,xMax);
  readFileMakeHisto("/data3/Higgs/170222/ggTo4e_Contin_MCFM701/ZZ4lAnalysis.root",hmass6,L,3);

  //*******  4mu  ************
  TH1F *hmass7=new TH1F("hmass7","ZZMass",nBin,xMin,xMax);
  readFileMakeHisto("/data3/Higgs/170222/ggTo4mu_Contin_MCFM701/ZZ4lAnalysis.root",hmass7,L,3);

  //*******  4tau  ************
  TH1F *hmass8=new TH1F("hmass8","ZZMass",nBin,xMin,xMax);
  readFileMakeHisto("/data3/Higgs/170222/ggTo4tau_Contin_MCFM701/ZZ4lAnalysis.root",hmass8,L,3);

  
  
  //******************************
  //sum and draw histo
  THStack *hs = new THStack("hs","");
  TCanvas *c2=new TCanvas("c2","c2");
  TH1F *hr=c2->DrawFrame(xMin,0.,xMax,110.);
  hr->SetXTitle("m_{4l} (GeV)");
  hr->SetYTitle("Events/4 GeV");
  //hr->SetTitle();
  hr->GetXaxis()->SetTitleOffset(1.2);
  hr->GetYaxis()->SetTitleOffset(1.2);
   
  gPad->SetLogx();
  hr->GetXaxis()->SetMoreLogLabels();
  hr->GetXaxis()->SetNoExponent();
  gPad->SetTickx();
  gPad->SetTicky();
 
  

  //Dark Blue histo
  TH1F*hmass10 = new TH1F(*hmass3+*hmass4+*hmass5+*hmass6+*hmass7+*hmass8);
  hmass10->SetLineColor(kBlue);
  hmass10->SetFillColor(TColor::GetColor("#3366ff"));
  hs->Add(hmass10);

  
  //Blue histo
  hmass0->SetLineColor(kBlue);
  hmass0->SetFillColor(TColor::GetColor("#8bc5ff"));
  hs->Add(hmass0);


  //Red histo
  hmass1->SetLineColor(kRed);
  hmass1->SetFillColor(TColor::GetColor("#ff9090"));
  hs->Add(hmass1);

  
  
  //All Data histogram
  hmass2->SetBinErrorOption(TH1::kPoisson);
 
  //call TGMaker Function
  TGraphAsymmErrors* tgmass2 = tgMaker(hmass2);

  tgmass2->SetName("tgmass2");
  tgmass2->SetMarkerStyle(20);
  tgmass2->SetMarkerSize(0.6);
  tgmass2->SetMarkerColor(1);
  tgmass2->SetLineColor(1);
    
  
  hs->Draw("hsame");
  tgmass2->Draw("P");
  


  //******************************
  //Legend
  TLegend *leg=new TLegend(0.63,0.68,0.83,0.82);
  leg->AddEntry("tgmass2","Data","pe");
  leg->AddEntry(hmass1,"H(125)","f");  
  leg->AddEntry(hmass0,"q#bar{q} #rightarrow ZZ, Z#gamma *","f");
  leg->AddEntry(hmass10,"gg #rightarrow ZZ, Z#gamma *","f");
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->Draw();
  //******************************

  c2->SaveAs("InvMass.png");

}// end plotInvMass

//******************************************************************

//Function(s) definition
//tgMaker Function: takes a histogram and returns a TGraphAsymmErrors  
TGraphAsymmErrors* tgMaker(const TH1F * h)
{

  Int_t nBins=h->GetNbinsX();
  
  vector <Double_t> vx;
  vector <Double_t> vy;
  vector <Double_t> vErrXUp;
  vector <Double_t> vErrXLow;
  vector <Double_t> vErrYUp;
  vector <Double_t> vErrYLow;
  Double_t xtemp, ytemp;
  Int_t a=0;
  Int_t n=0;

  
  for(Int_t i=1; i<=nBins; i++)
   {
     xtemp=h->GetBinCenter(i);
     ytemp=h->GetBinContent(i);

     
     if(ytemp>0.)
      {      
        vx.push_back(xtemp);
        vy.push_back(ytemp);

        vErrXUp.push_back(0.);
        vErrXLow.push_back(0.);

        vErrYUp.push_back(h->GetBinErrorUp(i));
        vErrYLow.push_back(h->GetBinErrorLow(i));
      }//end if
     else
      {    
        a++; 
      }//end else

   }//end for cycle


  n=nBins-a;
  
  Double_t x[n];
  Double_t y[n];
  Double_t ErrXUp[n];
  Double_t ErrXLow[n];
  Double_t ErrYUp[n];
  Double_t ErrYLow[n];
  
  for(Int_t j=1; j<=n; j++)
   {
     x[j]=vx[j];
     y[j]=vy[j];

     ErrXUp[j]=vErrXUp[j];
     ErrXLow[j]=vErrXLow[j];

     ErrYUp[j]=vErrYUp[j];
     ErrYLow[j]=vErrYLow[j];
   }

  return new TGraphAsymmErrors(n, x, y, ErrXLow, ErrXUp, ErrYLow, ErrYUp);
 

}//end tgMaker


//******************************************************************


//readFileMakehisto: takes a string (file path) and fills histogram with weight
void readFileMakeHisto(TString infile, TH1F *hist, Float_t L, Int_t k)
{ 

  Float_t ZZMass;
  Short_t ZZsel;
  Float_t xsec; //values are in [pb]!!
  Float_t PUWeight; 
  Float_t dataMCWeight;
  Float_t KFactor_QCD_ggZZ_Nominal; //ggZZ only
  Float_t KFactor_EW_qqZZ;          //qqZZ only
  Float_t KFactor_QCD_qqZZ_M;       //qqZZ only
  Double_t W=1.;

  //k=0 for qqZZ
  //k=1 for ggH125
  //k=2 for Data
  //k=3 for ggZZ
  
  //open input file
  TFile *file=new TFile(infile);

  //get Number of generated events
  TH1F *h=(TH1F*)file->Get("ZZTree/Counters");
  Double_t Ngen=h->GetBinContent(1);
   
  //read TTree and get ZZMass histogram
  TTree *tree = (TTree*)file->Get("ZZTree/candTree");
  //set branches address
  tree->SetBranchAddress("ZZMass",&ZZMass);
  tree->SetBranchAddress("ZZsel",&ZZsel);


  if(k!=2)
  { tree->SetBranchAddress("xsec",&xsec);
    tree->SetBranchAddress("PUWeight",&PUWeight);
    tree->SetBranchAddress("dataMCWeight",&dataMCWeight); }

  if(k==0)
  { tree->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ);
    tree->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M); }
  
  if(k==3) 
  { tree->SetBranchAddress("KFactor_QCD_ggZZ_Nominal",&KFactor_QCD_ggZZ_Nominal); }
  

  
  for(Int_t i=0; i<tree->GetEntries(); i++)
   {
     tree->GetEntry(i);

     if( !(ZZsel>=90) ) continue;

     //weight
     if(k==2)
     { hist->Fill(ZZMass); }
     else
     {
       if(k==0)
        { W=(xsec*1000.*L/Ngen)*PUWeight*dataMCWeight*KFactor_EW_qqZZ*KFactor_QCD_qqZZ_M; }
       
       if(k==1)
	    { W=(xsec*1000.*L/Ngen)*PUWeight*dataMCWeight; }
       
       if(k==3)
	    { W=(xsec*1000.*L/Ngen)*PUWeight*dataMCWeight*KFactor_QCD_ggZZ_Nominal; }

	   hist->Fill(ZZMass,W);
	
      }//end else
        
   }//end for
  file->Close();

}//end readFilemakeHisto

//******************************************************************
