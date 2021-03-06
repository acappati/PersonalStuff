#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Math/DistFunc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "plotUtils.C"

#define _DEBUG_ -1


string base_dir("../");

// string inclusive_input_file("../../Datacards13TeV/LegoCards/mlfitHZZ4L_ML_inclusive.root");
string inclusive_input_file(base_dir+"mlfitHZZ4L_ML_unblind_inclusive.root");
// string inclusive_input_file("../../Datacards13TeV/LegoCards/datacards_170306_7Cat_Mu/mlfitHZZ4L_ML_inclusive.root");

// by category
//string type("category");
// Bool_t tweakSpaces = false;
// Bool_t couplingModel = false;
//const Int_t nlines = 7;
string input_files[7] = {
  base_dir+"mlfitHZZ4L_ML_unblind_Untagged.root",
  base_dir+"mlfitHZZ4L_ML_unblind_VBF1JetTagged.root",
  base_dir+"mlfitHZZ4L_ML_unblind_VBF2JetTagged.root",
  base_dir+"mlfitHZZ4L_ML_unblind_VHhadrTagged.root",
  base_dir+"mlfitHZZ4L_ML_unblind_VHLeptTagged.root",
  base_dir+"mlfitHZZ4L_ML_unblind_VHMETtagged.root",
  base_dir+"mlfitHZZ4L_ML_unblind_ttHTagged.root"
};


//by process (5), VH, LepHad
string type("process");
string subtype("VH_LepHad");
Bool_t tweakSpaces = true;
Bool_t couplingModel = false;
const Int_t nlines = 5;
string process_input_file(base_dir+"higgsCombineMultiSignal_VH_LepHad_unblind.MultiDimFit.mH125.09.root");
string branches[nlines] = {
  "r_ggH",
  "r_VBF",
  "r_VH_had",
  "r_VH_lep",
  "r_ttH"
};
//string labels[nlines] = { "#mu_{ggH} ", "#mu_{VBF} ", "#mu_{VH} ", "#mu_{t#bar{t}H} ", };
string labels[nlines] = { "#mu_{ggH} = ", "#mu_{VBF} = ", "#mu_{VHhad} = ", "#mu_{VHlep} = ", "#mu_{t#bar{t}H} = ", };
Float_t leftMargin = 0.29;
Float_t yLabelSize = 0.055;
Int_t colorInclBand = TColor::GetColor("#C8F3C9");
Int_t colorInclLine = TColor::GetColor("#008A05");
Int_t colorErrorBars = TColor::GetColor("#FF9E00");
string labelCanvas = "muProcesses4_obs";
Float_t inclCentral = 0.991566;
Float_t inclDn      = -0.264199;
Float_t inclUp      = +0.329123;
Float_t central[nlines] = {  0.960,  0.670,  1.841,   8.410, };
Float_t dn     [nlines] = { -0.333, -0.670, -1.841,  -8.154, };
Float_t up     [nlines] = { +0.401, +1.613, +6.364, +13.075, };
Float_t xmax = 5.5;
// Float_t xmin = -4.;
Float_t xmin = 0;





void plotCategories(){

  setTDRStyle();
  TFile *incfitfile = TFile::Open(inclusive_input_file.c_str());
  TTree *inctreecomb = (TTree*)incfitfile->Get("tree_fit_sb");
  double mu = 0;
  double muLoErr = 0;
  double muHiErr = 0;
  inctreecomb->SetBranchAddress("mu", &mu);
  inctreecomb->SetBranchAddress("muLoErr", &muLoErr);
  inctreecomb->SetBranchAddress("muHiErr", &muHiErr);
  std::cout<<"reading inclusive \n";
  for(int i=0; i<(int)inctreecomb->GetEntries();++i) 
  {
    inctreecomb->GetEntry(i);
    if(i==0)
    {
       inclCentral = mu;
       // inclCentral = 1.;
       // inclDn = -muLoErr/mu;
       // inclUp = muHiErr/mu;
       inclDn = -muLoErr;
       inclUp = muHiErr;
    }
  }
  incfitfile->Close();
  if(type=="category")
  {
    for(int cha=0; cha<nlines; cha++)
    {
      TFile *fitfile = TFile::Open(input_files[cha].c_str());
      TTree *treecomb = (TTree*)fitfile->Get("tree_fit_sb");
      treecomb->SetBranchAddress("mu", &mu);
      treecomb->SetBranchAddress("muLoErr", &muLoErr);
      treecomb->SetBranchAddress("muHiErr", &muHiErr);
      std::cout<<"reading channel "<<cha<<"\n";
      for(int i=0; i<(int)treecomb->GetEntries();++i) 
      {
	treecomb->GetEntry(i);
	if(i==0)
	{
	  central[cha] = mu;
	  // central[cha] = 1.;
	  // dn[cha] = -muLoErr/mu;
	  // up[cha] = muHiErr/mu;
	  dn[cha] = -muLoErr;
	  up[cha] = muHiErr;
	}
      }
      fitfile->Close();
    }
  }
  else if(type=="process")
  {
    TFile *fitfile = TFile::Open(process_input_file.c_str());
    TTree *treecomb = (TTree*)fitfile->Get("limit");
    float vars[nlines];
    for(int proc=0;proc<nlines;proc++)
    {
      treecomb->SetBranchAddress(branches[proc].c_str(), &vars[proc]);
    }
    float tmp_central[nlines];
    // first loop on TTree entries to get the central value
    for(int i=0; i<(int)treecomb->GetEntries();++i) 
    {
      treecomb->GetEntry(i);

      if(i==10)
      {
        for(int proc=0;proc<nlines;proc++)
	{
          central[proc] = vars[proc];
          dn[proc] = 0;
          up[proc] = 0;
        }
      }
    }
    // second loop on TTree entries to assign error bars
    for(int i=0; i<(int)treecomb->GetEntries();++i) 
    {
      treecomb->GetEntry(i);

      if(i!=10)
      {
        for(int proc=0;proc<nlines;proc++)
	{
          if(vars[proc]-central[proc]>up[proc]) up[proc] = vars[proc] - central[proc];
          if(vars[proc]-central[proc]<dn[proc]) dn[proc] = vars[proc] - central[proc];
        }
      }

    }
      
    fitfile->Close();

  }
  gStyle->SetGridColor(kGray);

  TCanvas* c = new TCanvas(("c_"+labelCanvas).c_str(),"",500,500); 

  TStyle* style = (TStyle*)gROOT->GetStyle("tdrStyle")->Clone();
  style->SetPadLeftMargin(leftMargin);
  style->SetPadRightMargin(0.03);
  style->cd();
  c->cd();
  c->UseCurrentStyle();

  TH2F *h2 = new TH2F("h2",";#mu;",10,xmin,xmax,nlines,0.,nlines);
  h2->GetXaxis()->SetLabelSize(0.045);
  h2->GetXaxis()->SetTitleSize(0.055);
  h2->GetXaxis()->SetTitleOffset(1.05);
  //h2->GetXaxis()->SetNdivisions(-nlines);
  h2->SetStats(0);
  for(int i=1; i<=nlines; i++)
    h2->GetYaxis()->SetBinLabel(nlines+1-i,Form("%s%s%.2f^{#plus%.2f}_{#minus%.2f}%s",labels[i-1].c_str(),((tweakSpaces&&fabs(central[i-1])<10)?" ":""),central[i-1],up[i-1],fabs(dn[i-1]),((tweakSpaces&&fabs(central[i-1])<10)?" ":"")));
  h2->GetYaxis()->SetLabelSize(yLabelSize);
  h2->GetYaxis()->LabelsOption("h");
  h2->GetYaxis()->SetTickLength(0.);
  h2->GetYaxis()->SetNdivisions(-xmax);
  h2->Draw();

  TBox inclBand;
  inclBand.SetFillColor(colorInclBand);
  inclBand.DrawBox(inclCentral+inclDn,0.,inclCentral+inclUp,nlines);
  TLine* inclLine = new TLine(inclCentral,0.,inclCentral,nlines);
  inclLine->SetLineColor(colorInclLine);
  inclLine->SetLineWidth(2);
  inclLine->Draw();

  TLine* errorBar[nlines];
  TMarker* mark[nlines];
  for(int i=0; i<nlines; i++){
    float z = nlines-i-0.5;
    errorBar[i] = new TLine(central[i]+dn[i],z,central[i]+up[i],z);
    errorBar[i]->SetLineColor(colorErrorBars);
    errorBar[i]->SetLineWidth(3);
    errorBar[i]->Draw();
    mark[i] = new TMarker(central[i],z,21);
    mark[i]->SetMarkerColor(kBlack);
    mark[i]->SetMarkerSize(1.);
    std::cout<<central[i]<<"\n";
    // if(central[i]>1e-5) mark[i]->Draw();
    mark[i]->Draw();
  }

  TPaveText* pave = printInfo("H #rightarrow ZZ* #rightarrow 4#font[12]{l}",0.5,0.8,0.98,0.91);
  pave->SetTextSize(0.035);
  pave->SetTextAlign(31);
  pave->SetTextFont(42);
  pave->AddText("m_{H} = 125.09 GeV");

  TPaveText* paveComb = printInfo(Form("#mu_{comb.} = %.2f^{#plus%.2f}_{#minus%.2f}",inclCentral,inclUp,fabs(inclDn)),0.5,0.72,0.98,0.75);
  paveComb->SetTextSize(0.035);
  paveComb->SetTextAlign(31);
  paveComb->SetTextFont(42);
  TPad *p = new TPad("p","p",0.,0.,1.,1.); p->SetFillStyle(0); p->Draw(); p->cd();
  TBox inclBandLeg;
  inclBandLeg.SetFillColor(colorInclBand);
  inclBandLeg.DrawBox(0.68,0.73,0.72,0.76);
  TLine inclLineLeg;
  inclLineLeg.SetLineColor(colorInclLine);
  inclLineLeg.SetLineWidth(2);
  inclLineLeg.DrawLineNDC(0.7,0.73,0.7,0.76);


  //----- Print official CMS labels and lumi
  writeExtraText = false;
  extraText  = "Preliminary";
  string lumiText = "35.9 fb^{-1}";
  lumi_sqrtS = lumiText + " (13 TeV)";
  CMS_lumi( c, 0, 0 );

  c->Print("mu_perProdModes.png");
  c->Print("mu_perProdModes.pdf");

  c->RedrawAxis();
  c->Print(Form("plots/categories_%s.png",type.c_str()));
  c->Print(Form("plots/categories_%s.pdf",type.c_str()));
  c->Print(Form("plots/categories_%s.C",type.c_str()));
  // c->Print(Form("plots/categories_expected_%s.png",type.c_str()));
  // c->Print(Form("plots/categories_expected_%s.pdf",type.c_str()));
  // c->Print(Form("plots/categories_expected_%s.C",type.c_str()));


  //SaveCanvas("PlotsSignificance/",c);
  // SaveCanvas("$pl3/",c);

}

