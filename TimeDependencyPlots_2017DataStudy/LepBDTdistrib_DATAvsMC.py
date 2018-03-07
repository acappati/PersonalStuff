#!/usr/bin/en math

# *******************
# usage: 
#    python LepBDTdistrib_DATAvsMC.py
#
# ******************

import json
import ROOT, math, helper, CMSGraphics, CMS_lumi
from ROOT import TFile, TH1F, TCanvas, gSystem, TAttFill, TLegend, TRatioPlot, TPad, TStyle, THStack, TPaveText, gStyle
from helper import ReadJSON
from CMSGraphics import makeCMSCanvas, makeLegend
from ROOT import kBlue, kRed, kBlack, kWhite, kAzure, kOrange, kMagenta


# *****************************
# Declare all the variables

# data tree options 
ZZTree   = False
CRZLTree = False
ZTree    = True

# data periods options
# period = "data2016"
period = "data2017"
# *****************************
lumiText = '41.30 fb^{-1}'
#******************************


if(ZZTree):     treeText  = "ZZTree"
elif(CRZLTree): treeText  = "CRZLTree"
elif(ZTree):    treeText  = "ZTree"
else: print ("Error: wrong option!")

#******************************




# create output directory
OutputPath = "LepBDTdistrib_DATAvsMC_" + str(period) + "_" + str(treeText)
gSystem.Exec("mkdir -p " + OutputPath)
print "Output directory created!"




# *** DATA ***
#read histos from data file

histoDATA_input = TFile.Open("LepBDTdistrib_DATA_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoDATA_input.GetName(),'...'

LepBDTDATA = []

LepBDTDATA.append(histoDATA_input.Get('LepBDT leading ele'))
LepBDTDATA.append(histoDATA_input.Get('LepBDT leading ele in ECAL Barrel'))
LepBDTDATA.append(histoDATA_input.Get('LepBDT leading ele in ECAL Endcap'))

# LepBDTDATA.append(histoDATA_input.Get('LepBDT leading mu'))
# LepBDTDATA.append(histoDATA_input.Get('LepBDT leading mu in Muon Barrel'))
# LepBDTDATA.append(histoDATA_input.Get('LepBDT leading mu in Muon Endcap'))

if CRZLTree :
    LepBDTDATA.append(histoDATA_input.Get('LepBDT extraEl'))
    LepBDTDATA.append(histoDATA_input.Get('LepBDT extraEl in ECAL Barrel'))
    LepBDTDATA.append(histoDATA_input.Get('LepBDT extraEl in ECAL Endcap'))

    # LepBDTDATA.append(histoDATA_input.Get('LepBDT extraMu'))
    # LepBDTDATA.append(histoDATA_input.Get('LepBDT extraMu in Muon Barrel'))
    # LepBDTDATA.append(histoDATA_input.Get('LepBDT extraMu in Muon Endcap'))



# *** MC DY ***
#read histo from histoMC_DY.root

histoMC_input = TFile.Open("LepBDTdistrib_MC_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoMC_input.GetName(),'...'

LepBDTMC = []

LepBDTMC.append(histoMC_input.Get('LepBDT leading ele'))
LepBDTMC.append(histoMC_input.Get('LepBDT leading ele in ECAL Barrel'))
LepBDTMC.append(histoMC_input.Get('LepBDT leading ele in ECAL Endcap'))

# LepBDTMC.append(histoMC_input.Get('LepBDT leading mu'))
# LepBDTMC.append(histoMC_input.Get('LepBDT leading mu in Muon Barrel'))
# LepBDTMC.append(histoMC_input.Get('LepBDT leading mu in Muon Endcap'))

if CRZLTree :
    LepBDTMC.append(histoMC_input.Get('LepBDT extraEl'))
    LepBDTMC.append(histoMC_input.Get('LepBDT extraEl in ECAL Barrel'))
    LepBDTMC.append(histoMC_input.Get('LepBDT extraEl in ECAL Endcap'))

    # LepBDTMC.append(histoMC_input.Get('LepBDT extraMu'))
    # LepBDTMC.append(histoMC_input.Get('LepBDT extraMu in Muon Barrel'))
    # LepBDTMC.append(histoMC_input.Get('LepBDT extraMu in Muon Endcap'))



if CRZLTree :
    nameList = ['LepBDT_leading_ele',
                'LepBDT_leading_ele_InECALbarrel',
                'LepBDT_leading_ele_InECALendcap',
                # 'LepBDT_leading_mu',
                # 'LepBDT_leading_mu_InMuonBarrel',
                # 'LepBDT_leading_mu_InMuonEndcap',
                'LepBDT_extraEl',
                'LepBDT_extraEl_InECALbarrel',
                'LepBDT_extraEl_InECALendcap',
                # 'LepBDT_extraMu',
                # 'LepBDT_extraMu_InMuonBarrel',
                # 'LepBDT_extraMu_InMuonEndcap'
                ]
else : 
    nameList = ['LepBDT_leading_ele',
                'LepBDT_leading_ele_InECALbarrel',
                'LepBDT_leading_ele_InECALendcap',
                # 'LepBDT_leading_mu',
                # 'LepBDT_leading_mu_InMuonBarrel',
                # 'LepBDT_leading_mu_InMuonEndcap'
                ]
    
#canvas = []
                                   
# *** do the plots ***
for i in range(len(nameList)) :

    canvas = TCanvas("canvas","canvas",800,800)

    #normalization
    norm = 1 # norm MC to xsec
    # norm = LepBDTDATA[i].Integral() / LepBDTMC[i].Integral()   # norm MC to data

    #DATA hist
    LepBDTDATA[i].SetMarkerStyle(20)
    LepBDTDATA[i].SetMarkerSize(0.6)

    #MC hist
    LepBDTMC[i].Scale(norm) #normalize MC 
    LepBDTMC[i].SetFillColor(kMagenta-6)
    LepBDTMC[i].SetLineColor(kBlack)


    #upper plot pad
    pad1 = TPad("pad1","pad1", 0, 0.3, 1, 1.0)
    pad1.Draw()
    pad1.cd()


    LepBDTMC[i].Draw("histo")
    LepBDTDATA[i].Draw("sameEP") 


    LepBDTMC[i].GetXaxis().SetTitle("BDT score")
    LepBDTMC[i].GetXaxis().SetLabelFont(43)
    LepBDTMC[i].GetXaxis().SetLabelSize(15)
    LepBDTMC[i].GetYaxis().SetTitleSize(20)
    LepBDTMC[i].GetYaxis().SetTitleFont(43)
    LepBDTMC[i].GetYaxis().SetTitleOffset(1.8)
    LepBDTMC[i].GetYaxis().SetLabelFont(43)
    LepBDTMC[i].GetYaxis().SetLabelSize(15)
    LepBDTMC[i].GetYaxis().SetTitle("Events")

    LepBDTMC[i].SetTitle("")
    
    gStyle.SetOptStat(0)

    pad1.SetLogy() 


    # legend
    legend = TLegend(0.24,0.68,0.38,0.82)
    legend.AddEntry(LepBDTDATA[i],"Data", "p")
    legend.AddEntry(LepBDTMC[i],"DY MC","f")
    legend.SetFillColor(kWhite)
    legend.SetLineColor(kBlack)
    legend.SetTextFont(43)
    legend.SetTextSize(20)
    legend.Draw()
  

    canvas.Update()


    #lower plot pad
    canvas.cd()
    pad2 = TPad("pad2","pad2", 0, 0.05, 1, 0.3)
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()    #pad2 becomes the current pad


    #define ratio plot
    rp = TH1F(LepBDTDATA[i].Clone("rp"))
    rp.SetLineColor(kBlack)
    rp.SetMinimum(0.5)
    rp.SetMaximum(2.)
    rp.SetStats(0)
    rp.Divide(TH1F(LepBDTMC[i]))   #divide histo rp/MC
    rp.SetMarkerStyle(24)
    rp.SetTitle("") 
    
    rp.SetYTitle("Data/MC")
    rp.GetYaxis().SetNdivisions(505)
    rp.GetYaxis().SetTitleSize(20)
    rp.GetYaxis().SetTitleFont(43)
    rp.GetYaxis().SetTitleOffset(1.55)
    rp.GetYaxis().SetLabelFont(43)
    rp.GetYaxis().SetLabelSize(15)

    rp.GetXaxis().SetTitleSize(20)
    rp.GetXaxis().SetTitleFont(43)
    rp.GetXaxis().SetTitleOffset(4.)
    rp.GetXaxis().SetLabelFont(43)
    rp.GetXaxis().SetLabelSize(15)

    rp.Draw("ep")


    canvas.Update()


    #draw CMS and lumi text
    CMS_lumi.writeExtraText = True
    CMS_lumi.extraText      = "Preliminary"
    CMS_lumi.lumi_sqrtS     = lumiText + " (13 TeV)"
    CMS_lumi.cmsTextSize    = 0.6
    CMS_lumi.lumiTextSize   = 0.46
    CMS_lumi.extraOverCmsTextSize = 0.75
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 0, 0)
    
    
    canvas.Update()

     
    canvas.SaveAs(OutputPath + "/" + nameList[i] + ".pdf")
    canvas.SaveAs(OutputPath + "/" + nameList[i] + ".png")


print "plots done"
