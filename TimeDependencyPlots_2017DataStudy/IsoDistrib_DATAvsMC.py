#!/usr/bin/en math

# *******************
# usage: 
#    python IsoDistrib_DATAvsMC.py
#
# ******************

import json
import ROOT, math, helper, CMSGraphics, CMS_lumi
from ROOT import TFile, TH1F, TCanvas, gSystem, TAttFill, TLegend, TRatioPlot, TPad, TStyle, THStack, TPaveText, gStyle
from helper import ReadJSON
from CMSGraphics import makeCMSCanvas, makeLegend
from ROOT import kBlue, kRed, kBlack, kWhite, kAzure, kOrange


# *****************************
# Declare all the variables

# data tree options 
ZZTree    = False
CRZLLTree = True
CRZLTree  = False
ZTree     = False

# data periods options
# period = "data2016"
period = "data2017"
# *****************************
if period == "data2016" :
    lumiText = '35.9 fb^{-1}'
if period == "data2017" :
    lumiText = '41.86 fb^{-1}'
#******************************


if(ZZTree):      treeText  = "ZZTree"
elif(CRZLLTree): treeText  = "CRZLLTree"
elif(CRZLTree):  treeText  = "CRZLTree"
elif(ZTree):     treeText  = "ZTree"
else: print ("Error: wrong option!")

#******************************




# create output directory
OutputPath = "IsoDistrib_DATAvsMC_" + str(period) + "_" + str(treeText)
gSystem.Exec("mkdir -p " + OutputPath)
print "Output directory created!"




# *** DATA ***
#read histos from data file

histoDATA_input = TFile.Open("IsoDistrib_DATA_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoDATA_input.GetName(),'...'

IsoDATA = []

IsoDATA.append(histoDATA_input.Get('ISO leading ele'))
IsoDATA.append(histoDATA_input.Get('ISO leading ele in ECAL Barrel'))
IsoDATA.append(histoDATA_input.Get('ISO leading ele in ECAL Endcap'))

IsoDATA.append(histoDATA_input.Get('ISO max ele'))
IsoDATA.append(histoDATA_input.Get('ISO max ele in ECAL Barrel'))
IsoDATA.append(histoDATA_input.Get('ISO max ele in ECAL Endcap'))

IsoDATA.append(histoDATA_input.Get('ISO leading mu'))
IsoDATA.append(histoDATA_input.Get('ISO leading mu in Muon Barrel'))
IsoDATA.append(histoDATA_input.Get('ISO leading mu in Muon Endcap'))

IsoDATA.append(histoDATA_input.Get('ISO max mu'))
IsoDATA.append(histoDATA_input.Get('ISO max mu in Muon Barrel'))
IsoDATA.append(histoDATA_input.Get('ISO max mu in Muon Endcap'))

if not ZTree :
    IsoDATA.append(histoDATA_input.Get('ISO extraEl'))
    IsoDATA.append(histoDATA_input.Get('ISO extraEl in ECAL Barrel'))
    IsoDATA.append(histoDATA_input.Get('ISO extraEl in ECAL Endcap'))

    IsoDATA.append(histoDATA_input.Get('ISO extraMu'))
    IsoDATA.append(histoDATA_input.Get('ISO extraMu in Muon Barrel'))
    IsoDATA.append(histoDATA_input.Get('ISO extraMu in Muon Endcap'))



# *** MC DY ***
#read histo from histoMC_DY.root

histoMC_input = TFile.Open("IsoDistrib_MC_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoMC_input.GetName(),'...'

IsoMC = []

IsoMC.append(histoMC_input.Get('ISO leading ele'))
IsoMC.append(histoMC_input.Get('ISO leading ele in ECAL Barrel'))
IsoMC.append(histoMC_input.Get('ISO leading ele in ECAL Endcap'))

IsoMC.append(histoMC_input.Get('ISO max ele'))
IsoMC.append(histoMC_input.Get('ISO max ele in ECAL Barrel'))
IsoMC.append(histoMC_input.Get('ISO max ele in ECAL Endcap'))

IsoMC.append(histoMC_input.Get('ISO leading mu'))
IsoMC.append(histoMC_input.Get('ISO leading mu in Muon Barrel'))
IsoMC.append(histoMC_input.Get('ISO leading mu in Muon Endcap'))

IsoMC.append(histoMC_input.Get('ISO max mu'))
IsoMC.append(histoMC_input.Get('ISO max mu in Muon Barrel'))
IsoMC.append(histoMC_input.Get('ISO max mu in Muon Endcap'))

if not ZTree :
    IsoMC.append(histoMC_input.Get('ISO extraEl'))
    IsoMC.append(histoMC_input.Get('ISO extraEl in ECAL Barrel'))
    IsoMC.append(histoMC_input.Get('ISO extraEl in ECAL Endcap'))

    IsoMC.append(histoMC_input.Get('ISO extraMu'))
    IsoMC.append(histoMC_input.Get('ISO extraMu in Muon Barrel'))
    IsoMC.append(histoMC_input.Get('ISO extraMu in Muon Endcap'))



if ZTree :
    nameList = ['ISO_leading_ele',
                'ISO_leading_ele_InECALbarrel',
                'ISO_leading_ele_InECALendcap',
                'maxISO_ele',
                'maxISO_ele_InECALbarrel',
                'maxISO_ele_InECALendcap',
                'ISO_leading_mu',
                'ISO_leading_mu_InMuonBarrel',
                'ISO_leading_mu_InMuonEndcap',
                'maxISO_mu',
                'maxISO_mu_InMuonBarrel',
                'maxISO_mu_InMuonEndcap']

else :
    nameList = ['ISO_leading_ele',
                'ISO_leading_ele_InECALbarrel',
                'ISO_leading_ele_InECALendcap',
                'maxISO_ele',
                'maxISO_ele_InECALbarrel',
                'maxISO_ele_InECALendcap',
                'ISO_leading_mu',
                'ISO_leading_mu_InMuonBarrel',
                'ISO_leading_mu_InMuonEndcap',
                'maxISO_mu',
                'maxISO_mu_InMuonBarrel',
                'maxISO_mu_InMuonEndcap',
                'ISO_extraEl',
                'ISO_extraEl_InECALbarrel',
                'ISO_extraEl_InECALendcap',
                'ISO_extraMu',
                'ISO_extraMu_InMuonBarrel',
                'ISO_extraMu_InMuonEndcap']
    

                                   
# *** do the plots ***
for i in range(len(nameList)) :

    canvas = TCanvas("canvas","canvas",800,800)

    #normalize MC to data
    norm = IsoDATA[i].Integral() / IsoMC[i].Integral()

    #DATA hist
    IsoDATA[i].SetMarkerStyle(20)
    IsoDATA[i].SetMarkerSize(0.6)

    #MC hist
    IsoMC[i].Scale(norm) #normalize MC to data
    IsoMC[i].SetFillColor(kOrange-3)
    IsoMC[i].SetLineColor(kBlack)


    #upper plot pad
    pad1 = TPad("pad1","pad1", 0, 0.3, 1, 1.0)
    pad1.Draw()
    pad1.cd()


    IsoMC[i].Draw("histo")
    IsoDATA[i].Draw("sameEP") 


    IsoMC[i].GetXaxis().SetTitle("Lepton Isolation")
    IsoMC[i].GetXaxis().SetLabelFont(43)
    IsoMC[i].GetXaxis().SetLabelSize(15)
    IsoMC[i].GetYaxis().SetTitleSize(20)
    IsoMC[i].GetYaxis().SetTitleFont(43)
    IsoMC[i].GetYaxis().SetTitleOffset(1.8)
    IsoMC[i].GetYaxis().SetLabelFont(43)
    IsoMC[i].GetYaxis().SetLabelSize(15)
    IsoMC[i].GetYaxis().SetTitle("Events")

    IsoMC[i].SetTitle("")
    
    gStyle.SetOptStat(0)

    pad1.SetLogy() 


    # legend
    legend = TLegend(0.74,0.68,0.94,0.87)
    legend.AddEntry(IsoDATA[i],"Data", "p")
    legend.AddEntry(IsoMC[i],"Drell-Yan MC","f")
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
    rp = TH1F(IsoDATA[i].Clone("rp"))
    rp.SetLineColor(kBlack)
    rp.SetMinimum(0.5)
    rp.SetMaximum(2.)
    rp.SetStats(0)
    rp.Divide(TH1F(IsoMC[i]))   #divide histo rp/MC
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
