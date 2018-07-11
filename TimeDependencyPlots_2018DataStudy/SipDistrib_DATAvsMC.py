#!/usr/bin/en math

# *******************
# usage: 
#    python SipDistrib_DATAvsMC.py
#
# ******************

import json
import ROOT, math, helper, CMSGraphics, CMS_lumi
from ROOT import TFile, TH1F, TCanvas, gSystem, TAttFill, TLegend, TRatioPlot, TPad, TStyle, THStack, TPaveText, gStyle
from helper import ReadJSON
from CMSGraphics import makeCMSCanvas, makeLegend
from ROOT import kBlue, kRed, kBlack, kWhite, kAzure


# *****************************
# Declare all the variables

# data tree options 
ZZTree   = False
CRZLTree = False
ZTree    = True

# data periods options
# period = "data2016"
# period = "data2017"
period = "data2018"
# *****************************
lumiText = '16.6 fb^{-1}'
#******************************


if(ZZTree):     treeText  = "ZZTree"
elif(CRZLTree): treeText  = "CRZLTree"
elif(ZTree):    treeText  = "ZTree"
else: print ("Error: wrong option!")

#******************************




# create output directory
OutputPath = "SipDistrib_DATAvsMC_" + str(period) + "_" + str(treeText)
gSystem.Exec("mkdir -p " + OutputPath)
print "Output directory created!"




# *** DATA ***
#read histos from data file

histoDATA_input = TFile.Open("SipDistrib_DATA_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoDATA_input.GetName(),'...'

SipDATA = []

SipDATA.append(histoDATA_input.Get('SIP leading ele'))
SipDATA.append(histoDATA_input.Get('SIP leading ele in ECAL Barrel'))
SipDATA.append(histoDATA_input.Get('SIP leading ele in ECAL Endcap'))

SipDATA.append(histoDATA_input.Get('SIP max ele'))
SipDATA.append(histoDATA_input.Get('SIP max ele in ECAL Barrel'))
SipDATA.append(histoDATA_input.Get('SIP max ele in ECAL Endcap'))

SipDATA.append(histoDATA_input.Get('SIP leading mu'))
SipDATA.append(histoDATA_input.Get('SIP leading mu in Muon Barrel'))
SipDATA.append(histoDATA_input.Get('SIP leading mu in Muon Endcap'))

SipDATA.append(histoDATA_input.Get('SIP max mu'))
SipDATA.append(histoDATA_input.Get('SIP max mu in Muon Barrel'))
SipDATA.append(histoDATA_input.Get('SIP max mu in Muon Endcap'))

if not ZTree :
    SipDATA.append(histoDATA_input.Get('SIP extraEl'))
    SipDATA.append(histoDATA_input.Get('SIP extraEl in ECAL Barrel'))
    SipDATA.append(histoDATA_input.Get('SIP extraEl in ECAL Endcap'))

    SipDATA.append(histoDATA_input.Get('SIP extraMu'))
    SipDATA.append(histoDATA_input.Get('SIP extraMu in Muon Barrel'))
    SipDATA.append(histoDATA_input.Get('SIP extraMu in Muon Endcap'))



# *** MC DY ***
#read histo from histoMC_DY.root

histoMC_input = TFile.Open("SipDistrib_MC_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoMC_input.GetName(),'...'

SipMC = []

SipMC.append(histoMC_input.Get('SIP leading ele'))
SipMC.append(histoMC_input.Get('SIP leading ele in ECAL Barrel'))
SipMC.append(histoMC_input.Get('SIP leading ele in ECAL Endcap'))

SipMC.append(histoMC_input.Get('SIP max ele'))
SipMC.append(histoMC_input.Get('SIP max ele in ECAL Barrel'))
SipMC.append(histoMC_input.Get('SIP max ele in ECAL Endcap'))

SipMC.append(histoMC_input.Get('SIP leading mu'))
SipMC.append(histoMC_input.Get('SIP leading mu in Muon Barrel'))
SipMC.append(histoMC_input.Get('SIP leading mu in Muon Endcap'))

SipMC.append(histoMC_input.Get('SIP max mu'))
SipMC.append(histoMC_input.Get('SIP max mu in Muon Barrel'))
SipMC.append(histoMC_input.Get('SIP max mu in Muon Endcap'))

if not ZTree :
    SipMC.append(histoMC_input.Get('SIP extraEl'))
    SipMC.append(histoMC_input.Get('SIP extraEl in ECAL Barrel'))
    SipMC.append(histoMC_input.Get('SIP extraEl in ECAL Endcap'))

    SipMC.append(histoMC_input.Get('SIP extraMu'))
    SipMC.append(histoMC_input.Get('SIP extraMu in Muon Barrel'))
    SipMC.append(histoMC_input.Get('SIP extraMu in Muon Endcap'))



if ZTree :
    nameList = ['SIP_leading_ele',
                'SIP_leading_ele_InECALbarrel',
                'SIP_leading_ele_InECALendcap',
                'maxSIP_ele',
                'maxSIP_ele_InECALbarrel',
                'maxSIP_ele_InECALendcap',
                'SIP_leading_mu',
                'SIP_leading_mu_InMuonBarrel',
                'SIP_leading_mu_InMuonEndcap',
                'maxSIP_mu',
                'maxSIP_mu_InMuonBarrel',
                'maxSIP_mu_InMuonEndcap']

else : 
    nameList = ['SIP_leading_ele',
                'SIP_leading_ele_InECALbarrel',
                'SIP_leading_ele_InECALendcap',
                'maxSIP_ele',
                'maxSIP_ele_InECALbarrel',
                'maxSIP_ele_InECALendcap',
                'SIP_leading_mu',
                'SIP_leading_mu_InMuonBarrel',
                'SIP_leading_mu_InMuonEndcap',
                'maxSIP_mu',
                'maxSIP_mu_InMuonBarrel',
                'maxSIP_mu_InMuonEndcap',
                'SIP_extraEl',
                'SIP_extraEl_InECALbarrel',
                'SIP_extraEl_InECALendcap',
                'SIP_extraMu',
                'SIP_extraMu_InMuonBarrel',
                'SIP_extraMu_InMuonEndcap']
        

                                   
# *** do the plots ***
for i in range(len(nameList)) :

    canvas = TCanvas("canvas","canvas",800,800)

    #normalize MC 
    # norm = 1 # to lumi 
    norm = SipDATA[i].Integral() / SipMC[i].Integral() #to data

    #DATA hist
    SipDATA[i].SetMarkerStyle(20)
    SipDATA[i].SetMarkerSize(0.6)

    #MC hist
    SipMC[i].Scale(norm) #normalize MC 
    SipMC[i].SetFillColor(kAzure+6)
    SipMC[i].SetLineColor(kBlack)


    #upper plot pad
    pad1 = TPad("pad1","pad1", 0, 0.3, 1, 1.0)
    pad1.Draw()
    pad1.cd()


    SipMC[i].Draw("histo")
    SipDATA[i].Draw("sameEP") 


    SipMC[i].GetXaxis().SetTitle("SIP")
    SipMC[i].GetXaxis().SetLabelFont(43)
    SipMC[i].GetXaxis().SetLabelSize(15)
    SipMC[i].GetYaxis().SetTitleSize(20)
    SipMC[i].GetYaxis().SetTitleFont(43)
    SipMC[i].GetYaxis().SetTitleOffset(1.8)
    SipMC[i].GetYaxis().SetLabelFont(43)
    SipMC[i].GetYaxis().SetLabelSize(15)
    SipMC[i].GetYaxis().SetTitle("Events")

    SipMC[i].SetTitle("")
    
    gStyle.SetOptStat(0)


    # legend
    legend = TLegend(0.74,0.68,0.94,0.87)
    legend.AddEntry(SipDATA[i],"Data", "p")
    legend.AddEntry(SipMC[i],"Drell-Yan MC","f")
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
    rp = TH1F(SipDATA[i].Clone("rp"))
    rp.SetLineColor(kBlack)
    rp.SetMinimum(0.5)
    rp.SetMaximum(2.)
    rp.SetStats(0)
    rp.Divide(TH1F(SipMC[i]))   #divide histo rp/MC
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
