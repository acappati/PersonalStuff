#!/usr/bin/en math

# *******************
# usage: 
#    python yellowPlots.py
#
# structure: 
#    - read json files which contain DATA and MC DCB fit results obtained by running FitDATA.py and FitMC.py
#    - read root files where data and MC histos are stored 
#    - read 3 files: 2 MC and 1 DATA file (root files created by FitMC.py and FitDATA.py)
#    - append histos in 3 lists
#    - do plots with 2 pads: 1) with data vs MC histo, 2) with data/MC ratio  
#    - save plots in OutputPath
# ********************

import json
import ROOT, math, helper, CMSGraphics, CMS_lumi
from ROOT import TFile, TH1F, TCanvas, gSystem, TAttFill, TLegend, TRatioPlot, TPad, TStyle, THStack, TPaveText
from helper import ReadJSON
from CMSGraphics import makeCMSCanvas, makeLegend
from ROOT import kBlue, kRed, kBlack, kWhite



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
lumiText = '41.86 fb^{-1}'
#******************************


if(ZZTree):     treeText  = "ZZTree"
elif(CRZLTree): treeText  = "CRZLTree"
elif(ZTree):    treeText  = "ZTree"
else: print ("Error: wrong option!")


#******************************************************************************
# read fit values from output files obtained by running FitDATA.py and FitMC.py

with open("out_ZDBCmean_DATA_" + period + "_" + treeText  + ".json","r") as handle1 :
    massFitDATA_dict = json.load(handle1)

with open("out_ZDBCwidth_DATA_" + period + "_" + treeText  + ".json","r") as handle2 :
    widthFitDATA_dict = json.load(handle2)

with open("out_ZDBCmean_MC_DYJets_" + period + "_" + treeText  + ".json","r") as handle3 :
    massFitMC_dict = json.load(handle3)

with open("out_ZDBCwidth_MC_DYJets_" + period + "_" + treeText  + ".json","r") as handle4 :
    widthFitMC_dict = json.load(handle4)

print 'Input files with Fit values read!'
#****************************************************************************




# create output directory
OutputPath = "yellowPlots_" + str(period) + "_" + str(treeText)
gSystem.Exec("mkdir -p " + OutputPath)
print "Output directory created!"


# *** DATA ***
#read histos from data file
if(period == "data2016"):
    histoDATA_input = TFile.Open("histoDATA_data2016_" + str(treeText) + ".root")
elif(period == "data2017"):
    histoDATA_input = TFile.Open("histoDATA_data2017_" + str(treeText) + ".root")
else: 
    print "Wrong period!"
print 'Reading file', histoDATA_input.GetName(),'...'

ZMass_DATA = []

ZMass_DATA.append(histoDATA_input.Get('ZMass_ele'))              #ZMass , Z->ee    
if not ZTree :            
    ZMass_DATA.append(histoDATA_input.Get('ZMass_ele_extraMu'))  #ZMass , Z->ee + mu           
    ZMass_DATA.append(histoDATA_input.Get('ZMass_ele_extraEl'))  #ZMass , Z->ee + e            
ZMass_DATA.append(histoDATA_input.Get('ZMass_ele_EBEB'))         #ZMass , Z->ee Barrel-Barrel  
ZMass_DATA.append(histoDATA_input.Get('ZMass_ele_EBEE'))         #ZMass , Z->ee Barrel-Endcap  
ZMass_DATA.append(histoDATA_input.Get('ZMass_ele_EEEE'))         #ZMass , Z->ee Endcap-Endcap  
ZMass_DATA.append(histoDATA_input.Get('ZMass_mu'))               #ZMass , Z->mumu        
if not ZTree :   
    ZMass_DATA.append(histoDATA_input.Get('ZMass_mu_extraMu'))   #ZMass , Z->mumu + mu         
    ZMass_DATA.append(histoDATA_input.Get('ZMass_mu_extraEl'))   #ZMass , Z->mumu + e          
ZMass_DATA.append(histoDATA_input.Get('ZMass_mu_MBMB'))          #ZMass , Z->mumu Barrel-Barrel
ZMass_DATA.append(histoDATA_input.Get('ZMass_mu_MBME'))          #ZMass , Z->mumu Barrel-Endcap
ZMass_DATA.append(histoDATA_input.Get('ZMass_mu_MEME'))          #ZMass , Z->mumu Endcap-Endcap



# *** MC DY ***
#read histo from histoMC_DY.root
if(period == "data2016"):
    histoMC_DY_input = TFile.Open("histoMC_DYJets_data2016_" + str(treeText) + ".root")
elif(period == "data2017"):
    histoMC_DY_input = TFile.Open("histoMC_DYJets_data2017_" + str(treeText) + ".root")
else: 
    print "Wrong period!"
print 'Reading file', histoMC_DY_input.GetName(),'...'

ZMass_MC_DY = []

ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_ele'))              #ZMass , Z->ee
if not ZTree :
    ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_ele_extraMu'))  #ZMass , Z->ee + mu
    ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_ele_extraEl'))  #ZMass , Z->ee + e
ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_ele_EBEB'))         #ZMass , Z->ee Barrel-Barrel
ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_ele_EBEE'))         #ZMass , Z->ee Barrel-Endcap
ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_ele_EEEE'))         #ZMass , Z->ee Endcap-Endcap
ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_mu'))               #ZMass , Z->mumu
if not ZTree :
    ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_mu_extraMu'))   #ZMass , Z->mumu + mu
    ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_mu_extraEl'))   #ZMass , Z->mumu + e
ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_mu_MBMB'))          #ZMass , Z->mumu Barrel-Barrel
ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_mu_MBME'))          #ZMass , Z->mumu Barrel-Endcap
ZMass_MC_DY.append(histoMC_DY_input.Get('ZMass_mu_MEME'))          #ZMass , Z->mumu Endcap-Endcap


# *** MC TTJets ***
#read histo from histoMC_TTJets.root
if(period == "data2016"):
    histoMC_TTJets_input = TFile.Open("histoMC_TTJets_data2016_" + str(treeText) +".root")
elif(period == "data2017"):
    histoMC_TTJets_input = TFile.Open("histoMC_TTJets_data2017_" + str(treeText) +".root")
else: 
    print "Wrong period!"
print 'Reading file', histoMC_TTJets_input.GetName(),'...'

ZMass_MC_TTJets = []

ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_ele'))              #ZMass , Z->ee
if not ZTree :
    ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_ele_extraMu'))  #ZMass , Z->ee + mu
    ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_ele_extraEl'))  #ZMass , Z->ee + e
ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_ele_EBEB'))         #ZMass , Z->ee Barrel-Barrel
ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_ele_EBEE'))         #ZMass , Z->ee Barrel-Endcap
ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_ele_EEEE'))         #ZMass , Z->ee Endcap-Endcap
ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_mu'))               #ZMass , Z->mumu
if not ZTree :
    ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_mu_extraMu'))   #ZMass , Z->mumu + mu
    ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_mu_extraEl'))   #ZMass , Z->mumu + e
ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_mu_MBMB'))          #ZMass , Z->mumu Barrel-Barrel
ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_mu_MBME'))          #ZMass , Z->mumu Barrel-Endcap
ZMass_MC_TTJets.append(histoMC_TTJets_input.Get('ZMass_mu_MEME'))          #ZMass , Z->mumu Endcap-Endcap


if ZTree :
    nameList = ['ZMass_ele', 'ZMass_ele_EBEB', 'ZMass_ele_EBEE', 'ZMass_ele_EEEE', 'ZMass_mu', 'ZMass_mu_MBMB', 'ZMass_mu_MBME', 'ZMass_mu_MEME']
else :
    nameList = ['ZMass_ele', 'ZMass_ele_extraMu', 'ZMass_ele_extraEl', 'ZMass_ele_EBEB', 'ZMass_ele_EBEE', 'ZMass_ele_EEEE', 'ZMass_mu', 'ZMass_mu_extraMu', 'ZMass_mu_extraEl', 'ZMass_mu_MBMB', 'ZMass_mu_MBME', 'ZMass_mu_MEME']



# list of dictionary keys to loop over dictionaries in the correct order
keyList = ['Zee','Zee_extraMu','Zee_extraEl','Zee_EBEB','Zee_EBEE','Zee_EEEE','Zmumu','Zmumu_extraMu','Zmumu_extraEl','Zmumu_MBMB','Zmumu_MBME','Zmumu_MEME']


# *** do plots *** 
i = 0  #counter for histos lists

for key in keyList :

    # debug
    # print key, massFitDATA_dict[key]
    # print key, widthFitDATA_dict[key]                  
    # print key, massFitMC_dict[key]
    # print key, widthFitMC_dict[key]


    if (ZTree and massFitDATA_dict[key]==0) : continue

    canvas = TCanvas("canvas","canvas",800,800)
    
    hs = THStack("hs","")

    # norm = 1 # normalize to MC xsection 
    norm = ZMass_DATA[i].Integral() / (ZMass_MC_TTJets[i].Integral() + ZMass_MC_DY[i].Integral()) #normalize MC to data


    #DATA hist
    ZMass_DATA[i].SetMarkerStyle(20)
    ZMass_DATA[i].SetMarkerSize(0.6)
            
    #TTJets MC hist 
    ZMass_MC_TTJets[i].Scale(norm) #normalize MC to data
    ZMass_MC_TTJets[i].SetFillColor(3)
    ZMass_MC_TTJets[i].SetLineColor(kBlack)
    hs.Add(ZMass_MC_TTJets[i])

    # DY MC hist
    ZMass_MC_DY[i].Scale(norm) #normalize MC to data
    ZMass_MC_DY[i].SetFillColor(5)
    ZMass_MC_DY[i].SetLineColor(kBlack)
    hs.Add(ZMass_MC_DY[i])

    
    #upper plot pad
    pad1 = TPad("pad1","pad1", 0, 0.3, 1, 1.0)
    pad1.Draw()
    pad1.cd()

    
    hs.Draw("histo")
    ZMass_DATA[i].Draw("sameEP")    


    hs.GetXaxis().SetTitle("Mass [GeV/c^{2}]")
    hs.GetXaxis().SetLabelFont(43)
    hs.GetXaxis().SetLabelSize(15)
    hs.GetYaxis().SetTitleSize(20)
    hs.GetYaxis().SetTitleFont(43)
    hs.GetYaxis().SetTitleOffset(1.8)
    hs.GetYaxis().SetLabelFont(43)
    hs.GetYaxis().SetLabelSize(15)
    hs.GetYaxis().SetTitle("Events")

    
    # legend
    legend = TLegend(0.74,0.68,0.94,0.87)
    legend.AddEntry(ZMass_DATA[i],"Data", "p")
    legend.AddEntry(ZMass_MC_DY[i],"Drell-Yan","f")
    legend.AddEntry(ZMass_MC_TTJets[i],"t#bar{t}","f")
    legend.SetFillColor(kWhite)
    legend.SetLineColor(kBlack)
    legend.SetTextFont(43)
    legend.SetTextSize(20)
    legend.Draw()

    
    # box with fit results
    pv = TPaveText(0.64,0.35,0.94,0.65,"brNDC")
    pv.AddText("DATA:")
    pv.AddText("Z DCBmean = "  + str(round(massFitDATA_dict[key],2))  + " GeV")
    pv.AddText("Z DCBwidth = " + str(round(widthFitDATA_dict[key],2)) + " GeV")
    pv.AddText("MC:")         
    pv.AddText("Z DCBmean = "  + str(round(massFitMC_dict[key],2))  + " GeV")
    pv.AddText("Z DCBwidth = " + str(round(widthFitMC_dict[key],2)) + " GeV")
    pv.SetFillColor(kWhite)
    pv.SetBorderSize(1)
    pv.SetTextFont(40)
    pv.SetTextSize(0.037)
    pv.SetTextFont(42)
    pv.SetTextAlign(12) #text left aligned 
    # cange color of text 
    if "ele" in nameList[i] : 
        pv.SetTextColor(kBlue)
    elif "mu" in nameList[i] :
        pv.SetTextColor(kRed)
    else :
        pv.SetTextColor(kBlack)
    pv.Draw()


    # box with decay mode
    pv2 = TPaveText(0.19,0.75,0.35,0.85,"brNDC")
    if "ele" in nameList[i] : 
        pv2.AddText("Z #rightarrow e^{+} e^{-}")
        pv2.SetTextColor(kBlue)
    elif "mu" in nameList[i] :
        pv2.AddText("Z #rightarrow #mu^{+} #mu^{-}")
        pv2.SetTextColor(kRed)
    pv2.SetFillColor(kWhite)
    pv2.SetBorderSize(1)
    pv2.SetTextFont(40)
    pv2.SetTextSize(0.05)
    pv2.SetTextFont(42)
    pv2.SetTextAlign(22) #text centering
    pv2.Draw()

    
    canvas.Update()
    

    #lower plot pad
    canvas.cd()
    pad2 = TPad("pad2","pad2", 0, 0.05, 1, 0.3)
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()    #pad2 becomes the current pad

    #define ratio plot
    rp = TH1F(ZMass_DATA[i].Clone("rp"))
    rp.SetLineColor(kBlack)
    rp.SetMinimum(0.5)
    rp.SetMaximum(2.)
    rp.SetStats(0)
    rp.Divide(TH1F(ZMass_MC_DY[i]+ZMass_MC_TTJets[i]))   #divide histo rp/MC
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


    i = i +1   #counter for histos lists


print "plots done"



