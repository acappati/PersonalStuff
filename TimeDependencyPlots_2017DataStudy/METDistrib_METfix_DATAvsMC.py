#!/usr/bin/env python

# *******************
# usage: 
#    python METDistrib_METfix_DATAvsMC.py
#
# *******************


import json
import ROOT, math, helper, CMSGraphics, CMS_lumi
from ROOT import TFile, TH1, TH1F, TCanvas, gSystem, TAttFill, TLegend, TRatioPlot, TPad, TStyle, THStack, TPaveText, gStyle
from CMSGraphics import makeCMSCanvas, makeLegend
from helper import ReadJSON
from helper import DoSimpleFit, Result
from ROOT import kBlue, kRed, kBlack, kWhite, kAzure, kOrange


# *****************************
# Declare all the variables
# options
redoDATAHistos    = True
redoMCDYHistos    = True
redoMCTTbarHistos = True

# data tree options 
ZZTree    = False
CRZLLTree = False
CRZLTree  = False
ZTree     = True

# data periods options
#period = "data2017_preMETfix"
period = "data2017_withMETfix"
# *****************************


if(period == "data2017_preMETfix"):
    inputDATAtree    = TFile.Open("/data3/Higgs/180531_2017/AllData/ZZ4lAnalysis.root")        #2017 data (rereco json)   
    inputMCDYtree    = TFile.Open("/data3/Higgs/180531_2017/DYJetsToLL_M50/ZZ4lAnalysis.root") #DYJets 2017 MC
    inputMCTTbartree = TFile.Open("/data3/Higgs/180531_2017/TTTo2L2Nu/ZZ4lAnalysis.root")      #TTbarJets 2017 MC
    lumi     = 41.30   # fb-1
    lumiText = '41.30 fb^{-1}'
    if(ZZTree):
        treeDATA    = inputDATAtree.Get("ZZTree/candTree")
        treeMCDY    = inputMCDYtree.Get("ZZTree/candTree")
        treeMCTTbar = inputMCTTbartree.Get("ZZTree/candTree")
        treeText  = "ZZTree"
    elif(CRZLLTree):
        treeDATA    = inputDATAtree.Get("CRZLLTree/candTree")
        treeMCDY    = inputMCDYtree.Get("CRZLLTree/candTree")
        treeMCTTbar = inputMCTTbartree.Get("CRZLLTree/candTree")
        treeText  = "CRZLLTree"
    elif(CRZLTree):
        treeDATA    = inputDATAtree.Get("CRZLTree/candTree")
        treeMCDY    = inputMCDYtree.Get("CRZLTree/candTree")
        treeMCTTbar = inputMCTTbartree.Get("CRZLTree/candTree")
        treeText  = "CRZLTree"
    elif(ZTree):
        treeDATA    = inputDATAtree.Get("ZTree/candTree")
        treeMCDY    = inputMCDYtree.Get("ZTree/candTree")
        treeMCTTbar = inputMCTTbartree.Get("ZTree/candTree")
        treeText  = "ZTree"
    else:
        print ("Error: wrong option!")

elif(period == "data2017_withMETfix"):
    inputDATAtree    = TFile.Open("/data3/Higgs/181015_2017_MET/AllData/ZZ4lAnalysis.root")        #2017 data (rereco json)   
    inputMCDYtree    = TFile.Open("/data3/Higgs/181015_2017_MET/DYJetsToLL_M50/ZZ4lAnalysis.root") #DYJets 2017 MC 
    inputMCTTbartree = TFile.Open("/data3/Higgs/181015_2017_MET/TTTo2L2Nu/ZZ4lAnalysis.root")      #TTbarJets 2017 MC
    lumi     = 41.30   # fb-1
    lumiText = '41.30 fb^{-1}'
    if(ZZTree):
        treeDATA    = inputDATAtree.Get("ZZTree/candTree")
        treeMCDY    = inputMCDYtree.Get("ZZTree/candTree")
        treeMCTTbar = inputMCTTbartree.Get("ZZTree/candTree")
        treeText  = "ZZTree"
    elif(CRZLLTree):
        treeDATA    = inputDATAtree.Get("CRZLLTree/candTree")
        treeMCDY    = inputMCDYtree.Get("CRZLLTree/candTree")
        treeMCTTbar = inputMCTTbartree.Get("CRZLLTree/candTree")
        treeText  = "CRZLLTree"
    elif(CRZLTree):
        treeDATA    = inputDATAtree.Get("CRZLTree/candTree")
        treeMCDY    = inputMCDYtree.Get("CRZLTree/candTree")
        treeMCTTbar = inputMCTTbartree.Get("CRZLTree/candTree")
        treeText  = "CRZLTree"
    elif(ZTree):
        treeDATA    = inputDATAtree.Get("ZTree/candTree")
        treeMCDY    = inputMCDYtree.Get("ZTree/candTree")
        treeMCTTbar = inputMCTTbartree.Get("ZTree/candTree")
        treeText  = "ZTree"
    else:
        print ("Error: wrong option!")
else: 
    print ("Error: choose a period!")


# ********************
#  do data histos 
# ********************
if(redoDATAHistos) :

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    

    PFMET_hist_ele = TH1F('PFMET_ele','PFMET_ele',125,0,250) #Z->ee
    PFMET_hist_mu = TH1F('PFMET_mu','PFMET_mu',125,0,250)    #Z->mumu

    

    # read tree 
    print "reading tree", inputDATAtree.GetName(),treeText,treeDATA.GetName()  ,"..."
    
    treeDATA.SetBranchStatus("*",0)  # disable all branches
    if ZTree :
        treeDATA.SetBranchStatus("Zsel",1)
    else :
        treeDATA.SetBranchStatus("ZZsel",1)
    treeDATA.SetBranchStatus("LepLepId",1)
    treeDATA.SetBranchStatus("PFMET",1)
   
    

    for event in treeDATA:
        if ZTree :
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger
        else :
            if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger
  

        # Z->ee histos
        if(int(math.fabs(event.LepLepId[0])) == 11 ):
            
            PFMET_hist_ele.Fill(event.PFMET)


        # Z->mumu histos
        if(int(math.fabs(event.LepLepId[0])) == 13 ):    
            
            PFMET_hist_mu.Fill(event.PFMET)

        
    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_DATA = TFile.Open("METDistrib_DATA_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_DATA.cd()

    PFMET_hist_ele.Write()
    PFMET_hist_mu.Write()


    outFile_DATA.Close()
    print "DATA histo file created!"


# ********************
#  do MC DY histos
# ********************
if(redoMCDYHistos) :

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    PFMET_hist_MC_DY_ele = TH1F('PFMET_MC_DY_ele','PFMET_MC_DY_ele',125,0,250) #Z->ee
    PFMET_hist_MC_DY_mu = TH1F('PFMET_MC_DY_mu','PFMET_MC_DY_mu',125,0,250)    #Z->mumu


    # get partial event weight
    hcounters           = inputMCDYtree.Get("ZZTree/Counters")
    gen_sumWeights      = hcounters.GetBinContent(40)
    partialSampleWeight = lumi * 1000 / gen_sumWeights


    # read tree 
    print "reading tree", inputMCDYtree.GetName(),treeText,treeMCDY.GetName()  ,"..."
    for event in treeMCDY:
        if ZTree :
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger
        else :
            if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger


        weight = partialSampleWeight*event.xsec*event.overallEventWeight


        # Z->ee histos
        if(int(math.fabs(event.LepLepId[0])) == 11 ):           
            
            PFMET_hist_MC_DY_ele.Fill(event.PFMET,weight)

        # Z->mumu histos
        if(int(math.fabs(event.LepLepId[0])) == 13 ):

            PFMET_hist_MC_DY_mu.Fill(event.PFMET,weight)    

        
    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_MCDY = TFile.Open("METDistrib_MC_DY_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_MCDY.cd()

    PFMET_hist_MC_DY_ele.Write()
    PFMET_hist_MC_DY_mu.Write()


    outFile_MCDY.Close()
    print "MC DY histo file created!"
# ********************


# ********************
#  do MC TTbar histos
# ********************
if(redoMCTTbarHistos) :

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    PFMET_hist_MC_TTbar_ele = TH1F('PFMET_MC_TTbar_ele','PFMET_MC_TTbar_ele',125,0,250) #Z->ee
    PFMET_hist_MC_TTbar_mu = TH1F('PFMET_MC_TTbar_mu','PFMET_MC_TTbar_mu',125,0,250)    #Z->mumu


    # get partial event weight
    hcounters           = inputMCTTbartree.Get("ZZTree/Counters")
    gen_sumWeights      = hcounters.GetBinContent(40)
    partialSampleWeight = lumi * 1000 / gen_sumWeights


    # read tree 
    print "reading tree", inputMCTTbartree.GetName(),treeText,treeMCTTbar.GetName()  ,"..."
    for event in treeMCTTbar:
        if ZTree :
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger
        else :
            if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger


        weight = partialSampleWeight*event.xsec*event.overallEventWeight


        # Z->ee histos
        if(int(math.fabs(event.LepLepId[0])) == 11 ):

            PFMET_hist_MC_TTbar_ele.Fill(event.PFMET,weight)

        # Z->mumu histos
        if(int(math.fabs(event.LepLepId[0])) == 13 ):
            
            PFMET_hist_MC_TTbar_mu.Fill(event.PFMET,weight)


        
    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_MCTTbar = TFile.Open("METDistrib_MC_TTbar_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_MCTTbar.cd()

    PFMET_hist_MC_TTbar_ele.Write()
    PFMET_hist_MC_TTbar_mu.Write()


    outFile_MCTTbar.Close()
    print "MC TTbar histo file created!"
# ********************



# ****************************************
# create output directory 
outputDir = "METDistrib_DATAvsMC_" + str(period) + "_" + str(treeText)
gSystem.Exec("mkdir -p " + outputDir)
print "Output directory created!"


# **************************
# read data histos from file 
histoDATA_input = TFile.Open("METDistrib_DATA_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoDATA_input.GetName(),'...'

inDATA_list = []
inDATA_list.append(histoDATA_input.Get('PFMET_ele'))
inDATA_list.append(histoDATA_input.Get('PFMET_mu'))


# ****************************
# read DY MC histos from file 
histoMCDY_input = TFile.Open("METDistrib_MC_DY_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoMCDY_input.GetName(),'...'

inMCDY_list = []
inMCDY_list.append(histoMCDY_input.Get('PFMET_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('PFMET_MC_DY_mu'))



# ****************************
# read TTbar MC histos from file 
histoMCTTbar_input = TFile.Open("METDistrib_MC_TTbar_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoMCTTbar_input.GetName(),'...'

inMCTTbar_list = []
inMCTTbar_list.append(histoMCTTbar_input.Get('PFMET_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('PFMET_MC_TTbar_mu'))

       

# ******************************
# do DATA vs MC comparison plots  
for i in range(len(inDATA_list)) : 

    canvas = TCanvas("canvas","canvas",800,800)

    hs = THStack("hs","")

    norm = 1 # normalize to MC xsection 
    # norm = inDATA_list[i].Integral() / (inMCTTbar_list[i].Integral() + inMCDY_list[i].Integral()) #normalize MC to data


    #DATA hist
    inDATA_list[i].SetMarkerStyle(20)
    inDATA_list[i].SetMarkerSize(0.6)


    #MC TTbar hist
    inMCTTbar_list[i].Scale(norm) #normalize MC 
    inMCTTbar_list[i].SetFillColor(kAzure-2)
    inMCTTbar_list[i].SetLineColor(kBlack)
    hs.Add(inMCTTbar_list[i])

    #MC DY hist
    inMCDY_list[i].Scale(norm) #normalize MC
    inMCDY_list[i].SetFillColor(kOrange-3)
    inMCDY_list[i].SetLineColor(kBlack)
    hs.Add(inMCDY_list[i])


    #upper plot pad
    pad1 = TPad("pad1","pad1", 0, 0.3, 1, 1.0)
    pad1.Draw()
    pad1.cd()


    hs.SetMaximum(1.3*max(hs.GetMaximum(),inDATA_list[i].GetMaximum()))
    inDATA_list[i].SetMaximum(1.3*max(hs.GetMaximum(),inDATA_list[i].GetMaximum()))
    

    hs.Draw("histo") 
    inDATA_list[i].Draw("sameEP")
    
    
    hs.SetTitle("")
    hs.GetXaxis().SetTitle(inDATA_list[i].GetTitle())
    hs.GetXaxis().SetLabelFont(43)
    hs.GetXaxis().SetLabelSize(15)
    hs.GetYaxis().SetTitleSize(20)
    hs.GetYaxis().SetTitleFont(43)
    hs.GetYaxis().SetTitleOffset(1.8)
    hs.GetYaxis().SetLabelFont(43)
    hs.GetYaxis().SetLabelSize(15)
    hs.GetYaxis().SetTitle("Events")

    gStyle.SetOptStat(0)

    if "Pt" in inDATA_list[i].GetTitle() :
        pad1.SetLogy()


    # legend
    legend = TLegend(0.82,0.75,0.95,0.89)
    legend.AddEntry(inDATA_list[i],"Data", "p")
    legend.AddEntry(inMCDY_list[i],"DY MC","f")
    legend.AddEntry(inMCTTbar_list[i],"t#bar{t} MC","f")
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
    rp = TH1F(inDATA_list[i].Clone("rp"))
    rp.SetLineColor(kBlack)
    rp.SetMinimum(0.5)
    rp.SetMaximum(2.)
    rp.SetStats(0)
    rp.Divide(TH1F(inMCDY_list[i]+inMCTTbar_list[i]))   #divide histo rp/MC
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


    canvas.SaveAs(outputDir + "/" + inDATA_list[i].GetTitle() + ".pdf")
    canvas.SaveAs(outputDir + "/" + inDATA_list[i].GetTitle() + ".png")


print "plots done"
