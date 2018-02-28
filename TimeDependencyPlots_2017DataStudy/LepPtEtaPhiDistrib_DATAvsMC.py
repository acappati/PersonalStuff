#!/usr/bin/env python

# *******************
# usage: 
#    python LepPtEtaPhiDistrib_DATAvsMC.py
#
# *******************


import json
import ROOT, math, helper, CMSGraphics, CMS_lumi
from ROOT import TFile, TH1, TH1F, TCanvas, gSystem, TAttFill, TLegend, TRatioPlot, TPad, TStyle, THStack, TPaveText, gStyle
from CMSGraphics import makeCMSCanvas, makeLegend
from helper import ReadJSON
from helper import DoSimpleFit, Result
from ROOT import kBlue, kRed, kBlack, kWhite, kAzure, kOrange, kGreen



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
# period = "data2016"
period = "data2017"
# *****************************


#input file
if(period == "data2016"):
    inputDATAtree    = TFile.Open("/data3/Higgs/170222/AllData/ZZ4lAnalysis.root")        #2016 data
    inputMCDYtree    = TFile.Open("/data3/Higgs/170222/DYJetsToLL_M50/ZZ4lAnalysis.root") #DYJets 2016 MC
    inputMCTTbartree = TFile.Open("/data3/Higgs/170222/TTJets_DiLept/ZZ4lAnalysis.root")  #TTbarJets 2016 MC
    lumi     = 35.9  # fb-1
    lumiText = '35.9 fb^{-1}'
    if(ZZTree):
        treeDATA    = inputDATAtree.Get("ZZTree/candTree")
        treeMCDY    = inputMCDYtree.Get("ZZTree/candTree")
        treeMCTTbar = inputMCTTbartree.Get("ZZTree/candTree")
        treeText  = "ZZTree"
    elif(CRZLLTree):
        treeDATA    = inputDATAtree.Get("CRZLLTree/candTree")
        treeMCDY    = inputMCDYTtee.Get("CRZLLTree/candTree")
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

elif(period == "data2017"):
    inputDATAtree    = TFile.Open("/data3/Higgs/180218/AllData/ZZ4lAnalysis.root")        #2017 data 
    inputMCDYtree    = TFile.Open("/data3/Higgs/180218/DYJetsToLL_M50/ZZ4lAnalysis.root") #DYJets 2017 MC 
    inputMCTTbartree = TFile.Open("/data3/Higgs/180218/TTTo2L2Nu/ZZ4lAnalysis.root")      #TTbarJets 2017 MC 
    lumi     = 41.37   # fb-1
    lumiText = '41.37 fb^{-1}'
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

    TH1.SetDefaultSumw2() # set sumw2 = true for all the histograms created from now on

    # define data histograms Z->ee
    LepPt_hist_ele_1stlep        = TH1F('elePt_leadingEle',           'elePt_leadingEle',           100, 0, 300)
    LepPt_hist_ele_1stlep_barrel = TH1F('elePt_leadingEle_barrel',    'elePt_leadingEle_barrel',    100, 0, 300)
    LepPt_hist_ele_1stlep_endcap = TH1F('elePt_leadingEle_endcap',    'elePt_leadingEle_endcap',    100, 0, 300)
    LepPt_hist_ele_2ndlep        = TH1F('elePt_subleadingEle',        'elePt_subleadingEle',        100, 0, 300)
    LepPt_hist_ele_2ndlep_barrel = TH1F('elePt_subleadingEle_barrel', 'elePt_subleadingEle_barrel', 100, 0, 300)
    LepPt_hist_ele_2ndlep_endcap = TH1F('elePt_subleadingEle_endcap', 'elePt_subleadingEle_endcap', 100, 0, 300)

    LepEta_hist_ele_1stlep        = TH1F('eleEta_leadingEle',           'eleEta_leadingEle',           40, -3, 3)
    LepEta_hist_ele_1stlep_barrel = TH1F('eleEta_leadingEle_barrel',    'eleEta_leadingEle_barrel',    40, -3, 3)
    LepEta_hist_ele_1stlep_endcap = TH1F('eleEta_leadingEle_endcap',    'eleEta_leadingEle_endcap',    40, -3, 3)
    LepEta_hist_ele_2ndlep        = TH1F('eleEta_subleadingEle',        'eleEta_subleadingEle',        40, -3, 3)
    LepEta_hist_ele_2ndlep_barrel = TH1F('eleEta_subleadingEle_barrel', 'eleEta_subleadingEle_barrel', 40, -3, 3)
    LepEta_hist_ele_2ndlep_endcap = TH1F('eleEta_subleadingEle_endcap', 'eleEta_subleadingEle_endcap', 40, -3, 3)

    LepPhi_hist_ele_1stlep        = TH1F('elePhi_leadingEle',           'elePhi_leadingEle',           35, -3.3, 3.3)
    LepPhi_hist_ele_1stlep_barrel = TH1F('elePhi_leadingEle_barrel',    'elePhi_leadingEle_barrel',    35, -3.3, 3.3)
    LepPhi_hist_ele_1stlep_endcap = TH1F('elePhi_leadingEle_endcap',    'elePhi_leadingEle_endcap',    35, -3.3, 3.3)
    LepPhi_hist_ele_2ndlep        = TH1F('elePhi_subleadingEle',        'elePhi_subleadingEle',        35, -3.3, 3.3)
    LepPhi_hist_ele_2ndlep_barrel = TH1F('elePhi_subleadingEle_barrel', 'elePhi_subleadingEle_barrel', 35, -3.3, 3.3)
    LepPhi_hist_ele_2ndlep_endcap = TH1F('elePhi_subleadingEle_endcap', 'elePhi_subleadingEle_endcap', 35, -3.3, 3.3)


    # define data histograms Z->mumu
    LepPt_hist_mu_1stlep        = TH1F('muPt_leadingMu',           'muPt_leadingMu',           100, 0, 300)
    LepPt_hist_mu_1stlep_barrel = TH1F('muPt_leadingMu_barrel',    'muPt_leadingMu_barrel',    100, 0, 300)
    LepPt_hist_mu_1stlep_endcap = TH1F('muPt_leadingMu_endcap',    'muPt_leadingMu_endcap',    100, 0, 300)
    LepPt_hist_mu_2ndlep        = TH1F('muPt_subleadingMu',        'muPt_subleadingMu',        100, 0, 300)
    LepPt_hist_mu_2ndlep_barrel = TH1F('muPt_subleadingMu_barrel', 'muPt_subleadingMu_barrel', 100, 0, 300)
    LepPt_hist_mu_2ndlep_endcap = TH1F('muPt_subleadingMu_endcap', 'muPt_subleadingMu_endcap', 100, 0, 300)

    LepEta_hist_mu_1stlep        = TH1F('muEta_leadingMu',           'muEta_leadingMu',           40, -3, 3)
    LepEta_hist_mu_1stlep_barrel = TH1F('muEta_leadingMu_barrel',    'muEta_leadingMu_barrel',    40, -3, 3)
    LepEta_hist_mu_1stlep_endcap = TH1F('muEta_leadingMu_endcap',    'muEta_leadingMu_endcap',    40, -3, 3)
    LepEta_hist_mu_2ndlep        = TH1F('muEta_subleadingMu',        'muEta_subleadingMu',        40, -3, 3)
    LepEta_hist_mu_2ndlep_barrel = TH1F('muEta_subleadingMu_barrel', 'muEta_subleadingMu_barrel', 40, -3, 3)
    LepEta_hist_mu_2ndlep_endcap = TH1F('muEta_subleadingMu_endcap', 'muEta_subleadingMu_endcap', 40, -3, 3)

    LepPhi_hist_mu_1stlep        = TH1F('muPhi_leadingMu',           'muPhi_leadingMu',           35, -3.3, 3.3)
    LepPhi_hist_mu_1stlep_barrel = TH1F('muPhi_leadingMu_barrel',    'muPhi_leadingMu_barrel',    35, -3.3, 3.3)
    LepPhi_hist_mu_1stlep_endcap = TH1F('muPhi_leadingMu_endcap',    'muPhi_leadingMu_endcap',    35, -3.3, 3.3)
    LepPhi_hist_mu_2ndlep        = TH1F('muPhi_subleadingMu',        'muPhi_subleadingMu',        35, -3.3, 3.3)
    LepPhi_hist_mu_2ndlep_barrel = TH1F('muPhi_subleadingMu_barrel', 'muPhi_subleadingMu_barrel', 35, -3.3, 3.3)
    LepPhi_hist_mu_2ndlep_endcap = TH1F('muPhi_subleadingMu_endcap', 'muPhi_subleadingMu_endcap', 35, -3.3, 3.3)




    treeDATA.SetBranchStatus("*",0)  # disable all branches
    if ZTree :
        treeDATA.SetBranchStatus("Zsel",1)
    else : 
        treeDATA.SetBranchStatus("ZZsel",1)
    treeDATA.SetBranchStatus("LepLepId",1)
    treeDATA.SetBranchStatus("LepPt",1)
    treeDATA.SetBranchStatus("LepEta",1)
    treeDATA.SetBranchStatus("LepPhi",1)
    treeDATA.SetBranchStatus("ZMass",1)


    # read tree 
    print "reading tree", inputDATAtree.GetName(),treeText,treeDATA.GetName()  ,"..."
    for event in treeDATA:
        if ZTree :
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger
        else :
            if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger



        # if (max(event.LepPt[0],event.LepPt[1])<20 || min(event.LepPt[0],event.LepPt[1])<10) : continue    # cut as in the 4l sample
        # if max(event.LepPt[0],event.LepPt[1]) < 30. : continue  # cut on the leading lepton pt 

        if event.ZMass < 40. : continue    # cut in ZMass
        # if event.ZMass < 80. or event.ZMass > 100. : continue  # cut in ZMass


        
        # Z->ee histos
        if(int(math.fabs(event.LepLepId[0])) == 11 ):


            if event.LepPt[0] >= event.LepPt[1] :

                LepPt_hist_ele_1stlep.Fill(event.LepPt[0])
                LepPt_hist_ele_2ndlep.Fill(event.LepPt[1])
                LepEta_hist_ele_1stlep.Fill(event.LepEta[0])
                LepEta_hist_ele_2ndlep.Fill(event.LepEta[1])
                LepPhi_hist_ele_1stlep.Fill(event.LepPhi[0])
                LepPhi_hist_ele_2ndlep.Fill(event.LepPhi[1])

                
                if math.fabs(event.LepEta[0]) <= 1.479 :
                    LepPt_hist_ele_1stlep_barrel.Fill(event.LepPt[0])
                    LepEta_hist_ele_1stlep_barrel.Fill(event.LepEta[0])
                    LepPhi_hist_ele_1stlep_barrel.Fill(event.LepPhi[0])
                    
                elif math.fabs(event.LepEta[0]) > 1.479 :
                    LepPt_hist_ele_1stlep_endcap.Fill(event.LepPt[0])
                    LepEta_hist_ele_1stlep_endcap.Fill(event.LepEta[0])
                    LepPhi_hist_ele_1stlep_endcap.Fill(event.LepPhi[0])
                    

                if math.fabs(event.LepEta[1]) <= 1.479 :
                    LepPt_hist_ele_2ndlep_barrel.Fill(event.LepPt[1])
                    LepEta_hist_ele_2ndlep_barrel.Fill(event.LepEta[1])
                    LepPhi_hist_ele_2ndlep_barrel.Fill(event.LepPhi[1])

                elif math.fabs(event.LepEta[1]) > 1.479 :
                    LepPt_hist_ele_2ndlep_endcap.Fill(event.LepPt[1])
                    LepEta_hist_ele_2ndlep_endcap.Fill(event.LepEta[1])
                    LepPhi_hist_ele_2ndlep_endcap.Fill(event.LepPhi[1])
                    

            else :
                LepPt_hist_ele_1stlep.Fill(event.LepPt[1])
                LepPt_hist_ele_2ndlep.Fill(event.LepPt[0])
                LepEta_hist_ele_1stlep.Fill(event.LepEta[1])
                LepEta_hist_ele_2ndlep.Fill(event.LepEta[0])
                LepPhi_hist_ele_1stlep.Fill(event.LepPhi[1])
                LepPhi_hist_ele_2ndlep.Fill(event.LepPhi[0])

                
                if math.fabs(event.LepEta[1]) <= 1.479 :
                    LepPt_hist_ele_1stlep_barrel.Fill(event.LepPt[1])
                    LepEta_hist_ele_1stlep_barrel.Fill(event.LepEta[1])
                    LepPhi_hist_ele_1stlep_barrel.Fill(event.LepPhi[1])
                    
                elif math.fabs(event.LepEta[1]) > 1.479 :
                    LepPt_hist_ele_1stlep_endcap.Fill(event.LepPt[1])
                    LepEta_hist_ele_1stlep_endcap.Fill(event.LepEta[1])
                    LepPhi_hist_ele_1stlep_endcap.Fill(event.LepPhi[1])
                    

                if math.fabs(event.LepEta[0]) <= 1.479 :
                    LepPt_hist_ele_2ndlep_barrel.Fill(event.LepPt[0])
                    LepEta_hist_ele_2ndlep_barrel.Fill(event.LepEta[0])
                    LepPhi_hist_ele_2ndlep_barrel.Fill(event.LepPhi[0])

                elif math.fabs(event.LepEta[0]) > 1.479 :
                    LepPt_hist_ele_2ndlep_endcap.Fill(event.LepPt[0])
                    LepEta_hist_ele_2ndlep_endcap.Fill(event.LepEta[0])
                    LepPhi_hist_ele_2ndlep_endcap.Fill(event.LepPhi[0])
                


        # Z->mumu histos
        if(int(math.fabs(event.LepLepId[0])) == 13 ):

            
            if event.LepPt[0] >= event.LepPt[1] :

                LepPt_hist_mu_1stlep.Fill(event.LepPt[0])
                LepPt_hist_mu_2ndlep.Fill(event.LepPt[1])
                LepEta_hist_mu_1stlep.Fill(event.LepEta[0])
                LepEta_hist_mu_2ndlep.Fill(event.LepEta[1])
                LepPhi_hist_mu_1stlep.Fill(event.LepPhi[0])
                LepPhi_hist_mu_2ndlep.Fill(event.LepPhi[1])

                
                if math.fabs(event.LepEta[0]) <= 1. :
                    LepPt_hist_mu_1stlep_barrel.Fill(event.LepPt[0])
                    LepEta_hist_mu_1stlep_barrel.Fill(event.LepEta[0])
                    LepPhi_hist_mu_1stlep_barrel.Fill(event.LepPhi[0])
                    
                elif math.fabs(event.LepEta[0]) > 1. :
                    LepPt_hist_mu_1stlep_endcap.Fill(event.LepPt[0])
                    LepEta_hist_mu_1stlep_endcap.Fill(event.LepEta[0])
                    LepPhi_hist_mu_1stlep_endcap.Fill(event.LepPhi[0])
                    

                if math.fabs(event.LepEta[1]) <= 1. :
                    LepPt_hist_mu_2ndlep_barrel.Fill(event.LepPt[1])
                    LepEta_hist_mu_2ndlep_barrel.Fill(event.LepEta[1])
                    LepPhi_hist_mu_2ndlep_barrel.Fill(event.LepPhi[1])

                elif math.fabs(event.LepEta[1]) > 1. :
                    LepPt_hist_mu_2ndlep_endcap.Fill(event.LepPt[1])
                    LepEta_hist_mu_2ndlep_endcap.Fill(event.LepEta[1])
                    LepPhi_hist_mu_2ndlep_endcap.Fill(event.LepPhi[1])
                    

            else :
                LepPt_hist_mu_1stlep.Fill(event.LepPt[1])
                LepPt_hist_mu_2ndlep.Fill(event.LepPt[0])
                LepEta_hist_mu_1stlep.Fill(event.LepEta[1])
                LepEta_hist_mu_2ndlep.Fill(event.LepEta[0])
                LepPhi_hist_mu_1stlep.Fill(event.LepPhi[1])
                LepPhi_hist_mu_2ndlep.Fill(event.LepPhi[0])

                
                if math.fabs(event.LepEta[1]) <= 1. :
                    LepPt_hist_mu_1stlep_barrel.Fill(event.LepPt[1])
                    LepEta_hist_mu_1stlep_barrel.Fill(event.LepEta[1])
                    LepPhi_hist_mu_1stlep_barrel.Fill(event.LepPhi[1])
                    
                elif math.fabs(event.LepEta[1]) > 1. :
                    LepPt_hist_mu_1stlep_endcap.Fill(event.LepPt[1])
                    LepEta_hist_mu_1stlep_endcap.Fill(event.LepEta[1])
                    LepPhi_hist_mu_1stlep_endcap.Fill(event.LepPhi[1])
                    

                if math.fabs(event.LepEta[0]) <= 1. :
                    LepPt_hist_mu_2ndlep_barrel.Fill(event.LepPt[0])
                    LepEta_hist_mu_2ndlep_barrel.Fill(event.LepEta[0])
                    LepPhi_hist_mu_2ndlep_barrel.Fill(event.LepPhi[0])

                elif math.fabs(event.LepEta[0]) > 1. :
                    LepPt_hist_mu_2ndlep_endcap.Fill(event.LepPt[0])
                    LepEta_hist_mu_2ndlep_endcap.Fill(event.LepEta[0])
                    LepPhi_hist_mu_2ndlep_endcap.Fill(event.LepPhi[0])
                



    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_DATA = TFile.Open("LepPtEtaPhiDistrib_DATA_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_DATA.cd()


    # Z->ee
    LepPt_hist_ele_1stlep.Write()
    LepPt_hist_ele_1stlep_barrel.Write()
    LepPt_hist_ele_1stlep_endcap.Write()
    LepPt_hist_ele_2ndlep.Write()
    LepPt_hist_ele_2ndlep_barrel.Write()
    LepPt_hist_ele_2ndlep_endcap.Write()
                                     
    LepEta_hist_ele_1stlep.Write()
    LepEta_hist_ele_1stlep_barrel.Write()
    LepEta_hist_ele_1stlep_endcap.Write()
    LepEta_hist_ele_2ndlep.Write()
    LepEta_hist_ele_2ndlep_barrel.Write()
    LepEta_hist_ele_2ndlep_endcap.Write()
                                     
    LepPhi_hist_ele_1stlep.Write()
    LepPhi_hist_ele_1stlep_barrel.Write()
    LepPhi_hist_ele_1stlep_endcap.Write()
    LepPhi_hist_ele_2ndlep.Write()
    LepPhi_hist_ele_2ndlep_barrel.Write()
    LepPhi_hist_ele_2ndlep_endcap.Write()
                                     
                                     
    # Z->mumu 
    LepPt_hist_mu_1stlep.Write()
    LepPt_hist_mu_1stlep_barrel.Write()
    LepPt_hist_mu_1stlep_endcap.Write()
    LepPt_hist_mu_2ndlep.Write()
    LepPt_hist_mu_2ndlep_barrel.Write()
    LepPt_hist_mu_2ndlep_endcap.Write()
                                     
    LepEta_hist_mu_1stlep.Write()
    LepEta_hist_mu_1stlep_barrel.Write()
    LepEta_hist_mu_1stlep_endcap.Write()
    LepEta_hist_mu_2ndlep.Write()
    LepEta_hist_mu_2ndlep_barrel.Write()
    LepEta_hist_mu_2ndlep_endcap.Write()
                                     
    LepPhi_hist_mu_1stlep.Write()
    LepPhi_hist_mu_1stlep_barrel.Write()
    LepPhi_hist_mu_1stlep_endcap.Write()
    LepPhi_hist_mu_2ndlep.Write()
    LepPhi_hist_mu_2ndlep_barrel.Write()
    LepPhi_hist_mu_2ndlep_endcap.Write()


    outFile_DATA.Close()
    print "DATA histo file created!"



# ********************
#  do MC DY histos
# ********************
if(redoMCDYHistos) :

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on


    # define data histograms Z->ee
    LepPt_hist_ele_1stlep_MC_DY        = TH1F('elePt_leadingEle_MC_DY',           'elePt_leadingEle_MC_DY',           100, 0, 300)
    LepPt_hist_ele_1stlep_barrel_MC_DY = TH1F('elePt_leadingEle_barrel_MC_DY',    'elePt_leadingEle_barrel_MC_DY',    100, 0, 300)
    LepPt_hist_ele_1stlep_endcap_MC_DY = TH1F('elePt_leadingEle_endcap_MC_DY',    'elePt_leadingEle_endcap_MC_DY',    100, 0, 300)
    LepPt_hist_ele_2ndlep_MC_DY        = TH1F('elePt_subleadingEle_MC_DY',        'elePt_subleadingEle_MC_DY',        100, 0, 300)
    LepPt_hist_ele_2ndlep_barrel_MC_DY = TH1F('elePt_subleadingEle_barrel_MC_DY', 'elePt_subleadingEle_barrel_MC_DY', 100, 0, 300)
    LepPt_hist_ele_2ndlep_endcap_MC_DY = TH1F('elePt_subleadingEle_endcap_MC_DY', 'elePt_subleadingEle_endcap_MC_DY', 100, 0, 300)

    LepEta_hist_ele_1stlep_MC_DY        = TH1F('eleEta_leadingEle_MC_DY',           'eleEta_leadingEle_MC_DY',           40, -3, 3)
    LepEta_hist_ele_1stlep_barrel_MC_DY = TH1F('eleEta_leadingEle_barrel_MC_DY',    'eleEta_leadingEle_barrel_MC_DY',    40, -3, 3)
    LepEta_hist_ele_1stlep_endcap_MC_DY = TH1F('eleEta_leadingEle_endcap_MC_DY',    'eleEta_leadingEle_endcap_MC_DY',    40, -3, 3)
    LepEta_hist_ele_2ndlep_MC_DY        = TH1F('eleEta_subleadingEle_MC_DY',        'eleEta_subleadingEle_MC_DY',        40, -3, 3)
    LepEta_hist_ele_2ndlep_barrel_MC_DY = TH1F('eleEta_subleadingEle_barrel_MC_DY', 'eleEta_subleadingEle_barrel_MC_DY', 40, -3, 3)
    LepEta_hist_ele_2ndlep_endcap_MC_DY = TH1F('eleEta_subleadingEle_endcap_MC_DY', 'eleEta_subleadingEle_endcap_MC_DY', 40, -3, 3)

    LepPhi_hist_ele_1stlep_MC_DY        = TH1F('elePhi_leadingEle_MC_DY',           'elePhi_leadingEle_MC_DY',           35, -3.3, 3.3)
    LepPhi_hist_ele_1stlep_barrel_MC_DY = TH1F('elePhi_leadingEle_barrel_MC_DY',    'elePhi_leadingEle_barrel_MC_DY',    35, -3.3, 3.3)
    LepPhi_hist_ele_1stlep_endcap_MC_DY = TH1F('elePhi_leadingEle_endcap_MC_DY',    'elePhi_leadingEle_endcap_MC_DY',    35, -3.3, 3.3)
    LepPhi_hist_ele_2ndlep_MC_DY        = TH1F('elePhi_subleadingEle_MC_DY',        'elePhi_subleadingEle_MC_DY',        35, -3.3, 3.3)
    LepPhi_hist_ele_2ndlep_barrel_MC_DY = TH1F('elePhi_subleadingEle_barrel_MC_DY', 'elePhi_subleadingEle_barrel_MC_DY', 35, -3.3, 3.3)
    LepPhi_hist_ele_2ndlep_endcap_MC_DY = TH1F('elePhi_subleadingEle_endcap_MC_DY', 'elePhi_subleadingEle_endcap_MC_DY', 35, -3.3, 3.3)


    # define data histograms Z->mumu
    LepPt_hist_mu_1stlep_MC_DY        = TH1F('muPt_leadingMu_MC_DY',           'muPt_leadingMu_MC_DY',           100, 0, 300)
    LepPt_hist_mu_1stlep_barrel_MC_DY = TH1F('muPt_leadingMu_barrel_MC_DY',    'muPt_leadingMu_barrel_MC_DY',    100, 0, 300)
    LepPt_hist_mu_1stlep_endcap_MC_DY = TH1F('muPt_leadingMu_endcap_MC_DY',    'muPt_leadingMu_endcap_MC_DY',    100, 0, 300)
    LepPt_hist_mu_2ndlep_MC_DY        = TH1F('muPt_subleadingMu_MC_DY',        'muPt_subleadingMu_MC_DY',        100, 0, 300)
    LepPt_hist_mu_2ndlep_barrel_MC_DY = TH1F('muPt_subleadingMu_barrel_MC_DY', 'muPt_subleadingMu_barrel_MC_DY', 100, 0, 300)
    LepPt_hist_mu_2ndlep_endcap_MC_DY = TH1F('muPt_subleadingMu_endcap_MC_DY', 'muPt_subleadingMu_endcap_MC_DY', 100, 0, 300)

    LepEta_hist_mu_1stlep_MC_DY        = TH1F('muEta_leadingMu_MC_DY',           'muEta_leadingMu_MC_DY',           40, -3, 3)
    LepEta_hist_mu_1stlep_barrel_MC_DY = TH1F('muEta_leadingMu_barrel_MC_DY',    'muEta_leadingMu_barrel_MC_DY',    40, -3, 3)
    LepEta_hist_mu_1stlep_endcap_MC_DY = TH1F('muEta_leadingMu_endcap_MC_DY',    'muEta_leadingMu_endcap_MC_DY',    40, -3, 3)
    LepEta_hist_mu_2ndlep_MC_DY        = TH1F('muEta_subleadingMu_MC_DY',        'muEta_subleadingMu_MC_DY',        40, -3, 3)
    LepEta_hist_mu_2ndlep_barrel_MC_DY = TH1F('muEta_subleadingMu_barrel_MC_DY', 'muEta_subleadingMu_barrel_MC_DY', 40, -3, 3)
    LepEta_hist_mu_2ndlep_endcap_MC_DY = TH1F('muEta_subleadingMu_endcap_MC_DY', 'muEta_subleadingMu_endcap_MC_DY', 40, -3, 3)

    LepPhi_hist_mu_1stlep_MC_DY        = TH1F('muPhi_leadingMu_MC_DY',           'muPhi_leadingMu_MC_DY',           35, -3.3, 3.3)
    LepPhi_hist_mu_1stlep_barrel_MC_DY = TH1F('muPhi_leadingMu_barrel_MC_DY',    'muPhi_leadingMu_barrel_MC_DY',    35, -3.3, 3.3)
    LepPhi_hist_mu_1stlep_endcap_MC_DY = TH1F('muPhi_leadingMu_endcap_MC_DY',    'muPhi_leadingMu_endcap_MC_DY',    35, -3.3, 3.3)
    LepPhi_hist_mu_2ndlep_MC_DY        = TH1F('muPhi_subleadingMu_MC_DY',        'muPhi_subleadingMu_MC_DY',        35, -3.3, 3.3)
    LepPhi_hist_mu_2ndlep_barrel_MC_DY = TH1F('muPhi_subleadingMu_barrel_MC_DY', 'muPhi_subleadingMu_barrel_MC_DY', 35, -3.3, 3.3)
    LepPhi_hist_mu_2ndlep_endcap_MC_DY = TH1F('muPhi_subleadingMu_endcap_MC_DY', 'muPhi_subleadingMu_endcap_MC_DY', 35, -3.3, 3.3)



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



        # if (max(event.LepPt[0],event.LepPt[1])<20 || min(event.LepPt[0],event.LepPt[1])<10) : continue    # cut as in the 4l sample
        # if max(event.LepPt[0],event.LepPt[1]) < 30. : continue  # cut on the leading lepton pt 

        if event.ZMass < 40. : continue    # cut in ZMass
        # if event.ZMass < 80. or event.ZMass > 100. : continue  # cut in ZMass




        weight = partialSampleWeight*event.xsec*event.overallEventWeight


        # Z->ee histos
        if(int(math.fabs(event.LepLepId[0])) == 11 ):

            
            if event.LepPt[0] >= event.LepPt[1] :

                LepPt_hist_ele_1stlep_MC_DY.Fill(event.LepPt[0],weight)
                LepPt_hist_ele_2ndlep_MC_DY.Fill(event.LepPt[1],weight)
                LepEta_hist_ele_1stlep_MC_DY.Fill(event.LepEta[0],weight)
                LepEta_hist_ele_2ndlep_MC_DY.Fill(event.LepEta[1],weight)
                LepPhi_hist_ele_1stlep_MC_DY.Fill(event.LepPhi[0],weight)
                LepPhi_hist_ele_2ndlep_MC_DY.Fill(event.LepPhi[1],weight)

                
                if math.fabs(event.LepEta[0]) <= 1.479 :
                    LepPt_hist_ele_1stlep_barrel_MC_DY.Fill(event.LepPt[0],weight)
                    LepEta_hist_ele_1stlep_barrel_MC_DY.Fill(event.LepEta[0],weight)
                    LepPhi_hist_ele_1stlep_barrel_MC_DY.Fill(event.LepPhi[0],weight)
                    
                elif math.fabs(event.LepEta[0]) > 1.479 :
                    LepPt_hist_ele_1stlep_endcap_MC_DY.Fill(event.LepPt[0],weight)
                    LepEta_hist_ele_1stlep_endcap_MC_DY.Fill(event.LepEta[0],weight)
                    LepPhi_hist_ele_1stlep_endcap_MC_DY.Fill(event.LepPhi[0],weight)
                    

                if math.fabs(event.LepEta[1]) <= 1.479 :
                    LepPt_hist_ele_2ndlep_barrel_MC_DY.Fill(event.LepPt[1],weight)
                    LepEta_hist_ele_2ndlep_barrel_MC_DY.Fill(event.LepEta[1],weight)
                    LepPhi_hist_ele_2ndlep_barrel_MC_DY.Fill(event.LepPhi[1],weight)

                elif math.fabs(event.LepEta[1]) > 1.479 :
                    LepPt_hist_ele_2ndlep_endcap_MC_DY.Fill(event.LepPt[1],weight)
                    LepEta_hist_ele_2ndlep_endcap_MC_DY.Fill(event.LepEta[1],weight)
                    LepPhi_hist_ele_2ndlep_endcap_MC_DY.Fill(event.LepPhi[1],weight)
                    

            else :
                LepPt_hist_ele_1stlep_MC_DY.Fill(event.LepPt[1],weight)
                LepPt_hist_ele_2ndlep_MC_DY.Fill(event.LepPt[0],weight)
                LepEta_hist_ele_1stlep_MC_DY.Fill(event.LepEta[1],weight)
                LepEta_hist_ele_2ndlep_MC_DY.Fill(event.LepEta[0],weight)
                LepPhi_hist_ele_1stlep_MC_DY.Fill(event.LepPhi[1],weight)
                LepPhi_hist_ele_2ndlep_MC_DY.Fill(event.LepPhi[0],weight)

                
                if math.fabs(event.LepEta[1]) <= 1.479 :
                    LepPt_hist_ele_1stlep_barrel_MC_DY.Fill(event.LepPt[1],weight)
                    LepEta_hist_ele_1stlep_barrel_MC_DY.Fill(event.LepEta[1],weight)
                    LepPhi_hist_ele_1stlep_barrel_MC_DY.Fill(event.LepPhi[1],weight)
                    
                elif math.fabs(event.LepEta[1]) > 1.479 :
                    LepPt_hist_ele_1stlep_endcap_MC_DY.Fill(event.LepPt[1],weight)
                    LepEta_hist_ele_1stlep_endcap_MC_DY.Fill(event.LepEta[1],weight)
                    LepPhi_hist_ele_1stlep_endcap_MC_DY.Fill(event.LepPhi[1],weight)
                    

                if math.fabs(event.LepEta[0]) <= 1.479 :
                    LepPt_hist_ele_2ndlep_barrel_MC_DY.Fill(event.LepPt[0],weight)
                    LepEta_hist_ele_2ndlep_barrel_MC_DY.Fill(event.LepEta[0],weight)
                    LepPhi_hist_ele_2ndlep_barrel_MC_DY.Fill(event.LepPhi[0],weight)

                elif math.fabs(event.LepEta[0]) > 1.479 :
                    LepPt_hist_ele_2ndlep_endcap_MC_DY.Fill(event.LepPt[0],weight)
                    LepEta_hist_ele_2ndlep_endcap_MC_DY.Fill(event.LepEta[0],weight)
                    LepPhi_hist_ele_2ndlep_endcap_MC_DY.Fill(event.LepPhi[0],weight)
                


        # Z->mumu histos
        if(int(math.fabs(event.LepLepId[0])) == 13 ):

            
            if event.LepPt[0] >= event.LepPt[1] :

                LepPt_hist_mu_1stlep_MC_DY.Fill(event.LepPt[0],weight)
                LepPt_hist_mu_2ndlep_MC_DY.Fill(event.LepPt[1],weight)
                LepEta_hist_mu_1stlep_MC_DY.Fill(event.LepEta[0],weight)
                LepEta_hist_mu_2ndlep_MC_DY.Fill(event.LepEta[1],weight)
                LepPhi_hist_mu_1stlep_MC_DY.Fill(event.LepPhi[0],weight)
                LepPhi_hist_mu_2ndlep_MC_DY.Fill(event.LepPhi[1],weight)

                
                if math.fabs(event.LepEta[0]) <= 1. :
                    LepPt_hist_mu_1stlep_barrel_MC_DY.Fill(event.LepPt[0],weight)
                    LepEta_hist_mu_1stlep_barrel_MC_DY.Fill(event.LepEta[0],weight)
                    LepPhi_hist_mu_1stlep_barrel_MC_DY.Fill(event.LepPhi[0],weight)
                    
                elif math.fabs(event.LepEta[0]) > 1. :
                    LepPt_hist_mu_1stlep_endcap_MC_DY.Fill(event.LepPt[0],weight)
                    LepEta_hist_mu_1stlep_endcap_MC_DY.Fill(event.LepEta[0],weight)
                    LepPhi_hist_mu_1stlep_endcap_MC_DY.Fill(event.LepPhi[0],weight)
                    

                if math.fabs(event.LepEta[1]) <= 1. :
                    LepPt_hist_mu_2ndlep_barrel_MC_DY.Fill(event.LepPt[1],weight)
                    LepEta_hist_mu_2ndlep_barrel_MC_DY.Fill(event.LepEta[1],weight)
                    LepPhi_hist_mu_2ndlep_barrel_MC_DY.Fill(event.LepPhi[1],weight)

                elif math.fabs(event.LepEta[1]) > 1. :
                    LepPt_hist_mu_2ndlep_endcap_MC_DY.Fill(event.LepPt[1],weight)
                    LepEta_hist_mu_2ndlep_endcap_MC_DY.Fill(event.LepEta[1],weight)
                    LepPhi_hist_mu_2ndlep_endcap_MC_DY.Fill(event.LepPhi[1],weight)
                    

            else :
                LepPt_hist_mu_1stlep_MC_DY.Fill(event.LepPt[1],weight)
                LepPt_hist_mu_2ndlep_MC_DY.Fill(event.LepPt[0],weight)
                LepEta_hist_mu_1stlep_MC_DY.Fill(event.LepEta[1],weight)
                LepEta_hist_mu_2ndlep_MC_DY.Fill(event.LepEta[0],weight)
                LepPhi_hist_mu_1stlep_MC_DY.Fill(event.LepPhi[1],weight)
                LepPhi_hist_mu_2ndlep_MC_DY.Fill(event.LepPhi[0],weight)

                
                if math.fabs(event.LepEta[1]) <= 1. :
                    LepPt_hist_mu_1stlep_barrel_MC_DY.Fill(event.LepPt[1],weight)
                    LepEta_hist_mu_1stlep_barrel_MC_DY.Fill(event.LepEta[1],weight)
                    LepPhi_hist_mu_1stlep_barrel_MC_DY.Fill(event.LepPhi[1],weight)
                    
                elif math.fabs(event.LepEta[1]) > 1. :
                    LepPt_hist_mu_1stlep_endcap_MC_DY.Fill(event.LepPt[1],weight)
                    LepEta_hist_mu_1stlep_endcap_MC_DY.Fill(event.LepEta[1],weight)
                    LepPhi_hist_mu_1stlep_endcap_MC_DY.Fill(event.LepPhi[1],weight)
                    

                if math.fabs(event.LepEta[0]) <= 1. :
                    LepPt_hist_mu_2ndlep_barrel_MC_DY.Fill(event.LepPt[0],weight)
                    LepEta_hist_mu_2ndlep_barrel_MC_DY.Fill(event.LepEta[0],weight)
                    LepPhi_hist_mu_2ndlep_barrel_MC_DY.Fill(event.LepPhi[0],weight)

                elif math.fabs(event.LepEta[0]) > 1. :
                    LepPt_hist_mu_2ndlep_endcap_MC_DY.Fill(event.LepPt[0],weight)
                    LepEta_hist_mu_2ndlep_endcap_MC_DY.Fill(event.LepEta[0],weight)
                    LepPhi_hist_mu_2ndlep_endcap_MC_DY.Fill(event.LepPhi[0],weight)
                

    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_MCDY = TFile.Open("LepPtEtaPhiDistrib_MC_DY_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_MCDY.cd()


    # Zee histos
    LepPt_hist_ele_1stlep_MC_DY.Write()     
    LepPt_hist_ele_1stlep_barrel_MC_DY.Write()
    LepPt_hist_ele_1stlep_endcap_MC_DY.Write()  
    LepPt_hist_ele_2ndlep_MC_DY.Write()
    LepPt_hist_ele_2ndlep_barrel_MC_DY.Write() 
    LepPt_hist_ele_2ndlep_endcap_MC_DY.Write()
                                         
    LepEta_hist_ele_1stlep_MC_DY.Write()
    LepEta_hist_ele_1stlep_barrel_MC_DY.Write()
    LepEta_hist_ele_1stlep_endcap_MC_DY.Write()
    LepEta_hist_ele_2ndlep_MC_DY.Write()
    LepEta_hist_ele_2ndlep_barrel_MC_DY.Write()
    LepEta_hist_ele_2ndlep_endcap_MC_DY.Write()
                                         
    LepPhi_hist_ele_1stlep_MC_DY.Write()
    LepPhi_hist_ele_1stlep_barrel_MC_DY.Write()
    LepPhi_hist_ele_1stlep_endcap_MC_DY.Write()
    LepPhi_hist_ele_2ndlep_MC_DY.Write()
    LepPhi_hist_ele_2ndlep_barrel_MC_DY.Write()
    LepPhi_hist_ele_2ndlep_endcap_MC_DY.Write()
                                         
                                         
    # Zmumu histos    
    LepPt_hist_mu_1stlep_MC_DY.Write()
    LepPt_hist_mu_1stlep_barrel_MC_DY.Write()
    LepPt_hist_mu_1stlep_endcap_MC_DY.Write()
    LepPt_hist_mu_2ndlep_MC_DY.Write()
    LepPt_hist_mu_2ndlep_barrel_MC_DY.Write()
    LepPt_hist_mu_2ndlep_endcap_MC_DY.Write()
                                         
    LepEta_hist_mu_1stlep_MC_DY.Write()
    LepEta_hist_mu_1stlep_barrel_MC_DY.Write() 
    LepEta_hist_mu_1stlep_endcap_MC_DY.Write()
    LepEta_hist_mu_2ndlep_MC_DY.Write()
    LepEta_hist_mu_2ndlep_barrel_MC_DY.Write() 
    LepEta_hist_mu_2ndlep_endcap_MC_DY.Write()
                                         
    LepPhi_hist_mu_1stlep_MC_DY.Write()
    LepPhi_hist_mu_1stlep_barrel_MC_DY.Write() 
    LepPhi_hist_mu_1stlep_endcap_MC_DY.Write()
    LepPhi_hist_mu_2ndlep_MC_DY.Write()
    LepPhi_hist_mu_2ndlep_barrel_MC_DY.Write() 
    LepPhi_hist_mu_2ndlep_endcap_MC_DY.Write()



    outFile_MCDY.Close()
    print "MC DY histo file created!"
# ********************


# ********************
#  do MC TTbar histos
# ********************
if(redoMCTTbarHistos) :

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on


    # define data histograms Z->ee
    LepPt_hist_ele_1stlep_MC_TTbar        = TH1F('elePt_leadingEle_MC_TTbar',           'elePt_leadingEle_MC_TTbar',           100, 0, 300)
    LepPt_hist_ele_1stlep_barrel_MC_TTbar = TH1F('elePt_leadingEle_barrel_MC_TTbar',    'elePt_leadingEle_barrel_MC_TTbar',    100, 0, 300)
    LepPt_hist_ele_1stlep_endcap_MC_TTbar = TH1F('elePt_leadingEle_endcap_MC_TTbar',    'elePt_leadingEle_endcap_MC_TTbar',    100, 0, 300)
    LepPt_hist_ele_2ndlep_MC_TTbar        = TH1F('elePt_subleadingEle_MC_TTbar',        'elePt_subleadingEle_MC_TTbar',        100, 0, 300)
    LepPt_hist_ele_2ndlep_barrel_MC_TTbar = TH1F('elePt_subleadingEle_barrel_MC_TTbar', 'elePt_subleadingEle_barrel_MC_TTbar', 100, 0, 300)
    LepPt_hist_ele_2ndlep_endcap_MC_TTbar = TH1F('elePt_subleadingEle_endcap_MC_TTbar', 'elePt_subleadingEle_endcap_MC_TTbar', 100, 0, 300)

    LepEta_hist_ele_1stlep_MC_TTbar        = TH1F('eleEta_leadingEle_MC_TTbar',           'eleEta_leadingEle_MC_TTbar',           40, -3, 3)
    LepEta_hist_ele_1stlep_barrel_MC_TTbar = TH1F('eleEta_leadingEle_barrel_MC_TTbar',    'eleEta_leadingEle_barrel_MC_TTbar',    40, -3, 3)
    LepEta_hist_ele_1stlep_endcap_MC_TTbar = TH1F('eleEta_leadingEle_endcap_MC_TTbar',    'eleEta_leadingEle_endcap_MC_TTbar',    40, -3, 3)
    LepEta_hist_ele_2ndlep_MC_TTbar        = TH1F('eleEta_subleadingEle_MC_TTbar',        'eleEta_subleadingEle_MC_TTbar',        40, -3, 3)
    LepEta_hist_ele_2ndlep_barrel_MC_TTbar = TH1F('eleEta_subleadingEle_barrel_MC_TTbar', 'eleEta_subleadingEle_barrel_MC_TTbar', 40, -3, 3)
    LepEta_hist_ele_2ndlep_endcap_MC_TTbar = TH1F('eleEta_subleadingEle_endcap_MC_TTbar', 'eleEta_subleadingEle_endcap_MC_TTbar', 40, -3, 3)

    LepPhi_hist_ele_1stlep_MC_TTbar        = TH1F('elePhi_leadingEle_MC_TTbar',           'elePhi_leadingEle_MC_TTbar',           35, -3.3, 3.3)
    LepPhi_hist_ele_1stlep_barrel_MC_TTbar = TH1F('elePhi_leadingEle_barrel_MC_TTbar',    'elePhi_leadingEle_barrel_MC_TTbar',    35, -3.3, 3.3)
    LepPhi_hist_ele_1stlep_endcap_MC_TTbar = TH1F('elePhi_leadingEle_endcap_MC_TTbar',    'elePhi_leadingEle_endcap_MC_TTbar',    35, -3.3, 3.3)
    LepPhi_hist_ele_2ndlep_MC_TTbar        = TH1F('elePhi_subleadingEle_MC_TTbar',        'elePhi_subleadingEle_MC_TTbar',        35, -3.3, 3.3)
    LepPhi_hist_ele_2ndlep_barrel_MC_TTbar = TH1F('elePhi_subleadingEle_barrel_MC_TTbar', 'elePhi_subleadingEle_barrel_MC_TTbar', 35, -3.3, 3.3)
    LepPhi_hist_ele_2ndlep_endcap_MC_TTbar = TH1F('elePhi_subleadingEle_endcap_MC_TTbar', 'elePhi_subleadingEle_endcap_MC_TTbar', 35, -3.3, 3.3)


    # define data histograms Z->mumu
    LepPt_hist_mu_1stlep_MC_TTbar        = TH1F('muPt_leadingMu_MC_TTbar',           'muPt_leadingMu_MC_TTbar',           100, 0, 300)
    LepPt_hist_mu_1stlep_barrel_MC_TTbar = TH1F('muPt_leadingMu_barrel_MC_TTbar',    'muPt_leadingMu_barrel_MC_TTbar',    100, 0, 300)
    LepPt_hist_mu_1stlep_endcap_MC_TTbar = TH1F('muPt_leadingMu_endcap_MC_TTbar',    'muPt_leadingMu_endcap_MC_TTbar',    100, 0, 300)
    LepPt_hist_mu_2ndlep_MC_TTbar        = TH1F('muPt_subleadingMu_MC_TTbar',        'muPt_subleadingMu_MC_TTbar',        100, 0, 300)
    LepPt_hist_mu_2ndlep_barrel_MC_TTbar = TH1F('muPt_subleadingMu_barrel_MC_TTbar', 'muPt_subleadingMu_barrel_MC_TTbar', 100, 0, 300)
    LepPt_hist_mu_2ndlep_endcap_MC_TTbar = TH1F('muPt_subleadingMu_endcap_MC_TTbar', 'muPt_subleadingMu_endcap_MC_TTbar', 100, 0, 300)

    LepEta_hist_mu_1stlep_MC_TTbar        = TH1F('muEta_leadingMu_MC_TTbar',           'muEta_leadingMu_MC_TTbar',           40, -3, 3)
    LepEta_hist_mu_1stlep_barrel_MC_TTbar = TH1F('muEta_leadingMu_barrel_MC_TTbar',    'muEta_leadingMu_barrel_MC_TTbar',    40, -3, 3)
    LepEta_hist_mu_1stlep_endcap_MC_TTbar = TH1F('muEta_leadingMu_endcap_MC_TTbar',    'muEta_leadingMu_endcap_MC_TTbar',    40, -3, 3)
    LepEta_hist_mu_2ndlep_MC_TTbar        = TH1F('muEta_subleadingMu_MC_TTbar',        'muEta_subleadingMu_MC_TTbar',        40, -3, 3)
    LepEta_hist_mu_2ndlep_barrel_MC_TTbar = TH1F('muEta_subleadingMu_barrel_MC_TTbar', 'muEta_subleadingMu_barrel_MC_TTbar', 40, -3, 3)
    LepEta_hist_mu_2ndlep_endcap_MC_TTbar = TH1F('muEta_subleadingMu_endcap_MC_TTbar', 'muEta_subleadingMu_endcap_MC_TTbar', 40, -3, 3)

    LepPhi_hist_mu_1stlep_MC_TTbar        = TH1F('muPhi_leadingMu_MC_TTbar',           'muPhi_leadingMu_MC_TTbar',           35, -3.3, 3.3)
    LepPhi_hist_mu_1stlep_barrel_MC_TTbar = TH1F('muPhi_leadingMu_barrel_MC_TTbar',    'muPhi_leadingMu_barrel_MC_TTbar',    35, -3.3, 3.3)
    LepPhi_hist_mu_1stlep_endcap_MC_TTbar = TH1F('muPhi_leadingMu_endcap_MC_TTbar',    'muPhi_leadingMu_endcap_MC_TTbar',    35, -3.3, 3.3)
    LepPhi_hist_mu_2ndlep_MC_TTbar        = TH1F('muPhi_subleadingMu_MC_TTbar',        'muPhi_subleadingMu_MC_TTbar',        35, -3.3, 3.3)
    LepPhi_hist_mu_2ndlep_barrel_MC_TTbar = TH1F('muPhi_subleadingMu_barrel_MC_TTbar', 'muPhi_subleadingMu_barrel_MC_TTbar', 35, -3.3, 3.3)
    LepPhi_hist_mu_2ndlep_endcap_MC_TTbar = TH1F('muPhi_subleadingMu_endcap_MC_TTbar', 'muPhi_subleadingMu_endcap_MC_TTbar', 35, -3.3, 3.3)



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



        # if (max(event.LepPt[0],event.LepPt[1])<20 || min(event.LepPt[0],event.LepPt[1])<10) : continue    # cut as in the 4l sample
        # if max(event.LepPt[0],event.LepPt[1]) < 30. : continue  # cut on the leading lepton pt 

        if event.ZMass < 40. : continue    # cut in ZMass
        # if event.ZMass < 80. or event.ZMass > 100. : continue  # cut in ZMass



        weight = partialSampleWeight*event.xsec*event.overallEventWeight


        # Z->ee histos
        if(int(math.fabs(event.LepLepId[0])) == 11 ):

            
            if event.LepPt[0] >= event.LepPt[1] :

                LepPt_hist_ele_1stlep_MC_TTbar.Fill(event.LepPt[0],weight)
                LepPt_hist_ele_2ndlep_MC_TTbar.Fill(event.LepPt[1],weight)
                LepEta_hist_ele_1stlep_MC_TTbar.Fill(event.LepEta[0],weight)
                LepEta_hist_ele_2ndlep_MC_TTbar.Fill(event.LepEta[1],weight)
                LepPhi_hist_ele_1stlep_MC_TTbar.Fill(event.LepPhi[0],weight)
                LepPhi_hist_ele_2ndlep_MC_TTbar.Fill(event.LepPhi[1],weight)

                
                if math.fabs(event.LepEta[0]) <= 1.479 :
                    LepPt_hist_ele_1stlep_barrel_MC_TTbar.Fill(event.LepPt[0],weight)
                    LepEta_hist_ele_1stlep_barrel_MC_TTbar.Fill(event.LepEta[0],weight)
                    LepPhi_hist_ele_1stlep_barrel_MC_TTbar.Fill(event.LepPhi[0],weight)
                    
                elif math.fabs(event.LepEta[0]) > 1.479 :
                    LepPt_hist_ele_1stlep_endcap_MC_TTbar.Fill(event.LepPt[0],weight)
                    LepEta_hist_ele_1stlep_endcap_MC_TTbar.Fill(event.LepEta[0],weight)
                    LepPhi_hist_ele_1stlep_endcap_MC_TTbar.Fill(event.LepPhi[0],weight)
                    

                if math.fabs(event.LepEta[1]) <= 1.479 :
                    LepPt_hist_ele_2ndlep_barrel_MC_TTbar.Fill(event.LepPt[1],weight)
                    LepEta_hist_ele_2ndlep_barrel_MC_TTbar.Fill(event.LepEta[1],weight)
                    LepPhi_hist_ele_2ndlep_barrel_MC_TTbar.Fill(event.LepPhi[1],weight)

                elif math.fabs(event.LepEta[1]) > 1.479 :
                    LepPt_hist_ele_2ndlep_endcap_MC_TTbar.Fill(event.LepPt[1],weight)
                    LepEta_hist_ele_2ndlep_endcap_MC_TTbar.Fill(event.LepEta[1],weight)
                    LepPhi_hist_ele_2ndlep_endcap_MC_TTbar.Fill(event.LepPhi[1],weight)
                    

            else :
                LepPt_hist_ele_1stlep_MC_TTbar.Fill(event.LepPt[1],weight)
                LepPt_hist_ele_2ndlep_MC_TTbar.Fill(event.LepPt[0],weight)
                LepEta_hist_ele_1stlep_MC_TTbar.Fill(event.LepEta[1],weight)
                LepEta_hist_ele_2ndlep_MC_TTbar.Fill(event.LepEta[0],weight)
                LepPhi_hist_ele_1stlep_MC_TTbar.Fill(event.LepPhi[1],weight)
                LepPhi_hist_ele_2ndlep_MC_TTbar.Fill(event.LepPhi[0],weight)

                
                if math.fabs(event.LepEta[1]) <= 1.479 :
                    LepPt_hist_ele_1stlep_barrel_MC_TTbar.Fill(event.LepPt[1],weight)
                    LepEta_hist_ele_1stlep_barrel_MC_TTbar.Fill(event.LepEta[1],weight)
                    LepPhi_hist_ele_1stlep_barrel_MC_TTbar.Fill(event.LepPhi[1],weight)
                    
                elif math.fabs(event.LepEta[1]) > 1.479 :
                    LepPt_hist_ele_1stlep_endcap_MC_TTbar.Fill(event.LepPt[1],weight)
                    LepEta_hist_ele_1stlep_endcap_MC_TTbar.Fill(event.LepEta[1],weight)
                    LepPhi_hist_ele_1stlep_endcap_MC_TTbar.Fill(event.LepPhi[1],weight)
                    

                if math.fabs(event.LepEta[0]) <= 1.479 :
                    LepPt_hist_ele_2ndlep_barrel_MC_TTbar.Fill(event.LepPt[0],weight)
                    LepEta_hist_ele_2ndlep_barrel_MC_TTbar.Fill(event.LepEta[0],weight)
                    LepPhi_hist_ele_2ndlep_barrel_MC_TTbar.Fill(event.LepPhi[0],weight)

                elif math.fabs(event.LepEta[0]) > 1.479 :
                    LepPt_hist_ele_2ndlep_endcap_MC_TTbar.Fill(event.LepPt[0],weight)
                    LepEta_hist_ele_2ndlep_endcap_MC_TTbar.Fill(event.LepEta[0],weight)
                    LepPhi_hist_ele_2ndlep_endcap_MC_TTbar.Fill(event.LepPhi[0],weight)
                


        # Z->mumu histos
        if(int(math.fabs(event.LepLepId[0])) == 13 ):

            
            if event.LepPt[0] >= event.LepPt[1] :

                LepPt_hist_mu_1stlep_MC_TTbar.Fill(event.LepPt[0],weight)
                LepPt_hist_mu_2ndlep_MC_TTbar.Fill(event.LepPt[1],weight)
                LepEta_hist_mu_1stlep_MC_TTbar.Fill(event.LepEta[0],weight)
                LepEta_hist_mu_2ndlep_MC_TTbar.Fill(event.LepEta[1],weight)
                LepPhi_hist_mu_1stlep_MC_TTbar.Fill(event.LepPhi[0],weight)
                LepPhi_hist_mu_2ndlep_MC_TTbar.Fill(event.LepPhi[1],weight)

                
                if math.fabs(event.LepEta[0]) <= 1. :
                    LepPt_hist_mu_1stlep_barrel_MC_TTbar.Fill(event.LepPt[0],weight)
                    LepEta_hist_mu_1stlep_barrel_MC_TTbar.Fill(event.LepEta[0],weight)
                    LepPhi_hist_mu_1stlep_barrel_MC_TTbar.Fill(event.LepPhi[0],weight)
                    
                elif math.fabs(event.LepEta[0]) > 1. :
                    LepPt_hist_mu_1stlep_endcap_MC_TTbar.Fill(event.LepPt[0],weight)
                    LepEta_hist_mu_1stlep_endcap_MC_TTbar.Fill(event.LepEta[0],weight)
                    LepPhi_hist_mu_1stlep_endcap_MC_TTbar.Fill(event.LepPhi[0],weight)
                    

                if math.fabs(event.LepEta[1]) <= 1. :
                    LepPt_hist_mu_2ndlep_barrel_MC_TTbar.Fill(event.LepPt[1],weight)
                    LepEta_hist_mu_2ndlep_barrel_MC_TTbar.Fill(event.LepEta[1],weight)
                    LepPhi_hist_mu_2ndlep_barrel_MC_TTbar.Fill(event.LepPhi[1],weight)

                elif math.fabs(event.LepEta[1]) > 1. :
                    LepPt_hist_mu_2ndlep_endcap_MC_TTbar.Fill(event.LepPt[1],weight)
                    LepEta_hist_mu_2ndlep_endcap_MC_TTbar.Fill(event.LepEta[1],weight)
                    LepPhi_hist_mu_2ndlep_endcap_MC_TTbar.Fill(event.LepPhi[1],weight)
                    

            else :
                LepPt_hist_mu_1stlep_MC_TTbar.Fill(event.LepPt[1],weight)
                LepPt_hist_mu_2ndlep_MC_TTbar.Fill(event.LepPt[0],weight)
                LepEta_hist_mu_1stlep_MC_TTbar.Fill(event.LepEta[1],weight)
                LepEta_hist_mu_2ndlep_MC_TTbar.Fill(event.LepEta[0],weight)
                LepPhi_hist_mu_1stlep_MC_TTbar.Fill(event.LepPhi[1],weight)
                LepPhi_hist_mu_2ndlep_MC_TTbar.Fill(event.LepPhi[0],weight)

                
                if math.fabs(event.LepEta[1]) <= 1. :
                    LepPt_hist_mu_1stlep_barrel_MC_TTbar.Fill(event.LepPt[1],weight)
                    LepEta_hist_mu_1stlep_barrel_MC_TTbar.Fill(event.LepEta[1],weight)
                    LepPhi_hist_mu_1stlep_barrel_MC_TTbar.Fill(event.LepPhi[1],weight)
                    
                elif math.fabs(event.LepEta[1]) > 1. :
                    LepPt_hist_mu_1stlep_endcap_MC_TTbar.Fill(event.LepPt[1],weight)
                    LepEta_hist_mu_1stlep_endcap_MC_TTbar.Fill(event.LepEta[1],weight)
                    LepPhi_hist_mu_1stlep_endcap_MC_TTbar.Fill(event.LepPhi[1],weight)
                    

                if math.fabs(event.LepEta[0]) <= 1. :
                    LepPt_hist_mu_2ndlep_barrel_MC_TTbar.Fill(event.LepPt[0],weight)
                    LepEta_hist_mu_2ndlep_barrel_MC_TTbar.Fill(event.LepEta[0],weight)
                    LepPhi_hist_mu_2ndlep_barrel_MC_TTbar.Fill(event.LepPhi[0],weight)

                elif math.fabs(event.LepEta[0]) > 1. :
                    LepPt_hist_mu_2ndlep_endcap_MC_TTbar.Fill(event.LepPt[0],weight)
                    LepEta_hist_mu_2ndlep_endcap_MC_TTbar.Fill(event.LepEta[0],weight)
                    LepPhi_hist_mu_2ndlep_endcap_MC_TTbar.Fill(event.LepPhi[0],weight)
                

    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_MCTTbar = TFile.Open("LepPtEtaPhiDistrib_MC_TTbar_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_MCTTbar.cd()


    # Zee histos
    LepPt_hist_ele_1stlep_MC_TTbar.Write()     
    LepPt_hist_ele_1stlep_barrel_MC_TTbar.Write()
    LepPt_hist_ele_1stlep_endcap_MC_TTbar.Write()  
    LepPt_hist_ele_2ndlep_MC_TTbar.Write()
    LepPt_hist_ele_2ndlep_barrel_MC_TTbar.Write() 
    LepPt_hist_ele_2ndlep_endcap_MC_TTbar.Write()
                                         
    LepEta_hist_ele_1stlep_MC_TTbar.Write()
    LepEta_hist_ele_1stlep_barrel_MC_TTbar.Write()
    LepEta_hist_ele_1stlep_endcap_MC_TTbar.Write()
    LepEta_hist_ele_2ndlep_MC_TTbar.Write()
    LepEta_hist_ele_2ndlep_barrel_MC_TTbar.Write()
    LepEta_hist_ele_2ndlep_endcap_MC_TTbar.Write()
                                         
    LepPhi_hist_ele_1stlep_MC_TTbar.Write()
    LepPhi_hist_ele_1stlep_barrel_MC_TTbar.Write()
    LepPhi_hist_ele_1stlep_endcap_MC_TTbar.Write()
    LepPhi_hist_ele_2ndlep_MC_TTbar.Write()
    LepPhi_hist_ele_2ndlep_barrel_MC_TTbar.Write()
    LepPhi_hist_ele_2ndlep_endcap_MC_TTbar.Write()
                                         
                                         
    # Zmumu histos    
    LepPt_hist_mu_1stlep_MC_TTbar.Write()
    LepPt_hist_mu_1stlep_barrel_MC_TTbar.Write()
    LepPt_hist_mu_1stlep_endcap_MC_TTbar.Write()
    LepPt_hist_mu_2ndlep_MC_TTbar.Write()
    LepPt_hist_mu_2ndlep_barrel_MC_TTbar.Write()
    LepPt_hist_mu_2ndlep_endcap_MC_TTbar.Write()
                                         
    LepEta_hist_mu_1stlep_MC_TTbar.Write()
    LepEta_hist_mu_1stlep_barrel_MC_TTbar.Write() 
    LepEta_hist_mu_1stlep_endcap_MC_TTbar.Write()
    LepEta_hist_mu_2ndlep_MC_TTbar.Write()
    LepEta_hist_mu_2ndlep_barrel_MC_TTbar.Write() 
    LepEta_hist_mu_2ndlep_endcap_MC_TTbar.Write()
                                         
    LepPhi_hist_mu_1stlep_MC_TTbar.Write()
    LepPhi_hist_mu_1stlep_barrel_MC_TTbar.Write() 
    LepPhi_hist_mu_1stlep_endcap_MC_TTbar.Write()
    LepPhi_hist_mu_2ndlep_MC_TTbar.Write()
    LepPhi_hist_mu_2ndlep_barrel_MC_TTbar.Write() 
    LepPhi_hist_mu_2ndlep_endcap_MC_TTbar.Write()



    outFile_MCTTbar.Close()
    print "MC TTbar histo file created!"
# ********************




# ****************************************
# create output directory 
outputDir = "LepPtEtaPhiDistrib_DATAvsMC_" + str(period) + "_" + str(treeText)
gSystem.Exec("mkdir -p " + outputDir)
print "Output directory created!"



# **************************
# read data histos from file 
histoDATA_input = TFile.Open("LepPtEtaPhiDistrib_DATA_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoDATA_input.GetName(),'...'

inDATA_list = []


inDATA_list.append(histoDATA_input.Get('elePt_leadingEle'))
inDATA_list.append(histoDATA_input.Get('elePt_leadingEle_barrel'))
inDATA_list.append(histoDATA_input.Get('elePt_leadingEle_endcap'))
inDATA_list.append(histoDATA_input.Get('elePt_subleadingEle'))
inDATA_list.append(histoDATA_input.Get('elePt_subleadingEle_barrel'))
inDATA_list.append(histoDATA_input.Get('elePt_subleadingEle_endcap'))
inDATA_list.append(histoDATA_input.Get('eleEta_leadingEle'))
inDATA_list.append(histoDATA_input.Get('eleEta_leadingEle_barrel'))
inDATA_list.append(histoDATA_input.Get('eleEta_leadingEle_endcap'))
inDATA_list.append(histoDATA_input.Get('eleEta_subleadingEle'))
inDATA_list.append(histoDATA_input.Get('eleEta_subleadingEle_barrel'))
inDATA_list.append(histoDATA_input.Get('eleEta_subleadingEle_endcap'))
inDATA_list.append(histoDATA_input.Get('elePhi_leadingEle'))
inDATA_list.append(histoDATA_input.Get('elePhi_leadingEle_barrel'))
inDATA_list.append(histoDATA_input.Get('elePhi_leadingEle_endcap'))
inDATA_list.append(histoDATA_input.Get('elePhi_subleadingEle'))
inDATA_list.append(histoDATA_input.Get('elePhi_subleadingEle_barrel'))
inDATA_list.append(histoDATA_input.Get('elePhi_subleadingEle_endcap'))
inDATA_list.append(histoDATA_input.Get('muPt_leadingMu'))
inDATA_list.append(histoDATA_input.Get('muPt_leadingMu_barrel'))
inDATA_list.append(histoDATA_input.Get('muPt_leadingMu_endcap'))
inDATA_list.append(histoDATA_input.Get('muPt_subleadingMu'))
inDATA_list.append(histoDATA_input.Get('muPt_subleadingMu_barrel'))
inDATA_list.append(histoDATA_input.Get('muPt_subleadingMu_endcap'))
inDATA_list.append(histoDATA_input.Get('muEta_leadingMu'))
inDATA_list.append(histoDATA_input.Get('muEta_leadingMu_barrel'))
inDATA_list.append(histoDATA_input.Get('muEta_leadingMu_endcap'))
inDATA_list.append(histoDATA_input.Get('muEta_subleadingMu'))
inDATA_list.append(histoDATA_input.Get('muEta_subleadingMu_barrel'))
inDATA_list.append(histoDATA_input.Get('muEta_subleadingMu_endcap'))
inDATA_list.append(histoDATA_input.Get('muPhi_leadingMu'))
inDATA_list.append(histoDATA_input.Get('muPhi_leadingMu_barrel'))
inDATA_list.append(histoDATA_input.Get('muPhi_leadingMu_endcap'))
inDATA_list.append(histoDATA_input.Get('muPhi_subleadingMu'))
inDATA_list.append(histoDATA_input.Get('muPhi_subleadingMu_barrel'))
inDATA_list.append(histoDATA_input.Get('muPhi_subleadingMu_endcap'))



# ****************************
# read DY MC histos from file 
histoMCDY_input = TFile.Open("LepPtEtaPhiDistrib_MC_DY_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoMCDY_input.GetName(),'...'

inMCDY_list = []


inMCDY_list.append(histoMCDY_input.Get('elePt_leadingEle_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePt_leadingEle_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePt_leadingEle_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePt_subleadingEle_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePt_subleadingEle_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePt_subleadingEle_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('eleEta_leadingEle_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('eleEta_leadingEle_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('eleEta_leadingEle_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('eleEta_subleadingEle_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('eleEta_subleadingEle_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('eleEta_subleadingEle_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePhi_leadingEle_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePhi_leadingEle_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePhi_leadingEle_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePhi_subleadingEle_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePhi_subleadingEle_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('elePhi_subleadingEle_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPt_leadingMu_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPt_leadingMu_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPt_leadingMu_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPt_subleadingMu_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPt_subleadingMu_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPt_subleadingMu_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muEta_leadingMu_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muEta_leadingMu_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muEta_leadingMu_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muEta_subleadingMu_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muEta_subleadingMu_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muEta_subleadingMu_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPhi_leadingMu_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPhi_leadingMu_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPhi_leadingMu_endcap_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPhi_subleadingMu_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPhi_subleadingMu_barrel_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('muPhi_subleadingMu_endcap_MC_DY'))



# ****************************
# read TTbar MC histos from file 
histoMCTTbar_input = TFile.Open("LepPtEtaPhiDistrib_MC_TTbar_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoMCTTbar_input.GetName(),'...'

inMCTTbar_list = []


inMCTTbar_list.append(histoMCTTbar_input.Get('elePt_leadingEle_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePt_leadingEle_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePt_leadingEle_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePt_subleadingEle_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePt_subleadingEle_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePt_subleadingEle_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('eleEta_leadingEle_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('eleEta_leadingEle_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('eleEta_leadingEle_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('eleEta_subleadingEle_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('eleEta_subleadingEle_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('eleEta_subleadingEle_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePhi_leadingEle_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePhi_leadingEle_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePhi_leadingEle_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePhi_subleadingEle_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePhi_subleadingEle_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('elePhi_subleadingEle_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPt_leadingMu_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPt_leadingMu_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPt_leadingMu_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPt_subleadingMu_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPt_subleadingMu_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPt_subleadingMu_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muEta_leadingMu_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muEta_leadingMu_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muEta_leadingMu_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muEta_subleadingMu_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muEta_subleadingMu_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muEta_subleadingMu_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPhi_leadingMu_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPhi_leadingMu_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPhi_leadingMu_endcap_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPhi_subleadingMu_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPhi_subleadingMu_barrel_MC_TTbar'))
inMCTTbar_list.append(histoMCTTbar_input.Get('muPhi_subleadingMu_endcap_MC_TTbar'))




# ******************************
# do DATA vs MC comparison plots  
for i in range(len(inDATA_list)) : 

    canvas = TCanvas("canvas","canvas",800,800)


    hs = THStack("hs","")


    #DATA hist
    inDATA_list[i].SetMarkerStyle(20)
    inDATA_list[i].SetMarkerSize(0.6)


    #MC TTbar hist
    inMCTTbar_list[i].SetFillColor(kOrange-3)
    inMCTTbar_list[i].SetLineColor(kBlack)
    hs.Add(inMCTTbar_list[i])


    #MC DY hist
    inMCDY_list[i].SetFillColor(kGreen+1)
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
    legend.AddEntry(inMCTTbar_list[i],"t#bar{t}","f")
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
