#!/usr/bin/env python

# *******************
# usage: 
#    python jetsMetDistrib_DATAvsMC.py
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
    inputDATAtree    = TFile.Open("/data3/Higgs/180218/AllData/ZZ4lAnalysis.root")        #2017 data (rereco json 18/2/2018)   
    inputMCDYtree    = TFile.Open("/data3/Higgs/180218/DYJetsToLL_M50/ZZ4lAnalysis.root") #DYJets 2017 MC (rereco json 18/2/2018)
    inputMCTTbartree = TFile.Open("/data3/Higgs/180218/TTTo2L2Nu/ZZ4lAnalysis.root")      #TTbarJets 2017 MC (rereco json 18/2/2018)
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

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    # define data histograms Z->ee

    # nCleanedJets_hist_ele                   = TH1F('nCleanedJets_ele',                   'nCleanedJets_ele',                   30, 0, 30) 
    nCleanedJetsPt30_hist_ele               = TH1F('nCleanedJetsPt30_ele',               'nCleanedJetsPt30_ele',               15, 0, 15) 
    # nCleanedJetsPt30BTagged_hist_ele        = TH1F('nCleanedJetsPt30BTagged_ele',        'nCleanedJetsPt30BTagged_ele',        5,  0, 5)
    # nCleanedJetsPt30BTagged_bTagSF_hist_ele = TH1F('nCleanedJetsPt30BTagged_bTagSF_ele', 'nCleanedJetsPt30BTagged_bTagSF_ele', 5,  0, 5)

    JetPt_hist_ele               = TH1F('JetPt_inclusive_ele',           'JetPt_inclusive_ele',           120, 0, 1200)
    JetPt_hist_1stJet_ele        = TH1F('JetPt_leadingJet_ele',          'JetPt_leadingJet_ele',          120, 0, 1200)
    JetPt_hist_2ndjet_ele        = TH1F('JetPt_subLeadingJet_ele',       'JetPt_subLeadingJet_ele',       120, 0, 1200)
    JetPt_hist_1stJet_fwdeta_ele = TH1F('JetPt_leadingJet_fwdeta_ele',   'JetPt_leadingJet_fwdeta_ele',   120, 0, 1200)
    JetPt_hist_2ndjet_fwdeta_ele = TH1F('JetPt_subLeadingJet_fwdeta_ele','JetPt_subLeadingJet_fwdeta_ele',120, 0, 1200)

    JetEta_hist_ele               = TH1F('JetEta_inclusive_ele',           'JetEta_inclusive_ele',           50, -5, 5)
    JetEta_hist_1stjet_ele        = TH1F('JetEta_leadingJet_ele',          'JetEta_leadingJet_ele',          50, -5, 5)
    JetEta_hist_2ndjet_ele        = TH1F('JetEta_subLeadingJet_ele',       'JetEta_subLeadingJet_ele',       50, -5, 5)
    JetEta_hist_1stjet_fwdeta_ele = TH1F('JetEta_leadingJet_fwdeta_ele',   'JetEta_leadingJet_fwdeta_ele',   50, -5, 5)
    JetEta_hist_2ndjet_fwdeta_ele = TH1F('JetEta_subLeadingJet_fwdeta_ele','JetEta_subLeadingJet_fwdeta_ele',50, -5, 5)

    JetPhi_hist_ele               = TH1F('JetPhi_inclusive_ele',           'JetPhi_inclusive_ele',           33, -3.3, 3.3)
    JetPhi_hist_1stjet_ele        = TH1F('JetPhi_leadingJet_ele',          'JetPhi_leadingJet_ele',          33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_ele        = TH1F('JetPhi_subLeadingJet_ele',       'JetPhi_subLeadingJet_ele',       33, -3.3, 3.3) 
    JetPhi_hist_1stjet_fwdeta_ele = TH1F('JetPhi_leadingJet_fwdeta_ele',   'JetPhi_leadingJet_fwdeta_ele',   33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_fwdeta_ele = TH1F('JetPhi_subLeadingJet_fwdeta_ele','JetPhi_subLeadingJet_fwdeta_ele',33, -3.3, 3.3) 

    JetBTagger_hist_ele               = TH1F('JetBTagger_inclusive_ele',           'JetBTagger_inclusive_ele',           11, 0, 1.1)
    JetBTagger_hist_1stjet_ele        = TH1F('JetBTagger_leadingJet_ele',          'JetBTagger_leadingJet_ele',          11, 0, 1.1)
    JetBTagger_hist_2ndjet_ele        = TH1F('JetBTagger_subLeadingJet_ele',       'JetBTagger_subLeadingJet_ele',       11, 0, 1.1)
    JetBTagger_hist_1stjet_fwdeta_ele = TH1F('JetBTagger_leadingJet_fwdeta_ele',   'JetBTagger_leadingJet_fwdeta_ele',   11, 0, 1.1)
    JetBTagger_hist_2ndjet_fwdeta_ele = TH1F('JetBTagger_subLeadingJet_fwdeta_ele','JetBTagger_subLeadingJet_fwdeta_ele',11, 0, 1.1)

    JetIsBtagged_hist_ele               = TH1F('JetIsBtagged_inclusive_ele',           'JetIsBtagged_inclusive_ele',           11, 0, 1.1)
    JetIsBtagged_hist_1stjet_ele        = TH1F('JetIsBtagged_leadingJet_ele',          'JetIsBtagged_leadingJet_ele',          11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_ele        = TH1F('JetIsBtagged_subLeadingJet_ele',       'JetIsBtagged_subLeadingJet_ele',       11, 0, 1.1)
    JetIsBtagged_hist_1stjet_fwdeta_ele = TH1F('JetIsBtagged_leadingJet_fwdeta_ele',   'JetIsBtagged_leadingJet_fwdeta_ele',   11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_fwdeta_ele = TH1F('JetIsBtagged_subLeadingJet_fwdeta_ele','JetIsBtagged_subLeadingJet_fwdeta_ele',11, 0, 1.1)
    
    JetIsBtaggedWithSF_hist_ele               = TH1F('JetIsBtaggedWithSF_inclusive_ele',           'JetIsBtaggedWithSF_inclusive_ele',           11, 0, 1.1)
    JetIsBtaggedWithSF_hist_1stjet_ele        = TH1F('JetIsBtaggedWithSF_leadingJet_ele',          'JetIsBtaggedWithSF_leadingJet_ele',          11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_ele        = TH1F('JetIsBtaggedWithSF_subLeadingJet_ele',       'JetIsBtaggedWithSF_subLeadingJet_ele',       11, 0, 1.1) 
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_ele = TH1F('JetIsBtaggedWithSF_leadingJet_fwdeta_ele',   'JetIsBtaggedWithSF_leadingJet_fwdeta_ele',   11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_ele = TH1F('JetIsBtaggedWithSF_subLeadingJet_fwdeta_ele','JetIsBtaggedWithSF_subLeadingJet_fwdeta_ele',11, 0, 1.1)

    PFMET_hist_ele = TH1F('PFMET_ele','PFMET_ele',125,0,250)


    # define data histograms Z->mumu

    # nCleanedJets_hist_mu                   = TH1F('nCleanedJets_mu',                   'nCleanedJets_mu',                   30, 0, 30) 
    nCleanedJetsPt30_hist_mu               = TH1F('nCleanedJetsPt30_mu',               'nCleanedJetsPt30_mu',               15, 0, 15) 
    # nCleanedJetsPt30BTagged_hist_mu        = TH1F('nCleanedJetsPt30BTagged_mu',        'nCleanedJetsPt30BTagged_mu',        5,  0, 5)
    # nCleanedJetsPt30BTagged_bTagSF_hist_mu = TH1F('nCleanedJetsPt30BTagged_bTagSF_mu', 'nCleanedJetsPt30BTagged_bTagSF_mu', 5,  0, 5)

    JetPt_hist_mu               = TH1F('JetPt_inclusive_mu',           'JetPt_inclusive_mu',           120, 0, 1200)
    JetPt_hist_1stJet_mu        = TH1F('JetPt_leadingJet_mu',          'JetPt_leadingJet_mu',          120, 0, 1200)
    JetPt_hist_2ndjet_mu        = TH1F('JetPt_subLeadingJet_mu',       'JetPt_subLeadingJet_mu',       120, 0, 1200)
    JetPt_hist_1stJet_fwdeta_mu = TH1F('JetPt_leadingJet_fwdeta_mu',   'JetPt_leadingJet_fwdeta_mu',   120, 0, 1200)
    JetPt_hist_2ndjet_fwdeta_mu = TH1F('JetPt_subLeadingJet_fwdeta_mu','JetPt_subLeadingJet_fwdeta_mu',120, 0, 1200)

    JetEta_hist_mu               = TH1F('JetEta_inclusive_mu',           'JetEta_inclusive_mu',           50, -5, 5)
    JetEta_hist_1stjet_mu        = TH1F('JetEta_leadingJet_mu',          'JetEta_leadingJet_mu',          50, -5, 5)
    JetEta_hist_2ndjet_mu        = TH1F('JetEta_subLeadingJet_mu',       'JetEta_subLeadingJet_mu',       50, -5, 5)
    JetEta_hist_1stjet_fwdeta_mu = TH1F('JetEta_leadingJet_fwdeta_mu',   'JetEta_leadingJet_fwdeta_mu',   50, -5, 5)
    JetEta_hist_2ndjet_fwdeta_mu = TH1F('JetEta_subLeadingJet_fwdeta_mu','JetEta_subLeadingJet_fwdeta_mu',50, -5, 5)

    JetPhi_hist_mu               = TH1F('JetPhi_inclusive_mu',           'JetPhi_inclusive_mu',           33, -3.3, 3.3)
    JetPhi_hist_1stjet_mu        = TH1F('JetPhi_leadingJet_mu',          'JetPhi_leadingJet_mu',          33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_mu        = TH1F('JetPhi_subLeadingJet_mu',       'JetPhi_subLeadingJet_mu',       33, -3.3, 3.3) 
    JetPhi_hist_1stjet_fwdeta_mu = TH1F('JetPhi_leadingJet_fwdeta_mu',   'JetPhi_leadingJet_fwdeta_mu',   33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_fwdeta_mu = TH1F('JetPhi_subLeadingJet_fwdeta_mu','JetPhi_subLeadingJet_fwdeta_mu',33, -3.3, 3.3) 

    JetBTagger_hist_mu               = TH1F('JetBTagger_inclusive_mu',           'JetBTagger_inclusive_mu',           11, 0, 1.1)
    JetBTagger_hist_1stjet_mu        = TH1F('JetBTagger_leadingJet_mu',          'JetBTagger_leadingJet_mu',          11, 0, 1.1)
    JetBTagger_hist_2ndjet_mu        = TH1F('JetBTagger_subLeadingJet_mu',       'JetBTagger_subLeadingJet_mu',       11, 0, 1.1)
    JetBTagger_hist_1stjet_fwdeta_mu = TH1F('JetBTagger_leadingJet_fwdeta_mu',   'JetBTagger_leadingJet_fwdeta_mu',   11, 0, 1.1)
    JetBTagger_hist_2ndjet_fwdeta_mu = TH1F('JetBTagger_subLeadingJet_fwdeta_mu','JetBTagger_subLeadingJet_fwdeta_mu',11, 0, 1.1)

    JetIsBtagged_hist_mu               = TH1F('JetIsBtagged_inclusive_mu',           'JetIsBtagged_inclusive_mu',           11, 0, 1.1)
    JetIsBtagged_hist_1stjet_mu        = TH1F('JetIsBtagged_leadingJet_mu',          'JetIsBtagged_leadingJet_mu',          11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_mu        = TH1F('JetIsBtagged_subLeadingJet_mu',       'JetIsBtagged_subLeadingJet_mu',       11, 0, 1.1)
    JetIsBtagged_hist_1stjet_fwdeta_mu = TH1F('JetIsBtagged_leadingJet_fwdeta_mu',   'JetIsBtagged_leadingJet_fwdeta_mu',   11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_fwdeta_mu = TH1F('JetIsBtagged_subLeadingJet_fwdeta_mu','JetIsBtagged_subLeadingJet_fwdeta_mu',11, 0, 1.1)
    
    JetIsBtaggedWithSF_hist_mu               = TH1F('JetIsBtaggedWithSF_inclusive_mu',           'JetIsBtaggedWithSF_inclusive_mu',           11, 0, 1.1)
    JetIsBtaggedWithSF_hist_1stjet_mu        = TH1F('JetIsBtaggedWithSF_leadingJet_mu',          'JetIsBtaggedWithSF_leadingJet_mu',          11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_mu        = TH1F('JetIsBtaggedWithSF_subLeadingJet_mu',       'JetIsBtaggedWithSF_subLeadingJet_mu',       11, 0, 1.1) 
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_mu = TH1F('JetIsBtaggedWithSF_leadingJet_fwdeta_mu',   'JetIsBtaggedWithSF_leadingJet_fwdeta_mu',   11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_mu = TH1F('JetIsBtaggedWithSF_subLeadingJet_fwdeta_mu','JetIsBtaggedWithSF_subLeadingJet_fwdeta_mu',11, 0, 1.1)

    PFMET_hist_mu = TH1F('PFMET_mu','PFMET_mu',125,0,250)

    

    # read tree 
    print "reading tree", inputDATAtree.GetName(),treeText,treeDATA.GetName()  ,"..."
    for event in treeDATA:
        if ZTree :
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger
        else :
            if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger
  

        # Z->ee histos
        if(int(math.fabs(event.LepLepId[0])) == 11 ):

            # nCleanedJets_hist_ele.Fill(event.nCleanedJets)
            nCleanedJetsPt30_hist_ele.Fill(event.nCleanedJetsPt30)
            # nCleanedJetsPt30BTagged_hist_ele.Fill(event.nCleanedJetsPt30BTagged)
            # nCleanedJetsPt30BTagged_bTagSF_hist_ele.Fill(event.nCleanedJetsPt30BTagged_bTagSF)
            
            for i in range(len(event.JetPt)) :
                JetPt_hist_ele.Fill(event.JetPt[i])
                      
            for i in range(len(event.JetEta)) :  
                JetEta_hist_ele.Fill(event.JetEta[i])
                                            
            for i in range(len(event.JetPhi)) :
                JetPhi_hist_ele.Fill(event.JetPhi[i])
            

            if event.nCleanedJetsPt30 > 0 :

                for i in range(len(event.JetBTagger)) :
                    JetBTagger_hist_ele.Fill(event.JetBTagger[i])
                                            
                for i in range(len(event.JetIsBtagged)) :
                    JetIsBtagged_hist_ele.Fill(event.JetIsBtagged[i])
            
                for i in range(len(event.JetIsBtaggedWithSF)) :                
                    JetIsBtaggedWithSF_hist_ele.Fill(event.JetIsBtaggedWithSF[i])     
            
            
                JetPt_hist_1stJet_ele.Fill(event.JetPt[0])
                JetEta_hist_1stjet_ele.Fill(event.JetEta[0])
                JetPhi_hist_1stjet_ele.Fill(event.JetPhi[0])
                JetBTagger_hist_1stjet_ele.Fill(event.JetBTagger[0])
                JetIsBtagged_hist_1stjet_ele.Fill(event.JetIsBtagged[0])
                JetIsBtaggedWithSF_hist_1stjet_ele.Fill(event.JetIsBtaggedWithSF[0])
            
                if math.fabs(event.JetEta[0]) > 3.0 : 
                    JetPt_hist_1stJet_fwdeta_ele.Fill(event.JetPt[0])
                    JetEta_hist_1stjet_fwdeta_ele.Fill(event.JetEta[0])
                    JetPhi_hist_1stjet_fwdeta_ele.Fill(event.JetPhi[0])
                    JetBTagger_hist_1stjet_fwdeta_ele.Fill(event.JetBTagger[0])
                    JetIsBtagged_hist_1stjet_fwdeta_ele.Fill(event.JetIsBtagged[0])
                    JetIsBtaggedWithSF_hist_1stjet_fwdeta_ele.Fill(event.JetIsBtaggedWithSF[0])
            
            
                if event.nCleanedJetsPt30 > 1 :
                    JetPt_hist_2ndjet_ele.Fill(event.JetPt[1])
                    JetEta_hist_2ndjet_ele.Fill(event.JetEta[1])
                    JetPhi_hist_2ndjet_ele.Fill(event.JetPhi[1])
                    JetBTagger_hist_2ndjet_ele.Fill(event.JetBTagger[1])
                    JetIsBtagged_hist_2ndjet_ele.Fill(event.JetIsBtagged[1])
                    JetIsBtaggedWithSF_hist_2ndjet_ele.Fill(event.JetIsBtaggedWithSF[1])
            
                    if math.fabs(event.JetEta[1]) > 3.0 : 
                        JetPt_hist_2ndjet_fwdeta_ele.Fill(event.JetPt[1])
                        JetEta_hist_2ndjet_fwdeta_ele.Fill(event.JetEta[1])
                        JetPhi_hist_2ndjet_fwdeta_ele.Fill(event.JetPhi[1])
                        JetBTagger_hist_2ndjet_fwdeta_ele.Fill(event.JetBTagger[1])
                        JetIsBtagged_hist_2ndjet_fwdeta_ele.Fill(event.JetIsBtagged[1])
                        JetIsBtaggedWithSF_hist_2ndjet_fwdeta_ele.Fill(event.JetIsBtaggedWithSF[1])
                        
            
            PFMET_hist_ele.Fill(event.PFMET)


        # Z->mumu histos
        if(int(math.fabs(event.LepLepId[0])) == 13 ):

            # nCleanedJets_hist_mu.Fill(event.nCleanedJets)
            nCleanedJetsPt30_hist_mu.Fill(event.nCleanedJetsPt30)
            # nCleanedJetsPt30BTagged_hist_mu.Fill(event.nCleanedJetsPt30BTagged)
            # nCleanedJetsPt30BTagged_bTagSF_hist_mu.Fill(event.nCleanedJetsPt30BTagged_bTagSF)
            
            for i in range(len(event.JetPt)) :
                JetPt_hist_mu.Fill(event.JetPt[i])
                      
            for i in range(len(event.JetEta)) :  
                JetEta_hist_mu.Fill(event.JetEta[i])
                                            
            for i in range(len(event.JetPhi)) :
                JetPhi_hist_mu.Fill(event.JetPhi[i])
            

            if event.nCleanedJetsPt30 > 0 :

                for i in range(len(event.JetBTagger)) :
                    JetBTagger_hist_mu.Fill(event.JetBTagger[i])
                                            
                for i in range(len(event.JetIsBtagged)) :
                    JetIsBtagged_hist_mu.Fill(event.JetIsBtagged[i])
            
                for i in range(len(event.JetIsBtaggedWithSF)) :                
                    JetIsBtaggedWithSF_hist_mu.Fill(event.JetIsBtaggedWithSF[i])     
            
            
                JetPt_hist_1stJet_mu.Fill(event.JetPt[0])
                JetEta_hist_1stjet_mu.Fill(event.JetEta[0])
                JetPhi_hist_1stjet_mu.Fill(event.JetPhi[0])
                JetBTagger_hist_1stjet_mu.Fill(event.JetBTagger[0])
                JetIsBtagged_hist_1stjet_mu.Fill(event.JetIsBtagged[0])
                JetIsBtaggedWithSF_hist_1stjet_mu.Fill(event.JetIsBtaggedWithSF[0])
            
                if math.fabs(event.JetEta[0]) > 3.0 : 
                    JetPt_hist_1stJet_fwdeta_mu.Fill(event.JetPt[0])
                    JetEta_hist_1stjet_fwdeta_mu.Fill(event.JetEta[0])
                    JetPhi_hist_1stjet_fwdeta_mu.Fill(event.JetPhi[0])
                    JetBTagger_hist_1stjet_fwdeta_mu.Fill(event.JetBTagger[0])
                    JetIsBtagged_hist_1stjet_fwdeta_mu.Fill(event.JetIsBtagged[0])
                    JetIsBtaggedWithSF_hist_1stjet_fwdeta_mu.Fill(event.JetIsBtaggedWithSF[0])
            
            
                if event.nCleanedJetsPt30 > 1 :
                    JetPt_hist_2ndjet_mu.Fill(event.JetPt[1])
                    JetEta_hist_2ndjet_mu.Fill(event.JetEta[1])
                    JetPhi_hist_2ndjet_mu.Fill(event.JetPhi[1])
                    JetBTagger_hist_2ndjet_mu.Fill(event.JetBTagger[1])
                    JetIsBtagged_hist_2ndjet_mu.Fill(event.JetIsBtagged[1])
                    JetIsBtaggedWithSF_hist_2ndjet_mu.Fill(event.JetIsBtaggedWithSF[1])
            
                    if math.fabs(event.JetEta[1]) > 3.0 : 
                        JetPt_hist_2ndjet_fwdeta_mu.Fill(event.JetPt[1])
                        JetEta_hist_2ndjet_fwdeta_mu.Fill(event.JetEta[1])
                        JetPhi_hist_2ndjet_fwdeta_mu.Fill(event.JetPhi[1])
                        JetBTagger_hist_2ndjet_fwdeta_mu.Fill(event.JetBTagger[1])
                        JetIsBtagged_hist_2ndjet_fwdeta_mu.Fill(event.JetIsBtagged[1])
                        JetIsBtaggedWithSF_hist_2ndjet_fwdeta_mu.Fill(event.JetIsBtaggedWithSF[1])
                        
            
            PFMET_hist_mu.Fill(event.PFMET)

        
    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_DATA = TFile.Open("jetsMetDistrib_DATA_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_DATA.cd()

    # Zee histos
    # nCleanedJets_hist_ele.Write()
    nCleanedJetsPt30_hist_ele.Write()
    # nCleanedJetsPt30BTagged_hist_ele.Write()
    # nCleanedJetsPt30BTagged_bTagSF_hist_ele.Write()
    
    JetPt_hist_ele.Write()
    JetPt_hist_1stJet_ele.Write()
    JetPt_hist_2ndjet_ele.Write()
    JetPt_hist_1stJet_fwdeta_ele.Write()
    JetPt_hist_2ndjet_fwdeta_ele.Write()
     
    JetEta_hist_ele.Write()
    JetEta_hist_1stjet_ele.Write()
    JetEta_hist_2ndjet_ele.Write()
    JetEta_hist_1stjet_fwdeta_ele.Write()
    JetEta_hist_2ndjet_fwdeta_ele.Write()
                            
    JetPhi_hist_ele.Write()
    JetPhi_hist_1stjet_ele.Write()
    JetPhi_hist_2ndjet_ele.Write()
    JetPhi_hist_1stjet_fwdeta_ele.Write()
    JetPhi_hist_2ndjet_fwdeta_ele.Write()
                            
    JetBTagger_hist_ele.Write()
    JetBTagger_hist_1stjet_ele.Write()
    JetBTagger_hist_2ndjet_ele.Write()
    JetBTagger_hist_1stjet_fwdeta_ele.Write()
    JetBTagger_hist_2ndjet_fwdeta_ele.Write()
                            
    JetIsBtagged_hist_ele.Write()
    JetIsBtagged_hist_1stjet_ele.Write()
    JetIsBtagged_hist_2ndjet_ele.Write()
    JetIsBtagged_hist_1stjet_fwdeta_ele.Write()
    JetIsBtagged_hist_2ndjet_fwdeta_ele.Write()
                            
    JetIsBtaggedWithSF_hist_ele.Write()   
    JetIsBtaggedWithSF_hist_1stjet_ele.Write()
    JetIsBtaggedWithSF_hist_2ndjet_ele.Write()
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_ele.Write()
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_ele.Write()

    PFMET_hist_ele.Write()

    # Zmumu histos
    # nCleanedJets_hist_mu.Write()
    nCleanedJetsPt30_hist_mu.Write()
    # nCleanedJetsPt30BTagged_hist_mu.Write()
    # nCleanedJetsPt30BTagged_bTagSF_hist_mu.Write()
    
    JetPt_hist_mu.Write()
    JetPt_hist_1stJet_mu.Write()
    JetPt_hist_2ndjet_mu.Write()
    JetPt_hist_1stJet_fwdeta_mu.Write()
    JetPt_hist_2ndjet_fwdeta_mu.Write()
     
    JetEta_hist_mu.Write()
    JetEta_hist_1stjet_mu.Write()
    JetEta_hist_2ndjet_mu.Write()
    JetEta_hist_1stjet_fwdeta_mu.Write()
    JetEta_hist_2ndjet_fwdeta_mu.Write()
                            
    JetPhi_hist_mu.Write()
    JetPhi_hist_1stjet_mu.Write()
    JetPhi_hist_2ndjet_mu.Write()
    JetPhi_hist_1stjet_fwdeta_mu.Write()
    JetPhi_hist_2ndjet_fwdeta_mu.Write()
                            
    JetBTagger_hist_mu.Write()
    JetBTagger_hist_1stjet_mu.Write()
    JetBTagger_hist_2ndjet_mu.Write()
    JetBTagger_hist_1stjet_fwdeta_mu.Write()
    JetBTagger_hist_2ndjet_fwdeta_mu.Write()
                            
    JetIsBtagged_hist_mu.Write()
    JetIsBtagged_hist_1stjet_mu.Write()
    JetIsBtagged_hist_2ndjet_mu.Write()
    JetIsBtagged_hist_1stjet_fwdeta_mu.Write()
    JetIsBtagged_hist_2ndjet_fwdeta_mu.Write()
                            
    JetIsBtaggedWithSF_hist_mu.Write()   
    JetIsBtaggedWithSF_hist_1stjet_mu.Write()
    JetIsBtaggedWithSF_hist_2ndjet_mu.Write()
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_mu.Write()
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_mu.Write()

    PFMET_hist_mu.Write()


    outFile_DATA.Close()
    print "DATA histo file created!"


# ********************
#  do MC DY histos
# ********************
if(redoMCDYHistos) :

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    # define data histograms Z->ee

    # nCleanedJets_hist_MC_DY_ele                   = TH1F('nCleanedJets_MC_DY_ele',                   'nCleanedJets_MC_DY_ele',                   30, 0, 30) 
    nCleanedJetsPt30_hist_MC_DY_ele               = TH1F('nCleanedJetsPt30_MC_DY_ele',               'nCleanedJetsPt30_MC_DY_ele',               15, 0, 15) 
    # nCleanedJetsPt30BTagged_hist_MC_DY_ele        = TH1F('nCleanedJetsPt30BTagged_MC_DY_ele',        'nCleanedJetsPt30BTagged_MC_DY_ele',        5,  0, 5)
    # nCleanedJetsPt30BTagged_bTagSF_hist_MC_DY_ele = TH1F('nCleanedJetsPt30BTagged_bTagSF_MC_DY_ele', 'nCleanedJetsPt30BTagged_bTagSF_MC_DY_ele', 5,  0, 5)

    JetPt_hist_MC_DY_ele               = TH1F('JetPt_inclusive_MC_DY_ele',           'JetPt_inclusive_MC_DY_ele',           120, 0, 1200)
    JetPt_hist_1stJet_MC_DY_ele        = TH1F('JetPt_leadingJet_MC_DY_ele',          'JetPt_leadingJet_MC_DY_ele',          120, 0, 1200)
    JetPt_hist_2ndjet_MC_DY_ele        = TH1F('JetPt_subLeadingJet_MC_DY_ele',       'JetPt_subLeadingJet_MC_DY_ele',       120, 0, 1200)
    JetPt_hist_1stJet_fwdeta_MC_DY_ele = TH1F('JetPt_leadingJet_fwdeta_MC_DY_ele',   'JetPt_leadingJet_fwdeta_MC_DY_ele',   120, 0, 1200)
    JetPt_hist_2ndjet_fwdeta_MC_DY_ele = TH1F('JetPt_subLeadingJet_fwdeta_MC_DY_ele','JetPt_subLeadingJet_fwdeta_MC_DY_ele',120, 0, 1200)

    JetEta_hist_MC_DY_ele               = TH1F('JetEta_inclusive_MC_DY_ele',           'JetEta_inclusive_MC_DY_ele',           50, -5, 5)
    JetEta_hist_1stjet_MC_DY_ele        = TH1F('JetEta_leadingJet_MC_DY_ele',          'JetEta_leadingJet_MC_DY_ele',          50, -5, 5)
    JetEta_hist_2ndjet_MC_DY_ele        = TH1F('JetEta_subLeadingJet_MC_DY_ele',       'JetEta_subLeadingJet_MC_DY_ele',       50, -5, 5)
    JetEta_hist_1stjet_fwdeta_MC_DY_ele = TH1F('JetEta_leadingJet_fwdeta_MC_DY_ele',   'JetEta_leadingJet_fwdeta_MC_DY_ele',   50, -5, 5)
    JetEta_hist_2ndjet_fwdeta_MC_DY_ele = TH1F('JetEta_subLeadingJet_fwdeta_MC_DY_ele','JetEta_subLeadingJet_fwdeta_MC_DY_ele',50, -5, 5)

    JetPhi_hist_MC_DY_ele               = TH1F('JetPhi_inclusive_MC_DY_ele',           'JetPhi_inclusive_MC_DY_ele',           33, -3.3, 3.3)
    JetPhi_hist_1stjet_MC_DY_ele        = TH1F('JetPhi_leadingJet_MC_DY_ele',          'JetPhi_leadingJet_MC_DY_ele',          33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_MC_DY_ele        = TH1F('JetPhi_subLeadingJet_MC_DY_ele',       'JetPhi_subLeadingJet_MC_DY_ele',       33, -3.3, 3.3) 
    JetPhi_hist_1stjet_fwdeta_MC_DY_ele = TH1F('JetPhi_leadingJet_fwdeta_MC_DY_ele',   'JetPhi_leadingJet_fwdeta_MC_DY_ele',   33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_fwdeta_MC_DY_ele = TH1F('JetPhi_subLeadingJet_fwdeta_MC_DY_ele','JetPhi_subLeadingJet_fwdeta_MC_DY_ele',33, -3.3, 3.3) 

    JetBTagger_hist_MC_DY_ele               = TH1F('JetBTagger_inclusive_MC_DY_ele',           'JetBTagger_inclusive_MC_DY_ele',           11, 0, 1.1)
    JetBTagger_hist_1stjet_MC_DY_ele        = TH1F('JetBTagger_leadingJet_MC_DY_ele',          'JetBTagger_leadingJet_MC_DY_ele',          11, 0, 1.1)
    JetBTagger_hist_2ndjet_MC_DY_ele        = TH1F('JetBTagger_subLeadingJet_MC_DY_ele',       'JetBTagger_subLeadingJet_MC_DY_ele',       11, 0, 1.1)
    JetBTagger_hist_1stjet_fwdeta_MC_DY_ele = TH1F('JetBTagger_leadingJet_fwdeta_MC_DY_ele',   'JetBTagger_leadingJet_fwdeta_MC_DY_ele',   11, 0, 1.1)
    JetBTagger_hist_2ndjet_fwdeta_MC_DY_ele = TH1F('JetBTagger_subLeadingJet_fwdeta_MC_DY_ele','JetBTagger_subLeadingJet_fwdeta_MC_DY_ele',11, 0, 1.1)

    JetIsBtagged_hist_MC_DY_ele               = TH1F('JetIsBtagged_inclusive_MC_DY_ele',           'JetIsBtagged_inclusive_MC_DY_ele',           11, 0, 1.1)
    JetIsBtagged_hist_1stjet_MC_DY_ele        = TH1F('JetIsBtagged_leadingJet_MC_DY_ele',          'JetIsBtagged_leadingJet_MC_DY_ele',          11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_MC_DY_ele        = TH1F('JetIsBtagged_subLeadingJet_MC_DY_ele',       'JetIsBtagged_subLeadingJet_MC_DY_ele',       11, 0, 1.1)
    JetIsBtagged_hist_1stjet_fwdeta_MC_DY_ele = TH1F('JetIsBtagged_leadingJet_fwdeta_MC_DY_ele',   'JetIsBtagged_leadingJet_fwdeta_MC_DY_ele',   11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_fwdeta_MC_DY_ele = TH1F('JetIsBtagged_subLeadingJet_fwdeta_MC_DY_ele','JetIsBtagged_subLeadingJet_fwdeta_MC_DY_ele',11, 0, 1.1)
    
    JetIsBtaggedWithSF_hist_MC_DY_ele               = TH1F('JetIsBtaggedWithSF_inclusive_MC_DY_ele',           'JetIsBtaggedWithSF_inclusive_MC_DY_ele',           11, 0, 1.1)
    JetIsBtaggedWithSF_hist_1stjet_MC_DY_ele        = TH1F('JetIsBtaggedWithSF_leadingJet_MC_DY_ele',          'JetIsBtaggedWithSF_leadingJet_MC_DY_ele',          11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_MC_DY_ele        = TH1F('JetIsBtaggedWithSF_subLeadingJet_MC_DY_ele',       'JetIsBtaggedWithSF_subLeadingJet_MC_DY_ele',       11, 0, 1.1) 
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_DY_ele = TH1F('JetIsBtaggedWithSF_leadingJet_fwdeta_MC_DY_ele',   'JetIsBtaggedWithSF_leadingJet_fwdeta_MC_DY_ele',   11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_DY_ele = TH1F('JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_DY_ele','JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_DY_ele',11, 0, 1.1)

    PFMET_hist_MC_DY_ele = TH1F('PFMET_MC_DY_ele','PFMET_MC_DY_ele',125,0,250)


    # define data histograms Z->mumu

    # nCleanedJets_hist_MC_DY_mu                   = TH1F('nCleanedJets_MC_DY_mu',                   'nCleanedJets_MC_DY_mu',                   30, 0, 30) 
    nCleanedJetsPt30_hist_MC_DY_mu               = TH1F('nCleanedJetsPt30_MC_DY_mu',               'nCleanedJetsPt30_MC_DY_mu',               15, 0, 15) 
    # nCleanedJetsPt30BTagged_hist_MC_DY_mu        = TH1F('nCleanedJetsPt30BTagged_MC_DY_mu',        'nCleanedJetsPt30BTagged_MC_DY_mu',        5,  0, 5)
    # nCleanedJetsPt30BTagged_bTagSF_hist_MC_DY_mu = TH1F('nCleanedJetsPt30BTagged_bTagSF_MC_DY_mu', 'nCleanedJetsPt30BTagged_bTagSF_MC_DY_mu', 5,  0, 5)

    JetPt_hist_MC_DY_mu               = TH1F('JetPt_inclusive_MC_DY_mu',           'JetPt_inclusive_MC_DY_mu',           120, 0, 1200)
    JetPt_hist_1stJet_MC_DY_mu        = TH1F('JetPt_leadingJet_MC_DY_mu',          'JetPt_leadingJet_MC_DY_mu',          120, 0, 1200)
    JetPt_hist_2ndjet_MC_DY_mu        = TH1F('JetPt_subLeadingJet_MC_DY_mu',       'JetPt_subLeadingJet_MC_DY_mu',       120, 0, 1200)
    JetPt_hist_1stJet_fwdeta_MC_DY_mu = TH1F('JetPt_leadingJet_fwdeta_MC_DY_mu',   'JetPt_leadingJet_fwdeta_MC_DY_mu',   120, 0, 1200)
    JetPt_hist_2ndjet_fwdeta_MC_DY_mu = TH1F('JetPt_subLeadingJet_fwdeta_MC_DY_mu','JetPt_subLeadingJet_fwdeta_MC_DY_mu',120, 0, 1200)

    JetEta_hist_MC_DY_mu               = TH1F('JetEta_inclusive_MC_DY_mu',           'JetEta_inclusive_MC_DY_mu',           50, -5, 5)
    JetEta_hist_1stjet_MC_DY_mu        = TH1F('JetEta_leadingJet_MC_DY_mu',          'JetEta_leadingJet_MC_DY_mu',          50, -5, 5)
    JetEta_hist_2ndjet_MC_DY_mu        = TH1F('JetEta_subLeadingJet_MC_DY_mu',       'JetEta_subLeadingJet_MC_DY_mu',       50, -5, 5)
    JetEta_hist_1stjet_fwdeta_MC_DY_mu = TH1F('JetEta_leadingJet_fwdeta_MC_DY_mu',   'JetEta_leadingJet_fwdeta_MC_DY_mu',   50, -5, 5)
    JetEta_hist_2ndjet_fwdeta_MC_DY_mu = TH1F('JetEta_subLeadingJet_fwdeta_MC_DY_mu','JetEta_subLeadingJet_fwdeta_MC_DY_mu',50, -5, 5)

    JetPhi_hist_MC_DY_mu               = TH1F('JetPhi_inclusive_MC_DY_mu',           'JetPhi_inclusive_MC_DY_mu',           33, -3.3, 3.3)
    JetPhi_hist_1stjet_MC_DY_mu        = TH1F('JetPhi_leadingJet_MC_DY_mu',          'JetPhi_leadingJet_MC_DY_mu',          33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_MC_DY_mu        = TH1F('JetPhi_subLeadingJet_MC_DY_mu',       'JetPhi_subLeadingJet_MC_DY_mu',       33, -3.3, 3.3) 
    JetPhi_hist_1stjet_fwdeta_MC_DY_mu = TH1F('JetPhi_leadingJet_fwdeta_MC_DY_mu',   'JetPhi_leadingJet_fwdeta_MC_DY_mu',   33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_fwdeta_MC_DY_mu = TH1F('JetPhi_subLeadingJet_fwdeta_MC_DY_mu','JetPhi_subLeadingJet_fwdeta_MC_DY_mu',33, -3.3, 3.3) 

    JetBTagger_hist_MC_DY_mu               = TH1F('JetBTagger_inclusive_MC_DY_mu',           'JetBTagger_inclusive_MC_DY_mu',           11, 0, 1.1)
    JetBTagger_hist_1stjet_MC_DY_mu        = TH1F('JetBTagger_leadingJet_MC_DY_mu',          'JetBTagger_leadingJet_MC_DY_mu',          11, 0, 1.1)
    JetBTagger_hist_2ndjet_MC_DY_mu        = TH1F('JetBTagger_subLeadingJet_MC_DY_mu',       'JetBTagger_subLeadingJet_MC_DY_mu',       11, 0, 1.1)
    JetBTagger_hist_1stjet_fwdeta_MC_DY_mu = TH1F('JetBTagger_leadingJet_fwdeta_MC_DY_mu',   'JetBTagger_leadingJet_fwdeta_MC_DY_mu',   11, 0, 1.1)
    JetBTagger_hist_2ndjet_fwdeta_MC_DY_mu = TH1F('JetBTagger_subLeadingJet_fwdeta_MC_DY_mu','JetBTagger_subLeadingJet_fwdeta_MC_DY_mu',11, 0, 1.1)

    JetIsBtagged_hist_MC_DY_mu               = TH1F('JetIsBtagged_inclusive_MC_DY_mu',           'JetIsBtagged_inclusive_MC_DY_mu',           11, 0, 1.1)
    JetIsBtagged_hist_1stjet_MC_DY_mu        = TH1F('JetIsBtagged_leadingJet_MC_DY_mu',          'JetIsBtagged_leadingJet_MC_DY_mu',          11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_MC_DY_mu        = TH1F('JetIsBtagged_subLeadingJet_MC_DY_mu',       'JetIsBtagged_subLeadingJet_MC_DY_mu',       11, 0, 1.1)
    JetIsBtagged_hist_1stjet_fwdeta_MC_DY_mu = TH1F('JetIsBtagged_leadingJet_fwdeta_MC_DY_mu',   'JetIsBtagged_leadingJet_fwdeta_MC_DY_mu',   11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_fwdeta_MC_DY_mu = TH1F('JetIsBtagged_subLeadingJet_fwdeta_MC_DY_mu','JetIsBtagged_subLeadingJet_fwdeta_MC_DY_mu',11, 0, 1.1)
    
    JetIsBtaggedWithSF_hist_MC_DY_mu               = TH1F('JetIsBtaggedWithSF_inclusive_MC_DY_mu',           'JetIsBtaggedWithSF_inclusive_MC_DY_mu',           11, 0, 1.1)
    JetIsBtaggedWithSF_hist_1stjet_MC_DY_mu        = TH1F('JetIsBtaggedWithSF_leadingJet_MC_DY_mu',          'JetIsBtaggedWithSF_leadingJet_MC_DY_mu',          11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_MC_DY_mu        = TH1F('JetIsBtaggedWithSF_subLeadingJet_MC_DY_mu',       'JetIsBtaggedWithSF_subLeadingJet_MC_DY_mu',       11, 0, 1.1) 
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_DY_mu = TH1F('JetIsBtaggedWithSF_leadingJet_fwdeta_MC_DY_mu',   'JetIsBtaggedWithSF_leadingJet_fwdeta_MC_DY_mu',   11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_DY_mu = TH1F('JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_DY_mu','JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_DY_mu',11, 0, 1.1)

    PFMET_hist_MC_DY_mu = TH1F('PFMET_MC_DY_mu','PFMET_MC_DY_mu',125,0,250)




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

            # nCleanedJets_hist_MC_DY_ele.Fill(event.nCleanedJets,weight)
            nCleanedJetsPt30_hist_MC_DY_ele.Fill(event.nCleanedJetsPt30,weight)
            # nCleanedJetsPt30BTagged_hist_MC_DY_ele.Fill(event.nCleanedJetsPt30BTagged,weight)
            # nCleanedJetsPt30BTagged_bTagSF_hist_MC_DY_ele.Fill(event.nCleanedJetsPt30BTagged_bTagSF,weight)
            
            for i in range(len(event.JetPt)) :
                JetPt_hist_MC_DY_ele.Fill(event.JetPt[i],weight)
                      
            for i in range(len(event.JetEta)) :  
                JetEta_hist_MC_DY_ele.Fill(event.JetEta[i],weight)
                                            
            for i in range(len(event.JetPhi)) :
                JetPhi_hist_MC_DY_ele.Fill(event.JetPhi[i],weight)
              

            if event.nCleanedJetsPt30 > 0 :

                for i in range(len(event.JetBTagger)) :
                    JetBTagger_hist_MC_DY_ele.Fill(event.JetBTagger[i],weight)
                                            
                for i in range(len(event.JetIsBtagged)) :
                    JetIsBtagged_hist_MC_DY_ele.Fill(event.JetIsBtagged[i],weight)
            
                for i in range(len(event.JetIsBtaggedWithSF)) :                
                    JetIsBtaggedWithSF_hist_MC_DY_ele.Fill(event.JetIsBtaggedWithSF[i],weight)     
            
            
                JetPt_hist_1stJet_MC_DY_ele.Fill(event.JetPt[0],weight)
                JetEta_hist_1stjet_MC_DY_ele.Fill(event.JetEta[0],weight)
                JetPhi_hist_1stjet_MC_DY_ele.Fill(event.JetPhi[0],weight)
                JetBTagger_hist_1stjet_MC_DY_ele.Fill(event.JetBTagger[0],weight)
                JetIsBtagged_hist_1stjet_MC_DY_ele.Fill(event.JetIsBtagged[0],weight)
                JetIsBtaggedWithSF_hist_1stjet_MC_DY_ele.Fill(event.JetIsBtaggedWithSF[0],weight)
            
                if math.fabs(event.JetEta[0]) > 3.0 : 
                    JetPt_hist_1stJet_fwdeta_MC_DY_ele.Fill(event.JetPt[0],weight)
                    JetEta_hist_1stjet_fwdeta_MC_DY_ele.Fill(event.JetEta[0],weight)
                    JetPhi_hist_1stjet_fwdeta_MC_DY_ele.Fill(event.JetPhi[0],weight)
                    JetBTagger_hist_1stjet_fwdeta_MC_DY_ele.Fill(event.JetBTagger[0],weight)
                    JetIsBtagged_hist_1stjet_fwdeta_MC_DY_ele.Fill(event.JetIsBtagged[0],weight)
                    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_DY_ele.Fill(event.JetIsBtaggedWithSF[0],weight)
            
            
                if event.nCleanedJetsPt30 > 1 :
                    JetPt_hist_2ndjet_MC_DY_ele.Fill(event.JetPt[1],weight)
                    JetEta_hist_2ndjet_MC_DY_ele.Fill(event.JetEta[1],weight)
                    JetPhi_hist_2ndjet_MC_DY_ele.Fill(event.JetPhi[1],weight)
                    JetBTagger_hist_2ndjet_MC_DY_ele.Fill(event.JetBTagger[1],weight)
                    JetIsBtagged_hist_2ndjet_MC_DY_ele.Fill(event.JetIsBtagged[1],weight)
                    JetIsBtaggedWithSF_hist_2ndjet_MC_DY_ele.Fill(event.JetIsBtaggedWithSF[1],weight)
            
                    if math.fabs(event.JetEta[1]) > 3.0 : 
                        JetPt_hist_2ndjet_fwdeta_MC_DY_ele.Fill(event.JetPt[1],weight)
                        JetEta_hist_2ndjet_fwdeta_MC_DY_ele.Fill(event.JetEta[1],weight)
                        JetPhi_hist_2ndjet_fwdeta_MC_DY_ele.Fill(event.JetPhi[1],weight)
                        JetBTagger_hist_2ndjet_fwdeta_MC_DY_ele.Fill(event.JetBTagger[1],weight)
                        JetIsBtagged_hist_2ndjet_fwdeta_MC_DY_ele.Fill(event.JetIsBtagged[1],weight)
                        JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_DY_ele.Fill(event.JetIsBtaggedWithSF[1],weight)
                        
            
            PFMET_hist_MC_DY_ele.Fill(event.PFMET,weight)



        # Z->mumu histos
        if(int(math.fabs(event.LepLepId[0])) == 13 ):

            # nCleanedJets_hist_MC_DY_mu.Fill(event.nCleanedJets,weight)
            nCleanedJetsPt30_hist_MC_DY_mu.Fill(event.nCleanedJetsPt30,weight)
            # nCleanedJetsPt30BTagged_hist_MC_DY_mu.Fill(event.nCleanedJetsPt30BTagged,weight)
            # nCleanedJetsPt30BTagged_bTagSF_hist_MC_DY_mu.Fill(event.nCleanedJetsPt30BTagged_bTagSF,weight)
            
            for i in range(len(event.JetPt)) :
                JetPt_hist_MC_DY_mu.Fill(event.JetPt[i],weight)
                      
            for i in range(len(event.JetEta)) :  
                JetEta_hist_MC_DY_mu.Fill(event.JetEta[i],weight)
                                            
            for i in range(len(event.JetPhi)) :
                JetPhi_hist_MC_DY_mu.Fill(event.JetPhi[i],weight)

              
            if event.nCleanedJetsPt30 > 0 :

                for i in range(len(event.JetBTagger)) :
                    JetBTagger_hist_MC_DY_mu.Fill(event.JetBTagger[i],weight)
                                            
                for i in range(len(event.JetIsBtagged)) :
                    JetIsBtagged_hist_MC_DY_mu.Fill(event.JetIsBtagged[i],weight)
            
                for i in range(len(event.JetIsBtaggedWithSF)) :                
                    JetIsBtaggedWithSF_hist_MC_DY_mu.Fill(event.JetIsBtaggedWithSF[i],weight)     
            
            
                JetPt_hist_1stJet_MC_DY_mu.Fill(event.JetPt[0],weight)
                JetEta_hist_1stjet_MC_DY_mu.Fill(event.JetEta[0],weight)
                JetPhi_hist_1stjet_MC_DY_mu.Fill(event.JetPhi[0],weight)
                JetBTagger_hist_1stjet_MC_DY_mu.Fill(event.JetBTagger[0],weight)
                JetIsBtagged_hist_1stjet_MC_DY_mu.Fill(event.JetIsBtagged[0],weight)
                JetIsBtaggedWithSF_hist_1stjet_MC_DY_mu.Fill(event.JetIsBtaggedWithSF[0],weight)
            
                if math.fabs(event.JetEta[0]) > 3.0 : 
                    JetPt_hist_1stJet_fwdeta_MC_DY_mu.Fill(event.JetPt[0],weight)
                    JetEta_hist_1stjet_fwdeta_MC_DY_mu.Fill(event.JetEta[0],weight)
                    JetPhi_hist_1stjet_fwdeta_MC_DY_mu.Fill(event.JetPhi[0],weight)
                    JetBTagger_hist_1stjet_fwdeta_MC_DY_mu.Fill(event.JetBTagger[0],weight)
                    JetIsBtagged_hist_1stjet_fwdeta_MC_DY_mu.Fill(event.JetIsBtagged[0],weight)
                    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_DY_mu.Fill(event.JetIsBtaggedWithSF[0],weight)
            
            
                if event.nCleanedJetsPt30 > 1 :
                    JetPt_hist_2ndjet_MC_DY_mu.Fill(event.JetPt[1],weight)
                    JetEta_hist_2ndjet_MC_DY_mu.Fill(event.JetEta[1],weight)
                    JetPhi_hist_2ndjet_MC_DY_mu.Fill(event.JetPhi[1],weight)
                    JetBTagger_hist_2ndjet_MC_DY_mu.Fill(event.JetBTagger[1],weight)
                    JetIsBtagged_hist_2ndjet_MC_DY_mu.Fill(event.JetIsBtagged[1],weight)
                    JetIsBtaggedWithSF_hist_2ndjet_MC_DY_mu.Fill(event.JetIsBtaggedWithSF[1],weight)
            
                    if math.fabs(event.JetEta[1]) > 3.0 : 
                        JetPt_hist_2ndjet_fwdeta_MC_DY_mu.Fill(event.JetPt[1],weight)
                        JetEta_hist_2ndjet_fwdeta_MC_DY_mu.Fill(event.JetEta[1],weight)
                        JetPhi_hist_2ndjet_fwdeta_MC_DY_mu.Fill(event.JetPhi[1],weight)
                        JetBTagger_hist_2ndjet_fwdeta_MC_DY_mu.Fill(event.JetBTagger[1],weight)
                        JetIsBtagged_hist_2ndjet_fwdeta_MC_DY_mu.Fill(event.JetIsBtagged[1],weight)
                        JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_DY_mu.Fill(event.JetIsBtaggedWithSF[1],weight)
                        
            
            PFMET_hist_MC_DY_mu.Fill(event.PFMET,weight)


    

        
    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_MCDY = TFile.Open("jetsMetDistrib_MC_DY_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_MCDY.cd()

    # Zee histos
    # nCleanedJets_hist_MC_DY_ele.Write()
    nCleanedJetsPt30_hist_MC_DY_ele.Write()
    # nCleanedJetsPt30BTagged_hist_MC_DY_ele.Write()
    # nCleanedJetsPt30BTagged_bTagSF_hist_MC_DY_ele.Write()
    
    JetPt_hist_MC_DY_ele.Write()
    JetPt_hist_1stJet_MC_DY_ele.Write()
    JetPt_hist_2ndjet_MC_DY_ele.Write()
    JetPt_hist_1stJet_fwdeta_MC_DY_ele.Write()
    JetPt_hist_2ndjet_fwdeta_MC_DY_ele.Write()
                            
    JetEta_hist_MC_DY_ele.Write()
    JetEta_hist_1stjet_MC_DY_ele.Write()
    JetEta_hist_2ndjet_MC_DY_ele.Write()
    JetEta_hist_1stjet_fwdeta_MC_DY_ele.Write()
    JetEta_hist_2ndjet_fwdeta_MC_DY_ele.Write()
                            
    JetPhi_hist_MC_DY_ele.Write()
    JetPhi_hist_1stjet_MC_DY_ele.Write()
    JetPhi_hist_2ndjet_MC_DY_ele.Write()
    JetPhi_hist_1stjet_fwdeta_MC_DY_ele.Write()
    JetPhi_hist_2ndjet_fwdeta_MC_DY_ele.Write()
                            
    JetBTagger_hist_MC_DY_ele.Write()
    JetBTagger_hist_1stjet_MC_DY_ele.Write()
    JetBTagger_hist_2ndjet_MC_DY_ele.Write()
    JetBTagger_hist_1stjet_fwdeta_MC_DY_ele.Write()
    JetBTagger_hist_2ndjet_fwdeta_MC_DY_ele.Write()
                            
    JetIsBtagged_hist_MC_DY_ele.Write()
    JetIsBtagged_hist_1stjet_MC_DY_ele.Write()
    JetIsBtagged_hist_2ndjet_MC_DY_ele.Write()
    JetIsBtagged_hist_1stjet_fwdeta_MC_DY_ele.Write()
    JetIsBtagged_hist_2ndjet_fwdeta_MC_DY_ele.Write()
                            
    JetIsBtaggedWithSF_hist_MC_DY_ele.Write()   
    JetIsBtaggedWithSF_hist_1stjet_MC_DY_ele.Write()
    JetIsBtaggedWithSF_hist_2ndjet_MC_DY_ele.Write()
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_DY_ele.Write()
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_DY_ele.Write()

    PFMET_hist_MC_DY_ele.Write()

    
    # Zmumu histos
    # nCleanedJets_hist_MC_DY_mu.Write()
    nCleanedJetsPt30_hist_MC_DY_mu.Write()
    # nCleanedJetsPt30BTagged_hist_MC_DY_mu.Write()
    # nCleanedJetsPt30BTagged_bTagSF_hist_MC_DY_mu.Write()
    
    JetPt_hist_MC_DY_mu.Write()
    JetPt_hist_1stJet_MC_DY_mu.Write()
    JetPt_hist_2ndjet_MC_DY_mu.Write()
    JetPt_hist_1stJet_fwdeta_MC_DY_mu.Write()
    JetPt_hist_2ndjet_fwdeta_MC_DY_mu.Write()
                            
    JetEta_hist_MC_DY_mu.Write()
    JetEta_hist_1stjet_MC_DY_mu.Write()
    JetEta_hist_2ndjet_MC_DY_mu.Write()
    JetEta_hist_1stjet_fwdeta_MC_DY_mu.Write()
    JetEta_hist_2ndjet_fwdeta_MC_DY_mu.Write()
                            
    JetPhi_hist_MC_DY_mu.Write()
    JetPhi_hist_1stjet_MC_DY_mu.Write()
    JetPhi_hist_2ndjet_MC_DY_mu.Write()
    JetPhi_hist_1stjet_fwdeta_MC_DY_mu.Write()
    JetPhi_hist_2ndjet_fwdeta_MC_DY_mu.Write()
                            
    JetBTagger_hist_MC_DY_mu.Write()
    JetBTagger_hist_1stjet_MC_DY_mu.Write()
    JetBTagger_hist_2ndjet_MC_DY_mu.Write()
    JetBTagger_hist_1stjet_fwdeta_MC_DY_mu.Write()
    JetBTagger_hist_2ndjet_fwdeta_MC_DY_mu.Write()
                            
    JetIsBtagged_hist_MC_DY_mu.Write()
    JetIsBtagged_hist_1stjet_MC_DY_mu.Write()
    JetIsBtagged_hist_2ndjet_MC_DY_mu.Write()
    JetIsBtagged_hist_1stjet_fwdeta_MC_DY_mu.Write()
    JetIsBtagged_hist_2ndjet_fwdeta_MC_DY_mu.Write()
                            
    JetIsBtaggedWithSF_hist_MC_DY_mu.Write()   
    JetIsBtaggedWithSF_hist_1stjet_MC_DY_mu.Write()
    JetIsBtaggedWithSF_hist_2ndjet_MC_DY_mu.Write()
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_DY_mu.Write()
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_DY_mu.Write()

    PFMET_hist_MC_DY_mu.Write()


    outFile_MCDY.Close()
    print "MC DY histo file created!"
# ********************


# ********************
#  do MC TTbar histos
# ********************
if(redoMCTTbarHistos) :

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    # define data histograms Z->ee

    # nCleanedJets_hist_MC_TTbar_ele                   = TH1F('nCleanedJets_MC_TTbar_ele',                   'nCleanedJets_MC_TTbar_ele',                   30, 0, 30) 
    nCleanedJetsPt30_hist_MC_TTbar_ele               = TH1F('nCleanedJetsPt30_MC_TTbar_ele',               'nCleanedJetsPt30_MC_TTbar_ele',               15, 0, 15) 
    # nCleanedJetsPt30BTagged_hist_MC_TTbar_ele        = TH1F('nCleanedJetsPt30BTagged_MC_TTbar_ele',        'nCleanedJetsPt30BTagged_MC_TTbar_ele',        5,  0, 5)
    # nCleanedJetsPt30BTagged_bTagSF_hist_MC_TTbar_ele = TH1F('nCleanedJetsPt30BTagged_bTagSF_MC_TTbar_ele', 'nCleanedJetsPt30BTagged_bTagSF_MC_TTbar_ele', 5,  0, 5)

    JetPt_hist_MC_TTbar_ele               = TH1F('JetPt_inclusive_MC_TTbar_ele',           'JetPt_inclusive_MC_TTbar_ele',           120, 0, 1200)
    JetPt_hist_1stJet_MC_TTbar_ele        = TH1F('JetPt_leadingJet_MC_TTbar_ele',          'JetPt_leadingJet_MC_TTbar_ele',          120, 0, 1200)
    JetPt_hist_2ndjet_MC_TTbar_ele        = TH1F('JetPt_subLeadingJet_MC_TTbar_ele',       'JetPt_subLeadingJet_MC_TTbar_ele',       120, 0, 1200)
    JetPt_hist_1stJet_fwdeta_MC_TTbar_ele = TH1F('JetPt_leadingJet_fwdeta_MC_TTbar_ele',   'JetPt_leadingJet_fwdeta_MC_TTbar_ele',   120, 0, 1200)
    JetPt_hist_2ndjet_fwdeta_MC_TTbar_ele = TH1F('JetPt_subLeadingJet_fwdeta_MC_TTbar_ele','JetPt_subLeadingJet_fwdeta_MC_TTbar_ele',120, 0, 1200)

    JetEta_hist_MC_TTbar_ele               = TH1F('JetEta_inclusive_MC_TTbar_ele',           'JetEta_inclusive_MC_TTbar_ele',           50, -5, 5)
    JetEta_hist_1stjet_MC_TTbar_ele        = TH1F('JetEta_leadingJet_MC_TTbar_ele',          'JetEta_leadingJet_MC_TTbar_ele',          50, -5, 5)
    JetEta_hist_2ndjet_MC_TTbar_ele        = TH1F('JetEta_subLeadingJet_MC_TTbar_ele',       'JetEta_subLeadingJet_MC_TTbar_ele',       50, -5, 5)
    JetEta_hist_1stjet_fwdeta_MC_TTbar_ele = TH1F('JetEta_leadingJet_fwdeta_MC_TTbar_ele',   'JetEta_leadingJet_fwdeta_MC_TTbar_ele',   50, -5, 5)
    JetEta_hist_2ndjet_fwdeta_MC_TTbar_ele = TH1F('JetEta_subLeadingJet_fwdeta_MC_TTbar_ele','JetEta_subLeadingJet_fwdeta_MC_TTbar_ele',50, -5, 5)

    JetPhi_hist_MC_TTbar_ele               = TH1F('JetPhi_inclusive_MC_TTbar_ele',           'JetPhi_inclusive_MC_TTbar_ele',           33, -3.3, 3.3)
    JetPhi_hist_1stjet_MC_TTbar_ele        = TH1F('JetPhi_leadingJet_MC_TTbar_ele',          'JetPhi_leadingJet_MC_TTbar_ele',          33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_MC_TTbar_ele        = TH1F('JetPhi_subLeadingJet_MC_TTbar_ele',       'JetPhi_subLeadingJet_MC_TTbar_ele',       33, -3.3, 3.3) 
    JetPhi_hist_1stjet_fwdeta_MC_TTbar_ele = TH1F('JetPhi_leadingJet_fwdeta_MC_TTbar_ele',   'JetPhi_leadingJet_fwdeta_MC_TTbar_ele',   33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_fwdeta_MC_TTbar_ele = TH1F('JetPhi_subLeadingJet_fwdeta_MC_TTbar_ele','JetPhi_subLeadingJet_fwdeta_MC_TTbar_ele',33, -3.3, 3.3) 

    JetBTagger_hist_MC_TTbar_ele               = TH1F('JetBTagger_inclusive_MC_TTbar_ele',           'JetBTagger_inclusive_MC_TTbar_ele',           11, 0, 1.1)
    JetBTagger_hist_1stjet_MC_TTbar_ele        = TH1F('JetBTagger_leadingJet_MC_TTbar_ele',          'JetBTagger_leadingJet_MC_TTbar_ele',          11, 0, 1.1)
    JetBTagger_hist_2ndjet_MC_TTbar_ele        = TH1F('JetBTagger_subLeadingJet_MC_TTbar_ele',       'JetBTagger_subLeadingJet_MC_TTbar_ele',       11, 0, 1.1)
    JetBTagger_hist_1stjet_fwdeta_MC_TTbar_ele = TH1F('JetBTagger_leadingJet_fwdeta_MC_TTbar_ele',   'JetBTagger_leadingJet_fwdeta_MC_TTbar_ele',   11, 0, 1.1)
    JetBTagger_hist_2ndjet_fwdeta_MC_TTbar_ele = TH1F('JetBTagger_subLeadingJet_fwdeta_MC_TTbar_ele','JetBTagger_subLeadingJet_fwdeta_MC_TTbar_ele',11, 0, 1.1)

    JetIsBtagged_hist_MC_TTbar_ele               = TH1F('JetIsBtagged_inclusive_MC_TTbar_ele',           'JetIsBtagged_inclusive_MC_TTbar_ele',           11, 0, 1.1)
    JetIsBtagged_hist_1stjet_MC_TTbar_ele        = TH1F('JetIsBtagged_leadingJet_MC_TTbar_ele',          'JetIsBtagged_leadingJet_MC_TTbar_ele',          11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_MC_TTbar_ele        = TH1F('JetIsBtagged_subLeadingJet_MC_TTbar_ele',       'JetIsBtagged_subLeadingJet_MC_TTbar_ele',       11, 0, 1.1)
    JetIsBtagged_hist_1stjet_fwdeta_MC_TTbar_ele = TH1F('JetIsBtagged_leadingJet_fwdeta_MC_TTbar_ele',   'JetIsBtagged_leadingJet_fwdeta_MC_TTbar_ele',   11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_fwdeta_MC_TTbar_ele = TH1F('JetIsBtagged_subLeadingJet_fwdeta_MC_TTbar_ele','JetIsBtagged_subLeadingJet_fwdeta_MC_TTbar_ele',11, 0, 1.1)
    
    JetIsBtaggedWithSF_hist_MC_TTbar_ele               = TH1F('JetIsBtaggedWithSF_inclusive_MC_TTbar_ele',           'JetIsBtaggedWithSF_inclusive_MC_TTbar_ele',           11, 0, 1.1)
    JetIsBtaggedWithSF_hist_1stjet_MC_TTbar_ele        = TH1F('JetIsBtaggedWithSF_leadingJet_MC_TTbar_ele',          'JetIsBtaggedWithSF_leadingJet_MC_TTbar_ele',          11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_MC_TTbar_ele        = TH1F('JetIsBtaggedWithSF_subLeadingJet_MC_TTbar_ele',       'JetIsBtaggedWithSF_subLeadingJet_MC_TTbar_ele',       11, 0, 1.1) 
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_TTbar_ele = TH1F('JetIsBtaggedWithSF_leadingJet_fwdeta_MC_TTbar_ele',   'JetIsBtaggedWithSF_leadingJet_fwdeta_MC_TTbar_ele',   11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_TTbar_ele = TH1F('JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_TTbar_ele','JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_TTbar_ele',11, 0, 1.1)

    PFMET_hist_MC_TTbar_ele = TH1F('PFMET_MC_TTbar_ele','PFMET_MC_TTbar_ele',125,0,250)


    # define data histograms Z->mumu

    # nCleanedJets_hist_MC_TTbar_mu                   = TH1F('nCleanedJets_MC_TTbar_mu',                   'nCleanedJets_MC_TTbar_mu',                   30, 0, 30) 
    nCleanedJetsPt30_hist_MC_TTbar_mu               = TH1F('nCleanedJetsPt30_MC_TTbar_mu',               'nCleanedJetsPt30_MC_TTbar_mu',               15, 0, 15) 
    # nCleanedJetsPt30BTagged_hist_MC_TTbar_mu        = TH1F('nCleanedJetsPt30BTagged_MC_TTbar_mu',        'nCleanedJetsPt30BTagged_MC_TTbar_mu',        5,  0, 5)
    # nCleanedJetsPt30BTagged_bTagSF_hist_MC_TTbar_mu = TH1F('nCleanedJetsPt30BTagged_bTagSF_MC_TTbar_mu', 'nCleanedJetsPt30BTagged_bTagSF_MC_TTbar_mu', 5,  0, 5)

    JetPt_hist_MC_TTbar_mu               = TH1F('JetPt_inclusive_MC_TTbar_mu',           'JetPt_inclusive_MC_TTbar_mu',           120, 0, 1200)
    JetPt_hist_1stJet_MC_TTbar_mu        = TH1F('JetPt_leadingJet_MC_TTbar_mu',          'JetPt_leadingJet_MC_TTbar_mu',          120, 0, 1200)
    JetPt_hist_2ndjet_MC_TTbar_mu        = TH1F('JetPt_subLeadingJet_MC_TTbar_mu',       'JetPt_subLeadingJet_MC_TTbar_mu',       120, 0, 1200)
    JetPt_hist_1stJet_fwdeta_MC_TTbar_mu = TH1F('JetPt_leadingJet_fwdeta_MC_TTbar_mu',   'JetPt_leadingJet_fwdeta_MC_TTbar_mu',   120, 0, 1200)
    JetPt_hist_2ndjet_fwdeta_MC_TTbar_mu = TH1F('JetPt_subLeadingJet_fwdeta_MC_TTbar_mu','JetPt_subLeadingJet_fwdeta_MC_TTbar_mu',120, 0, 1200)

    JetEta_hist_MC_TTbar_mu               = TH1F('JetEta_inclusive_MC_TTbar_mu',           'JetEta_inclusive_MC_TTbar_mu',           50, -5, 5)
    JetEta_hist_1stjet_MC_TTbar_mu        = TH1F('JetEta_leadingJet_MC_TTbar_mu',          'JetEta_leadingJet_MC_TTbar_mu',          50, -5, 5)
    JetEta_hist_2ndjet_MC_TTbar_mu        = TH1F('JetEta_subLeadingJet_MC_TTbar_mu',       'JetEta_subLeadingJet_MC_TTbar_mu',       50, -5, 5)
    JetEta_hist_1stjet_fwdeta_MC_TTbar_mu = TH1F('JetEta_leadingJet_fwdeta_MC_TTbar_mu',   'JetEta_leadingJet_fwdeta_MC_TTbar_mu',   50, -5, 5)
    JetEta_hist_2ndjet_fwdeta_MC_TTbar_mu = TH1F('JetEta_subLeadingJet_fwdeta_MC_TTbar_mu','JetEta_subLeadingJet_fwdeta_MC_TTbar_mu',50, -5, 5)

    JetPhi_hist_MC_TTbar_mu               = TH1F('JetPhi_inclusive_MC_TTbar_mu',           'JetPhi_inclusive_MC_TTbar_mu',           33, -3.3, 3.3)
    JetPhi_hist_1stjet_MC_TTbar_mu        = TH1F('JetPhi_leadingJet_MC_TTbar_mu',          'JetPhi_leadingJet_MC_TTbar_mu',          33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_MC_TTbar_mu        = TH1F('JetPhi_subLeadingJet_MC_TTbar_mu',       'JetPhi_subLeadingJet_MC_TTbar_mu',       33, -3.3, 3.3) 
    JetPhi_hist_1stjet_fwdeta_MC_TTbar_mu = TH1F('JetPhi_leadingJet_fwdeta_MC_TTbar_mu',   'JetPhi_leadingJet_fwdeta_MC_TTbar_mu',   33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_fwdeta_MC_TTbar_mu = TH1F('JetPhi_subLeadingJet_fwdeta_MC_TTbar_mu','JetPhi_subLeadingJet_fwdeta_MC_TTbar_mu',33, -3.3, 3.3) 

    JetBTagger_hist_MC_TTbar_mu               = TH1F('JetBTagger_inclusive_MC_TTbar_mu',           'JetBTagger_inclusive_MC_TTbar_mu',           11, 0, 1.1)
    JetBTagger_hist_1stjet_MC_TTbar_mu        = TH1F('JetBTagger_leadingJet_MC_TTbar_mu',          'JetBTagger_leadingJet_MC_TTbar_mu',          11, 0, 1.1)
    JetBTagger_hist_2ndjet_MC_TTbar_mu        = TH1F('JetBTagger_subLeadingJet_MC_TTbar_mu',       'JetBTagger_subLeadingJet_MC_TTbar_mu',       11, 0, 1.1)
    JetBTagger_hist_1stjet_fwdeta_MC_TTbar_mu = TH1F('JetBTagger_leadingJet_fwdeta_MC_TTbar_mu',   'JetBTagger_leadingJet_fwdeta_MC_TTbar_mu',   11, 0, 1.1)
    JetBTagger_hist_2ndjet_fwdeta_MC_TTbar_mu = TH1F('JetBTagger_subLeadingJet_fwdeta_MC_TTbar_mu','JetBTagger_subLeadingJet_fwdeta_MC_TTbar_mu',11, 0, 1.1)

    JetIsBtagged_hist_MC_TTbar_mu               = TH1F('JetIsBtagged_inclusive_MC_TTbar_mu',           'JetIsBtagged_inclusive_MC_TTbar_mu',           11, 0, 1.1)
    JetIsBtagged_hist_1stjet_MC_TTbar_mu        = TH1F('JetIsBtagged_leadingJet_MC_TTbar_mu',          'JetIsBtagged_leadingJet_MC_TTbar_mu',          11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_MC_TTbar_mu        = TH1F('JetIsBtagged_subLeadingJet_MC_TTbar_mu',       'JetIsBtagged_subLeadingJet_MC_TTbar_mu',       11, 0, 1.1)
    JetIsBtagged_hist_1stjet_fwdeta_MC_TTbar_mu = TH1F('JetIsBtagged_leadingJet_fwdeta_MC_TTbar_mu',   'JetIsBtagged_leadingJet_fwdeta_MC_TTbar_mu',   11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_fwdeta_MC_TTbar_mu = TH1F('JetIsBtagged_subLeadingJet_fwdeta_MC_TTbar_mu','JetIsBtagged_subLeadingJet_fwdeta_MC_TTbar_mu',11, 0, 1.1)
    
    JetIsBtaggedWithSF_hist_MC_TTbar_mu               = TH1F('JetIsBtaggedWithSF_inclusive_MC_TTbar_mu',           'JetIsBtaggedWithSF_inclusive_MC_TTbar_mu',           11, 0, 1.1)
    JetIsBtaggedWithSF_hist_1stjet_MC_TTbar_mu        = TH1F('JetIsBtaggedWithSF_leadingJet_MC_TTbar_mu',          'JetIsBtaggedWithSF_leadingJet_MC_TTbar_mu',          11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_MC_TTbar_mu        = TH1F('JetIsBtaggedWithSF_subLeadingJet_MC_TTbar_mu',       'JetIsBtaggedWithSF_subLeadingJet_MC_TTbar_mu',       11, 0, 1.1) 
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_TTbar_mu = TH1F('JetIsBtaggedWithSF_leadingJet_fwdeta_MC_TTbar_mu',   'JetIsBtaggedWithSF_leadingJet_fwdeta_MC_TTbar_mu',   11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_TTbar_mu = TH1F('JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_TTbar_mu','JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_TTbar_mu',11, 0, 1.1)

    PFMET_hist_MC_TTbar_mu = TH1F('PFMET_MC_TTbar_mu','PFMET_MC_TTbar_mu',125,0,250)




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

            # nCleanedJets_hist_MC_TTbar_ele.Fill(event.nCleanedJets,weight)
            nCleanedJetsPt30_hist_MC_TTbar_ele.Fill(event.nCleanedJetsPt30,weight)
            # nCleanedJetsPt30BTagged_hist_MC_TTbar_ele.Fill(event.nCleanedJetsPt30BTagged,weight)
            # nCleanedJetsPt30BTagged_bTagSF_hist_MC_TTbar_ele.Fill(event.nCleanedJetsPt30BTagged_bTagSF,weight)
            
            for i in range(len(event.JetPt)) :
                JetPt_hist_MC_TTbar_ele.Fill(event.JetPt[i],weight)
                      
            for i in range(len(event.JetEta)) :  
                JetEta_hist_MC_TTbar_ele.Fill(event.JetEta[i],weight)
                                            
            for i in range(len(event.JetPhi)) :
                JetPhi_hist_MC_TTbar_ele.Fill(event.JetPhi[i],weight)
              

            if event.nCleanedJetsPt30 > 0 :

                for i in range(len(event.JetBTagger)) :
                    JetBTagger_hist_MC_TTbar_ele.Fill(event.JetBTagger[i],weight)
                                            
                for i in range(len(event.JetIsBtagged)) :
                    JetIsBtagged_hist_MC_TTbar_ele.Fill(event.JetIsBtagged[i],weight)
            
                for i in range(len(event.JetIsBtaggedWithSF)) :                
                    JetIsBtaggedWithSF_hist_MC_TTbar_ele.Fill(event.JetIsBtaggedWithSF[i],weight)     
            
            
                JetPt_hist_1stJet_MC_TTbar_ele.Fill(event.JetPt[0],weight)
                JetEta_hist_1stjet_MC_TTbar_ele.Fill(event.JetEta[0],weight)
                JetPhi_hist_1stjet_MC_TTbar_ele.Fill(event.JetPhi[0],weight)
                JetBTagger_hist_1stjet_MC_TTbar_ele.Fill(event.JetBTagger[0],weight)
                JetIsBtagged_hist_1stjet_MC_TTbar_ele.Fill(event.JetIsBtagged[0],weight)
                JetIsBtaggedWithSF_hist_1stjet_MC_TTbar_ele.Fill(event.JetIsBtaggedWithSF[0],weight)
            
                if math.fabs(event.JetEta[0]) > 3.0 : 
                    JetPt_hist_1stJet_fwdeta_MC_TTbar_ele.Fill(event.JetPt[0],weight)
                    JetEta_hist_1stjet_fwdeta_MC_TTbar_ele.Fill(event.JetEta[0],weight)
                    JetPhi_hist_1stjet_fwdeta_MC_TTbar_ele.Fill(event.JetPhi[0],weight)
                    JetBTagger_hist_1stjet_fwdeta_MC_TTbar_ele.Fill(event.JetBTagger[0],weight)
                    JetIsBtagged_hist_1stjet_fwdeta_MC_TTbar_ele.Fill(event.JetIsBtagged[0],weight)
                    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_TTbar_ele.Fill(event.JetIsBtaggedWithSF[0],weight)
            
            
                if event.nCleanedJetsPt30 > 1 :
                    JetPt_hist_2ndjet_MC_TTbar_ele.Fill(event.JetPt[1],weight)
                    JetEta_hist_2ndjet_MC_TTbar_ele.Fill(event.JetEta[1],weight)
                    JetPhi_hist_2ndjet_MC_TTbar_ele.Fill(event.JetPhi[1],weight)
                    JetBTagger_hist_2ndjet_MC_TTbar_ele.Fill(event.JetBTagger[1],weight)
                    JetIsBtagged_hist_2ndjet_MC_TTbar_ele.Fill(event.JetIsBtagged[1],weight)
                    JetIsBtaggedWithSF_hist_2ndjet_MC_TTbar_ele.Fill(event.JetIsBtaggedWithSF[1],weight)
            
                    if math.fabs(event.JetEta[1]) > 3.0 : 
                        JetPt_hist_2ndjet_fwdeta_MC_TTbar_ele.Fill(event.JetPt[1],weight)
                        JetEta_hist_2ndjet_fwdeta_MC_TTbar_ele.Fill(event.JetEta[1],weight)
                        JetPhi_hist_2ndjet_fwdeta_MC_TTbar_ele.Fill(event.JetPhi[1],weight)
                        JetBTagger_hist_2ndjet_fwdeta_MC_TTbar_ele.Fill(event.JetBTagger[1],weight)
                        JetIsBtagged_hist_2ndjet_fwdeta_MC_TTbar_ele.Fill(event.JetIsBtagged[1],weight)
                        JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_TTbar_ele.Fill(event.JetIsBtaggedWithSF[1],weight)
                        
            
            PFMET_hist_MC_TTbar_ele.Fill(event.PFMET,weight)



        # Z->mumu histos
        if(int(math.fabs(event.LepLepId[0])) == 13 ):

            # nCleanedJets_hist_MC_TTbar_mu.Fill(event.nCleanedJets,weight)
            nCleanedJetsPt30_hist_MC_TTbar_mu.Fill(event.nCleanedJetsPt30,weight)
            # nCleanedJetsPt30BTagged_hist_MC_TTbar_mu.Fill(event.nCleanedJetsPt30BTagged,weight)
            # nCleanedJetsPt30BTagged_bTagSF_hist_MC_TTbar_mu.Fill(event.nCleanedJetsPt30BTagged_bTagSF,weight)
            
            for i in range(len(event.JetPt)) :
                JetPt_hist_MC_TTbar_mu.Fill(event.JetPt[i],weight)
                      
            for i in range(len(event.JetEta)) :  
                JetEta_hist_MC_TTbar_mu.Fill(event.JetEta[i],weight)
                                            
            for i in range(len(event.JetPhi)) :
                JetPhi_hist_MC_TTbar_mu.Fill(event.JetPhi[i],weight)

              
            if event.nCleanedJetsPt30 > 0 :

                for i in range(len(event.JetBTagger)) :
                    JetBTagger_hist_MC_TTbar_mu.Fill(event.JetBTagger[i],weight)
                                            
                for i in range(len(event.JetIsBtagged)) :
                    JetIsBtagged_hist_MC_TTbar_mu.Fill(event.JetIsBtagged[i],weight)
            
                for i in range(len(event.JetIsBtaggedWithSF)) :                
                    JetIsBtaggedWithSF_hist_MC_TTbar_mu.Fill(event.JetIsBtaggedWithSF[i],weight)     
            
            
                JetPt_hist_1stJet_MC_TTbar_mu.Fill(event.JetPt[0],weight)
                JetEta_hist_1stjet_MC_TTbar_mu.Fill(event.JetEta[0],weight)
                JetPhi_hist_1stjet_MC_TTbar_mu.Fill(event.JetPhi[0],weight)
                JetBTagger_hist_1stjet_MC_TTbar_mu.Fill(event.JetBTagger[0],weight)
                JetIsBtagged_hist_1stjet_MC_TTbar_mu.Fill(event.JetIsBtagged[0],weight)
                JetIsBtaggedWithSF_hist_1stjet_MC_TTbar_mu.Fill(event.JetIsBtaggedWithSF[0],weight)
            
                if math.fabs(event.JetEta[0]) > 3.0 : 
                    JetPt_hist_1stJet_fwdeta_MC_TTbar_mu.Fill(event.JetPt[0],weight)
                    JetEta_hist_1stjet_fwdeta_MC_TTbar_mu.Fill(event.JetEta[0],weight)
                    JetPhi_hist_1stjet_fwdeta_MC_TTbar_mu.Fill(event.JetPhi[0],weight)
                    JetBTagger_hist_1stjet_fwdeta_MC_TTbar_mu.Fill(event.JetBTagger[0],weight)
                    JetIsBtagged_hist_1stjet_fwdeta_MC_TTbar_mu.Fill(event.JetIsBtagged[0],weight)
                    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_TTbar_mu.Fill(event.JetIsBtaggedWithSF[0],weight)
            
            
                if event.nCleanedJetsPt30 > 1 :
                    JetPt_hist_2ndjet_MC_TTbar_mu.Fill(event.JetPt[1],weight)
                    JetEta_hist_2ndjet_MC_TTbar_mu.Fill(event.JetEta[1],weight)
                    JetPhi_hist_2ndjet_MC_TTbar_mu.Fill(event.JetPhi[1],weight)
                    JetBTagger_hist_2ndjet_MC_TTbar_mu.Fill(event.JetBTagger[1],weight)
                    JetIsBtagged_hist_2ndjet_MC_TTbar_mu.Fill(event.JetIsBtagged[1],weight)
                    JetIsBtaggedWithSF_hist_2ndjet_MC_TTbar_mu.Fill(event.JetIsBtaggedWithSF[1],weight)
            
                    if math.fabs(event.JetEta[1]) > 3.0 : 
                        JetPt_hist_2ndjet_fwdeta_MC_TTbar_mu.Fill(event.JetPt[1],weight)
                        JetEta_hist_2ndjet_fwdeta_MC_TTbar_mu.Fill(event.JetEta[1],weight)
                        JetPhi_hist_2ndjet_fwdeta_MC_TTbar_mu.Fill(event.JetPhi[1],weight)
                        JetBTagger_hist_2ndjet_fwdeta_MC_TTbar_mu.Fill(event.JetBTagger[1],weight)
                        JetIsBtagged_hist_2ndjet_fwdeta_MC_TTbar_mu.Fill(event.JetIsBtagged[1],weight)
                        JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_TTbar_mu.Fill(event.JetIsBtaggedWithSF[1],weight)
                        
            
            PFMET_hist_MC_TTbar_mu.Fill(event.PFMET,weight)


    

        
    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_MCTTbar = TFile.Open("jetsMetDistrib_MC_TTbar_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_MCTTbar.cd()

    # Zee histos
    # nCleanedJets_hist_MC_TTbar_ele.Write()
    nCleanedJetsPt30_hist_MC_TTbar_ele.Write()
    # nCleanedJetsPt30BTagged_hist_MC_TTbar_ele.Write()
    # nCleanedJetsPt30BTagged_bTagSF_hist_MC_TTbar_ele.Write()
    
    JetPt_hist_MC_TTbar_ele.Write()
    JetPt_hist_1stJet_MC_TTbar_ele.Write()
    JetPt_hist_2ndjet_MC_TTbar_ele.Write()
    JetPt_hist_1stJet_fwdeta_MC_TTbar_ele.Write()
    JetPt_hist_2ndjet_fwdeta_MC_TTbar_ele.Write()
                            
    JetEta_hist_MC_TTbar_ele.Write()
    JetEta_hist_1stjet_MC_TTbar_ele.Write()
    JetEta_hist_2ndjet_MC_TTbar_ele.Write()
    JetEta_hist_1stjet_fwdeta_MC_TTbar_ele.Write()
    JetEta_hist_2ndjet_fwdeta_MC_TTbar_ele.Write()
                            
    JetPhi_hist_MC_TTbar_ele.Write()
    JetPhi_hist_1stjet_MC_TTbar_ele.Write()
    JetPhi_hist_2ndjet_MC_TTbar_ele.Write()
    JetPhi_hist_1stjet_fwdeta_MC_TTbar_ele.Write()
    JetPhi_hist_2ndjet_fwdeta_MC_TTbar_ele.Write()
                            
    JetBTagger_hist_MC_TTbar_ele.Write()
    JetBTagger_hist_1stjet_MC_TTbar_ele.Write()
    JetBTagger_hist_2ndjet_MC_TTbar_ele.Write()
    JetBTagger_hist_1stjet_fwdeta_MC_TTbar_ele.Write()
    JetBTagger_hist_2ndjet_fwdeta_MC_TTbar_ele.Write()
                            
    JetIsBtagged_hist_MC_TTbar_ele.Write()
    JetIsBtagged_hist_1stjet_MC_TTbar_ele.Write()
    JetIsBtagged_hist_2ndjet_MC_TTbar_ele.Write()
    JetIsBtagged_hist_1stjet_fwdeta_MC_TTbar_ele.Write()
    JetIsBtagged_hist_2ndjet_fwdeta_MC_TTbar_ele.Write()
                            
    JetIsBtaggedWithSF_hist_MC_TTbar_ele.Write()   
    JetIsBtaggedWithSF_hist_1stjet_MC_TTbar_ele.Write()
    JetIsBtaggedWithSF_hist_2ndjet_MC_TTbar_ele.Write()
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_TTbar_ele.Write()
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_TTbar_ele.Write()

    PFMET_hist_MC_TTbar_ele.Write()

    
    # Zmumu histos
    # nCleanedJets_hist_MC_TTbar_mu.Write()
    nCleanedJetsPt30_hist_MC_TTbar_mu.Write()
    # nCleanedJetsPt30BTagged_hist_MC_TTbar_mu.Write()
    # nCleanedJetsPt30BTagged_bTagSF_hist_MC_TTbar_mu.Write()
    
    JetPt_hist_MC_TTbar_mu.Write()
    JetPt_hist_1stJet_MC_TTbar_mu.Write()
    JetPt_hist_2ndjet_MC_TTbar_mu.Write()
    JetPt_hist_1stJet_fwdeta_MC_TTbar_mu.Write()
    JetPt_hist_2ndjet_fwdeta_MC_TTbar_mu.Write()
                            
    JetEta_hist_MC_TTbar_mu.Write()
    JetEta_hist_1stjet_MC_TTbar_mu.Write()
    JetEta_hist_2ndjet_MC_TTbar_mu.Write()
    JetEta_hist_1stjet_fwdeta_MC_TTbar_mu.Write()
    JetEta_hist_2ndjet_fwdeta_MC_TTbar_mu.Write()
                            
    JetPhi_hist_MC_TTbar_mu.Write()
    JetPhi_hist_1stjet_MC_TTbar_mu.Write()
    JetPhi_hist_2ndjet_MC_TTbar_mu.Write()
    JetPhi_hist_1stjet_fwdeta_MC_TTbar_mu.Write()
    JetPhi_hist_2ndjet_fwdeta_MC_TTbar_mu.Write()
                            
    JetBTagger_hist_MC_TTbar_mu.Write()
    JetBTagger_hist_1stjet_MC_TTbar_mu.Write()
    JetBTagger_hist_2ndjet_MC_TTbar_mu.Write()
    JetBTagger_hist_1stjet_fwdeta_MC_TTbar_mu.Write()
    JetBTagger_hist_2ndjet_fwdeta_MC_TTbar_mu.Write()
                            
    JetIsBtagged_hist_MC_TTbar_mu.Write()
    JetIsBtagged_hist_1stjet_MC_TTbar_mu.Write()
    JetIsBtagged_hist_2ndjet_MC_TTbar_mu.Write()
    JetIsBtagged_hist_1stjet_fwdeta_MC_TTbar_mu.Write()
    JetIsBtagged_hist_2ndjet_fwdeta_MC_TTbar_mu.Write()
                            
    JetIsBtaggedWithSF_hist_MC_TTbar_mu.Write()   
    JetIsBtaggedWithSF_hist_1stjet_MC_TTbar_mu.Write()
    JetIsBtaggedWithSF_hist_2ndjet_MC_TTbar_mu.Write()
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_TTbar_mu.Write()
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_TTbar_mu.Write()

    PFMET_hist_MC_TTbar_mu.Write()


    outFile_MCTTbar.Close()
    print "MC TTbar histo file created!"
# ********************



# ****************************************
# create output directory 
outputDir = "jetsMetDistrib_DATAvsMC_" + str(period) + "_" + str(treeText)
gSystem.Exec("mkdir -p " + outputDir)
print "Output directory created!"


# **************************
# read data histos from file 
histoDATA_input = TFile.Open("jetsMetDistrib_DATA_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoDATA_input.GetName(),'...'

inDATA_list = []

# inDATA_list.append(histoDATA_input.Get('nCleanedJets_ele'))
inDATA_list.append(histoDATA_input.Get('nCleanedJetsPt30_ele'))
# inDATA_list.append(histoDATA_input.Get('nCleanedJetsPt30BTagged_ele'))
# inDATA_list.append(histoDATA_input.Get('nCleanedJetsPt30BTagged_bTagSF_ele'))
inDATA_list.append(histoDATA_input.Get('JetPt_inclusive_ele'))
inDATA_list.append(histoDATA_input.Get('JetPt_leadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetPt_subLeadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetPt_leadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetPt_subLeadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetEta_inclusive_ele'))
inDATA_list.append(histoDATA_input.Get('JetEta_leadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetEta_subLeadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetEta_leadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetEta_subLeadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetPhi_inclusive_ele'))
inDATA_list.append(histoDATA_input.Get('JetPhi_leadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetPhi_subLeadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetPhi_leadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetPhi_subLeadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_inclusive_ele'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_leadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_subLeadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_leadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_subLeadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_inclusive_ele'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_leadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_subLeadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_leadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_subLeadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_inclusive_ele'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_leadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_subLeadingJet_ele'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_leadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_subLeadingJet_fwdeta_ele'))
inDATA_list.append(histoDATA_input.Get('PFMET_ele'))

# inDATA_list.append(histoDATA_input.Get('nCleanedJets_mu'))
inDATA_list.append(histoDATA_input.Get('nCleanedJetsPt30_mu'))
# inDATA_list.append(histoDATA_input.Get('nCleanedJetsPt30BTagged_mu'))
# inDATA_list.append(histoDATA_input.Get('nCleanedJetsPt30BTagged_bTagSF_mu'))
inDATA_list.append(histoDATA_input.Get('JetPt_inclusive_mu'))
inDATA_list.append(histoDATA_input.Get('JetPt_leadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetPt_subLeadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetPt_leadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetPt_subLeadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetEta_inclusive_mu'))
inDATA_list.append(histoDATA_input.Get('JetEta_leadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetEta_subLeadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetEta_leadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetEta_subLeadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetPhi_inclusive_mu'))
inDATA_list.append(histoDATA_input.Get('JetPhi_leadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetPhi_subLeadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetPhi_leadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetPhi_subLeadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_inclusive_mu'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_leadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_subLeadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_leadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_subLeadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_inclusive_mu'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_leadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_subLeadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_leadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_subLeadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_inclusive_mu'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_leadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_subLeadingJet_mu'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_leadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_subLeadingJet_fwdeta_mu'))
inDATA_list.append(histoDATA_input.Get('PFMET_mu'))


# ****************************
# read DY MC histos from file 
histoMCDY_input = TFile.Open("jetsMetDistrib_MC_DY_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoMCDY_input.GetName(),'...'

inMCDY_list = []

# inMCDY_list.append(histoMCDY_input.Get('nCleanedJets_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('nCleanedJetsPt30_MC_DY_ele'))
# inMCDY_list.append(histoMCDY_input.Get('nCleanedJetsPt30BTagged_MC_DY_ele'))
# inMCDY_list.append(histoMCDY_input.Get('nCleanedJetsPt30BTagged_bTagSF_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_inclusive_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_leadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_subLeadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_leadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_subLeadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_inclusive_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_leadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_subLeadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_leadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_subLeadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_inclusive_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_leadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_subLeadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_leadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_subLeadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_inclusive_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_leadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_subLeadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_leadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_subLeadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_inclusive_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_leadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_subLeadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_leadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_subLeadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_inclusive_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_leadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_subLeadingJet_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_leadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_DY_ele'))
inMCDY_list.append(histoMCDY_input.Get('PFMET_MC_DY_ele'))

# inMCDY_list.append(histoMCDY_input.Get('nCleanedJets_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('nCleanedJetsPt30_MC_DY_mu'))
# inMCDY_list.append(histoMCDY_input.Get('nCleanedJetsPt30BTagged_MC_DY_mu'))
# inMCDY_list.append(histoMCDY_input.Get('nCleanedJetsPt30BTagged_bTagSF_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_inclusive_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_leadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_subLeadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_leadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_subLeadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_inclusive_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_leadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_subLeadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_leadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_subLeadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_inclusive_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_leadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_subLeadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_leadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_subLeadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_inclusive_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_leadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_subLeadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_leadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_subLeadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_inclusive_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_leadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_subLeadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_leadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_subLeadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_inclusive_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_leadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_subLeadingJet_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_leadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_DY_mu'))
inMCDY_list.append(histoMCDY_input.Get('PFMET_MC_DY_mu'))



# ****************************
# read TTbar MC histos from file 
histoMCTTbar_input = TFile.Open("jetsMetDistrib_MC_TTbar_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoMCTTbar_input.GetName(),'...'

inMCTTbar_list = []

# inMCTTbar_list.append(histoMCTTbar_input.Get('nCleanedJets_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('nCleanedJetsPt30_MC_TTbar_ele'))
# inMCTTbar_list.append(histoMCTTbar_input.Get('nCleanedJetsPt30BTagged_MC_TTbar_ele'))
# inMCTTbar_list.append(histoMCTTbar_input.Get('nCleanedJetsPt30BTagged_bTagSF_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPt_inclusive_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPt_leadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPt_subLeadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPt_leadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPt_subLeadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetEta_inclusive_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetEta_leadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetEta_subLeadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetEta_leadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetEta_subLeadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPhi_inclusive_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPhi_leadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPhi_subLeadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPhi_leadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPhi_subLeadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetBTagger_inclusive_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetBTagger_leadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetBTagger_subLeadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetBTagger_leadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetBTagger_subLeadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtagged_inclusive_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtagged_leadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtagged_subLeadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtagged_leadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtagged_subLeadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtaggedWithSF_inclusive_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtaggedWithSF_leadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtaggedWithSF_subLeadingJet_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtaggedWithSF_leadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_TTbar_ele'))
inMCTTbar_list.append(histoMCTTbar_input.Get('PFMET_MC_TTbar_ele'))

# inMCTTbar_list.append(histoMCTTbar_input.Get('nCleanedJets_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('nCleanedJetsPt30_MC_TTbar_mu'))
# inMCTTbar_list.append(histoMCTTbar_input.Get('nCleanedJetsPt30BTagged_MC_TTbar_mu'))
# inMCTTbar_list.append(histoMCTTbar_input.Get('nCleanedJetsPt30BTagged_bTagSF_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPt_inclusive_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPt_leadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPt_subLeadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPt_leadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPt_subLeadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetEta_inclusive_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetEta_leadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetEta_subLeadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetEta_leadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetEta_subLeadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPhi_inclusive_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPhi_leadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPhi_subLeadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPhi_leadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetPhi_subLeadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetBTagger_inclusive_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetBTagger_leadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetBTagger_subLeadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetBTagger_leadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetBTagger_subLeadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtagged_inclusive_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtagged_leadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtagged_subLeadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtagged_leadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtagged_subLeadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtaggedWithSF_inclusive_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtaggedWithSF_leadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtaggedWithSF_subLeadingJet_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtaggedWithSF_leadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_TTbar_mu'))
inMCTTbar_list.append(histoMCTTbar_input.Get('PFMET_MC_TTbar_mu'))

       

# ******************************
# do DATA vs MC comparison plots  
for i in range(len(inDATA_list)) : 

    canvas = TCanvas("canvas","canvas",800,800)

    hs = THStack("hs","")

    #DATA hist
    inDATA_list[i].SetMarkerStyle(20)
    inDATA_list[i].SetMarkerSize(0.6)


    #MC TTbar hist
    inMCTTbar_list[i].SetFillColor(kAzure-2)
    inMCTTbar_list[i].SetLineColor(kBlack)
    hs.Add(inMCTTbar_list[i])

    #MC DY hist
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
