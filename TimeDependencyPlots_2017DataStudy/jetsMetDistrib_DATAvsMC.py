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
redoDATAHistos = True
redoMCDYHistos = True

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
    inputDATAtree = TFile.Open("/data3/Higgs/170222/AllData/ZZ4lAnalysis.root") #2016 data
    inputMCDYtree = TFile.Open("/data3/Higgs/170222/DYJetsToLL_M50/ZZ4lAnalysis.root") #DYJets 2016 MC
    lumi     = 35.9  # fb-1
    lumiText = '35.9 fb^{-1}'
    if(ZZTree):
        treeDATA  = inputDATAtree.Get("ZZTree/candTree")
        treeMCDY  = inputMCDYtree.Get("ZZTree/candTree")
        treeText  = "ZZTree"
    elif(CRZLLTree):
        treeDATA  = inputDATAtree.Get("CRZLLTree/candTree")
        treeMCDY  = inputMCDYTtee.Get("CRZLLTree/candTree")
        treeText  = "CRZLLTree"
    elif(CRZLTree):
        treeDATA  = inputDATAtree.Get("CRZLTree/candTree")
        treeMCDY  = inputMCDYtree.Get("CRZLTree/candTree")
        treeText  = "CRZLTree"
    elif(ZTree):
        treeDATA  = inputDATAtree.Get("ZTree/candTree")
        treeMCDY  = inputMCDYtree.Get("ZTree/candTree")
        treeText  = "ZTree"
    else:
        print ("Error: wrong option!")

elif(period == "data2017"):
    inputDATAtree = TFile.Open("/data3/Higgs/180122/AllData/ZZ4lAnalysis.root") #2017 data  
    inputMCDYtree = TFile.Open("/data3/Higgs/180122/DYJetsToLL_M50/ZZ4lAnalysis.root") #DYJets 2017 MC
    lumi     = 41.37   # fb-1
    lumiText = '41.37 fb^{-1}'
    if(ZZTree):
        treeDATA  = inputDATAtree.Get("ZZTree/candTree")
        treeMCDY  = inputMCDYtree.Get("ZZTree/candTree")
        treeText  = "ZZTree"
    elif(CRZLLTree):
        treeDATA  = inputDATAtree.Get("CRZLLTree/candTree")
        treeMCDY  = inputMCDYtree.Get("CRZLLTree/candTree")
        treeText  = "CRZLLTree"
    elif(CRZLTree):
        treeDATA  = inputDATAtree.Get("CRZLTree/candTree")
        treeMCDY  = inputMCDYtree.Get("CRZLTree/candTree")
        treeText  = "CRZLTree"
    elif(ZTree):
        treeDATA  = inputDATAtree.Get("ZTree/candTree")
        treeMCDY  = inputMCDYtree.Get("ZTree/candTree")
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

    # define data histograms 
    nCleanedJets_hist                   = TH1F('nCleanedJets',                   'nCleanedJets',                   30, 0, 30) 
    nCleanedJetsPt30_hist               = TH1F('nCleanedJetsPt30',               'nCleanedJetsPt30',               15, 0, 15) 
    nCleanedJetsPt30BTagged_hist        = TH1F('nCleanedJetsPt30BTagged',        'nCleanedJetsPt30BTagged',        5,  0, 5)
    nCleanedJetsPt30BTagged_bTagSF_hist = TH1F('nCleanedJetsPt30BTagged_bTagSF', 'nCleanedJetsPt30BTagged_bTagSF', 5,  0, 5)

    JetPt_hist               = TH1F('JetPt_inclusive',           'JetPt_inclusive',           120, 0, 1200)
    JetPt_hist_1stJet        = TH1F('JetPt_leadingJet',          'JetPt_leadingJet',          120, 0, 1200)
    JetPt_hist_2ndjet        = TH1F('JetPt_subLeadingJet',       'JetPt_subLeadingJet',       120, 0, 1200)
    JetPt_hist_1stJet_fwdeta = TH1F('JetPt_leadingJet_fwdeta',   'JetPt_leadingJet_fwdeta',   120, 0, 1200)
    JetPt_hist_2ndjet_fwdeta = TH1F('JetPt_subLeadingJet_fwdeta','JetPt_subLeadingJet_fwdeta',120, 0, 1200)

    JetEta_hist               = TH1F('JetEta_inclusive',           'JetEta_inclusive',           50, -5, 5)
    JetEta_hist_1stjet        = TH1F('JetEta_leadingJet',          'JetEta_leadingJet',          50, -5, 5)
    JetEta_hist_2ndjet        = TH1F('JetEta_subLeadingJet',       'JetEta_subLeadingJet',       50, -5, 5)
    JetEta_hist_1stjet_fwdeta = TH1F('JetEta_leadingJet_fwdeta',   'JetEta_leadingJet_fwdeta',   50, -5, 5)
    JetEta_hist_2ndjet_fwdeta = TH1F('JetEta_subLeadingJet_fwdeta','JetEta_subLeadingJet_fwdeta',50, -5, 5)

    JetPhi_hist               = TH1F('JetPhi_inclusive',           'JetPhi_inclusive',           33, -3.3, 3.3)
    JetPhi_hist_1stjet        = TH1F('JetPhi_leadingJet',          'JetPhi_leadingJet',          33, -3.3, 3.3) 
    JetPhi_hist_2ndjet        = TH1F('JetPhi_subLeadingJet',       'JetPhi_subLeadingJet',       33, -3.3, 3.3) 
    JetPhi_hist_1stjet_fwdeta = TH1F('JetPhi_leadingJet_fwdeta',   'JetPhi_leadingJet_fwdeta',   33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_fwdeta = TH1F('JetPhi_subLeadingJet_fwdeta','JetPhi_subLeadingJet_fwdeta',33, -3.3, 3.3) 

    JetBTagger_hist               = TH1F('JetBTagger_inclusive',           'JetBTagger_inclusive',           20, -11, 2)
    JetBTagger_hist_1stjet        = TH1F('JetBTagger_leadingJet',          'JetBTagger_leadingJet',          20, -11, 2)
    JetBTagger_hist_2ndjet        = TH1F('JetBTagger_subLeadingJet',       'JetBTagger_subLeadingJet',       20, -11, 2)
    JetBTagger_hist_1stjet_fwdeta = TH1F('JetBTagger_leadingJet_fwdeta',   'JetBTagger_leadingJet_fwdeta',   20, -11, 2)
    JetBTagger_hist_2ndjet_fwdeta = TH1F('JetBTagger_subLeadingJet_fwdeta','JetBTagger_subLeadingJet_fwdeta',20, -11, 2)

    JetIsBtagged_hist               = TH1F('JetIsBtagged_inclusive',           'JetIsBtagged_inclusive',           11, 0, 1.1)
    JetIsBtagged_hist_1stjet        = TH1F('JetIsBtagged_leadingJet',          'JetIsBtagged_leadingJet',          11, 0, 1.1)
    JetIsBtagged_hist_2ndjet        = TH1F('JetIsBtagged_subLeadingJet',       'JetIsBtagged_subLeadingJet',       11, 0, 1.1)
    JetIsBtagged_hist_1stjet_fwdeta = TH1F('JetIsBtagged_leadingJet_fwdeta',   'JetIsBtagged_leadingJet_fwdeta',   11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_fwdeta = TH1F('JetIsBtagged_subLeadingJet_fwdeta','JetIsBtagged_subLeadingJet_fwdeta',11, 0, 1.1)
    
    JetIsBtaggedWithSF_hist               = TH1F('JetIsBtaggedWithSF_inclusive',           'JetIsBtaggedWithSF_inclusive',           11, 0, 1.1)
    JetIsBtaggedWithSF_hist_1stjet        = TH1F('JetIsBtaggedWithSF_leadingJet',          'JetIsBtaggedWithSF_leadingJet',          11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet        = TH1F('JetIsBtaggedWithSF_subLeadingJet',       'JetIsBtaggedWithSF_subLeadingJet',       11, 0, 1.1) 
    JetIsBtaggedWithSF_hist_1stjet_fwdeta = TH1F('JetIsBtaggedWithSF_leadingJet_fwdeta',   'JetIsBtaggedWithSF_leadingJet_fwdeta',   11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta = TH1F('JetIsBtaggedWithSF_subLeadingJet_fwdeta','JetIsBtaggedWithSF_subLeadingJet_fwdeta',11, 0, 1.1)

    PFMET_hist = TH1F('PFMET','PFMET',250,0,250)

    

    # read tree 
    print "reading tree", inputDATAtree.GetName(),treeText,treeDATA.GetName()  ,"..."
    for event in treeDATA:
        if ZTree :
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger
        else :
            if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger
  

        nCleanedJets_hist.Fill(event.nCleanedJets)
        nCleanedJetsPt30_hist.Fill(event.nCleanedJetsPt30)
        nCleanedJetsPt30BTagged_hist.Fill(event.nCleanedJetsPt30BTagged)
        nCleanedJetsPt30BTagged_bTagSF_hist.Fill(event.nCleanedJetsPt30BTagged_bTagSF)
        
        for i in range(len(event.JetPt)) :
            JetPt_hist.Fill(event.JetPt[i])
                  
        for i in range(len(event.JetEta)) :  
            JetEta_hist.Fill(event.JetEta[i])
                                        
        for i in range(len(event.JetPhi)) :
            JetPhi_hist.Fill(event.JetPhi[i])
          
        for i in range(len(event.JetBTagger)) :
            JetBTagger_hist.Fill(event.JetBTagger[i])
                                        
        for i in range(len(event.JetIsBtagged)) :
            JetIsBtagged_hist.Fill(event.JetIsBtagged[i])
        
        for i in range(len(event.JetIsBtaggedWithSF)) :                
            JetIsBtaggedWithSF_hist.Fill(event.JetIsBtaggedWithSF[i])     
        

        if event.nCleanedJets > 0 :
            JetPt_hist_1stJet.Fill(event.JetPt[0])
            JetEta_hist_1stjet.Fill(event.JetEta[0])
            JetPhi_hist_1stjet.Fill(event.JetPhi[0])
            JetBTagger_hist_1stjet.Fill(event.JetBTagger[0])
            JetIsBtagged_hist_1stjet.Fill(event.JetIsBtagged[0])
            JetIsBtaggedWithSF_hist_1stjet.Fill(event.JetIsBtaggedWithSF[0])

            if math.fabs(event.JetEta[0]) > 3.0 : 
                JetPt_hist_1stJet_fwdeta.Fill(event.JetPt[0])
                JetEta_hist_1stjet_fwdeta.Fill(event.JetEta[0])
                JetPhi_hist_1stjet_fwdeta.Fill(event.JetPhi[0])
                JetBTagger_hist_1stjet_fwdeta.Fill(event.JetBTagger[0])
                JetIsBtagged_hist_1stjet_fwdeta.Fill(event.JetIsBtagged[0])
                JetIsBtaggedWithSF_hist_1stjet_fwdeta.Fill(event.JetIsBtaggedWithSF[0])
        
        
            if event.nCleanedJets > 1 :
                JetPt_hist_2ndjet.Fill(event.JetPt[1])
                JetEta_hist_2ndjet.Fill(event.JetEta[1])
                JetPhi_hist_2ndjet.Fill(event.JetPhi[1])
                JetBTagger_hist_2ndjet.Fill(event.JetBTagger[1])
                JetIsBtagged_hist_2ndjet.Fill(event.JetIsBtagged[1])
                JetIsBtaggedWithSF_hist_2ndjet.Fill(event.JetIsBtaggedWithSF[1])

                if math.fabs(event.JetEta[1]) > 3.0 : 
                    JetPt_hist_2ndjet_fwdeta.Fill(event.JetPt[1])
                    JetEta_hist_2ndjet_fwdeta.Fill(event.JetEta[1])
                    JetPhi_hist_2ndjet_fwdeta.Fill(event.JetPhi[1])
                    JetBTagger_hist_2ndjet_fwdeta.Fill(event.JetBTagger[1])
                    JetIsBtagged_hist_2ndjet_fwdeta.Fill(event.JetIsBtagged[1])
                    JetIsBtaggedWithSF_hist_2ndjet_fwdeta.Fill(event.JetIsBtaggedWithSF[1])
                    
            

        PFMET_hist.Fill(event.PFMET)

        
    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_DATA = TFile.Open("jetsMetDistrib_DATA_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_DATA.cd()

    
    nCleanedJets_hist.Write()
    nCleanedJetsPt30_hist.Write()
    nCleanedJetsPt30BTagged_hist.Write()
    nCleanedJetsPt30BTagged_bTagSF_hist.Write()
    
    JetPt_hist.Write()
    JetPt_hist_1stJet.Write()
    JetPt_hist_2ndjet.Write()
    JetPt_hist_1stJet_fwdeta.Write()
    JetPt_hist_2ndjet_fwdeta.Write()
     
    JetEta_hist.Write()
    JetEta_hist_1stjet.Write()
    JetEta_hist_2ndjet.Write()
    JetEta_hist_1stjet_fwdeta.Write()
    JetEta_hist_2ndjet_fwdeta.Write()
                            
    JetPhi_hist.Write()
    JetPhi_hist_1stjet.Write()
    JetPhi_hist_2ndjet.Write()
    JetPhi_hist_1stjet_fwdeta.Write()
    JetPhi_hist_2ndjet_fwdeta.Write()
                            
    JetBTagger_hist.Write()
    JetBTagger_hist_1stjet.Write()
    JetBTagger_hist_2ndjet.Write()
    JetBTagger_hist_1stjet_fwdeta.Write()
    JetBTagger_hist_2ndjet_fwdeta.Write()
                            
    JetIsBtagged_hist.Write()
    JetIsBtagged_hist_1stjet.Write()
    JetIsBtagged_hist_2ndjet.Write()
    JetIsBtagged_hist_1stjet_fwdeta.Write()
    JetIsBtagged_hist_2ndjet_fwdeta.Write()
                            
    JetIsBtaggedWithSF_hist.Write()   
    JetIsBtaggedWithSF_hist_1stjet.Write()
    JetIsBtaggedWithSF_hist_2ndjet.Write()
    JetIsBtaggedWithSF_hist_1stjet_fwdeta.Write()
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta.Write()

    PFMET_hist.Write()

    outFile_DATA.Close()
    print "DATA histo file created!"


# ********************
#  do MC DY histos
# ********************
if(redoMCDYHistos) :

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    # define data histograms 
    nCleanedJets_hist_MC_DY                   = TH1F('nCleanedJets_MC_DY',                   'nCleanedJets_MC_DY',                   30, 0, 30) 
    nCleanedJetsPt30_hist_MC_DY               = TH1F('nCleanedJetsPt30_MC_DY',               'nCleanedJetsPt30_MC_DY',               15, 0, 15) 
    nCleanedJetsPt30BTagged_hist_MC_DY        = TH1F('nCleanedJetsPt30BTagged_MC_DY',        'nCleanedJetsPt30BTagged_MC_DY',        5,  0, 5)
    nCleanedJetsPt30BTagged_bTagSF_hist_MC_DY = TH1F('nCleanedJetsPt30BTagged_bTagSF_MC_DY', 'nCleanedJetsPt30BTagged_bTagSF_MC_DY', 5,  0, 5)

    JetPt_hist_MC_DY               = TH1F('JetPt_inclusive_MC_DY',           'JetPt_inclusive_MC_DY',           120, 0, 1200)
    JetPt_hist_1stJet_MC_DY        = TH1F('JetPt_leadingJet_MC_DY',          'JetPt_leadingJet_MC_DY',          120, 0, 1200)
    JetPt_hist_2ndjet_MC_DY        = TH1F('JetPt_subLeadingJet_MC_DY',       'JetPt_subLeadingJet_MC_DY',       120, 0, 1200)
    JetPt_hist_1stJet_fwdeta_MC_DY = TH1F('JetPt_leadingJet_fwdeta_MC_DY',   'JetPt_leadingJet_fwdeta_MC_DY',   120, 0, 1200)
    JetPt_hist_2ndjet_fwdeta_MC_DY = TH1F('JetPt_subLeadingJet_fwdeta_MC_DY','JetPt_subLeadingJet_fwdeta_MC_DY',120, 0, 1200)

    JetEta_hist_MC_DY               = TH1F('JetEta_inclusive_MC_DY',           'JetEta_inclusive_MC_DY',           50, -5, 5)
    JetEta_hist_1stjet_MC_DY        = TH1F('JetEta_leadingJet_MC_DY',          'JetEta_leadingJet_MC_DY',          50, -5, 5)
    JetEta_hist_2ndjet_MC_DY        = TH1F('JetEta_subLeadingJet_MC_DY',       'JetEta_subLeadingJet_MC_DY',       50, -5, 5)
    JetEta_hist_1stjet_fwdeta_MC_DY = TH1F('JetEta_leadingJet_fwdeta_MC_DY',   'JetEta_leadingJet_fwdeta_MC_DY',   50, -5, 5)
    JetEta_hist_2ndjet_fwdeta_MC_DY = TH1F('JetEta_subLeadingJet_fwdeta_MC_DY','JetEta_subLeadingJet_fwdeta_MC_DY',50, -5, 5)

    JetPhi_hist_MC_DY               = TH1F('JetPhi_inclusive_MC_DY',           'JetPhi_inclusive_MC_DY',           33, -3.3, 3.3)
    JetPhi_hist_1stjet_MC_DY        = TH1F('JetPhi_leadingJet_MC_DY',          'JetPhi_leadingJet_MC_DY',          33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_MC_DY        = TH1F('JetPhi_subLeadingJet_MC_DY',       'JetPhi_subLeadingJet_MC_DY',       33, -3.3, 3.3) 
    JetPhi_hist_1stjet_fwdeta_MC_DY = TH1F('JetPhi_leadingJet_fwdeta_MC_DY',   'JetPhi_leadingJet_fwdeta_MC_DY',   33, -3.3, 3.3) 
    JetPhi_hist_2ndjet_fwdeta_MC_DY = TH1F('JetPhi_subLeadingJet_fwdeta_MC_DY','JetPhi_subLeadingJet_fwdeta_MC_DY',33, -3.3, 3.3) 

    JetBTagger_hist_MC_DY               = TH1F('JetBTagger_inclusive_MC_DY',           'JetBTagger_inclusive_MC_DY',           20, -11, 2)
    JetBTagger_hist_1stjet_MC_DY        = TH1F('JetBTagger_leadingJet_MC_DY',          'JetBTagger_leadingJet_MC_DY',          20, -11, 2)
    JetBTagger_hist_2ndjet_MC_DY        = TH1F('JetBTagger_subLeadingJet_MC_DY',       'JetBTagger_subLeadingJet_MC_DY',       20, -11, 2)
    JetBTagger_hist_1stjet_fwdeta_MC_DY = TH1F('JetBTagger_leadingJet_fwdeta_MC_DY',   'JetBTagger_leadingJet_fwdeta_MC_DY',   20, -11, 2)
    JetBTagger_hist_2ndjet_fwdeta_MC_DY = TH1F('JetBTagger_subLeadingJet_fwdeta_MC_DY','JetBTagger_subLeadingJet_fwdeta_MC_DY',20, -11, 2)

    JetIsBtagged_hist_MC_DY               = TH1F('JetIsBtagged_inclusive_MC_DY',           'JetIsBtagged_inclusive_MC_DY',           11, 0, 1.1)
    JetIsBtagged_hist_1stjet_MC_DY        = TH1F('JetIsBtagged_leadingJet_MC_DY',          'JetIsBtagged_leadingJet_MC_DY',          11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_MC_DY        = TH1F('JetIsBtagged_subLeadingJet_MC_DY',       'JetIsBtagged_subLeadingJet_MC_DY',       11, 0, 1.1)
    JetIsBtagged_hist_1stjet_fwdeta_MC_DY = TH1F('JetIsBtagged_leadingJet_fwdeta_MC_DY',   'JetIsBtagged_leadingJet_fwdeta_MC_DY',   11, 0, 1.1)
    JetIsBtagged_hist_2ndjet_fwdeta_MC_DY = TH1F('JetIsBtagged_subLeadingJet_fwdeta_MC_DY','JetIsBtagged_subLeadingJet_fwdeta_MC_DY',11, 0, 1.1)
    
    JetIsBtaggedWithSF_hist_MC_DY               = TH1F('JetIsBtaggedWithSF_inclusive_MC_DY',           'JetIsBtaggedWithSF_inclusive_MC_DY',           11, 0, 1.1)
    JetIsBtaggedWithSF_hist_1stjet_MC_DY        = TH1F('JetIsBtaggedWithSF_leadingJet_MC_DY',          'JetIsBtaggedWithSF_leadingJet_MC_DY',          11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_MC_DY        = TH1F('JetIsBtaggedWithSF_subLeadingJet_MC_DY',       'JetIsBtaggedWithSF_subLeadingJet_MC_DY',       11, 0, 1.1) 
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_DY = TH1F('JetIsBtaggedWithSF_leadingJet_fwdeta_MC_DY',   'JetIsBtaggedWithSF_leadingJet_fwdeta_MC_DY',   11, 0, 1.1)
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_DY = TH1F('JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_DY','JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_DY',11, 0, 1.1)

    PFMET_hist_MC_DY = TH1F('PFMET_MC_DY','PFMET_MC_DY',250,0,250)


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


        nCleanedJets_hist_MC_DY.Fill(event.nCleanedJets,weight)
        nCleanedJetsPt30_hist_MC_DY.Fill(event.nCleanedJetsPt30,weight)
        nCleanedJetsPt30BTagged_hist_MC_DY.Fill(event.nCleanedJetsPt30BTagged,weight)
        nCleanedJetsPt30BTagged_bTagSF_hist_MC_DY.Fill(event.nCleanedJetsPt30BTagged_bTagSF,weight)
        
        for i in range(len(event.JetPt)) :
            JetPt_hist_MC_DY.Fill(event.JetPt[i],weight)
                  
        for i in range(len(event.JetEta)) :  
            JetEta_hist_MC_DY.Fill(event.JetEta[i],weight)
                                        
        for i in range(len(event.JetPhi)) :
            JetPhi_hist_MC_DY.Fill(event.JetPhi[i],weight)
          
        for i in range(len(event.JetBTagger)) :
            JetBTagger_hist_MC_DY.Fill(event.JetBTagger[i],weight)
                                        
        for i in range(len(event.JetIsBtagged)) :
            JetIsBtagged_hist_MC_DY.Fill(event.JetIsBtagged[i],weight)
        
        for i in range(len(event.JetIsBtaggedWithSF)) :                
            JetIsBtaggedWithSF_hist_MC_DY.Fill(event.JetIsBtaggedWithSF[i],weight)     
        

        if event.nCleanedJets > 0 :
            JetPt_hist_1stJet_MC_DY.Fill(event.JetPt[0],weight)
            JetEta_hist_1stjet_MC_DY.Fill(event.JetEta[0],weight)
            JetPhi_hist_1stjet_MC_DY.Fill(event.JetPhi[0],weight)
            JetBTagger_hist_1stjet_MC_DY.Fill(event.JetBTagger[0],weight)
            JetIsBtagged_hist_1stjet_MC_DY.Fill(event.JetIsBtagged[0],weight)
            JetIsBtaggedWithSF_hist_1stjet_MC_DY.Fill(event.JetIsBtaggedWithSF[0],weight)

            if math.fabs(event.JetEta[0]) > 3.0 : 
                JetPt_hist_1stJet_fwdeta_MC_DY.Fill(event.JetPt[0],weight)
                JetEta_hist_1stjet_fwdeta_MC_DY.Fill(event.JetEta[0],weight)
                JetPhi_hist_1stjet_fwdeta_MC_DY.Fill(event.JetPhi[0],weight)
                JetBTagger_hist_1stjet_fwdeta_MC_DY.Fill(event.JetBTagger[0],weight)
                JetIsBtagged_hist_1stjet_fwdeta_MC_DY.Fill(event.JetIsBtagged[0],weight)
                JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_DY.Fill(event.JetIsBtaggedWithSF[0],weight)
        
        
            if event.nCleanedJets > 1 :
                JetPt_hist_2ndjet_MC_DY.Fill(event.JetPt[1],weight)
                JetEta_hist_2ndjet_MC_DY.Fill(event.JetEta[1],weight)
                JetPhi_hist_2ndjet_MC_DY.Fill(event.JetPhi[1],weight)
                JetBTagger_hist_2ndjet_MC_DY.Fill(event.JetBTagger[1],weight)
                JetIsBtagged_hist_2ndjet_MC_DY.Fill(event.JetIsBtagged[1],weight)
                JetIsBtaggedWithSF_hist_2ndjet_MC_DY.Fill(event.JetIsBtaggedWithSF[1],weight)

                if math.fabs(event.JetEta[1]) > 3.0 : 
                    JetPt_hist_2ndjet_fwdeta_MC_DY.Fill(event.JetPt[1],weight)
                    JetEta_hist_2ndjet_fwdeta_MC_DY.Fill(event.JetEta[1],weight)
                    JetPhi_hist_2ndjet_fwdeta_MC_DY.Fill(event.JetPhi[1],weight)
                    JetBTagger_hist_2ndjet_fwdeta_MC_DY.Fill(event.JetBTagger[1],weight)
                    JetIsBtagged_hist_2ndjet_fwdeta_MC_DY.Fill(event.JetIsBtagged[1],weight)
                    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_DY.Fill(event.JetIsBtaggedWithSF[1],weight)
                    
            

        PFMET_hist_MC_DY.Fill(event.PFMET,weight)

        
    #save histograms in a root file 
    print "saving histograms into root file ..."
    outFile_MCDY = TFile.Open("jetsMetDistrib_MC_DY_"+ period + "_" + treeText +".root", "RECREATE")
    outFile_MCDY.cd()

    
    nCleanedJets_hist_MC_DY.Write()
    nCleanedJetsPt30_hist_MC_DY.Write()
    nCleanedJetsPt30BTagged_hist_MC_DY.Write()
    nCleanedJetsPt30BTagged_bTagSF_hist_MC_DY.Write()
    
    JetPt_hist_MC_DY.Write()
    JetPt_hist_1stJet_MC_DY.Write()
    JetPt_hist_2ndjet_MC_DY.Write()
    JetPt_hist_1stJet_fwdeta_MC_DY.Write()
    JetPt_hist_2ndjet_fwdeta_MC_DY.Write()
                            
    JetEta_hist_MC_DY.Write()
    JetEta_hist_1stjet_MC_DY.Write()
    JetEta_hist_2ndjet_MC_DY.Write()
    JetEta_hist_1stjet_fwdeta_MC_DY.Write()
    JetEta_hist_2ndjet_fwdeta_MC_DY.Write()
                            
    JetPhi_hist_MC_DY.Write()
    JetPhi_hist_1stjet_MC_DY.Write()
    JetPhi_hist_2ndjet_MC_DY.Write()
    JetPhi_hist_1stjet_fwdeta_MC_DY.Write()
    JetPhi_hist_2ndjet_fwdeta_MC_DY.Write()
                            
    JetBTagger_hist_MC_DY.Write()
    JetBTagger_hist_1stjet_MC_DY.Write()
    JetBTagger_hist_2ndjet_MC_DY.Write()
    JetBTagger_hist_1stjet_fwdeta_MC_DY.Write()
    JetBTagger_hist_2ndjet_fwdeta_MC_DY.Write()
                            
    JetIsBtagged_hist_MC_DY.Write()
    JetIsBtagged_hist_1stjet_MC_DY.Write()
    JetIsBtagged_hist_2ndjet_MC_DY.Write()
    JetIsBtagged_hist_1stjet_fwdeta_MC_DY.Write()
    JetIsBtagged_hist_2ndjet_fwdeta_MC_DY.Write()
                            
    JetIsBtaggedWithSF_hist_MC_DY.Write()   
    JetIsBtaggedWithSF_hist_1stjet_MC_DY.Write()
    JetIsBtaggedWithSF_hist_2ndjet_MC_DY.Write()
    JetIsBtaggedWithSF_hist_1stjet_fwdeta_MC_DY.Write()
    JetIsBtaggedWithSF_hist_2ndjet_fwdeta_MC_DY.Write()

    PFMET_hist_MC_DY.Write()

    outFile_MCDY.Close()
    print "MC DY histo file created!"
# ********************



# create output directory 
outputDir = "jetsMetDistrib_DATAvsMC_" + str(period) + "_" + str(treeText)
gSystem.Exec("mkdir -p " + outputDir)
print "Output directory created!"


# **************************
# read data histos from file 
histoDATA_input = TFile.Open("jetsMetDistrib_DATA_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoDATA_input.GetName(),'...'

inDATA_list = []

inDATA_list.append(histoDATA_input.Get('nCleanedJets'))
inDATA_list.append(histoDATA_input.Get('nCleanedJetsPt30'))
inDATA_list.append(histoDATA_input.Get('nCleanedJetsPt30BTagged'))
inDATA_list.append(histoDATA_input.Get('nCleanedJetsPt30BTagged_bTagSF'))
inDATA_list.append(histoDATA_input.Get('JetPt_inclusive'))
inDATA_list.append(histoDATA_input.Get('JetPt_leadingJet'))
inDATA_list.append(histoDATA_input.Get('JetPt_subLeadingJet'))
inDATA_list.append(histoDATA_input.Get('JetPt_leadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetPt_subLeadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetEta_inclusive'))
inDATA_list.append(histoDATA_input.Get('JetEta_leadingJet'))
inDATA_list.append(histoDATA_input.Get('JetEta_subLeadingJet'))
inDATA_list.append(histoDATA_input.Get('JetEta_leadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetEta_subLeadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetPhi_inclusive'))
inDATA_list.append(histoDATA_input.Get('JetPhi_leadingJet'))
inDATA_list.append(histoDATA_input.Get('JetPhi_subLeadingJet'))
inDATA_list.append(histoDATA_input.Get('JetPhi_leadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetPhi_subLeadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_inclusive'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_leadingJet'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_subLeadingJet'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_leadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetBTagger_subLeadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_inclusive'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_leadingJet'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_subLeadingJet'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_leadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetIsBtagged_subLeadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_inclusive'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_leadingJet'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_subLeadingJet'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_leadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('JetIsBtaggedWithSF_subLeadingJet_fwdeta'))
inDATA_list.append(histoDATA_input.Get('PFMET'))


# ****************************
# read DY MC histos from file 
histoMCDY_input = TFile.Open("jetsMetDistrib_MC_DY_" + str(period) + "_" + str(treeText) + ".root")
print 'Reading file', histoMCDY_input.GetName(),'...'

inMCDY_list = []

inMCDY_list.append(histoMCDY_input.Get('nCleanedJets_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('nCleanedJetsPt30_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('nCleanedJetsPt30BTagged_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('nCleanedJetsPt30BTagged_bTagSF_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_inclusive_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_leadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_subLeadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_leadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetPt_subLeadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_inclusive_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_leadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_subLeadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_leadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetEta_subLeadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_inclusive_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_leadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_subLeadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_leadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetPhi_subLeadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_inclusive_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_leadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_subLeadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_leadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetBTagger_subLeadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_inclusive_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_leadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_subLeadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_leadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtagged_subLeadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_inclusive_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_leadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_subLeadingJet_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_leadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('JetIsBtaggedWithSF_subLeadingJet_fwdeta_MC_DY'))
inMCDY_list.append(histoMCDY_input.Get('PFMET_MC_DY'))

       

# ******************************
# do DATA vs MC comparison plots  
for i in range(len(inDATA_list)) : 

    canvas = TCanvas("canvas","canvas",800,800)

    #DATA hist
    inDATA_list[i].SetMarkerStyle(20)
    inDATA_list[i].SetMarkerSize(0.6)

    #MC hist
    inMCDY_list[i].SetFillColor(kOrange-3)
    inMCDY_list[i].SetLineColor(kBlack)


    #upper plot pad
    pad1 = TPad("pad1","pad1", 0, 0.3, 1, 1.0)
    pad1.Draw()
    pad1.cd()


    inMCDY_list[i].SetMaximum(1.3*max(inMCDY_list[i].GetMaximum(),inDATA_list[i].GetMaximum()))
    inDATA_list[i].SetMaximum(1.3*max(inMCDY_list[i].GetMaximum(),inDATA_list[i].GetMaximum()))
    

    inMCDY_list[i].Draw("histo") 
    inDATA_list[i].Draw("sameEP")
    
    
    inMCDY_list[i].SetTitle("")
    inMCDY_list[i].GetXaxis().SetTitle(inDATA_list[i].GetTitle())
    inMCDY_list[i].GetXaxis().SetLabelFont(43)
    inMCDY_list[i].GetXaxis().SetLabelSize(15)
    inMCDY_list[i].GetYaxis().SetTitleSize(20)
    inMCDY_list[i].GetYaxis().SetTitleFont(43)
    inMCDY_list[i].GetYaxis().SetTitleOffset(1.8)
    inMCDY_list[i].GetYaxis().SetLabelFont(43)
    inMCDY_list[i].GetYaxis().SetLabelSize(15)
    inMCDY_list[i].GetYaxis().SetTitle("Events")

    gStyle.SetOptStat(0)

    if "Pt" in inDATA_list[i].GetTitle() :
        pad1.SetLogy()


    # legend
    legend = TLegend(0.82,0.77,0.95,0.89)
    legend.AddEntry(inDATA_list[i],"Data", "p")
    legend.AddEntry(inMCDY_list[i],"DY MC","f")
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
    rp.Divide(TH1F(inMCDY_list[i]))   #divide histo rp/MC
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
