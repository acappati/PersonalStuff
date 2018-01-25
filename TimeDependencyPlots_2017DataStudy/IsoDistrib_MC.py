#!/usr/bin/env python

# *******************
# usage: 
#    python IsoDistrib_MC.py
#
# *******************


import json
import ROOT, helper, math
from ROOT import TFile, TH1, TH1F, TCanvas, gSystem
from helper import DoSimpleFit, Result, DoDCBunbinnedFit

# *****************************
# Declare all the variables
# options
redoHistos = True

# data tree options 
ZZTree   = False
CRZLTree = True
ZTree    = False

# data periods options
# period = "data2016"
period = "data2017"
# *****************************

# lumi 
if(period == "data2016"):
    lumi = 35.9     # fb-1
elif(period == "data2017"):
    lumi = 41.86    # fb-1


#input file (DY MC)
if(period == "data2017"):
    inputTree = TFile.Open("/data3/Higgs/180122/DYJetsToLL_M50/ZZ4lAnalysis.root") #DYJets 2017 MC
    if(ZZTree):
        tree      = inputTree.Get("ZZTree/candTree")
        treeText  = "ZZTree"
    elif(CRZLTree):
        tree      = inputTree.Get("CRZLTree/candTree")
        treeText  = "CRZLTree"
    elif(ZTree):
        tree      = inputTree.Get("ZTree/candTree")
        treeText  = "ZTree"
    else:
        print ("Error: wrong option!")
else: 
    print ("Error: choose a period!")




if(redoHistos) :

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    # define ISO histograms 
    # Z->ee
    Z_ele_1st_ISO_hist    = TH1F('ISO leading ele'               , 'ISO leading ele'               , 100, 0, 0.6) 
    Z_ele_1st_ISO_hist_EB = TH1F('ISO leading ele in ECAL Barrel', 'ISO leading ele in ECAL Barrel', 100, 0, 0.6) 
    Z_ele_1st_ISO_hist_EE = TH1F('ISO leading ele in ECAL Endcap', 'ISO leading ele in ECAL Endcap', 100, 0, 0.6) 

    Z_ele_max_ISO_hist    = TH1F('ISO max ele'               , 'ISO max ele'               , 100, 0, 0.6) 
    Z_ele_max_ISO_hist_EB = TH1F('ISO max ele in ECAL Barrel', 'ISO max ele in ECAL Barrel', 100, 0, 0.6) 
    Z_ele_max_ISO_hist_EE = TH1F('ISO max ele in ECAL Endcap', 'ISO max ele in ECAL Endcap', 100, 0, 0.6) 

    # Z->mumu
    Z_mu_1st_ISO_hist    = TH1F('ISO leading mu'               , 'ISO leading mu'               , 100, 0, 0.6) 
    Z_mu_1st_ISO_hist_MB = TH1F('ISO leading mu in Muon Barrel', 'ISO leading mu in Muon Barrel', 100, 0, 0.6) 
    Z_mu_1st_ISO_hist_ME = TH1F('ISO leading mu in Muon Endcap', 'ISO leading mu in Muon Endcap', 100, 0, 0.6) 

    Z_mu_max_ISO_hist    = TH1F('ISO max mu'               , 'ISO max mu'               , 100, 0, 0.6) 
    Z_mu_max_ISO_hist_MB = TH1F('ISO max mu in Muon Barrel', 'ISO max mu in Muon Barrel', 100, 0, 0.6) 
    Z_mu_max_ISO_hist_ME = TH1F('ISO max mu in Muon Endcap', 'ISO max mu in Muon Endcap', 100, 0, 0.6) 
   
    if CRZLTree :
        Z_ExtraEl_ISO_hist    = TH1F('ISO extraEl'               , 'ISO extraEl'               , 100, 0, 2.)
        Z_ExtraEl_ISO_hist_EB = TH1F('ISO extraEl in ECAL Barrel', 'ISO extraEl in ECAL Barrel', 100, 0, 2.)
        Z_ExtraEl_ISO_hist_EE = TH1F('ISO extraEl in ECAL Endcap', 'ISO extraEl in ECAL Endcap', 100, 0, 2.)

        Z_ExtraMu_ISO_hist    = TH1F('ISO extraMu'               , 'ISO extraMu'               , 100, 0, 2.)
        Z_ExtraMu_ISO_hist_MB = TH1F('ISO extraMu in Muon Barrel', 'ISO extraMu in Muon Barrel', 100, 0, 2.)
        Z_ExtraMu_ISO_hist_ME = TH1F('ISO extraMu in Muon Endcap', 'ISO extraMu in Muon Endcap', 100, 0, 2.)
    
    

    # get partial event weight
    hcounters           = inputTree.Get("ZZTree/Counters")
    gen_sumWeights      = hcounters.GetBinContent(40)
    partialSampleWeight = lumi * 1000 / gen_sumWeights


    # read TTree 
    print "reading tree", inputTree.GetName(),treeText,tree.GetName()  ,"..."
    for event in tree:
        if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger 


        weight = partialSampleWeight*event.xsec*event.overallEventWeight

        
        # Z->ee 
        if(int(math.fabs(event.LepLepId[0])) == 11 ) :
            
            if event.LepPt[0] >= event.LepPt[1] :
                Z_ele_1st_ISO_hist.Fill(event.LepCombRelIsoPF[0], weight)

                if math.fabs(event.LepEta[0]) <= 1.479 :
                    Z_ele_1st_ISO_hist_EB.Fill(event.LepCombRelIsoPF[0], weight)
                else :
                    Z_ele_1st_ISO_hist_EE.Fill(event.LepCombRelIsoPF[0], weight)

            else :
                Z_ele_1st_ISO_hist.Fill(event.LepCombRelIsoPF[1], weight)

                if math.fabs(event.LepEta[1]) <= 1.479 :
                    Z_ele_1st_ISO_hist_EB.Fill(event.LepCombRelIsoPF[1], weight)
                else :
                    Z_ele_1st_ISO_hist_EE.Fill(event.LepCombRelIsoPF[1], weight)

            
            if(event.LepCombRelIsoPF[0] >= event.LepCombRelIsoPF[1]):
                Z_ele_max_ISO_hist.Fill(event.LepCombRelIsoPF[0], weight)

                if math.fabs(event.LepEta[0]) <= 1.479 :
                    Z_ele_max_ISO_hist_EB.Fill(event.LepCombRelIsoPF[0], weight)
                else :
                    Z_ele_max_ISO_hist_EE.Fill(event.LepCombRelIsoPF[0], weight)

            else :
                Z_ele_max_ISO_hist.Fill(event.LepCombRelIsoPF[1], weight)
                        
                if math.fabs(event.LepEta[1]) <= 1.479 :
                    Z_ele_max_ISO_hist_EB.Fill(event.LepCombRelIsoPF[1], weight)
                else :
                    Z_ele_max_ISO_hist_EE.Fill(event.LepCombRelIsoPF[1], weight)

          

        # Z->mumu
        if(int(math.fabs(event.LepLepId[0])) == 13 ) :
            
            if event.LepPt[0] >= event.LepPt[1] :
                Z_mu_1st_ISO_hist.Fill(event.LepCombRelIsoPF[0], weight)

                if math.fabs(event.LepEta[0]) <= 1. :
                    Z_mu_1st_ISO_hist_MB.Fill(event.LepCombRelIsoPF[0], weight)
                else :
                    Z_mu_1st_ISO_hist_ME.Fill(event.LepCombRelIsoPF[0], weight)

            else :
                Z_mu_1st_ISO_hist.Fill(event.LepCombRelIsoPF[1], weight)

                if math.fabs(event.LepEta[1]) <= 1. :
                    Z_mu_1st_ISO_hist_MB.Fill(event.LepCombRelIsoPF[1], weight)
                else :
                    Z_mu_1st_ISO_hist_ME.Fill(event.LepCombRelIsoPF[1], weight)

            
            if(event.LepCombRelIsoPF[0] >= event.LepCombRelIsoPF[1]):
                Z_mu_max_ISO_hist.Fill(event.LepCombRelIsoPF[0], weight)

                if math.fabs(event.LepEta[0]) <= 1. :
                    Z_mu_max_ISO_hist_MB.Fill(event.LepCombRelIsoPF[0], weight)
                else :
                    Z_mu_max_ISO_hist_ME.Fill(event.LepCombRelIsoPF[0], weight)

            else :
                Z_mu_max_ISO_hist.Fill(event.LepCombRelIsoPF[1], weight)
                        
                if math.fabs(event.LepEta[1]) <= 1. :
                    Z_mu_max_ISO_hist_MB.Fill(event.LepCombRelIsoPF[1], weight)
                else :
                    Z_mu_max_ISO_hist_ME.Fill(event.LepCombRelIsoPF[1], weight)

        
        # extra lepton
        if CRZLTree :
                
                if(int(math.fabs(event.LepLepId[2])) == 11 ) :

                    Z_ExtraEl_ISO_hist.Fill(event.LepCombRelIsoPF[2], weight)
                    
                    if math.fabs(event.LepEta[2]) <= 1.479 :  
                        Z_ExtraEl_ISO_hist_EB.Fill(event.LepCombRelIsoPF[2], weight)
                    else :
                        Z_ExtraEl_ISO_hist_EE.Fill(event.LepCombRelIsoPF[2], weight)

                elif(int(math.fabs(event.LepLepId[2])) == 13 ) :
                    
                    Z_ExtraMu_ISO_hist.Fill(event.LepCombRelIsoPF[2], weight)

                    if math.fabs(event.LepEta[2]) <= 1. :
                        Z_ExtraMu_ISO_hist_MB.Fill(event.LepCombRelIsoPF[2], weight)
                    else :
                        Z_ExtraMu_ISO_hist_ME.Fill(event.LepCombRelIsoPF[2], weight)

                else :
                    print "Error: wrong particle ID!"


    #save histograms in a root file 
    print "saving histograms into root file ..."
    Iso_outFile = TFile.Open("IsoDistrib_MC_"+ period + "_" + treeText +".root", "RECREATE")
    Iso_outFile.cd()

    Z_ele_1st_ISO_hist.Write()
    Z_ele_1st_ISO_hist_EB.Write()
    Z_ele_1st_ISO_hist_EE.Write()
                                
    Z_ele_max_ISO_hist.Write()
    Z_ele_max_ISO_hist_EB.Write()
    Z_ele_max_ISO_hist_EE.Write()
                                                   
    Z_mu_1st_ISO_hist.Write()
    Z_mu_1st_ISO_hist_MB.Write()
    Z_mu_1st_ISO_hist_ME.Write()
                                
    Z_mu_max_ISO_hist.Write()
    Z_mu_max_ISO_hist_MB.Write()
    Z_mu_max_ISO_hist_ME.Write()
                                
    if CRZLTree :               
        Z_ExtraEl_ISO_hist.Write()
        Z_ExtraEl_ISO_hist_EB.Write()
        Z_ExtraEl_ISO_hist_EE.Write()
                                
        Z_ExtraMu_ISO_hist.Write()
        Z_ExtraMu_ISO_hist_MB.Write()
        Z_ExtraMu_ISO_hist_ME.Write()


    Iso_outFile.Close()
