#!/usr/bin/env python

# *******************
# usage: 
#    python LepBDTdistrib_DATA.py
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
CRZLTree = False
ZTree    = True

# data periods options
# period = "data2016"
period = "data2017"
# *****************************



#input file
if(period == "data2016"):
    data = TFile.Open("/data3/Higgs/170222/AllData/ZZ4lAnalysis.root") #2016 data
    lumi = 35.9  # fb-1
    if(ZZTree):
        tree      = data.Get("ZZTree/candTree")
        treeText  = "ZZTree"
    elif(CRZLTree):
        tree      = data.Get("CRZLTree/candTree")
        treeText  = "CRZLTree"
    elif(ZTree):
        tree      = data.Get("ZTree/candTree")
        treeText  = "ZTree"
    else:
        print ("Error: wrong option!")

elif(period == "data2017"):
    data = TFile.Open("/data3/Higgs/180416/AllData/ZZ4lAnalysis.root") #2017 data
    lumi = 41.30   # fb-1
    if(ZZTree):
        tree      = data.Get("ZZTree/candTree")
        treeText  = "ZZTree"
    elif(CRZLTree):
        tree      = data.Get("CRZLTree/candTree")
        treeText  = "CRZLTree"
    elif(ZTree):
        tree      = data.Get("ZTree/candTree")
        treeText  = "ZTree"
    else:
        print ("Error: wrong option!")
else: 
    print ("Error: choose a period!")




if(redoHistos) :

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    # define LepBDT histograms 
    # Z->ee
    Z_ele_1st_LepBDT_hist    = TH1F('LepBDT leading ele'               , 'LepBDT leading ele'               , 100, -1.1, 1.1) 
    Z_ele_1st_LepBDT_hist_EB = TH1F('LepBDT leading ele in ECAL Barrel', 'LepBDT leading ele in ECAL Barrel', 100, -1.1, 1.1) 
    Z_ele_1st_LepBDT_hist_EE = TH1F('LepBDT leading ele in ECAL Endcap', 'LepBDT leading ele in ECAL Endcap', 100, -1.1, 1.1) 

    # Z->mumu
    # Z_mu_1st_LepBDT_hist    = TH1F('LepBDT leading mu'               , 'LepBDT leading mu'               , 100, -1.1, 1.1) 
    # Z_mu_1st_LepBDT_hist_MB = TH1F('LepBDT leading mu in Muon Barrel', 'LepBDT leading mu in Muon Barrel', 100, -1.1, 1.1) 
    # Z_mu_1st_LepBDT_hist_ME = TH1F('LepBDT leading mu in Muon Endcap', 'LepBDT leading mu in Muon Endcap', 100, -1.1, 1.1) 

    if CRZLTree :
        Z_ExtraEl_LepBDT_hist    = TH1F('LepBDT extraEl'               , 'LepBDT extraEl'               , 100, -1.1, 1.1)
        Z_ExtraEl_LepBDT_hist_EB = TH1F('LepBDT extraEl in ECAL Barrel', 'LepBDT extraEl in ECAL Barrel', 100, -1.1, 1.1)
        Z_ExtraEl_LepBDT_hist_EE = TH1F('LepBDT extraEl in ECAL Endcap', 'LepBDT extraEl in ECAL Endcap', 100, -1.1, 1.1)

        # Z_ExtraMu_LepBDT_hist    = TH1F('LepBDT extraMu'               , 'LepBDT extraMu'               , 100, -1.1, 1.1)
        # Z_ExtraMu_LepBDT_hist_MB = TH1F('LepBDT extraMu in Muon Barrel', 'LepBDT extraMu in Muon Barrel', 100, -1.1, 1.1)
        # Z_ExtraMu_LepBDT_hist_ME = TH1F('LepBDT extraMu in Muon Endcap', 'LepBDT extraMu in Muon Endcap', 100, -1.1, 1.1)
    
    

    

    tree.SetBranchStatus("*",0)  # disable all branches
    if ZTree :
        tree.SetBranchStatus("Zsel",1)
    else :
        tree.SetBranchStatus("ZZsel",1)
    tree.SetBranchStatus("LepLepId",1)
    tree.SetBranchStatus("LepPt",1)
    tree.SetBranchStatus("LepEta",1)
    tree.SetBranchStatus("LepBDT",1)


    # read TTree 
    print "reading tree", data.GetName(),treeText,tree.GetName()  ,"..."
    for event in tree:
        if ZTree :
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger
        else :
            if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger
        
        
        # Z->ee 
        if(int(math.fabs(event.LepLepId[0])) == 11 ) :
            
            if event.LepPt[0] >= event.LepPt[1] :
                Z_ele_1st_LepBDT_hist.Fill(event.LepBDT[0])

                if math.fabs(event.LepEta[0]) <= 1.479 :
                    Z_ele_1st_LepBDT_hist_EB.Fill(event.LepBDT[0])
                else :
                    Z_ele_1st_LepBDT_hist_EE.Fill(event.LepBDT[0])

            else :
                Z_ele_1st_LepBDT_hist.Fill(event.LepBDT[1])

                if math.fabs(event.LepEta[1]) <= 1.479 :
                    Z_ele_1st_LepBDT_hist_EB.Fill(event.LepBDT[1])
                else :
                    Z_ele_1st_LepBDT_hist_EE.Fill(event.LepBDT[1])

                      

        # Z->mumu
        # if(int(math.fabs(event.LepLepId[0])) == 13 ) :
            
        #     if event.LepPt[0] >= event.LepPt[1] :
        #         Z_mu_1st_LepBDT_hist.Fill(event.LepBDT[0])

        #         if math.fabs(event.LepEta[0]) <= 1. :
        #             Z_mu_1st_LepBDT_hist_MB.Fill(event.LepBDT[0])
        #         else :
        #             Z_mu_1st_LepBDT_hist_ME.Fill(event.LepBDT[0])

        #     else :
        #         Z_mu_1st_LepBDT_hist.Fill(event.LepBDT[1])

        #         if math.fabs(event.LepEta[1]) <= 1. :
        #             Z_mu_1st_LepBDT_hist_MB.Fill(event.LepBDT[1])
        #         else :
        #             Z_mu_1st_LepBDT_hist_ME.Fill(event.LepBDT[1])

            
        
        # extra lepton
        if CRZLTree :
                
                if(int(math.fabs(event.LepLepId[2])) == 11 ) :

                    Z_ExtraEl_LepBDT_hist.Fill(event.LepBDT[2])
                    
                    if math.fabs(event.LepEta[2]) <= 1.479 :  
                        Z_ExtraEl_LepBDT_hist_EB.Fill(event.LepBDT[2])
                    else :
                        Z_ExtraEl_LepBDT_hist_EE.Fill(event.LepBDT[2])

                # elif(int(math.fabs(event.LepLepId[2])) == 13 ) :
                    
                #     Z_ExtraMu_LepBDT_hist.Fill(event.LepBDT[2])

                #     if math.fabs(event.LepEta[2]) <= 1. :
                #         Z_ExtraMu_LepBDT_hist_MB.Fill(event.LepBDT[2])
                #     else :
                #         Z_ExtraMu_LepBDT_hist_ME.Fill(event.LepBDT[2])

                else :
                    print "Error: wrong particle ID!"


    #save histograms in a root file 
    print "saving histograms into root file ..."
    LepBDT_outFile = TFile.Open("LepBDTdistrib_DATA_"+ period + "_" + treeText +".root", "RECREATE")
    LepBDT_outFile.cd()

    Z_ele_1st_LepBDT_hist.Write()
    Z_ele_1st_LepBDT_hist_EB.Write()
    Z_ele_1st_LepBDT_hist_EE.Write()
                                          
    # Z_mu_1st_LepBDT_hist.Write()
    # Z_mu_1st_LepBDT_hist_MB.Write()
    # Z_mu_1st_LepBDT_hist_ME.Write()
                                                                
    if CRZLTree :               
        Z_ExtraEl_LepBDT_hist.Write()
        Z_ExtraEl_LepBDT_hist_EB.Write()
        Z_ExtraEl_LepBDT_hist_EE.Write()
                                
        # Z_ExtraMu_LepBDT_hist.Write()
        # Z_ExtraMu_LepBDT_hist_MB.Write()
        # Z_ExtraMu_LepBDT_hist_ME.Write()


    LepBDT_outFile.Close()
