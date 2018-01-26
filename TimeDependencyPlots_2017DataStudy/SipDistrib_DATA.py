#!/usr/bin/env python

# *******************
# usage: 
#    python SipDistrib_DATA.py
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
    data = TFile.Open("/data3/Higgs/180122/AllData/ZZ4lAnalysis.root") #2017 data
    lumi = 41.86   # fb-1
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

    # define SIP histograms 
    # Z->ee
    Z_ele_1st_SIP_hist    = TH1F('SIP leading ele'               , 'SIP leading ele'               , 100, 0, 10) 
    Z_ele_1st_SIP_hist_EB = TH1F('SIP leading ele in ECAL Barrel', 'SIP leading ele in ECAL Barrel', 100, 0, 10) 
    Z_ele_1st_SIP_hist_EE = TH1F('SIP leading ele in ECAL Endcap', 'SIP leading ele in ECAL Endcap', 100, 0, 10) 

    Z_ele_max_SIP_hist    = TH1F('SIP max ele'               , 'SIP max ele'               , 100, 0, 10) 
    Z_ele_max_SIP_hist_EB = TH1F('SIP max ele in ECAL Barrel', 'SIP max ele in ECAL Barrel', 100, 0, 10) 
    Z_ele_max_SIP_hist_EE = TH1F('SIP max ele in ECAL Endcap', 'SIP max ele in ECAL Endcap', 100, 0, 10) 

    # Z->mumu
    Z_mu_1st_SIP_hist    = TH1F('SIP leading mu'               , 'SIP leading mu'               , 100, 0, 10) 
    Z_mu_1st_SIP_hist_MB = TH1F('SIP leading mu in Muon Barrel', 'SIP leading mu in Muon Barrel', 100, 0, 10) 
    Z_mu_1st_SIP_hist_ME = TH1F('SIP leading mu in Muon Endcap', 'SIP leading mu in Muon Endcap', 100, 0, 10) 

    Z_mu_max_SIP_hist    = TH1F('SIP max mu'               , 'SIP max mu'               , 100, 0, 10) 
    Z_mu_max_SIP_hist_MB = TH1F('SIP max mu in Muon Barrel', 'SIP max mu in Muon Barrel', 100, 0, 10) 
    Z_mu_max_SIP_hist_ME = TH1F('SIP max mu in Muon Endcap', 'SIP max mu in Muon Endcap', 100, 0, 10) 
   
    if not ZTree :
        Z_ExtraEl_SIP_hist    = TH1F('SIP extraEl'               , 'SIP extraEl'               , 100, 0, 10)
        Z_ExtraEl_SIP_hist_EB = TH1F('SIP extraEl in ECAL Barrel', 'SIP extraEl in ECAL Barrel', 100, 0, 10)
        Z_ExtraEl_SIP_hist_EE = TH1F('SIP extraEl in ECAL Endcap', 'SIP extraEl in ECAL Endcap', 100, 0, 10)

        Z_ExtraMu_SIP_hist    = TH1F('SIP extraMu'               , 'SIP extraMu'               , 100, 0, 10)
        Z_ExtraMu_SIP_hist_MB = TH1F('SIP extraMu in Muon Barrel', 'SIP extraMu in Muon Barrel', 100, 0, 10)
        Z_ExtraMu_SIP_hist_ME = TH1F('SIP extraMu in Muon Endcap', 'SIP extraMu in Muon Endcap', 100, 0, 10)
    
    

    



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
                Z_ele_1st_SIP_hist.Fill(event.LepSIP[0])

                if math.fabs(event.LepEta[0]) <= 1.479 :
                    Z_ele_1st_SIP_hist_EB.Fill(event.LepSIP[0])
                else :
                    Z_ele_1st_SIP_hist_EE.Fill(event.LepSIP[0])

            else :
                Z_ele_1st_SIP_hist.Fill(event.LepSIP[1])

                if math.fabs(event.LepEta[1]) <= 1.479 :
                    Z_ele_1st_SIP_hist_EB.Fill(event.LepSIP[1])
                else :
                    Z_ele_1st_SIP_hist_EE.Fill(event.LepSIP[1])

            
            if(event.LepSIP[0] >= event.LepSIP[1]):
                Z_ele_max_SIP_hist.Fill(event.LepSIP[0])

                if math.fabs(event.LepEta[0]) <= 1.479 :
                    Z_ele_max_SIP_hist_EB.Fill(event.LepSIP[0])
                else :
                    Z_ele_max_SIP_hist_EE.Fill(event.LepSIP[0])

            else :
                Z_ele_max_SIP_hist.Fill(event.LepSIP[1])
                        
                if math.fabs(event.LepEta[1]) <= 1.479 :
                    Z_ele_max_SIP_hist_EB.Fill(event.LepSIP[1])
                else :
                    Z_ele_max_SIP_hist_EE.Fill(event.LepSIP[1])

          

        # Z->mumu
        if(int(math.fabs(event.LepLepId[0])) == 13 ) :
            
            if event.LepPt[0] >= event.LepPt[1] :
                Z_mu_1st_SIP_hist.Fill(event.LepSIP[0])

                if math.fabs(event.LepEta[0]) <= 1. :
                    Z_mu_1st_SIP_hist_MB.Fill(event.LepSIP[0])
                else :
                    Z_mu_1st_SIP_hist_ME.Fill(event.LepSIP[0])

            else :
                Z_mu_1st_SIP_hist.Fill(event.LepSIP[1])

                if math.fabs(event.LepEta[1]) <= 1. :
                    Z_mu_1st_SIP_hist_MB.Fill(event.LepSIP[1])
                else :
                    Z_mu_1st_SIP_hist_ME.Fill(event.LepSIP[1])

            
            if(event.LepSIP[0] >= event.LepSIP[1]):
                Z_mu_max_SIP_hist.Fill(event.LepSIP[0])

                if math.fabs(event.LepEta[0]) <= 1. :
                    Z_mu_max_SIP_hist_MB.Fill(event.LepSIP[0])
                else :
                    Z_mu_max_SIP_hist_ME.Fill(event.LepSIP[0])

            else :
                Z_mu_max_SIP_hist.Fill(event.LepSIP[1])
                        
                if math.fabs(event.LepEta[1]) <= 1. :
                    Z_mu_max_SIP_hist_MB.Fill(event.LepSIP[1])
                else :
                    Z_mu_max_SIP_hist_ME.Fill(event.LepSIP[1])

        
        # extra lepton
        if not ZTree :
                
                if(int(math.fabs(event.LepLepId[2])) == 11 ) :

                    Z_ExtraEl_SIP_hist.Fill(event.LepSIP[2])
                    
                    if math.fabs(event.LepEta[2]) <= 1.479 :  
                        Z_ExtraEl_SIP_hist_EB.Fill(event.LepSIP[2])
                    else :
                        Z_ExtraEl_SIP_hist_EE.Fill(event.LepSIP[2])

                elif(int(math.fabs(event.LepLepId[2])) == 13 ) :
                    
                    Z_ExtraMu_SIP_hist.Fill(event.LepSIP[2])

                    if math.fabs(event.LepEta[2]) <= 1. :
                        Z_ExtraMu_SIP_hist_MB.Fill(event.LepSIP[2])
                    else :
                        Z_ExtraMu_SIP_hist_ME.Fill(event.LepSIP[2])

                else :
                    print "Error: wrong particle ID!"


    #save histograms in a root file 
    print "saving histograms into root file ..."
    Sip_outFile = TFile.Open("SipDistrib_DATA_"+ period + "_" + treeText +".root", "RECREATE")
    Sip_outFile.cd()

    Z_ele_1st_SIP_hist.Write()
    Z_ele_1st_SIP_hist_EB.Write()
    Z_ele_1st_SIP_hist_EE.Write()
                                
    Z_ele_max_SIP_hist.Write()
    Z_ele_max_SIP_hist_EB.Write()
    Z_ele_max_SIP_hist_EE.Write()
                                                   
    Z_mu_1st_SIP_hist.Write()
    Z_mu_1st_SIP_hist_MB.Write()
    Z_mu_1st_SIP_hist_ME.Write()
                                
    Z_mu_max_SIP_hist.Write()
    Z_mu_max_SIP_hist_MB.Write()
    Z_mu_max_SIP_hist_ME.Write()
                                
    if not ZTree :               
        Z_ExtraEl_SIP_hist.Write()
        Z_ExtraEl_SIP_hist_EB.Write()
        Z_ExtraEl_SIP_hist_EE.Write()
                                
        Z_ExtraMu_SIP_hist.Write()
        Z_ExtraMu_SIP_hist_MB.Write()
        Z_ExtraMu_SIP_hist_ME.Write()


    Sip_outFile.Close()
