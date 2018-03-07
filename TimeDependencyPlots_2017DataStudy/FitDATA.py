#!/usr/bin/env python


# *******************
# usage: 
#    python FitDATA.py
#
# structure:
#    - read Data file (inputTree) 
#    - store data events in histos and save them in a file.root (if redoHistos)
#    - read the root file 
#    - fit the histos with a function from helper.py and print on file fit results
# ********************

import json
import ROOT, helper, math
from ROOT import TFile, TH1, TH1F, TCanvas, gSystem
from helper import DoSimpleFit, Result, DoDCBunbinnedFit


# *****************************
# Declare all the variables
# fit options 
fitMC      = False  #true for fitting MC in FitMC.py 
fitDATA    = True   #true for fitting DATA in FitDATA.py 
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
    data = TFile.Open("/data3/Higgs/180305/AllData/ZZ4lAnalysis.root") #2017 data
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


# output directory 
outputDir = "FitResults_DATA_" + str(period) + "_" + str(treeText)
gSystem.Exec("mkdir -p " + outputDir)
print "Output directory created!"



if(redoHistos) : 

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    #define histogramms 
    #Z->ee histos
    ZMass_ele_hist      = TH1F( 'ZMass_ele'      , 'ZMass_ele'      , 120, 60, 120)  #ZMass , Z->ee
    ZMass_ele_hist_EBEB = TH1F( 'ZMass_ele_EBEB' , 'ZMass_ele_EBEB' , 120, 60, 120)  #ZMass , Z->ee , Barrel-Barrel
    ZMass_ele_hist_EBEE = TH1F( 'ZMass_ele_EBEE' , 'ZMass_ele_EBEE' , 120, 60, 120)  #ZMass , Z->ee , Barrel-Endcap
    ZMass_ele_hist_EEEE = TH1F( 'ZMass_ele_EEEE' , 'ZMass_ele_EEEE' , 120, 60, 120)  #ZMass , Z->ee , Endcap-Endcap
    if not ZTree :
        ZMass_ele_hist_extraMu = TH1F( 'ZMass_ele_extraMu', 'ZMass_ele_extraMu', 120, 60, 120)  #ZMass , Z->ee + Extra mu
        ZMass_ele_hist_extraEl = TH1F( 'ZMass_ele_extraEl', 'ZMass_ele_extraEl', 120, 60, 120)  #ZMass , Z->ee + Extra e
    #Z->mumu histos
    ZMass_mu_hist      = TH1F( 'ZMass_mu'      , 'ZMass_mu'      , 120, 60, 120)  #ZMass , Z->mumu
    ZMass_mu_hist_MBMB = TH1F( 'ZMass_mu_MBMB' , 'ZMass_mu_MBMB' , 120, 60, 120)  #ZMass , Z->mumu , Barrel-Barrel
    ZMass_mu_hist_MBME = TH1F( 'ZMass_mu_MBME' , 'ZMass_mu_MBME' , 120, 60, 120)  #ZMass , Z->mumu , Barrel-Endcap
    ZMass_mu_hist_MEME = TH1F( 'ZMass_mu_MEME' , 'ZMass_mu_MEME' , 120, 60, 120)  #ZMass , Z->mumu , Endcap-Endcap
    if not ZTree :
        ZMass_mu_hist_extraMu  = TH1F( 'ZMass_mu_extraMu' , 'ZMass_mu_extraMu' , 120, 60, 120)  #ZMass , Z->mumu + Extra mu
        ZMass_mu_hist_extraEl  = TH1F( 'ZMass_mu_extraEl' , 'ZMass_mu_extraEl' , 120, 60, 120)  #ZMass , Z->mumu + Extra e
    


    tree.SetBranchStatus("*",0)  # disable all branches
    if ZTree :
        tree.SetBranchStatus("Zsel",1)
    else : 
        tree.SetBranchStatus("ZZsel",1)
    tree.SetBranchStatus("LepLepId",1)
    if ZTree :
        tree.SetBranchStatus("ZMass",1)
    else :
        tree.SetBranchStatus("Z1Mass",1)
    tree.SetBranchStatus("LepEta",1)



    #read ttree ZTree
    if ZTree :
        print "reading tree", data.GetName(),treeText,tree.GetName()  ,"..."
        for event in tree:
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger 
            
            #Z->ee histos
            if(int(math.fabs(event.LepLepId[0])) == 11 ):
                ZMass_ele_hist.Fill(event.ZMass)             #ZMass , Z->ee
                if(math.fabs(event.LepEta[0]) < 1.479 and math.fabs(event.LepEta[1]) < 1.479):
                    ZMass_ele_hist_EBEB.Fill(event.ZMass)    #ZMass , Z->ee , Barrel-Barrel
                elif(math.fabs(event.LepEta[0]) > 1.479 and math.fabs(event.LepEta[1]) > 1.479):
                    ZMass_ele_hist_EEEE.Fill(event.ZMass)    #ZMass , Z->ee , Endcap-Endcap
                else:
                    ZMass_ele_hist_EBEE.Fill(event.ZMass)    #ZMass , Z->ee , Barrel-Endcap
            
            #Z->mumu histos
            elif(int(math.fabs(event.LepLepId[0])) == 13 ):
                ZMass_mu_hist.Fill(event.ZMass)              #ZMass , Z->mumu
                if(math.fabs(event.LepEta[0]) < 1. and math.fabs(event.LepEta[1]) < 1.):
                    ZMass_mu_hist_MBMB.Fill(event.ZMass)     #ZMass , Z->mumu , Barrel-Barrel
                elif(math.fabs(event.LepEta[0]) > 1. and math.fabs(event.LepEta[1]) > 1.):
                    ZMass_mu_hist_MEME.Fill(event.ZMass)     #ZMass , Z->mumu , Endcap-Endcap
                else:
                    ZMass_mu_hist_MBME.Fill(event.ZMass)     #ZMass , Z->mumu , Barrel-Endcap
        
    
    #read ttree ZZTree or CRZLTree
    else: 
        print "reading tree", data.GetName(),treeText,tree.GetName()  ,"..."
        for event in tree:
            if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger 
        
            #Z->ee histos
            if(int(math.fabs(event.LepLepId[0])) == 11 ):
                ZMass_ele_hist.Fill(event.Z1Mass)             #ZMass , Z->ee
                if(math.fabs(event.LepEta[0]) < 1.479 and math.fabs(event.LepEta[1]) < 1.479):
                    ZMass_ele_hist_EBEB.Fill(event.Z1Mass)    #ZMass , Z->ee , Barrel-Barrel
                elif(math.fabs(event.LepEta[0]) > 1.479 and math.fabs(event.LepEta[1]) > 1.479):
                    ZMass_ele_hist_EEEE.Fill(event.Z1Mass)    #ZMass , Z->ee , Endcap-Endcap
                else:
                    ZMass_ele_hist_EBEE.Fill(event.Z1Mass)    #ZMass , Z->ee , Barrel-Endcap
                if(int(math.fabs(event.LepLepId[2])) == 11 ):
                    ZMass_ele_hist_extraEl.Fill(event.Z1Mass) #ZMass , Z->ee   + Extra e
                else:
                    ZMass_ele_hist_extraMu.Fill(event.Z1Mass) #ZMass , Z->ee   + Extra mu
       
            #Z->mumu histos
            elif(int(math.fabs(event.LepLepId[0])) == 13 ):
                ZMass_mu_hist.Fill(event.Z1Mass)              #ZMass , Z->mumu
                if(math.fabs(event.LepEta[0]) < 1. and math.fabs(event.LepEta[1]) < 1.):
                    ZMass_mu_hist_MBMB.Fill(event.Z1Mass)     #ZMass , Z->mumu , Barrel-Barrel
                elif(math.fabs(event.LepEta[0]) > 1. and math.fabs(event.LepEta[1]) > 1.):
                    ZMass_mu_hist_MEME.Fill(event.Z1Mass)     #ZMass , Z->mumu , Endcap-Endcap
                else:
                    ZMass_mu_hist_MBME.Fill(event.Z1Mass)     #ZMass , Z->mumu , Barrel-Endcap
                if(int(math.fabs(event.LepLepId[2])) == 11 ):
                    ZMass_mu_hist_extraEl.Fill(event.Z1Mass)  #ZMass , Z->mumu + Extra e
                else:
                    ZMass_mu_hist_extraMu.Fill(event.Z1Mass)  #ZMass , Z->mumu + Extra mu


    #save histograms in a root file 
    print "saving histograms into root file ..."
    histoDATA = TFile.Open("histoDATA_"+ period + "_" + treeText +".root", "RECREATE")
    histoDATA.cd()
    ZMass_ele_hist.Write()  
    if not ZTree :
        ZMass_ele_hist_extraMu.Write() 
        ZMass_ele_hist_extraEl.Write() 
    ZMass_ele_hist_EBEB.Write()    
    ZMass_ele_hist_EBEE.Write()    
    ZMass_ele_hist_EEEE.Write()       
    ZMass_mu_hist.Write()      
    if not ZTree :    
        ZMass_mu_hist_extraMu.Write()  
        ZMass_mu_hist_extraEl.Write()  
    ZMass_mu_hist_MBMB.Write()     
    ZMass_mu_hist_MBME.Write()     
    ZMass_mu_hist_MEME.Write()     

    histoDATA.Close()

#read histo from histoDATA.root
histoDATA_input = TFile.Open("histoDATA_"+ period + "_" + treeText +".root")
print 'Reading file', histoDATA_input.GetName(),'...'

histogramList = []

histogramList.append(histoDATA_input.Get('ZMass_ele'))
if not ZTree :
    histogramList.append(histoDATA_input.Get('ZMass_ele_extraMu'))
    histogramList.append(histoDATA_input.Get('ZMass_ele_extraEl'))
histogramList.append(histoDATA_input.Get('ZMass_ele_EBEB'))
histogramList.append(histoDATA_input.Get('ZMass_ele_EBEE'))
histogramList.append(histoDATA_input.Get('ZMass_ele_EEEE'))
histogramList.append(histoDATA_input.Get('ZMass_mu'))
if not ZTree :
    histogramList.append(histoDATA_input.Get('ZMass_mu_extraMu'))
    histogramList.append(histoDATA_input.Get('ZMass_mu_extraEl'))
histogramList.append(histoDATA_input.Get('ZMass_mu_MBMB'))
histogramList.append(histoDATA_input.Get('ZMass_mu_MBME'))
histogramList.append(histoDATA_input.Get('ZMass_mu_MEME'))



#define lists for the fit function
if ZTree :
    histTitleList  = ['ZMass_ele_hist', 'ZMass_ele_hist_EBEB', 'ZMass_ele_hist_EBEE', 'ZMass_ele_hist_EEEE', 'ZMass_mu_hist', 'ZMass_mu_hist_MBMB', 'ZMass_mu_hist_MBME', 'ZMass_mu_hist_MEME']
else :
    histTitleList  = ['ZMass_ele_hist', 'ZMass_ele_hist_extraMu', 'ZMass_ele_hist_extraEl', 'ZMass_ele_hist_EBEB', 'ZMass_ele_hist_EBEE', 'ZMass_ele_hist_EEEE', 'ZMass_mu_hist', 'ZMass_mu_hist_extraMu', 'ZMass_mu_hist_extraEl', 'ZMass_mu_hist_MBMB', 'ZMass_mu_hist_MBME', 'ZMass_mu_hist_MEME']

luminosityList = [ lumi for i in range(len(histogramList))]

#do the fit 
fitResult = DoSimpleFit(histogramList, luminosityList, ZZTree, outputDir, histTitleList, fitMC, fitDATA)
#fitResult = DoDCBunbinnedFit(histogramList, luminosityList, ZZTree, outputDir, histTitleList, fitMC, fitDATA) #to be fixed!

print "Fit done!!"



# ************************************
# store fit results in dictionaries
if CRZLTree :

    massFitDATA_dict  = {'Zee': fitResult[0].mean, 'Zee_extraMu': fitResult[1].mean, 'Zee_extraEl': fitResult[2].mean, 'Zee_EBEB': fitResult[3].mean, 'Zee_EBEE': fitResult[4].mean, 'Zee_EEEE': fitResult[5].mean, 'Zmumu': fitResult[6].mean, 'Zmumu_extraMu': fitResult[7].mean, 'Zmumu_extraEl': fitResult[8].mean, 'Zmumu_MBMB': fitResult[9].mean, 'Zmumu_MBME': fitResult[10].mean, 'Zmumu_MEME': fitResult[11].mean}

    widthFitDATA_dict = {'Zee': fitResult[0].width, 'Zee_extraMu': fitResult[1].width, 'Zee_extraEl': fitResult[2].width, 'Zee_EBEB': fitResult[3].width, 'Zee_EBEE': fitResult[4].width, 'Zee_EEEE': fitResult[5].width, 'Zmumu': fitResult[6].width, 'Zmumu_extraMu': fitResult[7].width, 'Zmumu_extraEl': fitResult[8].width, 'Zmumu_MBMB': fitResult[9].width, 'Zmumu_MBME': fitResult[10].width, 'Zmumu_MEME': fitResult[11].width}


elif ZTree :

    massFitDATA_dict  = {'Zee': fitResult[0].mean, 'Zee_extraMu': 0, 'Zee_extraEl': 0, 'Zee_EBEB': fitResult[1].mean, 'Zee_EBEE': fitResult[2].mean, 'Zee_EEEE': fitResult[3].mean, 'Zmumu': fitResult[4].mean, 'Zmumu_extraMu': 0, 'Zmumu_extraEl': 0, 'Zmumu_MBMB': fitResult[5].mean, 'Zmumu_MBME': fitResult[6].mean, 'Zmumu_MEME': fitResult[7].mean }

    widthFitDATA_dict = {'Zee': fitResult[0].width, 'Zee_extraMu': 0, 'Zee_extraEl': 0, 'Zee_EBEB': fitResult[1].width, 'Zee_EBEE': fitResult[2].width, 'Zee_EEEE': fitResult[3].width, 'Zmumu': fitResult[4].width, 'Zmumu_extraMu': 0, 'Zmumu_extraEl': 0, 'Zmumu_MBMB': fitResult[5].width, 'Zmumu_MBME': fitResult[6].width, 'Zmumu_MEME': fitResult[7].width}




# ***********************************
# save result on output files
# save Z DCB mean result on file json

with open("out_ZDBCmean_DATA_" + period + "_" + treeText  + ".json","w") as handle1 :
    json.dump(massFitDATA_dict, handle1)

print "Z DCB mean output written on file"


# save Z DCB width result on file json
with open("out_ZDBCwidth_DATA_" + period + "_" + treeText  + ".json","w") as handle2 :
    json.dump(widthFitDATA_dict, handle2)

print "Z DCB width output written on file"
#************************************


