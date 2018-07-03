#!/usr/bin/env python


# *******************
# usage: 
#    - choose the options 
#    - python FitMC.py
#
# structure:
#    - read MC file (inputTree) (1 file by time)
#    - store weighted events in histos and save them in a file.root (if redoHistos)
#    - read the root file 
#    - fit the histos with a function from helper.py and print on file fit results
# ********************
 
import json
import ROOT, helper, math
from ROOT import TFile, TH1, TH1F, TCanvas, gSystem
from helper import DoSimpleFit, Result

# *****************************
# Declare all the variables
# fit options 
fitMC       = True    #true for fitting MC in FitMC.py 
fitDATA     = False   #true for fitting DATA in FitDATA.py 
redoHistos  = True
applyPU2018 = True    # true for prompt-reco 2018 data 

# data tree options 
ZZTree   = False
CRZLTree = False
ZTree    = True

# data periods options
period = "data2018"

# MC sample
MCsample = "DYJets"
# MCsample = "TTJets"
# *****************************


if(period == "data2018"):
    lumi = 16.6     # fb-1

#input file
if(MCsample == "DYJets"):
    inputTree = TFile.Open("/data3/Higgs/180416/MC_main/DYJetsToLL_M50/ZZ4lAnalysis.root") #DYJets
elif(MCsample == "TTJets"):
    inputTree = TFile.Open("/data3/Higgs/180416/MC_main/TTTo2L2Nu/ZZ4lAnalysis.root")  #TTJets 
else:
    print ("Error: wrong MC sample!")


#PU weights
if(applyPU2018):
    fPU = TFile.Open("../../data/PileUpWeights/puWeights_Fall2017MC_to_2018Data_314472-317591.root");
    hPUWeight= fPU.Get("weights")


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


#create output directory 
outputDir = "FitResults_MC_" + str(MCsample) + "_" + str(period) + "_" + str(treeText) 
gSystem.Exec("mkdir -p " + outputDir)
print "Output directories created!"




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
    


    # get partial event weight
    hcounters           = inputTree.Get("ZZTree/Counters")
    gen_sumWeights      = hcounters.GetBinContent(40)
    partialSampleWeight = lumi * 1000 / gen_sumWeights


    #read ttree ZTree
    if ZTree : 
        print "reading tree", inputTree.GetName(),treeText,tree.GetName()  ,"..."
        for event in tree:
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger 

            if(applyPU2018):
                if event.PUWeight == 0. : 
                    PUWeight_old = 0.00000001
                else :
                    PUWeight_old = event.PUWeight
                weight = partialSampleWeight*event.xsec*event.overallEventWeight/PUWeight_old*(hPUWeight.GetBinContent(hPUWeight.FindBin(event.NTrueInt)));
            else:
                weight = partialSampleWeight*event.xsec*event.overallEventWeight

            #Z->ee histos
            if(int(math.fabs(event.LepLepId[0])) == 11 ):
                ZMass_ele_hist.Fill(event.ZMass, weight)             #ZMass , Z->ee
                if(math.fabs(event.LepEta[0]) < 1.479 and math.fabs(event.LepEta[1]) < 1.479):
                    ZMass_ele_hist_EBEB.Fill(event.ZMass, weight)    #ZMass , Z->ee , Barrel-Barrel
                elif(math.fabs(event.LepEta[0]) > 1.479 and math.fabs(event.LepEta[1]) > 1.479):
                    ZMass_ele_hist_EEEE.Fill(event.ZMass, weight)    #ZMass , Z->ee , Endcap-Endcap
                else:
                    ZMass_ele_hist_EBEE.Fill(event.ZMass, weight)    #ZMass , Z->ee , Barrel-Endcap

            #Z->mumu histos
            elif(int(math.fabs(event.LepLepId[0])) == 13 ):
                ZMass_mu_hist.Fill(event.ZMass, weight)              #ZMass , Z->mumu
                if(math.fabs(event.LepEta[0]) < 1. and math.fabs(event.LepEta[1]) < 1.):
                    ZMass_mu_hist_MBMB.Fill(event.ZMass, weight)     #ZMass , Z->mumu , Barrel-Barrel
                elif(math.fabs(event.LepEta[0]) > 1. and math.fabs(event.LepEta[1]) > 1.):
                    ZMass_mu_hist_MEME.Fill(event.ZMass, weight)     #ZMass , Z->mumu , Endcap-Endcap
                else:
                    ZMass_mu_hist_MBME.Fill(event.ZMass, weight)     #ZMass , Z->mumu , Barrel-Endcap


    #read ttree ZZTree or CRZLTree
    else :    
        print "reading tree", inputTree.GetName(),treeText,tree.GetName()  ,"..."
        for event in tree:
            if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger 
        
            if(applyPU2018):
                weight = partialSampleWeight*event.xsec*event.overallEventWeight/event.PUWeight*(hPUWeight.GetBinContent(hPUWeight.FindBin(event.NTrueInt)));
            else:
                weight = partialSampleWeight*event.xsec*event.overallEventWeight
            
    
            #Z->ee histos
            if(int(math.fabs(event.LepLepId[0])) == 11 ):
                ZMass_ele_hist.Fill(event.Z1Mass, weight)             #ZMass , Z->ee
                if(math.fabs(event.LepEta[0]) < 1.479 and math.fabs(event.LepEta[1]) < 1.479):
                    ZMass_ele_hist_EBEB.Fill(event.Z1Mass, weight)    #ZMass , Z->ee , Barrel-Barrel
                elif(math.fabs(event.LepEta[0]) > 1.479 and math.fabs(event.LepEta[1]) > 1.479):
                    ZMass_ele_hist_EEEE.Fill(event.Z1Mass, weight)    #ZMass , Z->ee , Endcap-Endcap
                else:
                    ZMass_ele_hist_EBEE.Fill(event.Z1Mass, weight)    #ZMass , Z->ee , Barrel-Endcap
                if(int(math.fabs(event.LepLepId[2])) == 11 ):
                    ZMass_ele_hist_extraEl.Fill(event.Z1Mass, weight) #ZMass , Z->ee   + Extra e
                else:
                    ZMass_ele_hist_extraMu.Fill(event.Z1Mass, weight) #ZMass , Z->ee   + Extra mu
    
            #Z->mumu histos
            elif(int(math.fabs(event.LepLepId[0])) == 13 ):
                ZMass_mu_hist.Fill(event.Z1Mass, weight)              #ZMass , Z->mumu
                if(math.fabs(event.LepEta[0]) < 1. and math.fabs(event.LepEta[1]) < 1.):
                    ZMass_mu_hist_MBMB.Fill(event.Z1Mass, weight)     #ZMass , Z->mumu , Barrel-Barrel
                elif(math.fabs(event.LepEta[0]) > 1. and math.fabs(event.LepEta[1]) > 1.):
                    ZMass_mu_hist_MEME.Fill(event.Z1Mass, weight)     #ZMass , Z->mumu , Endcap-Endcap
                else:
                    ZMass_mu_hist_MBME.Fill(event.Z1Mass, weight)     #ZMass , Z->mumu , Barrel-Endcap
                if(int(math.fabs(event.LepLepId[2])) == 11 ):
                    ZMass_mu_hist_extraEl.Fill(event.Z1Mass, weight)  #ZMass , Z->mumu + Extra e
                else:
                    ZMass_mu_hist_extraMu.Fill(event.Z1Mass, weight)  #ZMass , Z->mumu + Extra mu

    #save histograms in a root file 
    print "saving histograms to histo MC ..."
    histoMC = TFile.Open("histoMC_" + MCsample + "_" + period + "_" + treeText + ".root", "RECREATE")
    histoMC.cd()
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

    histoMC.Close()


#read histo from histoMC.root
histoMC_input = TFile.Open("histoMC_" + MCsample + "_" + period + "_" + treeText + ".root")
print 'Reading file', histoMC_input.GetName(),'...'

histogramList = []

histogramList.append(histoMC_input.Get('ZMass_ele'))
if not ZTree :
    histogramList.append(histoMC_input.Get('ZMass_ele_extraMu'))
    histogramList.append(histoMC_input.Get('ZMass_ele_extraEl'))
histogramList.append(histoMC_input.Get('ZMass_ele_EBEB'))
histogramList.append(histoMC_input.Get('ZMass_ele_EBEE'))
histogramList.append(histoMC_input.Get('ZMass_ele_EEEE'))
histogramList.append(histoMC_input.Get('ZMass_mu'))
if not ZTree :
    histogramList.append(histoMC_input.Get('ZMass_mu_extraMu'))
    histogramList.append(histoMC_input.Get('ZMass_mu_extraEl'))
histogramList.append(histoMC_input.Get('ZMass_mu_MBMB'))
histogramList.append(histoMC_input.Get('ZMass_mu_MBME'))
histogramList.append(histoMC_input.Get('ZMass_mu_MEME'))



#define lists for the fit function
if ZTree :
    histTitleList  = ['ZMass_ele_hist', 'ZMass_ele_hist_EBEB', 'ZMass_ele_hist_EBEE', 'ZMass_ele_hist_EEEE', 'ZMass_mu_hist', 'ZMass_mu_hist_MBMB', 'ZMass_mu_hist_MBME', 'ZMass_mu_hist_MEME']
else : 
    histTitleList  = ['ZMass_ele_hist', 'ZMass_ele_hist_extraMu', 'ZMass_ele_hist_extraEl', 'ZMass_ele_hist_EBEB', 'ZMass_ele_hist_EBEE', 'ZMass_ele_hist_EEEE', 'ZMass_mu_hist', 'ZMass_mu_hist_extraMu', 'ZMass_mu_hist_extraEl', 'ZMass_mu_hist_MBMB', 'ZMass_mu_hist_MBME', 'ZMass_mu_hist_MEME']

luminosityList = [ lumi for i in range(len(histogramList))]


#do the fit
if(MCsample == "DYJets"):
    fitResult = DoSimpleFit(histogramList, luminosityList, ZZTree, outputDir, histTitleList, fitMC, fitDATA)

    print "Fit done!!"


    # ************************************
    # store fit results in dictionaries
    if CRZLTree :
        
        massFitMC_dict  = {'Zee': fitResult[0].mean, 'Zee_extraMu': fitResult[1].mean, 'Zee_extraEl': fitResult[2].mean, 'Zee_EBEB': fitResult[3].mean, 'Zee_EBEE': fitResult[4].mean, 'Zee_EEEE': fitResult[5].mean, 'Zmumu': fitResult[6].mean, 'Zmumu_extraMu': fitResult[7].mean, 'Zmumu_extraEl': fitResult[8].mean, 'Zmumu_MBMB': fitResult[9].mean, 'Zmumu_MBME': fitResult[10].mean, 'Zmumu_MEME': fitResult[11].mean}

        widthFitMC_dict = {'Zee': fitResult[0].width, 'Zee_extraMu': fitResult[1].width, 'Zee_extraEl': fitResult[2].width, 'Zee_EBEB': fitResult[3].width, 'Zee_EBEE': fitResult[4].width, 'Zee_EEEE': fitResult[5].width, 'Zmumu': fitResult[6].width, 'Zmumu_extraMu': fitResult[7].width, 'Zmumu_extraEl': fitResult[8].width, 'Zmumu_MBMB': fitResult[9].width, 'Zmumu_MBME': fitResult[10].width, 'Zmumu_MEME': fitResult[11].width}


    elif ZTree :

        massFitMC_dict  = {'Zee': fitResult[0].mean, 'Zee_extraMu': 0, 'Zee_extraEl': 0, 'Zee_EBEB': fitResult[1].mean, 'Zee_EBEE': fitResult[2].mean, 'Zee_EEEE': fitResult[3].mean, 'Zmumu': fitResult[4].mean, 'Zmumu_extraMu': 0, 'Zmumu_extraEl': 0, 'Zmumu_MBMB': fitResult[5].mean, 'Zmumu_MBME': fitResult[6].mean, 'Zmumu_MEME': fitResult[7].mean }

        widthFitMC_dict = {'Zee': fitResult[0].width, 'Zee_extraMu': 0, 'Zee_extraEl': 0, 'Zee_EBEB': fitResult[1].width, 'Zee_EBEE': fitResult[2].width, 'Zee_EEEE': fitResult[3].width, 'Zmumu': fitResult[4].width, 'Zmumu_extraMu': 0, 'Zmumu_extraEl': 0, 'Zmumu_MBMB': fitResult[5].width, 'Zmumu_MBME': fitResult[6].width, 'Zmumu_MEME': fitResult[7].width}


    # ***********************************
    # save result on output files
    # save Z DCB mean result on file json
    with open("out_ZDBCmean_MC_" + MCsample + "_" + period + "_" + treeText + ".json","w") as handle1 :
        json.dump(massFitMC_dict, handle1)

    print "Z DCB mean output written on file"


    # save Z DCB width result on file json
    with open("out_ZDBCwidth_MC_" + MCsample + "_" + period + "_" + treeText  + ".json","w") as handle2 :
        json.dump(widthFitMC_dict, handle2)

    print "Z DCB width output written on file"
    #************************************
    



    
    
    
    
    
    
