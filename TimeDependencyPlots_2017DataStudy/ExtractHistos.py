#!/usr/bin/env python


# *******************
# usage: 
#    python ExtractHistos.py
#
# structure:
#    - read data file and json splitted in lumiblocks (produced by JSON_calc/SplitPerLumi.py)
#    - Fill histos with data for each lumi block
#    - save histo in a root file
# ********************

# First import ROOT and helper functions
import ROOT, helper, math
ROOT.gROOT.SetBatch(True)
from helper import ReadJSON
from ROOT import TFile, TH1F

# data tree options 
ZZTree   = False
CRZLTree = False
ZTree    = True

# data periods options
# period = "data2016"
period = "data2017"



if(period == "data2016"):
    data     = TFile.Open("/data3/Higgs/170222/AllData/ZZ4lAnalysis.root")
    inputTXT = "JSON_calc/SplittedBlocks_2016data_0p5_new.txt"
    saveAs   = "2016data_0p5_reduced"

elif(period == "data2017"):
    data     = TFile.Open("/data3/Higgs/180305/AllData/ZZ4lAnalysis.root")
    inputTXT = "JSON_calc/SplittedBlocks_2017data_0p5_new.txt"
    saveAs   = "2017data_0p5_reduced"

else:
    print ("Error, wrong period chosen!!!")



# choose the data tree 
if(ZZTree):
    tree           = data.Get("ZZTree/candTree")
    outputFileRoot = saveAs + "_ZZTree_histos.root"
elif(CRZLTree):
    tree           = data.Get("CRZLTree/candTree")
    outputFileRoot = saveAs + "_CRZLTree_histos.root"
elif(ZTree):
    tree           = data.Get("ZTree/candTree")
    outputFileRoot = saveAs + "_ZTree_histos.root"
else: 
    print ("Error: wrong option!")



RunNum_B = []
LumiNum_B = []
RunNum_E = []
LumiNum_E = []
recorded = []

# Read the input JSON file and store it
ReadJSON(inputTXT, RunNum_B, LumiNum_B, RunNum_E, LumiNum_E, recorded)
# Z->ee histos
ZMass_ele_hist = []
ZMass_ele_hist_EBEB = []
ZMass_ele_hist_EBEE = []
ZMass_ele_hist_EEEE = []
if not ZTree :                   # extra lepton histos possible if number of leptons > 2
    ZMass_ele_hist_extraMu = []
    ZMass_ele_hist_extraEl = []

# Z->mumu histos
ZMass_mu_hist = []
ZMass_mu_hist_MBMB = []
ZMass_mu_hist_MBME = []
ZMass_mu_hist_MEME = []
if not ZTree :
    ZMass_mu_hist_extraMu = []
    ZMass_mu_hist_extraEl = []

# lepton histos
electronISO_Max_hist = []
electronISO_Min_hist = []
muonISO_Max_hist = []
muonISO_Min_hist = []
electronSIP_Max_hist = []
electronSIP_Min_hist = []
muonSIP_Max_hist = []
muonSIP_Min_hist = []



for i in range(0,len(RunNum_B)):
    
    histo2  = TH1F( 'ZMass_ele'         + str(i), 'ZMass_ele'         + str(i), 120, 60, 120)
    histo21 = TH1F( 'ZMass_ele_EBEB_'   + str(i), 'ZMass_ele_EBEB_'   + str(i), 120, 60, 120)
    histo22 = TH1F( 'ZMass_ele_EBEE_'   + str(i), 'ZMass_ele_EBEE_'   + str(i), 120, 60, 120)
    histo23 = TH1F( 'ZMass_ele_EEEE_'   + str(i), 'ZMass_ele_EEEE_'   + str(i), 120, 60, 120)
    if not ZTree :
        histo13 = TH1F( 'ZMass_ele_extraMu' + str(i), 'ZMass_ele_extraMu' + str(i), 120, 60, 120)
        histo16 = TH1F( 'ZMass_ele_extraEl' + str(i), 'ZMass_ele_extraEl' + str(i), 120, 60, 120)
   
    histo3  = TH1F( 'ZMass_mu'         + str(i), 'ZMass_mu'         + str(i), 120, 60, 120)
    histo24 = TH1F( 'ZMass_mu_MBMB_'   + str(i), 'ZMass_mu_MBMB_'   + str(i), 120, 60, 120)
    histo25 = TH1F( 'ZMass_mu_MBME_'   + str(i), 'ZMass_mu_MBME_'   + str(i), 120, 60, 120)
    histo26 = TH1F( 'ZMass_mu_MEME_'   + str(i), 'ZMass_mu_MEME_'   + str(i), 120, 60, 120)
    if not ZTree :
        histo14 = TH1F( 'ZMass_mu_extraMu' + str(i), 'ZMass_mu_extraMu' + str(i), 120, 60, 120)
        histo17 = TH1F( 'ZMass_mu_extraEl' + str(i), 'ZMass_mu_extraEl' + str(i), 120, 60, 120) 
    
    histo4  = TH1F( 'ElectronISO_Max' + str(i), 'ElectronISO_Max' + str(i), 100, 0, 0.5)
    histo5  = TH1F( 'ElectronISO_Min' + str(i), 'ElectronISO_Min' + str(i), 100, 0, 0.5)
    histo6  = TH1F( 'MuonISO_Max'     + str(i), 'MuonISO_Max'     + str(i), 100, 0, 0.5)
    histo7  = TH1F( 'MuonISO_Min'     + str(i), 'MuonISO_Min'     + str(i), 100, 0, 0.5)
    histo8  = TH1F( 'ElectronSIP_Max' + str(i), 'ElectronSIP_Max' + str(i), 100, 0, 5.0)
    histo9  = TH1F( 'ElectronSIP_Min' + str(i), 'ElectronSIP_Min' + str(i), 100, 0, 5.0)
    histo10 = TH1F( 'MuonSIP_Max'     + str(i), 'MuonSIP_Max'     + str(i), 100, 0, 5.0)
    histo11 = TH1F( 'MuonSIP_Min'     + str(i), 'MuonSIP_Min'     + str(i), 100, 0, 5.0)
    
       
    ZMass_ele_hist.append(histo2)
    ZMass_ele_hist_EBEB.append(histo21)
    ZMass_ele_hist_EBEE.append(histo22)
    ZMass_ele_hist_EEEE.append(histo23)
    if not ZTree :
        ZMass_ele_hist_extraMu.append(histo13)
        ZMass_ele_hist_extraEl.append(histo16)
    
    ZMass_mu_hist.append(histo3)
    ZMass_mu_hist_MBMB.append(histo24)
    ZMass_mu_hist_MBME.append(histo25)
    ZMass_mu_hist_MEME.append(histo26)
    if not ZTree :
        ZMass_mu_hist_extraMu.append(histo14)
        ZMass_mu_hist_extraEl.append(histo17)
    
    electronISO_Max_hist.append(histo4)
    electronISO_Min_hist.append(histo5)
    muonISO_Max_hist.append(histo6)
    muonISO_Min_hist.append(histo7)
    electronSIP_Max_hist.append(histo8)
    electronSIP_Min_hist.append(histo9)
    muonSIP_Max_hist.append(histo10)
    muonSIP_Min_hist.append(histo11)
    

print "Splitting data into lumi blocks."


tree.SetBranchStatus("*",0)  # disable all branches
if ZTree :
    tree.SetBranchStatus("Zsel",1)
else : 
    tree.SetBranchStatus("ZZsel",1)
tree.SetBranchStatus("RunNumber",1)
tree.SetBranchStatus("LepLepId",1)
tree.SetBranchStatus("ZMass",1)
tree.SetBranchStatus("LepEta",1)
tree.SetBranchStatus("LepCombRelIsoPF",1)
tree.SetBranchStatus("LepSIP",1)



if ZTree : 
    for i in range(0, len(RunNum_B)):
        print(str(i+1)+". lumi block...")
        for event in tree:
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger 
            
            if ( int(event.RunNumber) >= int(RunNum_B[i]) and int(event.RunNumber) <= int(RunNum_E[i])):
                
                #Z->ee histos
                if(int(math.fabs(event.LepLepId[0])) == 11 ):
                    ZMass_ele_hist[i].Fill(event.ZMass)
                    if(math.fabs(event.LepEta[0]) < 1.479 and math.fabs(event.LepEta[1]) < 1.479):
                        ZMass_ele_hist_EBEB[i].Fill(event.ZMass)
                    elif(math.fabs(event.LepEta[0]) > 1.479 and math.fabs(event.LepEta[1]) > 1.479):
                        ZMass_ele_hist_EEEE[i].Fill(event.ZMass)
                    else:
                        ZMass_ele_hist_EBEE[i].Fill(event.ZMass)
                    
                    if(event.LepCombRelIsoPF[0] > event.LepCombRelIsoPF[1]):
                        electronISO_Max_hist[i].Fill(event.LepCombRelIsoPF[0])
                        electronISO_Min_hist[i].Fill(event.LepCombRelIsoPF[1])
                    else:
                        electronISO_Max_hist[i].Fill(event.LepCombRelIsoPF[1])
                        electronISO_Min_hist[i].Fill(event.LepCombRelIsoPF[0])
                    
                    if(event.LepSIP[0] > event.LepSIP[1]):
                        electronSIP_Max_hist[i].Fill(event.LepSIP[0])
                        electronSIP_Min_hist[i].Fill(event.LepSIP[1])
                    else:
                        electronSIP_Max_hist[i].Fill(event.LepSIP[1])
                        electronSIP_Min_hist[i].Fill(event.LepSIP[0])

                #Z->mumu histos
                elif(int(math.fabs(event.LepLepId[0])) == 13 ):
                    ZMass_mu_hist[i].Fill(event.ZMass)
                    if(math.fabs(event.LepEta[0]) < 1. and math.fabs(event.LepEta[1]) < 1.):
                        ZMass_mu_hist_MBMB[i].Fill(event.ZMass)
                    elif(math.fabs(event.LepEta[0]) > 1. and math.fabs(event.LepEta[1]) > 1.):
                        ZMass_mu_hist_MEME[i].Fill(event.ZMass)
                    else:
                        ZMass_mu_hist_MBME[i].Fill(event.ZMass)

                    if(event.LepCombRelIsoPF[0] > event.LepCombRelIsoPF[1]):
                        muonISO_Max_hist[i].Fill(event.LepCombRelIsoPF[0])
                        muonISO_Min_hist[i].Fill(event.LepCombRelIsoPF[1])
                    else:
                        muonISO_Max_hist[i].Fill(event.LepCombRelIsoPF[1])
                        muonISO_Min_hist[i].Fill(event.LepCombRelIsoPF[0])
                
                    if(event.LepSIP[0] > event.LepSIP[1]):
                        muonSIP_Max_hist[i].Fill(event.LepSIP[0])
                        muonSIP_Min_hist[i].Fill(event.LepSIP[1])
                    else:
                        muonSIP_Max_hist[i].Fill(event.LepSIP[1])
                        muonSIP_Min_hist[i].Fill(event.LepSIP[0])


else : 

    for i in range(0, len(RunNum_B)):
        print(str(i+1)+". lumi block...")
        for event in tree:
            if ( event.ZZsel < 0 ) : continue # skip events that do not pass the trigger 

            if ( int(event.RunNumber) >= int(RunNum_B[i]) and int(event.RunNumber) <= int(RunNum_E[i])):

                #Z->ee histos
                if(int(math.fabs(event.LepLepId[0])) == 11 ):
                    ZMass_ele_hist[i].Fill(event.Z1Mass)
                    if(math.fabs(event.LepEta[0]) < 1.479 and math.fabs(event.LepEta[1]) < 1.479):
                        ZMass_ele_hist_EBEB[i].Fill(event.Z1Mass)
                    elif(math.fabs(event.LepEta[0]) > 1.479 and math.fabs(event.LepEta[1]) > 1.479):
                        ZMass_ele_hist_EEEE[i].Fill(event.Z1Mass)
                    else:
                        ZMass_ele_hist_EBEE[i].Fill(event.Z1Mass)
                    if(int(math.fabs(event.LepLepId[2])) == 11 ):
                        ZMass_ele_hist_extraEl[i].Fill(event.Z1Mass)
                    else:
                        ZMass_ele_hist_extraMu[i].Fill(event.Z1Mass)

                    if(event.LepCombRelIsoPF[0] > event.LepCombRelIsoPF[1]):
                        electronISO_Max_hist[i].Fill(event.LepCombRelIsoPF[0])
                        electronISO_Min_hist[i].Fill(event.LepCombRelIsoPF[1])
                    else:
                        electronISO_Max_hist[i].Fill(event.LepCombRelIsoPF[1])
                        electronISO_Min_hist[i].Fill(event.LepCombRelIsoPF[0])
                
                    if(event.LepSIP[0] > event.LepSIP[1]):
                        electronSIP_Max_hist[i].Fill(event.LepSIP[0])
                        electronSIP_Min_hist[i].Fill(event.LepSIP[1])
                    else:
                        electronSIP_Max_hist[i].Fill(event.LepSIP[1])
                        electronSIP_Min_hist[i].Fill(event.LepSIP[0])

                #Z->mumu histos
                elif(int(math.fabs(event.LepLepId[0])) == 13 ):
                    ZMass_mu_hist[i].Fill(event.Z1Mass)
                    if(math.fabs(event.LepEta[0]) < 1. and math.fabs(event.LepEta[1]) < 1.):
                        ZMass_mu_hist_MBMB[i].Fill(event.Z1Mass)
                    elif(math.fabs(event.LepEta[0]) > 1. and math.fabs(event.LepEta[1]) > 1.):
                        ZMass_mu_hist_MEME[i].Fill(event.Z1Mass)
                    else:
                        ZMass_mu_hist_MBME[i].Fill(event.Z1Mass)
                    if(int(math.fabs(event.LepLepId[2])) == 11 ):
                        ZMass_mu_hist_extraEl[i].Fill(event.Z1Mass)
                    else:
                        ZMass_mu_hist_extraMu[i].Fill(event.Z1Mass)
                    
                    if(event.LepCombRelIsoPF[0] > event.LepCombRelIsoPF[1]):
                        muonISO_Max_hist[i].Fill(event.LepCombRelIsoPF[0])
                        muonISO_Min_hist[i].Fill(event.LepCombRelIsoPF[1])
                    else:
                        muonISO_Max_hist[i].Fill(event.LepCombRelIsoPF[1])
                        muonISO_Min_hist[i].Fill(event.LepCombRelIsoPF[0])
                    
                    if(event.LepSIP[0] > event.LepSIP[1]):
                        muonSIP_Max_hist[i].Fill(event.LepSIP[0])
                        muonSIP_Min_hist[i].Fill(event.LepSIP[1])
                    else:
                        muonSIP_Max_hist[i].Fill(event.LepSIP[1])
                        muonSIP_Min_hist[i].Fill(event.LepSIP[0])



print "Writing output file ..."
output = TFile.Open( outputFileRoot, "RECREATE")
for i in range(0, len(RunNum_B)):
    ZMass_ele_hist[i].Write()
    ZMass_ele_hist_EBEB[i].Write()
    ZMass_ele_hist_EBEE[i].Write()
    ZMass_ele_hist_EEEE[i].Write()
    if not ZTree :
        ZMass_ele_hist_extraMu[i].Write()
        ZMass_ele_hist_extraEl[i].Write()
    
    ZMass_mu_hist[i].Write()
    ZMass_mu_hist_MBMB[i].Write()
    ZMass_mu_hist_MBME[i].Write()
    ZMass_mu_hist_MEME[i].Write()
    if not ZTree :
        ZMass_mu_hist_extraMu[i].Write()    
        ZMass_mu_hist_extraEl[i].Write()
    
    electronISO_Max_hist[i].Write()
    electronISO_Min_hist[i].Write()
    muonISO_Max_hist[i].Write()
    muonISO_Min_hist[i].Write()
    electronSIP_Max_hist[i].Write()
    electronSIP_Min_hist[i].Write()
    muonSIP_Max_hist[i].Write()
    muonSIP_Min_hist[i].Write()

    
output.Close()




