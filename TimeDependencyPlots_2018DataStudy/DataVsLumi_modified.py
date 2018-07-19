#!/usr/bin/env python


# *******************
# usage: 
#    python DATAVsLumi_modified.py (DATAVsLumi.py original verion from Toni)
#
# structure:
#    - read file with data histo divided per lumiblocks (created by ExtractHistos.py)
#    - fit mass and width histos with helper.py function
#    - do plots of Zmass, Zwidth, lepton ISO and SIP vs Lumi (using helper.py functions)
#    - save plots and fit results in 3 different directories: ZPlots, LeptonPlots, FitResults
# ********************


# First import ROOT and helper functions
import json
import ROOT, helper, math
ROOT.gROOT.SetBatch(True)
from helper import DoSimpleFit, DoDCBunbinnedFit, DoRooFit, ReadJSON, FillHisto, Result, GraphVsLumi, ZMultVsLumi, MeanRMSVsLumi, Ratio, SAME3VsLumi, SAME2VsLumi, TwoFileSAME3VsLumi,TwoFileSAME2VsLumi, PlotNpv
from ROOT import TFile, TH1F, TCanvas, gSystem, TLegend
from ROOT import kBlue, kRed, kBlack, kWhite

# *****************************
# Declare all the variables
# fit options 
fitMC       = False  #true for fitting MC in FitMC.py 
fitDATA     = False  #true for fitting DATA in FitDATA.py 
DoInclusive = False

# data tree options 
ZZTree   = False
CRZLTree = False
ZTree    = True

# data periods options
period = "data2018"
# *****************************


if(period == "data2018"):
    inputTXT  = "JSON_calc/SplittedBlocks_2018data_0p5_new.txt"
    if(ZZTree):
        inputROOT = "2018data_0p5_reduced_ZZTree_histos.root"
        treeText  = "ZZTree"
    elif(CRZLTree):
        inputROOT = "2018data_0p5_reduced_CRZLTree_histos.root"
        treeText  = "CRZLTree"
    elif(ZTree):
        inputROOT = "2018data_0p5_reduced_ZTree_histos.root"
        treeText  = "ZTree"
    else:
        print ("Error: wrong option!")

else: 
    print ("Error: choose a period!")



#******************************************************************************
# MC and DATA fit for orizontal lines in graphs
# read fit values from output files obtained by running FitDATA.py and FitMC.py
# 'Z_ele', 'Z_ele_extraMu', 'Z_ele_extraEl', 'Z_ele_EBEB', 'Z_ele_EBEE', 'Z_ele_EEEE', 'Z_mu', 'Z_mu_extraMu', 'Z_mu_extraEl', 'Z_mu_MBMB', 'Z_mu_MBME', 'Z_mu_MEME'

with open("out_ZDBCmean_DATA_" + period + "_" + treeText  + ".json","r") as handle1 :
    massFitDATA_dict = json.load(handle1)

with open("out_ZDBCwidth_DATA_" + period + "_" + treeText  + ".json","r") as handle2 :
    widthFitDATA_dict = json.load(handle2)

with open("out_ZDBCmean_MC_DYJets_" + period + "_" + treeText  + ".json","r") as handle3 :
    massFitMC_dict = json.load(handle3)

with open("out_ZDBCwidth_MC_DYJets_" + period + "_" + treeText  + ".json","r") as handle4 :
    widthFitMC_dict = json.load(handle4)

print 'Input files with Fit values read!'
#****************************************************************************




# make output directories
outDir_fit     = "FitResults_"  + str(period) + "_" + str(treeText)
outDir_ZPlots  = "ZPlots_" + str(period) + "_" + str(treeText)
outDir_LeptonPlots = "LeptonPlots_" + str(period) + "_" + str(treeText)

gSystem.Exec("mkdir -p " + outDir_fit)
gSystem.Exec("mkdir -p " + outDir_ZPlots)
gSystem.Exec("mkdir -p " + outDir_LeptonPlots)
print "Output directories created!"

data = TFile(inputROOT)

RunNum_B  = []
LumiNum_B = []
RunNum_E  = []
LumiNum_E = []
recorded  = []

# Read the input JSON file and store it
ReadJSON(inputTXT, RunNum_B, LumiNum_B, RunNum_E, LumiNum_E, recorded)

#define histos lists (each list element correspond to a different lumiblock)
ZMass_F1_ele_hist         = []
if not ZTree :
    ZMass_F1_ele_hist_extraMu = []
    ZMass_F1_ele_hist_extraEl = []
ZMass_F1_ele_hist_EBEB    = []
ZMass_F1_ele_hist_EBEE    = []
ZMass_F1_ele_hist_EEEE    = []

ZMass_F1_mu_hist         = []
if not ZTree :
    ZMass_F1_mu_hist_extraMu = []
    ZMass_F1_mu_hist_extraEl = []
ZMass_F1_mu_hist_MBMB    = []
ZMass_F1_mu_hist_MBME    = []
ZMass_F1_mu_hist_MEME    = []

ElectronISO_F1_Max_hist = []
ElectronISO_F1_Min_hist = []
MuonISO_F1_Max_hist     = []
MuonISO_F1_Min_hist     = []
ElectronSIP_F1_Max_hist = []
ElectronSIP_F1_Min_hist = []
MuonSIP_F1_Max_hist     = []
MuonSIP_F1_Min_hist     = []

# fill histos lists 
for i in range(0,len(recorded)):

    ZMass_F1_ele_hist.append(data.Get('ZMass_ele' + str(i)))
    ZMass_F1_ele_hist_EBEB.append(data.Get('ZMass_ele_EBEB_' + str(i)))
    ZMass_F1_ele_hist_EBEE.append(data.Get('ZMass_ele_EBEE_' + str(i)))
    ZMass_F1_ele_hist_EEEE.append(data.Get('ZMass_ele_EEEE_' + str(i)))
    if not ZTree :
        ZMass_F1_ele_hist_extraMu.append(data.Get('ZMass_ele_extraMu' + str(i)))
        ZMass_F1_ele_hist_extraEl.append(data.Get('ZMass_ele_extraEl' + str(i)))
    
    ZMass_F1_mu_hist.append(data.Get('ZMass_mu' + str(i)))
    ZMass_F1_mu_hist_MBMB.append(data.Get('ZMass_mu_MBMB_' + str(i)))
    ZMass_F1_mu_hist_MBME.append(data.Get('ZMass_mu_MBME_' + str(i)))
    ZMass_F1_mu_hist_MEME.append(data.Get('ZMass_mu_MEME_' + str(i)))
    if not ZTree :
        ZMass_F1_mu_hist_extraMu.append(data.Get('ZMass_mu_extraMu' + str(i)))
        ZMass_F1_mu_hist_extraEl.append(data.Get('ZMass_mu_extraEl' + str(i)))
    
    ElectronISO_F1_Max_hist.append(data.Get('ElectronISO_Max' + str(i)))
    ElectronISO_F1_Min_hist.append(data.Get('ElectronISO_Min' + str(i)))
    MuonISO_F1_Max_hist.append(data.Get('MuonISO_Max' + str(i)))
    MuonISO_F1_Min_hist.append(data.Get('MuonISO_Min' + str(i)))
    ElectronSIP_F1_Max_hist.append(data.Get('ElectronSIP_Max' + str(i)))
    ElectronSIP_F1_Min_hist.append(data.Get('ElectronSIP_Min' + str(i)))
    MuonSIP_F1_Max_hist.append(data.Get('MuonSIP_Max' + str(i)))
    MuonSIP_F1_Min_hist.append(data.Get('MuonSIP_Min' + str(i)))



#***************************************************
#inclusive ISO and SIP histos 
OutputPathISOSIP = "ISOSIP_Plots_" + period + "_" + treeText
gSystem.Exec("mkdir -p " + OutputPathISOSIP) # create output directory


for i in range(0,len(recorded)) :

     c1 = TCanvas("c1","c1",800,600)
     ElectronISO_F1_Max_hist[i].SetTitle(ElectronISO_F1_Max_hist[i].GetName())
     ElectronISO_F1_Max_hist[i].SetXTitle("ele Max ISO")
     ElectronISO_F1_Max_hist[i].SetYTitle("events")
     ElectronISO_F1_Max_hist[i].Draw("histo")
     c1.SaveAs(OutputPathISOSIP + "/" + ElectronISO_F1_Max_hist[i].GetName() + ".pdf")
     c1.SaveAs(OutputPathISOSIP + "/" + ElectronISO_F1_Max_hist[i].GetName() + ".png")
     c1.SaveAs(OutputPathISOSIP + "/" + ElectronISO_F1_Max_hist[i].GetName() + ".root")

     c2 = TCanvas("c2","c2",800,600)
     ElectronISO_F1_Min_hist[i].SetTitle(ElectronISO_F1_Min_hist[i].GetName())
     ElectronISO_F1_Min_hist[i].SetXTitle("ele Min ISO")
     ElectronISO_F1_Min_hist[i].SetYTitle("events")
     ElectronISO_F1_Min_hist[i].Draw("histo")
     c2.SaveAs(OutputPathISOSIP + "/" + ElectronISO_F1_Min_hist[i].GetName() + ".pdf")
     c2.SaveAs(OutputPathISOSIP + "/" + ElectronISO_F1_Min_hist[i].GetName() + ".png")
     c2.SaveAs(OutputPathISOSIP + "/" + ElectronISO_F1_Min_hist[i].GetName() + ".root")

     c3 = TCanvas("c3","c3",800,600)
     ElectronSIP_F1_Max_hist[i].SetTitle(ElectronSIP_F1_Max_hist[i].GetName())
     ElectronSIP_F1_Max_hist[i].SetXTitle("ele Max SIP")
     ElectronSIP_F1_Max_hist[i].SetYTitle("events")
     ElectronSIP_F1_Max_hist[i].Draw("histo")
     c3.SaveAs(OutputPathISOSIP + "/" + ElectronSIP_F1_Max_hist[i].GetName() + ".pdf")
     c3.SaveAs(OutputPathISOSIP + "/" + ElectronSIP_F1_Max_hist[i].GetName() + ".png")
     c3.SaveAs(OutputPathISOSIP + "/" + ElectronSIP_F1_Max_hist[i].GetName() + ".root")

     c4 = TCanvas("c3","c3",800,600)
     ElectronSIP_F1_Min_hist[i].SetTitle(ElectronSIP_F1_Min_hist[i].GetName())
     ElectronSIP_F1_Min_hist[i].SetXTitle("ele Min SIP")
     ElectronSIP_F1_Min_hist[i].SetYTitle("events")
     ElectronSIP_F1_Min_hist[i].Draw("histo")
     c4.SaveAs(OutputPathISOSIP + "/" + ElectronSIP_F1_Min_hist[i].GetName() + ".pdf")
     c4.SaveAs(OutputPathISOSIP + "/" + ElectronSIP_F1_Min_hist[i].GetName() + ".png")
     c4.SaveAs(OutputPathISOSIP + "/" + ElectronSIP_F1_Min_hist[i].GetName() + ".root")

     c5 = TCanvas("c1","c1",800,600)
     MuonISO_F1_Max_hist[i].SetTitle(MuonISO_F1_Max_hist[i].GetName())
     MuonISO_F1_Max_hist[i].SetXTitle("mu Max ISO")
     MuonISO_F1_Max_hist[i].SetYTitle("events")
     MuonISO_F1_Max_hist[i].Draw("histo")
     c5.SaveAs(OutputPathISOSIP + "/" + MuonISO_F1_Max_hist[i].GetName() + ".pdf")
     c5.SaveAs(OutputPathISOSIP + "/" + MuonISO_F1_Max_hist[i].GetName() + ".png")
     c5.SaveAs(OutputPathISOSIP + "/" + MuonISO_F1_Max_hist[i].GetName() + ".root")

     c6 = TCanvas("c2","c2",800,600)
     MuonISO_F1_Min_hist[i].SetTitle(MuonISO_F1_Min_hist[i].GetName())
     MuonISO_F1_Min_hist[i].SetXTitle("mu Min ISO")
     MuonISO_F1_Min_hist[i].SetYTitle("events")
     MuonISO_F1_Min_hist[i].Draw("histo")
     c6.SaveAs(OutputPathISOSIP + "/" + MuonISO_F1_Min_hist[i].GetName() + ".pdf")
     c6.SaveAs(OutputPathISOSIP + "/" + MuonISO_F1_Min_hist[i].GetName() + ".png")
     c6.SaveAs(OutputPathISOSIP + "/" + MuonISO_F1_Min_hist[i].GetName() + ".root")

     c7 = TCanvas("c3","c3",800,600)
     MuonSIP_F1_Max_hist[i].SetTitle(MuonSIP_F1_Max_hist[i].GetName())
     MuonSIP_F1_Max_hist[i].SetXTitle("mu Max SIP")
     MuonSIP_F1_Max_hist[i].SetYTitle("events")
     MuonSIP_F1_Max_hist[i].Draw("histo")
     c7.SaveAs(OutputPathISOSIP + "/" + MuonSIP_F1_Max_hist[i].GetName() + ".pdf")
     c7.SaveAs(OutputPathISOSIP + "/" + MuonSIP_F1_Max_hist[i].GetName() + ".png")
     c7.SaveAs(OutputPathISOSIP + "/" + MuonSIP_F1_Max_hist[i].GetName() + ".root")

     c8 = TCanvas("c3","c3",800,600)
     MuonSIP_F1_Min_hist[i].SetTitle(MuonSIP_F1_Min_hist[i].GetName())
     MuonSIP_F1_Min_hist[i].SetXTitle("mu Min SIP")
     MuonSIP_F1_Min_hist[i].SetYTitle("events")
     MuonSIP_F1_Min_hist[i].Draw("histo")
     c8.SaveAs(OutputPathISOSIP + "/" + MuonSIP_F1_Min_hist[i].GetName() + ".pdf")
     c8.SaveAs(OutputPathISOSIP + "/" + MuonSIP_F1_Min_hist[i].GetName() + ".png")
     c8.SaveAs(OutputPathISOSIP + "/" + MuonSIP_F1_Min_hist[i].GetName() + ".root")


#***************************************************



# DCB fit for ZPlots,

# binned fit
Data2018_result2         = DoSimpleFit(ZMass_F1_ele_hist,         recorded, ZZTree, outDir_fit, period + "_Z_ele",         fitMC, fitDATA) # Z->e, full acceptance
Data2018_result2_EBEB    = DoSimpleFit(ZMass_F1_ele_hist_EBEB,    recorded, ZZTree, outDir_fit, period + "_Z_ele_EBEB",    fitMC, fitDATA) # Z->e, BB
Data2018_result2_EBEE    = DoSimpleFit(ZMass_F1_ele_hist_EBEE,    recorded, ZZTree, outDir_fit, period + "_Z_ele_EBEE",    fitMC, fitDATA) # Z->e, BE
Data2018_result2_EEEE    = DoSimpleFit(ZMass_F1_ele_hist_EEEE,    recorded, ZZTree, outDir_fit, period + "_Z_ele_EEEE",    fitMC, fitDATA) # Z->e, EE
if not ZTree :
    Data2018_result2_extraMu = DoSimpleFit(ZMass_F1_ele_hist_extraMu, recorded, ZZTree, outDir_fit, period + "_Z_ele_extraMu", fitMC, fitDATA) # Z->e, extra mu
    Data2018_result2_extraEl = DoSimpleFit(ZMass_F1_ele_hist_extraEl, recorded, ZZTree, outDir_fit, period + "_Z_ele_extraEl", fitMC, fitDATA) # Z->e, extra e

Data2018_result3         = DoSimpleFit(ZMass_F1_mu_hist,          recorded, ZZTree, outDir_fit, period + "_Z_mu",          fitMC, fitDATA) # Z->mu, full acceptance
Data2018_result3_MBMB    = DoSimpleFit(ZMass_F1_mu_hist_MBMB,     recorded, ZZTree, outDir_fit, period + "_Z_mu_MBMB",     fitMC, fitDATA) # Z->mu, BB
Data2018_result3_MBME    = DoSimpleFit(ZMass_F1_mu_hist_MBME,     recorded, ZZTree, outDir_fit, period + "_Z_mu_MBME",     fitMC, fitDATA) # Z->mu, BE
Data2018_result3_MEME    = DoSimpleFit(ZMass_F1_mu_hist_MEME,     recorded, ZZTree, outDir_fit, period + "_Z_mu_MEME",     fitMC, fitDATA) # Z->mu, EE
if not ZTree :
    Data2018_result3_extraMu = DoSimpleFit(ZMass_F1_mu_hist_extraMu,  recorded, ZZTree, outDir_fit, period + "_Z_mu_extraMu",  fitMC, fitDATA) # Z->mu, extra mu
    Data2018_result3_extraEl = DoSimpleFit(ZMass_F1_mu_hist_extraEl,  recorded, ZZTree, outDir_fit, period + "_Z_mu_extraEl",  fitMC, fitDATA) # Z->mu, extra e

# unbinned fit - to be fixed
# Data2018_result2         = DoDCBunbinnedFit(ZMass_F1_ele_hist,         recorded, ZZTree, outDir_fit, period + "_Z_ele",         fitMC, fitDATA) # Z->e, full acceptance
# Data2018_result2_EBEB    = DoDCBunbinnedFit(ZMass_F1_ele_hist_EBEB,    recorded, ZZTree, outDir_fit, period + "_Z_ele_EBEB",    fitMC, fitDATA) # Z->e, BB
# Data2018_result2_EBEE    = DoDCBunbinnedFit(ZMass_F1_ele_hist_EBEE,    recorded, ZZTree, outDir_fit, period + "_Z_ele_EBEE",    fitMC, fitDATA) # Z->e, BE
# Data2018_result2_EEEE    = DoDCBunbinnedFit(ZMass_F1_ele_hist_EEEE,    recorded, ZZTree, outDir_fit, period + "_Z_ele_EEEE",    fitMC, fitDATA) # Z->e, EE
# if not ZTree :
#     Data2018_result2_extraMu = DoDCBunbinnedFit(ZMass_F1_ele_hist_extraMu, recorded, ZZTree, outDir_fit, period + "_Z_ele_extraMu", fitMC, fitDATA) # Z->e, extra mu
#     Data2018_result2_extraEl = DoDCBunbinnedFit(ZMass_F1_ele_hist_extraEl, recorded, ZZTree, outDir_fit, period + "_Z_ele_extraEl", fitMC, fitDATA) # Z->e, extra e

# Data2018_result3         = DoDCBunbinnedFit(ZMass_F1_mu_hist,          recorded, ZZTree, outDir_fit, period + "_Z_mu",          fitMC, fitDATA) # Z->mu, full acceptance
# Data2018_result3_MBMB    = DoDCBunbinnedFit(ZMass_F1_mu_hist_MBMB,     recorded, ZZTree, outDir_fit, period + "_Z_mu_MBMB",     fitMC, fitDATA) # Z->mu, BB
# Data2018_result3_MBME    = DoDCBunbinnedFit(ZMass_F1_mu_hist_MBME,     recorded, ZZTree, outDir_fit, period + "_Z_mu_MBME",     fitMC, fitDATA) # Z->mu, BE
# Data2018_result3_MEME    = DoDCBunbinnedFit(ZMass_F1_mu_hist_MEME,     recorded, ZZTree, outDir_fit, period + "_Z_mu_MEME",     fitMC, fitDATA) # Z->mu, EE
# if not ZTree :
#     Data2018_result3_extraMu = DoDCBunbinnedFit(ZMass_F1_mu_hist_extraMu,  recorded, ZZTree, outDir_fit, period + "_Z_mu_extraMu",  fitMC, fitDATA) # Z->mu, extra mu
#     Data2018_result3_extraEl = DoDCBunbinnedFit(ZMass_F1_mu_hist_extraEl,  recorded, ZZTree, outDir_fit, period + "_Z_mu_extraEl",  fitMC, fitDATA) # Z->mu, extra e





Data2018_ele1, Data2018_ele2                 = GraphVsLumi(Data2018_result2,         outDir_ZPlots, period + "_Z_ele")
Data2018_mu1,  Data2018_mu2                  = GraphVsLumi(Data2018_result3,         outDir_ZPlots, period + "_Z_mu")
Data2018_ele_EBEB1, Data2018_ele_EBEB2       = GraphVsLumi(Data2018_result2_EBEB,    outDir_ZPlots, period + "_Z_ele_EBEB")
Data2018_ele_EBEE1, Data2018_ele_EBEE2       = GraphVsLumi(Data2018_result2_EBEE,    outDir_ZPlots, period + "_Z_ele_EBEE")
Data2018_ele_EEEE1, Data2018_ele_EEEE2       = GraphVsLumi(Data2018_result2_EEEE,    outDir_ZPlots, period + "_Z_ele_EEEE")
Data2018_mu_MBMB1, Data2018_mu_MBMB2         = GraphVsLumi(Data2018_result3_MBMB,    outDir_ZPlots, period + "_Z_mu_MBMB")
Data2018_mu_MBME1, Data2018_mu_MBME2         = GraphVsLumi(Data2018_result3_MBME,    outDir_ZPlots, period + "_Z_mu_MBME")
Data2018_mu_MEME1, Data2018_mu_MEME2         = GraphVsLumi(Data2018_result3_MEME,    outDir_ZPlots, period + "_Z_mu_MEME")
if not ZTree :
    Data2018_ele1_extraMu, Data2018_ele2_extraMu = GraphVsLumi(Data2018_result2_extraMu, outDir_ZPlots, period + "_Z_ele_extraMu")
    Data2018_mu1_extraMu,  Data2018_mu2_extraMu  = GraphVsLumi(Data2018_result3_extraMu, outDir_ZPlots, period + "_Z_mu_extraMu")
    Data2018_ele1_extraEl, Data2018_ele2_extraEl = GraphVsLumi(Data2018_result2_extraEl, outDir_ZPlots, period + "_Z_ele_extraEl")
    Data2018_mu1_extraEl,  Data2018_mu2_extraEl  = GraphVsLumi(Data2018_result3_extraEl, outDir_ZPlots, period + "_Z_mu_extraEl")



Data2018_eleM         = ZMultVsLumi(ZMass_F1_ele_hist,         recorded, outDir_ZPlots, period + "_Z_ele")
Data2018_muM          = ZMultVsLumi(ZMass_F1_mu_hist,          recorded, outDir_ZPlots, period + "_Z_mu")
Data2018_eleM_EBEB    = ZMultVsLumi(ZMass_F1_ele_hist_EBEB,    recorded, outDir_ZPlots, period + "_Z_ele_EBEB")
Data2018_eleM_EBEE    = ZMultVsLumi(ZMass_F1_ele_hist_EBEE,    recorded, outDir_ZPlots, period + "_Z_ele_EBEE")
Data2018_eleM_EEEE    = ZMultVsLumi(ZMass_F1_ele_hist_EEEE,    recorded, outDir_ZPlots, period + "_Z_ele_EEEE")
Data2018_muM_MBMB     = ZMultVsLumi(ZMass_F1_mu_hist_MBMB,     recorded, outDir_ZPlots, period + "_Z_mu_MBMB")
Data2018_muM_MBME     = ZMultVsLumi(ZMass_F1_mu_hist_MBME,     recorded, outDir_ZPlots, period + "_Z_mu_MBME")
Data2018_muM_MEME     = ZMultVsLumi(ZMass_F1_mu_hist_MEME,     recorded, outDir_ZPlots, period + "_Z_mu_MEME")
if not ZTree :
    Data2018_eleM_extraMu = ZMultVsLumi(ZMass_F1_ele_hist_extraMu, recorded, outDir_ZPlots, period + "_ele_extraMu")
    Data2018_muM_extraMu  = ZMultVsLumi(ZMass_F1_mu_hist_extraMu,  recorded, outDir_ZPlots, period + "_mu_extraMu")
    Data2018_eleM_extraEl = ZMultVsLumi(ZMass_F1_ele_hist_extraEl, recorded, outDir_ZPlots, period + "_ele_extraEl")
    Data2018_muM_extraEl  = ZMultVsLumi(ZMass_F1_mu_hist_extraEl,  recorded, outDir_ZPlots, period + "_mu_extraEl")


Data2018_eleISOMax1,Data2018_eleISOMax2 = MeanRMSVsLumi(ElectronISO_F1_Max_hist, recorded, outDir_LeptonPlots, period + "_ElectronISO_F1Max")
Data2018_eleISOMin1,Data2018_eleISOMin2 = MeanRMSVsLumi(ElectronISO_F1_Min_hist, recorded, outDir_LeptonPlots, period + "_ElectronISO_F1Min")
Data2018_muISOMax1,Data2018_muISOMax2   = MeanRMSVsLumi(MuonISO_F1_Max_hist,     recorded, outDir_LeptonPlots, period + "_MuonISO_F1Max")
Data2018_muISOMin1,Data2018_muISOMin2   = MeanRMSVsLumi(MuonISO_F1_Min_hist,     recorded, outDir_LeptonPlots, period + "_MuonISO_F1Min")
Data2018_eleSIPMax1,Data2018_eleSIPMax2 = MeanRMSVsLumi(ElectronSIP_F1_Max_hist, recorded, outDir_LeptonPlots, period + "_ElectronSIP_F1Max")
Data2018_eleSIPMin1,Data2018_eleSIPMin2 = MeanRMSVsLumi(ElectronSIP_F1_Min_hist, recorded, outDir_LeptonPlots, period + "_ElectronSIP_F1Min")
Data2018_muSIPMax1,Data2018_muSIPMax2   = MeanRMSVsLumi(MuonSIP_F1_Max_hist,     recorded, outDir_LeptonPlots, period + "_MuonSIP_F1Max")
Data2018_muSIPMin1,Data2018_muSIPMin2   = MeanRMSVsLumi(MuonSIP_F1_Min_hist,     recorded, outDir_LeptonPlots, period + "_MuonSIP_F1Min")



#Ele vs mu
SAME3VsLumi(None,Data2018_ele1,Data2018_mu1, str(outDir_ZPlots) + "/" + period + "_ZMass_ele_mu_together", "Zmass",  0.,0.,massFitMC_dict['Zee'], massFitDATA_dict['Zee'], massFitMC_dict['Zmumu'], massFitDATA_dict['Zmumu'], False, period) 
SAME3VsLumi(None,Data2018_ele2,Data2018_mu2, str(outDir_ZPlots) + "/" + period + "_ZWidth_ele_mu_together","Zwidth", 0.,0.,widthFitMC_dict['Zee'],widthFitDATA_dict['Zee'],widthFitMC_dict['Zmumu'],widthFitDATA_dict['Zmumu'],False, period) 
SAME3VsLumi(None,Data2018_eleM,Data2018_muM, str(outDir_ZPlots) + "/" + period + "_ZMult_ele_mu_together", "Zmult",  0.,0.,0.,0.,0.,0.,False, period)

#Ele vs mu, extra mu    
if not ZTree :
    SAME3VsLumi(None,Data2018_ele1_extraMu,Data2018_mu1_extraMu, str(outDir_ZPlots) + "/" + period + "_ZMass_ele_mu_together_extraMu", "Zmass",  0.,0.,massFitMC_dict['Zee_extraMu'], massFitDATA_dict['Zee_extraMu'], massFitMC_dict['Zmumu_extraMu'], massFitDATA_dict['Zmumu_extraMu'], False, period) 
    SAME3VsLumi(None,Data2018_ele2_extraMu,Data2018_mu2_extraMu, str(outDir_ZPlots) + "/" + period + "_ZWidth_ele_mu_together_extraMu","Zwidth", 0.,0.,widthFitMC_dict['Zee_extraMu'],widthFitDATA_dict['Zee_extraMu'],widthFitMC_dict['Zmumu_extraMu'],widthFitDATA_dict['Zmumu_extraMu'],False, period) 
    SAME3VsLumi(None,Data2018_eleM_extraMu,Data2018_muM_extraMu, str(outDir_ZPlots) + "/" + period + "_ZMult_ele_mu_together_extraMu", "Zmult",  0.,0.,0.,0.,0.,0.,False, period) 

#Ele vs mu, extra ele
if not ZTree :
    SAME3VsLumi(None,Data2018_ele1_extraEl,Data2018_mu1_extraEl, str(outDir_ZPlots) + "/" + period + "_ZMass_ele_mu_together_extraEl", "Zmass", 0.,0.,massFitMC_dict['Zee_extraEl'], massFitDATA_dict['Zee_extraEl'], massFitMC_dict['Zmumu_extraEl'], massFitDATA_dict['Zmumu_extraEl'], False, period)  
    SAME3VsLumi(None,Data2018_ele2_extraEl,Data2018_mu2_extraEl, str(outDir_ZPlots) + "/" + period + "_ZWidth_ele_mu_together_extraEl","Zwidth",0.,0.,widthFitMC_dict['Zee_extraEl'],widthFitDATA_dict['Zee_extraEl'],widthFitMC_dict['Zmumu_extraEl'],widthFitDATA_dict['Zmumu_extraEl'],False, period)
    SAME3VsLumi(None,Data2018_eleM_extraEl,Data2018_muM_extraEl, str(outDir_ZPlots) + "/" + period + "_ZMult_ele_mu_together_extraEl", "Zmult", 0.,0.,0.,0.,0.,0.,False, period)

#Ele, EBEB vs EBEE vs EEEE        
SAME3VsLumi(Data2018_ele_EBEB1,Data2018_ele_EBEE1,Data2018_ele_EEEE1, str(outDir_ZPlots) + "/" + period + "_ZMass_ele_EBEB_EBEE_EEEE", "Zmass", massFitMC_dict['Zee_EBEB'], massFitDATA_dict['Zee_EBEB'], massFitMC_dict['Zee_EBEE'], massFitDATA_dict['Zee_EBEE'], massFitMC_dict['Zee_EEEE'], massFitDATA_dict['Zee_EEEE'], True, period)  
SAME3VsLumi(Data2018_ele_EBEB2,Data2018_ele_EBEE2,Data2018_ele_EEEE2, str(outDir_ZPlots) + "/" + period + "_ZWidth_ele_EBEB_EBEE_EEEE","Zwidth",widthFitMC_dict['Zee_EBEB'],widthFitDATA_dict['Zee_EBEB'],widthFitMC_dict['Zee_EBEE'],widthFitDATA_dict['Zee_EBEE'],widthFitMC_dict['Zee_EEEE'],widthFitDATA_dict['Zee_EEEE'],True, period)
SAME3VsLumi(Data2018_eleM_EBEB,Data2018_eleM_EBEE,Data2018_eleM_EEEE, str(outDir_ZPlots) + "/" + period + "_ZMult_ele_EBEB_EBEE_EEEE", "Zmult", 0.,0.,0.,0.,0.,0.,True, period)

#Mu, MBMB vs MBME vs MEME        
SAME3VsLumi(Data2018_mu_MBMB1,Data2018_mu_MBME1,Data2018_mu_MEME1, str(outDir_ZPlots) + "/" + period + "_ZMass_mu_MBMB_MBME_MEME", "Zmass", massFitMC_dict['Zmumu_MBMB'], massFitDATA_dict['Zmumu_MBMB'], massFitMC_dict['Zmumu_MBME'], massFitDATA_dict['Zmumu_MBME'], massFitMC_dict['Zmumu_MEME'], massFitDATA_dict['Zmumu_MEME'], True, period) 
SAME3VsLumi(Data2018_mu_MBMB2,Data2018_mu_MBME2,Data2018_mu_MEME2, str(outDir_ZPlots) + "/" + period + "_ZWidth_mu_MBMB_MBME_MEME","Zwidth",widthFitMC_dict['Zmumu_MBMB'],widthFitDATA_dict['Zmumu_MBMB'],widthFitMC_dict['Zmumu_MBME'],widthFitDATA_dict['Zmumu_MBME'],widthFitMC_dict['Zmumu_MEME'],widthFitDATA_dict['Zmumu_MEME'],True, period)
SAME3VsLumi(Data2018_muM_MBMB,Data2018_muM_MBME,Data2018_muM_MEME, str(outDir_ZPlots) + "/" + period + "_ZMult_mu_MBMB_MBME_MEME", "Zmult", 0.,0.,0.,0.,0.,0.,True, period)

#ISO, SIP
SAME2VsLumi(Data2018_eleISOMax1,Data2018_muISOMax1, str(outDir_LeptonPlots) + "/" + period + "_MaxISO_ele_mu_together","ISO", period)
SAME2VsLumi(Data2018_eleISOMin1,Data2018_muISOMin1, str(outDir_LeptonPlots) + "/" + period + "_MinISO_ele_mu_together","ISO", period)
SAME2VsLumi(Data2018_eleSIPMax1,Data2018_muSIPMax1, str(outDir_LeptonPlots) + "/" + period + "_MaxSIP_ele_mu_together","SIP", period)
SAME2VsLumi(Data2018_eleSIPMin1,Data2018_muSIPMin1, str(outDir_LeptonPlots) + "/" + period + "_MinSIP_ele_mu_together","SIP", period)



