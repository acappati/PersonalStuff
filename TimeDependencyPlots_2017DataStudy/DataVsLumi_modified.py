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
import ROOT, helper, math
ROOT.gROOT.SetBatch(True)
from helper import DoSimpleFit, DoRooFit, ReadJSON, FillHisto, Result, GraphVsLumi, ZMultVsLumi, MeanRMSVsLumi, Ratio, SAME3VsLumi, SAME2VsLumi, TwoFileSAME3VsLumi,TwoFileSAME2VsLumi, PlotNpv
from ROOT import TFile, TH1F, TCanvas, gSystem, TLegend
from ROOT import kBlue, kRed, kBlack, kWhite

# *****************************
# Declare all the variables
# fit options 
fitMC       = False  #true for fitting MC in FitMC.py 
fitDATA     = False  #true for fitting DATA in FitDATA.py 
DoInclusive = False
applyPU2017 = False

# data tree options 
ZZTree   = False
CRZLTree = True
ZTree    = False

# data periods options
# period = "data2016"
period = "data2017"
# *****************************


if(period == "data2016"):
    inputTXT  = "JSON_calc/SplittedBlocks_2016data_0p5_new.txt"
    if(ZZTree):
        inputROOT = "2016data_0p5_reduced_ZZTree_histos.root"
        treeText  = "ZZTree"
    elif(CRZLTree):
        inputROOT = "2016data_0p5_reduced_CRZLTree_histos.root"
        treeText  = "CRZLTree"
    elif(ZTree):
        inputROOT = "2016data_0p5_reduced_ZTree_histos.root"
        treeText  = "ZTree"
    else:
        print ("Error: wrong option!")

elif(period == "data2017"):
    inputTXT  = "JSON_calc/SplittedBlocks_2017data_0p5_new.txt"
    if(ZZTree):
        inputROOT = "2017data_0p5_reduced_ZZTree_histos.root"
        treeText  = "ZZTree"
    elif(CRZLTree):
        inputROOT = "2017data_0p5_reduced_CRZLTree_histos.root"
        treeText  = "CRZLTree"
    elif(ZTree):
        inputROOT = "2017data_0p5_reduced_ZTree_histos.root"
        treeText  = "ZTree"
    else:
        print ("Error: wrong option!")

else: 
    print ("Error: choose a period!")




# ***********************************************
# MC and DATA fit for orizontal lines in graphs
# values from FitMC.py and FitDATA.py scripts
# 'Z_ele', 'Z_ele_extraMu', 'Z_ele_extraEl', 'Z_ele_EBEB', 'Z_ele_EBEE', 'Z_ele_EEEE', 'Z_mu', 'Z_mu_extraMu', 'Z_mu_extraEl', 'Z_mu_MBMB', 'Z_mu_MBME', 'Z_mu_MEME'

# 2016 data
if(period == "data2016"):

    massFitDATA  = [90.59790612299334, 90.56586193497857, 90.61123369417315, 90.68748261478433, 90.4105096437472, 90.40042309419971, 90.95284851134053, 90.95161309749076, 90.95341231648246, 91.01802298154813, 90.94863739783477, 90.85522027213003]
    widthFitDATA = [2.7793220292028646, 2.7983184869956923, 2.770902871058797, 2.437303194241703, 3.345729864688405, 3.5543577358221774, 2.1039712874410155, 2.095408440226358, 2.1070842206971414, 1.782698407932242, 2.130510238786742, 2.5762502815831216]

    if(applyPU2017):
        # MC with 2017 PU weight
        massFitMC  = [90.65328942151869, 90.63221459025714, 90.66226591929488, 90.73905897160829, 90.49687522164919, 90.44664928472986, 90.96881109491301, 90.97561425677574, 90.96789317642569, 91.03452262626013, 90.95418168038236, 90.90261382484313]
        widthFitMC =  [2.708906632827382, 2.7409155191468426, 2.6970765254112408, 2.3881037584654785, 3.1840490179797953, 3.5311904444581, 2.058956442696573, 2.0325794026858595, 2.0661613654651156, 1.7393828049165871, 2.0939116562892806, 2.4692045993046148]
    else:
        massFitMC  =  [90.64052467378977, 90.62190039429665, 90.64751531089786, 90.72771221866611, 90.49328659226262, 90.40547950810351, 90.96597099539777, 90.96029768841535, 90.96955739282213, 91.02290971026528, 90.95416866647793, 90.91310312395267]
        widthFitMC = [2.698804434077731, 2.7003390875493243, 2.6980006774989267, 2.3785367213968605, 3.1622727628487697, 3.5180734187921154, 2.056758233455251, 2.033399531416025, 2.061778269002473, 1.7351410697864615, 2.0847651732174604, 2.4627711316803427]



# 2017 data
if(period == "data2017"):

    if(CRZLTree) :
        # 'Z_ele', 'Z_ele_extraMu' , 'Z_ele_extraEl' , 'Z_ele_EBEB', 'Z_ele_EBEE', 'Z_ele_EEEE', 'Z_mu', 'Z_mu_extraMu', 'Z_mu_extraEl' , 'Z_mu_MBMB', 'Z_mu_MBME', 'Z_mu_MEME'
        massFitDATA  =  [89.88965685750605, 89.83604009222441, 89.91744063612994, 90.10866981966532, 89.25015974265675, 89.17309896072628, 90.95801807564867, 90.96723524458726, 90.95314598861657, 91.00873418485791, 90.96065733772602, 90.86880022623268]
        widthFitDATA =  [3.0601841820913336, 3.0888748014814182, 3.0449043972409995, 2.7199786078358645, 3.846217064066446, 4.317488374749054, 2.120339359868412, 2.1156000446794594, 2.1228544933260625, 1.802036357361082, 2.1420573120444226, 2.624797465456659]       


        massFitMC  =  [90.65277015222286, 90.66368802692999, 90.64576457104178, 90.73902956981016, 90.51708924382893, 90.40247924553357, 90.9863966398518, 91.00730553897829, 90.97773097050066, 91.05584099471746, 90.95908670972341, 90.95491497422672]
        widthFitMC =  [2.7607404335865833, 2.800870093633652, 2.7554494130630625, 2.4549688506168317, 3.205917852505671, 3.562119743243076, 2.0490557722630482, 2.0255669688610363, 2.0601441072750966, 1.7366790997575798, 2.1148660467970832, 2.4145429630724697]


    if(ZTree) : 
        # 'Z_ele', 'Z_ele_extraMu' = 0, 'Z_ele_extraEl' = 0, 'Z_ele_EBEB', 'Z_ele_EBEE', 'Z_ele_EEEE', 'Z_mu', 'Z_mu_extraMu' = 0, 'Z_mu_extraEl' = 0, 'Z_mu_MBMB', 'Z_mu_MBME', 'Z_mu_MEME'
        massFitDATA  = [ 89.7259289974 , 0.0 , 0.0 , 89.9863243635 , 89.0330157062 , 88.9058383479 , 90.9586964084 , 0.0 , 0.0 , 91.0098081907 , 90.9602975573 , 90.8731034006 ]
        widthFitDATA = [ 3.08812216575 , 0.0 , 0.0 , 2.7286689019 , 3.86955462913 , 4.35953081802 , 2.03838726915 , 0.0 , 0.0 , 1.71251696754 , 2.0527381423 , 2.53118453796 ]

        
        massFitMC  = [ 90.5294468571 , 0.0 , 0.0 , 90.6310074073 , 90.3581021096 , 90.3065609158 , 90.9632405666 , 0.0 , 0.0 , 91.0188009053 , 90.9594607923 , 90.894544852 ]
        widthFitMC = [ 2.80998665972 , 0.0 , 0.0 , 2.44645216736 , 3.32158184654 , 3.57923319861 , 2.02814146029 , 0.0 , 0.0 , 1.68756636642 , 2.05136385477 , 2.47269415969 ]

   

# ***********************************************




# make output directories
outDir_fit         = "FitResults_"  + str(period) + "_" + str(treeText)
if(applyPU2017): 
    outDir_ZPlots  = "ZPlots_" + str(period) + "_" + str(treeText) + "_2017PU"
else:
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
OutputPathISOSIP               = "ISOSIP_Plots_" + period + "_" + treeText
NameList_ISOSIPMax_inclusive   = ["ISOMax_inclusive","SIPMax_inclusive"]
TitleList_ISOSIPMax_inclusive  = ["ISO","SIP"]
HistoList_ISOSIPMax_Ele_inclusive  = []
HistoList_ISOSIPMax_Mu_inclusive   = []
gSystem.Exec("mkdir -p " + OutputPathISOSIP) # create output directory

HistoList_ISOSIPMax_Ele_inclusive.append(ElectronISO_F1_Max_hist[0])
HistoList_ISOSIPMax_Ele_inclusive.append(ElectronSIP_F1_Max_hist[0])
HistoList_ISOSIPMax_Mu_inclusive.append(MuonISO_F1_Max_hist[0])
HistoList_ISOSIPMax_Mu_inclusive.append(MuonSIP_F1_Max_hist[0])

for i in range(1,len(recorded)) :

    HistoList_ISOSIPMax_Ele_inclusive[0] += ElectronISO_F1_Max_hist[i]
    HistoList_ISOSIPMax_Ele_inclusive[1] += ElectronSIP_F1_Max_hist[1]
    HistoList_ISOSIPMax_Mu_inclusive[0]  += MuonISO_F1_Max_hist[i]
    HistoList_ISOSIPMax_Mu_inclusive[1]  += MuonSIP_F1_Max_hist[i]

# draw ISO and SIP plots
for j in range(0,len(HistoList_ISOSIPMax_Ele_inclusive)) :
    
    canvas = TCanvas("canvas","canvas",800,600)
    HistoList_ISOSIPMax_Mu_inclusive[j].SetTitle(TitleList_ISOSIPMax_inclusive[j])
    HistoList_ISOSIPMax_Mu_inclusive[j].SetXTitle(TitleList_ISOSIPMax_inclusive[j])
    HistoList_ISOSIPMax_Mu_inclusive[j].SetYTitle("events")
    
    HistoList_ISOSIPMax_Ele_inclusive[j].SetLineColor(kBlue)
    HistoList_ISOSIPMax_Ele_inclusive[j].Scale(1./HistoList_ISOSIPMax_Ele_inclusive[j].Integral())
    HistoList_ISOSIPMax_Mu_inclusive[j].Scale(1./ HistoList_ISOSIPMax_Mu_inclusive[j].Integral())
    HistoList_ISOSIPMax_Mu_inclusive[j].SetLineColor(kRed)
    canvas.cd()
    HistoList_ISOSIPMax_Mu_inclusive[j].Draw("histo")
    HistoList_ISOSIPMax_Ele_inclusive[j].Draw("histosame")
        
    legend = TLegend(0.75,0.76,0.98,0.95)
    legend.AddEntry(HistoList_ISOSIPMax_Ele_inclusive[j],"Electrons", "f")
    legend.AddEntry(HistoList_ISOSIPMax_Mu_inclusive[j], "Muons","f")
    legend.SetFillColor(kWhite)
    legend.SetLineColor(kBlack)
    legend.Draw()    

    canvas.Update()
    canvas.SaveAs(OutputPathISOSIP + "/" + NameList_ISOSIPMax_inclusive[j] + ".pdf")
    canvas.SaveAs(OutputPathISOSIP + "/" + NameList_ISOSIPMax_inclusive[j] + ".png")

#***************************************************



# DCB fit for ZPlots
Data2017_result2         = DoSimpleFit(ZMass_F1_ele_hist,         recorded, ZZTree, outDir_fit, period + "_Z_ele",         fitMC, fitDATA) # Z->e, full acceptance
Data2017_result2_EBEB    = DoSimpleFit(ZMass_F1_ele_hist_EBEB,    recorded, ZZTree, outDir_fit, period + "_Z_ele_EBEB",    fitMC, fitDATA) # Z->e, BB
Data2017_result2_EBEE    = DoSimpleFit(ZMass_F1_ele_hist_EBEE,    recorded, ZZTree, outDir_fit, period + "_Z_ele_EBEE",    fitMC, fitDATA) # Z->e, BE
Data2017_result2_EEEE    = DoSimpleFit(ZMass_F1_ele_hist_EEEE,    recorded, ZZTree, outDir_fit, period + "_Z_ele_EEEE",    fitMC, fitDATA) # Z->e, EE
if not ZTree :
    Data2017_result2_extraMu = DoSimpleFit(ZMass_F1_ele_hist_extraMu, recorded, ZZTree, outDir_fit, period + "_Z_ele_extraMu", fitMC, fitDATA) # Z->e, extra mu
    Data2017_result2_extraEl = DoSimpleFit(ZMass_F1_ele_hist_extraEl, recorded, ZZTree, outDir_fit, period + "_Z_ele_extraEl", fitMC, fitDATA) # Z->e, extra e

Data2017_result3         = DoSimpleFit(ZMass_F1_mu_hist,          recorded, ZZTree, outDir_fit, period + "_Z_mu",          fitMC, fitDATA) # Z->mu, full acceptance
Data2017_result3_MBMB    = DoSimpleFit(ZMass_F1_mu_hist_MBMB,     recorded, ZZTree, outDir_fit, period + "_Z_mu_MBMB",     fitMC, fitDATA) # Z->mu, BB
Data2017_result3_MBME    = DoSimpleFit(ZMass_F1_mu_hist_MBME,     recorded, ZZTree, outDir_fit, period + "_Z_mu_MBME",     fitMC, fitDATA) # Z->mu, BE
Data2017_result3_MEME    = DoSimpleFit(ZMass_F1_mu_hist_MEME,     recorded, ZZTree, outDir_fit, period + "_Z_mu_MEME",     fitMC, fitDATA) # Z->mu, EE
if not ZTree :
    Data2017_result3_extraMu = DoSimpleFit(ZMass_F1_mu_hist_extraMu,  recorded, ZZTree, outDir_fit, period + "_Z_mu_extraMu",  fitMC, fitDATA) # Z->mu, extra mu
    Data2017_result3_extraEl = DoSimpleFit(ZMass_F1_mu_hist_extraEl,  recorded, ZZTree, outDir_fit, period + "_Z_mu_extraEl",  fitMC, fitDATA) # Z->mu, extra e




Data2017_ele1, Data2017_ele2                 = GraphVsLumi(Data2017_result2,         outDir_ZPlots, period + "_Z_ele")
Data2017_mu1,  Data2017_mu2                  = GraphVsLumi(Data2017_result3,         outDir_ZPlots, period + "_Z_mu")
Data2017_ele_EBEB1, Data2017_ele_EBEB2       = GraphVsLumi(Data2017_result2_EBEB,    outDir_ZPlots, period + "_Z_ele_EBEB")
Data2017_ele_EBEE1, Data2017_ele_EBEE2       = GraphVsLumi(Data2017_result2_EBEE,    outDir_ZPlots, period + "_Z_ele_EBEE")
Data2017_ele_EEEE1, Data2017_ele_EEEE2       = GraphVsLumi(Data2017_result2_EEEE,    outDir_ZPlots, period + "_Z_ele_EEEE")
Data2017_mu_MBMB1, Data2017_mu_MBMB2         = GraphVsLumi(Data2017_result3_MBMB,    outDir_ZPlots, period + "_Z_mu_MBMB")
Data2017_mu_MBME1, Data2017_mu_MBME2         = GraphVsLumi(Data2017_result3_MBME,    outDir_ZPlots, period + "_Z_mu_MBME")
Data2017_mu_MEME1, Data2017_mu_MEME2         = GraphVsLumi(Data2017_result3_MEME,    outDir_ZPlots, period + "_Z_mu_MEME")
if not ZTree :
    Data2017_ele1_extraMu, Data2017_ele2_extraMu = GraphVsLumi(Data2017_result2_extraMu, outDir_ZPlots, period + "_Z_ele_extraMu")
    Data2017_mu1_extraMu,  Data2017_mu2_extraMu  = GraphVsLumi(Data2017_result3_extraMu, outDir_ZPlots, period + "_Z_mu_extraMu")
    Data2017_ele1_extraEl, Data2017_ele2_extraEl = GraphVsLumi(Data2017_result2_extraEl, outDir_ZPlots, period + "_Z_ele_extraEl")
    Data2017_mu1_extraEl,  Data2017_mu2_extraEl  = GraphVsLumi(Data2017_result3_extraEl, outDir_ZPlots, period + "_Z_mu_extraEl")



Data2017_eleM         = ZMultVsLumi(ZMass_F1_ele_hist,         recorded, outDir_ZPlots, period + "_Z_ele")
Data2017_muM          = ZMultVsLumi(ZMass_F1_mu_hist,          recorded, outDir_ZPlots, period + "_Z_mu")
Data2017_eleM_EBEB    = ZMultVsLumi(ZMass_F1_ele_hist_EBEB,    recorded, outDir_ZPlots, period + "_Z_ele_EBEB")
Data2017_eleM_EBEE    = ZMultVsLumi(ZMass_F1_ele_hist_EBEE,    recorded, outDir_ZPlots, period + "_Z_ele_EBEE")
Data2017_eleM_EEEE    = ZMultVsLumi(ZMass_F1_ele_hist_EEEE,    recorded, outDir_ZPlots, period + "_Z_ele_EEEE")
Data2017_muM_MBMB     = ZMultVsLumi(ZMass_F1_mu_hist_MBMB,     recorded, outDir_ZPlots, period + "_Z_mu_MBMB")
Data2017_muM_MBME     = ZMultVsLumi(ZMass_F1_mu_hist_MBME,     recorded, outDir_ZPlots, period + "_Z_mu_MBME")
Data2017_muM_MEME     = ZMultVsLumi(ZMass_F1_mu_hist_MEME,     recorded, outDir_ZPlots, period + "_Z_mu_MEME")
if not ZTree :
    Data2017_eleM_extraMu = ZMultVsLumi(ZMass_F1_ele_hist_extraMu, recorded, outDir_ZPlots, period + "_ele_extraMu")
    Data2017_muM_extraMu  = ZMultVsLumi(ZMass_F1_mu_hist_extraMu,  recorded, outDir_ZPlots, period + "_mu_extraMu")
    Data2017_eleM_extraEl = ZMultVsLumi(ZMass_F1_ele_hist_extraEl, recorded, outDir_ZPlots, period + "_ele_extraEl")
    Data2017_muM_extraEl  = ZMultVsLumi(ZMass_F1_mu_hist_extraEl,  recorded, outDir_ZPlots, period + "_mu_extraEl")


Data2017_eleISOMax1,Data2017_eleISOMax2 = MeanRMSVsLumi(ElectronISO_F1_Max_hist, recorded, outDir_LeptonPlots, period + "_ElectronISO_F1Max")
Data2017_eleISOMin1,Data2017_eleISOMin2 = MeanRMSVsLumi(ElectronISO_F1_Min_hist, recorded, outDir_LeptonPlots, period + "_ElectronISO_F1Min")
Data2017_muISOMax1,Data2017_muISOMax2   = MeanRMSVsLumi(MuonISO_F1_Max_hist,     recorded, outDir_LeptonPlots, period + "_MuonISO_F1Max")
Data2017_muISOMin1,Data2017_muISOMin2   = MeanRMSVsLumi(MuonISO_F1_Min_hist,     recorded, outDir_LeptonPlots, period + "_MuonISO_F1Min")
Data2017_eleSIPMax1,Data2017_eleSIPMax2 = MeanRMSVsLumi(ElectronSIP_F1_Max_hist, recorded, outDir_LeptonPlots, period + "_ElectronSIP_F1Max")
Data2017_eleSIPMin1,Data2017_eleSIPMin2 = MeanRMSVsLumi(ElectronSIP_F1_Min_hist, recorded, outDir_LeptonPlots, period + "_ElectronSIP_F1Min")
Data2017_muSIPMax1,Data2017_muSIPMax2   = MeanRMSVsLumi(MuonSIP_F1_Max_hist,     recorded, outDir_LeptonPlots, period + "_MuonSIP_F1Max")
Data2017_muSIPMin1,Data2017_muSIPMin2   = MeanRMSVsLumi(MuonSIP_F1_Min_hist,     recorded, outDir_LeptonPlots, period + "_MuonSIP_F1Min")



#Ele vs mu
SAME3VsLumi(None,Data2017_ele1,Data2017_mu1, str(outDir_ZPlots) + "/" + period + "_ZMass_ele_mu_together", "Zmass",  0.,0.,massFitMC[0], massFitDATA[0], massFitMC[6], massFitDATA[6], False, period) 
SAME3VsLumi(None,Data2017_ele2,Data2017_mu2, str(outDir_ZPlots) + "/" + period + "_ZWidth_ele_mu_together","Zwidth", 0.,0.,widthFitMC[0],widthFitDATA[0],widthFitMC[6],widthFitDATA[6],False, period) 
SAME3VsLumi(None,Data2017_eleM,Data2017_muM, str(outDir_ZPlots) + "/" + period + "_ZMult_ele_mu_together", "Zmult",  0.,0.,0.,0.,0.,0.,False, period)

#Ele vs mu, extra mu    
if not ZTree :
    SAME3VsLumi(None,Data2017_ele1_extraMu,Data2017_mu1_extraMu, str(outDir_ZPlots) + "/" + period + "_ZMass_ele_mu_together_extraMu", "Zmass",  0.,0.,massFitMC[1], massFitDATA[1], massFitMC[7], massFitDATA[7], False, period) 
    SAME3VsLumi(None,Data2017_ele2_extraMu,Data2017_mu2_extraMu, str(outDir_ZPlots) + "/" + period + "_ZWidth_ele_mu_together_extraMu","Zwidth", 0.,0.,widthFitMC[1],widthFitDATA[1],widthFitMC[7],widthFitDATA[7],False, period) 
    SAME3VsLumi(None,Data2017_eleM_extraMu,Data2017_muM_extraMu, str(outDir_ZPlots) + "/" + period + "_ZMult_ele_mu_together_extraMu", "Zmult",  0.,0.,0.,0.,0.,0.,False, period) 

#Ele vs mu, extra ele
if not ZTree :
    SAME3VsLumi(None,Data2017_ele1_extraEl,Data2017_mu1_extraEl, str(outDir_ZPlots) + "/" + period + "_ZMass_ele_mu_together_extraEl", "Zmass", 0.,0.,massFitMC[2], massFitDATA[2], massFitMC[8], massFitDATA[8], False, period)  
    SAME3VsLumi(None,Data2017_ele2_extraEl,Data2017_mu2_extraEl, str(outDir_ZPlots) + "/" + period + "_ZWidth_ele_mu_together_extraEl","Zwidth",0.,0.,widthFitMC[2],widthFitDATA[2],widthFitMC[8],widthFitDATA[8],False, period)
    SAME3VsLumi(None,Data2017_eleM_extraEl,Data2017_muM_extraEl, str(outDir_ZPlots) + "/" + period + "_ZMult_ele_mu_together_extraEl", "Zmult", 0.,0.,0.,0.,0.,0.,False, period)

#Ele, EBEB vs EBEE vs EEEE        
SAME3VsLumi(Data2017_ele_EBEB1,Data2017_ele_EBEE1,Data2017_ele_EEEE1, str(outDir_ZPlots) + "/" + period + "_ZMass_ele_EBEB_EBEE_EEEE", "Zmass", massFitMC[3], massFitDATA[3], massFitMC[4], massFitDATA[4], massFitMC[5], massFitDATA[5], True, period)  
SAME3VsLumi(Data2017_ele_EBEB2,Data2017_ele_EBEE2,Data2017_ele_EEEE2, str(outDir_ZPlots) + "/" + period + "_ZWidth_ele_EBEB_EBEE_EEEE","Zwidth",widthFitMC[3],widthFitDATA[3],widthFitMC[4],widthFitDATA[4],widthFitMC[5],widthFitDATA[5],True, period)
SAME3VsLumi(Data2017_eleM_EBEB,Data2017_eleM_EBEE,Data2017_eleM_EEEE, str(outDir_ZPlots) + "/" + period + "_ZMult_ele_EBEB_EBEE_EEEE", "Zmult", 0.,0.,0.,0.,0.,0.,True, period)

#Mu, MBMB vs MBME vs MEME        
SAME3VsLumi(Data2017_mu_MBMB1,Data2017_mu_MBME1,Data2017_mu_MEME1, str(outDir_ZPlots) + "/" + period + "_ZMass_mu_MBMB_MBME_MEME", "Zmass", massFitMC[9], massFitDATA[9], massFitMC[10], massFitDATA[10], massFitMC[11], massFitDATA[11], True, period) 
SAME3VsLumi(Data2017_mu_MBMB2,Data2017_mu_MBME2,Data2017_mu_MEME2, str(outDir_ZPlots) + "/" + period + "_ZWidth_mu_MBMB_MBME_MEME","Zwidth",widthFitMC[9],widthFitDATA[9],widthFitMC[10],widthFitDATA[10],widthFitMC[11],widthFitDATA[11],True, period)
SAME3VsLumi(Data2017_muM_MBMB,Data2017_muM_MBME,Data2017_muM_MEME, str(outDir_ZPlots) + "/" + period + "_ZMult_mu_MBMB_MBME_MEME", "Zmult", 0.,0.,0.,0.,0.,0.,True, period)

#ISO, SIP
SAME2VsLumi(Data2017_eleISOMax1,Data2017_muISOMax1, str(outDir_LeptonPlots) + "/" + period + "_MaxISO_ele_mu_together","ISO", period)
SAME2VsLumi(Data2017_eleISOMin1,Data2017_muISOMin1, str(outDir_LeptonPlots) + "/" + period + "_MinISO_ele_mu_together","ISO", period)
SAME2VsLumi(Data2017_eleSIPMax1,Data2017_muSIPMax1, str(outDir_LeptonPlots) + "/" + period + "_MaxSIP_ele_mu_together","SIP", period)
SAME2VsLumi(Data2017_eleSIPMin1,Data2017_muSIPMin1, str(outDir_LeptonPlots) + "/" + period + "_MinSIP_ele_mu_together","SIP", period)



