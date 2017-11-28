#! /usr/bin/env python
import sys
import os
import re
import math
import yaml
from syncUtils import *
from scipy.special import erf
import ROOT
from ROOT import *
from ctypes import *
from array import *


#**********************************************
# usage:
#    ./plotWorkpacePDFs.py
#
# notes:
#    - first specify the input (for MC and workspaces) and output paths 
#    - and set redoHistos = True to re-do histograms from MC files
#
#    - this script works inside the test directory of the ZZAnalysis release
#
#**********************************************

lib = CDLL('libZZAnalysisAnalysisStep.so')

#---------------
# luminosity
lumi = 35.86706

# Mass interval
m4lMin = 105.
m4lMax = 140.

# Bin number
nBin = 100
binWidth = (m4lMax - m4lMin) / nBin

# Final States: 
nFinalStates = 3
sFinalStates = ["4mu", "4e", "2e2mu"]    
fs4mu   = 0 
fs4e    = 1 
fs2e2mu = 2

# Categories: 
nCategories = 7
sCategories = ['UnTagged','VBF1jTagged','VBF2jTagged','VHLeptTagged','VHHadrTagged','ttHTagged','VHMETTagged']

# Datasets
nDatasets = 13
sDatasets = ['ggH125', 'VBFH125', 'WplusH125', 'WminusH125', 'ZH125', 'ttH125', 'ZZTo4l', 'ggTo4e_Contin_MCFM701', 'ggTo4mu_Contin_MCFM701', 'ggTo4tau_Contin_MCFM701', 'ggTo2e2mu_Contin_MCFM701', 'ggTo2e2tau_Contin_MCFM701', 'ggTo2mu2tau_Contin_MCFM701']


# Processes
ggH    = 0
qqH    = 1
WH     = 2 
WH_lep = 3 
WH_had = 4 
ZH     = 5 
ZH_lep = 6 
ZH_had = 7 
ttH    = 8
qqZZ   = 9 
ggZZ   = 10

#-------------------------------
# Options 
redoHistos = False


#-------------------------------------------------
# input path where there are workspaces to be read 
#inputWSPath = "/afs/cern.ch/work/a/acappati/H4l/170404/Datacards13TeV_correctShapes/STXSCards/"      #correct shapes
#inputWSPath = '/afs/cern.ch/work/a/acappati/H4l/170404/Datacards13TeV_170626/STXSCards/'             #HIG-16-041 shapes
inputWSPath = '/afs/cern.ch/work/a/acappati/H4l/170404/Datacards13TeV/STXSCards/'

# input path where there are MC files
inputMCPath = "/data3/Higgs/170526_SystStudy/"

# output path where to save plots
outputPath1 = "ShapesFromWS"
gSystem.Exec("mkdir -p "+outputPath1)
outputPath2 = "ShapesFromWS/ShapesFromWS_singleShape"
gSystem.Exec("mkdir -p "+outputPath2)


#--------------------
# ggH, qqH, WH_had, WH_lep, ZH_had, ZH_lep histograms
histo_ggH    = [[TH1F('ggH_fs'+sFinalStates[fs]+'_cat'+str(cat),'ggH_fs'+sFinalStates[fs]+'_cat'+str(cat),nBin,m4lMin,m4lMax) for cat in range(nCategories)] for fs in range(nFinalStates)] 
histo_qqH    = [[TH1F('qqH_fs'+sFinalStates[fs]+'_cat'+str(cat),'qqH_fs'+sFinalStates[fs]+'_cat'+str(cat),nBin,m4lMin,m4lMax) for cat in range(nCategories)] for fs in range(nFinalStates)] 
histo_WH_had = [[TH1F('WH_had_fs'+sFinalStates[fs]+'_cat'+str(cat),'WH_had_fs'+sFinalStates[fs]+'_cat'+str(cat),nBin,m4lMin,m4lMax) for cat in range(nCategories)] for fs in range(nFinalStates)] 
histo_WH_lep = [[TH1F('WH_lep_fs'+sFinalStates[fs]+'_cat'+str(cat),'WH_lep_fs'+sFinalStates[fs]+'_cat'+str(cat),nBin,m4lMin,m4lMax) for cat in range(nCategories)] for fs in range(nFinalStates)]
histo_ZH_had = [[TH1F('ZH_had_fs'+sFinalStates[fs]+'_cat'+str(cat),'ZH_had_fs'+sFinalStates[fs]+'_cat'+str(cat),nBin,m4lMin,m4lMax) for cat in range(nCategories)] for fs in range(nFinalStates)] 
histo_ZH_lep = [[TH1F('ZH_lep_fs'+sFinalStates[fs]+'_cat'+str(cat),'ZH_lep_fs'+sFinalStates[fs]+'_cat'+str(cat),nBin,m4lMin,m4lMax) for cat in range(nCategories)] for fs in range(nFinalStates)] 
histo_ttH    = [[TH1F('ttH_fs'+sFinalStates[fs]+'_cat'+str(cat),'ttH_fs'+sFinalStates[fs]+'_cat'+str(cat),nBin,m4lMin,m4lMax) for cat in range(nCategories)] for fs in range(nFinalStates)] 
histo_qqZZ   = [[TH1F('qqZZ_fs'+sFinalStates[fs]+'_cat'+str(cat),'qqZZ_fs'+sFinalStates[fs]+'_cat'+str(cat),nBin,m4lMin,m4lMax) for cat in range(nCategories)] for fs in range(nFinalStates)] 
histo_ggZZ   = [[TH1F('ggZZ_fs'+sFinalStates[fs]+'_cat'+str(cat),'ggZZ_fs'+sFinalStates[fs]+'_cat'+str(cat),nBin,m4lMin,m4lMax) for cat in range(nCategories)] for fs in range(nFinalStates)] 

#------------------------
# bin width for rebinning
binWidth_ggH    = [[0 for cat in range(nCategories)] for fs in range(nFinalStates)] 
binWidth_qqH    = [[0 for cat in range(nCategories)] for fs in range(nFinalStates)] 
binWidth_WH_had = [[0 for cat in range(nCategories)] for fs in range(nFinalStates)] 
binWidth_WH_lep = [[0 for cat in range(nCategories)] for fs in range(nFinalStates)]
binWidth_ZH_had = [[0 for cat in range(nCategories)] for fs in range(nFinalStates)] 
binWidth_ZH_lep = [[0 for cat in range(nCategories)] for fs in range(nFinalStates)] 
binWidth_ttH    = [[0 for cat in range(nCategories)] for fs in range(nFinalStates)] 
binWidth_qqZZ   = [[0 for cat in range(nCategories)] for fs in range(nFinalStates)] 
binWidth_ggZZ   = [[0 for cat in range(nCategories)] for fs in range(nFinalStates)] 


#----------------
# read MC files
#----------------
if(redoHistos) :

    # Loop over datasets
    for d in range(nDatasets) :

        currentProcess = -1
        if sDatasets[d]=='ggH125'                     : currentProcess = ggH 
        if sDatasets[d]=='VBFH125'                    : currentProcess = qqH
        if sDatasets[d]=='WplusH125'                  : currentProcess = WH
        if sDatasets[d]=='WminusH125'                 : currentProcess = WH
        if sDatasets[d]=='ZH125'                      : currentProcess = ZH
        if sDatasets[d]=='ttH125'                     : currentProcess = ttH 
        if sDatasets[d]=='ZZTo4l'                     : currentProcess = qqZZ
        if sDatasets[d]=='ggTo4e_Contin_MCFM701'      : currentProcess = ggZZ
        if sDatasets[d]=='ggTo4mu_Contin_MCFM701'     : currentProcess = ggZZ
        if sDatasets[d]=='ggTo4tau_Contin_MCFM701'    : currentProcess = ggZZ
        if sDatasets[d]=='ggTo2e2mu_Contin_MCFM701'   : currentProcess = ggZZ
        if sDatasets[d]=='ggTo2e2tau_Contin_MCFM701'  : currentProcess = ggZZ
        if sDatasets[d]=='ggTo2mu2tau_Contin_MCFM701' : currentProcess = ggZZ

        initialCurrentProcess = currentProcess
        
        # Open and read MC file    
        inFileMC = ROOT.TFile.Open(inputMCPath + sDatasets[d] + "/ZZ4lAnalysis.root")
        inFileMC.cd()
        
        hCounters = inFileMC.Get('ZZTree/Counters') 
        gen_sumWeights = hCounters.GetBinContent(41) #use PUweight
        partialSampleWeight = lumi * 1000 / gen_sumWeights
        
        
        jets30pt = []
        jets30eta = []
        jets30phi = []
        jets30mass = []
        jets30QGLikelihood = []

        # Process tree 
        inTree = inFileMC.Get("ZZTree/candTree")
        
        print 'Processing dataset',sDatasets[d],'(',inTree.GetEntriesFast(),') ...'
        nEntries = 0
        for entry in inTree:
        
            nEntries += 1
        
            if nEntries % 10000 == 0: 
                print 'Processing entry :', nEntries 
        
            # event selection
            ZZsel =  entry.ZZsel
            if not(ZZsel>=90) : continue 
        
            ZZMass = entry.ZZMass
            if (ZZMass<m4lMin or ZZMass>m4lMax) : continue
        
            nRun           = entry.RunNumber
            nLumi          = entry.LumiNumber
            nEvent         = entry.EventNumber
            xsec           = entry.xsec
            genHEPMCweight = entry.genHEPMCweight
            dataMCWeight   = entry.dataMCWeight
            Z1Flav         = entry.Z1Flav
            Z2Flav         = entry.Z2Flav
            pfMet          = entry.PFMET
            genExtInfo     = entry.genExtInfo

            # kfactor
            kfactor = 1
            if currentProcess==qqZZ : kfactor = entry.KFactor_EW_qqZZ * entry.KFactor_QCD_qqZZ_M
            if currentProcess==ggZZ : kfactor = entry.KFactor_QCD_ggZZ_Nominal
            
            eventWeight = partialSampleWeight * xsec * kfactor * genHEPMCweight*dataMCWeight #use PUweight
            
            
            # fill arrays from TBranches
            for i in range(len(entry.JetPt)):
                if entry.JetPt[i]>30.:
                    jets30pt.append(entry.JetPt[i])
                    jets30eta.append(entry.JetEta[i])
                    jets30phi.append(entry.JetPhi[i])
                    jets30mass.append(entry.JetMass[i])
                    jets30QGLikelihood.append(entry.JetQGLikelihood[i])
            
            # find subprocess if needed
            if initialCurrentProcess==WH : 
                if genExtInfo>10 : currentProcess = WH_lep
                else : currentProcess = WH_had
            if initialCurrentProcess==ZH :
                if genExtInfo>10 : currentProcess = ZH_lep
                else : currentProcess = ZH_had
            
            # find final state 
            currentFinalState = -1
            if Z1Flav==-121 :
                if Z2Flav==-121 : currentFinalState = fs4e 
                elif Z2Flav==-169 : currentFinalState = fs2e2mu
                else : print 'error in event ',nRun,':',nLumi,':',nEvent,' , Z2Flav = ',Z2Flav
            elif Z1Flav==-169 : 
                if Z2Flav==-121 : currentFinalState = fs2e2mu
                elif Z2Flav==-169 : currentFinalState = fs4mu
                else : print 'error in event ',nRun,':',nLumi,':',nEvent,' , Z2Flav = ',Z2Flav
            else : print 'error in event ',nRun,':',nLumi,':',nEvent,' , Z1Flav = ',Z1Flav  

            # find category
            # Moriond2017 categories
            currentCategory = lib.categoryMor17(
                c_int(entry.nExtraLep),
                c_int(entry.nExtraZ),
                c_int(entry.nCleanedJetsPt30),
                c_int(entry.nCleanedJetsPt30BTagged),
                (c_float * len(jets30QGLikelihood)).from_buffer(array('f',jets30QGLikelihood)),
                c_float(entry.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal),
                c_float(entry.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal),
                c_float(entry.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal),
                c_float(entry.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal),
                c_float(entry.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal),
                c_float(entry.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal),
                c_float(entry.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal),
                (c_float * len(jets30phi)).from_buffer(array('f',jets30phi)),
	        c_float(ZZMass),
                c_float(pfMet),
                c_bool(True), #useVHMETTagged
                c_bool(False) #useQGTagging
                )

            # print nRun,':',nLumi,':',nEvent, 'fs',currentFinalState,'cat',currentCategory
        


            # fill histograms
            if currentProcess==ggH    : histo_ggH[currentFinalState][currentCategory].Fill(ZZMass,eventWeight)
            if currentProcess==qqH    : histo_qqH[currentFinalState][currentCategory].Fill(ZZMass,eventWeight)
            if currentProcess==WH_had : histo_WH_had[currentFinalState][currentCategory].Fill(ZZMass,eventWeight)
            if currentProcess==WH_lep : histo_WH_lep[currentFinalState][currentCategory].Fill(ZZMass,eventWeight)
            if currentProcess==ZH_had : histo_ZH_had[currentFinalState][currentCategory].Fill(ZZMass,eventWeight)
            if currentProcess==ZH_lep : histo_ZH_lep[currentFinalState][currentCategory].Fill(ZZMass,eventWeight)
            if currentProcess==ttH    : histo_ttH[currentFinalState][currentCategory].Fill(ZZMass,eventWeight)
            if currentProcess==qqZZ   : histo_qqZZ[currentFinalState][currentCategory].Fill(ZZMass,eventWeight)
            if currentProcess==ggZZ   : histo_ggZZ[currentFinalState][currentCategory].Fill(ZZMass,eventWeight)

    
        

        
    # Save histograms in a root file
    fOutHisto = ROOT.TFile('HistoMC.root','recreate')
    fOutHisto.cd()
    for fs in range(nFinalStates) :
        for cat in range(nCategories) :
            histo_ggH[fs][cat].Write()
            histo_qqH[fs][cat].Write()
            histo_WH_had[fs][cat].Write()
            histo_WH_lep[fs][cat].Write()
            histo_ZH_had[fs][cat].Write()
            histo_ZH_lep[fs][cat].Write()
            histo_ttH[fs][cat].Write()
            histo_qqZZ[fs][cat].Write()
            histo_ggZZ[fs][cat].Write()
    fOutHisto.Close()


#-----------------------------------------------------------
# read histograms from root file
fInHisto = ROOT.TFile.Open('HistoMC.root')
fInHisto.cd()
print 'Reading file', fInHisto.GetName(),'...'
for fs in range(nFinalStates) :
    for cat in range(nCategories) :
        #ggH
        histo_ggH[fs][cat]    = fInHisto.Get('ggH_fs'+sFinalStates[fs]+'_cat'+str(cat))
        ggH_entries = histo_ggH[fs][cat].GetEntries()
        if (ggH_entries<10000 and ggH_entries>=1000) : histo_ggH[fs][cat].Rebin(2)
        if (ggH_entries<500  and ggH_entries>=100)   : histo_ggH[fs][cat].Rebin(4)
        if ggH_entries<100                           : histo_ggH[fs][cat].Rebin(5)
        binWidth_ggH[fs][cat] = histo_ggH[fs][cat].GetBinCenter(2) - histo_ggH[fs][cat].GetBinCenter(1)

        #qqH
        histo_qqH[fs][cat]    = fInHisto.Get('qqH_fs'+sFinalStates[fs]+'_cat'+str(cat))
        qqH_entries = histo_qqH[fs][cat].GetEntries()
        if (qqH_entries<10000 and qqH_entries>=500) : histo_qqH[fs][cat].Rebin(2)
        if (qqH_entries<500   and qqH_entries>=100) : histo_qqH[fs][cat].Rebin(4)
        if qqH_entries<100                          : histo_qqH[fs][cat].Rebin(5)
        binWidth_qqH[fs][cat] = histo_qqH[fs][cat].GetBinCenter(2) - histo_qqH[fs][cat].GetBinCenter(1)
        #WH_had
        histo_WH_had[fs][cat] = fInHisto.Get('WH_had_fs'+sFinalStates[fs]+'_cat'+str(cat))
        WH_had_entries = histo_WH_had[fs][cat].GetEntries()
        if (WH_had_entries<10000 and WH_had_entries>=500) : histo_WH_had[fs][cat].Rebin(2)
        if (WH_had_entries<500   and WH_had_entries>=100) : histo_WH_had[fs][cat].Rebin(4)
        if WH_had_entries<100                             : histo_WH_had[fs][cat].Rebin(5)
        binWidth_WH_had[fs][cat] = histo_WH_had[fs][cat].GetBinCenter(2) - histo_WH_had[fs][cat].GetBinCenter(1)
        #WH_lep
        histo_WH_lep[fs][cat] = fInHisto.Get('WH_lep_fs'+sFinalStates[fs]+'_cat'+str(cat))
        WH_lep_entries = histo_WH_lep[fs][cat].GetEntries()
        if (WH_lep_entries<10000 and WH_lep_entries>=500) : histo_WH_lep[fs][cat].Rebin(2)
        if (WH_lep_entries<500   and WH_lep_entries>=100) : histo_WH_lep[fs][cat].Rebin(4)
        if WH_lep_entries<100                             : histo_WH_lep[fs][cat].Rebin(5)
        binWidth_WH_lep[fs][cat] = histo_WH_lep[fs][cat].GetBinCenter(2) - histo_WH_lep[fs][cat].GetBinCenter(1)
        #ZH_had
        histo_ZH_had[fs][cat] = fInHisto.Get('ZH_had_fs'+sFinalStates[fs]+'_cat'+str(cat))
        ZH_had_entries = histo_ZH_had[fs][cat].GetEntries()
        if (ZH_had_entries<10000 and ZH_had_entries>=500) : histo_ZH_had[fs][cat].Rebin(2)
        if (ZH_had_entries<500   and ZH_had_entries>=100) : histo_ZH_had[fs][cat].Rebin(4)
        if ZH_had_entries<100                             : histo_ZH_had[fs][cat].Rebin(5)
        binWidth_ZH_had[fs][cat] = histo_ZH_had[fs][cat].GetBinCenter(2) - histo_ZH_had[fs][cat].GetBinCenter(1)
        #ZH_lep
        histo_ZH_lep[fs][cat] = fInHisto.Get('ZH_lep_fs'+sFinalStates[fs]+'_cat'+str(cat))
        ZH_lep_entries = histo_ZH_lep[fs][cat].GetEntries()
        if (ZH_lep_entries<10000 and ZH_lep_entries>=500) : histo_ZH_lep[fs][cat].Rebin(2)
        if (ZH_lep_entries<500   and ZH_lep_entries>=100) : histo_ZH_lep[fs][cat].Rebin(4)
        if ZH_lep_entries<100                             : histo_ZH_lep[fs][cat].Rebin(5)
        binWidth_ZH_lep[fs][cat] = histo_ZH_lep[fs][cat].GetBinCenter(2) - histo_ZH_lep[fs][cat].GetBinCenter(1)
        #ttH
        histo_ttH[fs][cat]    = fInHisto.Get('ttH_fs'+sFinalStates[fs]+'_cat'+str(cat))
        ttH_entries = histo_ttH[fs][cat].GetEntries()
        if (ttH_entries<10000 and ttH_entries>=500) : histo_ttH[fs][cat].Rebin(2)
        if (ttH_entries<500   and ttH_entries>=100) : histo_ttH[fs][cat].Rebin(4)
        if ttH_entries<100                          : histo_ttH[fs][cat].Rebin(5)
        binWidth_ttH[fs][cat] = histo_ttH[fs][cat].GetBinCenter(2) - histo_ttH[fs][cat].GetBinCenter(1)
        #qqZZ
        histo_qqZZ[fs][cat]   = fInHisto.Get('qqZZ_fs'+sFinalStates[fs]+'_cat'+str(cat))
        qqZZ_entries = histo_qqZZ[fs][cat].GetEntries()
        if qqZZ_entries>=10000                        : histo_qqZZ[fs][cat].Rebin(2)
        if (qqZZ_entries<10000 and qqZZ_entries>=500) : histo_qqZZ[fs][cat].Rebin(4)
        if (qqZZ_entries<500   and qqZZ_entries>=100) : histo_qqZZ[fs][cat].Rebin(5)
        if qqZZ_entries<100                           : histo_qqZZ[fs][cat].Rebin(6)
        binWidth_qqZZ[fs][cat] = histo_qqZZ[fs][cat].GetBinCenter(2) - histo_qqZZ[fs][cat].GetBinCenter(1)
        #ggZZ
        histo_ggZZ[fs][cat]   = fInHisto.Get('ggZZ_fs'+sFinalStates[fs]+'_cat'+str(cat))
        ggZZ_entries = histo_ggZZ[fs][cat].GetEntries()
        if ggZZ_entries>=10000                        : histo_ggZZ[fs][cat].Rebin(2)
        if (ggZZ_entries<10000 and ggZZ_entries>=500) : histo_ggZZ[fs][cat].Rebin(4)
        if (ggZZ_entries<500   and ggZZ_entries>=100) : histo_ggZZ[fs][cat].Rebin(5)
        if ggZZ_entries<100                           : histo_ggZZ[fs][cat].Rebin(6)
        binWidth_ggZZ[fs][cat] = histo_ggZZ[fs][cat].GetBinCenter(2) - histo_ggZZ[fs][cat].GetBinCenter(1)

# ------------------------------------
#read workspaces and print shapes
#------------------------------------
for fs in range(nFinalStates) : 
    for cat in range(nCategories) :
        #open file
        fInput = ROOT.TFile(inputWSPath+"hzz4lcard_%s_%i.input.root" % (sFinalStates[fs], cat))
        fInput.cd()
        print 'Reading file',fInput.GetName(),'...'

        #access to workspace
        w = fInput.Get("w")
        w.Print() #print workspace informations

        # data = w.data("data_obs")
        # get PDFs from the WS
        ggHShape      = w.pdf("ggH_hzz_mass")
        qqHShape      = w.pdf("qqH_hzz_mass")
        WH_hadShape   = w.pdf("WH_had_hzz_mass") 
        #WH_had_Landau = w.pdf('WH_hadLandau')
        #WH_had_CB     = w.pdf('WH_hadCB') 
        WH_lepShape   = w.pdf("WH_lep_hzz_mass")
        WH_lep_Landau = w.pdf('WH_lepLandau')
        WH_lep_CB     = w.pdf('WH_lepCB')
        ZH_hadShape   = w.pdf("ZH_had_hzz_mass")
        #ZH_had_Landau = w.pdf('ZH_hadLandau')
        #ZH_had_CB     = w.pdf('ZH_hadCB')
        ZH_lepShape   = w.pdf("ZH_lep_hzz_mass")
        ZH_lep_Landau = w.pdf('ZH_lepLandau')
        ZH_lep_CB     = w.pdf('ZH_lepCB')
        ttHShape      = w.pdf("ttH_hzz_mass")
        ttH_Landau    = w.pdf('ttHLandau')
        ttH_CB        = w.pdf('ttHCB')
        qqZZShape     = w.pdf("qqZZ_hzz_mass")
        ggZZShape     = w.pdf("ggZZ_hzz_mass")
        zjetsShape    = w.pdf("zjets_hzz_mass")

        # get signal normalization from the WS
        ggH_norm    = w.function("ggH_hzz_norm").getVal() #evaluate the RooFormlaVar ggH_norm at the nominal value of the parameter MH=125
        qqH_norm    = w.function("qqH_hzz_norm").getVal()
        WH_had_norm = w.function("WH_had_hzz_norm").getVal()
        #WH_had_frac = w.function('fracWH_had_'+sFinalStates[fs]+'_13').getVal() #frac is the normalization for the 2 components: CB and Landau
        WH_lep_norm = w.function("WH_lep_hzz_norm").getVal()
        #WH_lep_frac = w.function('fracWH_lep_'+sFinalStates[fs]+'_13').getVal()
        WH_lep_frac = w.function('fracWH_lep').getVal()
        ZH_had_norm = w.function("ZH_had_hzz_norm").getVal()
        #ZH_had_frac = w.function('fracZH_had_'+sFinalStates[fs]+'_13').getVal()
        ZH_lep_norm = w.function("ZH_lep_hzz_norm").getVal()
        #ZH_lep_frac = w.function('fracZH_lep_'+sFinalStates[fs]+'_13').getVal()
        ZH_lep_frac = w.function('fracZH_lep').getVal()
        ttH_norm    = w.function("ttH_hzz_norm").getVal()
        #ttH_frac    = w.function('fracttH_'+sFinalStates[fs]+'_13').getVal()
        ttH_frac    = w.function('fracttH').getVal()

        # normalization from yaml files
        fileNormInput = open(inputWSPath+'configs/inputs/yields_per_tag_category_13TeV_'+ sFinalStates[fs] +'.yaml')
        normYamlDict  = yaml.load(fileNormInput) #create a dictionary with the content of the yaml file
        qqZZ_norm     = float(normYamlDict[sCategories[cat]]['qqZZ_hzz'])
        ggZZ_norm     = float(normYamlDict[sCategories[cat]]['ggZZ_hzz'])
        
        mass4l = w.var("mass4l")
        # kd = w.var("kd")
        
        #plots
        plotggH    = mass4l.frame()
        plotggH.SetTitle("ggH")
        plotqqH    = mass4l.frame()
        plotqqH.SetTitle("qqH")
        plotWH_had = mass4l.frame()
        plotWH_had.SetTitle("WH_had")
        plotWH_lep = mass4l.frame()
        plotWH_lep.SetTitle("WH_lep")
        plotZH_had = mass4l.frame()
        plotZH_had.SetTitle("ZH_had")
        plotZH_lep = mass4l.frame()
        plotZH_lep.SetTitle("ZH_lep")
        plotttH    = mass4l.frame()
        plotttH.SetTitle("ttH")
        plotqqZZ   = mass4l.frame()
        plotqqZZ.SetTitle("qqZZ")
        plotggZZ   = mass4l.frame()
        plotggZZ.SetTitle("ggZZ")
        plotzjets  = mass4l.frame()
        plotzjets.SetTitle("zjets")
        
        # data.plotOn(plot)
        ggHShape.plotOn(plotggH,         RooFit.Normalization(ggH_norm*binWidth_ggH[fs][cat], RooAbsReal.NumEvent), RooFit.LineColor(kOrange))

        qqHShape.plotOn(plotqqH,         RooFit.Normalization(qqH_norm*binWidth_qqH[fs][cat], RooAbsReal.NumEvent), RooFit.LineColor(kPink-6))

        #WH_had_Landau.plotOn(plotWH_had, RooFit.Normalization(WH_had_norm*(1-WH_had_frac)*binWidth_WH_had[fs][cat], RooAbsReal.NumEvent), RooFit.LineColor(kMagenta+2))
        #WH_had_CB.plotOn(plotWH_had,     RooFit.Normalization(WH_had_norm*WH_had_frac*binWidth_WH_had[fs][cat],     RooAbsReal.NumEvent), RooFit.LineColor(kCyan-6))
        WH_hadShape.plotOn(plotWH_had,   RooFit.Normalization(WH_had_norm*binWidth_WH_had[fs][cat],                 RooAbsReal.NumEvent), RooFit.LineColor(kGreen+2))

        WH_lep_Landau.plotOn(plotWH_lep, RooFit.Normalization(WH_lep_norm*(1-WH_lep_frac)*binWidth_WH_lep[fs][cat], RooAbsReal.NumEvent), RooFit.LineColor(kMagenta+2))
        WH_lep_CB.plotOn(plotWH_lep,     RooFit.Normalization(WH_lep_norm*WH_lep_frac*binWidth_WH_lep[fs][cat],     RooAbsReal.NumEvent), RooFit.LineColor(kCyan-6))
        WH_lepShape.plotOn(plotWH_lep,   RooFit.Normalization(WH_lep_norm*binWidth_WH_lep[fs][cat],                 RooAbsReal.NumEvent), RooFit.LineColor(kBlue-7))

        #ZH_had_Landau.plotOn(plotZH_had, RooFit.Normalization(ZH_had_norm*(1-ZH_had_frac)*binWidth_ZH_had[fs][cat], RooAbsReal.NumEvent), RooFit.LineColor(kMagenta+2))
        #ZH_had_CB.plotOn(plotZH_had,     RooFit.Normalization(ZH_had_norm*ZH_had_frac*binWidth_ZH_had[fs][cat],     RooAbsReal.NumEvent), RooFit.LineColor(kCyan-6))
        ZH_hadShape.plotOn(plotZH_had,   RooFit.Normalization(ZH_had_norm*binWidth_ZH_had[fs][cat],                 RooAbsReal.NumEvent), RooFit.LineColor(kOrange+1))

        ZH_lep_Landau.plotOn(plotZH_lep, RooFit.Normalization(ZH_lep_norm*(1-ZH_lep_frac)*binWidth_ZH_lep[fs][cat], RooAbsReal.NumEvent), RooFit.LineColor(kMagenta+2))
        ZH_lep_CB.plotOn(plotZH_lep,     RooFit.Normalization(ZH_lep_norm*ZH_lep_frac*binWidth_ZH_lep[fs][cat],     RooAbsReal.NumEvent), RooFit.LineColor(kCyan-6))
        ZH_lepShape.plotOn(plotZH_lep,   RooFit.Normalization(ZH_lep_norm*binWidth_ZH_lep[fs][cat],                 RooAbsReal.NumEvent), RooFit.LineColor(6))

        ttH_Landau.plotOn(plotttH,       RooFit.Normalization(ttH_norm*(1-ttH_frac)*binWidth_ttH[fs][cat],    RooAbsReal.NumEvent), RooFit.LineColor(kMagenta+2))
        ttH_CB.plotOn(plotttH,           RooFit.Normalization(ttH_norm*ttH_frac*binWidth_ttH[fs][cat],        RooAbsReal.NumEvent), RooFit.LineColor(kCyan-6))
        ttHShape.plotOn(plotttH,         RooFit.Normalization(ttH_norm*binWidth_ttH[fs][cat],                 RooAbsReal.NumEvent), RooFit.LineColor(kOrange-3))

        qqZZShape.plotOn(plotqqZZ,       RooFit.Normalization(qqZZ_norm*binWidth_qqZZ[fs][cat], RooAbsReal.NumEvent), RooFit.LineColor(kGreen-3))

        ggZZShape.plotOn(plotggZZ,       RooFit.Normalization(ggZZ_norm*binWidth_ggZZ[fs][cat],  RooAbsReal.NumEvent), RooFit.LineColor(kOrange+7))

        zjetsShape.plotOn(plotzjets, RooFit.LineColor(kMagenta+1))
        #kdplot = kd.frame()
        #data.plotOn(kdplot)
        

        # plot single canvas
        canvas1 = ROOT.TCanvas('canvas1','canvas1',600,600)
        canvas1.cd()
        histo_ggH[fs][cat].Draw()
        plotggH.Draw('same')
        canvas1.SaveAs(outputPath2+"/ggH_%s_%i.png" % (sFinalStates[fs], cat))
        canvas1.SaveAs(outputPath2+"/ggH_%s_%i.pdf" % (sFinalStates[fs], cat))

        canvas2 = ROOT.TCanvas('canvas2','canvas2',600,600)
        canvas2.cd()
        histo_qqH[fs][cat].Draw()
        plotqqH.Draw('same')
        canvas2.SaveAs(outputPath2+"/qqH_%s_%i.png" % (sFinalStates[fs], cat))
        canvas2.SaveAs(outputPath2+"/qqH_%s_%i.pdf" % (sFinalStates[fs], cat))

        canvas3 = ROOT.TCanvas('canvas3','canvas3',600,600)
        canvas3.cd()
        histo_WH_had[fs][cat].Draw()
        plotWH_had.Draw('same')
        canvas3.SaveAs(outputPath2+"/WH_had_%s_%i.png" % (sFinalStates[fs], cat))
        canvas3.SaveAs(outputPath2+"/WH_had_%s_%i.pdf" % (sFinalStates[fs], cat))

        canvas4 = ROOT.TCanvas('canvas4','canvas4',600,600)
        canvas4.cd()
        histo_WH_lep[fs][cat].Draw()
        plotWH_lep.Draw('same')
        canvas4.SaveAs(outputPath2+"/WH_lep_%s_%i.png" % (sFinalStates[fs], cat))
        canvas4.SaveAs(outputPath2+"/WH_lep_%s_%i.pdf" % (sFinalStates[fs], cat))

        canvas5 = ROOT.TCanvas('canvas5','canvas5',600,600)
        canvas5.cd()
        histo_ZH_had[fs][cat].Draw()
        plotZH_had.Draw('same')
        canvas5.SaveAs(outputPath2+"/ZH_had_%s_%i.png" % (sFinalStates[fs], cat))
        canvas5.SaveAs(outputPath2+"/ZH_had_%s_%i.pdf" % (sFinalStates[fs], cat))

        canvas6 = ROOT.TCanvas('canvas6','canvas6',600,600)
        canvas6.cd()
        histo_ZH_lep[fs][cat].Draw()
        plotZH_lep.Draw('same')
        canvas6.SaveAs(outputPath2+"/ZH_lep_%s_%i.png" % (sFinalStates[fs], cat))
        canvas6.SaveAs(outputPath2+"/ZH_lep_%s_%i.pdf" % (sFinalStates[fs], cat))

        canvas7 = ROOT.TCanvas('canvas7','canvas7',600,600)
        canvas7.cd()
        histo_ttH[fs][cat].Draw()
        plotttH.Draw('same')
        canvas7.SaveAs(outputPath2+"/ttH_%s_%i.png" % (sFinalStates[fs], cat))
        canvas7.SaveAs(outputPath2+"/ttH_%s_%i.pdf" % (sFinalStates[fs], cat))

        canvas8 = ROOT.TCanvas('canvas8','canvas8',600,600)
        canvas8.cd()
        histo_qqZZ[fs][cat].Draw()
        plotqqZZ.Draw('same')
        canvas8.SaveAs(outputPath2+"/qqZZ_%s_%i.png" % (sFinalStates[fs], cat))
        canvas8.SaveAs(outputPath2+"/qqZZ_%s_%i.pdf" % (sFinalStates[fs], cat))

        canvas9 = ROOT.TCanvas('canvas9','canvas9',600,600)
        canvas9.cd()
        histo_ggZZ[fs][cat].Draw()
        plotggZZ.Draw('same')
        canvas9.SaveAs(outputPath2+"/ggZZ_%s_%i.png" % (sFinalStates[fs], cat))
        canvas9.SaveAs(outputPath2+"/ggZZ_%s_%i.pdf" % (sFinalStates[fs], cat))

        canvas10 = ROOT.TCanvas('canvas10','canvas10',600,600)
        canvas10.cd()
        plotzjets.Draw()
        canvas10.SaveAs(outputPath2+"/zjets_%s_%i.png" % (sFinalStates[fs], cat))
        canvas10.SaveAs(outputPath2+"/zjets_%s_%i.pdf" % (sFinalStates[fs], cat))
        
        # plot multiple canvas
        canvas = ROOT.TCanvas('canvas','canvas',1300,700)
        canvas.Divide(5,2)
        canvas.cd(1)
        histo_ggH[fs][cat].Draw()
        plotggH.Draw('same')
        canvas.cd(2)
        histo_qqH[fs][cat].Draw()
        plotqqH.Draw('same')
        canvas.cd(3)
        histo_WH_had[fs][cat].Draw()
        plotWH_had.Draw('same')
        canvas.cd(4)
        histo_WH_lep[fs][cat].Draw()
        plotWH_lep.Draw('same')
        canvas.cd(5)
        histo_ZH_had[fs][cat].Draw()
        plotZH_had.Draw('same')
        canvas.cd(6)
        histo_ZH_lep[fs][cat].Draw()
        plotZH_lep.Draw('same')
        canvas.cd(7)
        histo_ttH[fs][cat].Draw()
        plotttH.Draw('same')
        canvas.cd(8)
        histo_qqZZ[fs][cat].Draw()
        plotqqZZ.Draw('same')
        canvas.cd(9)
        histo_ggZZ[fs][cat].Draw()
        plotggZZ.Draw('same')
        canvas.cd(10)
        plotzjets.Draw()
        

        canvas.SaveAs(outputPath1+"/pdfs_%s_%i.png" % (sFinalStates[fs], cat))
        canvas.SaveAs(outputPath1+"/pdfs_%s_%i.pdf" % (sFinalStates[fs], cat))
        
        
        
        
        
        
