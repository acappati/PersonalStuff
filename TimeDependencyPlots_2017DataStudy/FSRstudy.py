#!/usr/bin/env python


import ROOT, helper, math
from ROOT import TFile, TH1, TH1F, TCanvas, gSystem, TLegend, TAttFill, TColor, kBlue, kRed, kWhite, kBlack



# *****************************
# Declare all the variables

redoHistos = True

# data tree options 
ZZTree   = False
CRZLTree = False
ZTree    = True

# data periods options
period = "data2016"
# period = "data2017"
# *****************************


#input file
if(period == "data2016"):
    #data = TFile.Open("/data3/Higgs/170222/AllData/ZZ4lAnalysis.root") #2016 data
    data = TFile.Open("/data3/Higgs/170907_2016/AllData/ZZ4lAnalysis2016.root")
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
    #data = TFile.Open("/data3/Higgs/170907_Data2017_add/AllData/ZZ4lAnalysis.root")
    data = TFile.Open("/data3/Higgs/170907_Data2017_add/AllData/ZZ4lAnalysis_add170909.root")
    lumi = 13.88   # fb-1
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
outputDir = "studyFSR_" + str(period) + "_" + str(treeText)
gSystem.Exec("mkdir -p " + outputDir)
print "Output directory created!"



if(redoHistos):

    TH1.SetDefaultSumw2() # set sumw2 = true fro all the histograms created from now on

    # define histos
    ZMass_data_hist       = TH1F( 'ZMass_data',       'ZMass_data',       100, 20, 120)  # ZMass post FSR
    ZMassPreFSR_data_hist = TH1F( 'ZMassPreFSR_data', 'ZMassPreFSR_data', 100, 20, 120)  # ZMass pre FSR

    if ZTree: 
        #read ttree ZTree
        print "reading tree", data.GetName(),treeText,tree.GetName()  ,"..."

        nFSR = 0
        for event in tree :
            if ( event.Zsel < 0 ) : continue # skip events that do not pass the trigger 

            
            if event.ZMassPreFSR > 0 : # fsrPt.size() > 0
                nFSR += 1
                ZMass_data_hist.Fill(event.ZMass)  # ZMass post FSR
                ZMassPreFSR_data_hist.Fill(event.ZMassPreFSR)  # ZMass pre FSR

    
    #save histograms in a root file 
    print "saving histograms into root file ..."
    histoFSRstudy = TFile.Open("histoFSRstudy_"+ period + "_" + treeText +".root", "RECREATE")
    histoFSRstudy.cd()
    ZMass_data_hist.Write()
    ZMassPreFSR_data_hist.Write()

    histoFSRstudy.Close()




# read histo from FSRstudy root file
histoINPUT = TFile.Open("histoFSRstudy_"+ period + "_" + treeText +".root")
print 'Reading file', histoINPUT.GetName(),'...'

hist_ZMass_postFSR_data = histoINPUT.Get('ZMass_data')
hist_ZMass_preFSR_data  = histoINPUT.Get('ZMassPreFSR_data')



# *** do plots ***
canvas = TCanvas("canvas","canvas",800,800)


hist_ZMass_postFSR_data.SetLineColor(kBlue)
hist_ZMass_postFSR_data.SetLineWidth(2)

hist_ZMass_preFSR_data.SetLineColor(kRed)
hist_ZMass_preFSR_data.SetLineWidth(2)

hist_ZMass_postFSR_data.SetStats(0)
hist_ZMass_postFSR_data.SetTitle("")

hist_ZMass_postFSR_data.GetXaxis().SetTitleSize(20)
hist_ZMass_postFSR_data.GetXaxis().SetTitleFont(43)
hist_ZMass_postFSR_data.GetXaxis().SetTitleOffset(1.8)
hist_ZMass_postFSR_data.GetXaxis().SetTitle("Mass [GeV/c^{2}]")
hist_ZMass_postFSR_data.GetXaxis().SetLabelFont(43)
hist_ZMass_postFSR_data.GetXaxis().SetLabelSize(15)
hist_ZMass_postFSR_data.GetYaxis().SetTitleSize(20)
hist_ZMass_postFSR_data.GetYaxis().SetTitleFont(43)
hist_ZMass_postFSR_data.GetYaxis().SetTitleOffset(1.8)
hist_ZMass_postFSR_data.GetYaxis().SetLabelFont(43)
hist_ZMass_postFSR_data.GetYaxis().SetLabelSize(15)
hist_ZMass_postFSR_data.GetYaxis().SetTitle("Events")

canvas.cd()
hist_ZMass_postFSR_data.Draw("hist")
hist_ZMass_preFSR_data.Draw("histsame")



legend = TLegend(0.77,0.88,0.97,0.96)
legend.AddEntry(hist_ZMass_preFSR_data, "without FSR","l")
legend.AddEntry(hist_ZMass_postFSR_data,"with FSR ","l")
legend.SetFillColor(kWhite)
legend.SetLineColor(kBlack)
legend.SetTextFont(43)
legend.SetTextSize(20)
legend.Draw()

canvas.Update()

canvas.SaveAs(outputDir + "/ZMass_FSRstudy.pdf")
canvas.SaveAs(outputDir + "/ZMass_FSRstudy.png")
canvas.SaveAs(outputDir + "/ZMass_FSRstudy.root")

print "plots done"
