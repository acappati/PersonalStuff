#!/usr/bin/env python

import ROOT, array, DoubleCB, CMSGraphics, random, copy
ROOT.gROOT.SetBatch(True)
from ROOT import TCanvas, TFile, TH1F, TF1, gSystem
gSystem.Load('libRooFit')
from ROOT import RooFit, RooRealVar, RooGaussian, RooDataSet, RooArgList, RooTreeData, RooDataHist, RooCBShape, RooNumConvPdf, RooFFTConvPdf, RooGlobalFunc
from ROOT import RooCmdArg, RooArgSet, kFALSE, RooLinkedList, kBlue, kRed, kBlack, kOpenStar, kWhite, kGray
from ROOT import gStyle, TGraph, TGraphErrors, TMath, TMultiGraph, TLine, gPad, TGaxis, TLegend, TText, TLatex, TColor
from array import array
from DoubleCB import DoubleCB
from CMSGraphics import makeCMSCanvas, makeLegend, printLumiPrelLeft, printLumiPrelOut, printLumiLeft, printLumiOut
import math


class Result:
    def __init__(self, mean, width, lumi, mean_err, width_err):
        self.mean = mean
        self.width = width
        self.lumi = lumi
        self.mean_err = mean_err
        self.width_err = width_err



def GraphVsLumi(result, outputDir, title):
    #gSystem.Exec("mkdir -p ZPlots")
    can1 = makeCMSCanvas(str(random.random()),"mean vs lumi ",900,700)
    can2 = makeCMSCanvas(str(random.random()),"width vs lumi ",900,700)
    lumi = []
    lumi_err = []
    mean = []
    mean_err = []
    width = []
    width_err = []
    sumLumi = 0.
    for i in range(0,len(result)):
        sumLumi += float(result[i].lumi)
        lumi.append(sumLumi - float(result[i].lumi)/2)
        lumi_err.append(float(result[i].lumi)/2)
        mean.append(result[i].mean)
        mean_err.append(result[i].mean_err)
        width.append(result[i].width)
        width_err.append(result[i].width_err)
    graph1 = TGraphErrors(len(result),array('d',lumi),array('d',mean),array('d',lumi_err),array('d',mean_err))
    graph2 = TGraphErrors(len(result),array('d',lumi),array('d',width),array('d',lumi_err),array('d',width_err))
    can1.cd()
    graph1.SetTitle("")
    graph1.GetXaxis().SetTitle("Lumi [fb^{-1}]")
    graph1.GetYaxis().SetTitle("Mass [GeV]")
    graph1.SetMarkerStyle(20)
    graph1.SetMarkerSize(1)
    graph1.Draw("AP")
    printLumiPrelOut(can1)
    can1.SaveAs(str(outputDir) + "/"+ title +"_mean.pdf")
    can1.SaveAs(str(outputDir) + "/"+ title +"_mean.png")
    can2.cd()
    graph2.SetTitle("")
    graph2.GetXaxis().SetTitle("Lumi [fb^{-1}]")
    graph2.GetYaxis().SetTitle("Width [GeV]")
    graph2.SetMarkerStyle(20)
    graph2.SetMarkerSize(1)
    graph2.Draw("AP")
    printLumiPrelOut(can2)
    can2.SaveAs(str(outputDir) + "/"+ title +"_width.pdf")
    can2.SaveAs(str(outputDir) + "/"+ title +"_width.png")
    return graph1, graph2;


def ZMultVsLumi(histo, recorded, outputDir, title):
    #gSystem.Exec("mkdir -p ZPlots")
    can1 = makeCMSCanvas(str(random.random()),"mult vs lumi ",900,700)
    lumi = []
    lumi_err = []
    mult = []
    mult_err = []
    sumLumi = 0.
    for i in range(0,len(recorded)):
        sumLumi += float(recorded[i])
        lumi.append(sumLumi - float(recorded[i])/2)
        lumi_err.append(float(recorded[i])/2)
        mult.append(histo[i].GetEntries()/float(recorded[i]))
        mult_err.append(math.sqrt(histo[i].GetEntries()/float(recorded[i])))
    graph1 = TGraphErrors(len(lumi),array('d',lumi),array('d',mult),array('d',lumi_err),array('d',mult_err))
    can1.cd()
    graph1.SetTitle("")
    graph1.GetXaxis().SetTitle("Lumi [fb^{-1}]")
    graph1.GetYaxis().SetTitle("#Z / fb^{-1}")
    graph1.SetMarkerStyle(20)
    graph1.SetMarkerSize(1)
    graph1.Draw("AP")
    printLumiPrelOut(can1)
    can1.SaveAs(str(outputDir) + "/Z_multiplicity_"+title+".pdf")
    can1.SaveAs(str(outputDir) + "/Z_multiplicity_"+title+".png")
    return graph1;


def Ratio(histo1,histo2, recorded, title):
    #gSystem.Exec("mkdir -p ZPlots")
    can1 = makeCMSCanvas(str(random.random()),"mult vs lumi ",900,700)
    lumi = []
    lumi_err = []
    ratio = []
    ratio_err = []
    sumLumi = 0.
    for i in range(0,len(recorded)):
        sumLumi += float(recorded[i])
        lumi.append(sumLumi - float(recorded[i])/2)
        lumi_err.append(float(recorded[i])/2)
        ratio.append(histo1[i].GetEntries()/histo2[i].GetEntries())
        ratio_err.append(0)
    graph1 = TGraphErrors(len(recorded),array('d',lumi),array('d',ratio),array('d',lumi_err),array('d',ratio_err))
    can1.cd()
    graph1.SetTitle("")
    graph1.GetXaxis().SetTitle("Lumi [fb^{-1}]")
    graph1.GetYaxis().SetTitle("e^{+}e^{-}/#mu^{+}#mu^{-}")
    graph1.SetMarkerStyle(20)
    graph1.SetMarkerSize(1)
    graph1.Draw("AP")
    printLumiPrelOut(can1)
    can1.SaveAs("ZPlots/Z_ratio_"+title+".pdf")
    can1.SaveAs("ZPlots/Z_ratio_"+title+".png")
    return;


def MeanRMSVsLumi(histo, recorded, outputDir, title):
    can1 = makeCMSCanvas(str(random.random()),"mean vs lumi ",900,700)
    can2 = makeCMSCanvas(str(random.random()),"RMS vs lumi ",900,700)
    lumi = []
    lumi_err = []
    mean = []
    mean_err = []
    RMS = []
    RMS_err = []
    sumLumi = 0.
    for i in range(0,len(recorded)):
        sumLumi += float(recorded[i])
        lumi.append(sumLumi - float(recorded[i])/2)
        lumi_err.append(float(recorded[i])/2)
        mean.append(histo[i].GetMean())
        mean_err.append(0.)#dont put error on ISO and SIP histo[i].GetRMS()
        RMS.append(histo[i].GetRMS())
        RMS_err.append(0.)
    graph1 = TGraphErrors(len(recorded),array('d',lumi),array('d',mean),array('d',lumi_err),array('d',mean_err))
    graph2 = TGraphErrors(len(recorded),array('d',lumi),array('d',RMS),array('d',lumi_err),array('d',RMS_err))
    can1.cd()
    graph1.SetTitle("")
    graph1.GetXaxis().SetTitle("Lumi [fb^{-1}]")
    graph1.GetYaxis().SetTitle(" ")
    graph1.SetMarkerStyle(20)
    graph1.SetMarkerSize(1)
    graph1.Draw("AP")
    printLumiPrelOut(can1)
    can1.SaveAs(str(outputDir) + "/"+title+"_mean.pdf")
    can1.SaveAs(str(outputDir) + "/"+title+"_mean.png")
    can2.cd()
    graph2.SetTitle("")
    graph2.GetXaxis().SetTitle(" ")
    graph2.GetYaxis().SetTitle("Width [GeV]")
    graph2.SetMarkerStyle(20)
    graph2.SetMarkerSize(1)
    graph2.Draw("AP")
    printLumiPrelOut(can2)
    can2.SaveAs(str(outputDir) + "/"+title+"_width.pdf")
    can2.SaveAs(str(outputDir) + "/"+title+"_width.png")
    
    return graph1,graph2;


def DoSimpleFit(histo, lumi, ZZtree, outputDir, title, fitMC, fitDATA):  #fitMC = true for fitting MC in FitMC.py, fitDATA = true for fitting data in FitDATA.py
    
    can = []
    result = []
    for i in range(0,len(histo)):
        canv = makeCMSCanvas(str(random.random())+str(i),"Fit result "+str(i),900,700)
        can.append(canv)

    Z1_fitFunction = TF1(str(random.random()),DoubleCB(),60. ,120., 7)

    Z1_fitFunction.SetParLimits(0,80.,100.)#mean
    Z1_fitFunction.SetParLimits(1,0.,5.)#width
    Z1_fitFunction.SetParLimits(2,0.,2.)#alpha1
    Z1_fitFunction.SetParLimits(3,0.,10.)#n1
    Z1_fitFunction.SetParLimits(4,0.,2.)#alpha2
    Z1_fitFunction.SetParLimits(5,0.,10.)#n2
#    Z1_fitFunction.SetParLimits(6,0.,1000000.)#const

    Z1_fitFunction.SetParName(0,"Mean")#mean
    Z1_fitFunction.SetParName(1,"Sigma")#width
    Z1_fitFunction.SetParName(2,"#alpha_{1}")#alpha1
    Z1_fitFunction.SetParName(3,"n_{1}")#n1
    Z1_fitFunction.SetParName(4,"#alpha_{2}")#alpha2
    Z1_fitFunction.SetParName(5,"n_{2}")#n2
    Z1_fitFunction.SetParName(6,"C")#const

    for i in range(0,len(histo)):        
#       print "***", histo[i].GetTitle(), lumi[i], title, histo[i].GetEntries(), "***"
        # define fit function 
        if(fitMC) :          # fit function for fitting MC in FitMC.py
            Z1_fitFunction.SetParameters(91.,2.5,1.,1.,1.,1.,histo[i].Integral()/5.)
        elif(fitDATA):       # fit function for fitting DATA in FitDATA.py
            Z1_fitFunction.SetParameters(91.,2.5,1.,1.,1.,1.,histo[i].Integral()/5.)
        else : 
            if(ZZtree) :
                Z1_fitFunction.SetParameters(91.,2.5,1.,1.,1.,1.,10.)
            else:
                Z1_fitFunction.SetParameters(91.,2.5,1.,1.,1.,1.,histo[i].Integral()/5.)

        # do the fit 
        can[i].cd()        
        histo[i].Fit(Z1_fitFunction) # do the fit
        histo[i].Draw("E")
        gStyle.SetOptFit()
        mass = Z1_fitFunction.GetParameters()[0]
        width = Z1_fitFunction.GetParameters()[1]
        mass_err = Z1_fitFunction.GetParError(0)
        width_err = Z1_fitFunction.GetParError(1)
        result.append(Result(mass,width,lumi[i],mass_err,width_err))
        print ("Mass: " +str(mass) + " Width: " + str(width))
        printLumiPrelOut(can[i])
        
        if(fitMC) :
            can[i].SaveAs(str(outputDir) + "/"+ str(title[i]) +"_fit.pdf")
            can[i].SaveAs(str(outputDir) + "/"+ str(title[i]) +"_fit.png")
        elif(fitDATA) :
            can[i].SaveAs(str(outputDir) + "/"+ str(title[i]) +"_fit.pdf")
            can[i].SaveAs(str(outputDir) + "/"+ str(title[i]) +"_fit.png")
        else :
            can[i].SaveAs(str(outputDir) + "/"+ title +"_fit_"+str(i)+".pdf")
            can[i].SaveAs(str(outputDir) + "/"+ title +"_fit_"+str(i)+".png")
            
                
    return result;


#to be fixed!
def DoDCBunbinnedFit(histo, lumi, ZZtree, outputDir, title, fitMC, fitDATA):  #fitMC = true for fitting MC in FitMC.py, fitDATA = true for fitting data in FitDATA.py

    # create canvas list
    can = []
    datahist = []
    massplot = []
    result = []
    for i in range(0,len(histo)):
        canv = makeCMSCanvas(str(random.random())+str(i),"Fit result "+str(i),900,700)
        can.append(canv)

    # variable
    if "ele" in title:
      x1_zmass = RooRealVar("x1_zmass","m_{e^{+}e^{-}}",60,120)
    if "mu" in title:
      x1_zmass = RooRealVar("x1_zmass","m_{#mu^{+}#mu^{-}}",60,120)

    # define fit function (DCB)
    meanDCB      = RooRealVar("mean_{DCB}",    "mean of DCB",   60.,120. )
    sigmaDCB     = RooRealVar("#sigma_{DCB}",  "width of DCB",  0., 6.   )
    alphaDCB     = RooRealVar("#alpha_{DCB}",  "alpha of DCB",  0., 100. )
    nDCB         = RooRealVar("n_{DCB}",       "n of DCB",      0., 10.  )
    alpha2DCB    = RooRealVar("#alpha2_{DCB}", "alpha2 of DCB", 0., 100. )
    n2DCB        = RooRealVar("n2_{DCB}",      "n2 of DCB",     0., 10.  )
    DCB_function = DoubleCB( x1_zmass, meanDCB, sigmaDCB, alphaDCB, nDCB, alpha2DCB, n2DCB)


    # do the fit and plot
    for i in range(0,len(histo)):
         
        datahist[i] = RooDataHist("data","data",RooArgList(x1_zmass),RooFit.Import(histo))
        DCB_function.fitTo(datahist[i],ROOT.RooFit.Extended(1)) #do the fit

        massplot[i] = x1_zmass.frame()
        datahist[i].plotOn(massplot[i])
        DCB_function.plotOn(massplot[i])
        DCB_function.paramOn(massplot[i], RooFit.Layout(0.65,0.99,0.9))

        can[i].cd()
        gStyle.SetOptFit()
        massplot[i].Draw()
        
       

        mass      = meanDCB.getVal()
        width     = sigmaDCB.getVal()
        mass_err  = 0. #to check
        width_err = 0. #to check
        result.append(Result(mass,width,lumi[i],mass_err,width_err))
        print ("Mass: " +str(mass) + " Width: " + str(width))
        printLumiPrelOut(can[i])
 
        if(fitMC) :
            can[i].SaveAs(str(outputDir) + "/"+ str(title[i]) +"_fit.pdf")
            can[i].SaveAs(str(outputDir) + "/"+ str(title[i]) +"_fit.png")
        elif(fitDATA) :
            can[i].SaveAs(str(outputDir) + "/"+ str(title[i]) +"_fit.pdf")
            can[i].SaveAs(str(outputDir) + "/"+ str(title[i]) +"_fit.png")
        else :
            can[i].SaveAs(str(outputDir) + "/"+ title +"_fit_"+str(i)+".pdf")
            can[i].SaveAs(str(outputDir) + "/"+ title +"_fit_"+str(i)+".png")
            
                
    return result;

    

def DoLandauFit(histo, lumi, ZZtree, title):
    #gSystem.Exec("mkdir -p FitResults")
    can = []
    result = []
    for i in range(0,len(histo)):
        canv = makeCMSCanvas(str(random.random())+str(i),"Fit result "+str(i),900,700)
        can.append(canv)
    
    Z1_fitFunction = TF1('Z1_fitFunction',"[0]*TMath::Landau(x,[1],[2])",0,4)
    
    Z1_fitFunction.SetParameters(100,0.1,1.0)
    
    #Z1_fitFunction.SetParLimits(0,85.,95.)
    #Z1_fitFunction.SetParLimits(1,0.,5.)
    #Z1_fitFunction.SetParLimits(2,0.,2.)
    
    Z1_fitFunction.SetParName(0,"Const")#constant
    Z1_fitFunction.SetParName(1,"Position")#position
    Z1_fitFunction.SetParName(2,"Scale")#scale
    
    for i in range(0,len(histo)):
        can[i].cd()
        histo[i].Fit(Z1_fitFunction,"","",0.,2.)
        histo[i].Draw()
        gStyle.SetOptFit()
        position = Z1_fitFunction.GetParameters()[1]
        scale = Z1_fitFunction.GetParameters()[2]
        position_err = Z1_fitFunction.GetParError(1)
        scale_err = Z1_fitFunction.GetParError(2)
        result.append(Result(position,scale,lumi[i],position_err,scale_err))
        print ("Position: " +str(position) + " Scale: " + str(scale))
        printLumiPrelOut(can[i])
        can[i].SaveAs("FitResults/"+title+"_Landfit_"+str(i)+".pdf")
        can[i].SaveAs("FitResults/"+title+"_Landfit_"+str(i)+".png")
    return result;


def SAME3VsLumi(g1, g2, g3, title, ptype, lineMC1, lineDATA1, lineMC2, lineDATA2, lineMC3, lineDATA3, DoInclusive, dataPeriod):
    canvas = makeCMSCanvas(str(random.random()),"canvas",900,700)
    canvas.cd()
    graph2=copy.deepcopy(g2)
    graph3=copy.deepcopy(g3)
    graph2.SetMarkerColor(kBlue)#electrons
    graph3.SetMarkerColor(kRed)#muons
    multigraph = TMultiGraph()
    if(DoInclusive):
        graph1=copy.deepcopy(g1)
        multigraph.Add(graph1,"AP")
    multigraph.Add(graph2,"AP")
    multigraph.Add(graph3,"AP")
    multigraph.Draw("AP")
    TGaxis.SetMaxDigits(2)
    TGaxis.SetExponentOffset(-0.05, 0.02, "y")
    multigraph.GetXaxis().SetTitle("L [fb^{-1}]")
    multigraph.GetYaxis().SetTitleOffset(1.4)
    if(ptype == "Zmass"):
        multigraph.GetYaxis().SetTitle("M_{Z} [GeV]")
        # multigraph.GetYaxis().SetTitle("M_{l^{+}l^{-}} [GeV]")
        multigraph.SetMaximum(max(multigraph.GetHistogram().GetMaximum(),91.4))
        multigraph.SetMinimum(min(multigraph.GetHistogram().GetMinimum(),89.6))
    elif(ptype == "Zwidth"):
        multigraph.GetYaxis().SetTitle("#Gamma_{Z} [GeV]")
    elif(ptype == "Zmult"):
        multigraph.GetYaxis().SetTitle("#Z / fb^{-1}")
        if(not DoInclusive) :
            multigraph.SetMaximum(max(multigraph.GetHistogram().GetMaximum(),60000.)) # set y axis minimum at 60000.
            multigraph.SetMinimum(0.)     # set y axis minimum at 0. 
            # multigraph.SetMaximum(60000.)  #second type: vs 2016 plots 
            # multigraph.SetMinimum(25000.)      
    printLumiPrelOut(canvas)

    
    # Draw legend 
    legend = TLegend(0.93,0.84,0.99,0.93)
    if(DoInclusive):
        #legend.AddEntry(graph1,"inclusive","P")
        legend.AddEntry(graph1,"BB","P")
        legend.AddEntry(graph2,"BE","P")
        legend.AddEntry(graph3,"EE","P")
    else :
        legend.AddEntry(graph2,"e^{+}e^{-}","P")
        legend.AddEntry(graph3,"#mu^{+}#mu^{-}","P")
    legend.SetFillColor(kWhite)
    legend.SetLineColor(kBlack)
    legend.SetTextFont(43)
    legend.SetTextSize(20)
    legend.Draw()
    canvas.Update()

    
    # Draw letters for data-taking periods
    if(dataPeriod == "data2017"):
        textLetters = TLatex()
        textLetters.SetTextColor(kGray+1)
        textLetters.SetTextSize(0.03)
        if(ptype == "Zmass"):
            textLetters.DrawLatex(2.,  gPad.GetUymin()+0.2,"B")
            textLetters.DrawLatex(9.5, gPad.GetUymin()+0.2,"C")
            textLetters.DrawLatex(16., gPad.GetUymin()+0.2,"D")
            textLetters.DrawLatex(23., gPad.GetUymin()+0.2,"E")
            textLetters.DrawLatex(36., gPad.GetUymin()+0.2,"F")
        elif(ptype == "Zwidth"):
            textLetters.DrawLatex(2.,  gPad.GetUymin()+0.3,"B")
            textLetters.DrawLatex(9.5, gPad.GetUymin()+0.3,"C")
            textLetters.DrawLatex(16., gPad.GetUymin()+0.3,"D")
            textLetters.DrawLatex(23., gPad.GetUymin()+0.3,"E")
            textLetters.DrawLatex(36., gPad.GetUymin()+0.3,"F")
        elif(ptype == "Zmult") :
            textLetters.DrawLatex(2.,  260000,"B")
            textLetters.DrawLatex(9.5, 260000,"C")
            textLetters.DrawLatex(16., 260000,"D")
            textLetters.DrawLatex(23., 260000,"E")
            textLetters.DrawLatex(36., 260000,"F")

    if(dataPeriod == "data2018"):
        textLetters = TLatex()
        textLetters.SetTextColor(kGray+1)
        textLetters.SetTextSize(0.03)
        if(ptype == "Zmass"):
            textLetters.DrawLatex(2.,  gPad.GetUymin() + 0.6,"Av1")
            textLetters.DrawLatex(9.5, gPad.GetUymin() + 0.6,"Av2")
            textLetters.DrawLatex(12.5, gPad.GetUymin()+ 0.6,"Av3")
            textLetters.DrawLatex(15., gPad.GetUymin() + 0.6,"Bv1")
        elif(ptype == "Zwidth"):
            textLetters.DrawLatex(2.,  gPad.GetUymin() +0.3,"Av1")
            textLetters.DrawLatex(9.5, gPad.GetUymin() +0.3,"Av2")
            textLetters.DrawLatex(12.5, gPad.GetUymin()+0.3,"Av3")
            textLetters.DrawLatex(15., gPad.GetUymin() +0.3,"Bv1")
        elif(ptype == "Zmult") :
            textLetters.DrawLatex(2.,  260000,"Av1")
            textLetters.DrawLatex(9.5, 260000,"Av2")
            textLetters.DrawLatex(12.5, 260000,"Av3")
            textLetters.DrawLatex(15., 260000,"Bv1")

    
    # ****
    if(dataPeriod == "data2018"):
        # draw vertical lines that divide different data taking periods
        down    = gPad.GetUymin()
        up      = gPad.GetUymax()
        
        lineAv1 = TLine(7.643, down, 7.643, up) # Run2018Av1 up to 7.643 fb-1
        lineAv1.SetLineColor(kBlack)
        lineAv1.SetLineStyle(2)
        lineAv1.Draw()    
        
        lineAv2 = TLine(12.441, down, 12.441, up) # Run2018Av2 up to 12.441 fb-1
        lineAv2.SetLineColor(kBlack)
        lineAv2.SetLineStyle(2)
        lineAv2.Draw()
        
        lineAv3 = TLine(12.863, down, 12.863, up) # Run2018Av3 up to 12.863 fb-1
        lineAv3.SetLineColor(kBlack)
        lineAv3.SetLineStyle(2)
        lineAv3.Draw()
        

    if(dataPeriod == "data2017"):
        # draw vertical lines that divide different data taking periods
        down    = gPad.GetUymin()
        up      = gPad.GetUymax()
        
        lineB = TLine(4.793, down, 4.793, up) # Run2017B up to 4.793 fb-1
        lineB.SetLineColor(kBlack)
        lineB.SetLineStyle(2)
        lineB.Draw()    
        
        lineC = TLine(14.549, down, 14.549, up) # Run2017C up to 14.549 fb-1
        lineC.SetLineColor(kBlack)
        lineC.SetLineStyle(2)
        lineC.Draw()
        
        lineD = TLine(18.868, down, 18.868, up) # Run2017D up to 18.868 fb-1
        lineD.SetLineColor(kBlack)
        lineD.SetLineStyle(2)
        lineD.Draw()
        
        lineE = TLine(28.293, down, 28.293, up) # Run2017E up to 28.293 fb-1
        lineE.SetLineColor(kBlack)
        lineE.SetLineStyle(2)
        lineE.Draw()

    if(dataPeriod == "data2016"):
        # draw vertical lines that divide different data taking periods
        down    = gPad.GetUymin()
        up      = gPad.GetUymax()
        
        lineB = TLine(5.789, down, 5.789, up) # Run2016B up to 5.789 fb-1
        lineB.SetLineColor(kBlack)
        lineB.SetLineStyle(2)
        lineB.Draw()
        
        lineC = TLine(8.366, down, 8.366, up) # Run2016C up to 8.366 fb-1
        lineC.SetLineColor(kBlack)
        lineC.SetLineStyle(2)
        lineC.Draw() 
        
        lineD = TLine(12.616, down, 12.616, up) # Run2016D up to 12.616 fb-1
        lineD.SetLineColor(kBlack)
        lineD.SetLineStyle(2)
        lineD.Draw()    
        
        lineE = TLine(16.624, down, 16.624, up) # Run2016E up to 16.624 fb-1
        lineE.SetLineColor(kBlack)
        lineE.SetLineStyle(2)
        lineE.Draw()    
        
        lineF = TLine(19.725, down, 19.725, up) # Run2016F up to 19.725 fb-1
        lineF.SetLineColor(kBlack)
        lineF.SetLineStyle(2)
        lineF.Draw()    
        
        lineG = TLine(27.268, down, 27.268, up) # Run2016G up to 27.268 fb-1
        lineG.SetLineColor(kBlack)
        lineG.SetLineStyle(2)
        lineG.Draw()       
    # ****
    
    # draw orizontal lines for MC and DATA fit
    if(ptype == "Zmass" or ptype == "Zwidth") :

        leftEnd  = gPad.GetUxmin()
        rightEnd = gPad.GetUxmax()

        if(DoInclusive):
            line1 = TLine(leftEnd,lineMC1,rightEnd,lineMC1)
            line1.SetLineColor(kBlack)
            line1.SetLineStyle(1) 
            line1.Draw()

            line2 = TLine(leftEnd,lineDATA1,rightEnd,lineDATA1)
            line2.SetLineColor(kBlack)
            line2.SetLineStyle(2)
            line2.Draw()
        
        # line for graph 2: color blue
        line3 = TLine(leftEnd,lineMC2,rightEnd,lineMC2)
        line3.SetLineColor(kBlue)
        line3.SetLineStyle(1) 
        line3.Draw()

        line4 = TLine(leftEnd,lineDATA2,rightEnd,lineDATA2)
        line4.SetLineColor(kBlue)
        line4.SetLineStyle(2)
        line4.Draw()
    
        # line for graph 3: color red
        line5 = TLine(leftEnd,lineMC3,rightEnd,lineMC3)
        line5.SetLineColor(kRed)
        line5.SetLineStyle(1) 
        line5.Draw()

        line6 = TLine(leftEnd,lineDATA3,rightEnd,lineDATA3)
        line6.SetLineColor(kRed)
        line6.SetLineStyle(2)
        line6.Draw()
    # ***    

    canvas.SaveAs(title + ".root")
    canvas.SaveAs(title + ".pdf")
    canvas.SaveAs(title + ".png")
    return;


def SAME2VsLumi(g1, g2,title, ptype, dataPeriod):
    canvas = makeCMSCanvas(str(random.random()),"canvas",900,700)
    canvas.cd()
    graph1=copy.deepcopy(g1)
    graph2=copy.deepcopy(g2)
    graph1.SetMarkerColor(kBlue)#electrons
    graph2.SetMarkerColor(kRed)#muons
    multigraph = TMultiGraph()
    multigraph.Add(graph1,"AP")
    multigraph.Add(graph2,"AP")
    multigraph.Draw("AP")
    TGaxis.SetMaxDigits(2)
    TGaxis.SetExponentOffset(-0.06, 0.02, "y")
    multigraph.GetXaxis().SetTitle("L [fb^{-1}]")
    multigraph.GetYaxis().SetTitleOffset(1.4)
    if(ptype == "ISO"):
        multigraph.GetYaxis().SetTitle("Isolation")
        gPad.Modified()
        multigraph.SetMinimum(0.5*gPad.GetUymin())
        multigraph.SetMaximum(1.5*gPad.GetUymax())
    elif(ptype == "SIP"):
        multigraph.GetYaxis().SetTitle("SIP")
        multigraph.SetMinimum(min(multigraph.GetHistogram().GetMinimum(),1.))
        multigraph.SetMaximum(max(multigraph.GetHistogram().GetMinimum(),2.2))

        min(multigraph.GetHistogram().GetMinimum(),1.)
    printLumiPrelOut(canvas)

    
    # Draw legend 
    legend = TLegend(0.93,0.84,0.99,0.93)
    legend.AddEntry(graph1,"e^{+}e^{-}","P")
    legend.AddEntry(graph2,"#mu^{+}#mu^{-}","P")
    legend.SetFillColor(kWhite)
    legend.SetLineColor(kBlack)
    legend.SetTextFont(43)
    legend.SetTextSize(20)
    legend.Draw()
    canvas.Update()

    

    # Draw letters for data-taking periods
    if(dataPeriod == "data2017"):
        textLetters = TLatex()
        textLetters.SetTextColor(kGray+1)
        textLetters.SetTextSize(0.03)
        if(ptype == "ISO"):
            textLetters.DrawLatex(2.,  0.8*gPad.GetUymax(),"B")
            textLetters.DrawLatex(9.5, 0.8*gPad.GetUymax(),"C")
            textLetters.DrawLatex(16., 0.8*gPad.GetUymax(),"D")
            textLetters.DrawLatex(23., 0.8*gPad.GetUymax(),"E")
            textLetters.DrawLatex(36., 0.8*gPad.GetUymax(),"F")
        elif(ptype == "SIP"):
            textLetters.DrawLatex(2.,  1.5,"B")    
            textLetters.DrawLatex(9.5, 1.5,"C")
            textLetters.DrawLatex(16., 1.5,"D")
            textLetters.DrawLatex(23., 1.5,"E")
            textLetters.DrawLatex(36., 1.5,"F")

    if(dataPeriod == "data2018"):
        textLetters = TLatex()
        textLetters.SetTextColor(kGray+1)
        textLetters.SetTextSize(0.03)
        if(ptype == "ISO"):
            textLetters.DrawLatex(2.,  0.8*gPad.GetUymax(),"Av1")
            textLetters.DrawLatex(9.5, 0.8*gPad.GetUymax(),"Av2")
            textLetters.DrawLatex(12.5, 0.8*gPad.GetUymax(),"Av3")
            textLetters.DrawLatex(15., 0.8*gPad.GetUymax(),"Bv1")
        elif(ptype == "SIP"):
            textLetters.DrawLatex(2.,  1.5,"Av1")    
            textLetters.DrawLatex(9.5, 1.5,"Av2")
            textLetters.DrawLatex(12.5, 1.5,"Av3")
            textLetters.DrawLatex(15., 1.5,"Bv1")
         


    # ****
    if(dataPeriod == "data2018"):
        # draw vertical lines that divide different data taking periods
        down    = gPad.GetUymin()
        up      = gPad.GetUymax()
        
        lineAv1 = TLine(7.643, down, 7.643, up) # Run2018Av1 up to 7.643 fb-1
        lineAv1.SetLineColor(kBlack)
        lineAv1.SetLineStyle(2)
        lineAv1.Draw()
        
        lineAv2 = TLine(12.441, down, 12.441, up) # Run2018Av2 up to 12.441 fb-1
        lineAv2.SetLineColor(kBlack)
        lineAv2.SetLineStyle(2)
        lineAv2.Draw()
        
        lineAv3 = TLine(12.863, down, 12.863, up) # Run2018Av3 up to 12.863 fb-1
        lineAv3.SetLineColor(kBlack)
        lineAv3.SetLineStyle(2)
        lineAv3.Draw()
        
    # ****
    if(dataPeriod == "data2017"):
        # draw vertical lines that divide different data taking periods
        down    = gPad.GetUymin()
        up      = gPad.GetUymax()
        
        lineB = TLine(4.793, down, 4.793, up) # Run2017B up to 4.793 fb-1
        lineB.SetLineColor(kBlack)
        lineB.SetLineStyle(2)
        lineB.Draw()
        
        lineC = TLine(14.549, down, 14.549, up) # Run2017C up to 14.549 fb-1
        lineC.SetLineColor(kBlack)
        lineC.SetLineStyle(2)
        lineC.Draw()
        
        lineD = TLine(18.868, down, 18.868, up) # Run2017D up to 18.868 fb-1
        lineD.SetLineColor(kBlack)
        lineD.SetLineStyle(2)
        lineD.Draw()
        
        lineE = TLine(28.293, down, 28.293, up) # Run2017E up to 28.293 fb-1
        lineE.SetLineColor(kBlack)
        lineE.SetLineStyle(2)
        lineE.Draw()
    
    # ****
    if(dataPeriod == "data2016"):
        # draw vertical lines that divide different data taking periods
        down    = gPad.GetUymin()
        up      = gPad.GetUymax()
        
        lineB = TLine(5.789, down, 5.789, up) # Run2016B up to 5.789 fb-1
        lineB.SetLineColor(kBlack)
        lineB.SetLineStyle(2)
        lineB.Draw()
        
        lineC = TLine(8.366, down, 8.366, up) # Run2016C up to 8.366 fb-1
        lineC.SetLineColor(kBlack)
        lineC.SetLineStyle(2)
        lineC.Draw() 
        
        lineD = TLine(12.616, down, 12.616, up) # Run2016D up to 12.616 fb-1
        lineD.SetLineColor(kBlack)
        lineD.SetLineStyle(2)
        lineD.Draw()    
        
        lineE = TLine(16.624, down, 16.624, up) # Run2016E up to 16.624 fb-1
        lineE.SetLineColor(kBlack)
        lineE.SetLineStyle(2)
        lineE.Draw()    
        
        lineF = TLine(19.725, down, 19.725, up) # Run2016F up to 19.725 fb-1
        lineF.SetLineColor(kBlack)
        lineF.SetLineStyle(2)
        lineF.Draw()    
        
        lineG = TLine(27.268, down, 27.268, up) # Run2016G up to 27.268 fb-1
        lineG.SetLineColor(kBlack)
        lineG.SetLineStyle(2)
        lineG.Draw()       
    # ****
    

    canvas.SaveAs(title + ".root")
    canvas.SaveAs(title + ".pdf")
    canvas.SaveAs(title + ".png")
    return;


def TwoFileSAME3VsLumi(F1graph1, F1graph2, F1graph3, F2graph1, F2graph2, F2graph3, title, type, DoInclusive):
    canvas = makeCMSCanvas(str(random.random()),"canvas",900,700)
    canvas.cd()
    F1graph1.SetMarkerColor(kBlack)#electrons
    F1graph2.SetMarkerColor(kBlue)#electrons
    F1graph3.SetMarkerColor(kRed)#muons
    F2graph1.SetMarkerColor(kBlack)#electrons
    F2graph2.SetMarkerColor(kBlue)#electrons
    F2graph3.SetMarkerColor(kRed)#muons
    F2graph1.SetMarkerStyle(kOpenStar)#inclusive
    F2graph2.SetMarkerStyle(kOpenStar)#electrons
    F2graph3.SetMarkerStyle(kOpenStar)#muons
    multigraph = TMultiGraph()
    if(DoInclusive or type == "ZmassBARELL" or type == "ZwidthBARELL" or type == "ZmultBARELL"):
        multigraph.Add(F1graph1,"AP")
        multigraph.Add(F2graph1,"AP")
    multigraph.Add(F1graph2,"AP")
    multigraph.Add(F1graph3,"AP")
    multigraph.Add(F2graph2,"AP")
    multigraph.Add(F2graph3,"AP")
    multigraph.Draw("AP")
    multigraph.GetXaxis().SetTitle("L [fb^{-1}]")
    multigraph.GetXaxis().SetLimits(0.,40.)
    multigraph.GetYaxis().SetTitleOffset(1.4)
    gPad.Modified()

    if(type == "Zmass" or type == "ZmassBARELL"):
        multigraph.SetMaximum(1.008*gPad.GetUymax())
        multigraph.GetYaxis().SetTitle("M_{l^{+}l^{-}} [GeV]")
    elif(type == "Zwidth" or type == "ZwidthBARELL"):
        multigraph.SetMaximum(1.1*gPad.GetUymax())
        multigraph.GetYaxis().SetTitle("#Gamma [GeV]")
    elif(type == "Zmult" or type == "ZmultBARELL"):
        multigraph.SetMaximum(1.2*gPad.GetUymax())
        multigraph.GetYaxis().SetTitle("#Z")
    elif(type == "SIP"):
        multigraph.GetYaxis().SetTitle("SIP")
        multigraph.SetMaximum(1.1*gPad.GetUymax())
    printLumiPrelOut(canvas)
    canvas.Update()
    down = gPad.GetUymin()
    up = gPad.GetUymax()
    MC_muMass = 90.92
    MC_eleMass = 90.63
    Data_muMass = 90.95
    Data_eleMass = 90.68
    lineB = TLine(5.57,down,5.57, up)
    lineB.SetLineColor(kBlack)
    lineB.SetLineStyle(2)
    lineB.Draw()
    lineC = TLine(8.58,down,8.58,up)
    lineC.SetLineColor(kBlack)
    lineC.SetLineStyle(2)
    lineC.Draw()
    lineD = TLine(12.9,down,12.9,up)
    lineD.SetLineColor(kBlack)
    lineD.SetLineStyle(2)
    lineD.Draw()
    lineE = TLine(16.57,down,16.57,up)
    lineE.SetLineColor(kBlack)
    lineE.SetLineStyle(2)
    lineE.Draw()
    lineF = TLine(19.7,down,19.7,up)
    lineF.SetLineColor(kBlack)
    lineF.SetLineStyle(2)
    lineF.Draw()
    lineG = TLine(26.9,down,26.9,up)
    lineG.SetLineColor(kBlack)
    lineG.SetLineStyle(2)
    lineG.Draw()
    if(type == "Zmass"):
        lineMC_ele = TLine(0.,MC_eleMass,40., MC_eleMass)
        lineMC_ele.SetLineColor(kBlue)
        lineMC_ele.SetLineStyle(1)
        lineMC_ele.Draw()
        lineMC_mu = TLine(0.,MC_muMass,40., MC_muMass)
        lineMC_mu.SetLineColor(kRed)
        lineMC_mu.SetLineStyle(1)
        lineMC_mu.Draw()
        lineData_ele = TLine(0.,Data_eleMass,40., Data_eleMass)
        lineData_ele.SetLineColor(kBlue)
        lineData_ele.SetLineStyle(2)
        lineData_ele.Draw()
        lineData_mu = TLine(0.,Data_muMass,40., Data_muMass)
        lineData_mu.SetLineColor(kRed)
        lineData_mu.SetLineStyle(2)
        lineData_mu.Draw()
    legend = TLegend(0.80,0.75,0.965,0.93)
    if(type == "Zmass" or type == "Zwidth" or type == "Zmult"):
        legend.AddEntry(F1graph2,"e^{+}e^{-}","P")
        legend.AddEntry(F1graph3,"#mu^{+}#mu^{-}","P")
        legend.AddEntry(F2graph2,"e^{+}e^{-} ICHEP","P")
        legend.AddEntry(F2graph3,"#mu^{+}#mu^{-} ICHEP","P")
    else:
        legend.AddEntry(F1graph1,"EBEB","P")
        legend.AddEntry(F2graph1,"EBEB ICHEP","P")
        legend.AddEntry(F1graph2,"EBEE","P")
        legend.AddEntry(F1graph3,"EEEE","P")
        legend.AddEntry(F2graph2,"EBEE ICHEP","P")
        legend.AddEntry(F2graph3,"EEEE ICHEP","P")
    legend.SetTextFont(32)
    legend.Draw()
    canvas.SaveAs(title + ".pdf")
    canvas.SaveAs(title + ".png")
    return;


def TwoFileSAME2VsLumi(F1graph1, F1graph2, F2graph1, F2graph2, title, type):
    canvas = makeCMSCanvas(str(random.random()),"canvas",900,700)
    canvas.cd()
    F1graph1.SetMarkerColor(kBlue)#electrons
    F1graph2.SetMarkerColor(kRed)#muons
    F2graph1.SetMarkerColor(kBlue)#electrons
    F2graph2.SetMarkerColor(kRed)#muons
    F2graph1.SetMarkerStyle(kOpenStar)
    F2graph2.SetMarkerStyle(kOpenStar)
    multigraph = TMultiGraph()
    multigraph.Add(F1graph1,"AP")
    multigraph.Add(F1graph2,"AP")
    multigraph.Add(F2graph1,"AP")
    multigraph.Add(F2graph2,"AP")
    multigraph.Draw("AP")
    multigraph.GetXaxis().SetLimits(0.,40.)
    #TGaxis.SetMaxDigits(2)
    #TGaxis.SetExponentOffset(-0.06, 0.02, "y")
    multigraph.GetXaxis().SetTitle("L [fb^{-1}]")
    multigraph.GetYaxis().SetTitleOffset(1.4)
    if(type == "ISO"):
        multigraph.GetYaxis().SetTitle("Isolation")
        gPad.Modified()
        multigraph.SetMinimum(0.5*gPad.GetUymin())
        multigraph.SetMaximum(1.5*gPad.GetUymax())
    elif(type == "SIP"):
        multigraph.GetYaxis().SetTitle("SIP")
    printLumiPrelOut(canvas)
    canvas.Update()
    down = gPad.GetUymin()
    up = gPad.GetUymax()
    lineB = TLine(5.57,down,5.57, up)
    lineB.SetLineColor(kBlack)
    lineB.SetLineStyle(2)
    lineB.Draw()
    lineC = TLine(8.58,down,8.58,up)
    lineC.SetLineColor(kBlack)
    lineC.SetLineStyle(2)
    lineC.Draw()
    lineD = TLine(12.9,down,12.9,up)
    lineD.SetLineColor(kBlack)
    lineD.SetLineStyle(2)
    lineD.Draw()
    lineE = TLine(16.57,down,16.57,up)
    lineE.SetLineColor(kBlack)
    lineE.SetLineStyle(2)
    lineE.Draw()
    lineF = TLine(19.7,down,19.7,up)
    lineF.SetLineColor(kBlack)
    lineF.SetLineStyle(2)
    lineF.Draw()
    lineG = TLine(26.9,down,26.9,up)
    lineG.SetLineColor(kBlack)
    lineG.SetLineStyle(2)
    lineG.Draw()
    legend = TLegend(0.80,0.75,0.965,0.93)
    legend.AddEntry(F1graph1,"e^{+}e^{-}","P")
    legend.AddEntry(F1graph2,"#mu^{+}#mu^{-}","P")
    legend.AddEntry(F2graph1,"e^{+}e^{-} ICHEP","P")
    legend.AddEntry(F2graph2,"#mu^{+}#mu^{-} ICHEP","P")
    legend.SetTextFont(32)
    legend.Draw()

    canvas.SaveAs(title + ".pdf")
    canvas.SaveAs(title + ".png")
    return;


def ReadJSON(inputTXT, RunNum_B, LumiNum_B, RunNum_E, LumiNum_E, delivered):
    with open(inputTXT, "r") as Inp:
        print("Opening " + inputTXT + ".")
        for line in Inp:
            cleared_line = line.split(" ")
            RunNum_B.append(cleared_line[0])
            LumiNum_B.append(cleared_line[1])
            RunNum_E.append(cleared_line[2])
            LumiNum_E.append(cleared_line[3])
            delivered.append((cleared_line[4]))
    return;


def PlotNpv(inputTXT, lumiInp):
    canvas = makeCMSCanvas(str(random.random()),"Npv",900,700)
    canvas.Divide(1,2)
    canvas.cd(2)
    Npv = []
    error = []
    lumi = []
    lumi_err = []
    noerr = []
    with open(inputTXT, "r") as Inp:
        print("Opening " + inputTXT + ".")
        for line in Inp:
            cleared_line = line.split(" ")
            Npv.append(float(cleared_line[0]))
            error.append(float(cleared_line[1]))
    sumLumi = 0
    for i in range(0,len(lumiInp)):
        sumLumi += float(lumiInp[i])
        lumi.append(sumLumi - float(lumiInp[i])/2)
        lumi_err.append(float(lumiInp[i])/2)
        noerr.append(0.)

    graph1 = TGraphErrors(len(lumiInp),array('d',lumi),array('d',Npv),array('d',lumi_err),array('d',error))
    graph2 = TGraphErrors(len(lumiInp),array('d',lumi),array('d',Npv),array('d',lumi_err),array('d',noerr))
    graph1.SetTitle("")
    graph1.GetXaxis().SetTitle("Lumi [fb^{-1}]")
    graph1.GetXaxis().SetLimits(0.,40.)
    graph1.GetYaxis().SetTitle("Npv")
    graph1.SetMarkerStyle(20)
    graph1.SetMarkerSize(1)
    graph1.Draw("AP")
    canvas.cd(1)
    graph2.SetTitle("")
    graph2.GetXaxis().SetTitle("Lumi [fb^{-1}]")
    graph2.GetXaxis().SetLimits(0.,35.)
    graph2.GetYaxis().SetTitle("Npv")
    graph2.SetMarkerStyle(20)
    graph2.SetMarkerSize(1)
    graph2.Draw("AP")
    canvas.Update()
    canvas.cd(2)
    down = gPad.GetUymin()
    up = gPad.GetUymax()
    lineB = TLine(5.57,down,5.57, up)
    lineB.SetLineColor(kBlack)
    lineB.SetLineStyle(2)
    lineB.Draw()
    lineC = TLine(8.58,down,8.58,up)
    lineC.SetLineColor(kBlack)
    lineC.SetLineStyle(2)
    lineC.Draw()
    lineD = TLine(12.9,down,12.9,up)
    lineD.SetLineColor(kBlack)
    lineD.SetLineStyle(2)
    lineD.Draw()
    lineE = TLine(16.57,down,16.57,up)
    lineE.SetLineColor(kBlack)
    lineE.SetLineStyle(2)
    lineE.Draw()
    lineF = TLine(19.7,down,19.7,up)
    lineF.SetLineColor(kBlack)
    lineF.SetLineStyle(2)
    lineF.Draw()
    lineG = TLine(26.9,down,26.9,up)
    lineG.SetLineColor(kBlack)
    lineG.SetLineStyle(2)
    lineG.Draw()
    canvas.SaveAs("Npv.pdf")
    return;


def FillHisto(histoZ1Mass, tree):
    for event in tree:
        histoZ1Mass.Fill(event.Z1Mass)
        print(str(event.RunNumber) + " " + str(event.LumiNumber))
    return;


def DoCBFit(histo, ZZtree, title, p):
    canv = makeCMSCanvas(str(random.random()),"Fit result ",900,700)

    Z1_fitFunction = TF1('Z1_fitFunction',DoubleCB(),60. ,120., 7)

    Z1_fitFunction.SetParameters(p[0],p[1],p[2],p[3],p[4],p[5],p[6])

    Z1_fitFunction.SetParLimits(0,85.,95.)#mean
    Z1_fitFunction.SetParLimits(1,0.,5.)#width
    Z1_fitFunction.SetParLimits(2,0.,2.)#alpha1
    Z1_fitFunction.SetParLimits(3,0.,100.)#n1
    Z1_fitFunction.SetParLimits(4,0.,2.)#alpha2
    Z1_fitFunction.SetParLimits(5,0.,100.)#n2
    Z1_fitFunction.SetParLimits(6,0.,1000000.)#const

    Z1_fitFunction.SetParName(0,"Mean")#mean
    Z1_fitFunction.SetParName(1,"Sigma")#width
    Z1_fitFunction.SetParName(2,"#alpha_{1}")#alpha1
    Z1_fitFunction.SetParName(3,"n_{1}")#n1
    Z1_fitFunction.SetParName(4,"#alpha_{2}")#alpha2
    Z1_fitFunction.SetParName(5,"n_{2}")#n2
    Z1_fitFunction.SetParName(6,"C")#const

    canv.cd()
    histo.Fit(Z1_fitFunction)
    gStyle.SetOptFit()
    histo.SetTitle("")
    histo.Draw()
    print ("Mass: " +str(Z1_fitFunction.GetParameters()[0]) + " Width: " + str(Z1_fitFunction.GetParameters()[1]))
    printLumiPrelOut(canv)
    canv.SaveAs("DataVsMC/FitResults/"+title+"_fit.pdf")
    canv.SaveAs("DataVsMC/FitResults/"+title+"_fit.png")
    return Z1_fitFunction.GetParameters()[0], Z1_fitFunction.GetParError(0);


def DoRooFit(histo, title):
    can = makeCMSCanvas(str(random.random()),"Fit result ",900,700)
    
    #Varible
    if "ele" in title:
      x1 = RooRealVar("x1","m_{e^{+}e^{-}}",80,100)
    if "mu" in title:
      x1 = RooRealVar("x1","m_{#mu^{+}#mu^{-}}",80,100)

    #Define CB function
    m = RooRealVar("mean_{CB}","mean of gaussian",60,120)
    s = RooRealVar("#sigma_{CB}","width of gaussian",0,3)
    a = RooRealVar("#alpha_{CB}","mean of gaussian",0,100)
    n = RooRealVar("n_{CB}","width of gaussian",0,5)
    CB = RooCBShape("CB","CB PDF",x1, m, s, a, n)
    
    m.setConstant(kFALSE)
    s.setConstant(kFALSE)
    a.setConstant(kFALSE)
    n.setConstant(kFALSE)
    
    
    #Define Gaussian function
    mean1 = RooRealVar("mean_{G}","mean of gaussian",-60,60)
    sigma1 = RooRealVar("#sigma_{G}","width of gaussian",0,10)
    gauss1 = RooGaussian("gauss1","gaussian PDF",x1,mean1,sigma1)
    
    mean1.setConstant(kFALSE)
    sigma1.setConstant(kFALSE)
    
    #Starting values of the parameters
    mean1.setVal(1.0)
    sigma1.setVal(1.0)
    m.setVal(90.0)
    s.setVal(1.0)
    a.setVal(10.0)
    n.setVal(2.0)

    # Construct CB (x) gauss
    x1.setBins(10000, "cache")
    CBxG = RooFFTConvPdf("CBxG", "CB (X) gauss", x1, CB, gauss1)
    
    can.cd()
    d = RooDataHist("d","d",RooArgList(x1),RooFit.Import(histo))
    CBxG.fitTo(d, RooLinkedList())
   
    # Plot PDF and toy data overlaid
    xframe2 = x1.frame(RooFit.Name("xframe"),RooFit.Title("")) # RooPlot
    d.plotOn(xframe2, RooLinkedList() )
    CBxG.paramOn(xframe2, RooFit.Layout(0.65,0.99,0.9))
    xframe2.getAttText().SetTextSize(0.03)
    CBxG.plotOn(xframe2)
    xframe2.Draw()
    can.SaveAs("DataVsMC/FitResults/"+title+"_Roofit.pdf")
    can.SaveAs("DataVsMC/FitResults/"+title+"_Roofit.png")
    
    return;
