# **********************
#
# This macro:
#  - reads data from the input file DataSet.root (jPsi and Psi(2s) mass 
#  - fits the JPsi and Psi(2s) mass distrbution and plot the result
#  - 

import ROOT

fInput = ROOT.TFile("DataSet.root")
fInput.cd()
data = fInput.Get("data")

# variable
mass = ROOT.RooRealVar("mass","Invariant mass",2.,6.,"GeV/c^2")

# JPsi: CB
mean_JPsi  = ROOT.RooRealVar("mean_JPsi","mean_JPsi",3.1,2.8,3.5)
sigma_JPsi = ROOT.RooRealVar("sigma_JPsi","sigma_JPsi",0.2,0.0001,1.)
alpha_JPsi = ROOT.RooRealVar("alpha_JPsi","alpha_JPsi",1.,-5.,5.)
n_JPsi     = ROOT.RooRealVar("n_JPsi","n_JPsi",1.,-5.,5.)
JPsiPDF    = ROOT.RooCBShape("JPsiPDF","JPsiPDF",mass,mean_JPsi,sigma_JPsi,alpha_JPsi,n_JPsi)

# Psi(2s): CB
mean_Psi = ROOT.RooRealVar("mean_Psi","mean_Psi",3.7, 3.2,4.1)
PsiPDF   = ROOT.RooCBShape("PsiPDF","PsiPDF",mass,mean_Psi,sigma_JPsi,alpha_JPsi,n_JPsi)

# backgroud: Chebychev
a0 = ROOT.RooRealVar("a0","a0",-0.1,-2.,2.)
a1 = ROOT.RooRealVar("a1","a1",0.5,-2.,2.)
a2 = ROOT.RooRealVar("a2","a2",-0.1,-2.,2.)
bkgPDF = ROOT.RooChebychev("bkgPDF","bkgPDF",mass,ROOT.RooArgList(a0,a1,a2))

# pdf tot
NJPsi  = ROOT.RooRealVar("NJPsi","NJPsi",2000.,0.1,20000.)
NPsi   = ROOT.RooRealVar("NPsi","NPsi",100.,0.001,2000.)
Nbkg   = ROOT.RooRealVar("Nbkg","Nbkg",20000.,0.1,100000.)
totPDF = ROOT.RooAddPdf("totPDF","totPDF",ROOT.RooArgList(JPsiPDF,PsiPDF,bkgPDF), ROOT.RooArgList(NJPsi,NPsi,Nbkg) )

# fit
totPDF.fitTo(data, ROOT.RooFit.Extended(1))

# plot
massplot = mass.frame()
data.plotOn(massplot)
totPDF.plotOn(massplot)

canvas = ROOT.TCanvas()
canvas.cd()
massplot.Draw()
canvas.SaveAs("fit.png")

# save in workspace
workspace = ROOT.RooWorkspace("workspace")
getattr(workspace,'import')(totPDF)
getattr(workspace,'import')(data)

workspace.Print()

fOutput = ROOT.TFile("final_fit.root","RECREATE")
fOutput.cd()
workspace.Write()
fOutput.Close()






