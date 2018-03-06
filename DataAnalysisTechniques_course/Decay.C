// **************************
//
// purpose: compute the number of remaining nuclei after a radioactive decay
//
// usage:
//        root
//        .L Decay.C++
//        Decay(number of initial nuclei, decay constant, time interval in s, number of time intervals after which counting remaining nuclei, seed)
//
// **************************


#include <TRandom3.h>
#include <TH1I.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TAttFill.h>


using namespace std;



// function for the fit 
Double_t Exp_function(Double_t *variable,Double_t *parameter){

  Double_t x = variable[0];

  return parameter[0] * TMath::Exp(- parameter[1] * x);

}


// main function
void Decay(Int_t N0 = 100, Float_t alpha = 0.01, Float_t dt = 1, UInt_t seed = 2126)
{

  // prob of a nucleus to decay 
  Float_t decay_prob = alpha * dt; 

  // number of time intervals 
  Int_t nBins = 300;

  // total time 
  Float_t timeTot = dt * nBins;


  // random number for the loop 
  TRandom3 *rand = new TRandom3();
  //gRandom->SetSeed(seed);
  rand->SetSeed(seed);
 

  // histo for the number of remaining nuclei after a certain time 
  TH1I *Nremaining_hist = new TH1I("Nremaining_hist","Remaining Nuclei",nBins,0.,timeTot);
  Nremaining_hist->GetXaxis()->SetTitle(Form("Number of time intervals (dt = %.2f s)",dt));
  Nremaining_hist->GetYaxis()->SetTitle("N_{remaining nuclei}");
  Nremaining_hist->SetLineColor(kBlack);
  Nremaining_hist->SetFillColor(kGreen-9);


  // start from number of remaining nuclei: initially is N0 
  Int_t N = N0;
 
  for(Float_t time=0; time<timeTot; time += dt){      // loop from t=0 to t in step of dt (ndt=number of time delay)
    for(Int_t nuclei=0; nuclei<N; nuclei++){      // loop over each remaining parent nucleus 

      if(rand->Rndm() < decay_prob) N--; // gen random number, if it is lower than decay prob, nucleus decays

    }
    Nremaining_hist->Fill(time,N);
  }


  // plot the histo
  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
  canvas->Draw();
  canvas->cd();

  Nremaining_hist->Draw("histo");

  canvas->Update();


  // plot theoretical decay function  N = N0 * e^(-alpha*t)
  TF1 *functionTh = new TF1("functionTh",Exp_function,0,timeTot,2);
  functionTh->SetParameter(0,N0);
  functionTh->SetParameter(1,alpha);
  functionTh->SetLineColor(kRed);
  functionTh->Draw("same");

  canvas->Update();

  
  canvas->SaveAs("Remaining_nuclei.png");



  

}
