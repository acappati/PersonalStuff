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
#include <TAttFill.h>


using namespace std;


void Decay(Int_t N0 = 100, Float_t alpha = 0.01, Float_t dt = 1, Int_t ndt = 300, UInt_t seed = 2126)
{

  // prob of a nucleus to decay 
  Float_t decay_prob = alpha * dt; 


  // random number for the loop 
  TRandom3 *rand = new TRandom3();
  gRandom->SetSeed(seed);
 

  // histo for the number of remaining nuclei after a certain time 
  TH1I *Nremaining_hist = new TH1I("Nremaining_hist","Remaining Nuclei",ndt,0,ndt);
  Nremaining_hist->GetXaxis()->SetTitle(Form("Number of time intervals (dt = %.2f s)",dt));
  Nremaining_hist->GetYaxis()->SetTitle("N_{remaining nuclei}");
  Nremaining_hist->SetLineColor(kBlack);
  Nremaining_hist->SetFillColor(kGreen-9);


  // start from number of remaining nuclei: initially is N0 
  Int_t N = N0;
 
  for(Int_t i=0; i<ndt; i++){      // loop from t=0 to t in step of dt (ndt=number of time delay)
    for(Int_t j=0; j<N; j++){      // loop over each remaining parent nucleus 

      if(rand->Rndm() < decay_prob) N--; // gen random number, if it is lower than decay prob, nucleus decays

    }
    Nremaining_hist->Fill(i+0.5,N);
  }


  // plot the histo
  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
  canvas->Draw();
  canvas->cd();

  Nremaining_hist->Draw("histo");

  canvas->Update();
  
  canvas->SaveAs("Remaining_nuclei.pdf");



  

}
