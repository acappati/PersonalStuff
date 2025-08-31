// **************************
//
// purpose: compute the number of decays in a fixed time
//
// usage:
//        root
//        .L DecaysInFixedTime.C++
//        DecaysInFixedTime()
//
// **************************


#include <TRandom3.h>
#include <TH1I.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TAttFill.h>
#include <TTimeStamp.h>


using namespace std;


Float_t kNorm = 1000.; // number of experiments (= normalization for theoretical functions)


// *****
// function which return the number of nuclei which decay in time T
Float_t Decay(Int_t N0, Float_t alpha, Float_t dt, Float_t timeTot, bool changeSeed = kTRUE)
{

  // prob of a nucleus to decay 
  Float_t decay_prob = alpha * dt; 
  
  // random number for the loop 
  TRandom3 *rand = new TRandom3();           // use default TRandom3 seed = 4357
  if(changeSeed){
    TTimeStamp *tim = new TTimeStamp();
    rand->SetSeed(tim->GetTime());	     // can change the seed according to machine time
  }

  // initially the number of remaining nuclei is N0 
  Int_t N = N0;
 
  for(Float_t time=0; time<timeTot; time += dt){      // loop from t=0 to t in step of dt 
    for(Int_t nuclei=0; nuclei<N; nuclei++){          // loop over each remaining parent nucleus 

      if(rand->Rndm() < decay_prob) N--;   // gen random number, if it is lower than decay prob, nucleus decays 
    }
  }

  return N0-N; // return the number of decaded nuclei

}






// *****
// main function
void DecaysInFixedTime(Int_t N0 = 1000, Float_t alpha = 2.e-5, Float_t dt = 1., Float_t timeTot = 100.)
{

  //... slides 5 tans

}
