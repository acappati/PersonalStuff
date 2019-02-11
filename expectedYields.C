// ********************************
// 
// this script computes the expected yields for HH->4l+X
//
// run with:
//        root -l -b -q expectedYields.C++
//
// *******************************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

using namespace std;




void expectedYields(){


  double higgsMass = 125.09; //GeV

  double lumiRun2  = 140.;   //fb-1
  double xSecHH    = 33.49;  //fb

 
  // compute BR
  double brHtobb         = 5.809 * 0.1;
  double brHToWW         = 2.152 * 0.1;
  double brHtoTauTau     = 6.256 * 0.01; 
  double brHToGammaGamma = 2.270 * 0.001;
  double brHToZZ         = 2.641 * 0.01;

  double brZtoee   = 3.3632 * 0.01;
  double brZtomumu = 3.3662 * 0.01;

  double brWtoeNu  = 10.71 * 0.01;
  double brWtomuNu = 10.63 * 0.01;


  
  double brHto4l = brHToZZ *(brZtoee + brZtomumu)*(brZtoee + brZtomumu); 
  
  double brHtoWW2l2Nu = brHToWW *(brWtoeNu + brWtomuNu)*(brWtoeNu + brWtomuNu);

  cout<<"ZZ->4l  : "<<(brZtoee + brZtomumu)*(brZtoee + brZtomumu)<<endl;
  cout<<"brHto4l : "<<brHto4l<<endl;
  cout<<"brHtoWW2l2Nu : "<<brHtoWW2l2Nu<<endl;


  double brHHto4lbb         = brHto4l * brHtobb;
  double brHHto4ltautau     = brHto4l * brHtoTauTau;
  double brHHto4lGammaGamma = brHto4l * brHToGammaGamma;
  double brHHto4lWW2l2Nu    = brHto4l * brHtoWW2l2Nu;
  double brHHto8l           = brHto4l * brHto4l;

  cout<<"BR "<<endl;
  cout<<"brHHto4lbb         : "<<brHHto4lbb         <<endl;                
  cout<<"brHHto4ltautau     : "<<brHHto4ltautau     <<endl;
  cout<<"brHHto4lGammaGamma : "<<brHHto4lGammaGamma <<endl;
  cout<<"brHHto4lWW2l2Nu    : "<<brHHto4lWW2l2Nu    <<endl;   
  cout<<"brHHto8l           : "<<brHHto8l           <<endl;  
  cout<<"-----------------"<<endl;


  // compute expected yields 
  double yieldsHHto4lbb          = lumiRun2 * xSecHH * brHHto4lbb;
  double yieldsHHto4ltautau      = lumiRun2 * xSecHH * brHHto4ltautau;
  double yieldsHHto4lGammaGamma  = lumiRun2 * xSecHH * brHHto4lGammaGamma;
  double yieldsHHto4lWW2l2Nu     = lumiRun2 * xSecHH * brHHto4lWW2l2Nu;
  double yieldsHHto8l            = lumiRun2 * xSecHH * brHHto8l;

  cout<<"Expected Yields "<<endl;
  cout<<"yieldsHHto4lbb         : "<<yieldsHHto4lbb         <<endl;                
  cout<<"yieldsHHto4ltautau     : "<<yieldsHHto4ltautau     <<endl;
  cout<<"yieldsHHto4lGammaGamma : "<<yieldsHHto4lGammaGamma <<endl;
  cout<<"yieldsHHto4lWW2l2Nu    : "<<yieldsHHto4lWW2l2Nu    <<endl;   
  cout<<"yieldsHHto8l           : "<<yieldsHHto8l           <<endl;  
  cout<<"-----------------"<<endl;

}
