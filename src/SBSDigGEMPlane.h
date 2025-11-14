#ifndef SBSDIGGEMPLANE_H
#define SBSDIGGEMPLANE_H

#include <iostream>
#include <vector>
#include <map>
#include <TROOT.h>
//#include "TH1D.h"
#include "TRandom3.h"

class SBSDigGEMPlane {
 public:
  SBSDigGEMPlane();
  SBSDigGEMPlane(short layer, short mod, int nstrips, int nsamples = 6, double thr = 100, double offset = 0, double roangle = 0);
  virtual ~SBSDigGEMPlane();
  void Clear();
  
  //void SetStripThreshold(double thr){Striphr = ;};

  Short_t Layer(){return fLayer;};
  Short_t Module(){return fModule;};
  Double_t dX(){return fdX;};
  Double_t Xoffset(){return fXoffset;};
  Double_t ROangle(){return fROangle;};
  Int_t GetNStrips(){return fNStrips;};
  Short_t GetADC(int strip, int samp){return fStripADC[strip*fNSamples+samp];};
  Int_t GetADCSum(int strip){return fStripADCsum[strip];};
  void SetADC(int strip, int samp, int adc){
    if(strip<fNStrips){
      fStripADCsum[strip]+= adc-fStripADC[strip*fNSamples+samp];
      fStripADC[strip*fNSamples+samp] = adc;
    }
  };
  void AddADC(int strip, int samp, int adc){
    if(strip<fNStrips){
      fStripADC[strip*fNSamples+samp]+=adc;
      fStripADCsum[strip]+=adc;
    }//else{
    //printf("strip = %d / %d, sample %d /%d \n", strip, fNStrips, samp, fNSamples);
    //}
  };

  void SetPlaneStripPed(){ // This function will set the pedestal mean and RMS for each strip in the GEM plane.

    for (int istrip=0; istrip<fNStrips; istrip++){

      fStripPedOffset[istrip] = R->Gaus(0,50.); //using some hard-coded values here just for testing.
      fStripPedRMS[istrip] = R->Gaus(0, 9.);
    }
  };

  void SetPlaneAPVCM(){  // This function will set CM mean and std.dev values for each APV in the GEM plane.

    for (int istrip=0; istrip<fNStrips; istrip++){

      if ( istrip%128==0 ){
        int iAPV = istrip/128;
        fAPVCMMean[iAPV] = R->Gaus(1500, 200);
        fAPVCMSigma[iAPV] = 20;
      }
    }
  };

  Double_t GetStripPedMean(int strip) {return fStripPedOffset[strip];};
  Double_t GetStripPedRMS(int strip) {return fStripPedRMS[strip];};
  Double_t GetAPVCMMean(int apv) {return fAPVCMMean[apv];};
  Double_t GetAPVCMSigma(int apv) {return fAPVCMSigma[apv];};
  
 private:
  // ADC sampled value of strip array of each axis
  Short_t fLayer;
  Short_t fModule;
  Int_t fNStrips;
  Int_t fNSamples;
  Double_t fStripThr;//threshold for ADC sum
  Int_t* fStripADCsum;
  Short_t* fStripADC;
  double fdX;
  double fXoffset;
  double fROangle;
  //ClassDef(SBSDigGEMPlane, 1)

  // Ped and CM values.
  TRandom3* R = new TRandom3(0); //Random number generation for ped and CM testing purposes.
  Double_t* fStripPedOffset;
  Double_t* fStripPedRMS;
  Double_t* fAPVCMMean;
  Double_t* fAPVCMSigma;
  
};
#endif
