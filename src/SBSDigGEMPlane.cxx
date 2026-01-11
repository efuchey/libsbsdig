#include "SBSDigGEMPlane.h"
#include "TMath.h"

using namespace std;

//
// Class SBSDigGEMPlane
//
SBSDigGEMPlane::SBSDigGEMPlane() :
  fLayer(0), fNStrips(3840), fNSamples(6), fStripThr(100), fXoffset(0), fROangle(0)
{  
  // SetPlaneStripPed();
  // SetPlaneAPVCM();
}

SBSDigGEMPlane::SBSDigGEMPlane(short layer, short mod, int nstrips, int nsamples, double thr, double offset, double roangle) :
  fLayer(layer), fModule(mod), fNStrips(nstrips), fNSamples(nsamples), fStripThr(thr), fXoffset(offset), fROangle(roangle)
{
  //fModule = mod;
  fdX = fNStrips*4.e-4;
  
  fStripADCsum = new Int_t[fNStrips];
  fStripADC = new Short_t[fNStrips*fNSamples];
  fStripMipADC = new Short_t[fNStrips*fNSamples];
  fStripPedSubADC = new Short_t[fNStrips*fNSamples];
  fStripCMCorrADC = new Short_t[fNStrips*fNSamples];

  fNAPVs = int(fNStrips/fNChanAPV);


  Clear();
}

SBSDigGEMPlane::~SBSDigGEMPlane()
{
  Clear();
}


void SBSDigGEMPlane::Clear()
{
  memset(fStripADCsum, 0, fNStrips*sizeof(Int_t));
  memset(fStripADC, 0, fNStrips*fNSamples*sizeof(Short_t));
  memset(fStripPedSubADC, 0, fNStrips*fNSamples*sizeof(Short_t));
  memset(fStripCMCorrADC, 0, fNStrips*fNSamples*sizeof(Short_t));
  memset(fStripMipADC, 0, fNStrips*fNSamples*sizeof(Short_t));
}

//ClassImp(SBSDigGEMPlane);

// 'Online' CM and ZS.

void SBSDigGEMPlane::DoPedSub(){ // First subtract pedestals from each channel before calculating online CM.

  for ( int istrip=0; istrip < fNStrips; istrip++ ){
    for ( int isamp=0; isamp < fNSamples; isamp++ ){
      fStripPedSubADC[istrip*fNSamples+isamp] = fStripADC[istrip*fNSamples+isamp] - fPedestalMap[istrip].mean; 
    }    
  }
}

double SBSDigGEMPlane::GetOnlineCommonMode(const std::array <int,fNChanAPV>& adcarray, int apvnum){ // We are replicating 'falg 4' in SBSGEMModule::GetCommonMode() method here.
  
  double cm_mean = GetAPVCMMean(apvnum);
  double cm_rms = GetAPVCMSigma(apvnum); // Here I'm changing CM sigma to RMS to go with the SBS_OFFLINE convension. But shouldn't it be sigma as true CM values is not 0 ?

  double cm_temp = 0.;

  int n_keep = 0;

  for (int iter=0; iter<3; iter++){

    double cm_min = cm_mean - fCommonModeRange_nsigma*cm_rms;

    if (iter==0) cm_min = 0.;

    double cm_max = cm_mean + fCommonModeRange_nsigma*cm_rms;
    double sumADCinrange = 0.;
    n_keep = 0;

    for (int ichan=0; ichan<fNChanAPV; ichan++){//Size of the array will be 128.

      double ADCtemp = (double)adcarray[ichan];
      double pedrmstemp = fPedestalMap[apvnum*fNChanAPV+ichan].rms;

      if (iter != 0){
        cm_min = cm_temp - fCommonModeDanningMethod_NsigmaCut*2.5*pedrmstemp;
        cm_max = cm_temp + fCommonModeDanningMethod_NsigmaCut*2.5*pedrmstemp;
      }

      if ( ADCtemp >= cm_min && ADCtemp <= cm_max ){
        n_keep++;
        sumADCinrange += ADCtemp;
      }
    }

    if (n_keep == 0) return cm_mean;

    cm_temp = sumADCinrange / n_keep;
  }

  if( n_keep < fCommonModeMinStripsInRange ) return cm_mean;

  return cm_temp;

}

void SBSDigGEMPlane::ApplyOnlineCMCorr(){ // Calculate CM for per 128 strip-set (per APV) for each TS.
// Then apply the CM correction to the fStripADC array. 

  DoPedSub();// First have to subtract pedestal offsets.

  for ( int istrip=0; istrip < fNStrips; istrip += fNChanAPV ){

    int iAPV = int(istrip/fNChanAPV);

    for ( int isamp=0; isamp < fNSamples; isamp++ ){

      std::array <int, fNChanAPV> thisAPVthisSampPedSubADC{-1500};

      for ( int ichan=0; ichan < fNChanAPV; ichan++ ){
        thisAPVthisSampPedSubADC[ichan] = fStripPedSubADC[(istrip+ichan)*fNSamples+isamp];
      }

      double thisAPVthisSampOnineCM = GetOnlineCommonMode(thisAPVthisSampPedSubADC, iAPV);

      for ( int ichan=0; ichan < fNChanAPV; ichan++ ){
        fStripCMCorrADC[(istrip+ichan)*fNSamples+isamp] = TMath::Nint(fStripPedSubADC[(istrip+ichan)*fNSamples+isamp] - thisAPVthisSampOnineCM);
        //fStripADC[(istrip+ichan)*fNSamples+isamp] = fStripCMCorrADC[(istrip+ichan)*fNSamples+isamp];
        SetADC((istrip+ichan), isamp, fStripCMCorrADC[(istrip+ichan)*fNSamples+isamp]);
      }
    }
  }  
}

void SBSDigGEMPlane::ApplyOnlineZS(const double zs_thr_nsigma){

  for ( int istrip=0; istrip < fNStrips; istrip++ ){

    bool stripPassZS = GetADCSum(istrip)/double(fNSamples) > zs_thr_nsigma*fPedestalMap[istrip].rms;

    // zero-suppress the strips that did not pass the threshold.
    if (!stripPassZS){
      for ( int isamp=0; isamp < fNSamples; isamp++ ){
        SetADC(istrip, isamp, 0);
      }
    }
  }
}