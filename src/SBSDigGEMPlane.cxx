#include "SBSDigGEMPlane.h"
#include "TMath.h"

using namespace std;

//
// Class SBSDigGEMPlane
//
SBSDigGEMPlane::SBSDigGEMPlane() :
  fLayer(0), fNStrips(3840), fNSamples(6), fStripThr(100), fXoffset(0), fROangle(0)
{  
  SetPlaneStripPed();
  SetPlaneAPVCM();
}

SBSDigGEMPlane::SBSDigGEMPlane(short layer, short mod, int nstrips, int nsamples, double thr, double offset, double roangle) :
  fLayer(layer), fModule(mod), fNStrips(nstrips), fNSamples(nsamples), fStripThr(thr), fXoffset(offset), fROangle(roangle)
{
  //fModule = mod;
  fdX = fNStrips*4.e-4;
  
  fStripADCsum = new Int_t[fNStrips];
  fStripADC = new Short_t[fNStrips*fNSamples];
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
}

//ClassImp(SBSDigGEMPlane);

