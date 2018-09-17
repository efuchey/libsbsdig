#ifndef _TSBSSIMHCAL_H
#define _TSBSSIMHCAL_H

#include "TSBSSimDetector.h"
#include "TSBSSimAuxi.h"

// class TF1;
// class TF1Convolution;
// class TTree;
// class TFile;

class TSBSSimHCal : public TSBSSimDetector {
public:
  TSBSSimHCal(const char* name, short id);
  virtual ~TSBSSimHCal();
  // This loads the simulation event data
  virtual void LoadEventData(const std::vector<g4sbshitdata*> &evbuffer);
  //The following function is to accumulate data without reinitializing the event.
  virtual void LoadAccumulateData(const std::vector<g4sbshitdata*> &evbuffer);
  virtual void Digitize(TSBSSimEvent &event);

  virtual void Clear(Option_t *op = "");

  // Silence compiler warnings about Init from parent class
  using THaAnalysisObject::Init;
  // Initialize
  void Init();

  struct Signal {
    SimEncoder::fadc_data fadc;
    SimEncoder::tdc_data tdc;
    //std::vector<double> samples;
    std::vector<double> samples_raw;
    std::vector<double> times_histo;
    int nbins_times;
    double sumedep;
    double mint;
    double maxt;
    int nbins;
    int nbins_raw;
    int npe;
    double tdc_time;
    bool met_tdc_thresh;
    int sum;
    int dnraw;
    double dx_samples;
    double dx_raw;
    double dx_raw_time;
    Signal();
    void FillNPE(TSPEModel *model, double pulsenorm, double t, double toffset = 0.0);
    void Fill(double t);
    void Digitize(TSPEModel *model, double pulsenorm, double toffset,
        double max_val);
    void DigitizeOld();
    void Clear();
  };
private:
  //SPEModel *fSPE;
  TSPEModel *fSPE;
  std::vector<Signal> fSignals;
  bool fHasFADC;
  // TODO: Try to use the standard TPMTSignal class (but must have
  // the ability to provide samples)
  //std::vector<TPMTSignal> fSignals;

  ClassDef(TSBSSimHCal,1)
};

#endif //_TSBSSIMHCAL_H
