#ifndef SBSDIGBKGDGEN_H
#define SBSDIGBKGDGEN_H

#include <iostream>
#include <vector>
#include <map>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "SBSDigPMTDet.h"
#include "SBSDigGEMDet.h"


//________________________________
class SBSDigBkgdGen {

 public:
  SBSDigBkgdGen();
  SBSDigBkgdGen(TFile* f_bkgd, double timewindow, bool pmtbkgddig);
  ~SBSDigBkgdGen();
  void Initialize(TFile* f_bkgd);
  
  void GenerateBkgd(TRandom3* R, 
		    std::vector<SBSDigPMTDet*> pmtdets,
		    std::vector<int> detmap, 
		    std::vector<SBSDigGEMDet*> gemdets, 
		    std::vector<int> gemmap, 
		    double lumifrac);
  
  void WriteXCHistos();
  
 private:
  double fTimeWindow;
  
  Double_t* NhitsBBGEMs;
  TH1D* h_EdephitBBGEMs;
  TH1D** h_xhitBBGEMs;
  TH1D** h_yhitBBGEMs;
  TH1D** h_dxhitBBGEMs;
  TH1D** h_dyhitBBGEMs;
  
  Double_t* NhitsHCal;
  TH1D* h_EdephitHCal;
  TH1D* h_zhitHCal;
  
  Double_t* NhitsBBPS;
  TH1D* h_EdephitBBPS;
  
  Double_t* NhitsBBSH;
  TH1D* h_EdephitBBSH;
  
  Double_t* NhitsBBHodo;
  TH1D* h_EdephitBBHodo;
  TH1D* h_xhitBBHodo;
  
  Double_t* P1hitGRINCH;
  Double_t* P2hitsGRINCH;
  TH1D* h_NpeGRINCH;
  
  Double_t* NhitsECal;
  TH1D* h_EdephitECal;
  
  Double_t* NhitsCDet;
  TH1D* h_EdephitCDet;
  TH1D* h_xhitCDet;
  
  Double_t* NhitsFT;
  TH1D* h_EdephitFT;
  TH1D** h_xhitFT;
  TH1D** h_yhitFT;
  TH1D** h_dxhitFT;
  TH1D** h_dyhitFT;
  
  Double_t* NhitsFPP1;
  TH1D* h_EdephitFPP1;
  TH1D** h_xhitFPP1;
  TH1D** h_yhitFPP1;
  TH1D** h_dxhitFPP1;
  TH1D** h_dyhitFPP1;
  
  Double_t* NhitsFPP2;
  TH1D* h_EdephitFPP2;
  TH1D** h_xhitFPP2;
  TH1D** h_yhitFPP2;
  TH1D** h_dxhitFPP2;
  TH1D** h_dyhitFPP2;
  
  /*  
  //cross-check histos
  TH1D** h_NhitsBBGEMs_XC;
  TH1D* h_EdephitBBGEMs_XC;
  TH1D** h_xhitBBGEMs_XC;
  TH1D** h_yhitBBGEMs_XC;
  TH1D* h_modBBGEMs_XC;
  
  TH2D* h_NhitsHCal_XC;
  TH1D* h_EdephitHCal_XC;
  TH1D* h_zhitHCal_XC;
  
  TH2D* h_NhitsBBPS_XC;
  TH1D* h_EdephitBBPS_XC;
  
  TH2D* h_NhitsBBSH_XC;
  TH1D* h_EdephitBBSH_XC;
  
  TH2D* h_NhitsBBHodo_XC;
  TH1D* h_EdephitBBHodo_XC;
  TH1D* h_xhitBBHodo_XC;
  
  TH2D* h_NhitsGRINCH_XC;
  TH1D* h_NpeGRINCH_XC;
  */
  
  bool fPMTBkgdDig;
  
};

#endif // SBSDIGAUXI_H

