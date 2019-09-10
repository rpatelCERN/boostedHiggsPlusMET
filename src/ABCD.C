#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>


void makeABCDPlot(vector<TH1F*> dem_histos, TString type = "", TString tagType = "");
void makeRpfPlot(vector<TH1F*> dem_histos, TString bkgType = "", TString jetType = "", TString tagType = "");

void makeDoubleBPlots(vector<TH1F*> dem_histos, TString bkgType = "", TString which = "");
void makeStackPlot(vector<TH1F*> h_QCD,vector<TH1F*> h_TT,vector<TH1F*> h_WJets,vector<TH1F*> h_ZJets, vector<TH1F*> h_T5HH1300,vector<TH1F*> h_T5HH1700,vector<TH1F*> h_T5HH2100,TString which = "");
void makeSingleLeptStackPlot(vector<TH1F*> h_TT,vector<TH1F*> h_WJets,TString which = "");

void QuickROC(TH1F* signalHist, TH1F* bkgHist, TString which = "");

void makeMETStack(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets, TH1F* h_T5HH1300,TH1F* h_T5HH1700,TH1F* h_T5HH2100,TString which = "");
void makeDoubleBStack(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets, TH1F* h_T5HH1300,TH1F* h_T5HH1700,TH1F* h_T5HH2100, TString which = "");

void tableMCyields();

ofstream myfile;
bool getDeviation = false;
bool runABCDPlots = false;
bool runRPFPlots = false;
bool runDoubleBStack = false;
bool runCompareDoubleB = false;
bool runStacks = false;
bool runMETStacks = true;
bool runROC = false;

bool useDeepDoubleB = false;

TString directory = "lowMETbin/";
// TString directory = "medB/";


string whichRegion = "signal";
// string whichRegion = "singleLept";
// string whichRegion = "photon";


void ABCD() {
  TFile * f = TFile::Open("ALPHABEThistos_V12.root");
  // TFile * f = TFile::Open("ALPHABEThistos_V12_medB.root");


  // TFile * fSignal = TFile::Open("ALPHABEThistosMC2016_signalOnly.root");
  // TFile * fPhoton = TFile::Open("ALPHABEThistos_photon.root");
  // TFile * fSingleLept = TFile::Open("ALPHABEThistos_singleLept.root");

  TH1F * h_A_sum; TH1F * h_B_sum; TH1F * h_A1_sum; TH1F * h_B1_sum; TH1F *h_C_sum; TH1F *h_D_sum;
  TH1F * h_A_QCD; TH1F * h_B_QCD; TH1F * h_A1_QCD; TH1F * h_B1_QCD; TH1F * h_C_QCD; TH1F * h_D_QCD;
  TH1F * h_A_GJets; TH1F * h_B_GJets; TH1F * h_A1_GJets; TH1F * h_B1_GJets; TH1F * h_C_GJets; TH1F * h_D_GJets;
  TH1F * h_A_TT; TH1F * h_B_TT; TH1F * h_A1_TT; TH1F * h_B1_TT; TH1F * h_C_TT; TH1F * h_D_TT;
  TH1F * h_A_WJets; TH1F * h_B_WJets; TH1F * h_A1_WJets; TH1F * h_B1_WJets; TH1F * h_C_WJets; TH1F * h_D_WJets;
  TH1F * h_A_ZJets; TH1F * h_B_ZJets; TH1F * h_A1_ZJets; TH1F * h_B1_ZJets; TH1F * h_C_ZJets; TH1F * h_D_ZJets;

  TH1F * h_A_T5HH1300; TH1F * h_B_T5HH1300; TH1F * h_A1_T5HH1300; TH1F * h_B1_T5HH1300; TH1F * h_C_T5HH1300; TH1F * h_D_T5HH1300;
  TH1F * h_A_T5HH1700; TH1F * h_B_T5HH1700; TH1F * h_A1_T5HH1700; TH1F * h_B1_T5HH1700; TH1F * h_C_T5HH1700; TH1F * h_D_T5HH1700;
  TH1F * h_A_T5HH2100; TH1F * h_B_T5HH2100; TH1F * h_A1_T5HH2100; TH1F * h_B1_T5HH2100; TH1F * h_C_T5HH2100; TH1F * h_D_T5HH2100;


  TH1F * h_J1M_doubletagSR_sum; TH1F * h_J2M_doubletagSR_sum; TH1F * h_J1M_doubletagSB_sum; TH1F * h_J2M_doubletagSB_sum;
  TH1F * h_J1M_antitagSR_sum; TH1F * h_J2M_antitagSR_sum; TH1F * h_J1M_antitagSB_sum; TH1F * h_J2M_antitagSB_sum;
  TH1F * h_J1M_tagSR_sum; TH1F * h_J2M_tagSR_sum; TH1F * h_J1M_tagSB_sum; TH1F * h_J2M_tagSB_sum;

  TH1F * h_J1M_doubletagSR_QCD; TH1F * h_J2M_doubletagSR_QCD; TH1F * h_J1M_doubletagSB_QCD; TH1F * h_J2M_doubletagSB_QCD;
  TH1F * h_J1M_antitagSR_QCD; TH1F * h_J2M_antitagSR_QCD; TH1F * h_J1M_antitagSB_QCD; TH1F * h_J2M_antitagSB_QCD;
  TH1F * h_J1M_tagSR_QCD; TH1F * h_J2M_tagSR_QCD; TH1F * h_J1M_tagSB_QCD; TH1F * h_J2M_tagSB_QCD;

  TH1F * h_J1M_doubletagSR_GJets; TH1F * h_J2M_doubletagSR_GJets; TH1F * h_J1M_doubletagSB_GJets; TH1F * h_J2M_doubletagSB_GJets;
  TH1F * h_J1M_antitagSR_GJets; TH1F * h_J2M_antitagSR_GJets; TH1F * h_J1M_antitagSB_GJets; TH1F * h_J2M_antitagSB_GJets;
  TH1F * h_J1M_tagSR_GJets; TH1F * h_J2M_tagSR_GJets; TH1F * h_J1M_tagSB_GJets; TH1F * h_J2M_tagSB_GJets;

  TH1F * h_J1M_doubletagSR_TT; TH1F * h_J2M_doubletagSR_TT; TH1F * h_J1M_doubletagSB_TT; TH1F * h_J2M_doubletagSB_TT;
  TH1F * h_J1M_antitagSR_TT; TH1F * h_J2M_antitagSR_TT; TH1F * h_J1M_antitagSB_TT; TH1F * h_J2M_antitagSB_TT;
  TH1F * h_J1M_tagSR_TT; TH1F * h_J2M_tagSR_TT; TH1F * h_J1M_tagSB_TT; TH1F * h_J2M_tagSB_TT;

  TH1F * h_J1M_doubletagSR_WJets; TH1F * h_J2M_doubletagSR_WJets; TH1F * h_J1M_doubletagSB_WJets; TH1F * h_J2M_doubletagSB_WJets;
  TH1F * h_J1M_antitagSR_WJets; TH1F * h_J2M_antitagSR_WJets; TH1F * h_J1M_antitagSB_WJets; TH1F * h_J2M_antitagSB_WJets;
  TH1F * h_J1M_tagSR_WJets; TH1F * h_J2M_tagSR_WJets; TH1F * h_J1M_tagSB_WJets; TH1F * h_J2M_tagSB_WJets;

  TH1F * h_J1M_doubletagSR_ZJets; TH1F * h_J2M_doubletagSR_ZJets; TH1F * h_J1M_doubletagSB_ZJets; TH1F * h_J2M_doubletagSB_ZJets;
  TH1F * h_J1M_antitagSR_ZJets; TH1F * h_J2M_antitagSR_ZJets; TH1F * h_J1M_antitagSB_ZJets; TH1F * h_J2M_antitagSB_ZJets;
  TH1F * h_J1M_tagSR_ZJets; TH1F * h_J2M_tagSR_ZJets; TH1F * h_J1M_tagSB_ZJets; TH1F * h_J2M_tagSB_ZJets;

  TH1F * h_J1M_doubletagSR_T5HH1300; TH1F * h_J2M_doubletagSR_T5HH1300; TH1F * h_J1M_doubletagSB_T5HH1300; TH1F * h_J2M_doubletagSB_T5HH1300;
  TH1F * h_J1M_antitagSR_T5HH1300; TH1F * h_J2M_antitagSR_T5HH1300; TH1F * h_J1M_antitagSB_T5HH1300; TH1F * h_J2M_antitagSB_T5HH1300;
  TH1F * h_J1M_tagSR_T5HH1300; TH1F * h_J2M_tagSR_T5HH1300; TH1F * h_J1M_tagSB_T5HH1300; TH1F * h_J2M_tagSB_T5HH1300;

  TH1F * h_J1M_doubletagSR_T5HH1700; TH1F * h_J2M_doubletagSR_T5HH1700; TH1F * h_J1M_doubletagSB_T5HH1700; TH1F * h_J2M_doubletagSB_T5HH1700;
  TH1F * h_J1M_antitagSR_T5HH1700; TH1F * h_J2M_antitagSR_T5HH1700; TH1F * h_J1M_antitagSB_T5HH1700; TH1F * h_J2M_antitagSB_T5HH1700;
  TH1F * h_J1M_tagSR_T5HH1700; TH1F * h_J2M_tagSR_T5HH1700; TH1F * h_J1M_tagSB_T5HH1700; TH1F * h_J2M_tagSB_T5HH1700;

  TH1F * h_J1M_doubletagSR_T5HH2100; TH1F * h_J2M_doubletagSR_T5HH2100; TH1F * h_J1M_doubletagSB_T5HH2100; TH1F * h_J2M_doubletagSB_T5HH2100;
  TH1F * h_J1M_antitagSR_T5HH2100; TH1F * h_J2M_antitagSR_T5HH2100; TH1F * h_J1M_antitagSB_T5HH2100; TH1F * h_J2M_antitagSB_T5HH2100;
  TH1F * h_J1M_tagSR_T5HH2100; TH1F * h_J2M_tagSR_T5HH2100; TH1F * h_J1M_tagSB_T5HH2100; TH1F * h_J2M_tagSB_T5HH2100;

  //For pt
  TH1F * h_J1Pt_doubletagSR_sum; TH1F * h_J2Pt_doubletagSR_sum; TH1F * h_J1Pt_doubletagSB_sum; TH1F * h_J2Pt_doubletagSB_sum;
  TH1F * h_J1Pt_antitagSR_sum; TH1F * h_J2Pt_antitagSR_sum; TH1F * h_J1Pt_antitagSB_sum; TH1F * h_J2Pt_antitagSB_sum;
  TH1F * h_J1Pt_tagSR_sum; TH1F * h_J2Pt_tagSR_sum; TH1F * h_J1Pt_tagSB_sum; TH1F * h_J2Pt_tagSB_sum;

  TH1F * h_J1Pt_doubletagSR_QCD; TH1F * h_J2Pt_doubletagSR_QCD; TH1F * h_J1Pt_doubletagSB_QCD; TH1F * h_J2Pt_doubletagSB_QCD;
  TH1F * h_J1Pt_antitagSR_QCD; TH1F * h_J2Pt_antitagSR_QCD; TH1F * h_J1Pt_antitagSB_QCD; TH1F * h_J2Pt_antitagSB_QCD;
  TH1F * h_J1Pt_tagSR_QCD; TH1F * h_J2Pt_tagSR_QCD; TH1F * h_J1Pt_tagSB_QCD; TH1F * h_J2Pt_tagSB_QCD;

  TH1F * h_J1Pt_doubletagSR_GJets; TH1F * h_J2Pt_doubletagSR_GJets; TH1F * h_J1Pt_doubletagSB_GJets; TH1F * h_J2Pt_doubletagSB_GJets;
  TH1F * h_J1Pt_antitagSR_GJets; TH1F * h_J2Pt_antitagSR_GJets; TH1F * h_J1Pt_antitagSB_GJets; TH1F * h_J2Pt_antitagSB_GJets;
  TH1F * h_J1Pt_tagSR_GJets; TH1F * h_J2Pt_tagSR_GJets; TH1F * h_J1Pt_tagSB_GJets; TH1F * h_J2Pt_tagSB_GJets;

  TH1F * h_J1Pt_doubletagSR_TT; TH1F * h_J2Pt_doubletagSR_TT; TH1F * h_J1Pt_doubletagSB_TT; TH1F * h_J2Pt_doubletagSB_TT;
  TH1F * h_J1Pt_antitagSR_TT; TH1F * h_J2Pt_antitagSR_TT; TH1F * h_J1Pt_antitagSB_TT; TH1F * h_J2Pt_antitagSB_TT;
  TH1F * h_J1Pt_tagSR_TT; TH1F * h_J2Pt_tagSR_TT; TH1F * h_J1Pt_tagSB_TT; TH1F * h_J2Pt_tagSB_TT;

  TH1F * h_J1Pt_doubletagSR_WJets; TH1F * h_J2Pt_doubletagSR_WJets; TH1F * h_J1Pt_doubletagSB_WJets; TH1F * h_J2Pt_doubletagSB_WJets;
  TH1F * h_J1Pt_antitagSR_WJets; TH1F * h_J2Pt_antitagSR_WJets; TH1F * h_J1Pt_antitagSB_WJets; TH1F * h_J2Pt_antitagSB_WJets;
  TH1F * h_J1Pt_tagSR_WJets; TH1F * h_J2Pt_tagSR_WJets; TH1F * h_J1Pt_tagSB_WJets; TH1F * h_J2Pt_tagSB_WJets;

  TH1F * h_J1Pt_doubletagSR_ZJets; TH1F * h_J2Pt_doubletagSR_ZJets; TH1F * h_J1Pt_doubletagSB_ZJets; TH1F * h_J2Pt_doubletagSB_ZJets;
  TH1F * h_J1Pt_antitagSR_ZJets; TH1F * h_J2Pt_antitagSR_ZJets; TH1F * h_J1Pt_antitagSB_ZJets; TH1F * h_J2Pt_antitagSB_ZJets;
  TH1F * h_J1Pt_tagSR_ZJets; TH1F * h_J2Pt_tagSR_ZJets; TH1F * h_J1Pt_tagSB_ZJets; TH1F * h_J2Pt_tagSB_ZJets;

  TH1F * h_J1Pt_doubletagSR_T5HH1300; TH1F * h_J2Pt_doubletagSR_T5HH1300; TH1F * h_J1Pt_doubletagSB_T5HH1300; TH1F * h_J2Pt_doubletagSB_T5HH1300;
  TH1F * h_J1Pt_antitagSR_T5HH1300; TH1F * h_J2Pt_antitagSR_T5HH1300; TH1F * h_J1Pt_antitagSB_T5HH1300; TH1F * h_J2Pt_antitagSB_T5HH1300;
  TH1F * h_J1Pt_tagSR_T5HH1300; TH1F * h_J2Pt_tagSR_T5HH1300; TH1F * h_J1Pt_tagSB_T5HH1300; TH1F * h_J2Pt_tagSB_T5HH1300;

  TH1F * h_J1Pt_doubletagSR_T5HH1700; TH1F * h_J2Pt_doubletagSR_T5HH1700; TH1F * h_J1Pt_doubletagSB_T5HH1700; TH1F * h_J2Pt_doubletagSB_T5HH1700;
  TH1F * h_J1Pt_antitagSR_T5HH1700; TH1F * h_J2Pt_antitagSR_T5HH1700; TH1F * h_J1Pt_antitagSB_T5HH1700; TH1F * h_J2Pt_antitagSB_T5HH1700;
  TH1F * h_J1Pt_tagSR_T5HH1700; TH1F * h_J2Pt_tagSR_T5HH1700; TH1F * h_J1Pt_tagSB_T5HH1700; TH1F * h_J2Pt_tagSB_T5HH1700;

  TH1F * h_J1Pt_doubletagSR_T5HH2100; TH1F * h_J2Pt_doubletagSR_T5HH2100; TH1F * h_J1Pt_doubletagSB_T5HH2100; TH1F * h_J2Pt_doubletagSB_T5HH2100;
  TH1F * h_J1Pt_antitagSR_T5HH2100; TH1F * h_J2Pt_antitagSR_T5HH2100; TH1F * h_J1Pt_antitagSB_T5HH2100; TH1F * h_J2Pt_antitagSB_T5HH2100;
  TH1F * h_J1Pt_tagSR_T5HH2100; TH1F * h_J2Pt_tagSR_T5HH2100; TH1F * h_J1Pt_tagSB_T5HH2100; TH1F * h_J2Pt_tagSB_T5HH2100;

  //Double b
  TH1F * h_J1_DoubleBLowMET_doubletagSR_QCD; TH1F * h_J1_DoubleBHighMET_doubletagSR_QCD; TH1F * h_J1_DoubleBLowMET_doubletagSB_QCD; TH1F * h_J1_DoubleBHighMET_doubletagSB_QCD;
  TH1F * h_J1_DoubleBLowMET_tagSR_QCD; TH1F * h_J1_DoubleBHighMET_tagSR_QCD; TH1F * h_J1_DoubleBLowMET_tagSB_QCD; TH1F * h_J1_DoubleBHighMET_tagSB_QCD;
  TH1F * h_J1_DoubleBLowMET_antitagSR_QCD; TH1F * h_J1_DoubleBHighMET_antitagSR_QCD; TH1F * h_J1_DoubleBLowMET_antitagSB_QCD; TH1F * h_J1_DoubleBHighMET_antitagSB_QCD;

  TH1F * h_J1_DoubleBLowMET_doubletagSR_TT; TH1F * h_J1_DoubleBHighMET_doubletagSR_TT; TH1F * h_J1_DoubleBLowMET_doubletagSB_TT; TH1F * h_J1_DoubleBHighMET_doubletagSB_TT;
  TH1F * h_J1_DoubleBLowMET_tagSR_TT; TH1F * h_J1_DoubleBHighMET_tagSR_TT; TH1F * h_J1_DoubleBLowMET_tagSB_TT; TH1F * h_J1_DoubleBHighMET_tagSB_TT;
  TH1F * h_J1_DoubleBLowMET_antitagSR_TT; TH1F * h_J1_DoubleBHighMET_antitagSR_TT; TH1F * h_J1_DoubleBLowMET_antitagSB_TT; TH1F * h_J1_DoubleBHighMET_antitagSB_TT;

  TH1F * h_J1_DoubleBLowMET_doubletagSR_WJets; TH1F * h_J1_DoubleBHighMET_doubletagSR_WJets; TH1F * h_J1_DoubleBLowMET_doubletagSB_WJets; TH1F * h_J1_DoubleBHighMET_doubletagSB_WJets;
  TH1F * h_J1_DoubleBLowMET_tagSR_WJets; TH1F * h_J1_DoubleBHighMET_tagSR_WJets; TH1F * h_J1_DoubleBLowMET_tagSB_WJets; TH1F * h_J1_DoubleBHighMET_tagSB_WJets;
  TH1F * h_J1_DoubleBLowMET_antitagSR_WJets; TH1F * h_J1_DoubleBHighMET_antitagSR_WJets; TH1F * h_J1_DoubleBLowMET_antitagSB_WJets; TH1F * h_J1_DoubleBHighMET_antitagSB_WJets;

  TH1F * h_J1_DoubleBLowMET_doubletagSR_ZJets; TH1F * h_J1_DoubleBHighMET_doubletagSR_ZJets; TH1F * h_J1_DoubleBLowMET_doubletagSB_ZJets; TH1F * h_J1_DoubleBHighMET_doubletagSB_ZJets;
  TH1F * h_J1_DoubleBLowMET_tagSR_ZJets; TH1F * h_J1_DoubleBHighMET_tagSR_ZJets; TH1F * h_J1_DoubleBLowMET_tagSB_ZJets; TH1F * h_J1_DoubleBHighMET_tagSB_ZJets;
  TH1F * h_J1_DoubleBLowMET_antitagSR_ZJets; TH1F * h_J1_DoubleBHighMET_antitagSR_ZJets; TH1F * h_J1_DoubleBLowMET_antitagSB_ZJets; TH1F * h_J1_DoubleBHighMET_antitagSB_ZJets;

  TH1F * h_J1_DoubleBLowMET_doubletagSR_T5HH1300; TH1F * h_J1_DoubleBHighMET_doubletagSR_T5HH1300; TH1F * h_J1_DoubleBLowMET_doubletagSB_T5HH1300; TH1F * h_J1_DoubleBHighMET_doubletagSB_T5HH1300;
  TH1F * h_J1_DoubleBLowMET_tagSR_T5HH1300; TH1F * h_J1_DoubleBHighMET_tagSR_T5HH1300; TH1F * h_J1_DoubleBLowMET_tagSB_T5HH1300; TH1F * h_J1_DoubleBHighMET_tagSB_T5HH1300;
  TH1F * h_J1_DoubleBLowMET_antitagSR_T5HH1300; TH1F * h_J1_DoubleBHighMET_antitagSR_T5HH1300; TH1F * h_J1_DoubleBLowMET_antitagSB_T5HH1300; TH1F * h_J1_DoubleBHighMET_antitagSB_T5HH1300;

  TH1F * h_J1_DoubleBLowMET_doubletagSR_T5HH1700; TH1F * h_J1_DoubleBHighMET_doubletagSR_T5HH1700; TH1F * h_J1_DoubleBLowMET_doubletagSB_T5HH1700; TH1F * h_J1_DoubleBHighMET_doubletagSB_T5HH1700;
  TH1F * h_J1_DoubleBLowMET_tagSR_T5HH1700; TH1F * h_J1_DoubleBHighMET_tagSR_T5HH1700; TH1F * h_J1_DoubleBLowMET_tagSB_T5HH1700; TH1F * h_J1_DoubleBHighMET_tagSB_T5HH1700;
  TH1F * h_J1_DoubleBLowMET_antitagSR_T5HH1700; TH1F * h_J1_DoubleBHighMET_antitagSR_T5HH1700; TH1F * h_J1_DoubleBLowMET_antitagSB_T5HH1700; TH1F * h_J1_DoubleBHighMET_antitagSB_T5HH1700;

  TH1F * h_J1_DoubleBLowMET_doubletagSR_T5HH2100; TH1F * h_J1_DoubleBHighMET_doubletagSR_T5HH2100; TH1F * h_J1_DoubleBLowMET_doubletagSB_T5HH2100; TH1F * h_J1_DoubleBHighMET_doubletagSB_T5HH2100;
  TH1F * h_J1_DoubleBLowMET_tagSR_T5HH2100; TH1F * h_J1_DoubleBHighMET_tagSR_T5HH2100; TH1F * h_J1_DoubleBLowMET_tagSB_T5HH2100; TH1F * h_J1_DoubleBHighMET_tagSB_T5HH2100;
  TH1F * h_J1_DoubleBLowMET_antitagSR_T5HH2100; TH1F * h_J1_DoubleBHighMET_antitagSR_T5HH2100; TH1F * h_J1_DoubleBLowMET_antitagSB_T5HH2100; TH1F * h_J1_DoubleBHighMET_antitagSB_T5HH2100;


  TH1F * h_MET_all_ZJets; TH1F * h_MET_Single_ZJets; TH1F * h_MET_Double_ZJets;
  TH1F * h_MET_all_WJets; TH1F * h_MET_Single_WJets; TH1F * h_MET_Double_WJets;
  TH1F * h_MET_all_TT; TH1F * h_MET_Single_TT; TH1F * h_MET_Double_TT;
  TH1F * h_MET_all_QCD; TH1F * h_MET_Single_QCD; TH1F * h_MET_Double_QCD;
  TH1F * h_MET_all_T5HH1300; TH1F * h_MET_Single_T5HH1300; TH1F * h_MET_Double_T5HH1300;
  TH1F * h_MET_all_T5HH1700; TH1F * h_MET_Single_T5HH1700; TH1F * h_MET_Double_T5HH1700;
  TH1F * h_MET_all_T5HH2100; TH1F * h_MET_Single_T5HH2100; TH1F * h_MET_Double_T5HH2100;


  if (whichRegion=="signal") {
    h_A_sum = (TH1F*)f->Get("MET_doubletagSR_sum"); h_B_sum = (TH1F*)f->Get("MET_doubletagSB_sum");
    h_A1_sum = (TH1F*)f->Get("MET_tagSR_sum"); h_B1_sum = (TH1F*)f->Get("MET_tagSB_sum");
    h_C_sum = (TH1F*)f->Get("MET_antitagSR_sum"); h_D_sum = (TH1F*)f->Get("MET_antitagSB_sum");

    h_A_TT = (TH1F*)f->Get("MET_doubletagSR_TT"); h_B_TT = (TH1F*)f->Get("MET_doubletagSB_TT");
    h_A1_TT = (TH1F*)f->Get("MET_tagSR_TT"); h_B1_TT = (TH1F*)f->Get("MET_tagSB_TT");
    h_C_TT = (TH1F*)f->Get("MET_antitagSR_TT"); h_D_TT = (TH1F*)f->Get("MET_antitagSB_TT");

    h_A_WJets = (TH1F*)f->Get("MET_doubletagSR_WJets"); h_B_WJets = (TH1F*)f->Get("MET_doubletagSB_WJets");
    h_A1_WJets = (TH1F*)f->Get("MET_tagSR_WJets"); h_B1_WJets = (TH1F*)f->Get("MET_tagSB_WJets");
    h_C_WJets = (TH1F*)f->Get("MET_antitagSR_WJets");   h_D_WJets = (TH1F*)f->Get("MET_antitagSB_WJets");

    h_A_ZJets = (TH1F*)f->Get("MET_doubletagSR_ZJets"); h_B_ZJets = (TH1F*)f->Get("MET_doubletagSB_ZJets");
    h_A1_ZJets = (TH1F*)f->Get("MET_tagSR_ZJets"); h_B1_ZJets = (TH1F*)f->Get("MET_tagSB_ZJets");
    h_C_ZJets = (TH1F*)f->Get("MET_antitagSR_ZJets");   h_D_ZJets = (TH1F*)f->Get("MET_antitagSB_ZJets");

    h_A_QCD = (TH1F*)f->Get("MET_doubletagSR_QCD"); h_B_QCD = (TH1F*)f->Get("MET_doubletagSB_QCD");
    h_A1_QCD = (TH1F*)f->Get("MET_tagSR_QCD"); h_B1_QCD = (TH1F*)f->Get("MET_tagSB_QCD");
    h_C_QCD = (TH1F*)f->Get("MET_antitagSR_QCD"); h_D_QCD = (TH1F*)f->Get("MET_antitagSB_QCD");


    h_A_T5HH1300 = (TH1F*)f->Get("MET_doubletagSR_T5HH1300"); h_B_T5HH1300 = (TH1F*)f->Get("MET_doubletagSB_T5HH1300");
    h_A1_T5HH1300 = (TH1F*)f->Get("MET_tagSR_T5HH1300"); h_B1_T5HH1300 = (TH1F*)f->Get("MET_tagSB_T5HH1300");
    h_C_T5HH1300 = (TH1F*)f->Get("MET_antitagSR_T5HH1300"); h_D_T5HH1300 = (TH1F*)f->Get("MET_antitagSB_T5HH1300");

    h_A_T5HH1700 = (TH1F*)f->Get("MET_doubletagSR_T5HH1700"); h_B_T5HH1700 = (TH1F*)f->Get("MET_doubletagSB_T5HH1700");
    h_A1_T5HH1700 = (TH1F*)f->Get("MET_tagSR_T5HH1700"); h_B1_T5HH1700 = (TH1F*)f->Get("MET_tagSB_T5HH1700");
    h_C_T5HH1700 = (TH1F*)f->Get("MET_antitagSR_T5HH1700"); h_D_T5HH1700 = (TH1F*)f->Get("MET_antitagSB_T5HH1700");

    h_A_T5HH2100 = (TH1F*)f->Get("MET_doubletagSR_T5HH2100"); h_B_T5HH2100 = (TH1F*)f->Get("MET_doubletagSB_T5HH2100");
    h_A1_T5HH2100 = (TH1F*)f->Get("MET_tagSR_T5HH2100"); h_B1_T5HH2100 = (TH1F*)f->Get("MET_tagSB_T5HH2100");
    h_C_T5HH2100 = (TH1F*)f->Get("MET_antitagSR_T5HH2100"); h_D_T5HH2100 = (TH1F*)f->Get("MET_antitagSB_T5HH2100");



    h_J1M_doubletagSR_sum = (TH1F*)f->Get("J1_m_doubletagSR_sum"); h_J2M_doubletagSR_sum = (TH1F*)f->Get("J2_m_doubletagSR_sum");
    h_J1M_doubletagSB_sum = (TH1F*)f->Get("J1_m_doubletagSB_sum"); h_J2M_doubletagSB_sum = (TH1F*)f->Get("J2_m_doubletagSB_sum");
    h_J1M_antitagSR_sum = (TH1F*)f->Get("J1_m_antitagSR_sum"); h_J2M_antitagSR_sum = (TH1F*)f->Get("J2_m_antitagSR_sum");
    h_J1M_antitagSB_sum = (TH1F*)f->Get("J1_m_antitagSB_sum"); h_J2M_antitagSB_sum = (TH1F*)f->Get("J2_m_antitagSB_sum");
    h_J1M_tagSR_sum = (TH1F*)f->Get("J1_m_tagSR_sum"); h_J2M_tagSR_sum = (TH1F*)f->Get("J2_m_tagSR_sum");
    h_J1M_tagSB_sum = (TH1F*)f->Get("J1_m_tagSB_sum"); h_J2M_tagSB_sum = (TH1F*)f->Get("J2_m_tagSB_sum");

    h_J1M_doubletagSR_ZJets = (TH1F*)f->Get("J1_m_doubletagSR_ZJets"); h_J2M_doubletagSR_ZJets = (TH1F*)f->Get("J2_m_doubletagSR_ZJets");
    h_J1M_doubletagSB_ZJets = (TH1F*)f->Get("J1_m_doubletagSB_ZJets"); h_J2M_doubletagSB_ZJets = (TH1F*)f->Get("J2_m_doubletagSB_ZJets");
    h_J1M_antitagSR_ZJets = (TH1F*)f->Get("J1_m_antitagSR_ZJets"); h_J2M_antitagSR_ZJets = (TH1F*)f->Get("J2_m_antitagSR_ZJets");
    h_J1M_antitagSB_ZJets = (TH1F*)f->Get("J1_m_antitagSB_ZJets"); h_J2M_antitagSB_ZJets = (TH1F*)f->Get("J2_m_antitagSB_ZJets");
    h_J1M_tagSR_ZJets = (TH1F*)f->Get("J1_m_tagSR_ZJets"); h_J2M_tagSR_ZJets = (TH1F*)f->Get("J2_m_tagSR_ZJets");
    h_J1M_tagSB_ZJets = (TH1F*)f->Get("J1_m_tagSB_ZJets"); h_J2M_tagSB_ZJets = (TH1F*)f->Get("J2_m_tagSB_ZJets");

    h_J1M_doubletagSR_QCD = (TH1F*)f->Get("J1_m_doubletagSR_QCD"); h_J2M_doubletagSR_QCD = (TH1F*)f->Get("J2_m_doubletagSR_QCD");
    h_J1M_doubletagSB_QCD = (TH1F*)f->Get("J1_m_doubletagSB_QCD"); h_J2M_doubletagSB_QCD = (TH1F*)f->Get("J2_m_doubletagSB_QCD");
    h_J1M_antitagSR_QCD = (TH1F*)f->Get("J1_m_antitagSR_QCD"); h_J2M_antitagSR_QCD = (TH1F*)f->Get("J2_m_antitagSR_QCD");
    h_J1M_antitagSB_QCD = (TH1F*)f->Get("J1_m_antitagSB_QCD"); h_J2M_antitagSB_QCD = (TH1F*)f->Get("J2_m_antitagSB_QCD");
    h_J1M_tagSR_QCD = (TH1F*)f->Get("J1_m_tagSR_QCD"); h_J2M_tagSR_QCD = (TH1F*)f->Get("J2_m_tagSR_QCD");
    h_J1M_tagSB_QCD = (TH1F*)f->Get("J1_m_tagSB_QCD"); h_J2M_tagSB_QCD = (TH1F*)f->Get("J2_m_tagSB_QCD");

    h_J1M_doubletagSR_TT = (TH1F*)f->Get("J1_m_doubletagSR_TT"); h_J2M_doubletagSR_TT = (TH1F*)f->Get("J2_m_doubletagSR_TT");
    h_J1M_doubletagSB_TT = (TH1F*)f->Get("J1_m_doubletagSB_TT"); h_J2M_doubletagSB_TT = (TH1F*)f->Get("J2_m_doubletagSB_TT");
    h_J1M_antitagSR_TT = (TH1F*)f->Get("J1_m_antitagSR_TT"); h_J2M_antitagSR_TT = (TH1F*)f->Get("J2_m_antitagSR_TT");
    h_J1M_antitagSB_TT = (TH1F*)f->Get("J1_m_antitagSB_TT"); h_J2M_antitagSB_TT = (TH1F*)f->Get("J2_m_antitagSB_TT");
    h_J1M_tagSR_TT = (TH1F*)f->Get("J1_m_tagSR_TT"); h_J2M_tagSR_TT = (TH1F*)f->Get("J2_m_tagSR_TT");
    h_J1M_tagSB_TT = (TH1F*)f->Get("J1_m_tagSB_TT"); h_J2M_tagSB_TT = (TH1F*)f->Get("J2_m_tagSB_TT");

    h_J1M_doubletagSR_WJets = (TH1F*)f->Get("J1_m_doubletagSR_WJets"); h_J2M_doubletagSR_WJets = (TH1F*)f->Get("J2_m_doubletagSR_WJets");
    h_J1M_doubletagSB_WJets = (TH1F*)f->Get("J1_m_doubletagSB_WJets"); h_J2M_doubletagSB_WJets = (TH1F*)f->Get("J2_m_doubletagSB_WJets");
    h_J1M_antitagSR_WJets = (TH1F*)f->Get("J1_m_antitagSR_WJets"); h_J2M_antitagSR_WJets = (TH1F*)f->Get("J2_m_antitagSR_WJets");
    h_J1M_antitagSB_WJets = (TH1F*)f->Get("J1_m_antitagSB_WJets"); h_J2M_antitagSB_WJets = (TH1F*)f->Get("J2_m_antitagSB_WJets");
    h_J1M_tagSR_WJets = (TH1F*)f->Get("J1_m_tagSR_WJets"); h_J2M_tagSR_WJets = (TH1F*)f->Get("J2_m_tagSR_WJets");
    h_J1M_tagSB_WJets = (TH1F*)f->Get("J1_m_tagSB_WJets"); h_J2M_tagSB_WJets = (TH1F*)f->Get("J2_m_tagSB_WJets");


    //for signal
    h_J1M_doubletagSR_T5HH1300 = (TH1F*)f->Get("J1_m_doubletagSR_T5HH1300"); h_J2M_doubletagSR_T5HH1300 = (TH1F*)f->Get("J2_m_doubletagSR_T5HH1300");
    h_J1M_doubletagSB_T5HH1300 = (TH1F*)f->Get("J1_m_doubletagSB_T5HH1300"); h_J2M_doubletagSB_T5HH1300 = (TH1F*)f->Get("J2_m_doubletagSB_T5HH1300");
    h_J1M_antitagSR_T5HH1300 = (TH1F*)f->Get("J1_m_antitagSR_T5HH1300"); h_J2M_antitagSR_T5HH1300 = (TH1F*)f->Get("J2_m_antitagSR_T5HH1300");
    h_J1M_antitagSB_T5HH1300 = (TH1F*)f->Get("J1_m_antitagSB_T5HH1300"); h_J2M_antitagSB_T5HH1300 = (TH1F*)f->Get("J2_m_antitagSB_T5HH1300");
    h_J1M_tagSR_T5HH1300 = (TH1F*)f->Get("J1_m_tagSR_T5HH1300"); h_J2M_tagSR_T5HH1300 = (TH1F*)f->Get("J2_m_tagSR_T5HH1300");
    h_J1M_tagSB_T5HH1300 = (TH1F*)f->Get("J1_m_tagSB_T5HH1300"); h_J2M_tagSB_T5HH1300 = (TH1F*)f->Get("J2_m_tagSB_T5HH1300");

    h_J1M_doubletagSR_T5HH1700 = (TH1F*)f->Get("J1_m_doubletagSR_T5HH1700"); h_J2M_doubletagSR_T5HH1700 = (TH1F*)f->Get("J2_m_doubletagSR_T5HH1700");
    h_J1M_doubletagSB_T5HH1700 = (TH1F*)f->Get("J1_m_doubletagSB_T5HH1700"); h_J2M_doubletagSB_T5HH1700 = (TH1F*)f->Get("J2_m_doubletagSB_T5HH1700");
    h_J1M_antitagSR_T5HH1700 = (TH1F*)f->Get("J1_m_antitagSR_T5HH1700"); h_J2M_antitagSR_T5HH1700 = (TH1F*)f->Get("J2_m_antitagSR_T5HH1700");
    h_J1M_antitagSB_T5HH1700 = (TH1F*)f->Get("J1_m_antitagSB_T5HH1700"); h_J2M_antitagSB_T5HH1700 = (TH1F*)f->Get("J2_m_antitagSB_T5HH1700");
    h_J1M_tagSR_T5HH1700 = (TH1F*)f->Get("J1_m_tagSR_T5HH1700"); h_J2M_tagSR_T5HH1700 = (TH1F*)f->Get("J2_m_tagSR_T5HH1700");
    h_J1M_tagSB_T5HH1700 = (TH1F*)f->Get("J1_m_tagSB_T5HH1700"); h_J2M_tagSB_T5HH1700 = (TH1F*)f->Get("J2_m_tagSB_T5HH1700");

    h_J1M_doubletagSR_T5HH2100 = (TH1F*)f->Get("J1_m_doubletagSR_T5HH2100"); h_J2M_doubletagSR_T5HH2100 = (TH1F*)f->Get("J2_m_doubletagSR_T5HH2100");
    h_J1M_doubletagSB_T5HH2100 = (TH1F*)f->Get("J1_m_doubletagSB_T5HH2100"); h_J2M_doubletagSB_T5HH2100 = (TH1F*)f->Get("J2_m_doubletagSB_T5HH2100");
    h_J1M_antitagSR_T5HH2100 = (TH1F*)f->Get("J1_m_antitagSR_T5HH2100"); h_J2M_antitagSR_T5HH2100 = (TH1F*)f->Get("J2_m_antitagSR_T5HH2100");
    h_J1M_antitagSB_T5HH2100 = (TH1F*)f->Get("J1_m_antitagSB_T5HH2100"); h_J2M_antitagSB_T5HH2100 = (TH1F*)f->Get("J2_m_antitagSB_T5HH2100");
    h_J1M_tagSR_T5HH2100 = (TH1F*)f->Get("J1_m_tagSR_T5HH2100"); h_J2M_tagSR_T5HH2100 = (TH1F*)f->Get("J2_m_tagSR_T5HH2100");
    h_J1M_tagSB_T5HH2100 = (TH1F*)f->Get("J1_m_tagSB_T5HH2100"); h_J2M_tagSB_T5HH2100 = (TH1F*)f->Get("J2_m_tagSB_T5HH2100");



    //For pt
    h_J1Pt_doubletagSR_sum = (TH1F*)f->Get("J1_pt_doubletagSR_sum"); h_J2Pt_doubletagSR_sum = (TH1F*)f->Get("J2_pt_doubletagSR_sum");
    h_J1Pt_doubletagSB_sum = (TH1F*)f->Get("J1_pt_doubletagSB_sum"); h_J2Pt_doubletagSB_sum = (TH1F*)f->Get("J2_pt_doubletagSB_sum");
    h_J1Pt_antitagSR_sum = (TH1F*)f->Get("J1_pt_antitagSR_sum"); h_J2Pt_antitagSR_sum = (TH1F*)f->Get("J2_pt_antitagSR_sum");
    h_J1Pt_antitagSB_sum = (TH1F*)f->Get("J1_pt_antitagSB_sum"); h_J2Pt_antitagSB_sum = (TH1F*)f->Get("J2_pt_antitagSB_sum");
    h_J1Pt_tagSR_sum = (TH1F*)f->Get("J1_pt_tagSR_sum"); h_J2Pt_tagSR_sum = (TH1F*)f->Get("J2_pt_tagSR_sum");
    h_J1Pt_tagSB_sum = (TH1F*)f->Get("J1_pt_tagSB_sum"); h_J2Pt_tagSB_sum = (TH1F*)f->Get("J2_pt_tagSB_sum");

    h_J1Pt_doubletagSR_ZJets = (TH1F*)f->Get("J1_pt_doubletagSR_ZJets"); h_J2Pt_doubletagSR_ZJets = (TH1F*)f->Get("J2_pt_doubletagSR_ZJets");
    h_J1Pt_doubletagSB_ZJets = (TH1F*)f->Get("J1_pt_doubletagSB_ZJets"); h_J2Pt_doubletagSB_ZJets = (TH1F*)f->Get("J2_pt_doubletagSB_ZJets");
    h_J1Pt_antitagSR_ZJets = (TH1F*)f->Get("J1_pt_antitagSR_ZJets"); h_J2Pt_antitagSR_ZJets = (TH1F*)f->Get("J2_pt_antitagSR_ZJets");
    h_J1Pt_antitagSB_ZJets = (TH1F*)f->Get("J1_pt_antitagSB_ZJets"); h_J2Pt_antitagSB_ZJets = (TH1F*)f->Get("J2_pt_antitagSB_ZJets");
    h_J1Pt_tagSR_ZJets = (TH1F*)f->Get("J1_pt_tagSR_ZJets"); h_J2Pt_tagSR_ZJets = (TH1F*)f->Get("J2_pt_tagSR_ZJets");
    h_J1Pt_tagSB_ZJets = (TH1F*)f->Get("J1_pt_tagSB_ZJets"); h_J2Pt_tagSB_ZJets = (TH1F*)f->Get("J2_pt_tagSB_ZJets");

    h_J1Pt_doubletagSR_QCD = (TH1F*)f->Get("J1_pt_doubletagSR_QCD"); h_J2Pt_doubletagSR_QCD = (TH1F*)f->Get("J2_pt_doubletagSR_QCD");
    h_J1Pt_doubletagSB_QCD = (TH1F*)f->Get("J1_pt_doubletagSB_QCD"); h_J2Pt_doubletagSB_QCD = (TH1F*)f->Get("J2_pt_doubletagSB_QCD");
    h_J1Pt_antitagSR_QCD = (TH1F*)f->Get("J1_pt_antitagSR_QCD"); h_J2Pt_antitagSR_QCD = (TH1F*)f->Get("J2_pt_antitagSR_QCD");
    h_J1Pt_antitagSB_QCD = (TH1F*)f->Get("J1_pt_antitagSB_QCD"); h_J2Pt_antitagSB_QCD = (TH1F*)f->Get("J2_pt_antitagSB_QCD");
    h_J1Pt_tagSR_QCD = (TH1F*)f->Get("J1_pt_tagSR_QCD"); h_J2Pt_tagSR_QCD = (TH1F*)f->Get("J2_pt_tagSR_QCD");
    h_J1Pt_tagSB_QCD = (TH1F*)f->Get("J1_pt_tagSB_QCD"); h_J2Pt_tagSB_QCD = (TH1F*)f->Get("J2_pt_tagSB_QCD");

    h_J1Pt_doubletagSR_TT = (TH1F*)f->Get("J1_pt_doubletagSR_TT"); h_J2Pt_doubletagSR_TT = (TH1F*)f->Get("J2_pt_doubletagSR_TT");
    h_J1Pt_doubletagSB_TT = (TH1F*)f->Get("J1_pt_doubletagSB_TT"); h_J2Pt_doubletagSB_TT = (TH1F*)f->Get("J2_pt_doubletagSB_TT");
    h_J1Pt_antitagSR_TT = (TH1F*)f->Get("J1_pt_antitagSR_TT"); h_J2Pt_antitagSR_TT = (TH1F*)f->Get("J2_pt_antitagSR_TT");
    h_J1Pt_antitagSB_TT = (TH1F*)f->Get("J1_pt_antitagSB_TT"); h_J2Pt_antitagSB_TT = (TH1F*)f->Get("J2_pt_antitagSB_TT");
    h_J1Pt_tagSR_TT = (TH1F*)f->Get("J1_pt_tagSR_TT"); h_J2Pt_tagSR_TT = (TH1F*)f->Get("J2_pt_tagSR_TT");
    h_J1Pt_tagSB_TT = (TH1F*)f->Get("J1_pt_tagSB_TT"); h_J2Pt_tagSB_TT = (TH1F*)f->Get("J2_pt_tagSB_TT");

    h_J1Pt_doubletagSR_WJets = (TH1F*)f->Get("J1_pt_doubletagSR_WJets"); h_J2Pt_doubletagSR_WJets = (TH1F*)f->Get("J2_pt_doubletagSR_WJets");
    h_J1Pt_doubletagSB_WJets = (TH1F*)f->Get("J1_pt_doubletagSB_WJets"); h_J2Pt_doubletagSB_WJets = (TH1F*)f->Get("J2_pt_doubletagSB_WJets");
    h_J1Pt_antitagSR_WJets = (TH1F*)f->Get("J1_pt_antitagSR_WJets"); h_J2Pt_antitagSR_WJets = (TH1F*)f->Get("J2_pt_antitagSR_WJets");
    h_J1Pt_antitagSB_WJets = (TH1F*)f->Get("J1_pt_antitagSB_WJets"); h_J2Pt_antitagSB_WJets = (TH1F*)f->Get("J2_pt_antitagSB_WJets");
    h_J1Pt_tagSR_WJets = (TH1F*)f->Get("J1_pt_tagSR_WJets"); h_J2Pt_tagSR_WJets = (TH1F*)f->Get("J2_pt_tagSR_WJets");
    h_J1Pt_tagSB_WJets = (TH1F*)f->Get("J1_pt_tagSB_WJets"); h_J2Pt_tagSB_WJets = (TH1F*)f->Get("J2_pt_tagSB_WJets");

    //for signal
    h_J1Pt_doubletagSR_T5HH1300 = (TH1F*)f->Get("J1_pt_doubletagSR_T5HH1300"); h_J2Pt_doubletagSR_T5HH1300 = (TH1F*)f->Get("J2_pt_doubletagSR_T5HH1300");
    h_J1Pt_doubletagSB_T5HH1300 = (TH1F*)f->Get("J1_pt_doubletagSB_T5HH1300"); h_J2Pt_doubletagSB_T5HH1300 = (TH1F*)f->Get("J2_pt_doubletagSB_T5HH1300");
    h_J1Pt_antitagSR_T5HH1300 = (TH1F*)f->Get("J1_pt_antitagSR_T5HH1300"); h_J2Pt_antitagSR_T5HH1300 = (TH1F*)f->Get("J2_pt_antitagSR_T5HH1300");
    h_J1Pt_antitagSB_T5HH1300 = (TH1F*)f->Get("J1_pt_antitagSB_T5HH1300"); h_J2Pt_antitagSB_T5HH1300 = (TH1F*)f->Get("J2_pt_antitagSB_T5HH1300");
    h_J1Pt_tagSR_T5HH1300 = (TH1F*)f->Get("J1_pt_tagSR_T5HH1300"); h_J2Pt_tagSR_T5HH1300 = (TH1F*)f->Get("J2_pt_tagSR_T5HH1300");
    h_J1Pt_tagSB_T5HH1300 = (TH1F*)f->Get("J1_pt_tagSB_T5HH1300"); h_J2Pt_tagSB_T5HH1300 = (TH1F*)f->Get("J2_pt_tagSB_T5HH1300");

    h_J1Pt_doubletagSR_T5HH1700 = (TH1F*)f->Get("J1_pt_doubletagSR_T5HH1700"); h_J2Pt_doubletagSR_T5HH1700 = (TH1F*)f->Get("J2_pt_doubletagSR_T5HH1700");
    h_J1Pt_doubletagSB_T5HH1700 = (TH1F*)f->Get("J1_pt_doubletagSB_T5HH1700"); h_J2Pt_doubletagSB_T5HH1700 = (TH1F*)f->Get("J2_pt_doubletagSB_T5HH1700");
    h_J1Pt_antitagSR_T5HH1700 = (TH1F*)f->Get("J1_pt_antitagSR_T5HH1700"); h_J2Pt_antitagSR_T5HH1700 = (TH1F*)f->Get("J2_pt_antitagSR_T5HH1700");
    h_J1Pt_antitagSB_T5HH1700 = (TH1F*)f->Get("J1_pt_antitagSB_T5HH1700"); h_J2Pt_antitagSB_T5HH1700 = (TH1F*)f->Get("J2_pt_antitagSB_T5HH1700");
    h_J1Pt_tagSR_T5HH1700 = (TH1F*)f->Get("J1_pt_tagSR_T5HH1700"); h_J2Pt_tagSR_T5HH1700 = (TH1F*)f->Get("J2_pt_tagSR_T5HH1700");
    h_J1Pt_tagSB_T5HH1700 = (TH1F*)f->Get("J1_pt_tagSB_T5HH1700"); h_J2Pt_tagSB_T5HH1700 = (TH1F*)f->Get("J2_pt_tagSB_T5HH1700");

    h_J1Pt_doubletagSR_T5HH2100 = (TH1F*)f->Get("J1_pt_doubletagSR_T5HH2100"); h_J2Pt_doubletagSR_T5HH2100 = (TH1F*)f->Get("J2_pt_doubletagSR_T5HH2100");
    h_J1Pt_doubletagSB_T5HH2100 = (TH1F*)f->Get("J1_pt_doubletagSB_T5HH2100"); h_J2Pt_doubletagSB_T5HH2100 = (TH1F*)f->Get("J2_pt_doubletagSB_T5HH2100");
    h_J1Pt_antitagSR_T5HH2100 = (TH1F*)f->Get("J1_pt_antitagSR_T5HH2100"); h_J2Pt_antitagSR_T5HH2100 = (TH1F*)f->Get("J2_pt_antitagSR_T5HH2100");
    h_J1Pt_antitagSB_T5HH2100 = (TH1F*)f->Get("J1_pt_antitagSB_T5HH2100"); h_J2Pt_antitagSB_T5HH2100 = (TH1F*)f->Get("J2_pt_antitagSB_T5HH2100");
    h_J1Pt_tagSR_T5HH2100 = (TH1F*)f->Get("J1_pt_tagSR_T5HH2100"); h_J2Pt_tagSR_T5HH2100 = (TH1F*)f->Get("J2_pt_tagSR_T5HH2100");
    h_J1Pt_tagSB_T5HH2100 = (TH1F*)f->Get("J1_pt_tagSB_T5HH2100"); h_J2Pt_tagSB_T5HH2100 = (TH1F*)f->Get("J2_pt_tagSB_T5HH2100");



    h_J1_DoubleBLowMET_doubletagSR_QCD = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSR_QCD"); h_J1_DoubleBHighMET_doubletagSR_QCD = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSR_QCD");
    h_J1_DoubleBLowMET_doubletagSB_QCD = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSB_QCD"); h_J1_DoubleBHighMET_doubletagSB_QCD = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSB_QCD");
    h_J1_DoubleBLowMET_tagSR_QCD = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSR_QCD"); h_J1_DoubleBHighMET_tagSR_QCD = (TH1F*)f->Get("J1_DoubleB_highMET_tagSR_QCD");
    h_J1_DoubleBLowMET_tagSB_QCD = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSB_QCD"); h_J1_DoubleBHighMET_tagSB_QCD = (TH1F*)f->Get("J1_DoubleB_highMET_tagSB_QCD");
    h_J1_DoubleBLowMET_antitagSR_QCD = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSR_QCD"); h_J1_DoubleBHighMET_antitagSR_QCD = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSR_QCD");
    h_J1_DoubleBLowMET_antitagSB_QCD = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSB_QCD"); h_J1_DoubleBHighMET_antitagSB_QCD = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSB_QCD");


    h_J1_DoubleBLowMET_doubletagSR_TT = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSR_TT"); h_J1_DoubleBHighMET_doubletagSR_TT = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSR_TT");
    h_J1_DoubleBLowMET_doubletagSB_TT = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSB_TT"); h_J1_DoubleBHighMET_doubletagSB_TT = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSB_TT");
    h_J1_DoubleBLowMET_tagSR_TT = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSR_TT"); h_J1_DoubleBHighMET_tagSR_TT = (TH1F*)f->Get("J1_DoubleB_highMET_tagSR_TT");
    h_J1_DoubleBLowMET_tagSB_TT = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSB_TT"); h_J1_DoubleBHighMET_tagSB_TT = (TH1F*)f->Get("J1_DoubleB_highMET_tagSB_TT");
    h_J1_DoubleBLowMET_antitagSR_TT = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSR_TT"); h_J1_DoubleBHighMET_antitagSR_TT = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSR_TT");
    h_J1_DoubleBLowMET_antitagSB_TT = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSB_TT"); h_J1_DoubleBHighMET_antitagSB_TT = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSB_TT");


    h_J1_DoubleBLowMET_doubletagSR_WJets = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSR_WJets"); h_J1_DoubleBHighMET_doubletagSR_WJets = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSR_WJets");
    h_J1_DoubleBLowMET_doubletagSB_WJets = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSB_WJets"); h_J1_DoubleBHighMET_doubletagSB_WJets = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSB_WJets");
    h_J1_DoubleBLowMET_tagSR_WJets = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSR_WJets"); h_J1_DoubleBHighMET_tagSR_WJets = (TH1F*)f->Get("J1_DoubleB_highMET_tagSR_WJets");
    h_J1_DoubleBLowMET_tagSB_WJets = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSB_WJets"); h_J1_DoubleBHighMET_tagSB_WJets = (TH1F*)f->Get("J1_DoubleB_highMET_tagSB_WJets");
    h_J1_DoubleBLowMET_antitagSR_WJets = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSR_WJets"); h_J1_DoubleBHighMET_antitagSR_WJets = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSR_WJets");
    h_J1_DoubleBLowMET_antitagSB_WJets = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSB_WJets"); h_J1_DoubleBHighMET_antitagSB_WJets = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSB_WJets");


    h_J1_DoubleBLowMET_doubletagSR_ZJets = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSR_ZJets"); h_J1_DoubleBHighMET_doubletagSR_ZJets = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSR_ZJets");
    h_J1_DoubleBLowMET_doubletagSB_ZJets = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSB_ZJets"); h_J1_DoubleBHighMET_doubletagSB_ZJets = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSB_ZJets");
    h_J1_DoubleBLowMET_tagSR_ZJets = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSR_ZJets"); h_J1_DoubleBHighMET_tagSR_ZJets = (TH1F*)f->Get("J1_DoubleB_highMET_tagSR_ZJets");
    h_J1_DoubleBLowMET_tagSB_ZJets = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSB_ZJets"); h_J1_DoubleBHighMET_tagSB_ZJets = (TH1F*)f->Get("J1_DoubleB_highMET_tagSB_ZJets");
    h_J1_DoubleBLowMET_antitagSR_ZJets = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSR_ZJets"); h_J1_DoubleBHighMET_antitagSR_ZJets = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSR_ZJets");
    h_J1_DoubleBLowMET_antitagSB_ZJets = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSB_ZJets"); h_J1_DoubleBHighMET_antitagSB_ZJets = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSB_ZJets");


    h_J1_DoubleBLowMET_doubletagSR_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSR_T5HH1300"); h_J1_DoubleBHighMET_doubletagSR_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSR_T5HH1300");
    h_J1_DoubleBLowMET_doubletagSB_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSB_T5HH1300"); h_J1_DoubleBHighMET_doubletagSB_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSB_T5HH1300");
    h_J1_DoubleBLowMET_tagSR_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSR_T5HH1300"); h_J1_DoubleBHighMET_tagSR_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_highMET_tagSR_T5HH1300");
    h_J1_DoubleBLowMET_tagSB_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSB_T5HH1300"); h_J1_DoubleBHighMET_tagSB_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_highMET_tagSB_T5HH1300");
    h_J1_DoubleBLowMET_antitagSR_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSR_T5HH1300"); h_J1_DoubleBHighMET_antitagSR_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSR_T5HH1300");
    h_J1_DoubleBLowMET_antitagSB_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSB_T5HH1300"); h_J1_DoubleBHighMET_antitagSB_T5HH1300 = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSB_T5HH1300");


    h_J1_DoubleBLowMET_doubletagSR_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSR_T5HH1700"); h_J1_DoubleBHighMET_doubletagSR_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSR_T5HH1700");
    h_J1_DoubleBLowMET_doubletagSB_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSB_T5HH1700"); h_J1_DoubleBHighMET_doubletagSB_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSB_T5HH1700");
    h_J1_DoubleBLowMET_tagSR_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSR_T5HH1700"); h_J1_DoubleBHighMET_tagSR_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_highMET_tagSR_T5HH1700");
    h_J1_DoubleBLowMET_tagSB_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSB_T5HH1700"); h_J1_DoubleBHighMET_tagSB_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_highMET_tagSB_T5HH1700");
    h_J1_DoubleBLowMET_antitagSR_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSR_T5HH1700"); h_J1_DoubleBHighMET_antitagSR_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSR_T5HH1700");
    h_J1_DoubleBLowMET_antitagSB_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSB_T5HH1700"); h_J1_DoubleBHighMET_antitagSB_T5HH1700 = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSB_T5HH1700");


    h_J1_DoubleBLowMET_doubletagSR_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSR_T5HH2100"); h_J1_DoubleBHighMET_doubletagSR_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSR_T5HH2100");
    h_J1_DoubleBLowMET_doubletagSB_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_lowMET_doubletagSB_T5HH2100"); h_J1_DoubleBHighMET_doubletagSB_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_highMET_doubletagSB_T5HH2100");
    h_J1_DoubleBLowMET_tagSR_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSR_T5HH2100"); h_J1_DoubleBHighMET_tagSR_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_highMET_tagSR_T5HH2100");
    h_J1_DoubleBLowMET_tagSB_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_lowMET_tagSB_T5HH2100"); h_J1_DoubleBHighMET_tagSB_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_highMET_tagSB_T5HH2100");
    h_J1_DoubleBLowMET_antitagSR_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSR_T5HH2100"); h_J1_DoubleBHighMET_antitagSR_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSR_T5HH2100");
    h_J1_DoubleBLowMET_antitagSB_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_lowMET_antitagSB_T5HH2100"); h_J1_DoubleBHighMET_antitagSB_T5HH2100 = (TH1F*)f->Get("J1_DoubleB_highMET_antitagSB_T5HH2100");


    // //MET plots
    // h_MET_all_ZJets = (TH1F*)f->Get("MET_All_ZJets"); h_MET_Single_ZJets = (TH1F*)f->Get("MET_Single_ZJets"); h_MET_Double_ZJets = (TH1F*)f->Get("MET_Double_ZJets");
    // h_MET_all_WJets = (TH1F*)f->Get("MET_All_WJets"); h_MET_Single_WJets = (TH1F*)f->Get("MET_Single_WJets"); h_MET_Double_WJets = (TH1F*)f->Get("MET_Double_WJets");
    // h_MET_all_TT = (TH1F*)f->Get("MET_All_TT"); h_MET_Single_TT = (TH1F*)f->Get("MET_Single_TT"); h_MET_Double_TT = (TH1F*)f->Get("MET_Double_TT");
    // h_MET_all_QCD = (TH1F*)f->Get("MET_All_QCD"); h_MET_Single_QCD = (TH1F*)f->Get("MET_Single_QCD"); h_MET_Double_QCD = (TH1F*)f->Get("MET_Double_QCD");
    // h_MET_all_T5HH1300 = (TH1F*)f->Get("MET_All_T5HH1300"); h_MET_Single_T5HH1300 = (TH1F*)f->Get("MET_Single_T5HH1300"); h_MET_Double_T5HH1300 = (TH1F*)f->Get("MET_Double_T5HH1300");
    // h_MET_all_T5HH1700 = (TH1F*)f->Get("MET_All_T5HH1700"); h_MET_Single_T5HH1700 = (TH1F*)f->Get("MET_Single_T5HH1700"); h_MET_Double_T5HH1700 = (TH1F*)f->Get("MET_Double_T5HH1700");
    // h_MET_all_T5HH2100 = (TH1F*)f->Get("MET_All_T5HH2100"); h_MET_Single_T5HH2100 = (TH1F*)f->Get("MET_Single_T5HH2100"); h_MET_Double_T5HH2100 = (TH1F*)f->Get("MET_Double_T5HH2100");
  }



  vector<TH1F*> histos_ABCD_sum = {h_A_sum, h_B_sum, h_C_sum, h_D_sum};
  vector<TH1F*> histos_ABCD_QCD = {h_A_QCD, h_B_QCD, h_C_QCD, h_D_QCD};
  vector<TH1F*> histos_ABCD_GJets = {h_A_GJets, h_B_GJets, h_C_GJets, h_D_GJets};
  vector<TH1F*> histos_ABCD_TT = {h_A_TT, h_B_TT, h_C_TT, h_D_TT};
  vector<TH1F*> histos_ABCD_WJets = {h_A_WJets, h_B_WJets, h_C_WJets, h_D_WJets};
  vector<TH1F*> histos_ABCD_ZJets = {h_A_ZJets, h_B_ZJets, h_C_ZJets, h_D_ZJets};

  vector<TH1F*> histos_A1B1CD_sum = {h_A1_sum, h_B1_sum, h_C_sum, h_D_sum};
  vector<TH1F*> histos_A1B1CD_QCD = {h_A1_QCD, h_B1_QCD, h_C_QCD, h_D_QCD};
  vector<TH1F*> histos_A1B1CD_GJets = {h_A1_GJets, h_B1_GJets, h_C_GJets, h_D_GJets};
  vector<TH1F*> histos_A1B1CD_TT = {h_A1_TT, h_B1_TT, h_C_TT, h_D_TT};
  vector<TH1F*> histos_A1B1CD_WJets = {h_A1_WJets, h_B1_WJets, h_C_WJets, h_D_WJets};
  vector<TH1F*> histos_A1B1CD_ZJets = {h_A1_ZJets, h_B1_ZJets, h_C_ZJets, h_D_ZJets};


  vector<TH1F*> histos_Rpf_J1_sum = {h_J1M_doubletagSR_sum, h_J1M_doubletagSB_sum, h_J1M_antitagSR_sum, h_J1M_antitagSB_sum};
  vector<TH1F*> histos_Rpf_J2_sum = {h_J2M_doubletagSR_sum, h_J2M_doubletagSB_sum, h_J2M_antitagSR_sum, h_J2M_antitagSB_sum};
  vector<TH1F*> histos_Rpfsingle_J1_sum = {h_J1M_tagSR_sum, h_J1M_tagSB_sum, h_J1M_antitagSR_sum, h_J1M_antitagSB_sum};
  vector<TH1F*> histos_Rpfsingle_J2_sum = {h_J2M_tagSR_sum, h_J2M_tagSB_sum, h_J2M_antitagSR_sum, h_J2M_antitagSB_sum};

  vector<TH1F*> histos_Rpf_J1_QCD = {h_J1M_doubletagSR_QCD, h_J1M_doubletagSB_QCD, h_J1M_antitagSR_QCD, h_J1M_antitagSB_QCD};
  vector<TH1F*> histos_Rpf_J2_QCD = {h_J2M_doubletagSR_QCD, h_J2M_doubletagSB_QCD, h_J2M_antitagSR_QCD, h_J2M_antitagSB_QCD};
  vector<TH1F*> histos_Rpfsingle_J1_QCD = {h_J1M_tagSR_QCD, h_J1M_tagSB_QCD, h_J1M_antitagSR_QCD, h_J1M_antitagSB_QCD};
  vector<TH1F*> histos_Rpfsingle_J2_QCD = {h_J2M_tagSR_QCD, h_J2M_tagSB_QCD, h_J2M_antitagSR_QCD, h_J2M_antitagSB_QCD};

  vector<TH1F*> histos_Rpf_J1_GJets = {h_J1M_doubletagSR_GJets, h_J1M_doubletagSB_GJets, h_J1M_antitagSR_GJets, h_J1M_antitagSB_GJets};
  vector<TH1F*> histos_Rpf_J2_GJets = {h_J2M_doubletagSR_GJets, h_J2M_doubletagSB_GJets, h_J2M_antitagSR_GJets, h_J2M_antitagSB_GJets};
  vector<TH1F*> histos_Rpfsingle_J1_GJets = {h_J1M_tagSR_GJets, h_J1M_tagSB_GJets, h_J1M_antitagSR_GJets, h_J1M_antitagSB_GJets};
  vector<TH1F*> histos_Rpfsingle_J2_GJets = {h_J2M_tagSR_GJets, h_J2M_tagSB_GJets, h_J2M_antitagSR_GJets, h_J2M_antitagSB_GJets};

  vector<TH1F*> histos_Rpf_J1_TT = {h_J1M_doubletagSR_TT, h_J1M_doubletagSB_TT, h_J1M_antitagSR_TT, h_J1M_antitagSB_TT};
  vector<TH1F*> histos_Rpf_J2_TT = {h_J2M_doubletagSR_TT, h_J2M_doubletagSB_TT, h_J2M_antitagSR_TT, h_J2M_antitagSB_TT};
  vector<TH1F*> histos_Rpfsingle_J1_TT = {h_J1M_tagSR_TT, h_J1M_tagSB_TT, h_J1M_antitagSR_TT, h_J1M_antitagSB_TT};
  vector<TH1F*> histos_Rpfsingle_J2_TT = {h_J2M_tagSR_TT, h_J2M_tagSB_TT, h_J2M_antitagSR_TT, h_J2M_antitagSB_TT};

  vector<TH1F*> histos_Rpf_J1_WJets = {h_J1M_doubletagSR_WJets, h_J1M_doubletagSB_WJets, h_J1M_antitagSR_WJets, h_J1M_antitagSB_WJets};
  vector<TH1F*> histos_Rpf_J2_WJets = {h_J2M_doubletagSR_WJets, h_J2M_doubletagSB_WJets, h_J2M_antitagSR_WJets, h_J2M_antitagSB_WJets};
  vector<TH1F*> histos_Rpfsingle_J1_WJets = {h_J1M_tagSR_WJets, h_J1M_tagSB_WJets, h_J1M_antitagSR_WJets, h_J1M_antitagSB_WJets};
  vector<TH1F*> histos_Rpfsingle_J2_WJets = {h_J2M_tagSR_WJets, h_J2M_tagSB_WJets, h_J2M_antitagSR_WJets, h_J2M_antitagSB_WJets};

  vector<TH1F*> histos_Rpf_J1_ZJets = {h_J1M_doubletagSR_ZJets, h_J1M_doubletagSB_ZJets, h_J1M_antitagSR_ZJets, h_J1M_antitagSB_ZJets};
  vector<TH1F*> histos_Rpf_J2_ZJets = {h_J2M_doubletagSR_ZJets, h_J2M_doubletagSB_ZJets, h_J2M_antitagSR_ZJets, h_J2M_antitagSB_ZJets};
  vector<TH1F*> histos_Rpfsingle_J1_ZJets = {h_J1M_tagSR_ZJets, h_J1M_tagSB_ZJets, h_J1M_antitagSR_ZJets, h_J1M_antitagSB_ZJets};
  vector<TH1F*> histos_Rpfsingle_J2_ZJets = {h_J2M_tagSR_ZJets, h_J2M_tagSB_ZJets, h_J2M_antitagSR_ZJets, h_J2M_antitagSB_ZJets};



  vector<TH1F*> histos_j1mass_ZJets = {h_J1M_doubletagSR_ZJets,h_J1M_doubletagSB_ZJets,h_J1M_tagSR_ZJets,h_J1M_tagSB_ZJets, h_J1M_antitagSR_ZJets, h_J1M_antitagSB_ZJets};
  vector<TH1F*> histos_j1mass_WJets = {h_J1M_doubletagSR_WJets,h_J1M_doubletagSB_WJets,h_J1M_tagSR_WJets,h_J1M_tagSB_WJets, h_J1M_antitagSR_WJets, h_J1M_antitagSB_WJets};
  vector<TH1F*> histos_j1mass_TT    = {h_J1M_doubletagSR_TT,h_J1M_doubletagSB_TT,h_J1M_tagSR_TT,h_J1M_tagSB_TT, h_J1M_antitagSR_TT, h_J1M_antitagSB_TT};
  vector<TH1F*> histos_j1mass_QCD   = {h_J1M_doubletagSR_QCD,h_J1M_doubletagSB_QCD,h_J1M_tagSR_QCD,h_J1M_tagSB_QCD, h_J1M_antitagSR_QCD, h_J1M_antitagSB_QCD};
  vector<TH1F*> histos_j1mass_T5HH1300 = {h_J1M_doubletagSR_T5HH1300,h_J1M_doubletagSB_T5HH1300,h_J1M_tagSR_T5HH1300,h_J1M_tagSB_T5HH1300, h_J1M_antitagSR_T5HH1300, h_J1M_antitagSB_T5HH1300};
  vector<TH1F*> histos_j1mass_T5HH1700 = {h_J1M_doubletagSR_T5HH1700,h_J1M_doubletagSB_T5HH1700,h_J1M_tagSR_T5HH1700,h_J1M_tagSB_T5HH1700, h_J1M_antitagSR_T5HH1700, h_J1M_antitagSB_T5HH1700};
  vector<TH1F*> histos_j1mass_T5HH2100 = {h_J1M_doubletagSR_T5HH2100,h_J1M_doubletagSB_T5HH2100,h_J1M_tagSR_T5HH2100,h_J1M_tagSB_T5HH2100, h_J1M_antitagSR_T5HH2100, h_J1M_antitagSB_T5HH2100};

  vector<TH1F*> histos_j2mass_ZJets = {h_J2M_doubletagSR_ZJets,h_J2M_doubletagSB_ZJets,h_J2M_tagSR_ZJets,h_J2M_tagSB_ZJets, h_J2M_antitagSR_ZJets, h_J2M_antitagSB_ZJets};
  vector<TH1F*> histos_j2mass_WJets = {h_J2M_doubletagSR_WJets,h_J2M_doubletagSB_WJets,h_J2M_tagSR_WJets,h_J2M_tagSB_WJets, h_J2M_antitagSR_WJets, h_J2M_antitagSB_WJets};
  vector<TH1F*> histos_j2mass_TT    = {h_J2M_doubletagSR_TT,h_J2M_doubletagSB_TT,h_J2M_tagSR_TT,h_J2M_tagSB_TT, h_J2M_antitagSR_TT, h_J2M_antitagSB_TT};
  vector<TH1F*> histos_j2mass_QCD   = {h_J2M_doubletagSR_QCD,h_J2M_doubletagSB_QCD,h_J2M_tagSR_QCD,h_J2M_tagSB_QCD, h_J2M_antitagSR_QCD, h_J2M_antitagSB_QCD};
  vector<TH1F*> histos_j2mass_T5HH1300 = {h_J2M_doubletagSR_T5HH1300,h_J2M_doubletagSB_T5HH1300,h_J2M_tagSR_T5HH1300,h_J2M_tagSB_T5HH1300, h_J2M_antitagSR_T5HH1300, h_J2M_antitagSB_T5HH1300};
  vector<TH1F*> histos_j2mass_T5HH1700 = {h_J2M_doubletagSR_T5HH1700,h_J2M_doubletagSB_T5HH1700,h_J2M_tagSR_T5HH1700,h_J2M_tagSB_T5HH1700, h_J2M_antitagSR_T5HH1700, h_J2M_antitagSB_T5HH1700};
  vector<TH1F*> histos_j2mass_T5HH2100 = {h_J2M_doubletagSR_T5HH2100,h_J2M_doubletagSB_T5HH2100,h_J2M_tagSR_T5HH2100,h_J2M_tagSB_T5HH2100, h_J2M_antitagSR_T5HH2100, h_J2M_antitagSB_T5HH2100};

  vector<TH1F*> histos_j1pt_ZJets = {h_J1Pt_doubletagSR_ZJets,h_J1Pt_doubletagSB_ZJets,h_J1Pt_tagSR_ZJets,h_J1Pt_tagSB_ZJets, h_J1Pt_antitagSR_ZJets, h_J1Pt_antitagSB_ZJets};
  vector<TH1F*> histos_j1pt_WJets = {h_J1Pt_doubletagSR_WJets,h_J1Pt_doubletagSB_WJets,h_J1Pt_tagSR_WJets,h_J1Pt_tagSB_WJets, h_J1Pt_antitagSR_WJets, h_J1Pt_antitagSB_WJets};
  vector<TH1F*> histos_j1pt_TT    = {h_J1Pt_doubletagSR_TT,h_J1Pt_doubletagSB_TT,h_J1Pt_tagSR_TT,h_J1Pt_tagSB_TT, h_J1Pt_antitagSR_TT, h_J1Pt_antitagSB_TT};
  vector<TH1F*> histos_j1pt_QCD   = {h_J1Pt_doubletagSR_QCD,h_J1Pt_doubletagSB_QCD,h_J1Pt_tagSR_QCD,h_J1Pt_tagSB_QCD, h_J1Pt_antitagSR_QCD, h_J1Pt_antitagSB_QCD};
  vector<TH1F*> histos_j1pt_T5HH1300 = {h_J1Pt_doubletagSR_T5HH1300,h_J1Pt_doubletagSB_T5HH1300,h_J1Pt_tagSR_T5HH1300,h_J1Pt_tagSB_T5HH1300, h_J1Pt_antitagSR_T5HH1300, h_J1Pt_antitagSB_T5HH1300};
  vector<TH1F*> histos_j1pt_T5HH1700 = {h_J1Pt_doubletagSR_T5HH1700,h_J1Pt_doubletagSB_T5HH1700,h_J1Pt_tagSR_T5HH1700,h_J1Pt_tagSB_T5HH1700, h_J1Pt_antitagSR_T5HH1700, h_J1Pt_antitagSB_T5HH1700};
  vector<TH1F*> histos_j1pt_T5HH2100 = {h_J1Pt_doubletagSR_T5HH2100,h_J1Pt_doubletagSB_T5HH2100,h_J1Pt_tagSR_T5HH2100,h_J1Pt_tagSB_T5HH2100, h_J1Pt_antitagSR_T5HH2100, h_J1Pt_antitagSB_T5HH2100};

  vector<TH1F*> histos_j2pt_ZJets = {h_J2Pt_doubletagSR_ZJets,h_J2Pt_doubletagSB_ZJets,h_J2Pt_tagSR_ZJets,h_J2Pt_tagSB_ZJets, h_J2Pt_antitagSR_ZJets, h_J2Pt_antitagSB_ZJets};
  vector<TH1F*> histos_j2pt_WJets = {h_J2Pt_doubletagSR_WJets,h_J2Pt_doubletagSB_WJets,h_J2Pt_tagSR_WJets,h_J2Pt_tagSB_WJets, h_J2Pt_antitagSR_WJets, h_J2Pt_antitagSB_WJets};
  vector<TH1F*> histos_j2pt_TT    = {h_J2Pt_doubletagSR_TT,h_J2Pt_doubletagSB_TT,h_J2Pt_tagSR_TT,h_J2Pt_tagSB_TT, h_J2Pt_antitagSR_TT, h_J2Pt_antitagSB_TT};
  vector<TH1F*> histos_j2pt_QCD   = {h_J2Pt_doubletagSR_QCD,h_J2Pt_doubletagSB_QCD,h_J2Pt_tagSR_QCD,h_J2Pt_tagSB_QCD, h_J2Pt_antitagSR_QCD, h_J2Pt_antitagSB_QCD};
  vector<TH1F*> histos_j2pt_T5HH1300 = {h_J2Pt_doubletagSR_T5HH1300,h_J2Pt_doubletagSB_T5HH1300,h_J2Pt_tagSR_T5HH1300,h_J2Pt_tagSB_T5HH1300, h_J2Pt_antitagSR_T5HH1300, h_J2Pt_antitagSB_T5HH1300};
  vector<TH1F*> histos_j2pt_T5HH1700 = {h_J2Pt_doubletagSR_T5HH1700,h_J2Pt_doubletagSB_T5HH1700,h_J2Pt_tagSR_T5HH1700,h_J2Pt_tagSB_T5HH1700, h_J2Pt_antitagSR_T5HH1700, h_J2Pt_antitagSB_T5HH1700};
  vector<TH1F*> histos_j2pt_T5HH2100 = {h_J2Pt_doubletagSR_T5HH2100,h_J2Pt_doubletagSB_T5HH2100,h_J2Pt_tagSR_T5HH2100,h_J2Pt_tagSB_T5HH2100, h_J2Pt_antitagSR_T5HH2100, h_J2Pt_antitagSB_T5HH2100};

  //Combine lowMET and highMET into one TH1F
  TH1F* h_J1_DoubleBAllMET_doubletagSR_QCD   = (TH1F*)h_J1_DoubleBLowMET_doubletagSR_QCD->Clone("J1_DoubleBAllMET_doubletagSR_QCD"); h_J1_DoubleBAllMET_doubletagSR_QCD->Add(h_J1_DoubleBHighMET_doubletagSR_QCD);
  TH1F* h_J1_DoubleBAllMET_doubletagSR_TT    = (TH1F*)h_J1_DoubleBLowMET_doubletagSR_TT->Clone("J1_DoubleBAllMET_doubletagSR_TT");   h_J1_DoubleBAllMET_doubletagSR_TT->Add(h_J1_DoubleBHighMET_doubletagSR_TT);
  TH1F* h_J1_DoubleBAllMET_doubletagSR_WJets = (TH1F*)h_J1_DoubleBLowMET_doubletagSR_WJets->Clone("J1_DoubleBAllMET_doubletagSR_WJets"); h_J1_DoubleBAllMET_doubletagSR_WJets->Add(h_J1_DoubleBHighMET_doubletagSR_WJets);
  TH1F* h_J1_DoubleBAllMET_doubletagSR_ZJets = (TH1F*)h_J1_DoubleBLowMET_doubletagSR_ZJets->Clone("J1_DoubleBAllMET_doubletagSR_ZJets"); h_J1_DoubleBAllMET_doubletagSR_ZJets->Add(h_J1_DoubleBHighMET_doubletagSR_ZJets);
  TH1F* h_J1_DoubleBAllMET_doubletagSR_T5HH1300 = (TH1F*)h_J1_DoubleBLowMET_doubletagSR_T5HH1300->Clone("J1_DoubleBAllMET_doubletagSR_T5HH1300"); h_J1_DoubleBAllMET_doubletagSR_T5HH1300->Add(h_J1_DoubleBHighMET_doubletagSR_T5HH1300);
  TH1F* h_J1_DoubleBAllMET_doubletagSR_T5HH1700 = (TH1F*)h_J1_DoubleBLowMET_doubletagSR_T5HH1700->Clone("J1_DoubleBAllMET_doubletagSR_T5HH1700"); h_J1_DoubleBAllMET_doubletagSR_T5HH1700->Add(h_J1_DoubleBHighMET_doubletagSR_T5HH1700);
  TH1F* h_J1_DoubleBAllMET_doubletagSR_T5HH2100 = (TH1F*)h_J1_DoubleBLowMET_doubletagSR_T5HH2100->Clone("J1_DoubleBAllMET_doubletagSR_T5HH2100"); h_J1_DoubleBAllMET_doubletagSR_T5HH2100->Add(h_J1_DoubleBHighMET_doubletagSR_T5HH2100);

  TH1F* h_J1_DoubleBAllMET_tagSR_QCD   = (TH1F*)h_J1_DoubleBLowMET_tagSR_QCD->Clone("J1_DoubleBAllMET_tagSR_QCD"); h_J1_DoubleBAllMET_tagSR_QCD->Add(h_J1_DoubleBHighMET_tagSR_QCD);
  TH1F* h_J1_DoubleBAllMET_tagSR_TT    = (TH1F*)h_J1_DoubleBLowMET_tagSR_TT->Clone("J1_DoubleBAllMET_tagSR_TT"); h_J1_DoubleBAllMET_tagSR_TT->Add(h_J1_DoubleBHighMET_tagSR_TT);
  TH1F* h_J1_DoubleBAllMET_tagSR_WJets = (TH1F*)h_J1_DoubleBLowMET_tagSR_WJets->Clone("J1_DoubleBAllMET_tagSR_WJets"); h_J1_DoubleBAllMET_tagSR_WJets->Add(h_J1_DoubleBHighMET_tagSR_WJets);
  TH1F* h_J1_DoubleBAllMET_tagSR_ZJets = (TH1F*)h_J1_DoubleBLowMET_tagSR_ZJets->Clone("J1_DoubleBAllMET_tagSR_ZJets"); h_J1_DoubleBAllMET_tagSR_ZJets->Add(h_J1_DoubleBHighMET_tagSR_ZJets);
  TH1F* h_J1_DoubleBAllMET_tagSR_T5HH1300 = (TH1F*)h_J1_DoubleBLowMET_tagSR_T5HH1300->Clone("J1_DoubleBAllMET_tagSR_T5HH1300"); h_J1_DoubleBAllMET_tagSR_T5HH1300->Add(h_J1_DoubleBHighMET_tagSR_T5HH1300);
  TH1F* h_J1_DoubleBAllMET_tagSR_T5HH1700 = (TH1F*)h_J1_DoubleBLowMET_tagSR_T5HH1700->Clone("J1_DoubleBAllMET_tagSR_T5HH1700"); h_J1_DoubleBAllMET_tagSR_T5HH1700->Add(h_J1_DoubleBHighMET_tagSR_T5HH1700);
  TH1F* h_J1_DoubleBAllMET_tagSR_T5HH2100 = (TH1F*)h_J1_DoubleBLowMET_tagSR_T5HH2100->Clone("J1_DoubleBAllMET_tagSR_T5HH2100"); h_J1_DoubleBAllMET_tagSR_T5HH2100->Add(h_J1_DoubleBHighMET_tagSR_T5HH2100);


  if (runABCDPlots){
    if (whichRegion=="signal"){
      // makeABCDPlot(histos_ABCD_sum, "BkgSum", "Double");
      makeABCDPlot(histos_ABCD_QCD, "QCD", "Double");
      makeABCDPlot(histos_ABCD_TT, "TT", "Double");
      makeABCDPlot(histos_ABCD_WJets, "WJets", "Double");
      makeABCDPlot(histos_ABCD_ZJets, "ZJets", "Double");
      //
      // makeABCDPlot(histos_A1B1CD_sum, "BkgSum", "Single");
      makeABCDPlot(histos_A1B1CD_QCD, "QCD", "Single");
      makeABCDPlot(histos_A1B1CD_TT, "TT", "Single");
      makeABCDPlot(histos_A1B1CD_WJets, "WJets", "Single");
      makeABCDPlot(histos_A1B1CD_ZJets, "ZJets", "Single");
    }

    else if (whichRegion=="singleLept"){
      // makeABCDPlot(histos_ABCD_sum, "BkgSum", "Double");
      makeABCDPlot(histos_ABCD_TT, "TT", "Double");
      makeABCDPlot(histos_ABCD_WJets, "WJets", "Double");

      makeABCDPlot(histos_ABCD_sum, "BkgSum", "Single");
      makeABCDPlot(histos_A1B1CD_TT, "TT", "Single");
      makeABCDPlot(histos_A1B1CD_WJets, "WJets", "Single");
    }
    else if (whichRegion=="photon"){
      makeABCDPlot(histos_ABCD_QCD, "QCD", "Double");
      makeABCDPlot(histos_ABCD_GJets, "GJets", "Double");
      makeABCDPlot(histos_A1B1CD_QCD, "QCD", "Single");
      makeABCDPlot(histos_A1B1CD_GJets, "GJets", "Single");
    }


  }

  std::cout<<"Finished with ABCD plots\n";

  if (runRPFPlots){
    myfile.open("Rpf_deviations.txt");

    if (whichRegion=="signal"){
      //J1 Double tag region
      // makeRpfPlot(histos_Rpf_J1_sum, "BkgSum", "J1", "Double");
      makeRpfPlot(histos_Rpf_J1_QCD, "QCD", "J1", "Double");
      makeRpfPlot(histos_Rpf_J1_TT, "TT", "J1", "Double");
      makeRpfPlot(histos_Rpf_J1_WJets, "WJets", "J1", "Double");
      makeRpfPlot(histos_Rpf_J1_ZJets, "ZJets", "J1", "Double");

      //J1 Single tag region
      // makeRpfPlot(histos_Rpfsingle_J1_sum, "BkgSum", "J1", "Single");
      makeRpfPlot(histos_Rpfsingle_J1_QCD, "QCD", "J1", "Single");
      makeRpfPlot(histos_Rpfsingle_J1_TT, "TT", "J1", "Single");
      makeRpfPlot(histos_Rpfsingle_J1_WJets, "WJets", "J1", "Single");
      makeRpfPlot(histos_Rpfsingle_J1_ZJets, "ZJets", "J1", "Single");

      //J2 double tag region
      // makeRpfPlot(histos_Rpf_J2_sum, "BkgSum", "J2", "Double");
      makeRpfPlot(histos_Rpf_J2_QCD, "QCD", "J2", "Double");
      makeRpfPlot(histos_Rpf_J2_TT, "TT", "J2", "Double");
      makeRpfPlot(histos_Rpf_J2_WJets, "WJets", "J2", "Double");
      makeRpfPlot(histos_Rpf_J2_ZJets, "ZJets", "J2", "Double");

      //J2 single tag region
      // makeRpfPlot(histos_Rpfsingle_J2_sum, "BkgSum", "J2", "Single");
      makeRpfPlot(histos_Rpfsingle_J2_QCD, "QCD", "J2", "Single");
      makeRpfPlot(histos_Rpfsingle_J2_TT, "TT", "J2", "Single");
      makeRpfPlot(histos_Rpfsingle_J2_WJets, "WJets", "J2", "Single");
      makeRpfPlot(histos_Rpfsingle_J2_ZJets, "ZJets", "J2", "Single");
    }
    if (whichRegion=="singleLept"){
      //J1 Double tag region
      makeRpfPlot(histos_Rpf_J1_TT, "TT", "J1", "Double");
      makeRpfPlot(histos_Rpf_J1_WJets, "WJets", "J1", "Double");

      //J1 Single tag region
      makeRpfPlot(histos_Rpfsingle_J1_TT, "TT", "J1", "Single");
      makeRpfPlot(histos_Rpfsingle_J1_WJets, "WJets", "J1", "Single");

      //J2 double tag region
      makeRpfPlot(histos_Rpf_J2_TT, "TT", "J2", "Double");
      makeRpfPlot(histos_Rpf_J2_WJets, "WJets", "J2", "Double");

      //J2 single tag region
      makeRpfPlot(histos_Rpfsingle_J2_TT, "TT", "J2", "Single");
      makeRpfPlot(histos_Rpfsingle_J2_WJets, "WJets", "J2", "Single");
    }
    if (whichRegion=="photon"){
      //J1 Double tag region
      makeRpfPlot(histos_Rpf_J1_QCD, "QCD", "J1", "Double");
      makeRpfPlot(histos_Rpf_J1_GJets, "GJets", "J1", "Double");

      //J1 Single tag region
      makeRpfPlot(histos_Rpfsingle_J1_QCD, "QCD", "J1", "Single");
      makeRpfPlot(histos_Rpfsingle_J1_GJets, "GJets", "J1", "Single");

      //J2 double tag region
      makeRpfPlot(histos_Rpf_J2_QCD, "QCD", "J2", "Double");
      makeRpfPlot(histos_Rpf_J2_GJets, "GJets", "J2", "Double");

      //J2 single tag region
      makeRpfPlot(histos_Rpfsingle_J2_QCD, "QCD", "J2", "Single");
      makeRpfPlot(histos_Rpfsingle_J2_GJets, "GJets", "J2", "Single");
    }

    myfile.close();
  }

  if (runDoubleBStack && whichRegion=="signal") {
    // makeStackPlot(histos_doubleBJ1_QCD,histos_doubleBJ1_TT,histos_doubleBJ1_WJets,histos_doubleBJ1_ZJets,histos_doubleBJ1_T5HH1300,histos_doubleBJ1_T5HH1700,histos_doubleBJ1_T5HH2100,"DoubleB");
    makeDoubleBStack(h_J1_DoubleBAllMET_doubletagSR_QCD,h_J1_DoubleBAllMET_doubletagSR_TT,h_J1_DoubleBAllMET_doubletagSR_WJets,h_J1_DoubleBAllMET_doubletagSR_ZJets, h_J1_DoubleBAllMET_doubletagSR_T5HH1300,h_J1_DoubleBAllMET_doubletagSR_T5HH1700,h_J1_DoubleBAllMET_doubletagSR_T5HH2100,"doubleAll");
    makeDoubleBStack(h_J1_DoubleBLowMET_doubletagSR_QCD,h_J1_DoubleBLowMET_doubletagSR_TT,h_J1_DoubleBLowMET_doubletagSR_WJets,h_J1_DoubleBLowMET_doubletagSR_ZJets, h_J1_DoubleBLowMET_doubletagSR_T5HH1300,h_J1_DoubleBLowMET_doubletagSR_T5HH1700,h_J1_DoubleBLowMET_doubletagSR_T5HH2100,"doubleLow");
    makeDoubleBStack(h_J1_DoubleBHighMET_doubletagSR_QCD,h_J1_DoubleBHighMET_doubletagSR_TT,h_J1_DoubleBHighMET_doubletagSR_WJets,h_J1_DoubleBHighMET_doubletagSR_ZJets, h_J1_DoubleBHighMET_doubletagSR_T5HH1300,h_J1_DoubleBHighMET_doubletagSR_T5HH1700,h_J1_DoubleBHighMET_doubletagSR_T5HH2100,"doubleHigh");

    makeDoubleBStack(h_J1_DoubleBAllMET_tagSR_QCD,h_J1_DoubleBAllMET_tagSR_TT,h_J1_DoubleBAllMET_tagSR_WJets,h_J1_DoubleBAllMET_tagSR_ZJets, h_J1_DoubleBAllMET_tagSR_T5HH1300,h_J1_DoubleBAllMET_tagSR_T5HH1700,h_J1_DoubleBAllMET_tagSR_T5HH2100,"singleAll");
    makeDoubleBStack(h_J1_DoubleBLowMET_tagSR_QCD,h_J1_DoubleBLowMET_tagSR_TT,h_J1_DoubleBLowMET_tagSR_WJets,h_J1_DoubleBLowMET_tagSR_ZJets, h_J1_DoubleBLowMET_tagSR_T5HH1300,h_J1_DoubleBLowMET_tagSR_T5HH1700,h_J1_DoubleBLowMET_tagSR_T5HH2100,"singleLow");
    makeDoubleBStack(h_J1_DoubleBHighMET_tagSR_QCD,h_J1_DoubleBHighMET_tagSR_TT,h_J1_DoubleBHighMET_tagSR_WJets,h_J1_DoubleBHighMET_tagSR_ZJets, h_J1_DoubleBHighMET_tagSR_T5HH1300,h_J1_DoubleBHighMET_tagSR_T5HH1700,h_J1_DoubleBHighMET_tagSR_T5HH2100,"singleHigh");

  }

  if (runStacks && whichRegion=="signal") {
    makeStackPlot(histos_j1mass_QCD,histos_j1mass_TT,histos_j1mass_WJets,histos_j1mass_ZJets,histos_j1mass_T5HH1300,histos_j1mass_T5HH1700,histos_j1mass_T5HH2100,"leadmass");
    makeStackPlot(histos_j2mass_QCD,histos_j2mass_TT,histos_j2mass_WJets,histos_j2mass_ZJets,histos_j2mass_T5HH1300,histos_j2mass_T5HH1700,histos_j2mass_T5HH2100,"subleadmass");

    makeStackPlot(histos_j1pt_QCD,histos_j1pt_TT,histos_j1pt_WJets,histos_j1pt_ZJets,histos_j1pt_T5HH1300,histos_j1pt_T5HH1700,histos_j1pt_T5HH2100,"leadpt");
    makeStackPlot(histos_j2pt_QCD,histos_j2pt_TT,histos_j2pt_WJets,histos_j2pt_ZJets,histos_j2pt_T5HH1300,histos_j2pt_T5HH1700,histos_j2pt_T5HH2100,"subleadpt");

    // makeStackPlot(histos_doubleBJ1_QCD,histos_doubleBJ1_TT,histos_doubleBJ1_WJets,histos_doubleBJ1_ZJets,histos_doubleBJ1_T5HH1300,histos_doubleBJ1_T5HH1700,histos_doubleBJ1_T5HH2100,"DoubleB");
    // makeStackPlot(histos_doubleBJ2_QCD,histos_doubleBJ2_TT,histos_doubleBJ2_WJets,histos_doubleBJ2_ZJets,histos_doubleBJ2_T5HH1300,histos_doubleBJ2_T5HH1700,histos_doubleBJ2_T5HH2100,"DoubleB");
  }


  if (runStacks && whichRegion=="singleLept") {
    makeSingleLeptStackPlot(histos_j1mass_TT,histos_j1mass_WJets,"leadmass");
    makeSingleLeptStackPlot(histos_j2mass_TT,histos_j2mass_WJets,"subleadmass");

    makeSingleLeptStackPlot(histos_j1pt_TT,histos_j1pt_WJets,"leadpt");
    makeSingleLeptStackPlot(histos_j2pt_TT,histos_j2pt_WJets,"subleadpt");

    // makeSingleLeptStackPlot(histos_doubleBJ1_TT,histos_doubleBJ1_WJets,"DoubleB");
    // makeSingleLeptStackPlot(histos_doubleBHJ1_TT,histos_doubleBHJ1_WJets,"DoubleB-H");
  }

  if (runMETStacks && whichRegion=="signal") {
    // makeMETStack(h_MET_all_QCD,h_MET_all_TT,h_MET_all_WJets,h_MET_all_ZJets,h_MET_all_T5HH1300,h_MET_all_T5HH1700,h_MET_all_T5HH2100,"all");
    // makeMETStack(h_MET_Single_QCD,h_MET_Single_TT,h_MET_Single_WJets,h_MET_Single_ZJets,h_MET_Single_T5HH1300,h_MET_Single_T5HH1700,h_MET_Single_T5HH2100,"single");
    // makeMETStack(h_MET_Double_QCD,h_MET_Double_TT,h_MET_Double_WJets,h_MET_Double_ZJets,h_MET_Double_T5HH1300,h_MET_Double_T5HH1700,h_MET_Double_T5HH2100,"double");

    makeMETStack(h_A_QCD, h_A_TT, h_A_WJets, h_A_ZJets, h_A_T5HH1300, h_A_T5HH1700, h_A_T5HH2100, "double");
    makeMETStack(h_A1_QCD, h_A1_TT, h_A1_WJets, h_A1_ZJets, h_A1_T5HH1300, h_A1_T5HH1700, h_A1_T5HH2100, "single");
}


  // if (runCompareDoubleB){
  //   if (whichRegion=="signal"){
  //     makeDoubleBPlots(histos_doubleBJ1_ZJets, "ZJets");
  //     makeDoubleBPlots(histos_doubleBJ1_WJets, "WJets");
  //     makeDoubleBPlots(histos_doubleBJ1_TT, "TT");
  //     makeDoubleBPlots(histos_doubleBJ1_QCD, "QCD");
  //     makeDoubleBPlots(histos_doubleBJ1_T5HH1300, "T5HH1300");
  //     makeDoubleBPlots(histos_doubleBJ1_T5HH1700, "T5HH1700");
  //     makeDoubleBPlots(histos_doubleBJ1_T5HH2100, "T5HH2100");
  //
  //   }
  //   if (whichRegion=="singleLept"){
  //     makeDoubleBPlots(histos_doubleBJ1_WJets, "WJets");
  //     makeDoubleBPlots(histos_doubleBJ1_TT, "TT");
  //   }
  //
  // }

} //end method ABCD()

void makeABCDPlot(vector<TH1F*> dem_histos, TString bkgType = "", TString tagType = "") {
  TH1F *histo_A = dem_histos[0]; TH1F *histo_B = dem_histos[1];
  TH1F *histo_C = dem_histos[2]; TH1F *histo_D = dem_histos[3];
  histo_B->Multiply(histo_C);
  histo_B->Divide(histo_D);

  TH1F *histo_pred = histo_B;
  // TH1F histo_ratio = (*histo_pred)/(*histo_A);

  for (int i = 1; i<=3; i++){
    if (histo_A->GetBinContent(i)==0){
      std::cout<<"Attempting to solve problem"<<std::endl;
      histo_A->SetBinContent(i,0.001);
      if (histo_A->GetBinError(i)==0) histo_A->SetBinError(i,1.4);
    }
    if (histo_pred->GetBinContent(i)==0){
      std::cout<<"Attempting to solve problem"<<std::endl;
      histo_pred->SetBinContent(i,0.001);
      if (histo_pred->GetBinError(i)==0) histo_pred->SetBinError(i,1.4);
    }
  }
  // histo_A->Sumw2();
  // histo_pred->Sumw2();

  TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo_A, histo_pred, "pois");

  TString canvName = bkgType+tagType;

  double W = 800;
  double H = 600;
  double T = 0.08*H;
  double B = 0.12*H;
  double L = 0.08*W;
  double R = 0.08*W;
  TCanvas * can_h = new TCanvas(canvName,canvName, 50, 50, W, H);
  can_h->SetFillColor(0);
  can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0);
  can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W );
  can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H );
  can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0);
  can_h->SetTicky(0);

  double up_height     = 0.8;
  // double dw_correction = 1.30;
  double dw_correction = 1.18;
  double font_size_dw  = 0.1;
  double dw_height = (1. - up_height) * dw_correction;
  double dw_height_offset = 0.02;

  TPad *pad1 = new TPad("pad1", "top pad" , 0.0, 0.25, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.25);
  pad1->SetTickx(0);
  pad1->SetTicky(0);
  pad1->SetPad(0., 1 - up_height, 1., 1.00);
  pad1->SetFrameFillColor(0);
  pad1->SetFillColor(0);
  pad1->SetTopMargin(0.08);
  pad1->SetLeftMargin(0.08);
  pad1->SetRightMargin(0.06);
  // pad1->SetLogy();
  pad1->Draw();

  pad2->SetPad(0., 0., 1., dw_height+dw_height_offset);
  pad2->SetFillColor(0);
  pad2->SetFrameFillColor(0);
  pad2->SetBottomMargin(0.33);
  pad2->SetTopMargin(0.01);
  pad2->SetLeftMargin(0.08);
  pad2->SetRightMargin(0.06);
  pad2->Draw();
  pad1->cd();

  histo_pred->SetStats(0); histo_A->SetStats(0);
  histo_pred->SetFillColor(kBlue);
  histo_pred->SetFillStyle(3445);

  histo_pred->SetMarkerSize(0);
  histo_pred->SetLineWidth(1);
  histo_pred->SetLineColor(kBlue);

  histo_A->SetTitle(";;nEvents");
  TString title = bkgType+";;nEvents";
  histo_pred->SetTitle(title);
  histo_pred->SetMinimum(0);
  if (tagType == "Double" && bkgType == "TT") histo_pred->SetMaximum(30.0);
  //
  // if (bkgType == "WJets" && tagType == "Double") histo_pred->SetMaximum(0.5);
  // histo_pred->SetTitle(";;nEvents");
  histo_pred->GetYaxis()->SetTitleOffset(0.6);
  //histo_pred->GetYaxis()->SetTitleSize(0.14);
  histo_pred->GetXaxis()->SetLabelSize(0);
  // histo_pred->GetYaxis()->SetNdivisions(505);
  histo_A->SetLineColor(kBlack); histo_A->SetMarkerStyle(20);
  histo_A->SetMarkerColor(kBlack);

  TLegend* legend = new TLegend(0.60,0.7,0.94,0.92);
  //legend->SetHeader(bkgType,"C"); // option "C" allows to center the header
  if (tagType=="Double") {
    legend->AddEntry(histo_A,"Simulation: A2","lp");
    legend->AddEntry(histo_pred,"Prediction: B2*C/D","f");

  }
  else if (tagType=="Single") {
    legend->AddEntry(histo_A,"Simulation: A1","lp");
    legend->AddEntry(histo_pred,"Prediction: B1*C/D","f");
  }

  histo_pred->Draw("E2");
  histo_A->Draw("same");
  legend->Draw("same");

  pad2->cd();
  graph->SetTitle(";MET [GeV];sim/pred");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlack);
  graph->GetYaxis()->SetTitleOffset(0.25);
  graph->GetYaxis()->SetTitleSize(0.13);
  graph->GetXaxis()->SetTitleOffset(0.9);
  graph->GetXaxis()->SetTitleSize(0.14);
  graph->GetXaxis()->SetLabelSize(0.14);
  graph->GetYaxis()->SetLabelSize(0.14);
  graph->GetXaxis()->SetNdivisions(507);
  graph->GetYaxis()->SetNdivisions(505);
  // graph->GetXaxis()->SetRangeUser(300,1000);
  // graph->GetYaxis()->SetRangeUser(0.3,1.5);

  graph->Draw("APE");
  graph->GetXaxis()->SetRangeUser(150,1000);
  if (tagType == "Double" &&  !useDeepDoubleB) graph->GetYaxis()->SetRangeUser(0.0,4.0);
  else if (tagType == "Double") graph->GetYaxis()->SetRangeUser(0.0,2.0);
  else if (tagType == "Single" &&  bkgType == "ZJets" && useDeepDoubleB) graph->GetYaxis()->SetRangeUser(0.8,1.6);
  else if (tagType == "Single" &&  bkgType == "TT" && !useDeepDoubleB) graph->GetYaxis()->SetRangeUser(0.0,1.2);
  else graph->GetYaxis()->SetRangeUser(0.8,1.2);

  // if (bkgType == "TT") graph->GetYaxis()->SetRangeUser(0.2,4.5);
  // graph->Draw("APE");



  can_h->Modified();
  can_h->Update();



  TLine *line = new TLine(150,1.0,1000,1.0);
  line->SetLineColor(kRed);
  line->SetLineStyle(2);
  line->Draw("same");

  TString thisDirectory = "ABCDPlots/ABCD_";
  if (whichRegion == "singleLept") thisDirectory = "ABCDPlots/ABCD_SingleLept_";
  if (whichRegion == "photon") thisDirectory = "ABCDPlots/ABCD_Photon_";

  TString savename = directory+thisDirectory+tagType+"_"+bkgType;
  if (useDeepDoubleB) savename = directory+thisDirectory+tagType+"_"+bkgType+"_deepDB";
  can_h->SaveAs(savename+".root");
  can_h->SaveAs(savename+".pdf");
}

void makeRpfPlot(vector<TH1F*> dem_histos, TString bkgType = "", TString jetType = "", TString tagType = "") {
  TH1D * num   = (TH1D*)dem_histos[0]->Clone("num"+tagType+"_"+jetType+"_"+bkgType);
  TH1D * denom = (TH1D*)dem_histos[2]->Clone("denom"+tagType+"_"+jetType+"_"+bkgType);

  num->Add(dem_histos[1]); denom->Add(dem_histos[3]);
  num->Rebin(2); denom->Rebin(2);

  num->SetMarkerStyle(20); num->SetMarkerSize(0.85);
  num->SetMinimum(0.0);
  num->SetLineColor(kRed); num->SetMarkerColor(kRed);
  num->GetYaxis()->SetTitle("Normalized to unity");


  denom->SetMarkerStyle(20); denom->SetMarkerSize(0.85);
  denom->SetMinimum(0.0);
  denom->SetLineColor(kBlue); denom->SetMarkerColor(kBlue);

  TString thisType = tagType+"_"+jetType+"_"+bkgType;
  TString numTitle = "num"+thisType; TString denomTitle = "denom"+thisType;
  num->SetTitle(numTitle); denom->SetTitle(denomTitle);
  TString graphName = "ratio"+thisType;

  // num->Divide(denom);
  // num->SetName(graphName);
  // num->GetXaxis()->SetTitle("Soft-drop mass [GeV]");
  // num->GetYaxis()->SetTitle("R_{p/f}");
  // num->SetTitle(graphName);
  // num->GetYaxis()->SetTitleSize(0.055);
  // num->SetMarkerStyle(20);
  // num->SetMarkerSize(0.85);
  // num->SetMinimum(0.0);
  // num->SetLineColor(kRed);
  // num->SetMarkerColor(kRed);

  TGraphAsymmErrors * graph = new TGraphAsymmErrors(num, denom, "pois");
  graph->SetName(graphName);
  graph->GetXaxis()->SetTitle("Soft-drop mass [GeV]");
  graph->GetYaxis()->SetTitle("R_{p/f}");
  graph->SetTitle("");
  graph->GetYaxis()->SetTitleSize(0.055);
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.85);
  graph->SetMinimum(0.0);
  graph->SetLineColor(kBlack);
  graph->SetMarkerColor(kBlack);
  graph->GetXaxis()->SetRangeUser(50,250);

  double W = 800;
  double H = 600;
  double T = 0.08*H;
  double B = 0.12*H;
  double L = 0.15*W;
  double R = 0.02*W;
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0);
  can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0);
  can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W );
  can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H );
  can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0);
  can_h->SetTicky(0);

  double up_height     = 0.8;
  double dw_correction = 1.18;
  double font_size_dw  = 0.1;
  double dw_height = (1. - up_height) * dw_correction;
  double dw_height_offset = 0.02;

  TPad *pad1 = new TPad("pad1", "top pad" , 0.0, 0.25, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.25);
  pad1->SetTickx(0);
  pad1->SetTicky(0);
  pad1->SetPad(0., 1 - up_height, 1., 1.00);
  pad1->SetFrameFillColor(0);
  pad1->SetFillColor(0);
  pad1->SetTopMargin(0.08);
  pad1->SetLeftMargin(0.10);
  pad1->SetRightMargin(0.04);
  pad1->Draw();

  pad2->SetPad(0., 0., 1., dw_height+dw_height_offset);
  pad2->SetFillColor(0);
  pad2->SetFrameFillColor(0);
  pad2->SetBottomMargin(0.33);
  pad2->SetTopMargin(0.01);
  pad2->SetLeftMargin(0.10);
  pad2->SetRightMargin(0.04);
  pad2->Draw();
  pad1->cd();

  num->GetXaxis()->SetLabelSize(0);
  num->SetTitle(thisType);

  num->DrawNormalized();
  denom->DrawNormalized("same");

  TLegend* legend = new TLegend(0.5426065,0.7935484,0.9498246,0.9403226,NULL,"brNDC") ;
  if (tagType == "Single"){
    legend->AddEntry(num, "SingleTag (numerator)", "PL");
    legend->AddEntry(denom, "AntiTag (denominator)", "PL");
  }
  else if (tagType == "Double"){
    legend->AddEntry(num, "DoubleTag (numerator)", "PL");
    legend->AddEntry(denom, "AntiTag (denominator)", "PL");
  }
  legend->Draw();

  pad2->cd();
  graph->GetYaxis()->SetTitleOffset(0.25);
  graph->GetYaxis()->SetTitleSize(0.13);
  graph->GetXaxis()->SetTitleOffset(0.9);
  graph->GetXaxis()->SetTitleSize(0.14);
  graph->GetXaxis()->SetLabelSize(0.14);
  graph->GetYaxis()->SetLabelSize(0.10);
  // graph->GetXaxis()->SetNdivisions(507);
  graph->GetYaxis()->SetNdivisions(505);
  graph->Draw("APE");

  if (tagType=="Double"){
    if (bkgType=="QCD") graph->GetYaxis()->SetRangeUser(0,0.5);
    else if (bkgType=="TT") graph->GetYaxis()->SetRangeUser(0,1.2);
    else if (bkgType=="WJets" || bkgType=="ZJets") graph->GetYaxis()->SetRangeUser(0,0.03);
  }

  else if (tagType=="Single"){
    if (bkgType=="QCD") graph->GetYaxis()->SetRangeUser(0,1.0);
    else if (bkgType=="TT") graph->GetYaxis()->SetRangeUser(0,3.0);
    else if (bkgType=="WJets" || bkgType=="ZJets") graph->GetYaxis()->SetRangeUser(0,0.3);
  }
  graph->Draw("APE");
  can_h->Modified();
  can_h->Update();

  // TLine *line = new TLine(300,1.0,1000,1.0);
  // line->SetLineColor(kRed);
  // line->SetLineStyle(2);
  // line->Draw("same");

  if (getDeviation){
    TF1*f0=new TF1("f0","pol0", 50,250);
    graph->Fit("f0","Q","R",50,180);
    //TF1 * thisFit = graph->GetFunction("pol0");
    double p0 = f0->GetParameter(0); //this does, indeed, work
    graph->GetFunction("f0")->SetBit(TF1::kNotDraw);

    TLine *line = new TLine(50.0,p0,250.0,p0);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
    //line->SetLineStyle(2);
    line->Draw("same");


    int nbinsx = graph->GetN();
    double ax[nbinsx],ay[nbinsx];
    vector<double> deviation ={0};
    for (int i=0; i< nbinsx; i++){
      graph->GetPoint(i,ax[i],ay[i]);
      const double x0 = ax[i]; const double y0 = ay[i];
      const double errorY = graph->GetErrorYlow(i);
      double dev = (y0-p0)/errorY;
      deviation.push_back(dev);
    }
    myfile<<thisType<<": ";
    for (int i=0;i< deviation.size();i++) myfile<<deviation[i]<<", ";
    myfile<<"p0: "<<p0<<endl<<endl;

  }

  TString thisDirectory = "RPFPlots/";
  if (whichRegion == "singleLept") thisDirectory = "RPFPlots/SingleLept";
  if (whichRegion == "photon") thisDirectory = "RPFPlots/Photon";

  if (useDeepDoubleB){
    can_h->SaveAs(directory+thisDirectory+"Rpf_"+thisType+"_deepDB.pdf");
    num->SaveAs(directory+thisDirectory+"Num_"+thisType+"_deepDB.root");
    denom->SaveAs(directory+thisDirectory+"Denom_"+thisType+"_deepDB.root");
    graph->SaveAs(directory+thisDirectory+"Ratio_"+thisType+"_deepDB.root");
  }
  else {
    can_h->SaveAs(directory+thisDirectory+"Rpf_"+thisType+".pdf");
    num->SaveAs(directory+thisDirectory+"Num_"+thisType+".root");
    denom->SaveAs(directory+thisDirectory+"Denom_"+thisType+".root");
    graph->SaveAs(directory+thisDirectory+"Ratio_"+thisType+".root");
  }


}

void makeDoubleBPlots(vector<TH1F*> dem_histos, TString bkgType = "", TString which = ""){
  TH1F *h_SBLow = dem_histos[0]; TH1F *h_SR = dem_histos[1]; TH1F *h_SBHigh = dem_histos[2];

  h_SBLow->Rebin(2); h_SR->Rebin(2); h_SBHigh->Rebin(2);
  double W = 800;
  double H = 600;
  double T = 0.08*H;
  double B = 0.12*H;
  double L = 0.15*W;
  double R = 0.02*W;
  TString graphName = "";
  if (which =="") graphName = "DoubleB_"+bkgType;
  // else if (which =="DoubleB-H") graphName = "DoubleBH_"+bkgType;
  // else if (which =="tau32") graphName = "Tau32_"+bkgType;

  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0);
  can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0);
  can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W );
  can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H );
  can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0);
  can_h->SetTicky(0);
  can_h->SetLogy();

  h_SBLow->GetXaxis()->SetTitle("Double-b discriminator");
  if (useDeepDoubleB) h_SBLow->GetXaxis()->SetTitle("Deep double-b discriminator");
  // else if (which =="DoubleB-H") h_SBLow->GetXaxis()->SetTitle("Double-b discriminator, H");
  // else if (which =="tau32") {
  //   h_SBLow->GetXaxis()->SetTitle("Tau32");
  //   h_SBLow->GetXaxis()->SetRangeUser(0,1);
  // }

  h_SBLow->GetYaxis()->SetTitle("");
  h_SBLow->SetTitle(graphName);

  h_SBLow->SetMarkerStyle(20); //h_SBLow->SetMarkerSize(0.85);
  h_SBLow->SetLineColor(kBlack); h_SBLow->SetMarkerColor(kBlack);

  h_SR->SetMarkerStyle(20); //h_SR->SetMarkerSize(0.85);
  h_SR->SetLineColor(kRed); h_SR->SetMarkerColor(kRed);

  h_SBHigh->SetMarkerStyle(20); //h_SBHigh->SetMarkerSize(0.85);
  h_SBHigh->SetLineColor(kBlue); h_SBHigh->SetMarkerColor(kBlue);

  TLegend* legend = new TLegend(0.63,0.75,0.98,0.92) ;
  legend->AddEntry(h_SBLow, "AK8 Mass [50,85]", "PL");
  legend->AddEntry(h_SR, "AK8 Mass [85,135]", "PL");
  legend->AddEntry(h_SBHigh, "AK8 Mass [135,250]", "PL");
  h_SBLow->DrawNormalized();
  h_SR->DrawNormalized("same");
  h_SBHigh->DrawNormalized("same");
  legend->Draw("same");



  if (useDeepDoubleB) can_h->SaveAs(directory+"SupportingPlots/DoubleBCompare_"+bkgType+"_deepDB.pdf");
  else can_h->SaveAs(directory+"SupportingPlots/DoubleBCompare_"+bkgType+".pdf");

  // else if (which =="DoubleB-H") can_h->SaveAs("SupportingPlots/DoubleBHCompare_"+bkgType+".pdf");
  // else if (which =="tau32") can_h->SaveAs("SupportingPlots/Tau32Compare_"+bkgType+".pdf");

  // h_SBLow->ClearUnderflowAndOverflow();
  // double failIntegralSBLow = h_SBLow->Integral(0,16)/h_SBLow->Integral()*100;
  // double failIntegralSBHigh = h_SBHigh->Integral(0,16)/h_SBHigh->Integral()*100;
  // double failIntegralSR = h_SR->Integral(0,16)/h_SR->Integral()*100;
  //
  // cout<<bkgType<<": ";
  // cout<<"Percent failing, low mass: "<<failIntegralSBLow<<endl;
  // cout<<"Percent failing, SR: "<<failIntegralSR<<endl;
  // cout<<"Percent failing, high mass: "<<failIntegralSBHigh<<endl<<endl;

}

void makeStackPlot(vector<TH1F*> h_QCD,vector<TH1F*> h_TT,vector<TH1F*> h_WJets,vector<TH1F*> h_ZJets, vector<TH1F*> h_T5HH1300,vector<TH1F*> h_T5HH1700,vector<TH1F*> h_T5HH2100, TString which = ""){
  TH1D *h_QCD_stack = (TH1D*)h_QCD[0]->Clone("QCDAll");
  h_QCD_stack->Add(h_QCD[1]); h_QCD_stack->Add(h_QCD[2]);
  if (which!="DoubleB" && which!="DoubleB-H" && which!="tau32"){
    h_QCD_stack->Add(h_QCD[3]); h_QCD_stack->Add(h_QCD[4]); h_QCD_stack->Add(h_QCD[5]);
  }
  h_QCD_stack->SetFillColor(kGray);
  h_QCD_stack->SetMarkerStyle(21); h_QCD_stack->SetMarkerColor(kGray);

  TH1D *h_TT_stack = (TH1D*)h_TT[0]->Clone("TTAll");
  h_TT_stack->Add(h_TT[1]); h_TT_stack->Add(h_TT[2]);
  if (which!="DoubleB" && which!="DoubleB-H" && which!="tau32"){
    h_TT_stack->Add(h_TT[3]); h_TT_stack->Add(h_TT[4]); h_TT_stack->Add(h_TT[5]);
  }
  h_TT_stack->SetFillColor(kCyan);
  h_TT_stack->SetMarkerStyle(21); h_TT_stack->SetMarkerColor(kCyan);

  TH1D *h_WJets_stack = (TH1D*)h_WJets[0]->Clone("WJetsAll");
  h_WJets_stack->Add(h_WJets[1]); h_WJets_stack->Add(h_WJets[2]);
  if (which!="DoubleB" && which!="DoubleB-H" && which!="tau32"){
    h_WJets_stack->Add(h_WJets[3]); h_WJets_stack->Add(h_WJets[4]); h_WJets_stack->Add(h_WJets[5]);
  }
  h_WJets_stack->SetFillColor(kBlue);
  h_WJets_stack->SetMarkerStyle(21); h_WJets_stack->SetMarkerColor(kBlue);

  TH1D *h_ZJets_stack = (TH1D*)h_ZJets[0]->Clone("WZetsAll");
  h_ZJets_stack->Add(h_ZJets[1]); h_ZJets_stack->Add(h_ZJets[2]);
  if (which!="DoubleB" && which!="DoubleB-H" && which!="tau32"){
    h_ZJets_stack->Add(h_ZJets[3]); h_ZJets_stack->Add(h_ZJets[4]); h_ZJets_stack->Add(h_ZJets[5]);
  }
  h_ZJets_stack->SetFillColor(kGreen+2);
  h_ZJets_stack->SetMarkerStyle(21); h_ZJets_stack->SetMarkerColor(kGreen+2);

  TH1D *h_T5HH1300_stack = (TH1D*)h_T5HH1300[0]->Clone("T5HH1300");
  h_T5HH1300_stack->Add(h_T5HH1300[1]); h_T5HH1300_stack->Add(h_T5HH1300[2]);
  if (which!="DoubleB" && which!="DoubleB-H" && which!="tau32"){
    h_T5HH1300_stack->Add(h_T5HH1300[3]); h_T5HH1300_stack->Add(h_T5HH1300[4]); h_T5HH1300_stack->Add(h_T5HH1300[5]);
  }
  h_T5HH1300_stack->SetMarkerStyle(1); h_T5HH1300_stack->SetMarkerColor(kRed);
  h_T5HH1300_stack->SetLineColor(kRed); h_T5HH1300_stack->SetLineStyle(1); h_T5HH1300_stack->SetLineWidth(2);

  TH1D *h_T5HH1700_stack = (TH1D*)h_T5HH1700[0]->Clone("T5HH1700");
  h_T5HH1700_stack->Add(h_T5HH1700[1]); h_T5HH1700_stack->Add(h_T5HH1700[2]);
  if (which!="DoubleB" && which!="DoubleB-H" && which!="tau32"){
    h_T5HH1700_stack->Add(h_T5HH1700[3]); h_T5HH1700_stack->Add(h_T5HH1700[4]); h_T5HH1700_stack->Add(h_T5HH1700[5]);
  }
  h_T5HH1700_stack->SetMarkerStyle(1); h_T5HH1700_stack->SetMarkerColor(kYellow);
  h_T5HH1700_stack->SetLineColor(kYellow);

  TH1D *h_T5HH2100_stack = (TH1D*)h_T5HH2100[0]->Clone("T5HH2100");
  h_T5HH2100_stack->Add(h_T5HH2100[1]); h_T5HH2100_stack->Add(h_T5HH2100[2]);
  if (which!="DoubleB" && which!="DoubleB-H" && which!="tau32"){
    h_T5HH2100_stack->Add(h_T5HH2100[3]); h_T5HH2100_stack->Add(h_T5HH2100[4]); h_T5HH2100_stack->Add(h_T5HH2100[5]);
  }
  h_T5HH2100_stack->SetMarkerStyle(1); h_T5HH2100_stack->SetMarkerColor(kMagenta);
  h_T5HH2100_stack->SetLineColor(kMagenta);

  // h_ZJets_stack->GetXaxis()->SetTitle("double-B discriminator"); h_ZJets_stack->GetYaxis()->SetTitle("Events");
  // h_WJets_stack->GetXaxis()->SetTitle("double-B discriminator"); h_WJets_stack->GetYaxis()->SetTitle("Events");
  // h_TT_stack->GetXaxis()->SetTitle("double-B discriminator"); h_TT_stack->GetYaxis()->SetTitle("Events");
  // h_QCD_stack->GetXaxis()->SetTitle("double-B discriminator"); h_QCD_stack->GetYaxis()->SetTitle("Events");

  THStack * doubleBStack = new THStack("hs","");
  doubleBStack->Add(h_ZJets_stack);
  doubleBStack->Add(h_WJets_stack);
  doubleBStack->Add(h_TT_stack);
  doubleBStack->Add(h_QCD_stack);

  double W = 800;
  double H = 600;
  double T = 0.08*H;
  double B = 0.12*H;
  double L = 0.08*W;
  double R = 0.08*W;

  TString graphName = "";
  if (which =="doubleB" && useDeepDoubleB) graphName = "DeepDoubleB_Stack";
  else if (which == "doubleB") graphName = "DoubleB_Stack";
  else if (which == "tau32") graphName = "Tau32_Stack";
  else if (which == "leadmass") graphName = "LeadMass_Stack";
  else if (which == "subleadmass") graphName = "SubLeadMass_Stack";
  else if (which == "leadpt") graphName = "LeadPt_Stack";
  else if (which == "subleadpt") graphName = "SubLeadPt_Stack";

  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0);
  can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0);
  can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W );
  can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H );
  can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0);
  can_h->SetTicky(0);
  can_h->SetLogy();

  doubleBStack->Draw("hist");
  doubleBStack->GetYaxis()->SetTitle("Events");
  doubleBStack->SetMinimum(0.1);
  doubleBStack->SetMaximum(300);

  // if (which == "tau32") {
  //   doubleBStack->GetXaxis()->SetRangeUser(0,1);
  //   doubleBStack->SetMaximum(300);
  // }
  // if (which == "DoubleB-H") {
  //   doubleBStack->SetMaximum(300);
  // }


  if (which =="DoubleB" && useDeepDoubleB) doubleBStack->GetXaxis()->SetTitle("Deep double-B discriminator");
  else if (which == "DoubleB") doubleBStack->GetXaxis()->SetTitle("double-B discriminator");
  else if (which == "tau32") doubleBStack->GetXaxis()->SetTitle("Tau32");
  else if (which == "leadmass") doubleBStack->GetXaxis()->SetTitle("Lead jet soft-drop mass [GeV]");
  else if (which == "subleadmass") doubleBStack->GetXaxis()->SetTitle("Sub-lead jet soft-drop mass [GeV]");
  else if (which == "leadpt") doubleBStack->GetXaxis()->SetTitle("Lead jet p_{T} [GeV]");
  else if (which == "subleadpt") doubleBStack->GetXaxis()->SetTitle("Sub-lead jet p_{T} [GeV]");


  h_T5HH1300_stack->Draw("same l");
  h_T5HH1700_stack->Draw("same l");
  h_T5HH2100_stack->Draw("same l");

  TLegend* legend = new TLegend(0.35,0.8,0.9,0.9) ;
  legend->AddEntry(h_QCD_stack, "QCD", "f");
  legend->AddEntry(h_TT_stack, "TT", "f");
  legend->AddEntry(h_WJets_stack, "WJets", "f");
  legend->AddEntry(h_ZJets_stack, "ZJets", "f");
  legend->AddEntry(h_T5HH1300_stack, "T5HH1300", "l");
  legend->AddEntry(h_T5HH1700_stack, "T5HH1700", "l");
  legend->AddEntry(h_T5HH2100_stack, "T5HH2100", "l");
  legend->SetNColumns(4);
  legend->SetBorderSize(0);

  legend->Draw("same");
  can_h->Update();
  can_h->Modified();


  if (which =="DoubleB" && useDeepDoubleB) can_h->SaveAs(directory+"SupportingPlots/DeepDoubleBCompare_Stack.pdf");
  else if (which == "DoubleB") can_h->SaveAs(directory+"SupportingPlots/DoubleBCompare_Stack.pdf");
  else if (which == "tau32") can_h->SaveAs(directory+"SupportingPlots/Tau32Compare_Stack.pdf");
  else if (which == "leadmass") can_h->SaveAs(directory+"SupportingPlots/LeadJetMass_Stack.pdf");
  else if (which == "subleadmass") can_h->SaveAs(directory+"SupportingPlots/SubLeadJetMass_Stack.pdf");
  else if (which == "leadpt") can_h->SaveAs(directory+"SupportingPlots/LeadJetPt_Stack.pdf");
  else if (which == "subleadpt") can_h->SaveAs(directory+"SupportingPlots/SubLeadJetPt_Stack.pdf");
}

void makeSingleLeptStackPlot(vector<TH1F*> h_TT,vector<TH1F*> h_WJets, TString which = ""){
  TH1D *h_TT_stack = (TH1D*)h_TT[0]->Clone("TTAll");
  h_TT_stack->Add(h_TT[1]); h_TT_stack->Add(h_TT[2]);
  if (which!="DoubleB" && which!="DoubleB-H" && which!="tau32"){
    h_TT_stack->Add(h_TT[3]); h_TT_stack->Add(h_TT[4]); h_TT_stack->Add(h_TT[5]);
  }
  h_TT_stack->SetFillColor(kCyan);
  h_TT_stack->SetMarkerStyle(21); h_TT_stack->SetMarkerColor(kCyan);

  TH1D *h_WJets_stack = (TH1D*)h_WJets[0]->Clone("WJetsAll");
  h_WJets_stack->Add(h_WJets[1]); h_WJets_stack->Add(h_WJets[2]);
  if (which!="DoubleB" && which!="DoubleB-H" && which!="tau32"){
    h_WJets_stack->Add(h_WJets[3]); h_WJets_stack->Add(h_WJets[4]); h_WJets_stack->Add(h_WJets[5]);
  }
  h_WJets_stack->SetFillColor(kBlue);
  h_WJets_stack->SetMarkerStyle(21); h_WJets_stack->SetMarkerColor(kBlue);


  THStack * doubleBStack = new THStack("hs","");
  doubleBStack->Add(h_WJets_stack);
  doubleBStack->Add(h_TT_stack);

  double W = 800;
  double H = 600;
  double T = 0.08*H;
  double B = 0.12*H;
  double L = 0.08*W;
  double R = 0.08*W;

  TString graphName = "";
  if (which =="DoubleB" && useDeepDoubleB) graphName = "DeepDoubleB_Stack";
  else if (which == "DoubleB") graphName = "DoubleB_Stack";
  else if (which == "tau32") graphName = "Tau32_Stack";
  else if (which == "leadmass") graphName = "LeadMass_Stack";
  else if (which == "subleadmass") graphName = "SubLeadMass_Stack";
  else if (which == "leadpt") graphName = "LeadPt_Stack";
  else if (which == "subleadpt") graphName = "SubLeadPt_Stack";

  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0);
  can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0);
  can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W );
  can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H );
  can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0);
  can_h->SetTicky(0);
  can_h->SetLogy();

  doubleBStack->Draw("hist");
  doubleBStack->GetYaxis()->SetTitle("Events");
  doubleBStack->SetMinimum(0.1);
  doubleBStack->SetMaximum(300);

  // if (which == "tau32") {
  //   doubleBStack->GetXaxis()->SetRangeUser(0,1);
  //   doubleBStack->SetMaximum(300);
  // }
  // if (which == "DoubleB-H") {
  //   doubleBStack->SetMaximum(300);
  // }


  if (which =="DoubleB" && useDeepDoubleB) doubleBStack->GetXaxis()->SetTitle("Deep double-B discriminator");
  else if (which == "DoubleB") doubleBStack->GetXaxis()->SetTitle("double-B discriminator");
  else if (which == "tau32") doubleBStack->GetXaxis()->SetTitle("Tau32");
  else if (which == "leadmass") doubleBStack->GetXaxis()->SetTitle("Lead jet soft-drop mass [GeV]");
  else if (which == "subleadmass") doubleBStack->GetXaxis()->SetTitle("Sub-lead jet soft-drop mass [GeV]");
  else if (which == "leadpt") doubleBStack->GetXaxis()->SetTitle("Lead jet p_{T} [GeV]");
  else if (which == "subleadpt") doubleBStack->GetXaxis()->SetTitle("Sub-lead jet p_{T} [GeV]");

  TLegend* legend = new TLegend(0.35,0.8,0.9,0.9) ;
  legend->AddEntry(h_TT_stack, "TT", "f");
  legend->AddEntry(h_WJets_stack, "WJets", "f");
  legend->SetNColumns(4);
  legend->SetBorderSize(0);

  legend->Draw("same");
  can_h->Update();
  can_h->Modified();


  if (which =="DoubleB" && useDeepDoubleB) can_h->SaveAs(directory+"SupportingPlots/SingleLeptDeepDoubleBCompare_Stack.pdf");
  else if (which == "DoubleB") can_h->SaveAs(directory+"SupportingPlots/SingleLeptDoubleBCompare_Stack.pdf");
  else if (which == "tau32") can_h->SaveAs(directory+"SupportingPlots/SingleLeptTau32Compare_Stack.pdf");
  else if (which == "leadmass") can_h->SaveAs(directory+"SupportingPlots/SingleLeptLeadJetMass_Stack.pdf");
  else if (which == "subleadmass") can_h->SaveAs(directory+"SupportingPlots/SingleLeptSubLeadJetMass_Stack.pdf");
  else if (which == "leadpt") can_h->SaveAs(directory+"SupportingPlots/SingleLeptLeadJetPt_Stack.pdf");
  else if (which == "subleadpt") can_h->SaveAs(directory+"SupportingPlots/SingleLeptSubLeadJetPt_Stack.pdf");
}

void makeDoubleBStack(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets, TH1F* h_T5HH1300,TH1F* h_T5HH1700,TH1F* h_T5HH2100, TString which = ""){
  h_QCD->SetFillColor(kGray); h_QCD->SetMarkerStyle(21); h_QCD->SetMarkerColor(kGray);
  h_TT->SetFillColor(kCyan); h_TT->SetMarkerStyle(21); h_TT->SetMarkerColor(kCyan);
  h_WJets->SetFillColor(kBlue);h_WJets->SetMarkerStyle(21); h_WJets->SetMarkerColor(kBlue);
  h_ZJets->SetFillColor(kGreen+2);h_ZJets->SetMarkerStyle(21); h_ZJets->SetMarkerColor(kGreen+2);

  h_T5HH1300->SetMarkerStyle(1); h_T5HH1300->SetMarkerColor(kRed);
  h_T5HH1300->SetLineColor(kRed); h_T5HH1300->SetLineStyle(1); h_T5HH1300->SetLineWidth(2);
  h_T5HH1700->SetMarkerStyle(1); h_T5HH1700->SetMarkerColor(kYellow); h_T5HH1700->SetLineColor(kYellow);
  h_T5HH2100->SetMarkerStyle(1); h_T5HH2100->SetMarkerColor(kMagenta); h_T5HH2100->SetLineColor(kMagenta);


  THStack * Stack = new THStack("hs","");
  Stack->Add(h_ZJets);
  Stack->Add(h_WJets);
  Stack->Add(h_TT);
  Stack->Add(h_QCD);

  double W = 800;
  double H = 600;
  double T = 0.08*H;
  double B = 0.12*H;
  double L = 0.08*W;
  double R = 0.08*W;

  TString graphName = "";
  if (which =="doubleAll") graphName = "DoubleB_allMET_DoubleSR";
  else if (which =="doubleLow") graphName = "DoubleB_lowMET_DoubleSR";
  else if (which =="doubleHigh") graphName = "DoubleB_highMET_DoubleSR";
  else if (which == "singleAll") graphName = "DoubleB_allMET_SingleSR";
  else if (which == "singleLow") graphName = "DoubleB_lowMET_SingleSR";
  else if (which == "singleHigh") graphName = "DoubleB_highMET_SingleSR";


  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0);
  can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0);
  can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W );
  can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H );
  can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0);
  can_h->SetTicky(0);
  can_h->SetLogy();

  Stack->Draw("hist");
  Stack->GetYaxis()->SetTitle("Events");
  Stack->SetMinimum(0.1);
  Stack->SetMaximum(300);

  // if (which == "all") {
  //   METStack->SetMinimum(.5);
  //   METStack->SetMaximum(2000);
  // }

  Stack->GetXaxis()->SetTitle("Double-b discriminator");

  h_T5HH1300->Draw("same l");
  h_T5HH1700->Draw("same l");
  h_T5HH2100->Draw("same l");

  TLegend* legend = new TLegend(0.35,0.8,0.9,0.9) ;
  legend->AddEntry(h_QCD, "QCD", "f");
  legend->AddEntry(h_TT, "TT", "f");
  legend->AddEntry(h_WJets, "WJets", "f");
  legend->AddEntry(h_ZJets, "ZJets", "f");
  legend->AddEntry(h_T5HH1300, "T5HH1300", "l");
  legend->AddEntry(h_T5HH1700, "T5HH1700", "l");
  legend->AddEntry(h_T5HH2100, "T5HH2100", "l");
  legend->SetNColumns(4);
  legend->SetBorderSize(0);

  legend->Draw("same");
  can_h->Update();
  can_h->Modified();

  TString savename = directory+"SupportingPlots/DoubleB"+which+"_Stack.pdf";
  can_h->SaveAs(savename);

}

void makeMETStack(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets, TH1F* h_T5HH1300,TH1F* h_T5HH1700,TH1F* h_T5HH2100, TString which = ""){
  h_QCD->SetFillColor(kGray); h_QCD->SetMarkerStyle(21); h_QCD->SetMarkerColor(kGray);
  h_TT->SetFillColor(kCyan); h_TT->SetMarkerStyle(21); h_TT->SetMarkerColor(kCyan);
  h_WJets->SetFillColor(kBlue);h_WJets->SetMarkerStyle(21); h_WJets->SetMarkerColor(kBlue);
  h_ZJets->SetFillColor(kGreen+2);h_ZJets->SetMarkerStyle(21); h_ZJets->SetMarkerColor(kGreen+2);

  h_T5HH1300->SetMarkerStyle(1); h_T5HH1300->SetMarkerColor(kRed);
  h_T5HH1300->SetLineColor(kRed); h_T5HH1300->SetLineStyle(1); h_T5HH1300->SetLineWidth(2);
  h_T5HH1700->SetMarkerStyle(1); h_T5HH1700->SetMarkerColor(kYellow); h_T5HH1700->SetLineColor(kYellow);
  h_T5HH2100->SetMarkerStyle(1); h_T5HH2100->SetMarkerColor(kMagenta); h_T5HH2100->SetLineColor(kMagenta);


  THStack * METStack = new THStack("hs","");
  METStack->Add(h_ZJets);
  METStack->Add(h_WJets);
  METStack->Add(h_TT);
  METStack->Add(h_QCD);

  double W = 800;
  double H = 600;
  double T = 0.08*H;
  double B = 0.12*H;
  double L = 0.08*W;
  double R = 0.08*W;

  TString graphName = "";
  if (which =="all") graphName = "MET_Stack";
  else if (which == "single") graphName = "MET_Stack_Single";
  else if (which == "double") graphName = "MET_Stack_Double";



  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0);
  can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0);
  can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W );
  can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H );
  can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0);
  can_h->SetTicky(0);
  can_h->SetLogy();

  METStack->Draw("hist");
  METStack->GetYaxis()->SetTitle("Events");
  METStack->SetMinimum(0.1);
  METStack->SetMaximum(300);

  if (which == "single") {
    METStack->SetMaximum(2000);
  }

  // if (which == "tau32") {
  //   doubleBStack->GetXaxis()->SetRangeUser(0,1);
  //   doubleBStack->SetMaximum(300);
  // }
  // if (which == "DoubleB-H") {
  //   doubleBStack->SetMaximum(300);
  // }

  METStack->GetXaxis()->SetTitle("MET [GeV]");


  h_T5HH1300->Draw("same");
  h_T5HH1700->Draw("same");
  h_T5HH2100->Draw("same");

  TLegend* legend = new TLegend(0.35,0.8,0.9,0.9) ;
  legend->AddEntry(h_QCD, "QCD", "f");
  legend->AddEntry(h_TT, "TT", "f");
  legend->AddEntry(h_WJets, "WJets", "f");
  legend->AddEntry(h_ZJets, "ZJets", "f");
  legend->AddEntry(h_T5HH1300, "T5HH1300", "l");
  legend->AddEntry(h_T5HH1700, "T5HH1700", "l");
  legend->AddEntry(h_T5HH2100, "T5HH2100", "l");
  legend->SetNColumns(4);
  legend->SetBorderSize(0);

  legend->Draw("same");
  can_h->Update();
  can_h->Modified();

  TString savename = directory+"SupportingPlots/MET"+which+"_Stack.pdf";
  if (useDeepDoubleB) savename = directory+"SupportingPlots/MET"+which+"_deepDB_Stack.pdf";
  can_h->SaveAs(savename);

  // if (which =="all") can_h->SaveAs("SupportingPlots/MET_Stack.pdf");
  // else if (which == "single") can_h->SaveAs("SupportingPlots/METSingle_Stack.pdf");
  // else if (which == "double") can_h->SaveAs("SupportingPlots/METDouble_Stack.pdf");
}


void QuickROC(TH1F* signalHist, TH1F* bkgHist, TString which = ""){
  TString graphName = "";
  if (which =="") graphName = "DoubleB";
  else if (which == "DoubleB-H") graphName = "DoubleBH";
  else if (which == "tau32") graphName = "Tau32";

  double W = 800;
  double H = 600;
  double T = 0.08*H;
  double B = 0.12*H;
  double L = 0.08*W;
  double R = 0.08*W;

  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);

  can_h->SetFillColor(0);
  can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0);
  can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W );
  can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H );
  can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0);
  can_h->SetTicky(0);

  signalHist->SetLineColor(kRed);

  TLegend* leg2=new TLegend(0.6,0.65,0.88,0.85);
  leg2->AddEntry(signalHist, "T5HH", "F");
  leg2->AddEntry(bkgHist, "TT", "F");

  int numBins = signalHist->GetNbinsX();
  float totalNumSig = signalHist->Integral();
  float totalNumBkg = bkgHist->Integral();

  double x[100], y[100];
  // std::cout<<"num bins sig: "<< numBins<<", num bins bkg: "<< bkgHist->GetNbinsX()<<std::endl;
  for(int i=1; i<=numBins; ++i){
    int cutbin=bkgHist->FindBin(bkgHist->GetBinLowEdge(i));
    // std::cout<<"Cutbin: "<<cutbin<<std::endl;
    double bkg=bkgHist->Integral(1, cutbin);
    double sig=signalHist->Integral(1,cutbin);
    double sigEff = 1 - sig/totalNumSig;
    double bkgRej = bkg/totalNumBkg;

    std::cout<<" Cut Value "<<bkgHist->GetBinLowEdge(i+1)<<": ";
    // std::cout<<"     Total Bkg: "<<bkgHist->Integral(1, cutbin)<<" Total Sig: "<<signalHist->Integral(1,cutbin)<<std::endl;
    // std::cout<<"     Bkg Rejection: "<<bkgRej<<" Sig Eff: "<<sigEff<<std::endl;
    std::cout<<"     SigEff,BkgRej: "<<sigEff<<", "<<bkgRej<<endl;
    x[i] = sigEff;
    y[i] = bkgRej;
  }
  TGraph* gr = new TGraph(numBins,x,y);

  gr->SetLineColor(kRed);

  gr->Draw();
  if (which =="")  gr->SetTitle("MVA double-b discriminator");
  else if (which == "DoubleB-H")  gr->SetTitle("Deep double-b discriminator");
  gr->GetXaxis()->SetTitle("Signal efficiency (T5HH1300)");
  gr->GetYaxis()->SetTitle("Background rejection (TT)");
  gr->GetYaxis()->SetRangeUser(0.5,1.0);
  can_h->Modified();
  can_h->Update();


  if (which =="") can_h->SaveAs(directory+"SupportingPlots/DoubleB_ROC_TT.pdf");
  else if (which == "DoubleB-H") can_h->SaveAs(directory+"SupportingPlots/DoubleBH_ROC_TT.pdf");
  // else if (which == "tau32") can_h->SaveAs("SupportingPlots/Tau32_ROC.pdf");
}

void ABCD2(){
  TFile * f = TFile::Open("ALPHABEThistos_V12_METonly.root");
  TFile * fSignal = TFile::Open("ALPHABEThistos_V12_METonly_signalOnly.root");


  TH1F* h_MET_double_QCD = (TH1F*)f->Get("MET_doubletagSR_QCD");
  TH1F* h_MET_double_TT = (TH1F*)f->Get("MET_doubletagSR_TT");
  TH1F* h_MET_double_WJets = (TH1F*)f->Get("MET_doubletagSR_WJets");
  TH1F* h_MET_double_ZJets = (TH1F*)f->Get("MET_doubletagSR_ZJets");


  TH1F* h_MET_double_TChiHH400 = (TH1F*)f->Get("MET_doubletagSR_TChiHH400");

  // TH1F* h_MET_double_TChiHH425 = (TH1F*)f->Get("MET_doubletagSR_TChiHH425");
  // TH1F* h_MET_double_TChiHH450 = (TH1F*)f->Get("MET_doubletagSR_TChiHH450");
  // TH1F* h_MET_double_TChiHH475 = (TH1F*)f->Get("MET_doubletagSR_TChiHH475");
  // TH1F* h_MET_double_TChiHH500 = (TH1F*)f->Get("MET_doubletagSR_TChiHH500");
  // TH1F* h_MET_double_TChiHH525 = (TH1F*)f->Get("MET_doubletagSR_TChiHH525");
  // TH1F* h_MET_double_TChiHH550 = (TH1F*)f->Get("MET_doubletagSR_TChiHH550");
  // TH1F* h_MET_double_TChiHH575 = (TH1F*)f->Get("MET_doubletagSR_TChiHH575");
  // TH1F* h_MET_double_TChiHH600 = (TH1F*)f->Get("MET_doubletagSR_TChiHH600");
  // TH1F* h_MET_double_TChiHH625 = (TH1F*)f->Get("MET_doubletagSR_TChiHH625");
  // TH1F* h_MET_double_TChiHH650 = (TH1F*)f->Get("MET_doubletagSR_TChiHH650");
  // TH1F* h_MET_double_TChiHH675 = (TH1F*)f->Get("MET_doubletagSR_TChiHH675");
  TH1F* h_MET_double_TChiHH700 = (TH1F*)f->Get("MET_doubletagSR_TChiHH700");
  // TH1F* h_MET_double_TChiHH725 = (TH1F*)f->Get("MET_doubletagSR_TChiHH725");
  // TH1F* h_MET_double_TChiHH750 = (TH1F*)f->Get("MET_doubletagSR_TChiHH750");
  // TH1F* h_MET_double_TChiHH775 = (TH1F*)f->Get("MET_doubletagSR_TChiHH775");
  // TH1F* h_MET_double_TChiHH800 = (TH1F*)f->Get("MET_doubletagSR_TChiHH800");
  // TH1F* h_MET_double_TChiHH825 = (TH1F*)f->Get("MET_doubletagSR_TChiHH825");
  // TH1F* h_MET_double_TChiHH850 = (TH1F*)f->Get("MET_doubletagSR_TChiHH850");
  // TH1F* h_MET_double_TChiHH875 = (TH1F*)f->Get("MET_doubletagSR_TChiHH875");
  // TH1F* h_MET_double_TChiHH900 = (TH1F*)f->Get("MET_doubletagSR_TChiHH900");
  // TH1F* h_MET_double_TChiHH925 = (TH1F*)f->Get("MET_doubletagSR_TChiHH925");
  // TH1F* h_MET_double_TChiHH950 = (TH1F*)f->Get("MET_doubletagSR_TChiHH950");
  // TH1F* h_MET_double_TChiHH975 = (TH1F*)f->Get("MET_doubletagSR_TChiHH975");
  TH1F* h_MET_double_TChiHH1000 = (TH1F*)f->Get("MET_doubletagSR_TChiHH1000");

  TH1F* h_MET_single_QCD = (TH1F*)f->Get("MET_tagSR_QCD");
  TH1F* h_MET_single_TT = (TH1F*)f->Get("MET_tagSR_TT");
  TH1F* h_MET_single_WJets = (TH1F*)f->Get("MET_tagSR_WJets");
  TH1F* h_MET_single_ZJets = (TH1F*)f->Get("MET_tagSR_ZJets");

  TH1F* h_MET_single_TChiHH400 = (TH1F*)f->Get("MET_tagSR_TChiHH400");

  // TH1F* h_MET_single_TChiHH425 = (TH1F*)f->Get("MET_tagSR_TChiHH425");
  // TH1F* h_MET_single_TChiHH450 = (TH1F*)f->Get("MET_tagSR_TChiHH450");
  // TH1F* h_MET_single_TChiHH475 = (TH1F*)f->Get("MET_tagSR_TChiHH475");
  // TH1F* h_MET_single_TChiHH500 = (TH1F*)f->Get("MET_tagSR_TChiHH500");
  // TH1F* h_MET_single_TChiHH525 = (TH1F*)f->Get("MET_tagSR_TChiHH525");
  // TH1F* h_MET_single_TChiHH550 = (TH1F*)f->Get("MET_tagSR_TChiHH550");
  // TH1F* h_MET_single_TChiHH575 = (TH1F*)f->Get("MET_tagSR_TChiHH575");
  // TH1F* h_MET_single_TChiHH600 = (TH1F*)f->Get("MET_tagSR_TChiHH600");
  // TH1F* h_MET_single_TChiHH625 = (TH1F*)f->Get("MET_tagSR_TChiHH625");
  // TH1F* h_MET_single_TChiHH650 = (TH1F*)f->Get("MET_tagSR_TChiHH650");
  // TH1F* h_MET_single_TChiHH675 = (TH1F*)f->Get("MET_tagSR_TChiHH675");
  TH1F* h_MET_single_TChiHH700 = (TH1F*)f->Get("MET_tagSR_TChiHH700");
  // TH1F* h_MET_single_TChiHH725 = (TH1F*)f->Get("MET_tagSR_TChiHH725");
  // TH1F* h_MET_single_TChiHH750 = (TH1F*)f->Get("MET_tagSR_TChiHH750");
  // TH1F* h_MET_single_TChiHH775 = (TH1F*)f->Get("MET_tagSR_TChiHH775");
  // TH1F* h_MET_single_TChiHH800 = (TH1F*)f->Get("MET_tagSR_TChiHH800");
  // TH1F* h_MET_single_TChiHH825 = (TH1F*)f->Get("MET_tagSR_TChiHH825");
  // TH1F* h_MET_single_TChiHH850 = (TH1F*)f->Get("MET_tagSR_TChiHH850");
  // TH1F* h_MET_single_TChiHH875 = (TH1F*)f->Get("MET_tagSR_TChiHH875");
  // TH1F* h_MET_single_TChiHH900 = (TH1F*)f->Get("MET_tagSR_TChiHH900");
  // TH1F* h_MET_single_TChiHH925 = (TH1F*)f->Get("MET_tagSR_TChiHH925");
  // TH1F* h_MET_single_TChiHH950 = (TH1F*)f->Get("MET_tagSR_TChiHH950");
  // TH1F* h_MET_single_TChiHH975 = (TH1F*)f->Get("MET_tagSR_TChiHH975");
  TH1F* h_MET_single_TChiHH1000 = (TH1F*)f->Get("MET_tagSR_TChiHH1000");

  vector<float> double_MET1_yields; vector<float> double_MET2_yields; vector<float> double_MET3_yields; vector<float> double_MET4_yields; vector<float> double_MET5_yields;
  vector<float> single_MET1_yields; vector<float> single_MET2_yields; vector<float> single_MET3_yields; vector<float> single_MET4_yields; vector<float> single_MET5_yields;

  vector<TH1F*> double_hists = {h_MET_double_QCD, h_MET_double_TT, h_MET_double_WJets, h_MET_double_ZJets, h_MET_double_TChiHH400, h_MET_double_TChiHH700, h_MET_double_TChiHH1000};
  vector<TH1F*> single_hists = {h_MET_single_QCD, h_MET_single_TT, h_MET_single_WJets, h_MET_single_ZJets, h_MET_single_TChiHH400, h_MET_single_TChiHH700, h_MET_single_TChiHH1000};
  for (int i = 0; i<7; i++) { //loop through samples
    double_MET1_yields.push_back(double_hists[i]->GetBinContent(1)); single_MET1_yields.push_back(single_hists[i]->GetBinContent(1));
    double_MET2_yields.push_back(double_hists[i]->GetBinContent(2)); single_MET2_yields.push_back(single_hists[i]->GetBinContent(2));
    double_MET3_yields.push_back(double_hists[i]->GetBinContent(3)); single_MET3_yields.push_back(single_hists[i]->GetBinContent(3));
    double_MET4_yields.push_back(double_hists[i]->GetBinContent(4)); single_MET4_yields.push_back(single_hists[i]->GetBinContent(4));
    double_MET5_yields.push_back(double_hists[i]->GetBinContent(5)); single_MET5_yields.push_back(single_hists[i]->GetBinContent(5));
  }

  myfile.open("MCyields_V12.txt");
  myfile << "\\hline Double \\\\ \\hline \n MET bin [GeV] & QCD & TT & WJets & ZJets & TChiHH400 & TChiHH700 & TChiHH1000 \\\\ \n \\hline";
  myfile << "150-200 ";
  for (int i=0;i<7;i++) myfile <<" & "<< double_MET1_yields[i];
  myfile << " \\\\ \n 200-300 ";
  for (int i=0;i<7;i++) myfile <<" & "<< double_MET2_yields[i];
  myfile << " \\\\ \n 300-500 ";
  for (int i=0;i<7;i++) myfile <<" & "<< double_MET3_yields[i];
  myfile << " \\\\ \n 500-700 ";
  for (int i=0;i<7;i++) myfile <<" & "<< double_MET4_yields[i];
  myfile << " \\\\ \n $>$ 700 ";
  for (int i=0;i<7;i++) myfile <<" & "<< double_MET5_yields[i];
  myfile << " \\\\ \\hline \\hline \n Single \\\\ \\hline \n 150-200 ";


  for (int i=0;i<7;i++) myfile <<" & "<< single_MET1_yields[i];
  myfile << " \\\\ \n 200-300 ";
  for (int i=0;i<7;i++) myfile <<" & "<< single_MET2_yields[i];
  myfile << " \\\\ \n 300-500 ";
  for (int i=0;i<7;i++) myfile <<" & "<< single_MET3_yields[i];
  myfile << " \\\\ \n 500-700 ";
  for (int i=0;i<7;i++) myfile <<" & "<< single_MET4_yields[i];
  myfile << " \\\\ \n $>$ 700 ";
  for (int i=0;i<7;i++) myfile <<" & "<< single_MET5_yields[i];
  myfile << " \\\\ \\hline";

  myfile.close();

}
