#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
// #include "labelCMS.h"
// #include "analysisTools.h"
// #include <TGraphAsymmErrors.h>
//using namespace std;

void makeABCDPlot(vector<TH1F*> dem_histos, TString type = "", TString tagType = "");
void makeRpfPlot(vector<TH1F*> dem_histos, TString bkgType = "", TString jetType = "", TString tagType = "");

void makeDoubleBPlots(vector<TH1F*> dem_histos, TString bkgType = "");
void makeDoubleBStackPlot(vector<TH1F*> h_QCD,vector<TH1F*> h_TT,vector<TH1F*> h_WJets,vector<TH1F*> h_ZJets);

ofstream myfile;
bool getDeviation = false;
bool runABCDPlots = false;
bool runRPFPlots = false;
bool runCompareDoubleB = true;
bool isSignal = false;


void ABCD() {
  // TFile * f = TFile::Open("ALPHABEThistos_V17.root");
  // TFile * f = TFile::Open("ALPHABEThistos_signal_V17.root");
  TFile * f = TFile::Open("ALPHABEThistos_1l_V17.root");

  TH1F * h_A_sum; TH1F * h_B_sum; TH1F * h_A1_sum; TH1F * h_B1_sum; TH1F *h_C_sum; TH1F *h_D_sum;
  TH1F * h_A_QCD; TH1F * h_B_QCD; TH1F * h_A1_QCD; TH1F * h_B1_QCD; TH1F * h_C_QCD; TH1F * h_D_QCD;
  TH1F * h_A_TT; TH1F * h_B_TT; TH1F * h_A1_TT; TH1F * h_B1_TT; TH1F * h_C_TT; TH1F * h_D_TT;
  TH1F * h_A_WJets; TH1F * h_B_WJets; TH1F * h_A1_WJets; TH1F * h_B1_WJets; TH1F * h_C_WJets; TH1F * h_D_WJets;
  TH1F * h_A_ZJets; TH1F * h_B_ZJets; TH1F * h_A1_ZJets; TH1F * h_B1_ZJets; TH1F * h_C_ZJets; TH1F * h_D_ZJets;

  TH1F * h_J1M_doubletagSR_sum; TH1F * h_J2M_doubletagSR_sum; TH1F * h_J1M_doubletagSB_sum; TH1F * h_J2M_doubletagSB_sum;
  TH1F * h_J1M_antitagSR_sum; TH1F * h_J2M_antitagSR_sum; TH1F * h_J1M_antitagSB_sum; TH1F * h_J2M_antitagSB_sum;
  TH1F * h_J1M_tagSR_sum; TH1F * h_J2M_tagSR_sum; TH1F * h_J1M_tagSB_sum; TH1F * h_J2M_tagSB_sum;

  TH1F * h_J1M_doubletagSR_QCD; TH1F * h_J2M_doubletagSR_QCD; TH1F * h_J1M_doubletagSB_QCD; TH1F * h_J2M_doubletagSB_QCD;
  TH1F * h_J1M_antitagSR_QCD; TH1F * h_J2M_antitagSR_QCD; TH1F * h_J1M_antitagSB_QCD; TH1F * h_J2M_antitagSB_QCD;
  TH1F * h_J1M_tagSR_QCD; TH1F * h_J2M_tagSR_QCD; TH1F * h_J1M_tagSB_QCD; TH1F * h_J2M_tagSB_QCD;

  TH1F * h_J1M_doubletagSR_TT; TH1F * h_J2M_doubletagSR_TT; TH1F * h_J1M_doubletagSB_TT; TH1F * h_J2M_doubletagSB_TT;
  TH1F * h_J1M_antitagSR_TT; TH1F * h_J2M_antitagSR_TT; TH1F * h_J1M_antitagSB_TT; TH1F * h_J2M_antitagSB_TT;
  TH1F * h_J1M_tagSR_TT; TH1F * h_J2M_tagSR_TT; TH1F * h_J1M_tagSB_TT; TH1F * h_J2M_tagSB_TT;

  TH1F * h_J1M_doubletagSR_WJets; TH1F * h_J2M_doubletagSR_WJets; TH1F * h_J1M_doubletagSB_WJets; TH1F * h_J2M_doubletagSB_WJets;
  TH1F * h_J1M_antitagSR_WJets; TH1F * h_J2M_antitagSR_WJets; TH1F * h_J1M_antitagSB_WJets; TH1F * h_J2M_antitagSB_WJets;
  TH1F * h_J1M_tagSR_WJets; TH1F * h_J2M_tagSR_WJets; TH1F * h_J1M_tagSB_WJets; TH1F * h_J2M_tagSB_WJets;

  TH1F * h_J1M_doubletagSR_ZJets; TH1F * h_J2M_doubletagSR_ZJets; TH1F * h_J1M_doubletagSB_ZJets; TH1F * h_J2M_doubletagSB_ZJets;
  TH1F * h_J1M_antitagSR_ZJets; TH1F * h_J2M_antitagSR_ZJets; TH1F * h_J1M_antitagSB_ZJets; TH1F * h_J2M_antitagSB_ZJets;
  TH1F * h_J1M_tagSR_ZJets; TH1F * h_J2M_tagSR_ZJets; TH1F * h_J1M_tagSB_ZJets; TH1F * h_J2M_tagSB_ZJets;

  TH1F * h_doubleBSBLow_WJets; TH1F * h_doubleBSR_WJets; TH1F * h_doubleBSBHigh_WJets;
  TH1F * h_doubleBSBLow_ZJets; TH1F * h_doubleBSR_ZJets; TH1F * h_doubleBSBHigh_ZJets;
  TH1F * h_doubleBSBLow_TT; TH1F * h_doubleBSR_TT; TH1F * h_doubleBSBHigh_TT;
  TH1F * h_doubleBSBLow_QCD; TH1F * h_doubleBSR_QCD; TH1F * h_doubleBSBHigh_QCD;

  if (isSignal){
    h_A_sum = (TH1F*)f->Get("MET_doubletagSR_sum"); h_B_sum = (TH1F*)f->Get("MET_doubletagSB_sum");
    h_A1_sum = (TH1F*)f->Get("MET_tagSR_sum"); h_B1_sum = (TH1F*)f->Get("MET_tagSB_sum");
    h_C_sum = (TH1F*)f->Get("MET_antitagSR_sum"); h_D_sum = (TH1F*)f->Get("MET_antitagSB_sum");

    h_A_QCD = (TH1F*)f->Get("MET_doubletagSR_QCD"); h_B_QCD = (TH1F*)f->Get("MET_doubletagSB_QCD");
    h_A1_QCD = (TH1F*)f->Get("MET_tagSR_QCD"); h_B1_QCD = (TH1F*)f->Get("MET_tagSB_QCD");
    h_C_QCD = (TH1F*)f->Get("MET_antitagSR_QCD"); h_D_QCD = (TH1F*)f->Get("MET_antitagSB_QCD");

    h_A_TT = (TH1F*)f->Get("MET_doubletagSR_TT"); h_B_TT = (TH1F*)f->Get("MET_doubletagSB_TT");
    h_A1_TT = (TH1F*)f->Get("MET_tagSR_TT"); h_B1_TT = (TH1F*)f->Get("MET_tagSB_TT");
    h_C_TT = (TH1F*)f->Get("MET_antitagSR_TT"); h_D_TT = (TH1F*)f->Get("MET_antitagSB_TT");

    h_A_WJets = (TH1F*)f->Get("MET_doubletagSR_WJets"); h_B_WJets = (TH1F*)f->Get("MET_doubletagSB_WJets");
    h_A1_WJets = (TH1F*)f->Get("MET_tagSR_WJets"); h_B1_WJets = (TH1F*)f->Get("MET_tagSB_WJets");
    h_C_WJets = (TH1F*)f->Get("MET_antitagSR_WJets");   h_D_WJets = (TH1F*)f->Get("MET_antitagSB_WJets");

    h_A_ZJets = (TH1F*)f->Get("MET_doubletagSR_ZJets"); h_B_ZJets = (TH1F*)f->Get("MET_doubletagSB_ZJets");
    h_A1_ZJets = (TH1F*)f->Get("MET_tagSR_ZJets"); h_B1_ZJets = (TH1F*)f->Get("MET_tagSB_ZJets");
    h_C_ZJets = (TH1F*)f->Get("MET_antitagSR_ZJets");   h_D_ZJets = (TH1F*)f->Get("MET_antitagSB_ZJets");

    h_J1M_doubletagSR_sum = (TH1F*)f->Get("J1pt_M_doubletagSR_sum"); h_J2M_doubletagSR_sum = (TH1F*)f->Get("J2pt_M_doubletagSR_sum");
    h_J1M_doubletagSB_sum = (TH1F*)f->Get("J1pt_M_doubletagSB_sum"); h_J2M_doubletagSB_sum = (TH1F*)f->Get("J2pt_M_doubletagSB_sum");
    h_J1M_antitagSR_sum = (TH1F*)f->Get("J1pt_M_antitagSR_sum"); h_J2M_antitagSR_sum = (TH1F*)f->Get("J2pt_M_antitagSR_sum");
    h_J1M_antitagSB_sum = (TH1F*)f->Get("J1pt_M_antitagSB_sum"); h_J2M_antitagSB_sum = (TH1F*)f->Get("J2pt_M_antitagSB_sum");
    h_J1M_tagSR_sum = (TH1F*)f->Get("J1pt_M_tagSR_sum"); h_J2M_tagSR_sum = (TH1F*)f->Get("J2pt_M_tagSR_sum");
    h_J1M_tagSB_sum = (TH1F*)f->Get("J1pt_M_tagSB_sum"); h_J2M_tagSB_sum = (TH1F*)f->Get("J2pt_M_tagSB_sum");

    h_J1M_doubletagSR_QCD = (TH1F*)f->Get("J1pt_M_doubletagSR_QCD"); h_J2M_doubletagSR_QCD = (TH1F*)f->Get("J2pt_M_doubletagSR_QCD");
    h_J1M_doubletagSB_QCD = (TH1F*)f->Get("J1pt_M_doubletagSB_QCD"); h_J2M_doubletagSB_QCD = (TH1F*)f->Get("J2pt_M_doubletagSB_QCD");
    h_J1M_antitagSR_QCD = (TH1F*)f->Get("J1pt_M_antitagSR_QCD"); h_J2M_antitagSR_QCD = (TH1F*)f->Get("J2pt_M_antitagSR_QCD");
    h_J1M_antitagSB_QCD = (TH1F*)f->Get("J1pt_M_antitagSB_QCD"); h_J2M_antitagSB_QCD = (TH1F*)f->Get("J2pt_M_antitagSB_QCD");
    h_J1M_tagSR_QCD = (TH1F*)f->Get("J1pt_M_tagSR_QCD"); h_J2M_tagSR_QCD = (TH1F*)f->Get("J2pt_M_tagSR_QCD");
    h_J1M_tagSB_QCD = (TH1F*)f->Get("J1pt_M_tagSB_QCD"); h_J2M_tagSB_QCD = (TH1F*)f->Get("J2pt_M_tagSB_QCD");

    h_J1M_doubletagSR_TT = (TH1F*)f->Get("J1pt_M_doubletagSR_TT"); h_J2M_doubletagSR_TT = (TH1F*)f->Get("J2pt_M_doubletagSR_TT");
    h_J1M_doubletagSB_TT = (TH1F*)f->Get("J1pt_M_doubletagSB_TT"); h_J2M_doubletagSB_TT = (TH1F*)f->Get("J2pt_M_doubletagSB_TT");
    h_J1M_antitagSR_TT = (TH1F*)f->Get("J1pt_M_antitagSR_TT"); h_J2M_antitagSR_TT = (TH1F*)f->Get("J2pt_M_antitagSR_TT");
    h_J1M_antitagSB_TT = (TH1F*)f->Get("J1pt_M_antitagSB_TT"); h_J2M_antitagSB_TT = (TH1F*)f->Get("J2pt_M_antitagSB_TT");
    h_J1M_tagSR_TT = (TH1F*)f->Get("J1pt_M_tagSR_TT"); h_J2M_tagSR_TT = (TH1F*)f->Get("J2pt_M_tagSR_TT");
    h_J1M_tagSB_TT = (TH1F*)f->Get("J1pt_M_tagSB_TT"); h_J2M_tagSB_TT = (TH1F*)f->Get("J2pt_M_tagSB_TT");

    h_J1M_doubletagSR_WJets = (TH1F*)f->Get("J1pt_M_doubletagSR_WJets"); h_J2M_doubletagSR_WJets = (TH1F*)f->Get("J2pt_M_doubletagSR_WJets");
    h_J1M_doubletagSB_WJets = (TH1F*)f->Get("J1pt_M_doubletagSB_WJets"); h_J2M_doubletagSB_WJets = (TH1F*)f->Get("J2pt_M_doubletagSB_WJets");
    h_J1M_antitagSR_WJets = (TH1F*)f->Get("J1pt_M_antitagSR_WJets"); h_J2M_antitagSR_WJets = (TH1F*)f->Get("J2pt_M_antitagSR_WJets");
    h_J1M_antitagSB_WJets = (TH1F*)f->Get("J1pt_M_antitagSB_WJets"); h_J2M_antitagSB_WJets = (TH1F*)f->Get("J2pt_M_antitagSB_WJets");
    h_J1M_tagSR_WJets = (TH1F*)f->Get("J1pt_M_tagSR_WJets"); h_J2M_tagSR_WJets = (TH1F*)f->Get("J2pt_M_tagSR_WJets");
    h_J1M_tagSB_WJets = (TH1F*)f->Get("J1pt_M_tagSB_WJets"); h_J2M_tagSB_WJets = (TH1F*)f->Get("J2pt_M_tagSB_WJets");

    h_J1M_doubletagSR_ZJets = (TH1F*)f->Get("J1pt_M_doubletagSR_ZJets"); h_J2M_doubletagSR_ZJets = (TH1F*)f->Get("J2pt_M_doubletagSR_ZJets");
    h_J1M_doubletagSB_ZJets = (TH1F*)f->Get("J1pt_M_doubletagSB_ZJets"); h_J2M_doubletagSB_ZJets = (TH1F*)f->Get("J2pt_M_doubletagSB_ZJets");
    h_J1M_antitagSR_ZJets = (TH1F*)f->Get("J1pt_M_antitagSR_ZJets"); h_J2M_antitagSR_ZJets = (TH1F*)f->Get("J2pt_M_antitagSR_ZJets");
    h_J1M_antitagSB_ZJets = (TH1F*)f->Get("J1pt_M_antitagSB_ZJets"); h_J2M_antitagSB_ZJets = (TH1F*)f->Get("J2pt_M_antitagSB_ZJets");
    h_J1M_tagSR_ZJets = (TH1F*)f->Get("J1pt_M_tagSR_ZJets"); h_J2M_tagSR_ZJets = (TH1F*)f->Get("J2pt_M_tagSR_ZJets");
    h_J1M_tagSB_ZJets = (TH1F*)f->Get("J1pt_M_tagSB_ZJets"); h_J2M_tagSB_ZJets = (TH1F*)f->Get("J2pt_M_tagSB_ZJets");


    h_doubleBSBLow_ZJets = (TH1F*)f->Get("J1_DoubleB_doubleBSBLow_ZJets"); h_doubleBSR_ZJets = (TH1F*)f->Get("J1_DoubleB_doubleBSR_ZJets");
    h_doubleBSBHigh_ZJets = (TH1F*)f->Get("J1_DoubleB_doubleBSBHigh_ZJets");
    h_doubleBSBLow_QCD = (TH1F*)f->Get("J1_DoubleB_doubleBSBLow_QCD"); h_doubleBSR_QCD = (TH1F*)f->Get("J1_DoubleB_doubleBSR_QCD");
    h_doubleBSBHigh_QCD = (TH1F*)f->Get("J1_DoubleB_doubleBSBHigh_QCD");
  }

  h_doubleBSBLow_WJets = (TH1F*)f->Get("J1_DoubleB_doubleBSBLow_WJets"); h_doubleBSR_WJets = (TH1F*)f->Get("J1_DoubleB_doubleBSR_WJets");
  h_doubleBSBHigh_WJets = (TH1F*)f->Get("J1_DoubleB_doubleBSBHigh_WJets");
  h_doubleBSBLow_TT = (TH1F*)f->Get("J1_DoubleB_doubleBSBLow_TT"); h_doubleBSR_TT = (TH1F*)f->Get("J1_DoubleB_doubleBSR_TT");
  h_doubleBSBHigh_TT = (TH1F*)f->Get("J1_DoubleB_doubleBSBHigh_TT");

  vector<TH1F*> histos_ABCD_sum = {h_A_sum, h_B_sum, h_C_sum, h_D_sum};
  vector<TH1F*> histos_ABCD_QCD = {h_A_QCD, h_B_QCD, h_C_QCD, h_D_QCD};
  vector<TH1F*> histos_ABCD_TT = {h_A_TT, h_B_TT, h_C_TT, h_D_TT};
  vector<TH1F*> histos_ABCD_WJets = {h_A_WJets, h_B_WJets, h_C_WJets, h_D_WJets};
  vector<TH1F*> histos_ABCD_ZJets = {h_A_ZJets, h_B_ZJets, h_C_ZJets, h_D_ZJets};

  vector<TH1F*> histos_A1B1CD_sum = {h_A1_sum, h_B1_sum, h_C_sum, h_D_sum};
  vector<TH1F*> histos_A1B1CD_QCD = {h_A1_QCD, h_B1_QCD, h_C_QCD, h_D_QCD};
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


  vector<TH1F*> histos_doubleB_ZJets = {h_doubleBSBLow_ZJets,h_doubleBSR_ZJets,h_doubleBSBHigh_ZJets};
  vector<TH1F*> histos_doubleB_WJets = {h_doubleBSBLow_WJets,h_doubleBSR_WJets,h_doubleBSBHigh_WJets};
  vector<TH1F*> histos_doubleB_TT = {h_doubleBSBLow_TT,h_doubleBSR_TT,h_doubleBSBHigh_TT};
  vector<TH1F*> histos_doubleB_QCD = {h_doubleBSBLow_QCD,h_doubleBSR_QCD,h_doubleBSBHigh_QCD};


  if (runABCDPlots){
    makeABCDPlot(histos_ABCD_sum, "BkgSum", "Double");
    makeABCDPlot(histos_ABCD_QCD, "QCD", "Double");
    makeABCDPlot(histos_ABCD_TT, "TT", "Double");
    makeABCDPlot(histos_ABCD_WJets, "WJets", "Double");
    makeABCDPlot(histos_ABCD_ZJets, "ZJets", "Double");

    makeABCDPlot(histos_A1B1CD_sum, "BkgSum", "Single");
    makeABCDPlot(histos_A1B1CD_QCD, "QCD", "Single");
    makeABCDPlot(histos_A1B1CD_TT, "TT", "Single");
    makeABCDPlot(histos_A1B1CD_WJets, "WJets", "Single");
    makeABCDPlot(histos_A1B1CD_ZJets, "ZJets", "Single");
  }

  if (runRPFPlots){
    myfile.open("Rpf_deviations.txt");
    //J1 Double tag region
    makeRpfPlot(histos_Rpf_J1_sum, "BkgSum", "J1", "Double");
    makeRpfPlot(histos_Rpf_J1_QCD, "QCD", "J1", "Double");
    makeRpfPlot(histos_Rpf_J1_TT, "TT", "J1", "Double");
    makeRpfPlot(histos_Rpf_J1_WJets, "WJets", "J1", "Double");
    makeRpfPlot(histos_Rpf_J1_ZJets, "ZJets", "J1", "Double");

    //J1 Single tag region
    makeRpfPlot(histos_Rpfsingle_J1_sum, "BkgSum", "J1", "Single");
    makeRpfPlot(histos_Rpfsingle_J1_QCD, "QCD", "J1", "Single");
    makeRpfPlot(histos_Rpfsingle_J1_TT, "TT", "J1", "Single");
    makeRpfPlot(histos_Rpfsingle_J1_WJets, "WJets", "J1", "Single");
    makeRpfPlot(histos_Rpfsingle_J1_ZJets, "ZJets", "J1", "Single");

    //J2 double tag region
    makeRpfPlot(histos_Rpf_J2_sum, "BkgSum", "J2", "Double");
    makeRpfPlot(histos_Rpf_J2_QCD, "QCD", "J2", "Double");
    makeRpfPlot(histos_Rpf_J2_TT, "TT", "J2", "Double");
    makeRpfPlot(histos_Rpf_J2_WJets, "WJets", "J2", "Double");
    makeRpfPlot(histos_Rpf_J2_ZJets, "ZJets", "J2", "Double");

    //J2 single tag region
    makeRpfPlot(histos_Rpfsingle_J2_sum, "BkgSum", "J2", "Single");
    makeRpfPlot(histos_Rpfsingle_J2_QCD, "QCD", "J2", "Single");
    makeRpfPlot(histos_Rpfsingle_J2_TT, "TT", "J2", "Single");
    makeRpfPlot(histos_Rpfsingle_J2_WJets, "WJets", "J2", "Single");
    makeRpfPlot(histos_Rpfsingle_J2_ZJets, "ZJets", "J2", "Single");
    myfile.close();
  }


  if (runCompareDoubleB){
    makeDoubleBPlots(histos_doubleB_WJets, "WJets");
    makeDoubleBPlots(histos_doubleB_TT, "TT");
    if (isSignal){
      makeDoubleBPlots(histos_doubleB_ZJets, "ZJets");
      // makeDoubleBPlots(histos_doubleB_QCD, "QCD");
      // makeDoubleBStackPlot(histos_doubleB_QCD,histos_doubleB_TT,histos_doubleB_WJets,histos_doubleB_ZJets);

    }
  }
}

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

  TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo_pred, histo_A, "pois");

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
  // pad1->SetLogy(logy)
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
  if (bkgType == "WJets" && tagType == "Double") histo_pred->SetMaximum(0.5);
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
  graph->SetTitle(";MET [GeV];pred/sim");
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
  graph->GetXaxis()->SetRangeUser(300,1000);
  graph->GetYaxis()->SetRangeUser(0.2,1.8);
  if (bkgType == "TT") graph->GetYaxis()->SetRangeUser(0.2,4.5);
  // graph->Draw("APE");

  can_h->Modified();
  can_h->Update();



  TLine *line = new TLine(300,1.0,1000,1.0);
  line->SetLineColor(kRed);
  line->SetLineStyle(2);
  line->Draw("same");


  TString savename = "ABCDPlots/ABCD"+tagType+"_"+bkgType+"_V17";
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

  can_h->SaveAs("RPFPlots/Rpf_"+thisType+"_V17.pdf");
  num->SaveAs("RPFPlots/Num_"+thisType+"_V17.root");
  denom->SaveAs("RPFPlots/Denom_"+thisType+"_V17.root");
  graph->SaveAs("RPFPlots/Ratio_"+thisType+"_V17.root");

}

void makeDoubleBPlots(vector<TH1F*> dem_histos, TString bkgType = ""){
  TH1F *h_SBLow = dem_histos[0]; TH1F *h_SR = dem_histos[1]; TH1F *h_SBHigh = dem_histos[2];

  h_SBLow->Rebin(2); h_SR->Rebin(2); h_SBHigh->Rebin(2);
  double W = 800;
  double H = 600;
  double T = 0.08*H;
  double B = 0.12*H;
  double L = 0.15*W;
  double R = 0.02*W;
  TString graphName = "DoubleB_"+bkgType;

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

  h_SBLow->GetXaxis()->SetTitle("Double-b discriminator");
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
  can_h->SaveAs("SupportingPlots/DoubleBCompare_"+bkgType+".pdf");

  h_SBLow->ClearUnderflowAndOverflow();

  double failIntegralSBLow = h_SBLow->Integral(0,16)/h_SBLow->Integral()*100;
  double failIntegralSBHigh = h_SBHigh->Integral(0,16)/h_SBHigh->Integral()*100;
  double failIntegralSR = h_SR->Integral(0,16)/h_SR->Integral()*100;

  cout<<bkgType<<": ";
  cout<<"Percent failing, low mass: "<<failIntegralSBLow<<endl;
  cout<<"Percent failing, SR: "<<failIntegralSR<<endl;
  cout<<"Percent failing, high mass: "<<failIntegralSBHigh<<endl<<endl;

}

void makeDoubleBStackPlot(vector<TH1F*> h_QCD,vector<TH1F*> h_TT,vector<TH1F*> h_WJets,vector<TH1F*> h_ZJets){
  TH1D *h_QCD_stack = (TH1D*)h_QCD[0]->Clone("QCDAll");
  h_QCD_stack->Add(h_QCD[1]); h_QCD_stack->Add(h_QCD[2]);
  h_QCD_stack->SetFillColor(kGray);
  h_QCD_stack->SetMarkerStyle(21); h_QCD_stack->SetMarkerColor(kGray);

  TH1D *h_TT_stack = (TH1D*)h_TT[0]->Clone("TTAll");
  h_TT_stack->Add(h_TT[1]); h_TT_stack->Add(h_TT[2]);
  h_TT_stack->SetFillColor(kCyan);
  h_TT_stack->SetMarkerStyle(21); h_TT_stack->SetMarkerColor(kCyan);

  TH1D *h_WJets_stack = (TH1D*)h_WJets[0]->Clone("WJetsAll");
  h_WJets_stack->Add(h_WJets[1]); h_WJets_stack->Add(h_WJets[2]);
  h_WJets_stack->SetFillColor(kBlue);
  h_WJets_stack->SetMarkerStyle(21); h_WJets_stack->SetMarkerColor(kBlue);

  TH1D *h_ZJets_stack = (TH1D*)h_ZJets[0]->Clone("WZetsAll");
  h_ZJets_stack->Add(h_ZJets[1]); h_ZJets_stack->Add(h_ZJets[2]);
  h_ZJets_stack->SetFillColor(kGreen+2);
  h_ZJets_stack->SetMarkerStyle(21); h_ZJets_stack->SetMarkerColor(kGreen+2);

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
  TString graphName = "DoubleB_Stack";

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

  doubleBStack->Draw("hist");
  doubleBStack->GetXaxis()->SetTitle("double-B discriminator"); doubleBStack->GetYaxis()->SetTitle("Events");


  TLegend* legend = new TLegend(0.35,0.85,0.9,0.9) ;
  legend->AddEntry(h_QCD_stack, "QCD", "f");
  legend->AddEntry(h_TT_stack, "TT", "f");
  legend->AddEntry(h_WJets_stack, "WJets", "f");
  legend->AddEntry(h_ZJets_stack, "ZJets", "f");
  legend->SetNColumns(4);
  legend->SetBorderSize(0);

  legend->Draw("same");
  can_h->Update();
  can_h->Modified();
  can_h->SaveAs("SupportingPlots/DoubleBCompare_Stack.pdf");


}
