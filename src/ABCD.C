#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TF1.h>


void makeABCDPlot(vector<TH1F*> dem_histos, TString type, TString tagType, TString year);
void makeRpfPlot(vector<TH1F*> dem_histos, TString bkgType, TString jetType, TString tagType, TString year);
void makeStackPlot(vector<TH1F*> h_QCD,vector<TH1F*> h_TT,vector<TH1F*> h_WJets,vector<TH1F*> h_ZJets, vector<TH1F*> h_T5HH1300,vector<TH1F*> h_T5HH1700,vector<TH1F*> h_T5HH2100,TString which);
void makeSingleLeptStackPlot(vector<TH1F*> h_TT,vector<TH1F*> h_WJets,TString which);
void QuickROC(TH1F* signalHist, TH1F* bkgHist, TString which);
void makeMETStack(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets, TH1F* h_T5HH1300,TH1F* h_T5HH1700,TH1F* h_T5HH2100,TString which);
void makeClosureStackPlot(vector<TH1F*> QCD_histos, vector<TH1F*> TT_histos, vector<TH1F*> WJets_histos,vector<TH1F*> ZJets_histos, vector<TH1F*> sum_histos);
void makeFullBkgClosure(vector<TH1F*> dem_histos, TString bkgType, TString tagType, TString year);
void makeMETNorm(vector<TH1F*> dem_histos, TString tagType);


ofstream myfile;
bool getDeviation = false;
bool runABCDPlots = false;
bool runRPFPlots = true;
bool runStacks = false;
bool runROC = false;
bool runMETNorm = false;


string whichRegion = "signal";
// string whichRegion = "singleLept";
// string whichRegion = "photon";

// TString year = "2016";
TString year = "Run2";


TFile * f = TFile::Open("ALPHABET_V17_MET300_2BoostedH.root");
// TFile * f = TFile::Open("ALPHABETBoost_V17_MET250.root");


TFile * fSignal = TFile::Open("ALPHABETBoost_V17_signalOnly.root");
TFile * fPhoton;// = TFile::Open("ALPHABETBoostMC2016_V12like_photon.root");


// TFile * fSingleLept;// = TFile::Open("ALPHABETBoost_2017_1l.root");
TFile * fSingleLept = TFile::Open("ALPHABET_V17_1l_MET300_2BoostedH.root");
// TFile * fSingleLept = TFile::Open("ALPHABETdata_V17_1l_2BoostedH.root");


TFile* fout = new TFile("ABCDPlots_V17.root","recreate");
// TFile* fout = new TFile("ABCDPlots_V17_MET250.root","recreate");
// TFile* fout = new TFile("ABCDPlots_V17data_1l.root","recreate");
// TFile* fout = new TFile("ABCDPlots_V17_1l.root","recreate");
// TFile* fout = new TFile("ABCDPlots_V17_photon.root","recreate");

TDirectory *cdABCD  = fout->mkdir("ABCD");
TDirectory *cdClose = fout->mkdir("Closure");
TDirectory *cdAssum = fout->mkdir("Assumption");
TDirectory *cdRPFFinal  = fout->mkdir("RPF");
TDirectory *cdRPF  = fout->mkdir("RPF_Support");
TDirectory *cdOther  = fout->mkdir("OtherPlots");

void runABCD() {
  TH1F * h_2HLeadSR_sum; TH1F * h_0HLeadSR_sum;
  TH1F * h_2HLeadSB_sum; TH1F * h_0HLeadSB_sum;

  TH1F * h_A_data; TH1F * h_B_data; TH1F * h_A1_data; TH1F * h_B1_data; TH1F *h_C_data; TH1F *h_D_data;

  TH1F * h_A_sum; TH1F * h_B_sum; TH1F * h_A1_sum; TH1F * h_B1_sum; TH1F *h_C_sum; TH1F *h_D_sum;
  TH1F * h_A_QCD; TH1F * h_B_QCD; TH1F * h_A1_QCD; TH1F * h_B1_QCD; TH1F * h_C_QCD; TH1F * h_D_QCD;
  TH1F * h_A_GJets; TH1F * h_B_GJets; TH1F * h_A1_GJets; TH1F * h_B1_GJets; TH1F * h_C_GJets; TH1F * h_D_GJets;
  TH1F * h_A_TT; TH1F * h_B_TT; TH1F * h_A1_TT; TH1F * h_B1_TT; TH1F * h_C_TT; TH1F * h_D_TT;
  TH1F * h_A_TT_Di; TH1F * h_B_TT_Di; TH1F * h_A1_TT_Di; TH1F * h_B1_TT_Di; TH1F * h_C_TT_Di; TH1F * h_D_TT_Di;
  TH1F * h_A_TT_SL; TH1F * h_B_TT_SL; TH1F * h_A1_TT_SL; TH1F * h_B1_TT_SL; TH1F * h_C_TT_SL; TH1F * h_D_TT_SL;
  TH1F * h_A_WJets; TH1F * h_B_WJets; TH1F * h_A1_WJets; TH1F * h_B1_WJets; TH1F * h_C_WJets; TH1F * h_D_WJets;
  TH1F * h_A_ZJets; TH1F * h_B_ZJets; TH1F * h_A1_ZJets; TH1F * h_B1_ZJets; TH1F * h_C_ZJets; TH1F * h_D_ZJets;


  TH1F * h_A_TChiHH200; TH1F * h_B_TChiHH200; TH1F * h_A1_TChiHH200; TH1F * h_B1_TChiHH200; TH1F * h_C_TChiHH200; TH1F * h_D_TChiHH200;
  TH1F * h_A_TChiHH400; TH1F * h_B_TChiHH400; TH1F * h_A1_TChiHH400; TH1F * h_B1_TChiHH400; TH1F * h_C_TChiHH400; TH1F * h_D_TChiHH400;
  TH1F * h_A_TChiHH700; TH1F * h_B_TChiHH700; TH1F * h_A1_TChiHH700; TH1F * h_B1_TChiHH700; TH1F * h_C_TChiHH700; TH1F * h_D_TChiHH700;
  TH1F * h_A_TChiHH1000; TH1F * h_B_TChiHH1000; TH1F * h_A1_TChiHH1000; TH1F * h_B1_TChiHH1000; TH1F * h_C_TChiHH1000; TH1F * h_D_TChiHH1000;
  TH1F * h_A_TChiHH1300; TH1F * h_B_TChiHH1300; TH1F * h_A1_TChiHH1300; TH1F * h_B1_TChiHH1300; TH1F * h_C_TChiHH1300; TH1F * h_D_TChiHH1300;
  TH1F * h_A_TChiHH1500; TH1F * h_B_TChiHH1500; TH1F * h_A1_TChiHH1500; TH1F * h_B1_TChiHH1500; TH1F * h_C_TChiHH1500; TH1F * h_D_TChiHH1500;

  TH1F * h_J1M_doubletagSR_sum; TH1F * h_J2M_doubletagSR_sum; TH1F * h_J1M_doubletagSB_sum; TH1F * h_J2M_doubletagSB_sum;
  TH1F * h_J1M_antitagSR_sum; TH1F * h_J2M_antitagSR_sum; TH1F * h_J1M_antitagSB_sum; TH1F * h_J2M_antitagSB_sum;
  TH1F * h_J1M_tagSR_sum; TH1F * h_J2M_tagSR_sum; TH1F * h_J1M_tagSB_sum; TH1F * h_J2M_tagSB_sum;
  TH1F * h_J2M_mjBins_doubletagSR_sum; TH1F * h_J2M_mjBins_doubletagSB_sum;
  TH1F * h_J2M_mjBins_antitagSR_sum; TH1F * h_J2M_mjBins_antitagSB_sum;
  TH1F * h_J1M_doubletagSR_QCD; TH1F * h_J2M_doubletagSR_QCD; TH1F * h_J1M_doubletagSB_QCD; TH1F * h_J2M_doubletagSB_QCD;
  TH1F * h_J1M_antitagSR_QCD; TH1F * h_J2M_antitagSR_QCD; TH1F * h_J1M_antitagSB_QCD; TH1F * h_J2M_antitagSB_QCD;
  TH1F * h_J1M_tagSR_QCD; TH1F * h_J2M_tagSR_QCD; TH1F * h_J1M_tagSB_QCD; TH1F * h_J2M_tagSB_QCD;

  TH1F * h_J2M_mjBins_doubletagSR_QCD; TH1F * h_J2M_mjBins_doubletagSB_QCD;
  TH1F * h_J2M_mjBins_antitagSR_QCD; TH1F * h_J2M_mjBins_antitagSB_QCD;
  TH1F * h_J2M_mjBins_doubletagSR_GJets; TH1F * h_J2M_mjBins_doubletagSB_GJets;
  TH1F * h_J2M_mjBins_antitagSR_GJets; TH1F * h_J2M_mjBins_antitagSB_GJets;

  TH1F * h_J1M_doubletagSR_GJets; TH1F * h_J2M_doubletagSR_GJets; TH1F * h_J1M_doubletagSB_GJets; TH1F * h_J2M_doubletagSB_GJets;
  TH1F * h_J1M_antitagSR_GJets; TH1F * h_J2M_antitagSR_GJets; TH1F * h_J1M_antitagSB_GJets; TH1F * h_J2M_antitagSB_GJets;
  TH1F * h_J1M_tagSR_GJets; TH1F * h_J2M_tagSR_GJets; TH1F * h_J1M_tagSB_GJets; TH1F * h_J2M_tagSB_GJets;

  TH1F * h_J1M_doubletagSR_TT; TH1F * h_J2M_doubletagSR_TT; TH1F * h_J1M_doubletagSB_TT; TH1F * h_J2M_doubletagSB_TT;
  TH1F * h_J1M_antitagSR_TT; TH1F * h_J2M_antitagSR_TT; TH1F * h_J1M_antitagSB_TT; TH1F * h_J2M_antitagSB_TT;
  TH1F * h_J1M_tagSR_TT; TH1F * h_J2M_tagSR_TT; TH1F * h_J1M_tagSB_TT; TH1F * h_J2M_tagSB_TT;
  TH1F * h_J2M_mjBins_doubletagSR_TT; TH1F * h_J2M_mjBins_doubletagSB_TT;
  TH1F * h_J2M_mjBins_antitagSR_TT; TH1F * h_J2M_mjBins_antitagSB_TT;

  TH1F * h_J1M_doubletagSR_TT_SL; TH1F * h_J2M_doubletagSR_TT_SL; TH1F * h_J1M_doubletagSB_TT_SL; TH1F * h_J2M_doubletagSB_TT_SL;
  TH1F * h_J1M_antitagSR_TT_SL; TH1F * h_J2M_antitagSR_TT_SL; TH1F * h_J1M_antitagSB_TT_SL; TH1F * h_J2M_antitagSB_TT_SL;
  TH1F * h_J1M_tagSR_TT_SL; TH1F * h_J2M_tagSR_TT_SL; TH1F * h_J1M_tagSB_TT_SL; TH1F * h_J2M_tagSB_TT_SL;
  TH1F * h_J2M_mjBins_doubletagSR_TT_SL; TH1F * h_J2M_mjBins_doubletagSB_TT_SL;
  TH1F * h_J2M_mjBins_antitagSR_TT_SL; TH1F * h_J2M_mjBins_antitagSB_TT_SL;

  TH1F * h_J1M_doubletagSR_TT_Di; TH1F * h_J2M_doubletagSR_TT_Di; TH1F * h_J1M_doubletagSB_TT_Di; TH1F * h_J2M_doubletagSB_TT_Di;
  TH1F * h_J1M_antitagSR_TT_Di; TH1F * h_J2M_antitagSR_TT_Di; TH1F * h_J1M_antitagSB_TT_Di; TH1F * h_J2M_antitagSB_TT_Di;
  TH1F * h_J1M_tagSR_TT_Di; TH1F * h_J2M_tagSR_TT_Di; TH1F * h_J1M_tagSB_TT_Di; TH1F * h_J2M_tagSB_TT_Di;
  TH1F * h_J2M_mjBins_doubletagSR_TT_Di; TH1F * h_J2M_mjBins_doubletagSB_TT_Di;
  TH1F * h_J2M_mjBins_antitagSR_TT_Di; TH1F * h_J2M_mjBins_antitagSB_TT_Di;

  TH1F * h_J1M_doubletagSR_WJets; TH1F * h_J2M_doubletagSR_WJets; TH1F * h_J1M_doubletagSB_WJets; TH1F * h_J2M_doubletagSB_WJets;
  TH1F * h_J1M_antitagSR_WJets; TH1F * h_J2M_antitagSR_WJets; TH1F * h_J1M_antitagSB_WJets; TH1F * h_J2M_antitagSB_WJets;
  TH1F * h_J1M_tagSR_WJets; TH1F * h_J2M_tagSR_WJets; TH1F * h_J1M_tagSB_WJets; TH1F * h_J2M_tagSB_WJets;
  TH1F * h_J2M_mjBins_doubletagSR_WJets; TH1F * h_J2M_mjBins_doubletagSB_WJets;
  TH1F * h_J2M_mjBins_antitagSR_WJets; TH1F * h_J2M_mjBins_antitagSB_WJets;

  TH1F * h_J1M_doubletagSR_ZJets; TH1F * h_J2M_doubletagSR_ZJets; TH1F * h_J1M_doubletagSB_ZJets; TH1F * h_J2M_doubletagSB_ZJets;
  TH1F * h_J1M_antitagSR_ZJets; TH1F * h_J2M_antitagSR_ZJets; TH1F * h_J1M_antitagSB_ZJets; TH1F * h_J2M_antitagSB_ZJets;
  TH1F * h_J1M_tagSR_ZJets; TH1F * h_J2M_tagSR_ZJets; TH1F * h_J1M_tagSB_ZJets; TH1F * h_J2M_tagSB_ZJets;
  TH1F * h_J2M_mjBins_doubletagSR_ZJets; TH1F * h_J2M_mjBins_doubletagSB_ZJets;
  TH1F * h_J2M_mjBins_antitagSR_ZJets; TH1F * h_J2M_mjBins_antitagSB_ZJets;

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

  TH1F * h_J1Pt_doubletagSR_TT_Di; TH1F * h_J2Pt_doubletagSR_TT_Di; TH1F * h_J1Pt_doubletagSB_TT_Di; TH1F * h_J2Pt_doubletagSB_TT_Di;
  TH1F * h_J1Pt_antitagSR_TT_Di; TH1F * h_J2Pt_antitagSR_TT_Di; TH1F * h_J1Pt_antitagSB_TT_Di; TH1F * h_J2Pt_antitagSB_TT_Di;
  TH1F * h_J1Pt_tagSR_TT_Di; TH1F * h_J2Pt_tagSR_TT_Di; TH1F * h_J1Pt_tagSB_TT_Di; TH1F * h_J2Pt_tagSB_TT_Di;

  TH1F * h_J1Pt_doubletagSR_TT_SL; TH1F * h_J2Pt_doubletagSR_TT_SL; TH1F * h_J1Pt_doubletagSB_TT_SL; TH1F * h_J2Pt_doubletagSB_TT_SL;
  TH1F * h_J1Pt_antitagSR_TT_SL; TH1F * h_J2Pt_antitagSR_TT_SL; TH1F * h_J1Pt_antitagSB_TT_SL; TH1F * h_J2Pt_antitagSB_TT_SL;
  TH1F * h_J1Pt_tagSR_TT_SL; TH1F * h_J2Pt_tagSR_TT_SL; TH1F * h_J1Pt_tagSB_TT_SL; TH1F * h_J2Pt_tagSB_TT_SL;

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

  TH1F * h_notagSR_sum; TH1F * h_notagSB_sum; TH1F * h_notagSB2_sum;
  TH1F * h_notagSR_ZJets; TH1F * h_notagSB_ZJets;
  TH1F * h_notagSR_WJets; TH1F * h_notagSB_WJets;
  TH1F * h_notagSR_TT; TH1F * h_notagSB_TT;
  TH1F * h_notagSR_TT_Di; TH1F * h_notagSB_TT_Di;
  TH1F * h_notagSR_TT_SL; TH1F * h_notagSB_TT_SL;
  TH1F * h_notagSR_QCD; TH1F * h_notagSB_QCD;

  TH1F * h_0HSBOpt1; TH1F * h_0HSROpt1;
  TH1F * h_0HSBOpt2; TH1F * h_0HSROpt2;
  TH1F * h_0HSBOpt1_TT; TH1F * h_0HSBOpt1_WJets;
  TH1F * h_0HSROpt1_TT; TH1F * h_0HSROpt1_WJets;

  if (whichRegion=="signal") {
    h_2HLeadSR_sum = (TH1F*)f->Get("J2pt_M_doubletagSRLead_sum"); h_0HLeadSR_sum = (TH1F*)f->Get("J2pt_M_antitagSRLead_sum");
    h_2HLeadSB_sum = (TH1F*)f->Get("J1pt_M_doubletagSBLead_sum"); h_0HLeadSB_sum = (TH1F*)f->Get("J1pt_M_antitagSBLead_sum");

    h_A_sum = (TH1F*)f->Get("MET_doubletagSR_sum"); h_B_sum = (TH1F*)f->Get("MET_doubletagSB_sum");
    h_A1_sum = (TH1F*)f->Get("MET_tagSR_sum"); h_B1_sum = (TH1F*)f->Get("MET_tagSB_sum");
    h_C_sum = (TH1F*)f->Get("MET_antitagSR_sum"); h_D_sum = (TH1F*)f->Get("MET_antitagSB_sum");

    h_A_TT = (TH1F*)f->Get("MET_doubletagSR_TT"); h_B_TT = (TH1F*)f->Get("MET_doubletagSB_TT");
    h_A1_TT = (TH1F*)f->Get("MET_tagSR_TT"); h_B1_TT = (TH1F*)f->Get("MET_tagSB_TT");
    h_C_TT = (TH1F*)f->Get("MET_antitagSR_TT"); h_D_TT = (TH1F*)f->Get("MET_antitagSB_TT");

    // h_A_TT_Di = (TH1F*)f->Get("MET_doubletagSR_TT_Di"); h_B_TT_Di = (TH1F*)f->Get("MET_doubletagSB_TT_Di");
    // h_A1_TT_Di = (TH1F*)f->Get("MET_tagSR_TT_Di"); h_B1_TT_Di = (TH1F*)f->Get("MET_tagSB_TT_Di");
    // h_C_TT_Di = (TH1F*)f->Get("MET_antitagSR_TT_Di"); h_D_TT_Di = (TH1F*)f->Get("MET_antitagSB_TT_Di");
    // h_A_TT_SL = (TH1F*)f->Get("MET_doubletagSR_TT_SL"); h_B_TT_SL = (TH1F*)f->Get("MET_doubletagSB_TT_SL");
    // h_A1_TT_SL = (TH1F*)f->Get("MET_tagSR_TT_SL"); h_B1_TT_SL = (TH1F*)f->Get("MET_tagSB_TT_SL");
    // h_C_TT_SL = (TH1F*)f->Get("MET_antitagSR_TT_SL"); h_D_TT_SL = (TH1F*)f->Get("MET_antitagSB_TT_SL");

    h_A_WJets = (TH1F*)f->Get("MET_doubletagSR_WJets"); h_B_WJets = (TH1F*)f->Get("MET_doubletagSB_WJets");
    h_A1_WJets = (TH1F*)f->Get("MET_tagSR_WJets"); h_B1_WJets = (TH1F*)f->Get("MET_tagSB_WJets");
    h_C_WJets = (TH1F*)f->Get("MET_antitagSR_WJets");   h_D_WJets = (TH1F*)f->Get("MET_antitagSB_WJets");
    h_A_ZJets = (TH1F*)f->Get("MET_doubletagSR_ZJets"); h_B_ZJets = (TH1F*)f->Get("MET_doubletagSB_ZJets");
    h_A1_ZJets = (TH1F*)f->Get("MET_tagSR_ZJets"); h_B1_ZJets = (TH1F*)f->Get("MET_tagSB_ZJets");
    h_C_ZJets = (TH1F*)f->Get("MET_antitagSR_ZJets");   h_D_ZJets = (TH1F*)f->Get("MET_antitagSB_ZJets");
    h_A_QCD = (TH1F*)f->Get("MET_doubletagSR_QCD"); h_B_QCD = (TH1F*)f->Get("MET_doubletagSB_QCD");
    h_A1_QCD = (TH1F*)f->Get("MET_tagSR_QCD"); h_B1_QCD = (TH1F*)f->Get("MET_tagSB_QCD");
    h_C_QCD = (TH1F*)f->Get("MET_antitagSR_QCD"); h_D_QCD = (TH1F*)f->Get("MET_antitagSB_QCD");

    //new regions here
    h_notagSR_sum = (TH1F*)f->Get("MET_notagSR_sum"); h_notagSB_sum = (TH1F*)f->Get("MET_notagSB_sum");
    h_notagSR_ZJets = (TH1F*)f->Get("MET_notagSR_ZJets"); h_notagSB_ZJets = (TH1F*)f->Get("MET_notagSB_ZJets");
    h_notagSR_WJets = (TH1F*)f->Get("MET_notagSR_WJets"); h_notagSB_WJets = (TH1F*)f->Get("MET_notagSB_WJets");
    h_notagSR_TT = (TH1F*)f->Get("MET_notagSR_TT"); h_notagSB_TT = (TH1F*)f->Get("MET_notagSB_TT");
    // h_notagSR_TT_Di = (TH1F*)f->Get("MET_notagSR_TT_Di"); h_notagSB_TT_Di = (TH1F*)f->Get("MET_notagSB_TT_Di");
    // h_notagSR_TT_SL = (TH1F*)f->Get("MET_notagSR_TT_SL"); h_notagSB_TT_SL = (TH1F*)f->Get("MET_notagSB_TT_SL");
    h_notagSR_QCD = (TH1F*)f->Get("MET_notagSR_QCD"); h_notagSB_QCD = (TH1F*)f->Get("MET_notagSB_QCD");

    h_A_TChiHH200 = (TH1F*)fSignal->Get("MET_doubletagSR_TChiHH200"); h_B_TChiHH200 = (TH1F*)fSignal->Get("MET_doubletagSB_TChiHH200");
    h_A1_TChiHH200 = (TH1F*)fSignal->Get("MET_tagSR_TChiHH200"); h_B1_TChiHH200 = (TH1F*)fSignal->Get("MET_tagSB_TChiHH200");
    h_C_TChiHH200 = (TH1F*)fSignal->Get("MET_antitagSR_TChiHH200"); h_D_TChiHH200 = (TH1F*)fSignal->Get("MET_antitagSB_TChiHH200");
    h_A_TChiHH400 = (TH1F*)fSignal->Get("MET_doubletagSR_TChiHH400"); h_B_TChiHH400 = (TH1F*)fSignal->Get("MET_doubletagSB_TChiHH400");
    h_A1_TChiHH400 = (TH1F*)fSignal->Get("MET_tagSR_TChiHH400"); h_B1_TChiHH400 = (TH1F*)fSignal->Get("MET_tagSB_TChiHH400");
    h_C_TChiHH400 = (TH1F*)fSignal->Get("MET_antitagSR_TChiHH400"); h_D_TChiHH400 = (TH1F*)fSignal->Get("MET_antitagSB_TChiHH400");
    h_A_TChiHH700 = (TH1F*)fSignal->Get("MET_doubletagSR_TChiHH700"); h_B_TChiHH700 = (TH1F*)fSignal->Get("MET_doubletagSB_TChiHH700");
    h_A1_TChiHH700 = (TH1F*)fSignal->Get("MET_tagSR_TChiHH700"); h_B1_TChiHH700 = (TH1F*)fSignal->Get("MET_tagSB_TChiHH700");
    h_C_TChiHH700 = (TH1F*)fSignal->Get("MET_antitagSR_TChiHH700"); h_D_TChiHH700 = (TH1F*)fSignal->Get("MET_antitagSB_TChiHH700");
    h_A_TChiHH1000 = (TH1F*)fSignal->Get("MET_doubletagSR_TChiHH1000"); h_B_TChiHH1000 = (TH1F*)fSignal->Get("MET_doubletagSB_TChiHH1000");
    h_A1_TChiHH1000 = (TH1F*)fSignal->Get("MET_tagSR_TChiHH1000"); h_B1_TChiHH1000 = (TH1F*)fSignal->Get("MET_tagSB_TChiHH1000");
    h_C_TChiHH1000 = (TH1F*)fSignal->Get("MET_antitagSR_TChiHH1000"); h_D_TChiHH1000 = (TH1F*)fSignal->Get("MET_antitagSB_TChiHH1000");
    h_A_TChiHH1300 = (TH1F*)fSignal->Get("MET_doubletagSR_TChiHH1300"); h_B_TChiHH1300 = (TH1F*)fSignal->Get("MET_doubletagSB_TChiHH1300");
    h_A1_TChiHH1300 = (TH1F*)fSignal->Get("MET_tagSR_TChiHH1300"); h_B1_TChiHH1300 = (TH1F*)fSignal->Get("MET_tagSB_TChiHH1300");
    h_C_TChiHH1300 = (TH1F*)fSignal->Get("MET_antitagSR_TChiHH1300"); h_D_TChiHH1300 = (TH1F*)fSignal->Get("MET_antitagSB_TChiHH1300");
    h_A_TChiHH1500 = (TH1F*)fSignal->Get("MET_doubletagSR_TChiHH1500"); h_B_TChiHH1500 = (TH1F*)fSignal->Get("MET_doubletagSB_TChiHH1500");
    h_A1_TChiHH1500 = (TH1F*)fSignal->Get("MET_tagSR_TChiHH1500"); h_B1_TChiHH1500 = (TH1F*)fSignal->Get("MET_tagSB_TChiHH1500");
    h_C_TChiHH1500 = (TH1F*)fSignal->Get("MET_antitagSR_TChiHH1500"); h_D_TChiHH1500 = (TH1F*)fSignal->Get("MET_antitagSB_TChiHH1500");

    h_J1M_doubletagSR_sum = (TH1F*)f->Get("J1pt_M_doubletagSR_sum"); h_J2M_doubletagSR_sum = (TH1F*)f->Get("J2pt_M_doubletagSR_sum");
    h_J1M_doubletagSB_sum = (TH1F*)f->Get("J1pt_M_doubletagSB_sum"); h_J2M_doubletagSB_sum = (TH1F*)f->Get("J2pt_M_doubletagSB_sum");
    h_J1M_antitagSR_sum = (TH1F*)f->Get("J1pt_M_antitagSR_sum"); h_J2M_antitagSR_sum = (TH1F*)f->Get("J2pt_M_antitagSR_sum");
    h_J1M_antitagSB_sum = (TH1F*)f->Get("J1pt_M_antitagSB_sum"); h_J2M_antitagSB_sum = (TH1F*)f->Get("J2pt_M_antitagSB_sum");
    h_J1M_tagSR_sum = (TH1F*)f->Get("J1pt_M_tagSR_sum"); h_J2M_tagSR_sum = (TH1F*)f->Get("J2pt_M_tagSR_sum");
    h_J1M_tagSB_sum = (TH1F*)f->Get("J1pt_M_tagSB_sum"); h_J2M_tagSB_sum = (TH1F*)f->Get("J2pt_M_tagSB_sum");

    h_J1M_doubletagSR_ZJets = (TH1F*)f->Get("J1pt_M_doubletagSR_ZJets"); h_J2M_doubletagSR_ZJets = (TH1F*)f->Get("J2pt_M_doubletagSR_ZJets");
    h_J1M_doubletagSB_ZJets = (TH1F*)f->Get("J1pt_M_doubletagSB_ZJets"); h_J2M_doubletagSB_ZJets = (TH1F*)f->Get("J2pt_M_doubletagSB_ZJets");
    h_J1M_antitagSR_ZJets = (TH1F*)f->Get("J1pt_M_antitagSR_ZJets"); h_J2M_antitagSR_ZJets = (TH1F*)f->Get("J2pt_M_antitagSR_ZJets");
    h_J1M_antitagSB_ZJets = (TH1F*)f->Get("J1pt_M_antitagSB_ZJets"); h_J2M_antitagSB_ZJets = (TH1F*)f->Get("J2pt_M_antitagSB_ZJets");
    h_J1M_tagSR_ZJets = (TH1F*)f->Get("J1pt_M_tagSR_ZJets"); h_J2M_tagSR_ZJets = (TH1F*)f->Get("J2pt_M_tagSR_ZJets");
    h_J1M_tagSB_ZJets = (TH1F*)f->Get("J1pt_M_tagSB_ZJets"); h_J2M_tagSB_ZJets = (TH1F*)f->Get("J2pt_M_tagSB_ZJets");

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

    // h_J1M_doubletagSR_TT_SL = (TH1F*)f->Get("J1pt_M_doubletagSR_TT_SL"); h_J2M_doubletagSR_TT_SL = (TH1F*)f->Get("J2pt_M_doubletagSR_TT_SL");
    // h_J1M_doubletagSB_TT_SL = (TH1F*)f->Get("J1pt_M_doubletagSB_TT_SL"); h_J2M_doubletagSB_TT_SL = (TH1F*)f->Get("J2pt_M_doubletagSB_TT_SL");
    // h_J1M_antitagSR_TT_SL = (TH1F*)f->Get("J1pt_M_antitagSR_TT_SL"); h_J2M_antitagSR_TT_SL = (TH1F*)f->Get("J2pt_M_antitagSR_TT_SL");
    // h_J1M_antitagSB_TT_SL = (TH1F*)f->Get("J1pt_M_antitagSB_TT_SL"); h_J2M_antitagSB_TT_SL = (TH1F*)f->Get("J2pt_M_antitagSB_TT_SL");
    // h_J1M_tagSR_TT_SL = (TH1F*)f->Get("J1pt_M_tagSR_TT_SL"); h_J2M_tagSR_TT_SL = (TH1F*)f->Get("J2pt_M_tagSR_TT_SL");
    // h_J1M_tagSB_TT_SL = (TH1F*)f->Get("J1pt_M_tagSB_TT_SL"); h_J2M_tagSB_TT_SL = (TH1F*)f->Get("J2pt_M_tagSB_TT_SL");
    // h_J1M_doubletagSR_TT_Di = (TH1F*)f->Get("J1pt_M_doubletagSR_TT_Di"); h_J2M_doubletagSR_TT_Di = (TH1F*)f->Get("J2pt_M_doubletagSR_TT_Di");
    // h_J1M_doubletagSB_TT_Di = (TH1F*)f->Get("J1pt_M_doubletagSB_TT_Di"); h_J2M_doubletagSB_TT_Di = (TH1F*)f->Get("J2pt_M_doubletagSB_TT_Di");
    // h_J1M_antitagSR_TT_Di = (TH1F*)f->Get("J1pt_M_antitagSR_TT_Di"); h_J2M_antitagSR_TT_Di = (TH1F*)f->Get("J2pt_M_antitagSR_TT_Di");
    // h_J1M_antitagSB_TT_Di = (TH1F*)f->Get("J1pt_M_antitagSB_TT_Di"); h_J2M_antitagSB_TT_Di = (TH1F*)f->Get("J2pt_M_antitagSB_TT_Di");
    // h_J1M_tagSR_TT_Di = (TH1F*)f->Get("J1pt_M_tagSR_TT_Di"); h_J2M_tagSR_TT_Di = (TH1F*)f->Get("J2pt_M_tagSR_TT_Di");
    // h_J1M_tagSB_TT_Di = (TH1F*)f->Get("J1pt_M_tagSB_TT_Di"); h_J2M_tagSB_TT_Di = (TH1F*)f->Get("J2pt_M_tagSB_TT_Di");


    h_J1M_doubletagSR_WJets = (TH1F*)f->Get("J1pt_M_doubletagSR_WJets"); h_J2M_doubletagSR_WJets = (TH1F*)f->Get("J2pt_M_doubletagSR_WJets");
    h_J1M_doubletagSB_WJets = (TH1F*)f->Get("J1pt_M_doubletagSB_WJets"); h_J2M_doubletagSB_WJets = (TH1F*)f->Get("J2pt_M_doubletagSB_WJets");
    h_J1M_antitagSR_WJets = (TH1F*)f->Get("J1pt_M_antitagSR_WJets"); h_J2M_antitagSR_WJets = (TH1F*)f->Get("J2pt_M_antitagSR_WJets");
    h_J1M_antitagSB_WJets = (TH1F*)f->Get("J1pt_M_antitagSB_WJets"); h_J2M_antitagSB_WJets = (TH1F*)f->Get("J2pt_M_antitagSB_WJets");
    h_J1M_tagSR_WJets = (TH1F*)f->Get("J1pt_M_tagSR_WJets"); h_J2M_tagSR_WJets = (TH1F*)f->Get("J2pt_M_tagSR_WJets");
    h_J1M_tagSB_WJets = (TH1F*)f->Get("J1pt_M_tagSB_WJets"); h_J2M_tagSB_WJets = (TH1F*)f->Get("J2pt_M_tagSB_WJets");

    h_J2M_mjBins_doubletagSR_sum = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_sum"); h_J2M_mjBins_doubletagSB_sum = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_sum");
    h_J2M_mjBins_antitagSR_sum = (TH1F*)f->Get("J2_M_jetBins_antitagSR_sum"); h_J2M_mjBins_antitagSB_sum = (TH1F*)f->Get("J2_M_jetBins_antitagSB_sum");
    h_J2M_mjBins_doubletagSR_QCD = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_QCD"); h_J2M_mjBins_doubletagSB_QCD = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_QCD");
    h_J2M_mjBins_antitagSR_QCD = (TH1F*)f->Get("J2_M_jetBins_antitagSR_QCD"); h_J2M_mjBins_antitagSB_QCD = (TH1F*)f->Get("J2_M_jetBins_antitagSB_QCD");

    h_J2M_mjBins_doubletagSR_TT = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_TT"); h_J2M_mjBins_doubletagSB_TT = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_TT");
    h_J2M_mjBins_antitagSR_TT = (TH1F*)f->Get("J2_M_jetBins_antitagSR_TT"); h_J2M_mjBins_antitagSB_TT = (TH1F*)f->Get("J2_M_jetBins_antitagSB_TT");
    // h_J2M_mjBins_doubletagSR_TT_Di = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_TT_Di"); h_J2M_mjBins_doubletagSB_TT_Di = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_TT_Di");
    // h_J2M_mjBins_antitagSR_TT_Di = (TH1F*)f->Get("J2_M_jetBins_antitagSR_TT_Di"); h_J2M_mjBins_antitagSB_TT_Di = (TH1F*)f->Get("J2_M_jetBins_antitagSB_TT_Di");
    // h_J2M_mjBins_doubletagSR_TT_SL = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_TT_SL"); h_J2M_mjBins_doubletagSB_TT_SL = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_TT_SL");
    // h_J2M_mjBins_antitagSR_TT_SL = (TH1F*)f->Get("J2_M_jetBins_antitagSR_TT_SL"); h_J2M_mjBins_antitagSB_TT_SL = (TH1F*)f->Get("J2_M_jetBins_antitagSB_TT_SL");
    h_J2M_mjBins_doubletagSR_WJets = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_WJets"); h_J2M_mjBins_doubletagSB_WJets = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_WJets");
    h_J2M_mjBins_antitagSR_WJets = (TH1F*)f->Get("J2_M_jetBins_antitagSR_WJets"); h_J2M_mjBins_antitagSB_WJets = (TH1F*)f->Get("J2_M_jetBins_antitagSB_WJets");

    h_J2M_mjBins_doubletagSR_ZJets = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_ZJets"); h_J2M_mjBins_doubletagSB_ZJets = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_ZJets");
    h_J2M_mjBins_antitagSR_ZJets = (TH1F*)f->Get("J2_M_jetBins_antitagSR_ZJets"); h_J2M_mjBins_antitagSB_ZJets = (TH1F*)f->Get("J2_M_jetBins_antitagSB_ZJets");

    //for signal
    h_J1M_doubletagSR_T5HH1300 = (TH1F*)fSignal->Get("J1pt_M_doubletagSR_T5qqqqZH1300"); h_J2M_doubletagSR_T5HH1300 = (TH1F*)fSignal->Get("J2pt_M_doubletagSR_T5qqqqZH1300");
    h_J1M_doubletagSB_T5HH1300 = (TH1F*)fSignal->Get("J1pt_M_doubletagSB_T5qqqqZH1300"); h_J2M_doubletagSB_T5HH1300 = (TH1F*)fSignal->Get("J2pt_M_doubletagSB_T5qqqqZH1300");
    h_J1M_antitagSR_T5HH1300 = (TH1F*)fSignal->Get("J1pt_M_antitagSR_T5qqqqZH1300"); h_J2M_antitagSR_T5HH1300 = (TH1F*)fSignal->Get("J2pt_M_antitagSR_T5qqqqZH1300");
    h_J1M_antitagSB_T5HH1300 = (TH1F*)fSignal->Get("J1pt_M_antitagSB_T5qqqqZH1300"); h_J2M_antitagSB_T5HH1300 = (TH1F*)fSignal->Get("J2pt_M_antitagSB_T5qqqqZH1300");
    h_J1M_tagSR_T5HH1300 = (TH1F*)fSignal->Get("J1pt_M_tagSR_T5qqqqZH1300"); h_J2M_tagSR_T5HH1300 = (TH1F*)fSignal->Get("J2pt_M_tagSR_T5qqqqZH1300");
    h_J1M_tagSB_T5HH1300 = (TH1F*)fSignal->Get("J1pt_M_tagSB_T5qqqqZH1300"); h_J2M_tagSB_T5HH1300 = (TH1F*)fSignal->Get("J2pt_M_tagSB_T5qqqqZH1300");
    h_J1M_doubletagSR_T5HH1700 = (TH1F*)fSignal->Get("J1pt_M_doubletagSR_T5qqqqZH1700"); h_J2M_doubletagSR_T5HH1700 = (TH1F*)fSignal->Get("J2pt_M_doubletagSR_T5qqqqZH1700");
    h_J1M_doubletagSB_T5HH1700 = (TH1F*)fSignal->Get("J1pt_M_doubletagSB_T5qqqqZH1700"); h_J2M_doubletagSB_T5HH1700 = (TH1F*)fSignal->Get("J2pt_M_doubletagSB_T5qqqqZH1700");
    h_J1M_antitagSR_T5HH1700 = (TH1F*)fSignal->Get("J1pt_M_antitagSR_T5qqqqZH1700"); h_J2M_antitagSR_T5HH1700 = (TH1F*)fSignal->Get("J2pt_M_antitagSR_T5qqqqZH1700");
    h_J1M_antitagSB_T5HH1700 = (TH1F*)fSignal->Get("J1pt_M_antitagSB_T5qqqqZH1700"); h_J2M_antitagSB_T5HH1700 = (TH1F*)fSignal->Get("J2pt_M_antitagSB_T5qqqqZH1700");
    h_J1M_tagSR_T5HH1700 = (TH1F*)fSignal->Get("J1pt_M_tagSR_T5qqqqZH1700"); h_J2M_tagSR_T5HH1700 = (TH1F*)fSignal->Get("J2pt_M_tagSR_T5qqqqZH1700");
    h_J1M_tagSB_T5HH1700 = (TH1F*)fSignal->Get("J1pt_M_tagSB_T5qqqqZH1700"); h_J2M_tagSB_T5HH1700 = (TH1F*)fSignal->Get("J2pt_M_tagSB_T5qqqqZH1700");
    h_J1M_doubletagSR_T5HH2100 = (TH1F*)fSignal->Get("J1pt_M_doubletagSR_T5qqqqZH2100"); h_J2M_doubletagSR_T5HH2100 = (TH1F*)fSignal->Get("J2pt_M_doubletagSR_T5qqqqZH2100");
    h_J1M_doubletagSB_T5HH2100 = (TH1F*)fSignal->Get("J1pt_M_doubletagSB_T5qqqqZH2100"); h_J2M_doubletagSB_T5HH2100 = (TH1F*)fSignal->Get("J2pt_M_doubletagSB_T5qqqqZH2100");
    h_J1M_antitagSR_T5HH2100 = (TH1F*)fSignal->Get("J1pt_M_antitagSR_T5qqqqZH2100"); h_J2M_antitagSR_T5HH2100 = (TH1F*)fSignal->Get("J2pt_M_antitagSR_T5qqqqZH2100");
    h_J1M_antitagSB_T5HH2100 = (TH1F*)fSignal->Get("J1pt_M_antitagSB_T5qqqqZH2100"); h_J2M_antitagSB_T5HH2100 = (TH1F*)fSignal->Get("J2pt_M_antitagSB_T5qqqqZH2100");
    h_J1M_tagSR_T5HH2100 = (TH1F*)fSignal->Get("J1pt_M_tagSR_T5qqqqZH2100"); h_J2M_tagSR_T5HH2100 = (TH1F*)fSignal->Get("J2pt_M_tagSR_T5qqqqZH2100");
    h_J1M_tagSB_T5HH2100 = (TH1F*)fSignal->Get("J1pt_M_tagSB_T5qqqqZH2100"); h_J2M_tagSB_T5HH2100 = (TH1F*)fSignal->Get("J2pt_M_tagSB_T5qqqqZH2100");

    //For pt
    h_J1Pt_doubletagSR_sum = (TH1F*)f->Get("J1pt_Pt_doubletagSR_sum"); h_J2Pt_doubletagSR_sum = (TH1F*)f->Get("J2pt_Pt_doubletagSR_sum");
    h_J1Pt_doubletagSB_sum = (TH1F*)f->Get("J1pt_Pt_doubletagSB_sum"); h_J2Pt_doubletagSB_sum = (TH1F*)f->Get("J2pt_Pt_doubletagSB_sum");
    h_J1Pt_antitagSR_sum = (TH1F*)f->Get("J1pt_Pt_antitagSR_sum"); h_J2Pt_antitagSR_sum = (TH1F*)f->Get("J2pt_Pt_antitagSR_sum");
    h_J1Pt_antitagSB_sum = (TH1F*)f->Get("J1pt_Pt_antitagSB_sum"); h_J2Pt_antitagSB_sum = (TH1F*)f->Get("J2pt_Pt_antitagSB_sum");
    h_J1Pt_tagSR_sum = (TH1F*)f->Get("J1pt_Pt_tagSR_sum"); h_J2Pt_tagSR_sum = (TH1F*)f->Get("J2pt_Pt_tagSR_sum");
    h_J1Pt_tagSB_sum = (TH1F*)f->Get("J1pt_Pt_tagSB_sum"); h_J2Pt_tagSB_sum = (TH1F*)f->Get("J2pt_Pt_tagSB_sum");

    h_J1Pt_doubletagSR_ZJets = (TH1F*)f->Get("J1pt_Pt_doubletagSR_ZJets"); h_J2Pt_doubletagSR_ZJets = (TH1F*)f->Get("J2pt_Pt_doubletagSR_ZJets");
    h_J1Pt_doubletagSB_ZJets = (TH1F*)f->Get("J1pt_Pt_doubletagSB_ZJets"); h_J2Pt_doubletagSB_ZJets = (TH1F*)f->Get("J2pt_Pt_doubletagSB_ZJets");
    h_J1Pt_antitagSR_ZJets = (TH1F*)f->Get("J1pt_Pt_antitagSR_ZJets"); h_J2Pt_antitagSR_ZJets = (TH1F*)f->Get("J2pt_Pt_antitagSR_ZJets");
    h_J1Pt_antitagSB_ZJets = (TH1F*)f->Get("J1pt_Pt_antitagSB_ZJets"); h_J2Pt_antitagSB_ZJets = (TH1F*)f->Get("J2pt_Pt_antitagSB_ZJets");
    h_J1Pt_tagSR_ZJets = (TH1F*)f->Get("J1pt_Pt_tagSR_ZJets"); h_J2Pt_tagSR_ZJets = (TH1F*)f->Get("J2pt_Pt_tagSR_ZJets");
    h_J1Pt_tagSB_ZJets = (TH1F*)f->Get("J1pt_Pt_tagSB_ZJets"); h_J2Pt_tagSB_ZJets = (TH1F*)f->Get("J2pt_Pt_tagSB_ZJets");

    h_J1Pt_doubletagSR_QCD = (TH1F*)f->Get("J1pt_Pt_doubletagSR_QCD"); h_J2Pt_doubletagSR_QCD = (TH1F*)f->Get("J2pt_Pt_doubletagSR_QCD");
    h_J1Pt_doubletagSB_QCD = (TH1F*)f->Get("J1pt_Pt_doubletagSB_QCD"); h_J2Pt_doubletagSB_QCD = (TH1F*)f->Get("J2pt_Pt_doubletagSB_QCD");
    h_J1Pt_antitagSR_QCD = (TH1F*)f->Get("J1pt_Pt_antitagSR_QCD"); h_J2Pt_antitagSR_QCD = (TH1F*)f->Get("J2pt_Pt_antitagSR_QCD");
    h_J1Pt_antitagSB_QCD = (TH1F*)f->Get("J1pt_Pt_antitagSB_QCD"); h_J2Pt_antitagSB_QCD = (TH1F*)f->Get("J2pt_Pt_antitagSB_QCD");
    h_J1Pt_tagSR_QCD = (TH1F*)f->Get("J1pt_Pt_tagSR_QCD"); h_J2Pt_tagSR_QCD = (TH1F*)f->Get("J2pt_Pt_tagSR_QCD");
    h_J1Pt_tagSB_QCD = (TH1F*)f->Get("J1pt_Pt_tagSB_QCD"); h_J2Pt_tagSB_QCD = (TH1F*)f->Get("J2pt_Pt_tagSB_QCD");

    h_J1Pt_doubletagSR_TT = (TH1F*)f->Get("J1pt_Pt_doubletagSR_TT"); h_J2Pt_doubletagSR_TT = (TH1F*)f->Get("J2pt_Pt_doubletagSR_TT");
    h_J1Pt_doubletagSB_TT = (TH1F*)f->Get("J1pt_Pt_doubletagSB_TT"); h_J2Pt_doubletagSB_TT = (TH1F*)f->Get("J2pt_Pt_doubletagSB_TT");
    h_J1Pt_antitagSR_TT = (TH1F*)f->Get("J1pt_Pt_antitagSR_TT"); h_J2Pt_antitagSR_TT = (TH1F*)f->Get("J2pt_Pt_antitagSR_TT");
    h_J1Pt_antitagSB_TT = (TH1F*)f->Get("J1pt_Pt_antitagSB_TT"); h_J2Pt_antitagSB_TT = (TH1F*)f->Get("J2pt_Pt_antitagSB_TT");
    h_J1Pt_tagSR_TT = (TH1F*)f->Get("J1pt_Pt_tagSR_TT"); h_J2Pt_tagSR_TT = (TH1F*)f->Get("J2pt_Pt_tagSR_TT");
    h_J1Pt_tagSB_TT = (TH1F*)f->Get("J1pt_Pt_tagSB_TT"); h_J2Pt_tagSB_TT = (TH1F*)f->Get("J2pt_Pt_tagSB_TT");

    // h_J1Pt_doubletagSR_TT_SL = (TH1F*)f->Get("J1pt_Pt_doubletagSR_TT_SL"); h_J2Pt_doubletagSR_TT_SL = (TH1F*)f->Get("J2pt_Pt_doubletagSR_TT_SL");
    // h_J1Pt_doubletagSB_TT_SL = (TH1F*)f->Get("J1pt_Pt_doubletagSB_TT_SL"); h_J2Pt_doubletagSB_TT_SL = (TH1F*)f->Get("J2pt_Pt_doubletagSB_TT_SL");
    // h_J1Pt_antitagSR_TT_SL = (TH1F*)f->Get("J1pt_Pt_antitagSR_TT_SL"); h_J2Pt_antitagSR_TT_SL = (TH1F*)f->Get("J2pt_Pt_antitagSR_TT_SL");
    // h_J1Pt_antitagSB_TT_SL = (TH1F*)f->Get("J1pt_Pt_antitagSB_TT_SL"); h_J2Pt_antitagSB_TT_SL = (TH1F*)f->Get("J2pt_Pt_antitagSB_TT_SL");
    // h_J1Pt_tagSR_TT_SL = (TH1F*)f->Get("J1pt_Pt_tagSR_TT_SL"); h_J2Pt_tagSR_TT_SL = (TH1F*)f->Get("J2pt_Pt_tagSR_TT_SL");
    // h_J1Pt_tagSB_TT_SL = (TH1F*)f->Get("J1pt_Pt_tagSB_TT_SL"); h_J2Pt_tagSB_TT_SL = (TH1F*)f->Get("J2pt_Pt_tagSB_TT_SL");
    // h_J1Pt_doubletagSR_TT_Di = (TH1F*)f->Get("J1pt_Pt_doubletagSR_TT_Di"); h_J2Pt_doubletagSR_TT_Di = (TH1F*)f->Get("J2pt_Pt_doubletagSR_TT_Di");
    // h_J1Pt_doubletagSB_TT_Di = (TH1F*)f->Get("J1pt_Pt_doubletagSB_TT_Di"); h_J2Pt_doubletagSB_TT_Di = (TH1F*)f->Get("J2pt_Pt_doubletagSB_TT_Di");
    // h_J1Pt_antitagSR_TT_Di = (TH1F*)f->Get("J1pt_Pt_antitagSR_TT_Di"); h_J2Pt_antitagSR_TT_Di = (TH1F*)f->Get("J2pt_Pt_antitagSR_TT_Di");
    // h_J1Pt_antitagSB_TT_Di = (TH1F*)f->Get("J1pt_Pt_antitagSB_TT_Di"); h_J2Pt_antitagSB_TT_Di = (TH1F*)f->Get("J2pt_Pt_antitagSB_TT_Di");
    // h_J1Pt_tagSR_TT_Di = (TH1F*)f->Get("J1pt_Pt_tagSR_TT_Di"); h_J2Pt_tagSR_TT_Di = (TH1F*)f->Get("J2pt_Pt_tagSR_TT_Di");
    // h_J1Pt_tagSB_TT_Di = (TH1F*)f->Get("J1pt_Pt_tagSB_TT_Di"); h_J2Pt_tagSB_TT_Di = (TH1F*)f->Get("J2pt_Pt_tagSB_TT_Di");

    h_J1Pt_doubletagSR_WJets = (TH1F*)f->Get("J1pt_Pt_doubletagSR_WJets"); h_J2Pt_doubletagSR_WJets = (TH1F*)f->Get("J2pt_Pt_doubletagSR_WJets");
    h_J1Pt_doubletagSB_WJets = (TH1F*)f->Get("J1pt_Pt_doubletagSB_WJets"); h_J2Pt_doubletagSB_WJets = (TH1F*)f->Get("J2pt_Pt_doubletagSB_WJets");
    h_J1Pt_antitagSR_WJets = (TH1F*)f->Get("J1pt_Pt_antitagSR_WJets"); h_J2Pt_antitagSR_WJets = (TH1F*)f->Get("J2pt_Pt_antitagSR_WJets");
    h_J1Pt_antitagSB_WJets = (TH1F*)f->Get("J1pt_Pt_antitagSB_WJets"); h_J2Pt_antitagSB_WJets = (TH1F*)f->Get("J2pt_Pt_antitagSB_WJets");
    h_J1Pt_tagSR_WJets = (TH1F*)f->Get("J1pt_Pt_tagSR_WJets"); h_J2Pt_tagSR_WJets = (TH1F*)f->Get("J2pt_Pt_tagSR_WJets");
    h_J1Pt_tagSB_WJets = (TH1F*)f->Get("J1pt_Pt_tagSB_WJets"); h_J2Pt_tagSB_WJets = (TH1F*)f->Get("J2pt_Pt_tagSB_WJets");

    //for signal
    h_J1Pt_doubletagSR_T5HH1300 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSR_T5qqqqZH1300"); h_J2Pt_doubletagSR_T5HH1300 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSR_T5qqqqZH1300");
    h_J1Pt_doubletagSB_T5HH1300 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSB_T5qqqqZH1300"); h_J2Pt_doubletagSB_T5HH1300 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSB_T5qqqqZH1300");
    h_J1Pt_antitagSR_T5HH1300 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSR_T5qqqqZH1300"); h_J2Pt_antitagSR_T5HH1300 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSR_T5qqqqZH1300");
    h_J1Pt_antitagSB_T5HH1300 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSB_T5qqqqZH1300"); h_J2Pt_antitagSB_T5HH1300 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSB_T5qqqqZH1300");
    h_J1Pt_tagSR_T5HH1300 = (TH1F*)fSignal->Get("J1pt_Pt_tagSR_T5qqqqZH1300"); h_J2Pt_tagSR_T5HH1300 = (TH1F*)fSignal->Get("J2pt_Pt_tagSR_T5qqqqZH1300");
    h_J1Pt_tagSB_T5HH1300 = (TH1F*)fSignal->Get("J1pt_Pt_tagSB_T5qqqqZH1300"); h_J2Pt_tagSB_T5HH1300 = (TH1F*)fSignal->Get("J2pt_Pt_tagSB_T5qqqqZH1300");
    h_J1Pt_doubletagSR_T5HH1700 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSR_T5qqqqZH1700"); h_J2Pt_doubletagSR_T5HH1700 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSR_T5qqqqZH1700");
    h_J1Pt_doubletagSB_T5HH1700 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSB_T5qqqqZH1700"); h_J2Pt_doubletagSB_T5HH1700 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSB_T5qqqqZH1700");
    h_J1Pt_antitagSR_T5HH1700 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSR_T5qqqqZH1700"); h_J2Pt_antitagSR_T5HH1700 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSR_T5qqqqZH1700");
    h_J1Pt_antitagSB_T5HH1700 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSB_T5qqqqZH1700"); h_J2Pt_antitagSB_T5HH1700 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSB_T5qqqqZH1700");
    h_J1Pt_tagSR_T5HH1700 = (TH1F*)fSignal->Get("J1pt_Pt_tagSR_T5qqqqZH1700"); h_J2Pt_tagSR_T5HH1700 = (TH1F*)fSignal->Get("J2pt_Pt_tagSR_T5qqqqZH1700");
    h_J1Pt_tagSB_T5HH1700 = (TH1F*)fSignal->Get("J1pt_Pt_tagSB_T5qqqqZH1700"); h_J2Pt_tagSB_T5HH1700 = (TH1F*)fSignal->Get("J2pt_Pt_tagSB_T5qqqqZH1700");
    h_J1Pt_doubletagSR_T5HH2100 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSR_T5qqqqZH2100"); h_J2Pt_doubletagSR_T5HH2100 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSR_T5qqqqZH2100");
    h_J1Pt_doubletagSB_T5HH2100 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSB_T5qqqqZH2100"); h_J2Pt_doubletagSB_T5HH2100 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSB_T5qqqqZH2100");
    h_J1Pt_antitagSR_T5HH2100 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSR_T5qqqqZH2100"); h_J2Pt_antitagSR_T5HH2100 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSR_T5qqqqZH2100");
    h_J1Pt_antitagSB_T5HH2100 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSB_T5qqqqZH2100"); h_J2Pt_antitagSB_T5HH2100 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSB_T5qqqqZH2100");
    h_J1Pt_tagSR_T5HH2100 = (TH1F*)fSignal->Get("J1pt_Pt_tagSR_T5qqqqZH2100"); h_J2Pt_tagSR_T5HH2100 = (TH1F*)fSignal->Get("J2pt_Pt_tagSR_T5qqqqZH2100");
    h_J1Pt_tagSB_T5HH2100 = (TH1F*)fSignal->Get("J1pt_Pt_tagSB_T5qqqqZH2100"); h_J2Pt_tagSB_T5HH2100 = (TH1F*)fSignal->Get("J2pt_Pt_tagSB_T5qqqqZH2100");

    h_0HSBOpt1 = (TH1F*)f->Get("MET_antitagSBOpt1_sum"); h_0HSBOpt2 = (TH1F*)f->Get("MET_antitagSBOpt2_sum");
    h_0HSROpt1 = (TH1F*)f->Get("MET_antitagSROpt1_sum"); h_0HSROpt2 = (TH1F*)f->Get("MET_antitagSROpt2_sum");
  }

  if (whichRegion=="photon") {
    h_A_sum = (TH1F*)fPhoton->Get("MET_doubletagSR_sum"); h_B_sum = (TH1F*)fPhoton->Get("MET_doubletagSB_sum");
    h_A1_sum = (TH1F*)fPhoton->Get("MET_tagSR_sum"); h_B1_sum = (TH1F*)fPhoton->Get("MET_tagSB_sum");
    h_C_sum = (TH1F*)fPhoton->Get("MET_antitagSR_sum"); h_D_sum = (TH1F*)fPhoton->Get("MET_antitagSB_sum");

    h_A_QCD = (TH1F*)fPhoton->Get("MET_doubletagSR_QCD"); h_B_QCD = (TH1F*)fPhoton->Get("MET_doubletagSB_QCD");
    h_A1_QCD = (TH1F*)fPhoton->Get("MET_tagSR_QCD"); h_B1_QCD = (TH1F*)fPhoton->Get("MET_tagSB_QCD");
    h_C_QCD = (TH1F*)fPhoton->Get("MET_antitagSR_QCD"); h_D_QCD = (TH1F*)fPhoton->Get("MET_antitagSB_QCD");

    h_A_GJets = (TH1F*)fPhoton->Get("MET_doubletagSR_GJets"); h_B_GJets = (TH1F*)fPhoton->Get("MET_doubletagSB_GJets");
    h_A1_GJets = (TH1F*)fPhoton->Get("MET_tagSR_GJets"); h_B1_GJets = (TH1F*)fPhoton->Get("MET_tagSB_GJets");
    h_C_GJets = (TH1F*)fPhoton->Get("MET_antitagSR_GJets"); h_D_GJets = (TH1F*)fPhoton->Get("MET_antitagSB_GJets");

    h_J1M_doubletagSR_GJets = (TH1F*)fPhoton->Get("J1pt_M_doubletagSR_GJets"); h_J2M_doubletagSR_GJets = (TH1F*)fPhoton->Get("J2pt_M_doubletagSR_GJets");
    h_J1M_doubletagSB_GJets = (TH1F*)fPhoton->Get("J1pt_M_doubletagSB_GJets"); h_J2M_doubletagSB_GJets = (TH1F*)fPhoton->Get("J2pt_M_doubletagSB_GJets");
    h_J1M_antitagSR_GJets = (TH1F*)fPhoton->Get("J1pt_M_antitagSR_GJets"); h_J2M_antitagSR_GJets = (TH1F*)fPhoton->Get("J2pt_M_antitagSR_GJets");
    h_J1M_antitagSB_GJets = (TH1F*)fPhoton->Get("J1pt_M_antitagSB_GJets"); h_J2M_antitagSB_GJets = (TH1F*)fPhoton->Get("J2pt_M_antitagSB_GJets");
    h_J1M_tagSR_GJets = (TH1F*)fPhoton->Get("J1pt_M_tagSR_GJets"); h_J2M_tagSR_GJets = (TH1F*)fPhoton->Get("J2pt_M_tagSR_GJets");
    h_J1M_tagSB_GJets = (TH1F*)fPhoton->Get("J1pt_M_tagSB_GJets"); h_J2M_tagSB_GJets = (TH1F*)fPhoton->Get("J2pt_M_tagSB_GJets");

    h_J1M_doubletagSR_QCD = (TH1F*)fPhoton->Get("J1pt_M_doubletagSR_QCD"); h_J2M_doubletagSR_QCD = (TH1F*)fPhoton->Get("J2pt_M_doubletagSR_QCD");
    h_J1M_doubletagSB_QCD = (TH1F*)fPhoton->Get("J1pt_M_doubletagSB_QCD"); h_J2M_doubletagSB_QCD = (TH1F*)fPhoton->Get("J2pt_M_doubletagSB_QCD");
    h_J1M_antitagSR_QCD = (TH1F*)fPhoton->Get("J1pt_M_antitagSR_QCD"); h_J2M_antitagSR_QCD = (TH1F*)fPhoton->Get("J2pt_M_antitagSR_QCD");
    h_J1M_antitagSB_QCD = (TH1F*)fPhoton->Get("J1pt_M_antitagSB_QCD"); h_J2M_antitagSB_QCD = (TH1F*)fPhoton->Get("J2pt_M_antitagSB_QCD");
    h_J1M_tagSR_QCD = (TH1F*)fPhoton->Get("J1pt_M_tagSR_QCD"); h_J2M_tagSR_QCD = (TH1F*)fPhoton->Get("J2pt_M_tagSR_QCD");
    h_J1M_tagSB_QCD = (TH1F*)fPhoton->Get("J1pt_M_tagSB_QCD"); h_J2M_tagSB_QCD = (TH1F*)fPhoton->Get("J2pt_M_tagSB_QCD");

    h_J1M_doubletagSR_sum = (TH1F*)fPhoton->Get("J1pt_M_doubletagSR_sum"); h_J2M_doubletagSR_sum = (TH1F*)fPhoton->Get("J2pt_M_doubletagSR_sum");
    h_J1M_doubletagSB_sum = (TH1F*)fPhoton->Get("J1pt_M_doubletagSB_sum"); h_J2M_doubletagSB_sum = (TH1F*)fPhoton->Get("J2pt_M_doubletagSB_sum");
    h_J1M_antitagSR_sum = (TH1F*)fPhoton->Get("J1pt_M_antitagSR_sum"); h_J2M_antitagSR_sum = (TH1F*)fPhoton->Get("J2pt_M_antitagSR_sum");
    h_J1M_antitagSB_sum = (TH1F*)fPhoton->Get("J1pt_M_antitagSB_sum"); h_J2M_antitagSB_sum = (TH1F*)fPhoton->Get("J2pt_M_antitagSB_sum");
    h_J1M_tagSR_sum = (TH1F*)fPhoton->Get("J1pt_M_tagSR_sum"); h_J2M_tagSR_sum = (TH1F*)fPhoton->Get("J2pt_M_tagSR_sum");
    h_J1M_tagSB_sum = (TH1F*)fPhoton->Get("J1pt_M_tagSB_sum"); h_J2M_tagSB_sum = (TH1F*)fPhoton->Get("J2pt_M_tagSB_sum");

    h_J1Pt_doubletagSR_sum = (TH1F*)fPhoton->Get("J1pt_Pt_doubletagSR_sum"); h_J2Pt_doubletagSR_sum = (TH1F*)fPhoton->Get("J2pt_Pt_doubletagSR_sum");
    h_J1Pt_doubletagSB_sum = (TH1F*)fPhoton->Get("J1pt_Pt_doubletagSB_sum"); h_J2Pt_doubletagSB_sum = (TH1F*)fPhoton->Get("J2pt_Pt_doubletagSB_sum");
    h_J1Pt_antitagSR_sum = (TH1F*)fPhoton->Get("J1pt_Pt_antitagSR_sum"); h_J2Pt_antitagSR_sum = (TH1F*)fPhoton->Get("J2pt_Pt_antitagSR_sum");
    h_J1Pt_antitagSB_sum = (TH1F*)fPhoton->Get("J1pt_Pt_antitagSB_sum"); h_J2Pt_antitagSB_sum = (TH1F*)fPhoton->Get("J2pt_Pt_antitagSB_sum");
    h_J1Pt_tagSR_sum = (TH1F*)fPhoton->Get("J1pt_Pt_tagSR_sum"); h_J2Pt_tagSR_sum = (TH1F*)fPhoton->Get("J2pt_Pt_tagSR_sum");
    h_J1Pt_tagSB_sum = (TH1F*)fPhoton->Get("J1pt_Pt_tagSB_sum"); h_J2Pt_tagSB_sum = (TH1F*)fPhoton->Get("J2pt_Pt_tagSB_sum");

    h_J1Pt_doubletagSR_QCD = (TH1F*)fPhoton->Get("J1pt_Pt_doubletagSR_QCD"); h_J2Pt_doubletagSR_QCD = (TH1F*)fPhoton->Get("J2pt_Pt_doubletagSR_QCD");
    h_J1Pt_doubletagSB_QCD = (TH1F*)fPhoton->Get("J1pt_Pt_doubletagSB_QCD"); h_J2Pt_doubletagSB_QCD = (TH1F*)fPhoton->Get("J2pt_Pt_doubletagSB_QCD");
    h_J1Pt_antitagSR_QCD = (TH1F*)fPhoton->Get("J1pt_Pt_antitagSR_QCD"); h_J2Pt_antitagSR_QCD = (TH1F*)fPhoton->Get("J2pt_Pt_antitagSR_QCD");
    h_J1Pt_antitagSB_QCD = (TH1F*)fPhoton->Get("J1pt_Pt_antitagSB_QCD"); h_J2Pt_antitagSB_QCD = (TH1F*)fPhoton->Get("J2pt_Pt_antitagSB_QCD");
    h_J1Pt_tagSR_QCD = (TH1F*)fPhoton->Get("J1pt_Pt_tagSR_QCD"); h_J2Pt_tagSR_QCD = (TH1F*)fPhoton->Get("J2pt_Pt_tagSR_QCD");
    h_J1Pt_tagSB_QCD = (TH1F*)fPhoton->Get("J1pt_Pt_tagSB_QCD"); h_J2Pt_tagSB_QCD = (TH1F*)fPhoton->Get("J2pt_Pt_tagSB_QCD");

    h_J1Pt_doubletagSR_GJets = (TH1F*)fPhoton->Get("J1pt_Pt_doubletagSR_GJets"); h_J2Pt_doubletagSR_GJets = (TH1F*)fPhoton->Get("J2pt_Pt_doubletagSR_GJets");
    h_J1Pt_doubletagSB_GJets = (TH1F*)fPhoton->Get("J1pt_Pt_doubletagSB_GJets"); h_J2Pt_doubletagSB_GJets = (TH1F*)fPhoton->Get("J2pt_Pt_doubletagSB_GJets");
    h_J1Pt_antitagSR_GJets = (TH1F*)fPhoton->Get("J1pt_Pt_antitagSR_GJets"); h_J2Pt_antitagSR_GJets = (TH1F*)fPhoton->Get("J2pt_Pt_antitagSR_GJets");
    h_J1Pt_antitagSB_GJets = (TH1F*)fPhoton->Get("J1pt_Pt_antitagSB_GJets"); h_J2Pt_antitagSB_GJets = (TH1F*)fPhoton->Get("J2pt_Pt_antitagSB_GJets");
    h_J1Pt_tagSR_GJets = (TH1F*)fPhoton->Get("J1pt_Pt_tagSR_GJets"); h_J2Pt_tagSR_GJets = (TH1F*)fPhoton->Get("J2pt_Pt_tagSR_GJets");
    h_J1Pt_tagSB_GJets = (TH1F*)fPhoton->Get("J1pt_Pt_tagSB_GJets"); h_J2Pt_tagSB_GJets = (TH1F*)fPhoton->Get("J2pt_Pt_tagSB_GJets");

    h_J2M_mjBins_doubletagSR_sum = (TH1F*)fPhoton->Get("J2_M_jetBins_doubletagSR_sum"); h_J2M_mjBins_doubletagSB_sum = (TH1F*)fPhoton->Get("J2_M_jetBins_doubletagSB_sum");
    h_J2M_mjBins_antitagSR_sum = (TH1F*)fPhoton->Get("J2_M_jetBins_antitagSR_sum"); h_J2M_mjBins_antitagSB_sum = (TH1F*)fPhoton->Get("J2_M_jetBins_antitagSB_sum");
    h_J2M_mjBins_doubletagSR_QCD = (TH1F*)fPhoton->Get("J2_M_jetBins_doubletagSR_QCD"); h_J2M_mjBins_doubletagSB_QCD = (TH1F*)fPhoton->Get("J2_M_jetBins_doubletagSB_QCD");
    h_J2M_mjBins_antitagSR_QCD = (TH1F*)fPhoton->Get("J2_M_jetBins_antitagSR_QCD"); h_J2M_mjBins_antitagSB_QCD = (TH1F*)fPhoton->Get("J2_M_jetBins_antitagSB_QCD");
    h_J2M_mjBins_doubletagSR_GJets = (TH1F*)fPhoton->Get("J2_M_jetBins_doubletagSR_GJets"); h_J2M_mjBins_doubletagSB_GJets = (TH1F*)fPhoton->Get("J2_M_jetBins_doubletagSB_GJets");
    h_J2M_mjBins_antitagSR_GJets = (TH1F*)fPhoton->Get("J2_M_jetBins_antitagSR_GJets"); h_J2M_mjBins_antitagSB_GJets = (TH1F*)fPhoton->Get("J2_M_jetBins_antitagSB_GJets");

    h_0HSBOpt1 = (TH1F*)fPhoton->Get("MET_antitagSBOpt1_sum"); h_0HSBOpt2 = (TH1F*)fPhoton->Get("MET_antitagSBOpt2_sum");
    h_0HSROpt1 = (TH1F*)fPhoton->Get("MET_antitagSROpt1_sum"); h_0HSROpt2 = (TH1F*)fPhoton->Get("MET_antitagSROpt2_sum");
    h_2HLeadSR_sum = (TH1F*)fPhoton->Get("J2pt_M_doubletagSRLead_sum"); h_0HLeadSR_sum = (TH1F*)fPhoton->Get("J2pt_M_antitagSRLead_sum");
  }


  if (whichRegion=="singleLept") {
    // h_A_data = (TH1F*)fSingleLept->Get("MET_doubletagSR_data"); h_B_data = (TH1F*)fSingleLept->Get("MET_doubletagSB_data");
    // h_A1_data = (TH1F*)fSingleLept->Get("MET_tagSR_data"); h_B1_data = (TH1F*)fSingleLept->Get("MET_tagSB_data");
    // h_C_data = (TH1F*)fSingleLept->Get("MET_antitagSR_data"); h_D_data = (TH1F*)fSingleLept->Get("MET_antitagSB_data");

    // h_A_sum = (TH1F*)fSingleLept->Get("MET_doubletagSR_sum"); h_B_sum = (TH1F*)fSingleLept->Get("MET_doubletagSB_sum");
    // h_A1_sum = (TH1F*)fSingleLept->Get("MET_tagSR_sum"); h_B1_sum = (TH1F*)fSingleLept->Get("MET_tagSB_sum");
    // h_C_sum = (TH1F*)fSingleLept->Get("MET_antitagSR_sum"); h_D_sum = (TH1F*)fSingleLept->Get("MET_antitagSB_sum");
    h_A_TT = (TH1F*)fSingleLept->Get("MET_doubletagSR_TT"); h_B_TT = (TH1F*)fSingleLept->Get("MET_doubletagSB_TT");
    h_A1_TT = (TH1F*)fSingleLept->Get("MET_tagSR_TT"); h_B1_TT = (TH1F*)fSingleLept->Get("MET_tagSB_TT");
    h_C_TT = (TH1F*)fSingleLept->Get("MET_antitagSR_TT"); h_D_TT = (TH1F*)fSingleLept->Get("MET_antitagSB_TT");

    // h_A_TT_Di = (TH1F*)fSingleLept->Get("MET_doubletagSR_TT_Di"); h_B_TT_Di = (TH1F*)fSingleLept->Get("MET_doubletagSB_TT_Di");
    // h_A1_TT_Di = (TH1F*)fSingleLept->Get("MET_tagSR_TT_Di"); h_B1_TT_Di = (TH1F*)fSingleLept->Get("MET_tagSB_TT_Di");
    // h_C_TT_Di = (TH1F*)fSingleLept->Get("MET_antitagSR_TT_Di"); h_D_TT_Di = (TH1F*)fSingleLept->Get("MET_antitagSB_TT_Di");
    // h_A_TT_SL = (TH1F*)fSingleLept->Get("MET_doubletagSR_TT_SL"); h_B_TT_SL = (TH1F*)fSingleLept->Get("MET_doubletagSB_TT_SL");
    // h_A1_TT_SL = (TH1F*)fSingleLept->Get("MET_tagSR_TT_SL"); h_B1_TT_SL = (TH1F*)fSingleLept->Get("MET_tagSB_TT_SL");
    // h_C_TT_SL = (TH1F*)fSingleLept->Get("MET_antitagSR_TT_SL"); h_D_TT_SL = (TH1F*)fSingleLept->Get("MET_antitagSB_TT_SL");

    h_A_WJets = (TH1F*)fSingleLept->Get("MET_doubletagSR_WJets"); h_B_WJets = (TH1F*)fSingleLept->Get("MET_doubletagSB_WJets");
    h_A1_WJets = (TH1F*)fSingleLept->Get("MET_tagSR_WJets"); h_B1_WJets = (TH1F*)fSingleLept->Get("MET_tagSB_WJets");
    h_C_WJets = (TH1F*)fSingleLept->Get("MET_antitagSR_WJets");   h_D_WJets = (TH1F*)fSingleLept->Get("MET_antitagSB_WJets");

    h_A_sum = (TH1F*)h_A_TT->Clone("h_A_sum");  h_B_sum = (TH1F*)h_B_TT->Clone("h_B_sum");
    h_A1_sum = (TH1F*)h_A1_TT->Clone("h_A1_sum");  h_B1_sum = (TH1F*)h_B1_TT->Clone("h_B1_sum");
    h_C_sum = (TH1F*)h_C_TT->Clone("h_C_sum");  h_D_sum = (TH1F*)h_D_TT->Clone("h_D_sum");
    h_A_sum->Add(h_A_WJets); h_B_sum->Add(h_B_WJets);
    h_A1_sum->Add(h_A1_WJets); h_B1_sum->Add(h_B1_WJets);
    h_C_sum->Add(h_C_WJets); h_D_sum->Add(h_D_WJets);


    // h_J1M_doubletagSR_sum = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSR_sum"); h_J2M_doubletagSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSR_sum");
    // h_J1M_doubletagSB_sum = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSB_sum"); h_J2M_doubletagSB_sum = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSB_sum");
    // h_J1M_antitagSR_sum = (TH1F*)fSingleLept->Get("J1pt_M_antitagSR_sum"); h_J2M_antitagSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_antitagSR_sum");
    // h_J1M_antitagSB_sum = (TH1F*)fSingleLept->Get("J1pt_M_antitagSB_sum"); h_J2M_antitagSB_sum = (TH1F*)fSingleLept->Get("J2pt_M_antitagSB_sum");
    // h_J1M_tagSR_sum = (TH1F*)fSingleLept->Get("J1pt_M_tagSR_sum"); h_J2M_tagSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_tagSR_sum");
    // h_J1M_tagSB_sum = (TH1F*)fSingleLept->Get("J1pt_M_tagSB_sum"); h_J2M_tagSB_sum = (TH1F*)fSingleLept->Get("J2pt_M_tagSB_sum");

    h_J1M_doubletagSR_TT = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSR_TT"); h_J2M_doubletagSR_TT = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSR_TT");
    h_J1M_doubletagSB_TT = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSB_TT"); h_J2M_doubletagSB_TT = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSB_TT");
    h_J1M_antitagSR_TT = (TH1F*)fSingleLept->Get("J1pt_M_antitagSR_TT"); h_J2M_antitagSR_TT = (TH1F*)fSingleLept->Get("J2pt_M_antitagSR_TT");
    h_J1M_antitagSB_TT = (TH1F*)fSingleLept->Get("J1pt_M_antitagSB_TT"); h_J2M_antitagSB_TT = (TH1F*)fSingleLept->Get("J2pt_M_antitagSB_TT");
    h_J1M_tagSR_TT = (TH1F*)fSingleLept->Get("J1pt_M_tagSR_TT"); h_J2M_tagSR_TT = (TH1F*)fSingleLept->Get("J2pt_M_tagSR_TT");
    h_J1M_tagSB_TT = (TH1F*)fSingleLept->Get("J1pt_M_tagSB_TT"); h_J2M_tagSB_TT = (TH1F*)fSingleLept->Get("J2pt_M_tagSB_TT");

    // h_J1M_doubletagSR_TT_SL = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSR_TT_SL"); h_J2M_doubletagSR_TT_SL = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSR_TT_SL");
    // h_J1M_doubletagSB_TT_SL = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSB_TT_SL"); h_J2M_doubletagSB_TT_SL = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSB_TT_SL");
    // h_J1M_antitagSR_TT_SL = (TH1F*)fSingleLept->Get("J1pt_M_antitagSR_TT_SL"); h_J2M_antitagSR_TT_SL = (TH1F*)fSingleLept->Get("J2pt_M_antitagSR_TT_SL");
    // h_J1M_antitagSB_TT_SL = (TH1F*)fSingleLept->Get("J1pt_M_antitagSB_TT_SL"); h_J2M_antitagSB_TT_SL = (TH1F*)fSingleLept->Get("J2pt_M_antitagSB_TT_SL");
    // h_J1M_tagSR_TT_SL = (TH1F*)fSingleLept->Get("J1pt_M_tagSR_TT_SL"); h_J2M_tagSR_TT_SL = (TH1F*)fSingleLept->Get("J2pt_M_tagSR_TT_SL");
    // h_J1M_tagSB_TT_SL = (TH1F*)fSingleLept->Get("J1pt_M_tagSB_TT_SL"); h_J2M_tagSB_TT_SL = (TH1F*)fSingleLept->Get("J2pt_M_tagSB_TT_SL");
    // h_J1M_doubletagSR_TT_Di = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSR_TT_Di"); h_J2M_doubletagSR_TT_Di = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSR_TT_Di");
    // h_J1M_doubletagSB_TT_Di = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSB_TT_Di"); h_J2M_doubletagSB_TT_Di = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSB_TT_Di");
    // h_J1M_antitagSR_TT_Di = (TH1F*)fSingleLept->Get("J1pt_M_antitagSR_TT_Di"); h_J2M_antitagSR_TT_Di = (TH1F*)fSingleLept->Get("J2pt_M_antitagSR_TT_Di");
    // h_J1M_antitagSB_TT_Di = (TH1F*)fSingleLept->Get("J1pt_M_antitagSB_TT_Di"); h_J2M_antitagSB_TT_Di = (TH1F*)fSingleLept->Get("J2pt_M_antitagSB_TT_Di");
    // h_J1M_tagSR_TT_Di = (TH1F*)fSingleLept->Get("J1pt_M_tagSR_TT_Di"); h_J2M_tagSR_TT_Di = (TH1F*)fSingleLept->Get("J2pt_M_tagSR_TT_Di");
    // h_J1M_tagSB_TT_Di = (TH1F*)fSingleLept->Get("J1pt_M_tagSB_TT_Di"); h_J2M_tagSB_TT_Di = (TH1F*)fSingleLept->Get("J2pt_M_tagSB_TT_Di");

    h_J1M_doubletagSR_WJets = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSR_WJets"); h_J2M_doubletagSR_WJets = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSR_WJets");
    h_J1M_doubletagSB_WJets = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSB_WJets"); h_J2M_doubletagSB_WJets = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSB_WJets");
    h_J1M_antitagSR_WJets = (TH1F*)fSingleLept->Get("J1pt_M_antitagSR_WJets"); h_J2M_antitagSR_WJets = (TH1F*)fSingleLept->Get("J2pt_M_antitagSR_WJets");
    h_J1M_antitagSB_WJets = (TH1F*)fSingleLept->Get("J1pt_M_antitagSB_WJets"); h_J2M_antitagSB_WJets = (TH1F*)fSingleLept->Get("J2pt_M_antitagSB_WJets");
    h_J1M_tagSR_WJets = (TH1F*)fSingleLept->Get("J1pt_M_tagSR_WJets"); h_J2M_tagSR_WJets = (TH1F*)fSingleLept->Get("J2pt_M_tagSR_WJets");
    h_J1M_tagSB_WJets = (TH1F*)fSingleLept->Get("J1pt_M_tagSB_WJets"); h_J2M_tagSB_WJets = (TH1F*)fSingleLept->Get("J2pt_M_tagSB_WJets");

    h_J1Pt_doubletagSR_TT = (TH1F*)fSingleLept->Get("J1pt_Pt_doubletagSR_TT"); h_J2Pt_doubletagSR_TT = (TH1F*)fSingleLept->Get("J2pt_Pt_doubletagSR_TT");
    h_J1Pt_doubletagSB_TT = (TH1F*)fSingleLept->Get("J1pt_Pt_doubletagSB_TT"); h_J2Pt_doubletagSB_TT = (TH1F*)fSingleLept->Get("J2pt_Pt_doubletagSB_TT");
    h_J1Pt_antitagSR_TT = (TH1F*)fSingleLept->Get("J1pt_Pt_antitagSR_TT"); h_J2Pt_antitagSR_TT = (TH1F*)fSingleLept->Get("J2pt_Pt_antitagSR_TT");
    h_J1Pt_antitagSB_TT = (TH1F*)fSingleLept->Get("J1pt_Pt_antitagSB_TT"); h_J2Pt_antitagSB_TT = (TH1F*)fSingleLept->Get("J2pt_Pt_antitagSB_TT");
    h_J1Pt_tagSR_TT = (TH1F*)fSingleLept->Get("J1pt_Pt_tagSR_TT"); h_J2Pt_tagSR_TT = (TH1F*)fSingleLept->Get("J2pt_Pt_tagSR_TT");
    h_J1Pt_tagSB_TT = (TH1F*)fSingleLept->Get("J1pt_Pt_tagSB_TT"); h_J2Pt_tagSB_TT = (TH1F*)fSingleLept->Get("J2pt_Pt_tagSB_TT");

    // h_J1Pt_doubletagSR_TT_SL = (TH1F*)fSingleLept->Get("J1pt_Pt_doubletagSR_TT_SL"); h_J2Pt_doubletagSR_TT_SL = (TH1F*)fSingleLept->Get("J2pt_Pt_doubletagSR_TT_SL");
    // h_J1Pt_doubletagSB_TT_SL = (TH1F*)fSingleLept->Get("J1pt_Pt_doubletagSB_TT_SL"); h_J2Pt_doubletagSB_TT_SL = (TH1F*)fSingleLept->Get("J2pt_Pt_doubletagSB_TT_SL");
    // h_J1Pt_antitagSR_TT_SL = (TH1F*)fSingleLept->Get("J1pt_Pt_antitagSR_TT_SL"); h_J2Pt_antitagSR_TT_SL = (TH1F*)fSingleLept->Get("J2pt_Pt_antitagSR_TT_SL");
    // h_J1Pt_antitagSB_TT_SL = (TH1F*)fSingleLept->Get("J1pt_Pt_antitagSB_TT_SL"); h_J2Pt_antitagSB_TT_SL = (TH1F*)fSingleLept->Get("J2pt_Pt_antitagSB_TT_SL");
    // h_J1Pt_tagSR_TT_SL = (TH1F*)fSingleLept->Get("J1pt_Pt_tagSR_TT_SL"); h_J2Pt_tagSR_TT_SL = (TH1F*)fSingleLept->Get("J2pt_Pt_tagSR_TT_SL");
    // h_J1Pt_tagSB_TT_SL = (TH1F*)fSingleLept->Get("J1pt_Pt_tagSB_TT_SL"); h_J2Pt_tagSB_TT_SL = (TH1F*)fSingleLept->Get("J2pt_Pt_tagSB_TT_SL");
    // h_J1Pt_doubletagSR_TT_Di = (TH1F*)fSingleLept->Get("J1pt_Pt_doubletagSR_TT_Di"); h_J2Pt_doubletagSR_TT_Di = (TH1F*)fSingleLept->Get("J2pt_Pt_doubletagSR_TT_Di");
    // h_J1Pt_doubletagSB_TT_Di = (TH1F*)fSingleLept->Get("J1pt_Pt_doubletagSB_TT_Di"); h_J2Pt_doubletagSB_TT_Di = (TH1F*)fSingleLept->Get("J2pt_Pt_doubletagSB_TT_Di");
    // h_J1Pt_antitagSR_TT_Di = (TH1F*)fSingleLept->Get("J1pt_Pt_antitagSR_TT_Di"); h_J2Pt_antitagSR_TT_Di = (TH1F*)fSingleLept->Get("J2pt_Pt_antitagSR_TT_Di");
    // h_J1Pt_antitagSB_TT_Di = (TH1F*)fSingleLept->Get("J1pt_Pt_antitagSB_TT_Di"); h_J2Pt_antitagSB_TT_Di = (TH1F*)fSingleLept->Get("J2pt_Pt_antitagSB_TT_Di");
    // h_J1Pt_tagSR_TT_Di = (TH1F*)fSingleLept->Get("J1pt_Pt_tagSR_TT_Di"); h_J2Pt_tagSR_TT_Di = (TH1F*)fSingleLept->Get("J2pt_Pt_tagSR_TT_Di");
    // h_J1Pt_tagSB_TT_Di = (TH1F*)fSingleLept->Get("J1pt_Pt_tagSB_TT_Di"); h_J2Pt_tagSB_TT_Di = (TH1F*)fSingleLept->Get("J2pt_Pt_tagSB_TT_Di");

    h_J1Pt_doubletagSR_WJets = (TH1F*)fSingleLept->Get("J1pt_Pt_doubletagSR_WJets"); h_J2Pt_doubletagSR_WJets = (TH1F*)fSingleLept->Get("J2pt_Pt_doubletagSR_WJets");
    h_J1Pt_doubletagSB_WJets = (TH1F*)fSingleLept->Get("J1pt_Pt_doubletagSB_WJets"); h_J2Pt_doubletagSB_WJets = (TH1F*)fSingleLept->Get("J2pt_Pt_doubletagSB_WJets");
    h_J1Pt_antitagSR_WJets = (TH1F*)fSingleLept->Get("J1pt_Pt_antitagSR_WJets"); h_J2Pt_antitagSR_WJets = (TH1F*)fSingleLept->Get("J2pt_Pt_antitagSR_WJets");
    h_J1Pt_antitagSB_WJets = (TH1F*)fSingleLept->Get("J1pt_Pt_antitagSB_WJets"); h_J2Pt_antitagSB_WJets = (TH1F*)fSingleLept->Get("J2pt_Pt_antitagSB_WJets");
    h_J1Pt_tagSR_WJets = (TH1F*)fSingleLept->Get("J1pt_Pt_tagSR_WJets"); h_J2Pt_tagSR_WJets = (TH1F*)fSingleLept->Get("J2pt_Pt_tagSR_WJets");
    h_J1Pt_tagSB_WJets = (TH1F*)fSingleLept->Get("J1pt_Pt_tagSB_WJets"); h_J2Pt_tagSB_WJets = (TH1F*)fSingleLept->Get("J2pt_Pt_tagSB_WJets");

    // h_J2M_mjBins_doubletagSR_sum = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSR_sum"); h_J2M_mjBins_doubletagSB_sum = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSB_sum");
    // h_J2M_mjBins_antitagSR_sum = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSR_sum"); h_J2M_mjBins_antitagSB_sum = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSB_sum");
    h_J2M_mjBins_doubletagSR_TT = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSR_TT"); h_J2M_mjBins_doubletagSB_TT = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSB_TT");
    h_J2M_mjBins_antitagSR_TT = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSR_TT"); h_J2M_mjBins_antitagSB_TT = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSB_TT");

    // h_J2M_mjBins_doubletagSR_TT_Di = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_TT_Di"); h_J2M_mjBins_doubletagSB_TT_Di = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_TT_Di");
    // h_J2M_mjBins_antitagSR_TT_Di = (TH1F*)f->Get("J2_M_jetBins_antitagSR_TT_Di"); h_J2M_mjBins_antitagSB_TT_Di = (TH1F*)f->Get("J2_M_jetBins_antitagSB_TT_Di");
    // h_J2M_mjBins_doubletagSR_TT_SL = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_TT_SL"); h_J2M_mjBins_doubletagSB_TT_SL = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_TT_SL");
    // h_J2M_mjBins_antitagSR_TT_SL = (TH1F*)f->Get("J2_M_jetBins_antitagSR_TT_SL"); h_J2M_mjBins_antitagSB_TT_SL = (TH1F*)f->Get("J2_M_jetBins_antitagSB_TT_SL");

    h_J2M_mjBins_doubletagSR_WJets = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_WJets"); h_J2M_mjBins_doubletagSB_WJets = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_WJets");
    h_J2M_mjBins_antitagSR_WJets = (TH1F*)f->Get("J2_M_jetBins_antitagSR_WJets"); h_J2M_mjBins_antitagSB_WJets = (TH1F*)f->Get("J2_M_jetBins_antitagSB_WJets");
    //
    // h_0HSBOpt1 = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_sum"); //h_0HSBOpt2 = (TH1F*)fSingleLept->Get("MET_antitagSBOpt2_sum");
    // h_0HSROpt1 = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_sum"); //h_0HSROpt2 = (TH1F*)fSingleLept->Get("MET_antitagSROpt2_sum");

    h_0HSBOpt1_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_TT"); h_0HSBOpt1_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_WJets");
    h_0HSROpt1_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_TT"); h_0HSROpt1_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_WJets");

    // h_2HLeadSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSRLead_sum"); h_0HLeadSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_antitagSRLead_sum");

  }

  vector<TH1F*> histos_ABCD_data = {h_A_data, h_B_data, h_C_data, h_D_data};

  vector<TH1F*> histos_ABCD_sum = {h_A_sum, h_B_sum, h_C_sum, h_D_sum};
  vector<TH1F*> histos_ABCD_QCD = {h_A_QCD, h_B_QCD, h_C_QCD, h_D_QCD};
  vector<TH1F*> histos_ABCD_GJets = {h_A_GJets, h_B_GJets, h_C_GJets, h_D_GJets};
  vector<TH1F*> histos_ABCD_TT = {h_A_TT, h_B_TT, h_C_TT, h_D_TT};
  vector<TH1F*> histos_ABCD_TT_Di = {h_A_TT_Di, h_B_TT_Di, h_C_TT_Di, h_D_TT_Di};
  vector<TH1F*> histos_ABCD_TT_SL = {h_A_TT_SL, h_B_TT_SL, h_C_TT_SL, h_D_TT_SL};
  vector<TH1F*> histos_ABCD_WJets = {h_A_WJets, h_B_WJets, h_C_WJets, h_D_WJets};
  vector<TH1F*> histos_ABCD_ZJets = {h_A_ZJets, h_B_ZJets, h_C_ZJets, h_D_ZJets};

  vector<TH1F*> histos_ABCD_TChiHH200 = {h_A_TChiHH200, h_B_TChiHH200, h_C_TChiHH200, h_D_TChiHH200};
  vector<TH1F*> histos_ABCD_TChiHH400 = {h_A_TChiHH400, h_B_TChiHH400, h_C_TChiHH400, h_D_TChiHH400};
  vector<TH1F*> histos_ABCD_TChiHH700 = {h_A_TChiHH700, h_B_TChiHH700, h_C_TChiHH700, h_D_TChiHH700};
  vector<TH1F*> histos_ABCD_TChiHH1000 = {h_A_TChiHH1000, h_B_TChiHH1000, h_C_TChiHH1000, h_D_TChiHH1000};
  vector<TH1F*> histos_ABCD_TChiHH1300 = {h_A_TChiHH1300, h_B_TChiHH1300, h_C_TChiHH1300, h_D_TChiHH1300};
  vector<TH1F*> histos_ABCD_TChiHH1500 = {h_A_TChiHH1500, h_B_TChiHH1500, h_C_TChiHH1500, h_D_TChiHH1500};

  vector<TH1F*> histos_A1B1CD_sum = {h_A1_sum, h_B1_sum, h_C_sum, h_D_sum};
  vector<TH1F*> histos_A1B1CD_QCD = {h_A1_QCD, h_B1_QCD, h_C_QCD, h_D_QCD};
  vector<TH1F*> histos_A1B1CD_GJets = {h_A1_GJets, h_B1_GJets, h_C_GJets, h_D_GJets};
  vector<TH1F*> histos_A1B1CD_TT = {h_A1_TT, h_B1_TT, h_C_TT, h_D_TT};
  vector<TH1F*> histos_A1B1CD_TT_Di = {h_A1_TT_Di, h_B1_TT_Di, h_C_TT_Di, h_D_TT_Di};
  vector<TH1F*> histos_A1B1CD_TT_SL = {h_A1_TT_SL, h_B1_TT_SL, h_C_TT_SL, h_D_TT_SL};
  vector<TH1F*> histos_A1B1CD_WJets = {h_A1_WJets, h_B1_WJets, h_C_WJets, h_D_WJets};
  vector<TH1F*> histos_A1B1CD_ZJets = {h_A1_ZJets, h_B1_ZJets, h_C_ZJets, h_D_ZJets};

  vector<TH1F*> histos_A1B1CD_TChiHH200 = {h_A1_TChiHH200, h_B1_TChiHH200, h_C_TChiHH200, h_D_TChiHH200};
  vector<TH1F*> histos_A1B1CD_TChiHH400 = {h_A1_TChiHH400, h_B1_TChiHH400, h_C_TChiHH400, h_D_TChiHH400};
  vector<TH1F*> histos_A1B1CD_TChiHH700 = {h_A1_TChiHH700, h_B1_TChiHH700, h_C_TChiHH700, h_D_TChiHH700};
  vector<TH1F*> histos_A1B1CD_TChiHH1000 = {h_A1_TChiHH1000, h_B1_TChiHH1000, h_C_TChiHH1000, h_D_TChiHH1000};
  vector<TH1F*> histos_A1B1CD_TChiHH1300 = {h_A1_TChiHH1300, h_B1_TChiHH1300, h_C_TChiHH1300, h_D_TChiHH1300};
  vector<TH1F*> histos_A1B1CD_TChiHH1500 = {h_A1_TChiHH1500, h_B1_TChiHH1500, h_C_TChiHH1500, h_D_TChiHH1500};

  vector<TH1F*> histos_RPF_LeadSR_sum = {h_2HLeadSR_sum,h_0HLeadSR_sum};
  vector<TH1F*> histos_RPF_LeadSB_sum = {h_2HLeadSB_sum,h_0HLeadSB_sum};

  // vector<TH1F*> histos_RPF_LeadSR_sum = {h_J2M_doubletagSR_sum,h_2HLeadSR_sum,h_J2M_antitagSR_sum,h_0HLeadSR_sum};
  // vector<TH1F*> histos_RPF_LeadSB_sum = {h_J2M_doubletagSR_sum,h_2HLeadSB_sum,h_J2M_antitagSR_sum,h_0HLeadSB_sum};

  vector<TH1F*> histos_Rpf_J1_sum = {h_J1M_doubletagSR_sum, h_J1M_doubletagSB_sum, h_J1M_antitagSR_sum, h_J1M_antitagSB_sum};
  vector<TH1F*> histos_Rpf_J2_sum = {h_J2M_doubletagSR_sum, h_J2M_doubletagSB_sum, h_J2M_antitagSR_sum, h_J2M_antitagSB_sum};
  vector<TH1F*> histos_Rpf_J2_mJBins_sum = {h_J2M_mjBins_doubletagSR_sum, h_J2M_mjBins_doubletagSB_sum, h_J2M_mjBins_antitagSR_sum, h_J2M_mjBins_antitagSB_sum};
  vector<TH1F*> histos_Rpfsingle_J1_sum = {h_J1M_tagSR_sum, h_J1M_tagSB_sum, h_J1M_antitagSR_sum, h_J1M_antitagSB_sum};
  vector<TH1F*> histos_Rpfsingle_J2_sum = {h_J2M_tagSR_sum, h_J2M_tagSB_sum, h_J2M_antitagSR_sum, h_J2M_antitagSB_sum};

  vector<TH1F*> histos_Rpf_J1_QCD = {h_J1M_doubletagSR_QCD, h_J1M_doubletagSB_QCD, h_J1M_antitagSR_QCD, h_J1M_antitagSB_QCD};
  vector<TH1F*> histos_Rpf_J2_QCD = {h_J2M_doubletagSR_QCD, h_J2M_doubletagSB_QCD, h_J2M_antitagSR_QCD, h_J2M_antitagSB_QCD};
  vector<TH1F*> histos_Rpf_J2_mJBins_QCD = {h_J2M_mjBins_doubletagSR_QCD, h_J2M_mjBins_doubletagSB_QCD, h_J2M_mjBins_antitagSR_QCD, h_J2M_mjBins_antitagSB_QCD};
  vector<TH1F*> histos_Rpfsingle_J1_QCD = {h_J1M_tagSR_QCD, h_J1M_tagSB_QCD, h_J1M_antitagSR_QCD, h_J1M_antitagSB_QCD};
  vector<TH1F*> histos_Rpfsingle_J2_QCD = {h_J2M_tagSR_QCD, h_J2M_tagSB_QCD, h_J2M_antitagSR_QCD, h_J2M_antitagSB_QCD};

  vector<TH1F*> histos_Rpf_J1_GJets = {h_J1M_doubletagSR_GJets, h_J1M_doubletagSB_GJets, h_J1M_antitagSR_GJets, h_J1M_antitagSB_GJets};
  vector<TH1F*> histos_Rpf_J2_GJets = {h_J2M_doubletagSR_GJets, h_J2M_doubletagSB_GJets, h_J2M_antitagSR_GJets, h_J2M_antitagSB_GJets};
  vector<TH1F*> histos_Rpfsingle_J1_GJets = {h_J1M_tagSR_GJets, h_J1M_tagSB_GJets, h_J1M_antitagSR_GJets, h_J1M_antitagSB_GJets};
  vector<TH1F*> histos_Rpfsingle_J2_GJets = {h_J2M_tagSR_GJets, h_J2M_tagSB_GJets, h_J2M_antitagSR_GJets, h_J2M_antitagSB_GJets};

  vector<TH1F*> histos_Rpf_J1_TT = {h_J1M_doubletagSR_TT, h_J1M_doubletagSB_TT, h_J1M_antitagSR_TT, h_J1M_antitagSB_TT};
  vector<TH1F*> histos_Rpf_J2_TT = {h_J2M_doubletagSR_TT, h_J2M_doubletagSB_TT, h_J2M_antitagSR_TT, h_J2M_antitagSB_TT};
  vector<TH1F*> histos_Rpf_J2_mJBins_TT = {h_J2M_mjBins_doubletagSR_TT, h_J2M_mjBins_doubletagSB_TT, h_J2M_mjBins_antitagSR_TT, h_J2M_mjBins_antitagSB_TT};
  vector<TH1F*> histos_Rpfsingle_J1_TT = {h_J1M_tagSR_TT, h_J1M_tagSB_TT, h_J1M_antitagSR_TT, h_J1M_antitagSB_TT};
  vector<TH1F*> histos_Rpfsingle_J2_TT = {h_J2M_tagSR_TT, h_J2M_tagSB_TT, h_J2M_antitagSR_TT, h_J2M_antitagSB_TT};

  vector<TH1F*> histos_Rpf_J1_TT_Di = {h_J1M_doubletagSR_TT_Di, h_J1M_doubletagSB_TT_Di, h_J1M_antitagSR_TT_Di, h_J1M_antitagSB_TT_Di};
  vector<TH1F*> histos_Rpf_J2_TT_Di = {h_J2M_doubletagSR_TT_Di, h_J2M_doubletagSB_TT_Di, h_J2M_antitagSR_TT_Di, h_J2M_antitagSB_TT_Di};
  vector<TH1F*> histos_Rpf_J2_mJBins_TT_Di = {h_J2M_mjBins_doubletagSR_TT_Di, h_J2M_mjBins_doubletagSB_TT_Di, h_J2M_mjBins_antitagSR_TT_Di, h_J2M_mjBins_antitagSB_TT_Di};
  vector<TH1F*> histos_Rpfsingle_J1_TT_Di = {h_J1M_tagSR_TT_Di, h_J1M_tagSB_TT_Di, h_J1M_antitagSR_TT_Di, h_J1M_antitagSB_TT_Di};
  vector<TH1F*> histos_Rpfsingle_J2_TT_Di = {h_J2M_tagSR_TT_Di, h_J2M_tagSB_TT_Di, h_J2M_antitagSR_TT_Di, h_J2M_antitagSB_TT_Di};

  vector<TH1F*> histos_Rpf_J1_TT_SL = {h_J1M_doubletagSR_TT_SL, h_J1M_doubletagSB_TT_SL, h_J1M_antitagSR_TT_SL, h_J1M_antitagSB_TT_SL};
  vector<TH1F*> histos_Rpf_J2_TT_SL = {h_J2M_doubletagSR_TT_SL, h_J2M_doubletagSB_TT_SL, h_J2M_antitagSR_TT_SL, h_J2M_antitagSB_TT_SL};
  vector<TH1F*> histos_Rpf_J2_mJBins_TT_SL = {h_J2M_mjBins_doubletagSR_TT_SL, h_J2M_mjBins_doubletagSB_TT_SL, h_J2M_mjBins_antitagSR_TT_SL, h_J2M_mjBins_antitagSB_TT_SL};
  vector<TH1F*> histos_Rpfsingle_J1_TT_SL = {h_J1M_tagSR_TT_SL, h_J1M_tagSB_TT_SL, h_J1M_antitagSR_TT_SL, h_J1M_antitagSB_TT_SL};
  vector<TH1F*> histos_Rpfsingle_J2_TT_SL = {h_J2M_tagSR_TT_SL, h_J2M_tagSB_TT_SL, h_J2M_antitagSR_TT_SL, h_J2M_antitagSB_TT_SL};

  vector<TH1F*> histos_Rpf_J1_WJets = {h_J1M_doubletagSR_WJets, h_J1M_doubletagSB_WJets, h_J1M_antitagSR_WJets, h_J1M_antitagSB_WJets};
  vector<TH1F*> histos_Rpf_J2_WJets = {h_J2M_doubletagSR_WJets, h_J2M_doubletagSB_WJets, h_J2M_antitagSR_WJets, h_J2M_antitagSB_WJets};
  vector<TH1F*> histos_Rpf_J2_mJBins_WJets = {h_J2M_mjBins_doubletagSR_WJets, h_J2M_mjBins_doubletagSB_WJets, h_J2M_mjBins_antitagSR_WJets, h_J2M_mjBins_antitagSB_WJets};
  vector<TH1F*> histos_Rpfsingle_J1_WJets = {h_J1M_tagSR_WJets, h_J1M_tagSB_WJets, h_J1M_antitagSR_WJets, h_J1M_antitagSB_WJets};
  vector<TH1F*> histos_Rpfsingle_J2_WJets = {h_J2M_tagSR_WJets, h_J2M_tagSB_WJets, h_J2M_antitagSR_WJets, h_J2M_antitagSB_WJets};

  vector<TH1F*> histos_Rpf_J1_ZJets = {h_J1M_doubletagSR_ZJets, h_J1M_doubletagSB_ZJets, h_J1M_antitagSR_ZJets, h_J1M_antitagSB_ZJets};
  vector<TH1F*> histos_Rpf_J2_ZJets = {h_J2M_doubletagSR_ZJets, h_J2M_doubletagSB_ZJets, h_J2M_antitagSR_ZJets, h_J2M_antitagSB_ZJets};
  vector<TH1F*> histos_Rpf_J2_mJBins_ZJets = {h_J2M_mjBins_doubletagSR_ZJets, h_J2M_mjBins_doubletagSB_ZJets, h_J2M_mjBins_antitagSR_ZJets, h_J2M_mjBins_antitagSB_ZJets};
  vector<TH1F*> histos_Rpfsingle_J1_ZJets = {h_J1M_tagSR_ZJets, h_J1M_tagSB_ZJets, h_J1M_antitagSR_ZJets, h_J1M_antitagSB_ZJets};
  vector<TH1F*> histos_Rpfsingle_J2_ZJets = {h_J2M_tagSR_ZJets, h_J2M_tagSB_ZJets, h_J2M_antitagSR_ZJets, h_J2M_antitagSB_ZJets};

  vector<TH1F*> histos_j1mass_ZJets = {h_J1M_doubletagSR_ZJets,h_J1M_doubletagSB_ZJets,h_J1M_tagSR_ZJets,h_J1M_tagSB_ZJets, h_J1M_antitagSR_ZJets, h_J1M_antitagSB_ZJets};
  vector<TH1F*> histos_j1mass_WJets = {h_J1M_doubletagSR_WJets,h_J1M_doubletagSB_WJets,h_J1M_tagSR_WJets,h_J1M_tagSB_WJets, h_J1M_antitagSR_WJets, h_J1M_antitagSB_WJets};
  vector<TH1F*> histos_j1mass_TT    = {h_J1M_doubletagSR_TT,h_J1M_doubletagSB_TT,h_J1M_tagSR_TT,h_J1M_tagSB_TT, h_J1M_antitagSR_TT, h_J1M_antitagSB_TT};
  vector<TH1F*> histos_j1mass_TT_Di    = {h_J1M_doubletagSR_TT_Di,h_J1M_doubletagSB_TT_Di,h_J1M_tagSR_TT_Di,h_J1M_tagSB_TT_Di, h_J1M_antitagSR_TT_Di, h_J1M_antitagSB_TT_Di};
  vector<TH1F*> histos_j1mass_TT_SL    = {h_J1M_doubletagSR_TT_SL,h_J1M_doubletagSB_TT_SL,h_J1M_tagSR_TT_SL,h_J1M_tagSB_TT_SL, h_J1M_antitagSR_TT_SL, h_J1M_antitagSB_TT_SL};
  vector<TH1F*> histos_j1mass_QCD   = {h_J1M_doubletagSR_QCD,h_J1M_doubletagSB_QCD,h_J1M_tagSR_QCD,h_J1M_tagSB_QCD, h_J1M_antitagSR_QCD, h_J1M_antitagSB_QCD};
  vector<TH1F*> histos_j1mass_T5HH1300 = {h_J1M_doubletagSR_T5HH1300,h_J1M_doubletagSB_T5HH1300,h_J1M_tagSR_T5HH1300,h_J1M_tagSB_T5HH1300, h_J1M_antitagSR_T5HH1300, h_J1M_antitagSB_T5HH1300};
  vector<TH1F*> histos_j1mass_T5HH1700 = {h_J1M_doubletagSR_T5HH1700,h_J1M_doubletagSB_T5HH1700,h_J1M_tagSR_T5HH1700,h_J1M_tagSB_T5HH1700, h_J1M_antitagSR_T5HH1700, h_J1M_antitagSB_T5HH1700};
  vector<TH1F*> histos_j1mass_T5HH2100 = {h_J1M_doubletagSR_T5HH2100,h_J1M_doubletagSB_T5HH2100,h_J1M_tagSR_T5HH2100,h_J1M_tagSB_T5HH2100, h_J1M_antitagSR_T5HH2100, h_J1M_antitagSB_T5HH2100};

  vector<TH1F*> histos_j2mass_ZJets = {h_J2M_doubletagSR_ZJets,h_J2M_doubletagSB_ZJets,h_J2M_tagSR_ZJets,h_J2M_tagSB_ZJets, h_J2M_antitagSR_ZJets, h_J2M_antitagSB_ZJets};
  vector<TH1F*> histos_j2mass_WJets = {h_J2M_doubletagSR_WJets,h_J2M_doubletagSB_WJets,h_J2M_tagSR_WJets,h_J2M_tagSB_WJets, h_J2M_antitagSR_WJets, h_J2M_antitagSB_WJets};
  vector<TH1F*> histos_j2mass_TT    = {h_J2M_doubletagSR_TT,h_J2M_doubletagSB_TT,h_J2M_tagSR_TT,h_J2M_tagSB_TT, h_J2M_antitagSR_TT, h_J2M_antitagSB_TT};
  vector<TH1F*> histos_j2mass_TT_Di    = {h_J2M_doubletagSR_TT_Di,h_J2M_doubletagSB_TT_Di,h_J2M_tagSR_TT_Di,h_J2M_tagSB_TT_Di, h_J2M_antitagSR_TT_Di, h_J2M_antitagSB_TT_Di};
  vector<TH1F*> histos_j2mass_TT_SL    = {h_J2M_doubletagSR_TT_SL,h_J2M_doubletagSB_TT_SL,h_J2M_tagSR_TT_SL,h_J2M_tagSB_TT_SL, h_J2M_antitagSR_TT_SL, h_J2M_antitagSB_TT_SL};
  vector<TH1F*> histos_j2mass_QCD   = {h_J2M_doubletagSR_QCD,h_J2M_doubletagSB_QCD,h_J2M_tagSR_QCD,h_J2M_tagSB_QCD, h_J2M_antitagSR_QCD, h_J2M_antitagSB_QCD};
  vector<TH1F*> histos_j2mass_T5HH1300 = {h_J2M_doubletagSR_T5HH1300,h_J2M_doubletagSB_T5HH1300,h_J2M_tagSR_T5HH1300,h_J2M_tagSB_T5HH1300, h_J2M_antitagSR_T5HH1300, h_J2M_antitagSB_T5HH1300};
  vector<TH1F*> histos_j2mass_T5HH1700 = {h_J2M_doubletagSR_T5HH1700,h_J2M_doubletagSB_T5HH1700,h_J2M_tagSR_T5HH1700,h_J2M_tagSB_T5HH1700, h_J2M_antitagSR_T5HH1700, h_J2M_antitagSB_T5HH1700};
  vector<TH1F*> histos_j2mass_T5HH2100 = {h_J2M_doubletagSR_T5HH2100,h_J2M_doubletagSB_T5HH2100,h_J2M_tagSR_T5HH2100,h_J2M_tagSB_T5HH2100, h_J2M_antitagSR_T5HH2100, h_J2M_antitagSB_T5HH2100};

  vector<TH1F*> histos_j1pt_ZJets = {h_J1Pt_doubletagSR_ZJets,h_J1Pt_doubletagSB_ZJets,h_J1Pt_tagSR_ZJets,h_J1Pt_tagSB_ZJets, h_J1Pt_antitagSR_ZJets, h_J1Pt_antitagSB_ZJets};
  vector<TH1F*> histos_j1pt_WJets = {h_J1Pt_doubletagSR_WJets,h_J1Pt_doubletagSB_WJets,h_J1Pt_tagSR_WJets,h_J1Pt_tagSB_WJets, h_J1Pt_antitagSR_WJets, h_J1Pt_antitagSB_WJets};
  vector<TH1F*> histos_j1pt_TT    = {h_J1Pt_doubletagSR_TT,h_J1Pt_doubletagSB_TT,h_J1Pt_tagSR_TT,h_J1Pt_tagSB_TT, h_J1Pt_antitagSR_TT, h_J1Pt_antitagSB_TT};
  vector<TH1F*> histos_j1pt_TT_Di    = {h_J1Pt_doubletagSR_TT_Di,h_J1Pt_doubletagSB_TT_Di,h_J1Pt_tagSR_TT_Di,h_J1Pt_tagSB_TT_Di, h_J1Pt_antitagSR_TT_Di, h_J1Pt_antitagSB_TT_Di};
  vector<TH1F*> histos_j1pt_TT_SL    = {h_J1Pt_doubletagSR_TT_SL,h_J1Pt_doubletagSB_TT_SL,h_J1Pt_tagSR_TT_SL,h_J1Pt_tagSB_TT_SL, h_J1Pt_antitagSR_TT_SL, h_J1Pt_antitagSB_TT_SL};
  vector<TH1F*> histos_j1pt_QCD   = {h_J1Pt_doubletagSR_QCD,h_J1Pt_doubletagSB_QCD,h_J1Pt_tagSR_QCD,h_J1Pt_tagSB_QCD, h_J1Pt_antitagSR_QCD, h_J1Pt_antitagSB_QCD};
  vector<TH1F*> histos_j1pt_T5HH1300 = {h_J1Pt_doubletagSR_T5HH1300,h_J1Pt_doubletagSB_T5HH1300,h_J1Pt_tagSR_T5HH1300,h_J1Pt_tagSB_T5HH1300, h_J1Pt_antitagSR_T5HH1300, h_J1Pt_antitagSB_T5HH1300};
  vector<TH1F*> histos_j1pt_T5HH1700 = {h_J1Pt_doubletagSR_T5HH1700,h_J1Pt_doubletagSB_T5HH1700,h_J1Pt_tagSR_T5HH1700,h_J1Pt_tagSB_T5HH1700, h_J1Pt_antitagSR_T5HH1700, h_J1Pt_antitagSB_T5HH1700};
  vector<TH1F*> histos_j1pt_T5HH2100 = {h_J1Pt_doubletagSR_T5HH2100,h_J1Pt_doubletagSB_T5HH2100,h_J1Pt_tagSR_T5HH2100,h_J1Pt_tagSB_T5HH2100, h_J1Pt_antitagSR_T5HH2100, h_J1Pt_antitagSB_T5HH2100};

  vector<TH1F*> histos_j2pt_ZJets = {h_J2Pt_doubletagSR_ZJets,h_J2Pt_doubletagSB_ZJets,h_J2Pt_tagSR_ZJets,h_J2Pt_tagSB_ZJets, h_J2Pt_antitagSR_ZJets, h_J2Pt_antitagSB_ZJets};
  vector<TH1F*> histos_j2pt_WJets = {h_J2Pt_doubletagSR_WJets,h_J2Pt_doubletagSB_WJets,h_J2Pt_tagSR_WJets,h_J2Pt_tagSB_WJets, h_J2Pt_antitagSR_WJets, h_J2Pt_antitagSB_WJets};
  vector<TH1F*> histos_j2pt_TT    = {h_J2Pt_doubletagSR_TT,h_J2Pt_doubletagSB_TT,h_J2Pt_tagSR_TT,h_J2Pt_tagSB_TT, h_J2Pt_antitagSR_TT, h_J2Pt_antitagSB_TT};
  vector<TH1F*> histos_j2pt_TT_Di    = {h_J2Pt_doubletagSR_TT_Di,h_J2Pt_doubletagSB_TT_Di,h_J2Pt_tagSR_TT_Di,h_J2Pt_tagSB_TT_Di, h_J2Pt_antitagSR_TT_Di, h_J2Pt_antitagSB_TT_Di};
  vector<TH1F*> histos_j2pt_TT_SL    = {h_J2Pt_doubletagSR_TT_SL,h_J2Pt_doubletagSB_TT_SL,h_J2Pt_tagSR_TT_SL,h_J2Pt_tagSB_TT_SL, h_J2Pt_antitagSR_TT_SL, h_J2Pt_antitagSB_TT_SL};
  vector<TH1F*> histos_j2pt_QCD   = {h_J2Pt_doubletagSR_QCD,h_J2Pt_doubletagSB_QCD,h_J2Pt_tagSR_QCD,h_J2Pt_tagSB_QCD, h_J2Pt_antitagSR_QCD, h_J2Pt_antitagSB_QCD};
  vector<TH1F*> histos_j2pt_T5HH1300 = {h_J2Pt_doubletagSR_T5HH1300,h_J2Pt_doubletagSB_T5HH1300,h_J2Pt_tagSR_T5HH1300,h_J2Pt_tagSB_T5HH1300, h_J2Pt_antitagSR_T5HH1300, h_J2Pt_antitagSB_T5HH1300};
  vector<TH1F*> histos_j2pt_T5HH1700 = {h_J2Pt_doubletagSR_T5HH1700,h_J2Pt_doubletagSB_T5HH1700,h_J2Pt_tagSR_T5HH1700,h_J2Pt_tagSB_T5HH1700, h_J2Pt_antitagSR_T5HH1700, h_J2Pt_antitagSB_T5HH1700};
  vector<TH1F*> histos_j2pt_T5HH2100 = {h_J2Pt_doubletagSR_T5HH2100,h_J2Pt_doubletagSB_T5HH2100,h_J2Pt_tagSR_T5HH2100,h_J2Pt_tagSB_T5HH2100, h_J2Pt_antitagSR_T5HH2100, h_J2Pt_antitagSB_T5HH2100};

  vector<TH1F*> histos_assump_sum = {h_notagSR_sum,h_notagSB_sum};
  vector<TH1F*> histos_assump2_sum = {h_notagSR_sum,h_notagSB2_sum};
  vector<TH1F*> histos_assump_ZJets = {h_notagSR_ZJets,h_notagSB_ZJets};
  vector<TH1F*> histos_assump_WJets = {h_notagSR_WJets,h_notagSB_WJets};
  vector<TH1F*> histos_assump_TT = {h_notagSR_TT,h_notagSB_TT};
  vector<TH1F*> histos_assump_TT_Di = {h_notagSR_TT_Di,h_notagSB_TT_Di};
  vector<TH1F*> histos_assump_TT_SL = {h_notagSR_TT_SL,h_notagSB_TT_SL};
  vector<TH1F*> histos_assump_QCD = {h_notagSR_QCD,h_notagSB_QCD};

  // vector<TH1F*> h_METShape_2H = {h_A_sum,h_B_sum,h_C_sum,h_D_sum};
  // vector<TH1F*> h_METShape_2H = {h_A_TT,h_A_WJets,h_0HSROpt1_TT,h_0HSROpt1_WJets,h_0HSBOpt1_TT,h_0HSBOpt1_WJets};
  vector<TH1F*> h_METShape_2H = {h_A_sum,h_D_sum,h_0HSBOpt1,h_C_sum,h_0HSROpt1};
  vector<TH1F*> h_METShape_1H = {h_A1_sum,h_D_sum,h_0HSBOpt1,h_C_sum,h_0HSROpt1};

  //Revisit optional cuts when V18 is out
  // vector<TH1F*> h_METShape_2H = {h_A_sum,h_D_sum,h_0HSBOpt2,h_C_sum,h_0HSROpt2};
  // vector<TH1F*> h_METShape_1H = {h_A1_sum,h_D_sum,h_0HSBOpt2,h_C_sum,h_0HSROpt2};


  if (runABCDPlots){
    if (whichRegion=="signal"){
      // makeFullBkgClosure(histos_ABCD_sum, "BkgSum", "Double",year);
      makeABCDPlot(histos_ABCD_sum, "BkgSum", "Double",year);
      makeABCDPlot(histos_ABCD_QCD, "QCD", "Double",year);
      makeABCDPlot(histos_ABCD_TT, "TT", "Double",year);
      // makeABCDPlot(histos_ABCD_TT_Di, "TT_Di", "Double",year);
      // makeABCDPlot(histos_ABCD_TT_SL, "TT_SL", "Double",year);
      makeABCDPlot(histos_ABCD_WJets, "WJets", "Double",year);
      makeABCDPlot(histos_ABCD_ZJets, "ZJets", "Double",year);
      makeABCDPlot(histos_A1B1CD_sum, "BkgSum", "Single",year);
      makeABCDPlot(histos_A1B1CD_QCD, "QCD", "Single",year);
      makeABCDPlot(histos_A1B1CD_TT, "TT", "Single",year);
      // makeABCDPlot(histos_A1B1CD_TT_Di, "TT_Di", "Single",year);
      // makeABCDPlot(histos_A1B1CD_TT_SL, "TT_SL", "Single",year);
      makeABCDPlot(histos_A1B1CD_WJets, "WJets", "Single",year);
      makeABCDPlot(histos_A1B1CD_ZJets, "ZJets", "Single",year);
    }

    else if (whichRegion=="singleLept"){
      makeFullBkgClosure(histos_ABCD_sum, "BkgSum", "Double",year);

      // makeABCDPlot(histos_ABCD_data, "Data", "Double",year);

      // makeABCDPlot(histos_ABCD_sum, "BkgSum", "Double",year);
      makeABCDPlot(histos_ABCD_TT, "TT", "Double",year);
      // makeABCDPlot(histos_ABCD_TT_Di, "TT_Di", "Double",year);
      // makeABCDPlot(histos_ABCD_TT_SL, "TT_SL", "Double",year);
      makeABCDPlot(histos_ABCD_WJets, "WJets", "Double",year);

      // makeABCDPlot(histos_ABCD_sum, "BkgSum", "Single",year);
      makeABCDPlot(histos_A1B1CD_TT, "TT", "Single",year);
      // makeABCDPlot(histos_A1B1CD_TT_Di, "TT_Di", "Single",year);
      // makeABCDPlot(histos_A1B1CD_TT_SL, "TT_SL", "Single",year);
      makeABCDPlot(histos_A1B1CD_WJets, "WJets", "Single",year);

    }
    else if (whichRegion=="photon"){
      makeABCDPlot(histos_ABCD_QCD, "QCD", "Double",year);
      makeABCDPlot(histos_ABCD_GJets, "GJets", "Double",year);
      makeABCDPlot(histos_A1B1CD_QCD, "QCD", "Single",year);
      makeABCDPlot(histos_A1B1CD_GJets, "GJets", "Single",year);
    }
  }

   // makeRpfPlot(histos_RPF_LeadSR_sum, "BkgSum", "J2", "Double",year);
   // makeRpfPlot(histos_RPF_LeadSB_sum, "BkgSum", "J2", "Double",year);


  if (runRPFPlots){
    // myfile.open("Rpf_deviations.txt");
    if (whichRegion=="signal"){
      //J1 Double tag region
      makeRpfPlot(histos_Rpf_J1_sum, "BkgSum", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_QCD, "QCD", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_TT, "TT", "J1", "Double",year);
      // makeRpfPlot(histos_Rpf_J1_TT_Di, "TT_Di", "J1", "Double",year);
      // makeRpfPlot(histos_Rpf_J1_TT_SL, "TT_SL", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_WJets, "WJets", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_ZJets, "ZJets", "J1", "Double",year);

      //J1 Single tag region
      makeRpfPlot(histos_Rpfsingle_J1_sum, "BkgSum", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_QCD, "QCD", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_TT, "TT", "J1", "Single",year);
      // makeRpfPlot(histos_Rpfsingle_J1_TT_Di, "TT_Di", "J1", "Single",year);
      // makeRpfPlot(histos_Rpfsingle_J1_TT_SL, "TT_SL", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_WJets, "WJets", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_ZJets, "ZJets", "J1", "Single",year);

      //J2 double tag region
      makeRpfPlot(histos_Rpf_J2_sum, "BkgSum", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_QCD, "QCD", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_TT, "TT", "J2", "Double",year);
      // makeRpfPlot(histos_Rpf_J2_TT_Di, "TT_Di", "J2", "Double",year);
      // makeRpfPlot(histos_Rpf_J2_TT_SL, "TT_SL", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_WJets, "WJets", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_ZJets, "ZJets", "J2", "Double",year);

      //J2 double tag region, with three mass bins
      makeRpfPlot(histos_Rpf_J2_mJBins_sum, "BkgSum_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_QCD, "QCD_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_TT, "TT_mJ", "J2", "Double",year);
      // makeRpfPlot(histos_Rpf_J2_mJBins_TT_Di, "TT_Di_mJ", "J2", "Double",year);
      // makeRpfPlot(histos_Rpf_J2_mJBins_TT_SL, "TT_SL_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_WJets, "WJets_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_ZJets, "ZJets_mJ", "J2", "Double",year);


      //J2 single tag region
      makeRpfPlot(histos_Rpfsingle_J2_sum, "BkgSum", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_QCD, "QCD", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_TT, "TT", "J2", "Single",year);
      // makeRpfPlot(histos_Rpfsingle_J2_TT_Di, "TT_Di", "J2", "Single",year);
      // makeRpfPlot(histos_Rpfsingle_J2_TT_SL, "TT_SL", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_WJets, "WJets", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_ZJets, "ZJets", "J2", "Single",year);
    }
    if (whichRegion=="singleLept"){
      //J1 Double tag region
      // makeRpfPlot(histos_Rpf_J1_sum, "BkgSum", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_TT, "TT", "J1", "Double",year);
      // makeRpfPlot(histos_Rpf_J1_TT_Di, "TT_Di", "J1", "Double",year);
      // makeRpfPlot(histos_Rpf_J1_TT_SL, "TT_SL", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_WJets, "WJets", "J1", "Double",year);


      //J1 Single tag region
      // makeRpfPlot(histos_Rpfsingle_J1_sum, "BkgSum", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_TT, "TT", "J1", "Single",year);
      // makeRpfPlot(histos_Rpfsingle_J1_TT_Di, "TT_Di", "J1", "Single",year);
      // makeRpfPlot(histos_Rpfsingle_J1_TT_SL, "TT_SL", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_WJets, "WJets", "J1", "Single",year);


      //J2 double tag region
      // makeRpfPlot(histos_Rpf_J2_sum, "BkgSum", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_TT, "TT", "J2", "Double",year);
      // makeRpfPlot(histos_Rpf_J2_TT_Di, "TT_Di", "J2", "Double",year);
      // makeRpfPlot(histos_Rpf_J2_TT_SL, "TT_SL", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_WJets, "WJets", "J2", "Double",year);


      //J2 single tag region
      // makeRpfPlot(histos_Rpfsingle_J2_sum, "BkgSum", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_TT, "TT", "J2", "Single",year);
      // makeRpfPlot(histos_Rpfsingle_J2_TT_Di, "TT_Di", "J2", "Single",year);
      // makeRpfPlot(histos_Rpfsingle_J2_TT_SL, "TT_SL", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_WJets, "WJets", "J2", "Single",year);

      //J2 double tag region, with three mass bins
      // makeRpfPlot(histos_Rpf_J2_mJBins_sum, "BkgSum_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_TT, "TT_mJ", "J2", "Double",year);
      // makeRpfPlot(histos_Rpf_J2_mJBins_TT_Di, "TT_Di_mJ", "J2", "Double",year);
      // makeRpfPlot(histos_Rpf_J2_mJBins_TT_SL, "TT_SL_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_WJets, "WJets_mJ", "J2", "Double",year);
    }
    if (whichRegion=="photon"){
      //J1 Double tag region
      makeRpfPlot(histos_Rpf_J1_sum, "BkgSum", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_QCD, "QCD", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_GJets, "GJets", "J1", "Double",year);

      //J1 Single tag region
      makeRpfPlot(histos_Rpfsingle_J1_sum, "BkgSum", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_QCD, "QCD", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_GJets, "GJets", "J1", "Single",year);

      //J2 double tag region
      makeRpfPlot(histos_Rpf_J2_sum, "BkgSum", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_QCD, "QCD", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_GJets, "GJets", "J2", "Double",year);

      //J2 single tag region
      makeRpfPlot(histos_Rpfsingle_J2_sum, "BkgSum", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_QCD, "QCD", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_GJets, "GJets", "J2", "Single",year);

      // makeRpfPlot(histos_RPF_LeadSR_sum, "BkgSum", "J2", "Double",year);
    }

    // myfile.close();
  }


  // if (runStacks && whichRegion=="signal") {
  //   makeStackPlot(histos_j1mass_QCD,histos_j1mass_TT,histos_j1mass_WJets,histos_j1mass_ZJets,histos_j1mass_T5HH1300,histos_j1mass_T5HH1700,histos_j1mass_T5HH2100,"leadmass");
  //   makeStackPlot(histos_j2mass_QCD,histos_j2mass_TT,histos_j2mass_WJets,histos_j2mass_ZJets,histos_j2mass_T5HH1300,histos_j2mass_T5HH1700,histos_j2mass_T5HH2100,"subleadmass");
  //   makeStackPlot(histos_j1pt_QCD,histos_j1pt_TT,histos_j1pt_WJets,histos_j1pt_ZJets,histos_j1pt_T5HH1300,histos_j1pt_T5HH1700,histos_j1pt_T5HH2100,"leadpt");
  //   makeStackPlot(histos_j2pt_QCD,histos_j2pt_TT,histos_j2pt_WJets,histos_j2pt_ZJets,histos_j2pt_T5HH1300,histos_j2pt_T5HH1700,histos_j2pt_T5HH2100,"subleadpt");
  // }
  //
  // if (runStacks && whichRegion=="singleLept") {
  //   makeSingleLeptStackPlot(histos_j1mass_TT,histos_j1mass_WJets,"leadmass");
  //   makeSingleLeptStackPlot(histos_j2mass_TT,histos_j2mass_WJets,"subleadmass");
  //   makeSingleLeptStackPlot(histos_j1pt_TT,histos_j1pt_WJets,"leadpt");
  //   makeSingleLeptStackPlot(histos_j2pt_TT,histos_j2pt_WJets,"subleadpt");
  // }
  //


  if (runMETNorm){
    makeMETNorm(h_METShape_2H, "2H");
    // makeMETNorm(h_METShape_1H, "1H");
  }

  // if (runROC){
  //   TH1F * WJetsH_Added = (TH1F*)histos_doubleBHJ1_WJets[0]->Clone("WJetsH");
  //   WJetsH_Added->Add(histos_doubleBHJ1_WJets[1]); WJetsH_Added->Add(histos_doubleBHJ1_WJets[2]);
  //
  //   TH1F * TTH_Added = (TH1F*)histos_doubleBHJ1_TT[0]->Clone("TTH");
  //   TTH_Added->Add(histos_doubleBHJ1_TT[1]); TTH_Added->Add(histos_doubleBHJ1_TT[2]);
  //
  //   TH1F * T5HH1300H_Added = (TH1F*)histos_doubleBHJ1_T5HH1300[0]->Clone("T5HH1300H");
  //   T5HH1300H_Added->Add(histos_doubleBHJ1_T5HH1300[1]); T5HH1300H_Added->Add(histos_doubleBHJ1_T5HH1300[2]);
  //
  //   TH1F * T5HH1700H_Added = (TH1F*)histos_doubleBHJ1_T5HH1700[0]->Clone("T5HH1700H");
  //   T5HH1700H_Added->Add(histos_doubleBHJ1_T5HH1700[1]); T5HH1700H_Added->Add(histos_doubleBHJ1_T5HH1700[2]);
  //
  //   TH1F * T5HH2100H_Added = (TH1F*)histos_doubleBHJ1_T5HH2100[0]->Clone("T5HH2100H");
  //   T5HH2100H_Added->Add(histos_doubleBHJ1_T5HH2100[1]); T5HH2100H_Added->Add(histos_doubleBHJ1_T5HH2100[2]);
  //
  //   QuickROC(T5HH1300H_Added, WJetsH_Added, "DoubleB-H");
  //   QuickROC(T5HH1300H_Added, TTH_Added, "DoubleB-H");
  //
  //   TH1F * WJets_Added = (TH1F*)histos_doubleBJ1_WJets[0]->Clone("WJets");
  //   WJets_Added->Add(histos_doubleBJ1_WJets[1]); WJets_Added->Add(histos_doubleBJ1_WJets[2]);
  //
  //   TH1F * TT_Added = (TH1F*)histos_doubleBJ1_TT[0]->Clone("TT");
  //   TT_Added->Add(histos_doubleBJ1_TT[1]); TT_Added->Add(histos_doubleBJ1_TT[2]);
  //
  //   TH1F * T5HH1300_Added = (TH1F*)histos_doubleBJ1_T5HH1300[0]->Clone("T5HH1300");
  //   T5HH1300_Added->Add(histos_doubleBJ1_T5HH1300[1]); T5HH1300_Added->Add(histos_doubleBJ1_T5HH1300[2]);
  //
  //   QuickROC(T5HH1300_Added, WJets_Added, "");
  //   QuickROC(T5HH1300_Added, TT_Added, "");
  // }

  fout->Close();
  f->Close();
  // fSignal->Close();
} //end method ABCD()


void makeABCDPlot(vector<TH1F*> dem_histos, TString bkgType = "", TString tagType = "", TString year = "") {
  TH1F *histo_A = (TH1F*)dem_histos[0]->Clone("histo_A"); TH1F *histo_B = (TH1F*)dem_histos[1]->Clone("histo_B");
  TH1F *histo_C = (TH1F*)dem_histos[2]->Clone("histo_C"); TH1F *histo_D = (TH1F*)dem_histos[3]->Clone("histo_D");
  TH1F *histo_pred = (TH1F*)histo_B->Clone("histo_pred");
  bool isData = false;
  TString histoName = dem_histos[0]->GetName();
  if (histoName.Contains("data")) isData=true;

  histo_pred->Multiply(histo_C);histo_pred->Divide(histo_D);

  if (!isData){
    int numMETBins = histo_pred->GetNbinsX();
    for (int i = 1; i<=numMETBins; i++){
      if (histo_A->GetBinContent(i)<0.001){
        std::cout<<"Attempting to solve problem"<<std::endl;
        histo_A->SetBinContent(i,0.001);
        if (histo_A->GetBinError(i)<0.001) histo_A->SetBinError(i,1.4);
      }
      if (histo_pred->GetBinContent(i)<0.001){
        std::cout<<"Attempting to solve problem"<<std::endl;
        histo_pred->SetBinContent(i,0.001);
        if (histo_pred->GetBinError(i)<0.001) histo_pred->SetBinError(i,1.4);
      }
    }
  }
  else {
    for (int i = 1; i<=3; i++){
      if (histo_A->GetBinContent(i)<0.001){
        std::cout<<"Attempting to solve problem"<<std::endl;
        histo_A->SetBinContent(i,0.001);
        if (histo_A->GetBinError(i)<0.001) histo_A->SetBinError(i,1.4);
      }
      if (histo_pred->GetBinContent(i)<0.001){
        std::cout<<"Attempting to solve problem"<<std::endl;
        histo_pred->SetBinContent(i,0.001);
        if (histo_pred->GetBinError(i)<0.001) histo_pred->SetBinError(i,1.4);
      }
    }
  }

  // histo_A->SetBinError(3,1.0);
  // histo_pred->SetBinError(3,1.0);
  //
  // cout<<"Will work!"<<endl;
  // cout<<"A2: "<<histo_A->GetBinContent(1)<<"+/-"<<histo_A->GetBinError(1);
  // cout<<", "<<histo_A->GetBinContent(2)<<"+/-"<<histo_A->GetBinError(2);
  // cout<<", "<<histo_A->GetBinContent(3)<<"+/-"<<histo_A->GetBinError(3)<<endl;
  //
  // cout<<"Pred: "<<histo_pred->GetBinContent(1)<<"+/-"<<histo_pred->GetBinError(1);
  // cout<<", "<<histo_pred->GetBinContent(2)<<"+/-"<<histo_pred->GetBinError(2);
  // cout<<", "<<histo_pred->GetBinContent(3)<<"+/-"<<histo_pred->GetBinError(3)<<endl;

  //apply kappa correction
  // if (!isData) histo_pred->Scale(0.66);

  TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo_A, histo_pred, "pois");
  TString canvName = bkgType+tagType;
  double W = 800; double H = 600;
  double T = 0.08*H; double B = 0.12*H;
  double L = 0.12*W; double R = 0.04*W;
  TCanvas * can_h = new TCanvas(canvName,canvName, 50, 50, W, H);
  can_h->SetFillColor(0); can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W ); can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H ); can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0); can_h->SetTicky(0);

  double up_height = 0.8; // double dw_correction = 1.30;
  double dw_correction = 1.18; double font_size_dw  = 0.1;
  double dw_height = (1. - up_height) * dw_correction;
  double dw_height_offset = 0.02;

  TPad *pad1 = new TPad("pad1", "top pad" , 0.0, 0.25, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.25);
  pad1->SetTickx(0); pad1->SetTicky(0);
  pad1->SetPad(0., 1 - up_height, 1., 1.00);
  pad1->SetFrameFillColor(0); pad1->SetFillColor(0);
  pad1->SetTopMargin(0.08); pad1->SetLeftMargin(0.08); pad1->SetRightMargin(0.06);
  pad1->Draw();

  pad2->SetPad(0., 0., 1., dw_height+dw_height_offset);
  pad2->SetFillColor(0); pad2->SetFrameFillColor(0);
  pad2->SetBottomMargin(0.33); pad2->SetTopMargin(0.01);
  pad2->SetLeftMargin(0.08); pad2->SetRightMargin(0.06);
  pad2->Draw(); pad1->cd();

  histo_pred->SetStats(0); histo_A->SetStats(0);
  histo_pred->SetFillColor(kBlue); histo_pred->SetFillStyle(3445);

  histo_pred->SetMarkerSize(0); histo_pred->SetLineWidth(1);
  histo_pred->SetLineColor(kBlue);

  TString title = bkgType+"_"+year+";;nEvents";
  histo_A->SetTitle(title); histo_pred->SetTitle(title);
  histo_pred->SetMinimum(0); histo_A->SetMinimum(0);

  histo_A->SetLineColor(kBlack); histo_A->SetMarkerStyle(20);
  histo_A->SetMarkerColor(kBlack);

  TLegend* legend = new TLegend(0.60,0.7,0.93,0.91);
  legend->SetBorderSize(0);
  if (!isData){
    if (tagType=="Double") {
      legend->AddEntry(histo_A,"MC: A2","lp");
      legend->AddEntry(histo_pred,"Pred: B2*C/D","f");
    }
    else if (tagType=="Single") {
      legend->AddEntry(histo_A,"MC: A1","lp");
      legend->AddEntry(histo_pred,"Pred: B1*C/D","f");
    }
  }

  else {
    if (tagType=="Double") {
      legend->AddEntry(histo_A,"Data A2","lp");
      legend->AddEntry(histo_pred,"Pred: B2*C/D","f");
    }
    else if (tagType=="Single") {
      legend->AddEntry(histo_A,"Data A1","lp");
      legend->AddEntry(histo_pred,"Pred: B1*C/D","f");
    }
  }

  float pred_bin1 = histo_pred->GetBinContent(1);
  float sig_bin1 = histo_A->GetBinContent(1);

  if (pred_bin1>sig_bin1) {
    histo_pred->GetYaxis()->SetTitleOffset(0.6);
    histo_pred->GetXaxis()->SetLabelSize(0);
    histo_pred->Draw("E2");
    histo_A->Draw("same");
  }
  else {
    histo_A->GetYaxis()->SetTitleOffset(0.6);
    histo_A->GetXaxis()->SetLabelSize(0);
    histo_A->Draw("same");
    histo_pred->Draw("E2 same");
  }
  legend->Draw("same");

  pad2->cd();
  graph->SetTitle(";MET [GeV]; #kappa (Sim/Pred)");
  graph->SetMarkerStyle(20);
  graph->SetLineColor(kBlack); graph->SetMarkerColor(kBlack);
  graph->GetYaxis()->SetTitleOffset(0.25);
  graph->GetYaxis()->SetTitleSize(0.13);
  graph->GetXaxis()->SetTitleOffset(0.9);
  graph->GetXaxis()->SetTitleSize(0.14);
  graph->GetXaxis()->SetLabelSize(0.14);
  graph->GetYaxis()->SetLabelSize(0.14);
  graph->GetXaxis()->SetNdivisions(507);
  graph->GetYaxis()->SetNdivisions(505);
  // graph->GetYaxis()->SetRangeUser(0.3,1.5);

  graph->Draw("APE");
  graph->GetXaxis()->SetRangeUser(300,1400);
  graph->GetYaxis()->SetRangeUser(-1.0,5.0);
  if (tagType=="Single") graph->GetYaxis()->SetRangeUser(0.0,2.0);
  graph->Draw("APE");
  can_h->Modified(); can_h->Update();
  graph->GetXaxis()->SetRangeUser(300,1400);
  can_h->Modified(); can_h->Update();


  TLine *line = new TLine(300,1.0,1400,1.0);
  line->SetLineColor(kRed); line->SetLineStyle(2);
  line->Draw("same");

  TString savename = "ABCD_"+tagType+"_"+bkgType+"_"+year;
  cdABCD->cd();
  can_h->Write(savename);
}

void makeFullBkgClosure(vector<TH1F*> dem_histos, TString bkgType = "", TString tagType = "", TString year = "") {
  TH1F *histo_A = (TH1F*)dem_histos[0]->Clone("histo_A"); TH1F *histo_B = (TH1F*)dem_histos[1]->Clone("histo_B");
  TH1F *histo_C = (TH1F*)dem_histos[2]->Clone("histo_C"); TH1F *histo_D = (TH1F*)dem_histos[3]->Clone("histo_D");
  TH1F *histo_pred = (TH1F*)histo_B->Clone("histo_pred");

  TH1F *h_finalBkg = (TH1F*)histo_pred->Clone("h_finalBkg");
  histo_pred->Multiply(histo_C);histo_pred->Divide(histo_D);

  float bkgNorm = histo_pred->Integral();
  // float kappa = 0.66; //SL
  // float bkgFrac1 = 0.88; //SL
  // float bkgFrac2 = 0.098; //SL
  // float bkgFrac3 = 0.019; //SL

  float kappa = 1.1; //Signal
  float bkgFrac1 = 0.86; //Signal
  float bkgFrac2 = 0.118; //Signal
  float bkgFrac3 = 0.0251; //Signal

  h_finalBkg->SetBinContent(1,bkgFrac1*bkgNorm*kappa);
  h_finalBkg->SetBinContent(2,bkgFrac2*bkgNorm*kappa);
  h_finalBkg->SetBinContent(3,bkgFrac3*bkgNorm*kappa);

  // int numMETBins = h_finalBkg->GetNbinsX();
  // for (int i = 1; i<=numMETBins; i++){
  //     h_finalBkg->SetBinContent(i,0.001);
  //     if (histo_pred->GetBinError(i)<0.001) histo_pred->SetBinError(i,1.4);
  // }
  // histo_A->SetBinError(3,1.0);
  // histo_pred->SetBinError(3,1.0);

  TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo_A, h_finalBkg, "pois");
  TString canvName = bkgType+tagType;
  double W = 800; double H = 600;
  double T = 0.08*H; double B = 0.12*H;
  double L = 0.12*W; double R = 0.04*W;
  TCanvas * can_h = new TCanvas(canvName,canvName, 50, 50, W, H);
  can_h->SetFillColor(0); can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W ); can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H ); can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0); can_h->SetTicky(0);

  double up_height = 0.8; // double dw_correction = 1.30;
  double dw_correction = 1.18; double font_size_dw  = 0.1;
  double dw_height = (1. - up_height) * dw_correction;
  double dw_height_offset = 0.02;

  TPad *pad1 = new TPad("pad1", "top pad" , 0.0, 0.25, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.25);
  pad1->SetTickx(0); pad1->SetTicky(0);
  pad1->SetPad(0., 1 - up_height, 1., 1.00);
  pad1->SetFrameFillColor(0); pad1->SetFillColor(0);
  pad1->SetTopMargin(0.08); pad1->SetLeftMargin(0.08); pad1->SetRightMargin(0.06);
  pad1->Draw();

  pad2->SetPad(0., 0., 1., dw_height+dw_height_offset);
  pad2->SetFillColor(0); pad2->SetFrameFillColor(0);
  pad2->SetBottomMargin(0.33); pad2->SetTopMargin(0.01);
  pad2->SetLeftMargin(0.08); pad2->SetRightMargin(0.06);
  pad2->Draw(); pad1->cd();

  h_finalBkg->SetStats(0); histo_A->SetStats(0);
  h_finalBkg->SetFillColor(kBlue); h_finalBkg->SetFillStyle(3445);

  h_finalBkg->SetMarkerSize(0); histo_pred->SetLineWidth(1);
  h_finalBkg->SetLineColor(kBlue);

  TString title = bkgType+"_"+year+";;nEvents";
  histo_A->SetTitle(title); h_finalBkg->SetTitle(title);
  h_finalBkg->SetMinimum(0); histo_A->SetMinimum(0);

  histo_A->SetLineColor(kBlack); histo_A->SetMarkerStyle(20);
  histo_A->SetMarkerColor(kBlack);

  TLegend* legend = new TLegend(0.60,0.7,0.93,0.91);
  legend->SetBorderSize(0);
  if (tagType=="Double") {
    legend->AddEntry(histo_A,"MC: A2","lp");
    legend->AddEntry(h_finalBkg,"Pred: (B2*C/D)*kappa*bkgFrac","f");
  }
  else if (tagType=="Single") {
    legend->AddEntry(histo_A,"MC: A1","lp");
    legend->AddEntry(h_finalBkg,"Pred: (B1*C/D)*kappa*bkgFrac","f");
  }

  float pred_bin1 = h_finalBkg->GetBinContent(1);
  float sig_bin1 = histo_A->GetBinContent(1);

  if (pred_bin1>sig_bin1) {
    h_finalBkg->GetYaxis()->SetTitleOffset(0.6);
    h_finalBkg->GetXaxis()->SetLabelSize(0);
    h_finalBkg->Draw("E2");
    histo_A->Draw("same");
  }
  else {
    histo_A->GetYaxis()->SetTitleOffset(0.6);
    histo_A->GetXaxis()->SetLabelSize(0);
    histo_A->Draw("same");
    h_finalBkg->Draw("E2 same");
  }
  legend->Draw("same");

  pad2->cd();
  graph->SetTitle(";MET [GeV]; MC/Pred");
  graph->SetMarkerStyle(20);
  graph->SetLineColor(kBlack); graph->SetMarkerColor(kBlack);
  graph->GetYaxis()->SetTitleOffset(0.25);
  graph->GetYaxis()->SetTitleSize(0.13);
  graph->GetXaxis()->SetTitleOffset(0.9);
  graph->GetXaxis()->SetTitleSize(0.14);
  graph->GetXaxis()->SetLabelSize(0.14);
  graph->GetYaxis()->SetLabelSize(0.14);
  graph->GetXaxis()->SetNdivisions(507);
  graph->GetYaxis()->SetNdivisions(505);
  // graph->GetYaxis()->SetRangeUser(0.3,1.5);

  graph->Draw("APE");
  graph->GetXaxis()->SetRangeUser(300,1400);
  graph->GetYaxis()->SetRangeUser(-1.0,5.0);
  if (tagType=="Single") graph->GetYaxis()->SetRangeUser(0.0,2.0);
  graph->Draw("APE");
  can_h->Modified(); can_h->Update();
  graph->GetXaxis()->SetRangeUser(300,1400);
  can_h->Modified(); can_h->Update();


  TLine *line = new TLine(300,1.0,1400,1.0);
  line->SetLineColor(kRed); line->SetLineStyle(2);
  line->Draw("same");

  cdClose->cd();
  can_h->Write("ABCD_Closure");
}

void makeClosureStackPlot(vector<TH1F*> QCD_histos, vector<TH1F*> TT_histos, vector<TH1F*> WJets_histos,vector<TH1F*> ZJets_histos, vector<TH1F*> sum_histos){
  TH1F *QCD_histo_A = (TH1F*)QCD_histos[0]->Clone("QCD_histo_A");
  TH1F *TT_histo_A = (TH1F*)TT_histos[0]->Clone("TT_histo_A");
  TH1F *WJets_histo_A = (TH1F*)WJets_histos[0]->Clone("WJets_histo_A");
  TH1F *ZJets_histo_A = (TH1F*)ZJets_histos[0]->Clone("ZJets_histo_A");

  TH1F *sum_B = (TH1F*)sum_histos[1]->Clone("sum_B");
  TH1F *sum_C = (TH1F*)sum_histos[2]->Clone("sum_C");
  TH1F *sum_D = (TH1F*)sum_histos[3]->Clone("sum_D");
  sum_B->Multiply(sum_C); sum_B->Divide(sum_D);
  TH1F *histo_pred = (TH1F*)sum_B->Clone("histo_pred");

  TH1F *histo_MC  = (TH1F*)QCD_histo_A->Clone("histo_MC");
  histo_MC->Add(TT_histo_A); histo_MC->Add(WJets_histo_A); histo_MC->Add(ZJets_histo_A);

  histo_pred->SetFillColor(kRed);
  histo_pred->SetFillStyle(3445);
  QCD_histo_A->SetLineColor(kGray);QCD_histo_A->SetFillColor(kGray);
  TT_histo_A->SetLineColor(kCyan+1);TT_histo_A->SetFillColor(kCyan+1);
  WJets_histo_A->SetLineColor(kBlue);WJets_histo_A->SetFillColor(kBlue);
  ZJets_histo_A->SetLineColor(kGreen+2);ZJets_histo_A->SetFillColor(kGreen+2);

  histo_MC->SetMarkerStyle(20);
  histo_MC->SetLineColor(kBlack);histo_MC->SetMarkerColor(kBlack);

  TLegend* legend = new TLegend(0.55,0.65,0.85,0.85,NULL,"brNDC");
  legend->SetBorderSize(0);
  legend->SetNColumns(2);
  legend->AddEntry(QCD_histo_A, "QCD", "F");
  legend->AddEntry(TT_histo_A, "TT", "F");
  legend->AddEntry(WJets_histo_A, "WJets", "F");
  legend->AddEntry(ZJets_histo_A, "ZJets", "F");
  legend->AddEntry(histo_pred, "Prediction", "F");
  legend->AddEntry(histo_MC, "MC Truth", "PL");

  THStack* stack_closure_MC = new THStack("stack_closure_MC","stack_closure_MC");
  stack_closure_MC->Add(QCD_histo_A);
  stack_closure_MC->Add(TT_histo_A);
  stack_closure_MC->Add(WJets_histo_A);
  stack_closure_MC->Add(ZJets_histo_A);

  TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo_MC, histo_pred, "pois");
  graph->GetXaxis()->SetTitle("MET [GeV]");
  graph->GetYaxis()->SetTitle("#kappa (Sim/Pred)");
  graph->SetTitle("");
  graph->GetYaxis()->SetTitleSize(0.055);
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.85);
  // graph->SetMinimum(0.0);
  graph->SetLineColor(kBlack);
  graph->SetMarkerColor(kBlack);
  // graph->GetXaxis()->SetRangeUser(50,250);

  TCanvas *c1 = new TCanvas("c1","MC Closure",10,10,1000,800);
  double up_height     = 0.8; double dw_correction = 1.18;
  double font_size_dw  = 0.1;
  double dw_height = (1. - up_height) * dw_correction;
  double dw_height_offset = 0.02;

  TPad *pad1 = new TPad("pad1", "top pad" , 0.0, 0.25, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.25);
  pad1->SetTickx(0); pad1->SetTicky(0);
  pad1->SetFrameFillColor(0); pad1->SetFillColor(0);
  pad1->SetTopMargin(0.08); pad1->SetLeftMargin(0.10);
  pad1->SetRightMargin(0.04);
  pad1->Draw();

  pad2->SetFillColor(0); pad2->SetFrameFillColor(0);
  pad2->SetBottomMargin(0.33); pad2->SetTopMargin(0.03);
  pad2->SetLeftMargin(0.10); pad2->SetRightMargin(0.04);
  pad2->Draw(); pad1->cd();

  histo_pred->SetTitle(";;Events (137 fb^{-1})");
  histo_pred->GetYaxis()->SetTitleOffset(0.7);
  histo_pred->Draw("E2");
  stack_closure_MC->Draw("hist same");
  histo_pred->Draw("E2 same");
  histo_MC->Draw("same");
  legend->Draw("same");

  pad2->cd();
  graph->GetYaxis()->SetTitleOffset(0.2);
  graph->GetYaxis()->SetTitleSize(0.13);
  graph->GetXaxis()->SetTitleOffset(0.9);
  graph->GetXaxis()->SetTitleSize(0.14);
  graph->GetXaxis()->SetLabelSize(0.14);
  graph->GetYaxis()->SetLabelSize(0.10);
  // graph->GetXaxis()->SetNdivisions(507);
  graph->GetYaxis()->SetNdivisions(505);
  graph->Draw("APE");
  graph->GetXaxis()->SetRangeUser(250,1400);
  graph->Draw("APE");
  TLine *line = new TLine(250,1.0,1400,1.0);
  line->SetLineColor(kRed); line->SetLineStyle(2);
  line->Draw("same");
  c1->Modified();c1->Update();

  cdClose->cd();
  c1->Write("ABCD_Closure");
}

void makeRpfPlot(vector<TH1F*> dem_histos, TString bkgType = "", TString jetType = "", TString tagType = "", TString year = "") {
  TH1D * num   = (TH1D*)dem_histos[1]->Clone("num"+tagType+"_"+jetType+"_"+bkgType); //this should be 2H SB
  TH1D * denom = (TH1D*)dem_histos[3]->Clone("denom"+tagType+"_"+jetType+"_"+bkgType);  //this should be 0H SB
  // num->Add(dem_histos[1]); denom->Add(dem_histos[3]);
  // num->Rebin(8);

  if (!bkgType.Contains("_mJ")){
    num->Rebin(2);
    denom->Rebin(2);
  }


  num->SetMarkerStyle(20); num->SetMarkerSize(0.85);
  num->SetLineColor(kRed); num->SetMarkerColor(kRed);
  denom->SetMarkerStyle(20); denom->SetMarkerSize(0.85);
  denom->SetLineColor(kBlue); denom->SetMarkerColor(kBlue);
  num->SetMinimum(0.0);denom->SetMinimum(0.0);

  TString thisType = tagType+"_"+jetType+"_"+bkgType+"_"+year;
  num->SetTitle(thisType+";;Normalized to area");
  num->GetYaxis()->SetTitleOffset(0.77);

  TString graphName = "ratio"+thisType;
  // num->Sumw2(); denom->Sumw2(); //already done
  TGraphAsymmErrors * graph = new TGraphAsymmErrors(num, denom, "pois");
  graph->SetName(graphName);
  graph->SetMarkerStyle(20); graph->SetMarkerSize(0.85);
  graph->SetMinimum(-0.1); graph->GetXaxis()->SetRangeUser(50,250);
  graph->SetLineColor(kBlack); graph->SetMarkerColor(kBlack);

  double W = 800; double H = 600;
  double T = 0.08*H; double B = 0.12*H;
  double L = 0.15*W; double R = 0.02*W;
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0); can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W ); can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H ); can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0); can_h->SetTicky(0);

  double up_height = 0.8; double dw_correction = 1.18;
  double font_size_dw  = 0.1; double dw_height_offset = 0.02;
  double dw_height = (1. - up_height) * dw_correction;

  TPad *pad1 = new TPad("pad1", "top pad" , 0.0, 0.25, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.25);
  pad1->SetTickx(0); pad1->SetTicky(0);
  pad1->SetPad(0., 1 - up_height, 1., 1.00);
  pad1->SetFrameFillColor(0); pad1->SetFillColor(0);
  pad1->SetTopMargin(0.08); pad1->SetLeftMargin(0.10); pad1->SetRightMargin(0.04);
  pad1->Draw();

  pad2->SetPad(0., 0., 1., dw_height+dw_height_offset);
  pad2->SetFillColor(0); pad2->SetFrameFillColor(0);
  pad2->SetBottomMargin(0.33); pad2->SetTopMargin(0.01);
  pad2->SetLeftMargin(0.10); pad2->SetRightMargin(0.04);
  pad2->Draw(); pad1->cd();

  num->GetXaxis()->SetLabelSize(0);
  // num->Rebin(2);
  // denom->Rebin(2);
  num->DrawNormalized();
  denom->DrawNormalized("same");

  TLegend* legend = new TLegend(0.5,0.65,0.95,0.88,NULL,"brNDC") ;
  legend->SetBorderSize(0);
  // legend->SetHeader(year);
  if (tagType == "Single"){
    legend->AddEntry(num, "1H SB (numerator)", "PL");
    legend->AddEntry(denom, "0H SB (denominator)", "PL");
    graph->SetTitle(";Soft-drop mass [GeV];R(1H/0H)");
  }
  else if (tagType == "Double"){
    legend->AddEntry(num, "2H SB (numerator)", "PL");
    legend->AddEntry(denom, "0H SB (denominator)", "PL");
    graph->SetTitle(";Soft-drop mass [GeV];R(2H/0H)");
  }

  legend->Draw();

  pad2->cd();
  graph->Draw("APE");
  graph->GetYaxis()->SetTitleOffset(0.25); graph->GetYaxis()->SetTitleSize(0.13);
  graph->GetXaxis()->SetTitleOffset(0.9); graph->GetXaxis()->SetTitleSize(0.14);
  graph->GetXaxis()->SetLabelSize(0.14); graph->GetYaxis()->SetLabelSize(0.10);
  // graph->GetXaxis()->SetNdivisions(507);
  graph->GetYaxis()->SetNdivisions(505);


  // graph->GetYaxis()->SetTitleSize(0.055);
  graph->Draw("APE");
  graph->GetYaxis()->SetRangeUser(-0.1,0.5);

  if (tagType=="Double"){
    if (bkgType=="WJets" || bkgType=="ZJets") graph->GetYaxis()->SetRangeUser(-0.005,0.015);
    else if (bkgType=="BkgSum") graph->GetYaxis()->SetRangeUser(-0.04,0.09);
  }

  else if (tagType=="Single"){
    if (bkgType=="TT") graph->GetYaxis()->SetRangeUser(0,2.5);
    else if (bkgType=="WJets" || bkgType=="ZJets") graph->GetYaxis()->SetRangeUser(0.07,0.15);
  }

  // graph->Draw("APE");
  can_h->Modified();
  can_h->Update();

  TF1*f0=new TF1("f0","pol0", 50,150);
  graph->Fit("f0","Q","R",50,150);
  double p0 = f0->GetParameter(0); //this does, indeed, work
  double error = f0->GetParError(0);
  graph->GetFunction("f0")->SetBit(TF1::kNotDraw);

  TLine *line = new TLine(50.0,p0,250.0,p0);
  line->SetLineColor(kRed);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw("same");
  std::string lineConst = to_string(p0); std::string lineErr = to_string(error);
  TString constString = "Const: "+lineConst+" #pm "+ lineErr;
  TLatex *t = new TLatex(0.2,0.85,constString);
  t->SetNDC(); t->SetTextFont(52);
  t->SetTextSize(0.13);
  t->Draw();

  // TF1*f0=new TF1("f0","pol0", 50,250);
  // // graph->Fit("f0","Q","R",50,250);
  // graph->Fit("f0","Q","R",50,180);
  // gStyle->SetOptFit();

  // if (getDeviation){
  //   TF1*f0=new TF1("f0","pol0", 50,250);
  //   graph->Fit("f0","Q","R",50,180);
  //   //TF1 * thisFit = graph->GetFunction("pol0");
  //   double p0 = f0->GetParameter(0); //this does, indeed, work
  //   graph->GetFunction("f0")->SetBit(TF1::kNotDraw);
  //
  //   TLine *line = new TLine(50.0,p0,250.0,p0);
  //   line->SetLineColor(kRed);
  //   line->SetLineWidth(3);
  //   //line->SetLineStyle(2);
  //   line->Draw("same");
  //
  //
  //   int nbinsx = graph->GetN();
  //   double ax[nbinsx],ay[nbinsx];
  //   vector<double> deviation ={0};
  //   for (int i=0; i< nbinsx; i++){
  //     graph->GetPoint(i,ax[i],ay[i]);
  //     const double x0 = ax[i]; const double y0 = ay[i];
  //     const double errorY = graph->GetErrorYlow(i);
  //     double dev = (y0-p0)/errorY;
  //     deviation.push_back(dev);
  //   }
  //   myfile<<thisType<<": ";
  //   for (int i=0;i< deviation.size();i++) myfile<<deviation[i]<<", ";
  //   myfile<<"p0: "<<p0<<endl<<endl;
  //
  // }
    cdRPFFinal->cd();
    can_h->Write(thisType);

    cdRPF->cd();
    num->Write("Num_"+thisType);
    denom->Write("Denom_"+thisType);
    graph->Write("RPFRatio_"+thisType);
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

  THStack * doubleBStack = new THStack("hs","");
  doubleBStack->Add(h_ZJets_stack);
  doubleBStack->Add(h_WJets_stack);
  doubleBStack->Add(h_TT_stack);
  doubleBStack->Add(h_QCD_stack);

  double W = 800; double H = 600;
  double T = 0.08*H; double B = 0.12*H;
  double L = 0.08*W; double R = 0.08*W;
  TString graphName = "";
  if (which == "doubleB") graphName = "DoubleB_Stack";
  else if (which == "tau32") graphName = "Tau32_Stack";
  else if (which == "leadmass") graphName = "LeadMass_Stack";
  else if (which == "subleadmass") graphName = "SubLeadMass_Stack";
  else if (which == "leadpt") graphName = "LeadPt_Stack";
  else if (which == "subleadpt") graphName = "SubLeadPt_Stack";

  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0); can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W ); can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H ); can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0); can_h->SetTicky(0);
  can_h->SetLogy();

  doubleBStack->Draw("hist");
  doubleBStack->GetYaxis()->SetTitle("Events");
  doubleBStack->SetMinimum(0.1); doubleBStack->SetMaximum(300);

  if (which == "DoubleB") doubleBStack->GetXaxis()->SetTitle("double-B discriminator");
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
  can_h->Update();can_h->Modified();

  TString savename = which+"_Stack";
  cdOther->cd();
  can_h->Write(savename);
}


void makeSingleLeptStackPlot(vector<TH1F*> h_TT,vector<TH1F*> h_WJets, TString which = "") {
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

  double W = 800; double H = 600;
  double T = 0.08*H; double B = 0.12*H;
  double L = 0.08*W; double R = 0.08*W;

  TString graphName = which+"_Stack";
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0); can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W ); can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H ); can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0); can_h->SetTicky(0);
  can_h->SetLogy();

  doubleBStack->Draw("hist");
  doubleBStack->GetYaxis()->SetTitle("Events");
  doubleBStack->SetMinimum(0.1); doubleBStack->SetMaximum(300);

  if (which == "DoubleB") doubleBStack->GetXaxis()->SetTitle("double-B discriminator");
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
  can_h->Update(); can_h->Modified();

  TString savename = which+"_SingleLept_Stack";
  cdOther->cd();
  can_h->Write(savename);
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

  double W = 800; double H = 600;
  double T = 0.08*H; double B = 0.12*H;
  double L = 0.08*W; double R = 0.08*W;

  TString graphName = "METStack_"+which;
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0); can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W ); can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H ); can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0); can_h->SetTicky(0);
  can_h->SetLogy();

  METStack->Draw("hist");
  METStack->SetTitle(";MET [GeV];Events");
  METStack->SetMinimum(0.1); METStack->SetMaximum(300);

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
  can_h->Update(); can_h->Modified();

  TString savename = "MET"+which+"_Stack";
  cdOther->cd();
  can_h->Write(savename);
}


void QuickROC(TH1F* signalHist, TH1F* bkgHist, TString which = ""){
  TString graphName = which;
  double W = 800; double H = 600;
  double T = 0.08*H; double B = 0.12*H;
  double L = 0.08*W; double R = 0.08*W;
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);

  can_h->SetFillColor(0); can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W ); can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H ); can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0); can_h->SetTicky(0);

  signalHist->SetLineColor(kRed);
  TLegend* leg2=new TLegend(0.6,0.65,0.88,0.85);
  leg2->AddEntry(signalHist, "T5HH", "F");
  leg2->AddEntry(bkgHist, "TT", "F");

  int numBins = signalHist->GetNbinsX();
  float totalNumSig = signalHist->Integral();
  float totalNumBkg = bkgHist->Integral();

  double x[100], y[100];
  for(int i=1; i<=numBins; ++i){
    int cutbin=bkgHist->FindBin(bkgHist->GetBinLowEdge(i));
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
  can_h->Modified(); can_h->Update();

  TString savename = "ROC_doubleB";
  cdOther->cd();
  can_h->Write(savename);
}

void makeMETNorm(vector<TH1F*> dem_histos, TString tagType = "") {
  //dem_histos will be [SR,0HSB,0HSBCuts,0HSR,0HSRCuts]
  TH1D * h_SR   = (TH1D*)dem_histos[0]->Clone("h_SR");
  // TH1D * h_0HSB_noCuts = (TH1D*)dem_histos[1]->Clone("h_0HSB_noCuts");
  TH1D * h_0HSB_withCuts = (TH1D*)dem_histos[2]->Clone("h_0HSB_withCuts");
  // TH1D * h_0HSR_noCuts = (TH1D*)dem_histos[3]->Clone("h_0HSR_noCuts");
  TH1D * h_0HSR_withCuts = (TH1D*)dem_histos[4]->Clone("h_0HSR_withCuts");

  // TH1D * h_SR_TT   = (TH1D*)dem_histos[0]->Clone("h_SR_TT"); //this will be 2HSR or 1HSR
  // TH1D * h_SR_WJets   = (TH1D*)dem_histos[1]->Clone("h_SR_WJets");
  // TH1D * h_0HSR_TT= (TH1D*)dem_histos[2]->Clone("h_0HSR_TT");
  // TH1D * h_0HSR_WJets= (TH1D*)dem_histos[3]->Clone("h_0HSR_WJets");
  // TH1D * h_0HSB_TT = (TH1D*)dem_histos[4]->Clone("h_0HSB_TT");
  // TH1D * h_0HSB_WJets = (TH1D*)dem_histos[5]->Clone("h_0HSB_WJets");

  //Adding SR and SB together
  // h_0HSB_noCuts->Add(h_0HSR_noCuts);
  h_0HSB_withCuts->Add(h_0HSR_withCuts);


  //try normalizing to 1?
  h_SR->Scale(1/h_SR->Integral());
  h_0HSB_withCuts->Scale(1/h_0HSB_withCuts->Integral());

  h_SR->SetMarkerStyle(20); h_SR->SetMarkerSize(0.85);
  h_SR->SetLineColor(kRed); h_SR->SetMarkerColor(kRed);


  h_0HSB_withCuts->SetMarkerStyle(20); h_0HSB_withCuts->SetMarkerSize(0.85);
  h_0HSB_withCuts->SetLineColor(kBlack); h_0HSB_withCuts->SetMarkerColor(kBlack);

  h_SR->SetTitle("BTagsM>1;MET [GeV];Normalized to area");

  TString graphName = "MET";
  TGraphAsymmErrors * graph = new TGraphAsymmErrors(h_SR, h_0HSB_withCuts, "pois");
  graph->SetName(graphName);
  graph->SetMarkerStyle(20); graph->SetMarkerSize(0.85);
  graph->SetMinimum(-0.1); graph->GetXaxis()->SetRangeUser(300,1400);
  graph->SetLineColor(kBlack); graph->SetMarkerColor(kBlack);
  graph->SetTitle(";MET [GeV];2HSR / 0HCuts");

  double W = 800; double H = 600;
  double T = 0.08*H; double B = 0.12*H;
  double L = 0.12*W; double R = 0.06*W;
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0); can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W ); can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H ); can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0); can_h->SetTicky(0);

  double up_height = 0.8; double dw_correction = 1.18;
  double font_size_dw  = 0.1; double dw_height_offset = 0.02;
  double dw_height = (1. - up_height) * dw_correction;
  TPad *pad1 = new TPad("pad1", "top pad" , 0.0, 0.25, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.25);
  pad1->SetTickx(0); pad1->SetTicky(0);
  pad1->SetPad(0., 1 - up_height, 1., 1.00);
  pad1->SetFrameFillColor(0); pad1->SetFillColor(0);
  pad1->SetTopMargin(0.08); pad1->SetLeftMargin(0.10); pad1->SetRightMargin(0.04);
  pad1->Draw();

  pad2->SetPad(0., 0., 1., dw_height+dw_height_offset);
  pad2->SetFillColor(0); pad2->SetFrameFillColor(0);
  pad2->SetBottomMargin(0.33); pad2->SetTopMargin(0.01);
  pad2->SetLeftMargin(0.10); pad2->SetRightMargin(0.04);
  pad2->Draw(); pad1->cd();
  h_SR->GetXaxis()->SetLabelSize(0);

  h_SR->Draw();
  h_0HSB_withCuts->Draw("same");

  TLegend* legend = new TLegend(0.5,0.65,0.95,0.88,NULL,"brNDC") ;
  legend->SetBorderSize(0);
  TString legendEntryName = tagType+" SR";
  legend->AddEntry(h_SR, legendEntryName, "PL");
  legend->AddEntry(h_0HSB_withCuts, "0H SB+SR,BTagsM>1", "PL");
  legend->Draw();

  pad2->cd();
  graph->Draw("APE");
  graph->GetYaxis()->SetTitleOffset(0.25); graph->GetYaxis()->SetTitleSize(0.13);
  graph->GetXaxis()->SetTitleOffset(0.9); graph->GetXaxis()->SetTitleSize(0.14);
  graph->GetXaxis()->SetLabelSize(0.14); graph->GetYaxis()->SetLabelSize(0.10);
  graph->GetYaxis()->SetNdivisions(505);
  graph->Draw("APE");
  double *this_y = graph->GetY();
  cout<<"MET bin1: "<<this_y[0];
  cout<<", MET bin2: "<<this_y[1];
  cout<<", MET bin3: "<<this_y[2]<<endl;


  TF1*f0=new TF1("f0","pol0", 300,1400);
  graph->Fit("f0","Q","R",300,1400);
  double p0 = f0->GetParameter(0); //this does, indeed, work
  double error = f0->GetParError(0);
  graph->GetFunction("f0")->SetBit(TF1::kNotDraw);

  TLine *line = new TLine(300.0,p0,1400.0,p0);
  line->SetLineColor(kRed);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw("same");
  std::string lineConst = to_string(p0); std::string lineErr = to_string(error);
  TString constString = "Const: "+lineConst+" #pm "+ lineErr;
  TLatex *t = new TLatex(0.2,0.85,constString);
  t->SetNDC(); t->SetTextFont(52);
  t->SetTextSize(0.13);
  t->Draw();

  cdOther->cd();
  can_h->Write(tagType+"MET");
}
