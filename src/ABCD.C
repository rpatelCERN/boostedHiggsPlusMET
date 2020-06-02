#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPie.h>
#include <TF1.h>


void makeABCDPlot(vector<TH1F*> dem_histos, TString type, TString tagType, TString year);
void makeRpfPlot(vector<TH1F*> dem_histos, TString bkgType, TString jetType, TString tagType, TString year);
void makeStackPlot(vector<TH1F*> h_QCD,vector<TH1F*> h_TT,vector<TH1F*> h_WJets,vector<TH1F*> h_ZJets, vector<TH1F*> h_T5HH1300,vector<TH1F*> h_T5HH1700,vector<TH1F*> h_T5HH2100,TString which);
void makeSingleLeptStackPlot(vector<TH1F*> h_TT,vector<TH1F*> h_WJets,TString which);
void makeMETStack(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets, TH1F* h_T5HH1300,TH1F* h_T5HH1700,TH1F* h_T5HH2100,TString which);
void makeFullBkgClosure(vector<TH1F*> dem_histos, TString bkgType, TString tagType, TString year);
void makeMETNorm(vector<TH1F*> dem_histos, TString tagType);
void makeMETNormCompare(vector<TH1F*> dem_histos, TString tagType);
TH1F* make0EventUncSum(vector<TH1F*> dem_histos); //in this order: QCD, TT, WJets,ZJets
TH1F* make0EventUncSum_1l(vector<TH1F*> dem_histos); //in this order: TT, WJets
void styleCanvas(TCanvas *can_h);
void styleGraph(TGraph *graph);
void tableOfYields(vector<TH1F*> dem_histos, TString bkgType, TString tagType); //A, B, C, D
void tableOfMETNorm(vector<TH1F*> dem_histos, TString bkgType, TString tagType); //A, D, DOpt, C, COpt
void pieChart(vector<TH1F*> h_QCD, vector<TH1F*> h_WJets, vector<TH1F*> h_ZJets, vector<TH1F*> h_TT, TString regionLabel);
void pieChart1l(vector<TH1F*> h_WJets, vector<TH1F*> h_TT, TString regionLabel);


//for style
double W = 800; double H = 600;
double T = 0.08*H; double B = 0.12*H;
double L = 0.12*W; double R = 0.04*W;
double up_height = 0.8; // double dw_correction = 1.30;
double dw_correction = 1.18; double font_size_dw  = 0.1;
double dw_height = (1. - up_height) * dw_correction;
double dw_height_offset = 0.02;

//for running
ofstream myfile;
bool runABCDPlots = true;
bool runFullBkg = true;
bool runRPFPlots = true;
bool runStacks = false;
bool runMETNorm = true;
bool runTableOfYields = true;
bool runPies = true;


// string whichRegion = "signal";
TString whichRegion = "singleLept";
// string whichRegion = "photon";

TString year = "Run2";


TFile * f = TFile::Open("test.root");
TFile * fSignal = TFile::Open("ALPHABETBoost_V17_signalOnly.root");
TFile * fPhoton;// = TFile::Open("ALPHABETBoostMC2016_V12like_photon.root");
TFile * fSingleLept = TFile::Open("ALPHABET_V17_1l_MET300_2BoostedH.root");

TString fout_name = "ABCDPlots_"+whichRegion+"_test.root";
TFile* fout = new TFile(fout_name,"recreate");


TDirectory *cdABCD  = fout->mkdir("ABCD");
TDirectory *cdClose = fout->mkdir("Closure");
TDirectory *cdAssum = fout->mkdir("Assumption");
TDirectory *cdRPFFinal  = fout->mkdir("RPF");
TDirectory *cdRPF  = fout->mkdir("RPF_Support");
TDirectory *cdPie  = fout->mkdir("PieCharts");
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
  TH1F * h_0HSBOpt3; TH1F * h_0HSROpt3;
  TH1F * h_0HSBOpt4; TH1F * h_0HSROpt4;

  TH1F * h_0HSBOpt1_TT; TH1F * h_0HSBOpt1_WJets; TH1F * h_0HSBOpt1_ZJets; TH1F * h_0HSBOpt1_QCD;
  TH1F * h_0HSROpt1_TT; TH1F * h_0HSROpt1_WJets; TH1F * h_0HSROpt1_ZJets; TH1F * h_0HSROpt1_QCD;
  TH1F * h_0HSBOpt2_TT; TH1F * h_0HSBOpt2_WJets; TH1F * h_0HSBOpt2_ZJets; TH1F * h_0HSBOpt2_QCD;
  TH1F * h_0HSROpt2_TT; TH1F * h_0HSROpt2_WJets; TH1F * h_0HSROpt2_ZJets; TH1F * h_0HSROpt2_QCD;
  TH1F * h_0HSBOpt3_TT; TH1F * h_0HSBOpt3_WJets; TH1F * h_0HSBOpt3_ZJets; TH1F * h_0HSBOpt3_QCD;
  TH1F * h_0HSROpt3_TT; TH1F * h_0HSROpt3_WJets; TH1F * h_0HSROpt3_ZJets; TH1F * h_0HSROpt3_QCD;
  TH1F * h_0HSBOpt4_TT; TH1F * h_0HSBOpt4_WJets; TH1F * h_0HSBOpt4_ZJets; TH1F * h_0HSBOpt4_QCD;
  TH1F * h_0HSROpt4_TT; TH1F * h_0HSROpt4_WJets; TH1F * h_0HSROpt4_ZJets; TH1F * h_0HSROpt4_QCD;

  if (whichRegion=="signal") {
    h_2HLeadSR_sum = (TH1F*)f->Get("J2pt_M_doubletagSRLead_sum"); h_0HLeadSR_sum = (TH1F*)f->Get("J2pt_M_antitagSRLead_sum");
    h_2HLeadSB_sum = (TH1F*)f->Get("J1pt_M_doubletagSBLead_sum"); h_0HLeadSB_sum = (TH1F*)f->Get("J1pt_M_antitagSBLead_sum");

    // h_A_sum = (TH1F*)f->Get("MET_doubletagSR_sum"); h_B_sum = (TH1F*)f->Get("MET_doubletagSB_sum");
    // h_A1_sum = (TH1F*)f->Get("MET_tagSR_sum"); h_B1_sum = (TH1F*)f->Get("MET_tagSB_sum");
    // h_C_sum = (TH1F*)f->Get("MET_antitagSR_sum"); h_D_sum = (TH1F*)f->Get("MET_antitagSB_sum");

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

    h_0HSROpt1_QCD = (TH1F*)f->Get("MET_antitagSROpt1_QCD"); h_0HSBOpt1_QCD = (TH1F*)f->Get("MET_antitagSBOpt1_QCD");
    h_0HSROpt1_TT = (TH1F*)f->Get("MET_antitagSROpt1_TT"); h_0HSBOpt1_TT = (TH1F*)f->Get("MET_antitagSBOpt1_TT");
    h_0HSROpt1_WJets = (TH1F*)f->Get("MET_antitagSROpt1_WJets"); h_0HSBOpt1_WJets = (TH1F*)f->Get("MET_antitagSBOpt1_WJets");
    h_0HSROpt1_ZJets = (TH1F*)f->Get("MET_antitagSROpt1_ZJets"); h_0HSBOpt1_ZJets = (TH1F*)f->Get("MET_antitagSBOpt1_ZJets");

    h_0HSROpt2_QCD = (TH1F*)f->Get("MET_antitagSROpt2_QCD"); h_0HSBOpt2_QCD = (TH1F*)f->Get("MET_antitagSBOpt2_QCD");
    h_0HSROpt2_TT = (TH1F*)f->Get("MET_antitagSROpt2_TT"); h_0HSBOpt2_TT = (TH1F*)f->Get("MET_antitagSBOpt2_TT");
    h_0HSROpt2_WJets = (TH1F*)f->Get("MET_antitagSROpt2_WJets"); h_0HSBOpt2_WJets = (TH1F*)f->Get("MET_antitagSBOpt2_WJets");
    h_0HSROpt2_ZJets = (TH1F*)f->Get("MET_antitagSROpt2_ZJets"); h_0HSBOpt2_ZJets = (TH1F*)f->Get("MET_antitagSBOpt2_ZJets");

    h_0HSROpt3_QCD = (TH1F*)f->Get("MET_antitagSROpt3_QCD"); h_0HSBOpt3_QCD = (TH1F*)f->Get("MET_antitagSBOpt3_QCD");
    h_0HSROpt3_TT = (TH1F*)f->Get("MET_antitagSROpt3_TT"); h_0HSBOpt3_TT = (TH1F*)f->Get("MET_antitagSBOpt3_TT");
    h_0HSROpt3_WJets = (TH1F*)f->Get("MET_antitagSROpt3_WJets"); h_0HSBOpt3_WJets = (TH1F*)f->Get("MET_antitagSBOpt3_WJets");
    h_0HSROpt3_ZJets = (TH1F*)f->Get("MET_antitagSROpt3_ZJets"); h_0HSBOpt3_ZJets = (TH1F*)f->Get("MET_antitagSBOpt3_ZJets");

    h_0HSROpt4_QCD = (TH1F*)f->Get("MET_antitagSROpt4_QCD"); h_0HSBOpt4_QCD = (TH1F*)f->Get("MET_antitagSBOpt4_QCD");
    h_0HSROpt4_TT = (TH1F*)f->Get("MET_antitagSROpt4_TT"); h_0HSBOpt4_TT = (TH1F*)f->Get("MET_antitagSBOpt4_TT");
    h_0HSROpt4_WJets = (TH1F*)f->Get("MET_antitagSROpt4_WJets"); h_0HSBOpt4_WJets = (TH1F*)f->Get("MET_antitagSBOpt4_WJets");
    h_0HSROpt4_ZJets = (TH1F*)f->Get("MET_antitagSROpt4_ZJets"); h_0HSBOpt4_ZJets = (TH1F*)f->Get("MET_antitagSBOpt4_ZJets");

    //Adding uncertainties to 0-event bins
    //TH1F make0EventUncSum(vector<TH1F*> dem_histos); //in this order: QCD, TT, WJets,ZJets
    h_A_sum = make0EventUncSum({h_A_QCD,h_A_TT,h_A_WJets,h_A_ZJets});
    h_B_sum = make0EventUncSum({h_B_QCD,h_B_TT,h_B_WJets,h_B_ZJets});
    h_A1_sum = make0EventUncSum({h_A1_QCD,h_A1_TT,h_A1_WJets,h_A1_ZJets});
    h_B1_sum = make0EventUncSum({h_B1_QCD,h_B1_TT,h_B1_WJets,h_B1_ZJets});
    h_C_sum = make0EventUncSum({h_C_QCD,h_C_TT,h_C_WJets,h_C_ZJets});
    h_D_sum = make0EventUncSum({h_D_QCD,h_D_TT,h_D_WJets,h_D_ZJets});

    h_0HSROpt1 = make0EventUncSum({h_0HSROpt1_QCD,h_0HSROpt1_TT,h_0HSROpt1_WJets,h_0HSROpt1_ZJets});
    h_0HSBOpt1 = make0EventUncSum({h_0HSBOpt1_QCD,h_0HSBOpt1_TT,h_0HSBOpt1_WJets,h_0HSBOpt1_ZJets});
    h_0HSROpt2 = make0EventUncSum({h_0HSROpt2_QCD,h_0HSROpt2_TT,h_0HSROpt2_WJets,h_0HSROpt2_ZJets});
    h_0HSBOpt2 = make0EventUncSum({h_0HSBOpt2_QCD,h_0HSBOpt2_TT,h_0HSBOpt2_WJets,h_0HSBOpt2_ZJets});
    h_0HSROpt3 = make0EventUncSum({h_0HSROpt3_QCD,h_0HSROpt3_TT,h_0HSROpt3_WJets,h_0HSROpt3_ZJets});
    h_0HSBOpt3 = make0EventUncSum({h_0HSBOpt3_QCD,h_0HSBOpt3_TT,h_0HSBOpt3_WJets,h_0HSBOpt3_ZJets});
    h_0HSROpt4 = make0EventUncSum({h_0HSROpt4_QCD,h_0HSROpt4_TT,h_0HSROpt4_WJets,h_0HSROpt4_ZJets});
    h_0HSBOpt4 = make0EventUncSum({h_0HSBOpt4_QCD,h_0HSBOpt4_TT,h_0HSBOpt4_WJets,h_0HSBOpt4_ZJets});

    h_A_sum->Write("h_A_sum");
    h_A1_sum->Write("h_A1_sum");
    h_B_sum->Write("h_B_sum");
    h_B1_sum->Write("h_B1_sum");
    h_C_sum->Write("h_C_sum");
    h_D_sum->Write("h_D_sum");
    h_0HSROpt1->Write("h_COpt1_sum");
    h_0HSROpt2->Write("h_COpt2_sum");
    h_0HSROpt3->Write("h_COpt3_sum");
    h_0HSROpt4->Write("h_COpt4_sum");
    h_0HSBOpt1->Write("h_DOpt1_sum");
    h_0HSBOpt2->Write("h_DOpt2_sum");
    h_0HSBOpt3->Write("h_DOpt3_sum");
    h_0HSBOpt4->Write("h_DOpt4_sum");


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

    // h_0HSROpt1 = (TH1F*)f->Get("MET_antitagSROpt1_sum"); h_0HSBOpt1 = (TH1F*)f->Get("MET_antitagSBOpt1_sum");
    // h_0HSROpt2 = (TH1F*)f->Get("MET_antitagSROpt2_sum"); h_0HSBOpt2 = (TH1F*)f->Get("MET_antitagSBOpt2_sum");
    // h_0HSROpt3 = (TH1F*)f->Get("MET_antitagSROpt3_sum"); h_0HSBOpt3 = (TH1F*)f->Get("MET_antitagSBOpt3_sum");
    // h_0HSROpt4 = (TH1F*)f->Get("MET_antitagSROpt4_sum"); h_0HSBOpt4 = (TH1F*)f->Get("MET_antitagSBOpt4_sum");
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

    // h_0HSBOpt1 = (TH1F*)fPhoton->Get("MET_antitagSBOpt1_sum"); h_0HSBOpt2 = (TH1F*)fPhoton->Get("MET_antitagSBOpt2_sum");
    // h_0HSROpt1 = (TH1F*)fPhoton->Get("MET_antitagSROpt1_sum"); h_0HSROpt2 = (TH1F*)fPhoton->Get("MET_antitagSROpt2_sum");
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

    h_0HSROpt1_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_TT"); h_0HSBOpt1_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_TT");
    h_0HSROpt1_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_WJets"); h_0HSBOpt1_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_WJets");
    h_0HSROpt2_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt2_TT"); h_0HSBOpt2_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt2_TT");
    h_0HSROpt2_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt2_WJets"); h_0HSBOpt2_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt2_WJets");
    h_0HSROpt3_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt3_TT"); h_0HSBOpt3_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt3_TT");
    h_0HSROpt3_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt3_WJets"); h_0HSBOpt3_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt3_WJets");
    h_0HSROpt4_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt4_TT"); h_0HSBOpt4_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt4_TT");
    h_0HSROpt4_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt4_WJets"); h_0HSBOpt4_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt4_WJets");

    // h_A_sum = (TH1F*)h_A_TT->Clone("h_A_sum");  h_B_sum = (TH1F*)h_B_TT->Clone("h_B_sum");
    // h_A1_sum = (TH1F*)h_A1_TT->Clone("h_A1_sum");  h_B1_sum = (TH1F*)h_B1_TT->Clone("h_B1_sum");
    // h_C_sum = (TH1F*)h_C_TT->Clone("h_C_sum");  h_D_sum = (TH1F*)h_D_TT->Clone("h_D_sum");
    // h_A_sum->Add(h_A_WJets); h_B_sum->Add(h_B_WJets);
    // h_A1_sum->Add(h_A1_WJets); h_B1_sum->Add(h_B1_WJets);
    // h_C_sum->Add(h_C_WJets); h_D_sum->Add(h_D_WJets);
    //Adding uncertainties to 0-event bins

    //TH1F make0EventUncSum_1l(vector<TH1F*> dem_histos); //in this order: TT, WJets
    h_A_sum = make0EventUncSum_1l({h_A_TT,h_A_WJets});
    h_B_sum = make0EventUncSum_1l({h_B_TT,h_B_WJets});
    h_A1_sum = make0EventUncSum_1l({h_A1_TT,h_A1_WJets});
    h_B1_sum = make0EventUncSum_1l({h_B1_TT,h_B1_WJets});
    h_C_sum = make0EventUncSum_1l({h_C_TT,h_C_WJets});
    h_D_sum = make0EventUncSum_1l({h_D_TT,h_D_WJets,});

    h_A_sum->Write("h_A_sum");
    h_A1_sum->Write("h_A1_sum");
    h_B_sum->Write("h_B_sum");
    h_B1_sum->Write("h_B1_sum");
    h_C_sum->Write("h_C_sum");
    h_D_sum->Write("h_D_sum");

    h_0HSROpt1 = make0EventUncSum_1l({h_0HSROpt1_TT,h_0HSROpt1_WJets});
    h_0HSBOpt1 = make0EventUncSum_1l({h_0HSBOpt1_TT,h_0HSBOpt1_WJets});
    h_0HSROpt2 = make0EventUncSum_1l({h_0HSROpt2_TT,h_0HSROpt2_WJets});
    h_0HSBOpt2 = make0EventUncSum_1l({h_0HSBOpt2_TT,h_0HSBOpt2_WJets});
    h_0HSROpt3 = make0EventUncSum_1l({h_0HSROpt3_TT,h_0HSROpt3_WJets});
    h_0HSBOpt3 = make0EventUncSum_1l({h_0HSBOpt3_TT,h_0HSBOpt3_WJets});
    h_0HSROpt4 = make0EventUncSum_1l({h_0HSROpt4_TT,h_0HSROpt4_WJets});
    h_0HSBOpt4 = make0EventUncSum_1l({h_0HSBOpt4_TT,h_0HSBOpt4_WJets});

    h_0HSROpt1->Write("h_COpt1_sum");
    h_0HSROpt2->Write("h_COpt2_sum");
    h_0HSROpt3->Write("h_COpt3_sum");
    h_0HSROpt4->Write("h_COpt4_sum");
    h_0HSBOpt1->Write("h_DOpt1_sum");
    h_0HSBOpt2->Write("h_DOpt2_sum");
    h_0HSBOpt3->Write("h_DOpt3_sum");
    h_0HSBOpt4->Write("h_DOpt4_sum");


    h_J1M_doubletagSR_sum = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSR_sum"); h_J2M_doubletagSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSR_sum");
    h_J1M_doubletagSB_sum = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSB_sum"); h_J2M_doubletagSB_sum = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSB_sum");
    h_J1M_antitagSR_sum = (TH1F*)fSingleLept->Get("J1pt_M_antitagSR_sum"); h_J2M_antitagSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_antitagSR_sum");
    h_J1M_antitagSB_sum = (TH1F*)fSingleLept->Get("J1pt_M_antitagSB_sum"); h_J2M_antitagSB_sum = (TH1F*)fSingleLept->Get("J2pt_M_antitagSB_sum");
    h_J1M_tagSR_sum = (TH1F*)fSingleLept->Get("J1pt_M_tagSR_sum"); h_J2M_tagSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_tagSR_sum");
    h_J1M_tagSB_sum = (TH1F*)fSingleLept->Get("J1pt_M_tagSB_sum"); h_J2M_tagSB_sum = (TH1F*)fSingleLept->Get("J2pt_M_tagSB_sum");

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

    h_J2M_mjBins_doubletagSR_sum = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSR_sum"); h_J2M_mjBins_doubletagSB_sum = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSB_sum");
    h_J2M_mjBins_antitagSR_sum = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSR_sum"); h_J2M_mjBins_antitagSB_sum = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSB_sum");
    h_J2M_mjBins_doubletagSR_TT = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSR_TT"); h_J2M_mjBins_doubletagSB_TT = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSB_TT");
    h_J2M_mjBins_antitagSR_TT = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSR_TT"); h_J2M_mjBins_antitagSB_TT = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSB_TT");

    // h_J2M_mjBins_doubletagSR_TT_Di = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_TT_Di"); h_J2M_mjBins_doubletagSB_TT_Di = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_TT_Di");
    // h_J2M_mjBins_antitagSR_TT_Di = (TH1F*)f->Get("J2_M_jetBins_antitagSR_TT_Di"); h_J2M_mjBins_antitagSB_TT_Di = (TH1F*)f->Get("J2_M_jetBins_antitagSB_TT_Di");
    // h_J2M_mjBins_doubletagSR_TT_SL = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_TT_SL"); h_J2M_mjBins_doubletagSB_TT_SL = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_TT_SL");
    // h_J2M_mjBins_antitagSR_TT_SL = (TH1F*)f->Get("J2_M_jetBins_antitagSR_TT_SL"); h_J2M_mjBins_antitagSB_TT_SL = (TH1F*)f->Get("J2_M_jetBins_antitagSB_TT_SL");

    h_J2M_mjBins_doubletagSR_WJets = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSR_WJets"); h_J2M_mjBins_doubletagSB_WJets = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSB_WJets");
    h_J2M_mjBins_antitagSR_WJets = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSR_WJets"); h_J2M_mjBins_antitagSB_WJets = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSB_WJets");
    //
    // h_0HSBOpt1 = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_sum"); h_0HSROpt1 = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_sum");
    // h_0HSBOpt2 = (TH1F*)fSingleLept->Get("MET_antitagSBOpt2_sum"); h_0HSROpt2 = (TH1F*)fSingleLept->Get("MET_antitagSROpt2_sum");
    // h_0HSBOpt3 = (TH1F*)fSingleLept->Get("MET_antitagSBOpt3_sum"); h_0HSROpt3 = (TH1F*)fSingleLept->Get("MET_antitagSROpt3_sum");
    // h_0HSBOpt4 = (TH1F*)fSingleLept->Get("MET_antitagSBOpt4_sum"); h_0HSROpt4 = (TH1F*)fSingleLept->Get("MET_antitagSROpt4_sum");

    // h_0HSBOpt1_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_TT"); h_0HSBOpt1_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_WJets");
    // h_0HSROpt1_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_TT"); h_0HSROpt1_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_WJets");
    //
    // //set up sum for opt cuts
    // h_0HSBOpt1 = (TH1F*)h_0HSBOpt1_TT->Clone("h_0HSBOpt1_sum");
    // h_0HSROpt1 = (TH1F*)h_0HSROpt1_TT->Clone("h_0HSROpt1_sum");
    // h_0HSBOpt1->Add(h_0HSBOpt1_WJets);
    // h_0HSROpt1->Add(h_0HSROpt1_WJets);


    // h_2HLeadSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSRLead_sum"); h_0HLeadSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_antitagSRLead_sum");

  }

  vector<TH1F*> histos_ABCD_data = {h_A_data, h_B_data, h_C_data, h_D_data};
  vector<TH1F*> histos_A1B1CD_data = {h_A1_data, h_B1_data, h_C_data, h_D_data};

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

  // vector<TH1F*> h_METShape_1H = {h_A1_sum,h_D_sum,h_0HSBOpt4,h_C_sum,h_0HSROpt4};
  vector<TH1F*> h_METShape_1H = {h_A1_sum,h_D_sum,h_0HSBOpt2,h_C_sum,h_0HSROpt2};

  if (runPies){
    if (whichRegion=="signal"){
      pieChart({h_A_QCD,h_B_QCD}, {h_A_WJets,h_B_WJets}, {h_A_ZJets,h_B_ZJets}, {h_A_TT,h_B_TT}, "2H");
      pieChart({h_A_QCD}, {h_A_WJets}, {h_A_ZJets}, {h_A_TT}, "2HSR");
      pieChart({h_B_QCD}, {h_B_WJets}, {h_B_ZJets}, {h_B_TT}, "2HSB");

      pieChart({h_A1_QCD,h_B1_QCD}, {h_A1_WJets,h_B1_WJets}, {h_A1_ZJets,h_B1_ZJets}, {h_A1_TT,h_B1_TT}, "1H");
      pieChart({h_A1_QCD}, {h_A1_WJets}, {h_A1_ZJets}, {h_A1_TT}, "1HSR");
      pieChart({h_B1_QCD}, {h_B1_WJets}, {h_B1_ZJets}, {h_B1_TT}, "1HSB");

      pieChart({h_C_QCD,h_D_QCD}, {h_C_WJets,h_D_WJets}, {h_C_ZJets,h_D_ZJets}, {h_C_TT,h_D_TT}, "0H");
      pieChart({h_C_QCD}, {h_C_WJets}, {h_C_ZJets}, {h_C_TT}, "0HSR");
      pieChart({h_D_QCD}, {h_D_WJets}, {h_D_ZJets}, {h_D_TT}, "0HSB");

      pieChart({h_0HSROpt2_QCD,h_0HSBOpt2_QCD}, {h_0HSROpt2_WJets,h_0HSBOpt2_WJets}, {h_0HSROpt2_ZJets,h_0HSBOpt2_ZJets}, {h_0HSROpt2_TT,h_0HSBOpt2_TT}, "0H,BTagsM>0");
      pieChart({h_0HSROpt1_QCD,h_0HSBOpt1_QCD}, {h_0HSROpt1_WJets,h_0HSBOpt1_WJets}, {h_0HSROpt1_ZJets,h_0HSBOpt1_ZJets}, {h_0HSROpt1_TT,h_0HSBOpt1_TT}, "0H,BTagsM>1");
    }
    else if (whichRegion=="singleLept"){
      pieChart1l({h_A_WJets,h_B_WJets}, {h_A_TT,h_B_TT}, "2H");
      pieChart1l({h_A_WJets}, {h_A_TT}, "2HSR");
      pieChart1l({h_B_WJets},{h_B_TT}, "2HSB");

      pieChart1l({h_A1_WJets,h_B1_WJets}, {h_A1_TT,h_B1_TT}, "1H");
      pieChart1l({h_A1_WJets}, {h_A1_TT}, "1HSR");
      pieChart1l({h_B1_WJets}, {h_B1_TT}, "1HSB");

      pieChart1l({h_C_WJets,h_D_WJets}, {h_C_TT,h_D_TT}, "0H");
      pieChart1l({h_C_WJets}, {h_C_TT}, "0HSR");
      pieChart1l({h_D_WJets}, {h_D_TT}, "0HSB");

      pieChart1l({h_0HSROpt2_WJets,h_0HSBOpt2_WJets}, {h_0HSROpt2_TT,h_0HSBOpt2_TT}, "0H,BTagsM>0");
      pieChart1l({h_0HSROpt1_WJets,h_0HSBOpt1_WJets}, {h_0HSROpt1_TT,h_0HSBOpt1_TT}, "0H,BTagsM>1");
    }
  }


  if (runFullBkg) {
    makeFullBkgClosure(histos_ABCD_sum, "BkgSum", "Double",year);
    makeFullBkgClosure(histos_A1B1CD_sum, "BkgSum", "Single",year);
  }

  if (runABCDPlots){
    if (whichRegion=="signal"){
      makeABCDPlot(histos_ABCD_sum, "BkgSum", "Double",year);
      makeABCDPlot(histos_A1B1CD_sum, "BkgSum", "Single",year);

      makeABCDPlot(histos_ABCD_QCD, "QCD", "Double",year);
      makeABCDPlot(histos_ABCD_TT, "TT", "Double",year);
      makeABCDPlot(histos_ABCD_WJets, "WJets", "Double",year);
      makeABCDPlot(histos_ABCD_ZJets, "ZJets", "Double",year);

      makeABCDPlot(histos_A1B1CD_QCD, "QCD", "Single",year);
      makeABCDPlot(histos_A1B1CD_TT, "TT", "Single",year);
      makeABCDPlot(histos_A1B1CD_WJets, "WJets", "Single",year);
      makeABCDPlot(histos_A1B1CD_ZJets, "ZJets", "Single",year);
    }

    else if (whichRegion=="singleLept"){
      makeABCDPlot(histos_ABCD_sum, "BkgSum", "Double",year);
      makeABCDPlot(histos_A1B1CD_sum, "BkgSum", "Single",year);

      // makeABCDPlot(histos_ABCD_data, "Data", "Double",year);
      // makeABCDPlot(histos_A1B1CD_data, "Data", "Single",year);

      makeABCDPlot(histos_ABCD_TT, "TT", "Double",year);
      makeABCDPlot(histos_ABCD_WJets, "WJets", "Double",year);

      makeABCDPlot(histos_A1B1CD_TT, "TT", "Single",year);
      makeABCDPlot(histos_A1B1CD_WJets, "WJets", "Single",year);
    }
    else if (whichRegion=="photon"){
      makeABCDPlot(histos_ABCD_QCD, "QCD", "Double",year);
      makeABCDPlot(histos_ABCD_GJets, "GJets", "Double",year);
      makeABCDPlot(histos_A1B1CD_QCD, "QCD", "Single",year);
      makeABCDPlot(histos_A1B1CD_GJets, "GJets", "Single",year);
    }
  }

  if (runRPFPlots){
    if (whichRegion=="signal"){
      //J1 Double tag region
      makeRpfPlot(histos_Rpf_J1_sum, "BkgSum", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_QCD, "QCD", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_TT, "TT", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_WJets, "WJets", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_ZJets, "ZJets", "J1", "Double",year);

      //J1 Single tag region
      makeRpfPlot(histos_Rpfsingle_J1_sum, "BkgSum", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_QCD, "QCD", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_TT, "TT", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_WJets, "WJets", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_ZJets, "ZJets", "J1", "Single",year);

      //J2 double tag region
      makeRpfPlot(histos_Rpf_J2_sum, "BkgSum", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_QCD, "QCD", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_TT, "TT", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_WJets, "WJets", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_ZJets, "ZJets", "J2", "Double",year);

      //J2 double tag region, with three mass bins
      makeRpfPlot(histos_Rpf_J2_mJBins_sum, "BkgSum_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_QCD, "QCD_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_TT, "TT_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_WJets, "WJets_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_ZJets, "ZJets_mJ", "J2", "Double",year);


      //J2 single tag region
      makeRpfPlot(histos_Rpfsingle_J2_sum, "BkgSum", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_QCD, "QCD", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_TT, "TT", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_WJets, "WJets", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_ZJets, "ZJets", "J2", "Single",year);
    }
    if (whichRegion=="singleLept"){
      //J1 Double tag region
      makeRpfPlot(histos_Rpf_J1_sum, "BkgSum", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_TT, "TT", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_WJets, "WJets", "J1", "Double",year);


      //J1 Single tag region
      makeRpfPlot(histos_Rpfsingle_J1_sum, "BkgSum", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_TT, "TT", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_WJets, "WJets", "J1", "Single",year);

      //J2 double tag region
      makeRpfPlot(histos_Rpf_J2_sum, "BkgSum", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_TT, "TT", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_WJets, "WJets", "J2", "Double",year);

      //J2 single tag region
      makeRpfPlot(histos_Rpfsingle_J2_sum, "BkgSum", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_TT, "TT", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_WJets, "WJets", "J2", "Single",year);

      //J2 double tag region, with three mass bins
      makeRpfPlot(histos_Rpf_J2_mJBins_sum, "BkgSum_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_TT, "TT_mJ", "J2", "Double",year);
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
    }
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
    makeMETNorm(h_METShape_2H, "Double");
    makeMETNorm(h_METShape_1H, "Single");

    tableOfMETNorm(h_METShape_2H, "BkgSum", "Double");
    tableOfMETNorm(h_METShape_1H, "BkgSum", "Single");
  }

  if (runTableOfYields){
    if (whichRegion=="signal"){
      tableOfYields(histos_ABCD_sum, "BkgSum", "Double");
      tableOfYields(histos_A1B1CD_sum, "BkgSum", "Single");

      tableOfYields(histos_ABCD_QCD, "QCD", "Double");
      tableOfYields(histos_ABCD_TT, "TT", "Double");
      tableOfYields(histos_ABCD_WJets, "WJets", "Double");
      tableOfYields(histos_ABCD_ZJets, "ZJets", "Double");

      tableOfYields(histos_A1B1CD_QCD, "QCD", "Single");
      tableOfYields(histos_A1B1CD_TT, "TT", "Single");
      tableOfYields(histos_A1B1CD_WJets, "WJets", "Single");
      tableOfYields(histos_A1B1CD_ZJets, "ZJets", "Single");
    }

    else if (whichRegion=="singleLept"){
      tableOfYields(histos_ABCD_sum, "BkgSum", "Double");
      tableOfYields(histos_A1B1CD_sum, "BkgSum", "Single");

      tableOfYields(histos_ABCD_TT, "TT", "Double");
      tableOfYields(histos_ABCD_WJets, "WJets", "Double");

      tableOfYields(histos_A1B1CD_TT, "TT", "Single");
      tableOfYields(histos_A1B1CD_WJets, "WJets", "Single");
    }
    else if (whichRegion=="photon"){
      tableOfYields(histos_ABCD_QCD, "QCD", "Double");
      tableOfYields(histos_ABCD_GJets, "GJets", "Double");
      tableOfYields(histos_A1B1CD_QCD, "QCD", "Single");
      tableOfYields(histos_A1B1CD_GJets, "GJets", "Single");
    }
  }

  // fout->Close();
  // f->Close();
  // fSignal->Close();
} //end method ABCD()


void makeABCDPlot(vector<TH1F*> dem_histos, TString bkgType = "", TString tagType = "", TString year = "") {
  TH1F *histo_A = (TH1F*)dem_histos[0]->Clone("histo_A"); TH1F *histo_B = (TH1F*)dem_histos[1]->Clone("histo_B");
  TH1F *histo_C = (TH1F*)dem_histos[2]->Clone("histo_C"); TH1F *histo_D = (TH1F*)dem_histos[3]->Clone("histo_D");

  bool isData = false;
  TString histoName = dem_histos[0]->GetName();
  if (histoName.Contains("data")) isData=true;


  int numMETBins = histo_A->GetNbinsX();
  for (int i = 1; i<=numMETBins; i++){
    if (histo_A->GetBinContent(i)<0.001){
      if (histoName.Contains("TT")) histo_A->SetBinError(i,0.795);
      if (histoName.Contains("WJets")) histo_A->SetBinError(i,0.210);
      if (histoName.Contains("ZJets")) histo_A->SetBinError(i,0.328);
    }
    if (histo_B->GetBinContent(i)<0.001){
      if (histoName.Contains("TT")) histo_B->SetBinError(i,0.795);
      if (histoName.Contains("WJets")) histo_B->SetBinError(i,0.210);
      if (histoName.Contains("ZJets")) histo_B->SetBinError(i,0.328);
    }
    if (histo_C->GetBinContent(i)<0.001){
      if (histoName.Contains("TT")) histo_C->SetBinError(i,0.795);
      if (histoName.Contains("WJets")) histo_C->SetBinError(i,0.210);
      if (histoName.Contains("ZJets")) histo_C->SetBinError(i,0.328);
    }
    if (histo_D->GetBinContent(i)<0.001){
      if (histoName.Contains("TT")) histo_D->SetBinError(i,0.795);
      if (histoName.Contains("WJets")) histo_D->SetBinError(i,0.210);
      if (histoName.Contains("ZJets")) histo_D->SetBinError(i,0.328);
    }
  }

  TH1F *histo_pred = (TH1F*)histo_B->Clone("histo_pred");
  histo_pred->Multiply(histo_C);histo_pred->Divide(histo_D);

  TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo_A, histo_pred, "pois");
  graph->SetTitle(";MET [GeV]; #kappa (Sim/Pred)");
  styleGraph(graph);

  TString canvName = bkgType+tagType;
  TCanvas * can_h = new TCanvas(canvName,canvName, 50, 50, W, H);
  styleCanvas(can_h);


  TPad *pad1 = new TPad("pad1", "top pad" , 0.0, 0.25, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.25);
  pad1->SetTickx(0); pad1->SetTicky(0);
  pad1->SetPad(0., 1 - up_height, 1., 1.00);
  pad1->SetFrameFillColor(0); pad1->SetFillColor(0);
  pad1->SetTopMargin(0.08); pad1->SetLeftMargin(0.10); pad1->SetRightMargin(0.06);
  pad1->Draw();

  pad2->SetPad(0., 0., 1., dw_height+dw_height_offset);
  pad2->SetFillColor(0); pad2->SetFrameFillColor(0);
  pad2->SetBottomMargin(0.33); pad2->SetTopMargin(0.01);
  pad2->SetLeftMargin(0.10); pad2->SetRightMargin(0.06);
  pad2->Draw();
  pad1->cd();

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
    histo_pred->GetYaxis()->SetTitleOffset(0.8);
    histo_pred->GetXaxis()->SetLabelSize(0);
    histo_pred->Draw("E2");
    histo_A->Draw("same");
  }
  else {
    histo_A->GetYaxis()->SetTitleOffset(0.8);
    histo_A->GetXaxis()->SetLabelSize(0);
    histo_A->Draw("same");
    histo_pred->Draw("E2 same");
  }
  legend->Draw("same");
  pad2->cd();

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
  delete can_h;
}

void makeFullBkgClosure(vector<TH1F*> dem_histos, TString bkgType = "", TString tagType = "", TString year = "") {
  TH1F *histo_A = (TH1F*)dem_histos[0]->Clone("histo_A"); TH1F *histo_B = (TH1F*)dem_histos[1]->Clone("histo_B");
  TH1F *histo_C = (TH1F*)dem_histos[2]->Clone("histo_C"); TH1F *histo_D = (TH1F*)dem_histos[3]->Clone("histo_D");

  double b_int; double b_error; b_int = histo_B->IntegralAndError(0,5,b_error,"");
  double c_int; double c_error; c_int = histo_C->IntegralAndError(0,5,c_error,"");
  double d_int; double d_error; d_int = histo_D->IntegralAndError(0,5,d_error,"");

  TH1F *histo_pred = (TH1F*)histo_B->Clone("histo_pred");
  histo_pred->Multiply(histo_C);histo_pred->Divide(histo_D);

  double bkgNorm; double bkgNorm_error;
  bkgNorm = histo_pred->IntegralAndError(0,5,bkgNorm_error,"");

  TH1F *h_finalBkg = (TH1F*)histo_pred->Clone("h_finalBkg");

  float kappa; float kappa_error;
  float bkgFrac1; float bkgFrac1_error;
  float bkgFrac2; float bkgFrac2_error;
  float bkgFrac3; float bkgFrac3_error;

  if (whichRegion=="signal"){
    if (tagType=="Double"){
      kappa = 1.1; kappa_error = 0.18;
      bkgFrac1 = 0.857;   bkgFrac1_error = 0.023;
      bkgFrac2 = 0.118;  bkgFrac2_error = 0.0079;
      bkgFrac3 = 0.0251; bkgFrac3_error = 0.0033;
    }
    else if (tagType=="Single"){
      kappa = 0.95; kappa_error = 0.05;
      bkgFrac1 = 0.821;   bkgFrac1_error = 0.0111;
      bkgFrac2 = 0.136;  bkgFrac2_error = 0.0040;
      bkgFrac3 = 0.044; bkgFrac3_error = 0.0020;
    }
  }

  if (whichRegion=="singleLept"){
    if (tagType=="Double"){
      kappa = 0.66; kappa_error = 0.137;
      bkgFrac1 = 0.923;  bkgFrac1_error = 0.0149;
      bkgFrac2 = 0.067; bkgFrac2_error = 0.0038;
      bkgFrac3 = 0.009; bkgFrac3_error = 0.0013;
    }
    else if (tagType=="Single"){
      kappa = 0.91; kappa_error = 0.0424;
      bkgFrac1 = 0.911;  bkgFrac1_error = 0.0089;
      bkgFrac2 = 0.077; bkgFrac2_error = 0.0024;
      bkgFrac3 = 0.012; bkgFrac3_error = 0.0008;
    }
  }

  float bin1Content = bkgNorm*kappa*bkgFrac1;
  float bin2Content = bkgNorm*kappa*bkgFrac2;
  float bin3Content = bkgNorm*kappa*bkgFrac3;

  h_finalBkg->SetBinContent(1,bin1Content);
  h_finalBkg->SetBinContent(2,bin2Content);
  h_finalBkg->SetBinContent(3,bin3Content);

  float bin1Error = bin1Content*sqrt( TMath::Power(bkgNorm_error/bkgNorm,2) + TMath::Power(kappa_error/kappa,2) + TMath::Power(bkgFrac1_error/bkgFrac1,2) );
  float bin2Error = bin2Content*sqrt( TMath::Power(bkgNorm_error/bkgNorm,2) + TMath::Power(kappa_error/kappa,2) + TMath::Power(bkgFrac2_error/bkgFrac2,2) );
  float bin3Error = bin3Content*sqrt( TMath::Power(bkgNorm_error/bkgNorm,2) + TMath::Power(kappa_error/kappa,2) + TMath::Power(bkgFrac3_error/bkgFrac3,2) );

  h_finalBkg->SetBinError(1,bin1Error);
  h_finalBkg->SetBinError(2,bin2Error);
  h_finalBkg->SetBinError(3,bin3Error);

  TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo_A, h_finalBkg, "pois");
  styleGraph(graph);
  graph->SetTitle(";MET [GeV]; MC/Pred");

  TString canvName = bkgType+tagType;
  TCanvas * can_h = new TCanvas(canvName,canvName, 50, 50, W, H);
  styleCanvas(can_h);


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

  TLegend* legend = new TLegend(0.40,0.7,0.93,0.91);
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
  can_h->Write("FullBkg_"+tagType);
  delete can_h;
}

void makeRpfPlot(vector<TH1F*> dem_histos, TString bkgType = "", TString jetType = "", TString tagType = "", TString year = "") {
  TH1D * num   = (TH1D*)dem_histos[1]->Clone("num"+tagType+"_"+jetType+"_"+bkgType); //this should be 2H SB
  TH1D * denom = (TH1D*)dem_histos[3]->Clone("denom"+tagType+"_"+jetType+"_"+bkgType);  //this should be 0H SB
  num->Add(dem_histos[0]); denom->Add(dem_histos[2]); //Add 2H SR to num, and 0H SR to denom

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
  TGraphAsymmErrors * graph = new TGraphAsymmErrors(num, denom, "pois");
  styleGraph(graph);
  graph->SetName(graphName);

  graph->SetMinimum(-0.1); graph->GetXaxis()->SetRangeUser(50,250);
  graph->SetLineColor(kBlack); graph->SetMarkerColor(kBlack);

  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  styleCanvas(can_h);

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
  num->DrawNormalized();
  denom->DrawNormalized("same");

  TLegend* legend = new TLegend(0.5,0.65,0.95,0.88,NULL,"brNDC") ;
  legend->SetBorderSize(0);
  if (tagType == "Single"){
    // legend->AddEntry(num, "1H SB (numerator)", "PL");
    // legend->AddEntry(denom, "0H SB (denominator)", "PL");
    legend->AddEntry(num, "1H SR+SB (numerator)", "PL");
    legend->AddEntry(denom, "0H SR+SB (denominator)", "PL");
    graph->SetTitle(";Soft-drop mass [GeV];R(1H/0H)");
  }
  else if (tagType == "Double"){
    // legend->AddEntry(num, "2H SB (numerator)", "PL");
    // legend->AddEntry(denom, "0H SB (denominator)", "PL");
    legend->AddEntry(num, "2H SR+SB (numerator)", "PL");
    legend->AddEntry(denom, "0H SR+SB (denominator)", "PL");
    graph->SetTitle(";Soft-drop mass [GeV];R(2H/0H)");
  }

  legend->Draw();

  pad2->cd();
  graph->Draw("APE");
  graph->GetYaxis()->SetRangeUser(-0.1,0.5);
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

  cdRPFFinal->cd();
  can_h->Write(thisType);

  cdRPF->cd();
  num->Write("Num_"+thisType);
  denom->Write("Denom_"+thisType);
  graph->Write("RPFRatio_"+thisType);
  delete can_h;
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

  TString graphName = "";
  if (which == "doubleB") graphName = "DoubleB_Stack";
  else if (which == "tau32") graphName = "Tau32_Stack";
  else if (which == "leadmass") graphName = "LeadMass_Stack";
  else if (which == "subleadmass") graphName = "SubLeadMass_Stack";
  else if (which == "leadpt") graphName = "LeadPt_Stack";
  else if (which == "subleadpt") graphName = "SubLeadPt_Stack";

  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  styleCanvas(can_h);
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
  delete can_h;
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

  TString graphName = which+"_Stack";
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  styleCanvas(can_h);
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
  delete can_h;
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


  TString graphName = "METStack_"+which;
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  styleCanvas(can_h);
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
  delete can_h;
}

void makeMETNorm(vector<TH1F*> dem_histos, TString tagType = "") {
  //dem_histos will be [SR,0HSB,0HSBCuts,0HSR,0HSRCuts]
  TH1D * h_SR   = (TH1D*)dem_histos[0]->Clone("h_SR");
  TH1D * h_0HSB_withCuts = (TH1D*)dem_histos[2]->Clone("h_0HSB_withCuts");
  TH1D * h_0HSR_withCuts = (TH1D*)dem_histos[4]->Clone("h_0HSR_withCuts");

  h_0HSB_withCuts->Add(h_0HSR_withCuts);

  h_SR->SetMarkerStyle(20); h_SR->SetMarkerSize(0.85); h_SR->SetLineColor(kRed); h_SR->SetMarkerColor(kRed);
  h_0HSB_withCuts->SetMarkerStyle(20); h_0HSB_withCuts->SetMarkerSize(0.85); h_0HSB_withCuts->SetLineColor(kBlack); h_0HSB_withCuts->SetMarkerColor(kBlack);

  if (tagType=="Single" && whichRegion=="signal") h_SR->SetTitle("0l, BTagsM>0;MET [GeV];Normalized to area");
  else if (tagType=="Double" && whichRegion=="signal") h_SR->SetTitle("0l, BTagsM>1;MET [GeV];Normalized to area");
  else if (tagType=="Single" && whichRegion == "singleLept") h_SR->SetTitle("1l, BTagsM>0;MET [GeV];Normalized to area");
  else if (tagType=="Double" && whichRegion=="singleLept") h_SR->SetTitle("1l, BTagsM>1;MET [GeV];Normalized to area");

  TString graphName = "MET";
  TGraphAsymmErrors * graph = new TGraphAsymmErrors(h_SR, h_0HSB_withCuts, "pois");
  styleGraph(graph);
  graph->SetMinimum(-0.1); graph->GetXaxis()->SetRangeUser(300,1400);
  graph->SetTitle(";MET [GeV];2HSR / 0HCuts");
  if (tagType=="Single") graph->SetTitle(";MET [GeV];1HSR / 0HCuts");

  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  styleCanvas(can_h);

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
  h_SR->GetYaxis()->SetTitleOffset(0.8);

  //try normalizing to 1?
  h_SR->Scale(1/h_SR->Integral());
  h_0HSB_withCuts->Scale(1/h_0HSB_withCuts->Integral());

  h_SR->Draw();
  h_0HSB_withCuts->Draw("same");

  TLegend* legend = new TLegend(0.5,0.65,0.95,0.88,NULL,"brNDC") ;
  legend->SetBorderSize(0);
  TString legendEntryName = tagType+" SR";
  if (tagType=="Single") {
    legend->AddEntry(h_SR, "1H SR", "PL");
    legend->AddEntry(h_0HSB_withCuts, "0H SB+SR,BTagsM>0", "PL");
  }
  else if (tagType=="Double") {
    legend->AddEntry(h_SR, "2H SR", "PL");
    legend->AddEntry(h_0HSB_withCuts, "0H SB+SR,BTagsM>1", "PL");
  }
  legend->Draw();

  pad2->cd();
  graph->Draw("APE");
  graph->GetYaxis()->SetTitleOffset(0.25); graph->GetYaxis()->SetTitleSize(0.13);
  graph->GetXaxis()->SetTitleOffset(0.9); graph->GetXaxis()->SetTitleSize(0.14);
  graph->GetXaxis()->SetLabelSize(0.14); graph->GetYaxis()->SetLabelSize(0.10);
  graph->GetYaxis()->SetNdivisions(505);
  graph->Draw("APE");


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
  delete can_h;
}

void makeMETNormCompare(vector<TH1F*> dem_histos, TString tagType = "") {
  //dem_histos will be [SR,0HSB,0HSBCuts,0HSR,0HSRCuts]
  TH1D * h_SR   = (TH1D*)dem_histos[0]->Clone("h_SR");
  TH1D * h_0HSB_noCuts = (TH1D*)dem_histos[1]->Clone("h_0HSB_noCuts");
  TH1D * h_0HSB_withCuts = (TH1D*)dem_histos[2]->Clone("h_0HSB_withCuts");
  TH1D * h_0HSR_noCuts = (TH1D*)dem_histos[3]->Clone("h_0HSR_noCuts");
  TH1D * h_0HSR_withCuts = (TH1D*)dem_histos[4]->Clone("h_0HSR_withCuts");

  //Adding SR and SB together
  h_0HSB_noCuts->Add(h_0HSR_noCuts);
  h_0HSB_withCuts->Add(h_0HSR_withCuts);

  h_SR->SetMarkerStyle(20); h_SR->SetMarkerSize(0.85); h_SR->SetLineColor(kRed); h_SR->SetMarkerColor(kRed);
  h_0HSB_withCuts->SetMarkerStyle(20); h_0HSB_withCuts->SetMarkerSize(0.85); h_0HSB_withCuts->SetLineColor(kBlack); h_0HSB_withCuts->SetMarkerColor(kBlack);
  h_0HSB_noCuts->SetMarkerStyle(20); h_0HSB_noCuts->SetMarkerSize(0.85); h_0HSB_noCuts->SetLineColor(kBlue); h_0HSB_noCuts->SetMarkerColor(kBlue);

  if (tagType=="Single") h_SR->SetTitle("BTagsM>0;MET [GeV];Normalized to area");
  else h_SR->SetTitle("BTagsM>1;MET [GeV];Normalized to area");

  if (whichRegion == "singleLept" && tagType=="Single") h_SR->SetTitle("1-l, BTagsM>0;MET [GeV];Normalized to area");
  else if (whichRegion == "singleLept" && tagType=="Double") h_SR->SetTitle("1-l, BTagsM>1;MET [GeV];Normalized to area");

  TString graphName = "MET";
  TGraphAsymmErrors * graph = new TGraphAsymmErrors(h_SR, h_0HSB_withCuts, "pois");
  styleGraph(graph);
  graph->SetMinimum(-0.1); graph->GetXaxis()->SetRangeUser(300,1400);
  graph->SetTitle(";MET [GeV];2HSR / 0HCuts");

  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  styleCanvas(can_h);

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

  //try normalizing to 1?
  h_SR->Scale(1/h_SR->Integral());
  h_0HSB_withCuts->Scale(1/h_0HSB_withCuts->Integral());
  h_SR->Draw();
  h_0HSB_withCuts->Draw("same");

  TLegend* legend = new TLegend(0.5,0.65,0.95,0.88,NULL,"brNDC") ;
  legend->SetBorderSize(0);
  TString legendEntryName = tagType+" SR";
  legend->AddEntry(h_SR, legendEntryName, "PL");
  legend->AddEntry(h_0HSB_noCuts, "0H SB+SR, no cuts", "PL");
  if (tagType=="Single") legend->AddEntry(h_0HSB_withCuts, "0H SB+SR,BTagsM>0", "PL");
  else if (tagType=="Double") legend->AddEntry(h_0HSB_withCuts, "0H SB+SR,BTagsM>1", "PL");
  legend->Draw();

  pad2->cd();
  graph->Draw("APE");
  graph->GetYaxis()->SetTitleOffset(0.25); graph->GetYaxis()->SetTitleSize(0.13);
  graph->GetXaxis()->SetTitleOffset(0.9); graph->GetXaxis()->SetTitleSize(0.14);
  graph->GetXaxis()->SetLabelSize(0.14); graph->GetYaxis()->SetLabelSize(0.10);
  graph->GetYaxis()->SetNdivisions(505);
  graph->Draw("APE");

  // TGraphAsymmErrors * graph = new TGraphAsymmErrors(h_SR, h_0HSB_withCuts, "pois");
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
  delete can_h;
}

void tableOfYields(vector<TH1F*> dem_histos, TString bkgType, TString tagType){  //A, B, C, D
  TH1F * histo_A = (TH1F*)dem_histos[0]->Clone("histo_A");
  TH1F * histo_B = (TH1F*)dem_histos[1]->Clone("histo_B");
  TH1F * histo_C = (TH1F*)dem_histos[2]->Clone("histo_C");
  TH1F * histo_D = (TH1F*)dem_histos[3]->Clone("histo_D");

  double a_int; double a_error; a_int = histo_A->IntegralAndError(0,5,a_error,"");
  double b_int; double b_error; b_int = histo_B->IntegralAndError(0,5,b_error,"");
  double c_int; double c_error; c_int = histo_C->IntegralAndError(0,5,c_error,"");
  double d_int; double d_error; d_int = histo_D->IntegralAndError(0,5,d_error,"");

  TH1F *histo_pred = (TH1F*)histo_B->Clone("histo_pred");
  histo_pred->Multiply(histo_C);histo_pred->Divide(histo_D);

  double bkgNorm; double bkgNorm_error;
  bkgNorm = histo_pred->IntegralAndError(0,5,bkgNorm_error,"");

  double kappa = a_int/bkgNorm;
  double kappa_error = sqrt(a_error*a_error + bkgNorm_error*bkgNorm_error)/a_int;

  TString yieldsFileName;
  if (tagType=="Double") yieldsFileName = "Yields/2HYields_"+whichRegion+"_"+bkgType+".txt";
  else if (tagType=="Single") yieldsFileName = "Yields/1HYields_"+whichRegion+"_"+bkgType+".txt";

  ofstream yields;
  yields.open(yieldsFileName);

  yields<<"% "<< bkgType <<endl;
  yields <<"\\documentclass[11pt, oneside]{article}\n\n";
  yields<<"\\usepackage{numprint}\n"<<"\\begin{document}\n"<<"\\begin{table}\n";
  yields<<"\\npdecimalsign{.}\n"<<"\\nprounddigits{1}\n"<<"\\centering\n";
  yields<<"\\begin{tabular}{ |c|n{4}{1}n{2}{3}| }\n"<<"\\hline\n";
  yields<<"Region &\\multicolumn{2}{c|}{MC Yields} \\\\ [0.5ex]  \\hline \\hline\n";
  if (tagType=="Double") {
    yields<<"2H SR & "<< a_int <<" & \\pm"<< a_error <<"\\\\ \n";
    yields<<"2H SB & "<< b_int <<" & \\pm"<< b_error <<"\\\\ \n";
  }
  else if (tagType=="Single") {
    yields<<"1H SR & "<< a_int <<" & \\pm"<< a_error <<"\\\\ \n";
    yields<<"1H SB & "<< b_int <<" & \\pm"<< b_error <<"\\\\ \n";
  }
  yields<<"0H SR & "<< c_int <<" & \\pm"<< c_error <<"\\\\ \n";
  yields<<"0H SB & "<< d_int <<" & \\pm"<< d_error <<"\\\\ \n";
  yields<<"Prediction & "<< bkgNorm <<" & \\pm"<< bkgNorm_error <<"\\\\ \\hline\n";
  yields<<"\\hline\n";
  yields<<"Kappa & "<< kappa << " & \\pm"<< kappa_error <<"\\\\ \\hline\n";
  yields<<"\\end{tabular}\n"<<"\\npnoround\n"<<"\\end{table}\n"<<"\\end{document}\n";
  yields.close();
}

void tableOfMETNorm(vector<TH1F*> dem_histos, TString bkgType, TString tagType){  //A, D, DOpt, C, COpt
  TH1F * histo_A = (TH1F*)dem_histos[0]->Clone("histo_A");
  TH1F * histo_DNoCuts = (TH1F*)dem_histos[1]->Clone("histo_DNoCuts");
  TH1F * histo_DWithCuts = (TH1F*)dem_histos[2]->Clone("histo_DWithCuts");
  TH1F * histo_CNoCuts = (TH1F*)dem_histos[3]->Clone("histo_CNoCuts");
  TH1F * histo_CWithCuts = (TH1F*)dem_histos[4]->Clone("histo_CWithCuts");

  TH1F *h_2H_Norm = (TH1F*) histo_A->Clone("h_2H_Norm");
  h_2H_Norm->Scale(1/h_2H_Norm->Integral(1,4));
  TH1F *h_0H_Norm = (TH1F*) histo_DWithCuts->Clone("h_0H_Norm");
  h_0H_Norm->Add(histo_CWithCuts);
  h_0H_Norm->Scale(1/h_0H_Norm->Integral(1,4));

  histo_DWithCuts->Add(histo_CWithCuts);

  TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo_A, histo_DWithCuts, "pois");
  TF1*f0=new TF1("f0","pol0", 300,1400);
  graph->Fit("f0","Q","R",300,1400);
  double p0 = f0->GetParameter(0); //this does, indeed, work
  double error = f0->GetParError(0);
  // graph->GetFunction("f0")->SetBit(TF1::kNotDraw);

  double bin1Obs = histo_A->GetBinContent(1)/histo_DWithCuts->GetBinContent(1);
  double bin2Obs = histo_A->GetBinContent(2)/histo_DWithCuts->GetBinContent(2);
  double bin3Obs = histo_A->GetBinContent(3)/histo_DWithCuts->GetBinContent(3);

  double bin1Dev = 1.0 + (bin1Obs-p0)/p0;
  double bin2Dev = 1.0 + (bin2Obs-p0)/p0;
  double bin3Dev = 1.0 + (bin3Obs-p0)/p0;


  TString yieldsFileName;
  if (tagType=="Double") yieldsFileName = "Yields/2HMETShape_"+whichRegion+"_"+bkgType+".txt";
  else if (tagType=="Single") yieldsFileName = "Yields/1HMETShape_"+whichRegion+"_"+bkgType+".txt";

  ofstream yields;
  yields.open(yieldsFileName);

  yields<<"% Background: "<< bkgType <<endl;
  yields<<"% Region: "<< whichRegion <<endl;
  yields<<"% Best fit: "<< p0 <<" +/- "<<error <<endl;

  yields <<"\\documentclass[12pt]{article}\n\n";
  yields<<"\\usepackage{numprint}\n"<<"\\begin{document}\n"<<"\\begin{table}\n";
  yields<<"\\npdecimalsign{.}\n"<<"\\nprounddigits{3}\n"<<"\\centering\n";
  yields<<"\\begin{tabular}{ c|n{4}{3}n{2}{3}|n{3}{3}n{2}{3}|n{4}{3} }\n"<<"\\hline\n";
  yields<<"MET region & \\multicolumn{2}{c|}{MC Yields} & \\multicolumn{2}{c|}{Bkg Frac.} & \\multicolumn{1}{c}{Deviation}  \\\\ \n";
  yields<<"\\hline \\hline\n";

  if (tagType=="Double") yields<<"\\multicolumn{6}{c}{2H SR} \\\\ \n";
  else if (tagType=="Single") yields<<"\\multicolumn{6}{c}{1H SR} \\\\ \n";
  yields<<"\\hline\n";

  yields<<"$[300,500]$ GeV & "<<  histo_A->GetBinContent(1)<<"& \\pm "<< histo_A->GetBinError(1)<< "& " << h_2H_Norm->GetBinContent(1)<<" & \\pm "<<h_2H_Norm->GetBinError(1)<<" & "<< " \\\\ \n";
  yields<<"$[500,700]$ GeV & "<<  histo_A->GetBinContent(2)<<"& \\pm "<< histo_A->GetBinError(2)<< "& " << h_2H_Norm->GetBinContent(2)<<" & \\pm "<<h_2H_Norm->GetBinError(2)<<" & "<< " \\\\ \n";
  yields<<"$[700+]$ GeV & "<<  histo_A->GetBinContent(3)<<"& \\pm "<< histo_A->GetBinError(3)<< "& " << h_2H_Norm->GetBinContent(3)<<" & \\pm "<<h_2H_Norm->GetBinError(3)<<" & "<< " \\\\ \n";

  yields<<"\\hline\n";
  if (tagType=="Double") yields<<"\\multicolumn{6}{c}{0H, $B_{M}>$1} \\\\ \n";
  else if (tagType=="Single") yields<<"\\multicolumn{6}{c}{0H, $B_{M}>$0} \\\\ \n";
  yields<<"\\hline\n";
  yields<<"$[300,500]$ GeV & "<<  histo_DWithCuts->GetBinContent(1)<<" & \\pm "<< histo_DWithCuts->GetBinError(1)<< "& " << h_0H_Norm->GetBinContent(1)<<" & \\pm "<<h_0H_Norm->GetBinError(1)<<" & "<< bin1Dev<<" \\\\ \n";
  yields<<"$[500,700]$ GeV & "<<  histo_DWithCuts->GetBinContent(2)<<" & \\pm "<< histo_DWithCuts->GetBinError(2)<< "& " << h_0H_Norm->GetBinContent(2)<<" & \\pm "<<h_0H_Norm->GetBinError(2)<<" & "<< bin2Dev<<" \\\\ \n";
  yields<<"$[700+]$ GeV & "<<  histo_DWithCuts->GetBinContent(3)<<" & \\pm "<< histo_DWithCuts->GetBinError(3)<< "& " << h_0H_Norm->GetBinContent(3)<<" & \\pm "<<h_0H_Norm->GetBinError(3)<<" & "<< bin3Dev<< " \\\\ \n";
  yields<<"\\hline\n";


  yields<<"\\end{tabular}\n"<<"\\npnoround\n"<<"\\end{table}\n"<<"\\end{document}\n";
  yields.close();
}


void pieChart(vector<TH1F*> h_QCD, vector<TH1F*> h_WJets, vector<TH1F*> h_ZJets, vector<TH1F*> h_TT, TString regionLabel){
  TH1F * h_Q = (TH1F*)h_QCD[0]->Clone("h_Q"); for (int i=1;i<h_QCD.size();i++) h_Q->Add(h_QCD[i]);
  TH1F * h_W = (TH1F*)h_WJets[0]->Clone("h_W"); for (int i=1;i<h_WJets.size();i++) h_W->Add(h_WJets[i]);
  TH1F * h_Z = (TH1F*)h_ZJets[0]->Clone("h_Z"); for (int i=1;i<h_ZJets.size();i++) h_Z->Add(h_ZJets[i]);
  TH1F * h_T = (TH1F*)h_TT[0]->Clone("h_T");  for (int i=1;i<h_TT.size();i++) h_T->Add(h_TT[i]);

  double Q_yields = h_Q->Integral();
  double W_yields = h_W->Integral();
  double Z_yields = h_Z->Integral();
  double T_yields = h_T->Integral();

  TString yieldsFileName;
  yieldsFileName = whichRegion+", "+regionLabel;
  TCanvas * can_h = new TCanvas(yieldsFileName,yieldsFileName, 50, 50, 1200, 1200);

  double vals[] = {Q_yields,W_yields,Z_yields,T_yields};
  int colors[] = {kSpring-5,kAzure+1,kOrange-2,kGreen+3};
  const char *labels[] = {"QCD","WJets","ZJets","TT"};
  int nvals = sizeof(vals)/sizeof(vals[0]);
  TPie *pie1 = new TPie("pie1", yieldsFileName,nvals,vals,colors, labels);
  pie1->SetAngularOffset(355.);
  pie1->SetAngle3D(45.);
  pie1->SetHeight(0.04);
  pie1->SetCircle(0.5,0.5,.25);
  pie1->SetLabelsOffset(.03);
  pie1->SetLabelFormat("%perc");
  pie1->Draw("3d");

  TString savename = whichRegion+"_"+regionLabel;
  cdPie->cd();
  if (regionLabel=="2H"){
    TCanvas * can2 = new TCanvas("PieChartLegend","PieChartLegend", 50, 50, 1200, 1200);
    TLegend *leg = pie1->MakeLegend();
    leg->SetBorderSize(0);
    leg->Draw();
    can2->Write("Legend");
    delete can2;
  }
  can_h->Write(savename);

  delete pie1;
  delete can_h;
}

void pieChart1l(vector<TH1F*> h_WJets, vector<TH1F*> h_TT, TString regionLabel){
  TH1F * h_W = (TH1F*)h_WJets[0]->Clone("h_W"); for (int i=1;i<h_WJets.size();i++) h_W->Add(h_WJets[i]);
  TH1F * h_T = (TH1F*)h_TT[0]->Clone("h_T");  for (int i=1;i<h_TT.size();i++) h_T->Add(h_TT[i]);

  double W_yields = h_W->Integral();
  double T_yields = h_T->Integral();

  TString yieldsFileName;
  yieldsFileName = whichRegion+", "+regionLabel;
  TCanvas * can_h = new TCanvas(yieldsFileName,yieldsFileName, 50, 50, 1200, 1200);

  double vals[] = {W_yields,T_yields};
  int colors[] = {kAzure+1,kGreen+3};
  const char *labels[] = {"WJets","TT"};
  int nvals = sizeof(vals)/sizeof(vals[0]);
  TPie *pie1 = new TPie("pie1", yieldsFileName,nvals,vals,colors, labels);
  pie1->SetAngularOffset(355.);
  pie1->SetAngle3D(45.);
  pie1->SetHeight(0.04);
  pie1->SetCircle(0.5,0.5,.25);
  pie1->SetLabelsOffset(.03);
  pie1->SetLabelFormat("%perc");
  pie1->Draw("3d");

  TString savename = whichRegion+"_"+regionLabel;
  cdPie->cd();
  if (regionLabel=="2H"){
    TLegend *leg = pie1->MakeLegend();
    TCanvas * can2 = new TCanvas("PieChartLegend","PieChartLegend", 50, 50, 1200, 1200);
    leg->SetBorderSize(0);
    leg->Draw();
    can2->Write("Legend");
    delete can2;
  }
  can_h->Write(savename);

  delete pie1;
  delete can_h;
}


TH1F *make0EventUncSum(vector<TH1F*> dem_histos) {
  //Histograms should be in this order: QCD, TT, WJets,ZJets
  TH1F * h_QCD   = (TH1F*)dem_histos[0]->Clone("h_QCD");
  TH1F * h_TT    = (TH1F*)dem_histos[1]->Clone("h_TT");
  TH1F * h_WJets = (TH1F*)dem_histos[2]->Clone("h_WJets");
  TH1F * h_ZJets = (TH1F*)dem_histos[3]->Clone("h_ZJets");

  for (int i=1;i<=3;i++){
    // if (h_QCD->GetBinContent(i)<0.001) h_QCD->SetBinError(i,19.55);
    if (h_TT->GetBinContent(i)<0.001) h_TT->SetBinError(i,0.795);
    if (h_WJets->GetBinContent(i)<0.001) h_WJets->SetBinError(i,0.210);
    if (h_ZJets->GetBinContent(i)<0.001) h_ZJets->SetBinError(i,0.328);
  }
  TH1F *h_sum = (TH1F*)h_QCD->Clone("h_sum");
  h_sum->Add(h_TT);
  h_sum->Add(h_WJets);
  h_sum->Add(h_ZJets);
  return h_sum;
}

TH1F *make0EventUncSum_1l(vector<TH1F*> dem_histos) {
  //Histograms should be in this order: TT, WJets
  TH1F * h_TT    = (TH1F*)dem_histos[0]->Clone("h_TT");
  TH1F * h_WJets = (TH1F*)dem_histos[1]->Clone("h_WJets");

  for (int i=1;i<=3;i++){
    if (h_TT->GetBinContent(i)<0.001) h_TT->SetBinError(i,0.795);
    if (h_WJets->GetBinContent(i)<0.001) h_WJets->SetBinError(i,0.210);
  }
  TH1F *h_sum = (TH1F*)h_TT->Clone("h_sum");
  h_sum->Add(h_WJets);
  return h_sum;
}

void styleCanvas(TCanvas * can_h){
  can_h->SetFillColor(0); can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W ); can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H ); can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0); can_h->SetTicky(0);
}

void styleGraph(TGraph *graph){
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
}
