#include <TH1F.h>
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
void makeStackPlot(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets,TH1F* h_SnglT, TH1F* h_T5HH1600,TH1F* h_T5HH2000,TH1F* h_T5HH2200,TH1F* h_TChiHH600,TH1F* h_TChiHH800,TH1F* h_TChiHH1000,TString which);
void makeSingleLeptStackPlot(vector<TH1F*> h_TT,vector<TH1F*> h_WJets,vector<TH1F*> h_SnglT,TString which);
void makeMETStack(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets, TH1F* h_SnglT, TH1F* h_T5HH1600,TH1F* h_T5HH2000,TH1F* h_T5HH2200,TH1F* h_TChiHH600,TH1F* h_TChiHH800,TH1F* h_TChiHH1000,TString tagType);
void makeFullBkgClosure(vector<TH1F*> dem_histos, TString bkgType, TString tagType, TString year);
void makePureBkgClosure(vector<TH1F*> dem_histos, TString bkgType, TString tagType, TString year);
void makeMETNorm(vector<TH1F*> dem_histos, TString bkgType, TString tagType);
void makeMETNormCompare(vector<TH1F*> dem_histos, TString bkgType, TString tagType);
TH1F* make0EventUncSum(vector<TH1F*> dem_histos); //in this order: QCD, TT, WJets,ZJets,SnglT
TH1F* make0EventUncSum_1l(vector<TH1F*> dem_histos); //in this order: SnglT, TT, WJets
void styleCanvas(TCanvas *can_h);
void styleGraph(TGraph *graph);
void DrawOverflow(TH1F* &h);
void tableOfYields(vector<TH1F*> dem_histos, TString bkgType, TString tagType); //A, B, C, D
void tableOfMETNorm(vector<TH1F*> dem_histos, TString bkgType, TString tagType); //A, D, DOpt, C, COpt
void pieChart(vector<TH1F*> h_QCD, vector<TH1F*> h_WJets, vector<TH1F*> h_ZJets, vector<TH1F*> h_TT, vector<TH1F*> h_SnglT, TString regionLabel, TString bin);
void pieChart1l(vector<TH1F*> h_WJets, vector<TH1F*> h_TT, vector<TH1F*> h_SnglT, TString regionLabel, TString bin);
void pieChartPhoton(vector<TH1F*> h_GJets, vector<TH1F*> h_QCD, TString regionLabel, TString bin);
void tableOfYieldsPure(vector<TH1F*> dem_histos, TString bkgType, TString tagType); //A, B, C, D
void massCorrelations(vector<TH1F*> dem_histos, TString bkgType); //A, B, C, D



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
bool runFullBkg = false;
bool runRPFPlots = true;
bool runStacks = false;
bool runMETNorm = true;
bool runTableOfYields = true;
bool runPies = false;
bool showNoCutLine = false;
bool runPure = false;
bool runDM = false;
bool runMassCorrelations = true;

bool savePDFs = true;

// string whichRegion = "signal";
// string whichRegion = "singleLept";
string whichRegion = "photon";

TString year = "Run2";

// TFile * f = TFile::Open("ALPHABET_V18_2BoostedH.root");
TFile * f = TFile::Open("/eos/uscms/store/user/emacdona/boostedHiggsPlusMET/ALPHABET_0lLoose_noPU.root"); //to check if photon stats look OK



TFile * fSignal;// = TFile::Open("ALPHABETBoost_V17_signalOnly.root");
// TFile * fPhoton = TFile::Open("ALPHABET_V18_photon.root");
TFile * fPhoton = TFile::Open("/eos/uscms/store/user/emacdona/boostedHiggsPlusMET/ALPHABET_photonTight_noPU.root");


TFile * fSingleLept = TFile::Open("ALPHABET_V18_1l_2BoostedH.root");

// TString fout_name = "ABCDPlots_"+whichRegion+"_test.root";
TString fout_name = "ABCDPlots_"+whichRegion+"_testLosoeCuts.root";
TFile* fout = new TFile(fout_name,"recreate");


TDirectory *cdABCD  = fout->mkdir("ABCD");
TDirectory *cdClose = fout->mkdir("Closure");
TDirectory *cdAssum = fout->mkdir("Assumption");
TDirectory *cdRPFFinal  = fout->mkdir("RPF");
TDirectory *cdRPF  = fout->mkdir("RPF_Support");
TDirectory *cdPie  = fout->mkdir("PieCharts");
TDirectory *cdOther  = fout->mkdir("OtherPlots");

void runABCD() {
  TH1F * h_avgM_A_sum; TH1F * h_avgM_B_sum; TH1F * h_avgM_A1_sum; TH1F * h_avgM_B1_sum; TH1F * h_avgM_C_sum; TH1F * h_avgM_D_sum;
  TH1F * h_A_data; TH1F * h_B_data; TH1F * h_A1_data; TH1F * h_B1_data; TH1F *h_C_data; TH1F *h_D_data;

  TH1F * h_A_sum; TH1F * h_B_sum; TH1F * h_A1_sum; TH1F * h_B1_sum; TH1F *h_C_sum; TH1F *h_D_sum;
  TH1F * hP_A_sum; TH1F * hP_B_sum; TH1F * hP_A1_sum; TH1F * hP_B1_sum; TH1F *hP_C_sum; TH1F *hP_D_sum;
  TH1F * h_A5_sum; TH1F * h_B5_sum; TH1F * h_A15_sum; TH1F * h_B15_sum; TH1F *h_C5_sum; TH1F *h_D5_sum;
  TH1F * h_A_QCD; TH1F * h_B_QCD; TH1F * h_A1_QCD; TH1F * h_B1_QCD; TH1F * h_C_QCD; TH1F * h_D_QCD;
  TH1F * h_A_SnglT; TH1F * h_B_SnglT; TH1F * h_A1_SnglT; TH1F * h_B1_SnglT; TH1F * h_C_SnglT; TH1F * h_D_SnglT;
  TH1F * h_A_GJets; TH1F * h_B_GJets; TH1F * h_A1_GJets; TH1F * h_B1_GJets; TH1F * h_C_GJets; TH1F * h_D_GJets;
  TH1F * h_A_TT; TH1F * h_B_TT; TH1F * h_A1_TT; TH1F * h_B1_TT; TH1F * h_C_TT; TH1F * h_D_TT;
  TH1F * h_A_TT_Di; TH1F * h_B_TT_Di; TH1F * h_A1_TT_Di; TH1F * h_B1_TT_Di; TH1F * h_C_TT_Di; TH1F * h_D_TT_Di;
  TH1F * h_A_TT_SL; TH1F * h_B_TT_SL; TH1F * h_A1_TT_SL; TH1F * h_B1_TT_SL; TH1F * h_C_TT_SL; TH1F * h_D_TT_SL;
  TH1F * h_A_WJets; TH1F * h_B_WJets; TH1F * h_A1_WJets; TH1F * h_B1_WJets; TH1F * h_C_WJets; TH1F * h_D_WJets;
  TH1F * h_A_ZJets; TH1F * h_B_ZJets; TH1F * h_A1_ZJets; TH1F * h_B1_ZJets; TH1F * h_C_ZJets; TH1F * h_D_ZJets;

  //Checking low vs high PU
  TH1F * h_AlowPU_sum; TH1F * h_BlowPU_sum; TH1F * h_A1lowPU_sum; TH1F * h_B1lowPU_sum; TH1F *h_ClowPU_sum; TH1F *h_DlowPU_sum;
  TH1F * h_AhighPU_sum; TH1F * h_BhighPU_sum; TH1F * h_A1highPU_sum; TH1F * h_B1highPU_sum; TH1F *h_ChighPU_sum; TH1F *h_DhighPU_sum;


  TH1F * h_A_TChiHH200; TH1F * h_B_TChiHH200; TH1F * h_A1_TChiHH200; TH1F * h_B1_TChiHH200; TH1F * h_C_TChiHH200; TH1F * h_D_TChiHH200;
  TH1F * h_A_TChiHH400; TH1F * h_B_TChiHH400; TH1F * h_A1_TChiHH400; TH1F * h_B1_TChiHH400; TH1F * h_C_TChiHH400; TH1F * h_D_TChiHH400;
  TH1F * h_A_TChiHH700; TH1F * h_B_TChiHH700; TH1F * h_A1_TChiHH700; TH1F * h_B1_TChiHH700; TH1F * h_C_TChiHH700; TH1F * h_D_TChiHH700;
  TH1F * h_A_TChiHH1000; TH1F * h_B_TChiHH1000; TH1F * h_A1_TChiHH1000; TH1F * h_B1_TChiHH1000; TH1F * h_C_TChiHH1000; TH1F * h_D_TChiHH1000;
  TH1F * h_A_TChiHH1300; TH1F * h_B_TChiHH1300; TH1F * h_A1_TChiHH1300; TH1F * h_B1_TChiHH1300; TH1F * h_C_TChiHH1300; TH1F * h_D_TChiHH1300;
  TH1F * h_A_TChiHH1500; TH1F * h_B_TChiHH1500; TH1F * h_A1_TChiHH1500; TH1F * h_B1_TChiHH1500; TH1F * h_C_TChiHH1500; TH1F * h_D_TChiHH1500;

  TH1F * h_J1M_doubletagSR_sum; TH1F * h_J1M_doubletagSB_sum; TH1F * h_J1M_tagSR_sum; TH1F * h_J1M_tagSB_sum; TH1F * h_J1M_antitagSR_sum; TH1F * h_J1M_antitagSB_sum;
  TH1F * h_J2M_doubletagSR_sum; TH1F * h_J2M_doubletagSB_sum; TH1F * h_J2M_tagSR_sum; TH1F * h_J2M_tagSB_sum; TH1F * h_J2M_antitagSR_sum; TH1F * h_J2M_antitagSB_sum;
  TH1F * h_J1J2M_doubletagSR_sum; TH1F * h_J1J2M_doubletagSB_sum; TH1F * h_J1J2M_tagSR_sum; TH1F * h_J1J2M_tagSB_sum;
  TH1F * h_J1J2M_antitagSR_sum; TH1F * h_J1J2M_antitagSB_sum;

  TH1F * h_J1M_mjBins_doubletagSR_sum; TH1F * h_J1M_mjBins_doubletagSB_sum; TH1F * h_J1M_mjBins_tagSR_sum; TH1F * h_J1M_mjBins_tagSB_sum; TH1F * h_J1M_mjBins_antitagSR_sum; TH1F * h_J1M_mjBins_antitagSB_sum;
  TH1F * h_J2M_mjBins_doubletagSR_sum; TH1F * h_J2M_mjBins_doubletagSB_sum; TH1F * h_J2M_mjBins_tagSR_sum; TH1F * h_J2M_mjBins_tagSB_sum; TH1F * h_J2M_mjBins_antitagSR_sum; TH1F * h_J2M_mjBins_antitagSB_sum;

  TH1F * h_baseline_j1mass_TT; TH1F * h_J1M_doubletagSR_TT; TH1F * h_J1M_doubletagSB_TT; TH1F * h_J1M_tagSR_TT; TH1F * h_J1M_tagSB_TT; TH1F * h_J1M_antitagSR_TT; TH1F * h_J1M_antitagSB_TT;
  TH1F * h_baseline_j2mass_TT; TH1F * h_J2M_doubletagSR_TT; TH1F * h_J2M_doubletagSB_TT; TH1F * h_J2M_tagSR_TT; TH1F * h_J2M_tagSB_TT; TH1F * h_J2M_antitagSR_TT; TH1F * h_J2M_antitagSB_TT;
  TH1F * h_J1M_mjBins_doubletagSR_TT; TH1F * h_J1M_mjBins_doubletagSB_TT; TH1F * h_J1M_mjBins_tagSR_TT; TH1F * h_J1M_mjBins_tagSB_TT; TH1F * h_J1M_mjBins_antitagSR_TT; TH1F * h_J1M_mjBins_antitagSB_TT;
  TH1F * h_J2M_mjBins_doubletagSR_TT; TH1F * h_J2M_mjBins_doubletagSB_TT; TH1F * h_J2M_mjBins_tagSR_TT; TH1F * h_J2M_mjBins_tagSB_TT; TH1F * h_J2M_mjBins_antitagSR_TT; TH1F * h_J2M_mjBins_antitagSB_TT;

  TH1F * h_J1M_doubletagSR_TT_Di; TH1F * h_J1M_doubletagSB_TT_Di; TH1F * h_J1M_tagSR_TT_Di; TH1F * h_J1M_tagSB_TT_Di; TH1F * h_J1M_antitagSR_TT_Di; TH1F * h_J1M_antitagSB_TT_Di;
  TH1F * h_J2M_doubletagSR_TT_Di; TH1F * h_J2M_doubletagSB_TT_Di; TH1F * h_J2M_tagSR_TT_Di; TH1F * h_J2M_tagSB_TT_Di; TH1F * h_J2M_antitagSR_TT_Di; TH1F * h_J2M_antitagSB_TT_Di;
  TH1F * h_J1M_mjBins_doubletagSR_TT_Di; TH1F * h_J1M_mjBins_doubletagSB_TT_Di; TH1F * h_J1M_mjBins_tagSR_TT_Di; TH1F * h_J1M_mjBins_tagSB_TT_Di; TH1F * h_J1M_mjBins_antitagSR_TT_Di; TH1F * h_J1M_mjBins_antitagSB_TT_Di;
  TH1F * h_J2M_mjBins_doubletagSR_TT_Di; TH1F * h_J2M_mjBins_doubletagSB_TT_Di; TH1F * h_J2M_mjBins_tagSR_TT_Di; TH1F * h_J2M_mjBins_tagSB_TT_Di; TH1F * h_J2M_mjBins_antitagSR_TT_Di; TH1F * h_J2M_mjBins_antitagSB_TT_Di;

  TH1F * h_J1M_doubletagSR_TT_SL; TH1F * h_J1M_doubletagSB_TT_SL; TH1F * h_J1M_tagSR_TT_SL; TH1F * h_J1M_tagSB_TT_SL; TH1F * h_J1M_antitagSR_TT_SL; TH1F * h_J1M_antitagSB_TT_SL;
  TH1F * h_J2M_doubletagSR_TT_SL; TH1F * h_J2M_doubletagSB_TT_SL; TH1F * h_J2M_tagSR_TT_SL; TH1F * h_J2M_tagSB_TT_SL; TH1F * h_J2M_antitagSR_TT_SL; TH1F * h_J2M_antitagSB_TT_SL;
  TH1F * h_J1M_mjBins_doubletagSR_TT_SL; TH1F * h_J1M_mjBins_doubletagSB_TT_SL; TH1F * h_J1M_mjBins_tagSR_TT_SL; TH1F * h_J1M_mjBins_tagSB_TT_SL; TH1F * h_J1M_mjBins_antitagSR_TT_SL; TH1F * h_J1M_mjBins_antitagSB_TT_SL;
  TH1F * h_J2M_mjBins_doubletagSR_TT_SL; TH1F * h_J2M_mjBins_doubletagSB_TT_SL; TH1F * h_J2M_mjBins_tagSR_TT_SL; TH1F * h_J2M_mjBins_tagSB_TT_SL; TH1F * h_J2M_mjBins_antitagSR_TT_SL; TH1F * h_J2M_mjBins_antitagSB_TT_SL;

  TH1F * h_baseline_j1mass_SnglT; TH1F * h_J1M_doubletagSR_SnglT; TH1F * h_J1M_doubletagSB_SnglT; TH1F * h_J1M_tagSR_SnglT; TH1F * h_J1M_tagSB_SnglT; TH1F * h_J1M_antitagSR_SnglT; TH1F * h_J1M_antitagSB_SnglT;
  TH1F * h_baseline_j2mass_SnglT; TH1F * h_J2M_doubletagSR_SnglT; TH1F * h_J2M_doubletagSB_SnglT; TH1F * h_J2M_tagSR_SnglT; TH1F * h_J2M_tagSB_SnglT; TH1F * h_J2M_antitagSR_SnglT; TH1F * h_J2M_antitagSB_SnglT;
  TH1F * h_J1M_mjBins_doubletagSR_SnglT; TH1F * h_J1M_mjBins_doubletagSB_SnglT; TH1F * h_J1M_mjBins_tagSR_SnglT; TH1F * h_J1M_mjBins_tagSB_SnglT; TH1F * h_J1M_mjBins_antitagSR_SnglT; TH1F * h_J1M_mjBins_antitagSB_SnglT;
  TH1F * h_J2M_mjBins_doubletagSR_SnglT; TH1F * h_J2M_mjBins_doubletagSB_SnglT; TH1F * h_J2M_mjBins_tagSR_SnglT; TH1F * h_J2M_mjBins_tagSB_SnglT; TH1F * h_J2M_mjBins_antitagSR_SnglT; TH1F * h_J2M_mjBins_antitagSB_SnglT;

  TH1F * h_baseline_j1mass_QCD; TH1F * h_J1M_doubletagSR_QCD; TH1F * h_J1M_doubletagSB_QCD; TH1F * h_J1M_tagSR_QCD; TH1F * h_J1M_tagSB_QCD; TH1F * h_J1M_antitagSR_QCD; TH1F * h_J1M_antitagSB_QCD;
  TH1F * h_baseline_j2mass_QCD; TH1F * h_J2M_doubletagSR_QCD; TH1F * h_J2M_doubletagSB_QCD; TH1F * h_J2M_tagSR_QCD; TH1F * h_J2M_tagSB_QCD; TH1F * h_J2M_antitagSR_QCD; TH1F * h_J2M_antitagSB_QCD;
  TH1F * h_J1M_mjBins_doubletagSR_QCD; TH1F * h_J1M_mjBins_doubletagSB_QCD; TH1F * h_J1M_mjBins_tagSR_QCD; TH1F * h_J1M_mjBins_tagSB_QCD; TH1F * h_J1M_mjBins_antitagSR_QCD; TH1F * h_J1M_mjBins_antitagSB_QCD;
  TH1F * h_J2M_mjBins_doubletagSR_QCD; TH1F * h_J2M_mjBins_doubletagSB_QCD; TH1F * h_J2M_mjBins_tagSR_QCD; TH1F * h_J2M_mjBins_tagSB_QCD; TH1F * h_J2M_mjBins_antitagSR_QCD; TH1F * h_J2M_mjBins_antitagSB_QCD;

  TH1F * h_baseline_j1mass_WJets; TH1F * h_J1M_doubletagSR_WJets; TH1F * h_J1M_doubletagSB_WJets; TH1F * h_J1M_tagSR_WJets; TH1F * h_J1M_tagSB_WJets; TH1F * h_J1M_antitagSR_WJets; TH1F * h_J1M_antitagSB_WJets;
  TH1F * h_baseline_j2mass_WJets; TH1F * h_J2M_doubletagSR_WJets; TH1F * h_J2M_doubletagSB_WJets; TH1F * h_J2M_tagSR_WJets; TH1F * h_J2M_tagSB_WJets; TH1F * h_J2M_antitagSR_WJets; TH1F * h_J2M_antitagSB_WJets;
  TH1F * h_J1M_mjBins_doubletagSR_WJets; TH1F * h_J1M_mjBins_doubletagSB_WJets; TH1F * h_J1M_mjBins_tagSR_WJets; TH1F * h_J1M_mjBins_tagSB_WJets; TH1F * h_J1M_mjBins_antitagSR_WJets; TH1F * h_J1M_mjBins_antitagSB_WJets;
  TH1F * h_J2M_mjBins_doubletagSR_WJets; TH1F * h_J2M_mjBins_doubletagSB_WJets; TH1F * h_J2M_mjBins_tagSR_WJets; TH1F * h_J2M_mjBins_tagSB_WJets; TH1F * h_J2M_mjBins_antitagSR_WJets; TH1F * h_J2M_mjBins_antitagSB_WJets;

  TH1F * h_baseline_j1mass_ZJets; TH1F * h_J1M_doubletagSR_ZJets; TH1F * h_J1M_doubletagSB_ZJets; TH1F * h_J1M_tagSR_ZJets; TH1F * h_J1M_tagSB_ZJets; TH1F * h_J1M_antitagSR_ZJets; TH1F * h_J1M_antitagSB_ZJets;
  TH1F * h_baseline_j2mass_ZJets; TH1F * h_J2M_doubletagSR_ZJets; TH1F * h_J2M_doubletagSB_ZJets; TH1F * h_J2M_tagSR_ZJets; TH1F * h_J2M_tagSB_ZJets; TH1F * h_J2M_antitagSR_ZJets; TH1F * h_J2M_antitagSB_ZJets;
  TH1F * h_J1M_mjBins_doubletagSR_ZJets; TH1F * h_J1M_mjBins_doubletagSB_ZJets; TH1F * h_J1M_mjBins_tagSR_ZJets; TH1F * h_J1M_mjBins_tagSB_ZJets; TH1F * h_J1M_mjBins_antitagSR_ZJets; TH1F * h_J1M_mjBins_antitagSB_ZJets;
  TH1F * h_J2M_mjBins_doubletagSR_ZJets; TH1F * h_J2M_mjBins_doubletagSB_ZJets; TH1F * h_J2M_mjBins_tagSR_ZJets; TH1F * h_J2M_mjBins_tagSB_ZJets; TH1F * h_J2M_mjBins_antitagSR_ZJets; TH1F * h_J2M_mjBins_antitagSB_ZJets;

  TH1F * h_J1M_doubletagSR_GJets; TH1F * h_J1M_doubletagSB_GJets; TH1F * h_J1M_tagSR_GJets; TH1F * h_J1M_tagSB_GJets; TH1F * h_J1M_antitagSR_GJets; TH1F * h_J1M_antitagSB_GJets;
  TH1F * h_J2M_doubletagSR_GJets; TH1F * h_J2M_doubletagSB_GJets; TH1F * h_J2M_tagSR_GJets; TH1F * h_J2M_tagSB_GJets; TH1F * h_J2M_antitagSR_GJets; TH1F * h_J2M_antitagSB_GJets;
  TH1F * h_J1M_mjBins_doubletagSR_GJets; TH1F * h_J1M_mjBins_doubletagSB_GJets; TH1F * h_J1M_mjBins_tagSR_GJets; TH1F * h_J1M_mjBins_tagSB_GJets; TH1F * h_J1M_mjBins_antitagSR_GJets; TH1F * h_J1M_mjBins_antitagSB_GJets;
  TH1F * h_J2M_mjBins_doubletagSR_GJets; TH1F * h_J2M_mjBins_doubletagSB_GJets; TH1F * h_J2M_mjBins_tagSR_GJets; TH1F * h_J2M_mjBins_tagSB_GJets; TH1F * h_J2M_mjBins_antitagSR_GJets; TH1F * h_J2M_mjBins_antitagSB_GJets;


  TH1F * h_baseline_j1mass_T5HH1600; TH1F * h_baseline_j2mass_T5HH1600; TH1F * h_J1M_doubletagSR_T5HH1600; TH1F * h_J2M_doubletagSR_T5HH1600; TH1F * h_J1M_doubletagSB_T5HH1600; TH1F * h_J2M_doubletagSB_T5HH1600;
  TH1F * h_J1M_antitagSR_T5HH1600; TH1F * h_J2M_antitagSR_T5HH1600; TH1F * h_J1M_antitagSB_T5HH1600; TH1F * h_J2M_antitagSB_T5HH1600;
  TH1F * h_J1M_tagSR_T5HH1600; TH1F * h_J2M_tagSR_T5HH1600; TH1F * h_J1M_tagSB_T5HH1600; TH1F * h_J2M_tagSB_T5HH1600;
  TH1F * h_baseline_j1mass_T5HH2000; TH1F * h_baseline_j2mass_T5HH2000; TH1F * h_J1M_doubletagSR_T5HH2000; TH1F * h_J2M_doubletagSR_T5HH2000; TH1F * h_J1M_doubletagSB_T5HH2000; TH1F * h_J2M_doubletagSB_T5HH2000;
  TH1F * h_J1M_antitagSR_T5HH2000; TH1F * h_J2M_antitagSR_T5HH2000; TH1F * h_J1M_antitagSB_T5HH2000; TH1F * h_J2M_antitagSB_T5HH2000;
  TH1F * h_J1M_tagSR_T5HH2000; TH1F * h_J2M_tagSR_T5HH2000; TH1F * h_J1M_tagSB_T5HH2000; TH1F * h_J2M_tagSB_T5HH2000;
  TH1F * h_baseline_j1mass_T5HH2200; TH1F * h_baseline_j2mass_T5HH2200; TH1F * h_J1M_doubletagSR_T5HH2200; TH1F * h_J2M_doubletagSR_T5HH2200; TH1F * h_J1M_doubletagSB_T5HH2200; TH1F * h_J2M_doubletagSB_T5HH2200;
  TH1F * h_J1M_antitagSR_T5HH2200; TH1F * h_J2M_antitagSR_T5HH2200; TH1F * h_J1M_antitagSB_T5HH2200; TH1F * h_J2M_antitagSB_T5HH2200;
  TH1F * h_J1M_tagSR_T5HH2200; TH1F * h_J2M_tagSR_T5HH2200; TH1F * h_J1M_tagSB_T5HH2200; TH1F * h_J2M_tagSB_T5HH2200;

  TH1F * h_baseline_j1mass_TChiHH600; TH1F * h_baseline_j2mass_TChiHH600; TH1F * h_J1M_doubletagSR_TChiHH600; TH1F * h_J2M_doubletagSR_TChiHH600; TH1F * h_J1M_doubletagSB_TChiHH600; TH1F * h_J2M_doubletagSB_TChiHH600;
  TH1F * h_J1M_antitagSR_TChiHH600; TH1F * h_J2M_antitagSR_TChiHH600; TH1F * h_J1M_antitagSB_TChiHH600; TH1F * h_J2M_antitagSB_TChiHH600;
  TH1F * h_J1M_tagSR_TChiHH600; TH1F * h_J2M_tagSR_TChiHH600; TH1F * h_J1M_tagSB_TChiHH600; TH1F * h_J2M_tagSB_TChiHH600;
  TH1F * h_baseline_j1mass_TChiHH800; TH1F * h_baseline_j2mass_TChiHH800; TH1F * h_J1M_doubletagSR_TChiHH800; TH1F * h_J2M_doubletagSR_TChiHH800; TH1F * h_J1M_doubletagSB_TChiHH800; TH1F * h_J2M_doubletagSB_TChiHH800;
  TH1F * h_J1M_antitagSR_TChiHH800; TH1F * h_J2M_antitagSR_TChiHH800; TH1F * h_J1M_antitagSB_TChiHH800; TH1F * h_J2M_antitagSB_TChiHH800;
  TH1F * h_J1M_tagSR_TChiHH800; TH1F * h_J2M_tagSR_TChiHH800; TH1F * h_J1M_tagSB_TChiHH800; TH1F * h_J2M_tagSB_TChiHH800;
  TH1F * h_baseline_j1mass_TChiHH1000; TH1F * h_baseline_j2mass_TChiHH1000; TH1F * h_J1M_doubletagSR_TChiHH1000; TH1F * h_J2M_doubletagSR_TChiHH1000; TH1F * h_J1M_doubletagSB_TChiHH1000; TH1F * h_J2M_doubletagSB_TChiHH1000;
  TH1F * h_J1M_antitagSR_TChiHH1000; TH1F * h_J2M_antitagSR_TChiHH1000; TH1F * h_J1M_antitagSB_TChiHH1000; TH1F * h_J2M_antitagSB_TChiHH1000;
  TH1F * h_J1M_tagSR_TChiHH1000; TH1F * h_J2M_tagSR_TChiHH1000; TH1F * h_J1M_tagSB_TChiHH1000; TH1F * h_J2M_tagSB_TChiHH1000;

  //For pt
  TH1F * h_J1Pt_doubletagSR_sum; TH1F * h_J2Pt_doubletagSR_sum; TH1F * h_J1Pt_doubletagSB_sum; TH1F * h_J2Pt_doubletagSB_sum;
  TH1F * h_J1Pt_antitagSR_sum; TH1F * h_J2Pt_antitagSR_sum; TH1F * h_J1Pt_antitagSB_sum; TH1F * h_J2Pt_antitagSB_sum;
  TH1F * h_J1Pt_tagSR_sum; TH1F * h_J2Pt_tagSR_sum; TH1F * h_J1Pt_tagSB_sum; TH1F * h_J2Pt_tagSB_sum;

  TH1F * h_baseline_j1pt_QCD; TH1F * h_baseline_j2pt_QCD; TH1F * h_J1Pt_doubletagSR_QCD; TH1F * h_J2Pt_doubletagSR_QCD; TH1F * h_J1Pt_doubletagSB_QCD; TH1F * h_J2Pt_doubletagSB_QCD;
  TH1F * h_J1Pt_antitagSR_QCD; TH1F * h_J2Pt_antitagSR_QCD; TH1F * h_J1Pt_antitagSB_QCD; TH1F * h_J2Pt_antitagSB_QCD;
  TH1F * h_J1Pt_tagSR_QCD; TH1F * h_J2Pt_tagSR_QCD; TH1F * h_J1Pt_tagSB_QCD; TH1F * h_J2Pt_tagSB_QCD;

  TH1F * h_baseline_j1pt_SnglT;  TH1F * h_baseline_j2pt_SnglT; TH1F * h_J1Pt_doubletagSR_SnglT; TH1F * h_J2Pt_doubletagSR_SnglT; TH1F * h_J1Pt_doubletagSB_SnglT; TH1F * h_J2Pt_doubletagSB_SnglT;
  TH1F * h_J1Pt_antitagSR_SnglT; TH1F * h_J2Pt_antitagSR_SnglT; TH1F * h_J1Pt_antitagSB_SnglT; TH1F * h_J2Pt_antitagSB_SnglT;
  TH1F * h_J1Pt_tagSR_SnglT; TH1F * h_J2Pt_tagSR_SnglT; TH1F * h_J1Pt_tagSB_SnglT; TH1F * h_J2Pt_tagSB_SnglT;

  TH1F * h_J1Pt_doubletagSR_GJets; TH1F * h_J2Pt_doubletagSR_GJets; TH1F * h_J1Pt_doubletagSB_GJets; TH1F * h_J2Pt_doubletagSB_GJets;
  TH1F * h_J1Pt_antitagSR_GJets; TH1F * h_J2Pt_antitagSR_GJets; TH1F * h_J1Pt_antitagSB_GJets; TH1F * h_J2Pt_antitagSB_GJets;
  TH1F * h_J1Pt_tagSR_GJets; TH1F * h_J2Pt_tagSR_GJets; TH1F * h_J1Pt_tagSB_GJets; TH1F * h_J2Pt_tagSB_GJets;

  TH1F * h_baseline_j1pt_TT;  TH1F * h_baseline_j2pt_TT;  TH1F * h_J1Pt_doubletagSR_TT; TH1F * h_J2Pt_doubletagSR_TT; TH1F * h_J1Pt_doubletagSB_TT; TH1F * h_J2Pt_doubletagSB_TT;
  TH1F * h_J1Pt_antitagSR_TT; TH1F * h_J2Pt_antitagSR_TT; TH1F * h_J1Pt_antitagSB_TT; TH1F * h_J2Pt_antitagSB_TT;
  TH1F * h_J1Pt_tagSR_TT; TH1F * h_J2Pt_tagSR_TT; TH1F * h_J1Pt_tagSB_TT; TH1F * h_J2Pt_tagSB_TT;

  TH1F * h_J1Pt_doubletagSR_TT_Di; TH1F * h_J2Pt_doubletagSR_TT_Di; TH1F * h_J1Pt_doubletagSB_TT_Di; TH1F * h_J2Pt_doubletagSB_TT_Di;
  TH1F * h_J1Pt_antitagSR_TT_Di; TH1F * h_J2Pt_antitagSR_TT_Di; TH1F * h_J1Pt_antitagSB_TT_Di; TH1F * h_J2Pt_antitagSB_TT_Di;
  TH1F * h_J1Pt_tagSR_TT_Di; TH1F * h_J2Pt_tagSR_TT_Di; TH1F * h_J1Pt_tagSB_TT_Di; TH1F * h_J2Pt_tagSB_TT_Di;

  TH1F * h_J1Pt_doubletagSR_TT_SL; TH1F * h_J2Pt_doubletagSR_TT_SL; TH1F * h_J1Pt_doubletagSB_TT_SL; TH1F * h_J2Pt_doubletagSB_TT_SL;
  TH1F * h_J1Pt_antitagSR_TT_SL; TH1F * h_J2Pt_antitagSR_TT_SL; TH1F * h_J1Pt_antitagSB_TT_SL; TH1F * h_J2Pt_antitagSB_TT_SL;
  TH1F * h_J1Pt_tagSR_TT_SL; TH1F * h_J2Pt_tagSR_TT_SL; TH1F * h_J1Pt_tagSB_TT_SL; TH1F * h_J2Pt_tagSB_TT_SL;

  TH1F * h_baseline_j1pt_WJets;  TH1F * h_baseline_j2pt_WJets; TH1F * h_J1Pt_doubletagSR_WJets; TH1F * h_J2Pt_doubletagSR_WJets; TH1F * h_J1Pt_doubletagSB_WJets; TH1F * h_J2Pt_doubletagSB_WJets;
  TH1F * h_J1Pt_antitagSR_WJets; TH1F * h_J2Pt_antitagSR_WJets; TH1F * h_J1Pt_antitagSB_WJets; TH1F * h_J2Pt_antitagSB_WJets;
  TH1F * h_J1Pt_tagSR_WJets; TH1F * h_J2Pt_tagSR_WJets; TH1F * h_J1Pt_tagSB_WJets; TH1F * h_J2Pt_tagSB_WJets;

  TH1F * h_baseline_j1pt_ZJets; TH1F * h_baseline_j2pt_ZJets; TH1F * h_J1Pt_doubletagSR_ZJets; TH1F * h_J2Pt_doubletagSR_ZJets; TH1F * h_J1Pt_doubletagSB_ZJets; TH1F * h_J2Pt_doubletagSB_ZJets;
  TH1F * h_J1Pt_antitagSR_ZJets; TH1F * h_J2Pt_antitagSR_ZJets; TH1F * h_J1Pt_antitagSB_ZJets; TH1F * h_J2Pt_antitagSB_ZJets;
  TH1F * h_J1Pt_tagSR_ZJets; TH1F * h_J2Pt_tagSR_ZJets; TH1F * h_J1Pt_tagSB_ZJets; TH1F * h_J2Pt_tagSB_ZJets;

  TH1F * h_baseline_j1pt_T5HH1600; TH1F * h_baseline_j2pt_T5HH1600; TH1F * h_J1Pt_doubletagSR_T5HH1600; TH1F * h_J2Pt_doubletagSR_T5HH1600; TH1F * h_J1Pt_doubletagSB_T5HH1600; TH1F * h_J2Pt_doubletagSB_T5HH1600;
  TH1F * h_J1Pt_antitagSR_T5HH1600; TH1F * h_J2Pt_antitagSR_T5HH1600; TH1F * h_J1Pt_antitagSB_T5HH1600; TH1F * h_J2Pt_antitagSB_T5HH1600;
  TH1F * h_J1Pt_tagSR_T5HH1600; TH1F * h_J2Pt_tagSR_T5HH1600; TH1F * h_J1Pt_tagSB_T5HH1600; TH1F * h_J2Pt_tagSB_T5HH1600;

  TH1F * h_baseline_j1pt_T5HH2000; TH1F * h_baseline_j2pt_T5HH2000; TH1F * h_J1Pt_doubletagSR_T5HH2000; TH1F * h_J2Pt_doubletagSR_T5HH2000; TH1F * h_J1Pt_doubletagSB_T5HH2000; TH1F * h_J2Pt_doubletagSB_T5HH2000;
  TH1F * h_J1Pt_antitagSR_T5HH2000; TH1F * h_J2Pt_antitagSR_T5HH2000; TH1F * h_J1Pt_antitagSB_T5HH2000; TH1F * h_J2Pt_antitagSB_T5HH2000;
  TH1F * h_J1Pt_tagSR_T5HH2000; TH1F * h_J2Pt_tagSR_T5HH2000; TH1F * h_J1Pt_tagSB_T5HH2000; TH1F * h_J2Pt_tagSB_T5HH2000;

  TH1F * h_baseline_j1pt_T5HH2200; TH1F * h_baseline_j2pt_T5HH2200; TH1F * h_J1Pt_doubletagSR_T5HH2200; TH1F * h_J2Pt_doubletagSR_T5HH2200; TH1F * h_J1Pt_doubletagSB_T5HH2200; TH1F * h_J2Pt_doubletagSB_T5HH2200;
  TH1F * h_J1Pt_antitagSR_T5HH2200; TH1F * h_J2Pt_antitagSR_T5HH2200; TH1F * h_J1Pt_antitagSB_T5HH2200; TH1F * h_J2Pt_antitagSB_T5HH2200;
  TH1F * h_J1Pt_tagSR_T5HH2200; TH1F * h_J2Pt_tagSR_T5HH2200; TH1F * h_J1Pt_tagSB_T5HH2200; TH1F * h_J2Pt_tagSB_T5HH2200;

  TH1F * h_baseline_j1pt_TChiHH600; TH1F * h_baseline_j2pt_TChiHH600; TH1F * h_J1Pt_doubletagSR_TChiHH600; TH1F * h_J2Pt_doubletagSR_TChiHH600; TH1F * h_J1Pt_doubletagSB_TChiHH600; TH1F * h_J2Pt_doubletagSB_TChiHH600;
  TH1F * h_J1Pt_antitagSR_TChiHH600; TH1F * h_J2Pt_antitagSR_TChiHH600; TH1F * h_J1Pt_antitagSB_TChiHH600; TH1F * h_J2Pt_antitagSB_TChiHH600;
  TH1F * h_J1Pt_tagSR_TChiHH600; TH1F * h_J2Pt_tagSR_TChiHH600; TH1F * h_J1Pt_tagSB_TChiHH600; TH1F * h_J2Pt_tagSB_TChiHH600;

  TH1F * h_baseline_j1pt_TChiHH800; TH1F * h_baseline_j2pt_TChiHH800; TH1F * h_J1Pt_doubletagSR_TChiHH800; TH1F * h_J2Pt_doubletagSR_TChiHH800; TH1F * h_J1Pt_doubletagSB_TChiHH800; TH1F * h_J2Pt_doubletagSB_TChiHH800;
  TH1F * h_J1Pt_antitagSR_TChiHH800; TH1F * h_J2Pt_antitagSR_TChiHH800; TH1F * h_J1Pt_antitagSB_TChiHH800; TH1F * h_J2Pt_antitagSB_TChiHH800;
  TH1F * h_J1Pt_tagSR_TChiHH800; TH1F * h_J2Pt_tagSR_TChiHH800; TH1F * h_J1Pt_tagSB_TChiHH800; TH1F * h_J2Pt_tagSB_TChiHH800;

  TH1F * h_baseline_j1pt_TChiHH1000; TH1F * h_baseline_j2pt_TChiHH1000; TH1F * h_J1Pt_doubletagSR_TChiHH1000; TH1F * h_J2Pt_doubletagSR_TChiHH1000; TH1F * h_J1Pt_doubletagSB_TChiHH1000; TH1F * h_J2Pt_doubletagSB_TChiHH1000;
  TH1F * h_J1Pt_antitagSR_TChiHH1000; TH1F * h_J2Pt_antitagSR_TChiHH1000; TH1F * h_J1Pt_antitagSB_TChiHH1000; TH1F * h_J2Pt_antitagSB_TChiHH1000;
  TH1F * h_J1Pt_tagSR_TChiHH1000; TH1F * h_J2Pt_tagSR_TChiHH1000; TH1F * h_J1Pt_tagSB_TChiHH1000; TH1F * h_J2Pt_tagSB_TChiHH1000;


  TH1F * h_DOpt1_sum; TH1F * h_COpt1_sum;
  TH1F * h_DOpt2_sum; TH1F * h_COpt2_sum;
  TH1F * h_DOpt3_sum; TH1F * h_COpt3_sum;
  TH1F * h_DOpt4_sum; TH1F * h_COpt4_sum;
  TH1F * h_DOpt5_sum; TH1F * h_COpt5_sum;

  TH1F * hP_DOpt1_sum; TH1F * hP_COpt1_sum;
  TH1F * hP_DOpt2_sum; TH1F * hP_COpt2_sum;
  TH1F * hP_DOpt3_sum; TH1F * hP_COpt3_sum;
  TH1F * hP_DOpt4_sum; TH1F * hP_COpt4_sum;
  TH1F * hP_DOpt5_sum; TH1F * hP_COpt5_sum;

  TH1F * h_D5Opt1_sum; TH1F * h_C5Opt1_sum;
  TH1F * h_D5Opt2_sum; TH1F * h_C5Opt2_sum;
  TH1F * h_D5Opt3_sum; TH1F * h_C5Opt3_sum;
  TH1F * h_D5Opt4_sum; TH1F * h_C5Opt4_sum;


  TH1F * h_COpt1_TT;    TH1F * h_DOpt1_TT;
  TH1F * h_COpt1_GJets; TH1F * h_DOpt1_GJets;
  TH1F * h_COpt1_WJets; TH1F * h_DOpt1_WJets;
  TH1F * h_COpt1_ZJets; TH1F * h_DOpt1_ZJets;
  TH1F * h_COpt1_QCD; TH1F * h_DOpt1_QCD;
  TH1F * h_COpt1_SnglT;TH1F * h_DOpt1_SnglT;

  TH1F * h_COpt2_TT;    TH1F * h_DOpt2_TT;
  TH1F * h_COpt2_GJets; TH1F * h_DOpt2_GJets;
  TH1F * h_COpt2_WJets; TH1F * h_DOpt2_WJets;
  TH1F * h_COpt2_ZJets; TH1F * h_DOpt2_ZJets;
  TH1F * h_COpt2_QCD; TH1F * h_DOpt2_QCD;
  TH1F * h_COpt2_SnglT;TH1F * h_DOpt2_SnglT;

  TH1F * h_COpt3_TT;    TH1F * h_DOpt3_TT;
  TH1F * h_COpt3_GJets; TH1F * h_DOpt3_GJets;
  TH1F * h_COpt3_WJets; TH1F * h_DOpt3_WJets;
  TH1F * h_COpt3_ZJets; TH1F * h_DOpt3_ZJets;
  TH1F * h_COpt3_QCD; TH1F * h_DOpt3_QCD;
  TH1F * h_COpt3_SnglT;TH1F * h_DOpt3_SnglT;

  TH1F * h_COpt4_TT;    TH1F * h_DOpt4_TT;
  TH1F * h_COpt4_GJets; TH1F * h_DOpt4_GJets;
  TH1F * h_COpt4_WJets; TH1F * h_DOpt4_WJets;
  TH1F * h_COpt4_ZJets; TH1F * h_DOpt4_ZJets;
  TH1F * h_COpt4_QCD; TH1F * h_DOpt4_QCD;
  TH1F * h_COpt4_SnglT;TH1F * h_DOpt4_SnglT;

  TH1F * h_COpt5_TT;    TH1F * h_DOpt5_TT;
  TH1F * h_COpt5_GJets; TH1F * h_DOpt5_GJets;
  TH1F * h_COpt5_WJets; TH1F * h_DOpt5_WJets;
  TH1F * h_COpt5_ZJets; TH1F * h_DOpt5_ZJets;
  TH1F * h_COpt5_QCD; TH1F * h_DOpt5_QCD;
  TH1F * h_COpt5_SnglT;TH1F * h_DOpt5_SnglT;

  if (whichRegion=="signal") {
    // h_A_sum = (TH1F*)f->Get("MET_doubletagSR_sum"); h_B_sum = (TH1F*)f->Get("MET_doubletagSB_sum");
    // h_A1_sum = (TH1F*)f->Get("MET_tagSR_sum"); h_B1_sum = (TH1F*)f->Get("MET_tagSB_sum");
    // h_C_sum = (TH1F*)f->Get("MET_antitagSR_sum"); h_D_sum = (TH1F*)f->Get("MET_antitagSB_sum");

    h_AlowPU_sum = (TH1F*)f->Get("MET_doubletagSRlowPU_sum"); h_BlowPU_sum = (TH1F*)f->Get("MET_doubletagSBlowPU_sum");
    h_A1lowPU_sum = (TH1F*)f->Get("MET_tagSRlowPU_sum"); h_B1lowPU_sum = (TH1F*)f->Get("MET_tagSBlowPU_sum");
    h_ClowPU_sum = (TH1F*)f->Get("MET_antitagSRlowPU_sum"); h_DlowPU_sum = (TH1F*)f->Get("MET_antitagSBlowPU_sum");

    h_AhighPU_sum = (TH1F*)f->Get("MET_doubletagSRhighPU_sum"); h_BhighPU_sum = (TH1F*)f->Get("MET_doubletagSBhighPU_sum");
    h_A1highPU_sum = (TH1F*)f->Get("MET_tagSRhighPU_sum"); h_B1highPU_sum = (TH1F*)f->Get("MET_tagSBhighPU_sum");
    h_ChighPU_sum = (TH1F*)f->Get("MET_antitagSRhighPU_sum"); h_DhighPU_sum = (TH1F*)f->Get("MET_antitagSBhighPU_sum");

    if (runDM) {
      h_avgM_A_sum = (TH1F*)f->Get("avgM_doubletagSR_sum");
      h_avgM_B_sum = (TH1F*)f->Get("avgM_doubletagSB_sum");
      h_avgM_A1_sum = (TH1F*)f->Get("avgM_tagSR_sum");
      h_avgM_B1_sum = (TH1F*)f->Get("avgM_tagSB_sum");
      h_avgM_C_sum = (TH1F*)f->Get("avgM_antitagSR_sum");
      h_avgM_D_sum = (TH1F*)f->Get("avgM_antitagSB_sum");
    }

    h_A5_sum = (TH1F*)f->Get("MET5_doubletagSR_sum"); h_B5_sum = (TH1F*)f->Get("MET5_doubletagSB_sum");
    h_A15_sum = (TH1F*)f->Get("MET5_tagSR_sum"); h_B15_sum = (TH1F*)f->Get("MET5_tagSB_sum");
    h_C5_sum = (TH1F*)f->Get("MET5_antitagSR_sum"); h_D5_sum = (TH1F*)f->Get("MET5_antitagSB_sum");

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
    h_A_SnglT = (TH1F*)f->Get("MET_doubletagSR_SnglT"); h_B_SnglT = (TH1F*)f->Get("MET_doubletagSB_SnglT");
    h_A1_SnglT = (TH1F*)f->Get("MET_tagSR_SnglT"); h_B1_SnglT = (TH1F*)f->Get("MET_tagSB_SnglT");
    h_C_SnglT = (TH1F*)f->Get("MET_antitagSR_SnglT"); h_D_SnglT = (TH1F*)f->Get("MET_antitagSB_SnglT");

    h_COpt1_QCD = (TH1F*)f->Get("MET_antitagSROpt1_QCD"); h_DOpt1_QCD = (TH1F*)f->Get("MET_antitagSBOpt1_QCD");
    h_COpt1_SnglT = (TH1F*)f->Get("MET_antitagSROpt1_SnglT"); h_DOpt1_SnglT = (TH1F*)f->Get("MET_antitagSBOpt1_SnglT");
    h_COpt1_TT = (TH1F*)f->Get("MET_antitagSROpt1_TT"); h_DOpt1_TT = (TH1F*)f->Get("MET_antitagSBOpt1_TT");
    h_COpt1_WJets = (TH1F*)f->Get("MET_antitagSROpt1_WJets"); h_DOpt1_WJets = (TH1F*)f->Get("MET_antitagSBOpt1_WJets");
    h_COpt1_ZJets = (TH1F*)f->Get("MET_antitagSROpt1_ZJets"); h_DOpt1_ZJets = (TH1F*)f->Get("MET_antitagSBOpt1_ZJets");

    h_C5Opt1_sum = (TH1F*)f->Get("MET5_antitagSROpt1_sum"); h_D5Opt1_sum = (TH1F*)f->Get("MET5_antitagSBOpt1_sum");
    h_C5Opt2_sum = (TH1F*)f->Get("MET5_antitagSROpt2_sum"); h_D5Opt2_sum = (TH1F*)f->Get("MET5_antitagSBOpt2_sum");
    h_C5Opt3_sum = (TH1F*)f->Get("MET5_antitagSROpt3_sum"); h_D5Opt3_sum = (TH1F*)f->Get("MET5_antitagSBOpt3_sum");
    h_C5Opt4_sum = (TH1F*)f->Get("MET5_antitagSROpt4_sum"); h_D5Opt4_sum = (TH1F*)f->Get("MET5_antitagSBOpt4_sum");

    h_COpt2_QCD = (TH1F*)f->Get("MET_antitagSROpt2_QCD"); h_DOpt2_QCD = (TH1F*)f->Get("MET_antitagSBOpt2_QCD");
    h_COpt2_SnglT = (TH1F*)f->Get("MET_antitagSROpt2_SnglT"); h_DOpt2_SnglT = (TH1F*)f->Get("MET_antitagSBOpt2_SnglT");
    h_COpt2_TT = (TH1F*)f->Get("MET_antitagSROpt2_TT"); h_DOpt2_TT = (TH1F*)f->Get("MET_antitagSBOpt2_TT");
    h_COpt2_WJets = (TH1F*)f->Get("MET_antitagSROpt2_WJets"); h_DOpt2_WJets = (TH1F*)f->Get("MET_antitagSBOpt2_WJets");
    h_COpt2_ZJets = (TH1F*)f->Get("MET_antitagSROpt2_ZJets"); h_DOpt2_ZJets = (TH1F*)f->Get("MET_antitagSBOpt2_ZJets");

    h_COpt3_QCD = (TH1F*)f->Get("MET_antitagSROpt3_QCD"); h_DOpt3_QCD = (TH1F*)f->Get("MET_antitagSBOpt3_QCD");
    h_COpt3_SnglT = (TH1F*)f->Get("MET_antitagSROpt3_SnglT"); h_DOpt3_SnglT = (TH1F*)f->Get("MET_antitagSBOpt3_SnglT");
    h_COpt3_TT = (TH1F*)f->Get("MET_antitagSROpt3_TT"); h_DOpt3_TT = (TH1F*)f->Get("MET_antitagSBOpt3_TT");
    h_COpt3_WJets = (TH1F*)f->Get("MET_antitagSROpt3_WJets"); h_DOpt3_WJets = (TH1F*)f->Get("MET_antitagSBOpt3_WJets");
    h_COpt3_ZJets = (TH1F*)f->Get("MET_antitagSROpt3_ZJets"); h_DOpt3_ZJets = (TH1F*)f->Get("MET_antitagSBOpt3_ZJets");

    h_COpt4_QCD = (TH1F*)f->Get("MET_antitagSROpt4_QCD"); h_DOpt4_QCD = (TH1F*)f->Get("MET_antitagSBOpt4_QCD");
    h_COpt4_SnglT = (TH1F*)f->Get("MET_antitagSROpt4_SnglT"); h_DOpt4_SnglT = (TH1F*)f->Get("MET_antitagSBOpt4_SnglT");
    h_COpt4_TT = (TH1F*)f->Get("MET_antitagSROpt4_TT"); h_DOpt4_TT = (TH1F*)f->Get("MET_antitagSBOpt4_TT");
    h_COpt4_WJets = (TH1F*)f->Get("MET_antitagSROpt4_WJets"); h_DOpt4_WJets = (TH1F*)f->Get("MET_antitagSBOpt4_WJets");
    h_COpt4_ZJets = (TH1F*)f->Get("MET_antitagSROpt4_ZJets"); h_DOpt4_ZJets = (TH1F*)f->Get("MET_antitagSBOpt4_ZJets");

    h_COpt5_QCD = (TH1F*)f->Get("MET_antitagSROpt5_QCD"); h_DOpt5_QCD = (TH1F*)f->Get("MET_antitagSBOpt5_QCD");
    h_COpt5_SnglT = (TH1F*)f->Get("MET_antitagSROpt5_SnglT"); h_DOpt5_SnglT = (TH1F*)f->Get("MET_antitagSBOpt5_SnglT");
    h_COpt5_TT = (TH1F*)f->Get("MET_antitagSROpt5_TT"); h_DOpt5_TT = (TH1F*)f->Get("MET_antitagSBOpt5_TT");
    h_COpt5_WJets = (TH1F*)f->Get("MET_antitagSROpt5_WJets"); h_DOpt5_WJets = (TH1F*)f->Get("MET_antitagSBOpt5_WJets");
    h_COpt5_ZJets = (TH1F*)f->Get("MET_antitagSROpt5_ZJets"); h_DOpt5_ZJets = (TH1F*)f->Get("MET_antitagSBOpt5_ZJets");

    //Adding uncertainties to 0-event bins
    h_A_sum = make0EventUncSum({h_A_QCD,h_A_TT,h_A_WJets,h_A_ZJets,h_A_SnglT});
    h_B_sum = make0EventUncSum({h_B_QCD,h_B_TT,h_B_WJets,h_B_ZJets,h_B_SnglT});
    h_A1_sum = make0EventUncSum({h_A1_QCD,h_A1_TT,h_A1_WJets,h_A1_ZJets,h_A1_SnglT});
    h_B1_sum = make0EventUncSum({h_B1_QCD,h_B1_TT,h_B1_WJets,h_B1_ZJets,h_B1_SnglT});
    h_C_sum = make0EventUncSum({h_C_QCD,h_C_TT,h_C_WJets,h_C_ZJets,h_C_SnglT});
    h_D_sum = make0EventUncSum({h_D_QCD,h_D_TT,h_D_WJets,h_D_ZJets,h_D_SnglT});

    h_COpt1_sum = make0EventUncSum({h_COpt1_QCD,h_COpt1_TT,h_COpt1_WJets,h_COpt1_ZJets,h_COpt1_SnglT});
    h_DOpt1_sum = make0EventUncSum({h_DOpt1_QCD,h_DOpt1_TT,h_DOpt1_WJets,h_DOpt1_ZJets,h_DOpt1_SnglT});
    h_COpt2_sum = make0EventUncSum({h_COpt2_QCD,h_COpt2_TT,h_COpt2_WJets,h_COpt2_ZJets,h_COpt2_SnglT});
    h_DOpt2_sum = make0EventUncSum({h_DOpt2_QCD,h_DOpt2_TT,h_DOpt2_WJets,h_DOpt2_ZJets,h_DOpt2_SnglT});
    h_COpt3_sum = make0EventUncSum({h_COpt3_QCD,h_COpt3_TT,h_COpt3_WJets,h_COpt3_ZJets,h_COpt3_SnglT});
    h_DOpt3_sum = make0EventUncSum({h_DOpt3_QCD,h_DOpt3_TT,h_DOpt3_WJets,h_DOpt3_ZJets,h_DOpt3_SnglT});
    h_COpt4_sum = make0EventUncSum({h_COpt4_QCD,h_COpt4_TT,h_COpt4_WJets,h_COpt4_ZJets,h_COpt4_SnglT});
    h_DOpt4_sum = make0EventUncSum({h_DOpt4_QCD,h_DOpt4_TT,h_DOpt4_WJets,h_DOpt4_ZJets,h_DOpt4_SnglT});
    h_COpt5_sum = make0EventUncSum({h_COpt5_QCD,h_COpt5_TT,h_COpt5_WJets,h_COpt5_ZJets,h_COpt5_SnglT});
    h_DOpt5_sum = make0EventUncSum({h_DOpt5_QCD,h_DOpt5_TT,h_DOpt5_WJets,h_DOpt5_ZJets,h_DOpt5_SnglT});

    h_A_sum->Write("ABCDRegions/h_A_sum");
    h_A1_sum->Write("ABCDRegions/h_A1_sum");
    h_B_sum->Write("ABCDRegions/h_B_sum");
    h_B1_sum->Write("ABCDRegions/h_B1_sum");
    h_C_sum->Write("ABCDRegions/h_C_sum");
    h_D_sum->Write("ABCDRegions/h_D_sum");


    h_COpt1_sum->Write("ABCDRegions/h_C_M1_sum");
    h_COpt2_sum->Write("ABCDRegions/h_C_M0_sum");
    h_COpt3_sum->Write("ABCDRegions/h_C_T0_sum");
    h_COpt4_sum->Write("ABCDRegions/h_C_L1M0_sum");
    h_COpt5_sum->Write("ABCDRegions/h_C_M1T0_sum");
    h_DOpt1_sum->Write("ABCDRegions/h_D_M1_sum");
    h_DOpt2_sum->Write("ABCDRegions/h_D_M0_sum");
    h_DOpt3_sum->Write("ABCDRegions/h_D_T0_sum");
    h_DOpt4_sum->Write("ABCDRegions/h_D_L1M0_sum");
    h_DOpt5_sum->Write("ABCDRegions/h_D_M1T0_sum");


    if (0){
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
    }
    h_J1M_doubletagSR_sum = (TH1F*)f->Get("J1pt_M_doubletagSR_sum"); h_J2M_doubletagSR_sum = (TH1F*)f->Get("J2pt_M_doubletagSR_sum");
    h_J1M_doubletagSB_sum = (TH1F*)f->Get("J1pt_M_doubletagSB_sum"); h_J2M_doubletagSB_sum = (TH1F*)f->Get("J2pt_M_doubletagSB_sum");
    h_J1M_antitagSR_sum = (TH1F*)f->Get("J1pt_M_antitagSR_sum"); h_J2M_antitagSR_sum = (TH1F*)f->Get("J2pt_M_antitagSR_sum");
    h_J1M_antitagSB_sum = (TH1F*)f->Get("J1pt_M_antitagSB_sum"); h_J2M_antitagSB_sum = (TH1F*)f->Get("J2pt_M_antitagSB_sum");
    h_J1M_tagSR_sum = (TH1F*)f->Get("J1pt_M_tagSR_sum"); h_J2M_tagSR_sum = (TH1F*)f->Get("J2pt_M_tagSR_sum");
    h_J1M_tagSB_sum = (TH1F*)f->Get("J1pt_M_tagSB_sum"); h_J2M_tagSB_sum = (TH1F*)f->Get("J2pt_M_tagSB_sum");

    h_J1J2M_doubletagSR_sum = (TH1F*)h_J1M_doubletagSR_sum->Clone("J1J2M_doubletagSR_sum"); h_J1J2M_doubletagSR_sum->Add(h_J2M_doubletagSR_sum);
    h_J1J2M_doubletagSB_sum = (TH1F*)h_J1M_doubletagSB_sum->Clone("J1J2M_doubletagSB_sum"); h_J1J2M_doubletagSB_sum->Add(h_J2M_doubletagSB_sum);
    h_J1J2M_tagSR_sum = (TH1F*)h_J1M_tagSR_sum->Clone("J1J2M_tagSR_sum"); h_J1J2M_tagSR_sum->Add(h_J2M_tagSR_sum);
    h_J1J2M_tagSB_sum = (TH1F*)h_J1M_tagSB_sum->Clone("J1J2M_tagSB_sum"); h_J1J2M_tagSB_sum->Add(h_J2M_tagSB_sum);
    h_J1J2M_antitagSR_sum = (TH1F*)h_J1M_antitagSR_sum->Clone("J1J2M_antitagSR_sum"); h_J1J2M_antitagSR_sum->Add(h_J2M_antitagSR_sum);
    h_J1J2M_antitagSB_sum = (TH1F*)h_J1M_antitagSB_sum->Clone("J1J2M_antitagSB_sum"); h_J1J2M_antitagSB_sum->Add(h_J2M_antitagSB_sum);

    h_J1M_doubletagSR_ZJets = (TH1F*)f->Get("J1pt_M_doubletagSR_ZJets"); h_J2M_doubletagSR_ZJets = (TH1F*)f->Get("J2pt_M_doubletagSR_ZJets");
    h_J1M_doubletagSB_ZJets = (TH1F*)f->Get("J1pt_M_doubletagSB_ZJets"); h_J2M_doubletagSB_ZJets = (TH1F*)f->Get("J2pt_M_doubletagSB_ZJets");
    h_J1M_antitagSR_ZJets = (TH1F*)f->Get("J1pt_M_antitagSR_ZJets"); h_J2M_antitagSR_ZJets = (TH1F*)f->Get("J2pt_M_antitagSR_ZJets");
    h_J1M_antitagSB_ZJets = (TH1F*)f->Get("J1pt_M_antitagSB_ZJets"); h_J2M_antitagSB_ZJets = (TH1F*)f->Get("J2pt_M_antitagSB_ZJets");
    h_J1M_tagSR_ZJets = (TH1F*)f->Get("J1pt_M_tagSR_ZJets"); h_J2M_tagSR_ZJets = (TH1F*)f->Get("J2pt_M_tagSR_ZJets");
    h_J1M_tagSB_ZJets = (TH1F*)f->Get("J1pt_M_tagSB_ZJets"); h_J2M_tagSB_ZJets = (TH1F*)f->Get("J2pt_M_tagSB_ZJets");
    h_baseline_j1mass_ZJets = (TH1F*)f->Get("J1pt_M_baseline_ZJets"); h_baseline_j2mass_ZJets = (TH1F*)f->Get("J2pt_M_baseline_ZJets");

    h_J1M_doubletagSR_QCD = (TH1F*)f->Get("J1pt_M_doubletagSR_QCD"); h_J2M_doubletagSR_QCD = (TH1F*)f->Get("J2pt_M_doubletagSR_QCD");
    h_J1M_doubletagSB_QCD = (TH1F*)f->Get("J1pt_M_doubletagSB_QCD"); h_J2M_doubletagSB_QCD = (TH1F*)f->Get("J2pt_M_doubletagSB_QCD");
    h_J1M_antitagSR_QCD = (TH1F*)f->Get("J1pt_M_antitagSR_QCD"); h_J2M_antitagSR_QCD = (TH1F*)f->Get("J2pt_M_antitagSR_QCD");
    h_J1M_antitagSB_QCD = (TH1F*)f->Get("J1pt_M_antitagSB_QCD"); h_J2M_antitagSB_QCD = (TH1F*)f->Get("J2pt_M_antitagSB_QCD");
    h_J1M_tagSR_QCD = (TH1F*)f->Get("J1pt_M_tagSR_QCD"); h_J2M_tagSR_QCD = (TH1F*)f->Get("J2pt_M_tagSR_QCD");
    h_J1M_tagSB_QCD = (TH1F*)f->Get("J1pt_M_tagSB_QCD"); h_J2M_tagSB_QCD = (TH1F*)f->Get("J2pt_M_tagSB_QCD");
    h_baseline_j1mass_QCD = (TH1F*)f->Get("J1pt_M_baseline_QCD"); h_baseline_j2mass_QCD = (TH1F*)f->Get("J2pt_M_baseline_QCD");


    h_J1M_doubletagSR_SnglT = (TH1F*)f->Get("J1pt_M_doubletagSR_SnglT"); h_J2M_doubletagSR_SnglT = (TH1F*)f->Get("J2pt_M_doubletagSR_SnglT");
    h_J1M_doubletagSB_SnglT = (TH1F*)f->Get("J1pt_M_doubletagSB_SnglT"); h_J2M_doubletagSB_SnglT = (TH1F*)f->Get("J2pt_M_doubletagSB_SnglT");
    h_J1M_antitagSR_SnglT = (TH1F*)f->Get("J1pt_M_antitagSR_SnglT"); h_J2M_antitagSR_SnglT = (TH1F*)f->Get("J2pt_M_antitagSR_SnglT");
    h_J1M_antitagSB_SnglT = (TH1F*)f->Get("J1pt_M_antitagSB_SnglT"); h_J2M_antitagSB_SnglT = (TH1F*)f->Get("J2pt_M_antitagSB_SnglT");
    h_J1M_tagSR_SnglT = (TH1F*)f->Get("J1pt_M_tagSR_SnglT"); h_J2M_tagSR_SnglT = (TH1F*)f->Get("J2pt_M_tagSR_SnglT");
    h_J1M_tagSB_SnglT = (TH1F*)f->Get("J1pt_M_tagSB_SnglT"); h_J2M_tagSB_SnglT = (TH1F*)f->Get("J2pt_M_tagSB_SnglT");
    h_baseline_j1mass_SnglT = (TH1F*)f->Get("J1pt_M_baseline_SnglT"); h_baseline_j2mass_SnglT = (TH1F*)f->Get("J2pt_M_baseline_SnglT");


    h_J1M_doubletagSR_TT = (TH1F*)f->Get("J1pt_M_doubletagSR_TT"); h_J2M_doubletagSR_TT = (TH1F*)f->Get("J2pt_M_doubletagSR_TT");
    h_J1M_doubletagSB_TT = (TH1F*)f->Get("J1pt_M_doubletagSB_TT"); h_J2M_doubletagSB_TT = (TH1F*)f->Get("J2pt_M_doubletagSB_TT");
    h_J1M_antitagSR_TT = (TH1F*)f->Get("J1pt_M_antitagSR_TT"); h_J2M_antitagSR_TT = (TH1F*)f->Get("J2pt_M_antitagSR_TT");
    h_J1M_antitagSB_TT = (TH1F*)f->Get("J1pt_M_antitagSB_TT"); h_J2M_antitagSB_TT = (TH1F*)f->Get("J2pt_M_antitagSB_TT");
    h_J1M_tagSR_TT = (TH1F*)f->Get("J1pt_M_tagSR_TT"); h_J2M_tagSR_TT = (TH1F*)f->Get("J2pt_M_tagSR_TT");
    h_J1M_tagSB_TT = (TH1F*)f->Get("J1pt_M_tagSB_TT"); h_J2M_tagSB_TT = (TH1F*)f->Get("J2pt_M_tagSB_TT");
    h_baseline_j1mass_TT = (TH1F*)f->Get("J1pt_M_baseline_TT"); h_baseline_j2mass_TT = (TH1F*)f->Get("J2pt_M_baseline_TT");


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
    h_baseline_j1mass_WJets = (TH1F*)f->Get("J1pt_M_baseline_WJets"); h_baseline_j2mass_WJets = (TH1F*)f->Get("J2pt_M_baseline_WJets");


    h_J2M_mjBins_doubletagSR_sum = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_sum"); h_J2M_mjBins_doubletagSB_sum = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_sum");
    h_J2M_mjBins_antitagSR_sum = (TH1F*)f->Get("J2_M_jetBins_antitagSR_sum"); h_J2M_mjBins_antitagSB_sum = (TH1F*)f->Get("J2_M_jetBins_antitagSB_sum");
    h_J2M_mjBins_doubletagSR_QCD = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_QCD"); h_J2M_mjBins_doubletagSB_QCD = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_QCD");
    h_J2M_mjBins_antitagSR_QCD = (TH1F*)f->Get("J2_M_jetBins_antitagSR_QCD"); h_J2M_mjBins_antitagSB_QCD = (TH1F*)f->Get("J2_M_jetBins_antitagSB_QCD");
    h_J2M_mjBins_doubletagSR_SnglT = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_SnglT"); h_J2M_mjBins_doubletagSB_SnglT = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_SnglT");
    h_J2M_mjBins_antitagSR_SnglT = (TH1F*)f->Get("J2_M_jetBins_antitagSR_SnglT"); h_J2M_mjBins_antitagSB_SnglT = (TH1F*)f->Get("J2_M_jetBins_antitagSB_SnglT");

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
    if (0){
      h_J1M_doubletagSR_T5HH1600 = (TH1F*)fSignal->Get("J1pt_M_doubletagSR_T5qqqqZH1600"); h_J2M_doubletagSR_T5HH1600 = (TH1F*)fSignal->Get("J2pt_M_doubletagSR_T5qqqqZH1600");
      h_J1M_doubletagSB_T5HH1600 = (TH1F*)fSignal->Get("J1pt_M_doubletagSB_T5qqqqZH1600"); h_J2M_doubletagSB_T5HH1600 = (TH1F*)fSignal->Get("J2pt_M_doubletagSB_T5qqqqZH1600");
      h_J1M_antitagSR_T5HH1600 = (TH1F*)fSignal->Get("J1pt_M_antitagSR_T5qqqqZH1600"); h_J2M_antitagSR_T5HH1600 = (TH1F*)fSignal->Get("J2pt_M_antitagSR_T5qqqqZH1600");
      h_J1M_antitagSB_T5HH1600 = (TH1F*)fSignal->Get("J1pt_M_antitagSB_T5qqqqZH1600"); h_J2M_antitagSB_T5HH1600 = (TH1F*)fSignal->Get("J2pt_M_antitagSB_T5qqqqZH1600");
      h_J1M_tagSR_T5HH1600 = (TH1F*)fSignal->Get("J1pt_M_tagSR_T5qqqqZH1600"); h_J2M_tagSR_T5HH1600 = (TH1F*)fSignal->Get("J2pt_M_tagSR_T5qqqqZH1600");
      h_J1M_tagSB_T5HH1600 = (TH1F*)fSignal->Get("J1pt_M_tagSB_T5qqqqZH1600"); h_J2M_tagSB_T5HH1600 = (TH1F*)fSignal->Get("J2pt_M_tagSB_T5qqqqZH1600");
      h_baseline_j1mass_T5HH1600 = (TH1F*)fSignal->Get("J1pt_M_baseline_T5qqqqZH1600"); h_baseline_j2mass_T5HH1600 = (TH1F*)fSignal->Get("J2pt_M_baseline_T5qqqqZH1600");

      h_J1M_doubletagSR_T5HH2000 = (TH1F*)fSignal->Get("J1pt_M_doubletagSR_T5qqqqZH2000"); h_J2M_doubletagSR_T5HH2000 = (TH1F*)fSignal->Get("J2pt_M_doubletagSR_T5qqqqZH2000");
      h_J1M_doubletagSB_T5HH2000 = (TH1F*)fSignal->Get("J1pt_M_doubletagSB_T5qqqqZH2000"); h_J2M_doubletagSB_T5HH2000 = (TH1F*)fSignal->Get("J2pt_M_doubletagSB_T5qqqqZH2000");
      h_J1M_antitagSR_T5HH2000 = (TH1F*)fSignal->Get("J1pt_M_antitagSR_T5qqqqZH2000"); h_J2M_antitagSR_T5HH2000 = (TH1F*)fSignal->Get("J2pt_M_antitagSR_T5qqqqZH2000");
      h_J1M_antitagSB_T5HH2000 = (TH1F*)fSignal->Get("J1pt_M_antitagSB_T5qqqqZH2000"); h_J2M_antitagSB_T5HH2000 = (TH1F*)fSignal->Get("J2pt_M_antitagSB_T5qqqqZH2000");
      h_J1M_tagSR_T5HH2000 = (TH1F*)fSignal->Get("J1pt_M_tagSR_T5qqqqZH2000"); h_J2M_tagSR_T5HH2000 = (TH1F*)fSignal->Get("J2pt_M_tagSR_T5qqqqZH2000");
      h_J1M_tagSB_T5HH2000 = (TH1F*)fSignal->Get("J1pt_M_tagSB_T5qqqqZH2000"); h_J2M_tagSB_T5HH2000 = (TH1F*)fSignal->Get("J2pt_M_tagSB_T5qqqqZH2000");
      h_baseline_j1mass_T5HH2000 = (TH1F*)fSignal->Get("J1pt_M_baseline_T5qqqqZH2000"); h_baseline_j2mass_T5HH2000 = (TH1F*)fSignal->Get("J2pt_M_baseline_T5qqqqZH2000");

      h_J1M_doubletagSR_T5HH2200 = (TH1F*)fSignal->Get("J1pt_M_doubletagSR_T5qqqqZH2200"); h_J2M_doubletagSR_T5HH2200 = (TH1F*)fSignal->Get("J2pt_M_doubletagSR_T5qqqqZH2200");
      h_J1M_doubletagSB_T5HH2200 = (TH1F*)fSignal->Get("J1pt_M_doubletagSB_T5qqqqZH2200"); h_J2M_doubletagSB_T5HH2200 = (TH1F*)fSignal->Get("J2pt_M_doubletagSB_T5qqqqZH2200");
      h_J1M_antitagSR_T5HH2200 = (TH1F*)fSignal->Get("J1pt_M_antitagSR_T5qqqqZH2200"); h_J2M_antitagSR_T5HH2200 = (TH1F*)fSignal->Get("J2pt_M_antitagSR_T5qqqqZH2200");
      h_J1M_antitagSB_T5HH2200 = (TH1F*)fSignal->Get("J1pt_M_antitagSB_T5qqqqZH2200"); h_J2M_antitagSB_T5HH2200 = (TH1F*)fSignal->Get("J2pt_M_antitagSB_T5qqqqZH2200");
      h_J1M_tagSR_T5HH2200 = (TH1F*)fSignal->Get("J1pt_M_tagSR_T5qqqqZH2200"); h_J2M_tagSR_T5HH2200 = (TH1F*)fSignal->Get("J2pt_M_tagSR_T5qqqqZH2200");
      h_J1M_tagSB_T5HH2200 = (TH1F*)fSignal->Get("J1pt_M_tagSB_T5qqqqZH2200"); h_J2M_tagSB_T5HH2200 = (TH1F*)fSignal->Get("J2pt_M_tagSB_T5qqqqZH2200");
      h_baseline_j1mass_T5HH2200 = (TH1F*)fSignal->Get("J1pt_M_baseline_T5qqqqZH2200"); h_baseline_j2mass_T5HH2200 = (TH1F*)fSignal->Get("J2pt_M_baseline_T5qqqqZH2200");


      h_J1M_doubletagSR_TChiHH600 = (TH1F*)fSignal->Get("J1pt_M_doubletagSR_TChiHH600"); h_J2M_doubletagSR_TChiHH600 = (TH1F*)fSignal->Get("J2pt_M_doubletagSR_TChiHH600");
      h_J1M_doubletagSB_TChiHH600 = (TH1F*)fSignal->Get("J1pt_M_doubletagSB_TChiHH600"); h_J2M_doubletagSB_TChiHH600 = (TH1F*)fSignal->Get("J2pt_M_doubletagSB_TChiHH600");
      h_J1M_antitagSR_TChiHH600 = (TH1F*)fSignal->Get("J1pt_M_antitagSR_TChiHH600"); h_J2M_antitagSR_TChiHH600 = (TH1F*)fSignal->Get("J2pt_M_antitagSR_TChiHH600");
      h_J1M_antitagSB_TChiHH600 = (TH1F*)fSignal->Get("J1pt_M_antitagSB_TChiHH600"); h_J2M_antitagSB_TChiHH600 = (TH1F*)fSignal->Get("J2pt_M_antitagSB_TChiHH600");
      h_J1M_tagSR_TChiHH600 = (TH1F*)fSignal->Get("J1pt_M_tagSR_TChiHH600"); h_J2M_tagSR_TChiHH600 = (TH1F*)fSignal->Get("J2pt_M_tagSR_TChiHH600");
      h_J1M_tagSB_TChiHH600 = (TH1F*)fSignal->Get("J1pt_M_tagSB_TChiHH600"); h_J2M_tagSB_TChiHH600 = (TH1F*)fSignal->Get("J2pt_M_tagSB_TChiHH600");
      h_baseline_j1mass_TChiHH600 = (TH1F*)fSignal->Get("J1pt_M_baseline_TChiHH600"); h_baseline_j2mass_TChiHH600 = (TH1F*)fSignal->Get("J2pt_M_baseline_TChiHH600");

      h_J1M_doubletagSR_TChiHH800 = (TH1F*)fSignal->Get("J1pt_M_doubletagSR_TChiHH800"); h_J2M_doubletagSR_TChiHH800 = (TH1F*)fSignal->Get("J2pt_M_doubletagSR_TChiHH800");
      h_J1M_doubletagSB_TChiHH800 = (TH1F*)fSignal->Get("J1pt_M_doubletagSB_TChiHH800"); h_J2M_doubletagSB_TChiHH800 = (TH1F*)fSignal->Get("J2pt_M_doubletagSB_TChiHH800");
      h_J1M_antitagSR_TChiHH800 = (TH1F*)fSignal->Get("J1pt_M_antitagSR_TChiHH800"); h_J2M_antitagSR_TChiHH800 = (TH1F*)fSignal->Get("J2pt_M_antitagSR_TChiHH800");
      h_J1M_antitagSB_TChiHH800 = (TH1F*)fSignal->Get("J1pt_M_antitagSB_TChiHH800"); h_J2M_antitagSB_TChiHH800 = (TH1F*)fSignal->Get("J2pt_M_antitagSB_TChiHH800");
      h_J1M_tagSR_TChiHH800 = (TH1F*)fSignal->Get("J1pt_M_tagSR_TChiHH800"); h_J2M_tagSR_TChiHH800 = (TH1F*)fSignal->Get("J2pt_M_tagSR_TChiHH800");
      h_J1M_tagSB_TChiHH800 = (TH1F*)fSignal->Get("J1pt_M_tagSB_TChiHH800"); h_J2M_tagSB_TChiHH800 = (TH1F*)fSignal->Get("J2pt_M_tagSB_TChiHH800");
      h_baseline_j1mass_TChiHH800 = (TH1F*)fSignal->Get("J1pt_M_baseline_TChiHH800"); h_baseline_j2mass_TChiHH800 = (TH1F*)fSignal->Get("J2pt_M_baseline_TChiHH800");

      h_J1M_doubletagSR_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_M_doubletagSR_TChiHH1000"); h_J2M_doubletagSR_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_M_doubletagSR_TChiHH1000");
      h_J1M_doubletagSB_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_M_doubletagSB_TChiHH1000"); h_J2M_doubletagSB_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_M_doubletagSB_TChiHH1000");
      h_J1M_antitagSR_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_M_antitagSR_TChiHH1000"); h_J2M_antitagSR_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_M_antitagSR_TChiHH1000");
      h_J1M_antitagSB_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_M_antitagSB_TChiHH1000"); h_J2M_antitagSB_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_M_antitagSB_TChiHH1000");
      h_J1M_tagSR_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_M_tagSR_TChiHH1000"); h_J2M_tagSR_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_M_tagSR_TChiHH1000");
      h_J1M_tagSB_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_M_tagSB_TChiHH1000"); h_J2M_tagSB_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_M_tagSB_TChiHH1000");
      h_baseline_j1mass_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_M_baseline_TChiHH1000"); h_baseline_j2mass_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_M_baseline_TChiHH1000");
    }
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
    h_baseline_j1pt_ZJets = (TH1F*)f->Get("J1pt_Pt_baseline_ZJets"); h_baseline_j2pt_ZJets = (TH1F*)f->Get("J2pt_Pt_baseline_ZJets");

    h_J1Pt_doubletagSR_QCD = (TH1F*)f->Get("J1pt_Pt_doubletagSR_QCD"); h_J2Pt_doubletagSR_QCD = (TH1F*)f->Get("J2pt_Pt_doubletagSR_QCD");
    h_J1Pt_doubletagSB_QCD = (TH1F*)f->Get("J1pt_Pt_doubletagSB_QCD"); h_J2Pt_doubletagSB_QCD = (TH1F*)f->Get("J2pt_Pt_doubletagSB_QCD");
    h_J1Pt_antitagSR_QCD = (TH1F*)f->Get("J1pt_Pt_antitagSR_QCD"); h_J2Pt_antitagSR_QCD = (TH1F*)f->Get("J2pt_Pt_antitagSR_QCD");
    h_J1Pt_antitagSB_QCD = (TH1F*)f->Get("J1pt_Pt_antitagSB_QCD"); h_J2Pt_antitagSB_QCD = (TH1F*)f->Get("J2pt_Pt_antitagSB_QCD");
    h_J1Pt_tagSR_QCD = (TH1F*)f->Get("J1pt_Pt_tagSR_QCD"); h_J2Pt_tagSR_QCD = (TH1F*)f->Get("J2pt_Pt_tagSR_QCD");
    h_J1Pt_tagSB_QCD = (TH1F*)f->Get("J1pt_Pt_tagSB_QCD"); h_J2Pt_tagSB_QCD = (TH1F*)f->Get("J2pt_Pt_tagSB_QCD");
    h_baseline_j1pt_QCD = (TH1F*)f->Get("J1pt_Pt_baseline_QCD"); h_baseline_j2pt_QCD = (TH1F*)f->Get("J2pt_Pt_baseline_QCD");


    h_J1Pt_doubletagSR_SnglT = (TH1F*)f->Get("J1pt_Pt_doubletagSR_SnglT"); h_J2Pt_doubletagSR_SnglT = (TH1F*)f->Get("J2pt_Pt_doubletagSR_SnglT");
    h_J1Pt_doubletagSB_SnglT = (TH1F*)f->Get("J1pt_Pt_doubletagSB_SnglT"); h_J2Pt_doubletagSB_SnglT = (TH1F*)f->Get("J2pt_Pt_doubletagSB_SnglT");
    h_J1Pt_antitagSR_SnglT = (TH1F*)f->Get("J1pt_Pt_antitagSR_SnglT"); h_J2Pt_antitagSR_SnglT = (TH1F*)f->Get("J2pt_Pt_antitagSR_SnglT");
    h_J1Pt_antitagSB_SnglT = (TH1F*)f->Get("J1pt_Pt_antitagSB_SnglT"); h_J2Pt_antitagSB_SnglT = (TH1F*)f->Get("J2pt_Pt_antitagSB_SnglT");
    h_J1Pt_tagSR_SnglT = (TH1F*)f->Get("J1pt_Pt_tagSR_SnglT"); h_J2Pt_tagSR_SnglT = (TH1F*)f->Get("J2pt_Pt_tagSR_SnglT");
    h_J1Pt_tagSB_SnglT = (TH1F*)f->Get("J1pt_Pt_tagSB_SnglT"); h_J2Pt_tagSB_SnglT = (TH1F*)f->Get("J2pt_Pt_tagSB_SnglT");
    h_baseline_j1pt_SnglT = (TH1F*)f->Get("J1pt_Pt_baseline_SnglT"); h_baseline_j2pt_SnglT = (TH1F*)f->Get("J2pt_Pt_baseline_SnglT");


    h_J1Pt_doubletagSR_TT = (TH1F*)f->Get("J1pt_Pt_doubletagSR_TT"); h_J2Pt_doubletagSR_TT = (TH1F*)f->Get("J2pt_Pt_doubletagSR_TT");
    h_J1Pt_doubletagSB_TT = (TH1F*)f->Get("J1pt_Pt_doubletagSB_TT"); h_J2Pt_doubletagSB_TT = (TH1F*)f->Get("J2pt_Pt_doubletagSB_TT");
    h_J1Pt_antitagSR_TT = (TH1F*)f->Get("J1pt_Pt_antitagSR_TT"); h_J2Pt_antitagSR_TT = (TH1F*)f->Get("J2pt_Pt_antitagSR_TT");
    h_J1Pt_antitagSB_TT = (TH1F*)f->Get("J1pt_Pt_antitagSB_TT"); h_J2Pt_antitagSB_TT = (TH1F*)f->Get("J2pt_Pt_antitagSB_TT");
    h_J1Pt_tagSR_TT = (TH1F*)f->Get("J1pt_Pt_tagSR_TT"); h_J2Pt_tagSR_TT = (TH1F*)f->Get("J2pt_Pt_tagSR_TT");
    h_J1Pt_tagSB_TT = (TH1F*)f->Get("J1pt_Pt_tagSB_TT"); h_J2Pt_tagSB_TT = (TH1F*)f->Get("J2pt_Pt_tagSB_TT");
    h_baseline_j1pt_TT = (TH1F*)f->Get("J1pt_Pt_baseline_TT"); h_baseline_j2pt_TT = (TH1F*)f->Get("J2pt_Pt_baseline_TT");


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
    h_baseline_j1pt_WJets = (TH1F*)f->Get("J1pt_Pt_baseline_WJets"); h_baseline_j2pt_WJets = (TH1F*)f->Get("J2pt_Pt_baseline_WJets");


    //for signal
    if (0){
      h_J1Pt_doubletagSR_T5HH1600 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSR_T5qqqqZH1600"); h_J2Pt_doubletagSR_T5HH1600 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSR_T5qqqqZH1600");
      h_J1Pt_doubletagSB_T5HH1600 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSB_T5qqqqZH1600"); h_J2Pt_doubletagSB_T5HH1600 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSB_T5qqqqZH1600");
      h_J1Pt_antitagSR_T5HH1600 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSR_T5qqqqZH1600"); h_J2Pt_antitagSR_T5HH1600 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSR_T5qqqqZH1600");
      h_J1Pt_antitagSB_T5HH1600 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSB_T5qqqqZH1600"); h_J2Pt_antitagSB_T5HH1600 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSB_T5qqqqZH1600");
      h_J1Pt_tagSR_T5HH1600 = (TH1F*)fSignal->Get("J1pt_Pt_tagSR_T5qqqqZH1600"); h_J2Pt_tagSR_T5HH1600 = (TH1F*)fSignal->Get("J2pt_Pt_tagSR_T5qqqqZH1600");
      h_J1Pt_tagSB_T5HH1600 = (TH1F*)fSignal->Get("J1pt_Pt_tagSB_T5qqqqZH1600"); h_J2Pt_tagSB_T5HH1600 = (TH1F*)fSignal->Get("J2pt_Pt_tagSB_T5qqqqZH1600");
      h_baseline_j1pt_T5HH1600 = (TH1F*)fSignal->Get("J1pt_Pt_baseline_T5qqqqZH1600"); h_baseline_j2pt_T5HH1600 = (TH1F*)fSignal->Get("J2pt_Pt_baseline_T5qqqqZH1600");

      h_J1Pt_doubletagSR_T5HH2000 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSR_T5qqqqZH2000"); h_J2Pt_doubletagSR_T5HH2000 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSR_T5qqqqZH2000");
      h_J1Pt_doubletagSB_T5HH2000 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSB_T5qqqqZH2000"); h_J2Pt_doubletagSB_T5HH2000 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSB_T5qqqqZH2000");
      h_J1Pt_antitagSR_T5HH2000 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSR_T5qqqqZH2000"); h_J2Pt_antitagSR_T5HH2000 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSR_T5qqqqZH2000");
      h_J1Pt_antitagSB_T5HH2000 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSB_T5qqqqZH2000"); h_J2Pt_antitagSB_T5HH2000 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSB_T5qqqqZH2000");
      h_J1Pt_tagSR_T5HH2000 = (TH1F*)fSignal->Get("J1pt_Pt_tagSR_T5qqqqZH2000"); h_J2Pt_tagSR_T5HH2000 = (TH1F*)fSignal->Get("J2pt_Pt_tagSR_T5qqqqZH2000");
      h_J1Pt_tagSB_T5HH2000 = (TH1F*)fSignal->Get("J1pt_Pt_tagSB_T5qqqqZH2000"); h_J2Pt_tagSB_T5HH2000 = (TH1F*)fSignal->Get("J2pt_Pt_tagSB_T5qqqqZH2000");
      h_baseline_j1pt_T5HH2000 = (TH1F*)fSignal->Get("J1pt_Pt_baseline_T5qqqqZH2000"); h_baseline_j2pt_T5HH2000 = (TH1F*)fSignal->Get("J2pt_Pt_baseline_T5qqqqZH2000");

      h_J1Pt_doubletagSR_T5HH2200 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSR_T5qqqqZH2200"); h_J2Pt_doubletagSR_T5HH2200 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSR_T5qqqqZH2200");
      h_J1Pt_doubletagSB_T5HH2200 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSB_T5qqqqZH2200"); h_J2Pt_doubletagSB_T5HH2200 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSB_T5qqqqZH2200");
      h_J1Pt_antitagSR_T5HH2200 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSR_T5qqqqZH2200"); h_J2Pt_antitagSR_T5HH2200 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSR_T5qqqqZH2200");
      h_J1Pt_antitagSB_T5HH2200 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSB_T5qqqqZH2200"); h_J2Pt_antitagSB_T5HH2200 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSB_T5qqqqZH2200");
      h_J1Pt_tagSR_T5HH2200 = (TH1F*)fSignal->Get("J1pt_Pt_tagSR_T5qqqqZH2200"); h_J2Pt_tagSR_T5HH2200 = (TH1F*)fSignal->Get("J2pt_Pt_tagSR_T5qqqqZH2200");
      h_J1Pt_tagSB_T5HH2200 = (TH1F*)fSignal->Get("J1pt_Pt_tagSB_T5qqqqZH2200"); h_J2Pt_tagSB_T5HH2200 = (TH1F*)fSignal->Get("J2pt_Pt_tagSB_T5qqqqZH2200");
      h_baseline_j1pt_T5HH2200 = (TH1F*)fSignal->Get("J1pt_Pt_baseline_T5qqqqZH2200"); h_baseline_j2pt_T5HH2200 = (TH1F*)fSignal->Get("J2pt_Pt_baseline_T5qqqqZH2200");

      h_J1Pt_doubletagSR_TChiHH600 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSR_TChiHH600"); h_J2Pt_doubletagSR_TChiHH600 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSR_TChiHH600");
      h_J1Pt_doubletagSB_TChiHH600 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSB_TChiHH600"); h_J2Pt_doubletagSB_TChiHH600 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSB_TChiHH600");
      h_J1Pt_antitagSR_TChiHH600 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSR_TChiHH600"); h_J2Pt_antitagSR_TChiHH600 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSR_TChiHH600");
      h_J1Pt_antitagSB_TChiHH600 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSB_TChiHH600"); h_J2Pt_antitagSB_TChiHH600 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSB_TChiHH600");
      h_J1Pt_tagSR_TChiHH600 = (TH1F*)fSignal->Get("J1pt_Pt_tagSR_TChiHH600"); h_J2Pt_tagSR_TChiHH600 = (TH1F*)fSignal->Get("J2pt_Pt_tagSR_TChiHH600");
      h_J1Pt_tagSB_TChiHH600 = (TH1F*)fSignal->Get("J1pt_Pt_tagSB_TChiHH600"); h_J2Pt_tagSB_TChiHH600 = (TH1F*)fSignal->Get("J2pt_Pt_tagSB_TChiHH600");
      h_baseline_j1pt_TChiHH600 = (TH1F*)fSignal->Get("J1pt_Pt_baseline_TChiHH600"); h_baseline_j2pt_TChiHH600 = (TH1F*)fSignal->Get("J2pt_Pt_baseline_TChiHH600");

      h_J1Pt_doubletagSR_TChiHH800 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSR_TChiHH800"); h_J2Pt_doubletagSR_TChiHH800 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSR_TChiHH800");
      h_J1Pt_doubletagSB_TChiHH800 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSB_TChiHH800"); h_J2Pt_doubletagSB_TChiHH800 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSB_TChiHH800");
      h_J1Pt_antitagSR_TChiHH800 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSR_TChiHH800"); h_J2Pt_antitagSR_TChiHH800 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSR_TChiHH800");
      h_J1Pt_antitagSB_TChiHH800 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSB_TChiHH800"); h_J2Pt_antitagSB_TChiHH800 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSB_TChiHH800");
      h_J1Pt_tagSR_TChiHH800 = (TH1F*)fSignal->Get("J1pt_Pt_tagSR_TChiHH800"); h_J2Pt_tagSR_TChiHH800 = (TH1F*)fSignal->Get("J2pt_Pt_tagSR_TChiHH800");
      h_J1Pt_tagSB_TChiHH800 = (TH1F*)fSignal->Get("J1pt_Pt_tagSB_TChiHH800"); h_J2Pt_tagSB_TChiHH800 = (TH1F*)fSignal->Get("J2pt_Pt_tagSB_TChiHH800");
      h_baseline_j1pt_TChiHH800 = (TH1F*)fSignal->Get("J1pt_Pt_baseline_TChiHH800"); h_baseline_j2pt_TChiHH800 = (TH1F*)fSignal->Get("J2pt_Pt_baseline_TChiHH800");

      h_J1Pt_doubletagSR_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSR_TChiHH1000"); h_J2Pt_doubletagSR_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSR_TChiHH1000");
      h_J1Pt_doubletagSB_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_Pt_doubletagSB_TChiHH1000"); h_J2Pt_doubletagSB_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_Pt_doubletagSB_TChiHH1000");
      h_J1Pt_antitagSR_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSR_TChiHH1000"); h_J2Pt_antitagSR_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSR_TChiHH1000");
      h_J1Pt_antitagSB_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_Pt_antitagSB_TChiHH1000"); h_J2Pt_antitagSB_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_Pt_antitagSB_TChiHH1000");
      h_J1Pt_tagSR_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_Pt_tagSR_TChiHH1000"); h_J2Pt_tagSR_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_Pt_tagSR_TChiHH1000");
      h_J1Pt_tagSB_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_Pt_tagSB_TChiHH1000"); h_J2Pt_tagSB_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_Pt_tagSB_TChiHH1000");
      h_baseline_j1pt_TChiHH1000 = (TH1F*)fSignal->Get("J1pt_Pt_baseline_TChiHH1000"); h_baseline_j2pt_TChiHH1000 = (TH1F*)fSignal->Get("J2pt_Pt_baseline_TChiHH1000");
    }
    // h_COpt1_sum = (TH1F*)f->Get("MET_antitagSROpt1_sum"); h_DOpt1_sum = (TH1F*)f->Get("MET_antitagSBOpt1_sum");
    // h_COpt2_sum = (TH1F*)f->Get("MET_antitagSROpt2_sum"); h_DOpt2_sum = (TH1F*)f->Get("MET_antitagSBOpt2_sum");
    // h_COpt3_sum = (TH1F*)f->Get("MET_antitagSROpt3_sum"); h_DOpt3_sum = (TH1F*)f->Get("MET_antitagSBOpt3_sum");
    // h_COpt4_sum = (TH1F*)f->Get("MET_antitagSROpt4_sum"); h_DOpt4_sum = (TH1F*)f->Get("MET_antitagSBOpt4_sum");
  }

  if (whichRegion=="photon") {
    h_A_sum = (TH1F*)fPhoton->Get("MET_doubletagSR_sum"); h_B_sum = (TH1F*)fPhoton->Get("MET_doubletagSB_sum");
    h_A1_sum = (TH1F*)fPhoton->Get("MET_tagSR_sum"); h_B1_sum = (TH1F*)fPhoton->Get("MET_tagSB_sum");
    h_C_sum = (TH1F*)fPhoton->Get("MET_antitagSR_sum"); h_D_sum = (TH1F*)fPhoton->Get("MET_antitagSB_sum");

    hP_A_sum = (TH1F*)fPhoton->Get("METPhoton_doubletagSR_sum"); hP_B_sum = (TH1F*)fPhoton->Get("METPhoton_doubletagSB_sum");
    hP_A1_sum = (TH1F*)fPhoton->Get("METPhoton_tagSR_sum"); hP_B1_sum = (TH1F*)fPhoton->Get("METPhoton_tagSB_sum");
    hP_C_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSR_sum"); hP_D_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSB_sum");

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

    h_J1J2M_doubletagSR_sum = (TH1F*)h_J1M_doubletagSR_sum->Clone("J1J2M_doubletagSR_sum"); h_J1J2M_doubletagSR_sum->Add(h_J2M_doubletagSR_sum);
    h_J1J2M_doubletagSB_sum = (TH1F*)h_J1M_doubletagSB_sum->Clone("J1J2M_doubletagSB_sum"); h_J1J2M_doubletagSB_sum->Add(h_J2M_doubletagSB_sum);
    h_J1J2M_tagSR_sum = (TH1F*)h_J1M_tagSR_sum->Clone("J1J2M_tagSR_sum"); h_J1J2M_tagSR_sum->Add(h_J2M_tagSR_sum);
    h_J1J2M_tagSB_sum = (TH1F*)h_J1M_tagSB_sum->Clone("J1J2M_tagSB_sum"); h_J1J2M_tagSB_sum->Add(h_J2M_tagSB_sum);
    h_J1J2M_antitagSR_sum = (TH1F*)h_J1M_antitagSR_sum->Clone("J1J2M_antitagSR_sum"); h_J1J2M_antitagSR_sum->Add(h_J2M_antitagSR_sum);
    h_J1J2M_antitagSB_sum = (TH1F*)h_J1M_antitagSB_sum->Clone("J1J2M_antitagSB_sum"); h_J1J2M_antitagSB_sum->Add(h_J2M_antitagSB_sum);

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

    // h_DOpt1_sum = (TH1F*)fPhoton->Get("MET_antitagSBOpt1_sum"); h_DOpt2_sum = (TH1F*)fPhoton->Get("MET_antitagSBOpt2_sum");
    // h_COpt1_sum = (TH1F*)fPhoton->Get("MET_antitagSROpt1_sum"); h_COpt2_sum = (TH1F*)fPhoton->Get("MET_antitagSROpt2_sum");
    h_COpt1_QCD = (TH1F*)fPhoton->Get("MET_antitagSROpt1_QCD"); h_DOpt1_QCD = (TH1F*)fPhoton->Get("MET_antitagSBOpt1_QCD");
    h_COpt1_GJets = (TH1F*)fPhoton->Get("MET_antitagSROpt1_GJets"); h_DOpt1_GJets = (TH1F*)fPhoton->Get("MET_antitagSBOpt1_GJets");
    h_COpt1_sum = (TH1F*)fPhoton->Get("MET_antitagSROpt1_sum"); h_DOpt1_sum = (TH1F*)fPhoton->Get("MET_antitagSBOpt1_sum");
    hP_COpt1_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSROpt1_sum"); hP_DOpt1_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSBOpt1_sum");

    h_C5Opt1_sum = (TH1F*)fPhoton->Get("MET5_antitagSROpt1_sum"); h_D5Opt1_sum = (TH1F*)fPhoton->Get("MET5_antitagSBOpt1_sum");
    h_C5Opt2_sum = (TH1F*)fPhoton->Get("MET5_antitagSROpt2_sum"); h_D5Opt2_sum = (TH1F*)fPhoton->Get("MET5_antitagSBOpt2_sum");
    h_C5Opt3_sum = (TH1F*)fPhoton->Get("MET5_antitagSROpt3_sum"); h_D5Opt3_sum = (TH1F*)fPhoton->Get("MET5_antitagSBOpt3_sum");
    h_C5Opt4_sum = (TH1F*)fPhoton->Get("MET5_antitagSROpt4_sum"); h_D5Opt4_sum = (TH1F*)fPhoton->Get("MET5_antitagSBOpt4_sum");

    h_COpt2_QCD = (TH1F*)fPhoton->Get("MET_antitagSROpt2_QCD"); h_DOpt2_QCD = (TH1F*)fPhoton->Get("MET_antitagSBOpt2_QCD");
    h_COpt2_GJets = (TH1F*)fPhoton->Get("MET_antitagSROpt2_GJets"); h_DOpt2_GJets = (TH1F*)fPhoton->Get("MET_antitagSBOpt2_GJets");
    h_COpt2_sum = (TH1F*)fPhoton->Get("MET_antitagSROpt2_sum"); h_DOpt2_sum = (TH1F*)fPhoton->Get("MET_antitagSBOpt2_sum");
    hP_COpt2_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSROpt2_sum"); hP_DOpt2_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSBOpt2_sum");

    h_COpt3_QCD = (TH1F*)fPhoton->Get("MET_antitagSROpt3_QCD"); h_DOpt3_QCD = (TH1F*)fPhoton->Get("MET_antitagSBOpt3_QCD");
    h_COpt3_GJets = (TH1F*)fPhoton->Get("MET_antitagSROpt3_GJets"); h_DOpt3_GJets = (TH1F*)fPhoton->Get("MET_antitagSBOpt3_GJets");
    h_COpt3_sum = (TH1F*)fPhoton->Get("MET_antitagSROpt3_sum"); h_DOpt3_sum = (TH1F*)fPhoton->Get("MET_antitagSBOpt3_sum");
    hP_COpt3_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSROpt3_sum"); hP_DOpt3_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSBOpt3_sum");

    h_COpt4_QCD = (TH1F*)fPhoton->Get("MET_antitagSROpt4_QCD"); h_DOpt4_QCD = (TH1F*)fPhoton->Get("MET_antitagSBOpt4_QCD");
    h_COpt4_GJets = (TH1F*)fPhoton->Get("MET_antitagSROpt4_GJets"); h_DOpt4_GJets = (TH1F*)fPhoton->Get("MET_antitagSBOpt4_GJets");
    h_COpt4_sum = (TH1F*)fPhoton->Get("MET_antitagSROpt4_sum"); h_DOpt4_sum = (TH1F*)fPhoton->Get("MET_antitagSBOpt4_sum");
    hP_COpt4_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSROpt4_sum"); hP_DOpt4_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSBOpt4_sum");

    h_COpt5_QCD = (TH1F*)fPhoton->Get("MET_antitagSROpt5_QCD"); h_DOpt5_QCD = (TH1F*)fPhoton->Get("MET_antitagSBOpt5_QCD");
    h_COpt5_GJets = (TH1F*)fPhoton->Get("MET_antitagSROpt5_GJets"); h_DOpt5_GJets = (TH1F*)fPhoton->Get("MET_antitagSBOpt5_GJets");
    h_COpt5_sum = (TH1F*)fPhoton->Get("MET_antitagSROpt5_sum"); h_DOpt5_sum = (TH1F*)fPhoton->Get("MET_antitagSBOpt5_sum");
    hP_COpt5_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSROpt5_sum"); hP_DOpt5_sum = (TH1F*)fPhoton->Get("METPhoton_antitagSBOpt5_sum");
  }


  if (whichRegion=="singleLept") {
    if (runDM) {
      h_avgM_A_sum = (TH1F*)fSingleLept->Get("avgM_doubletagSR_sum");
      h_avgM_B_sum = (TH1F*)fSingleLept->Get("avgM_doubletagSB_sum");
      h_avgM_A1_sum = (TH1F*)fSingleLept->Get("avgM_tagSR_sum");
      h_avgM_B1_sum = (TH1F*)fSingleLept->Get("avgM_tagSB_sum");
      h_avgM_C_sum = (TH1F*)fSingleLept->Get("avgM_antitagSR_sum");
      h_avgM_D_sum = (TH1F*)fSingleLept->Get("avgM_antitagSB_sum");
    }

    // h_A_data = (TH1F*)fSingleLept->Get("MET_doubletagSR_data"); h_B_data = (TH1F*)fSingleLept->Get("MET_doubletagSB_data");
    // h_A1_data = (TH1F*)fSingleLept->Get("MET_tagSR_data"); h_B1_data = (TH1F*)fSingleLept->Get("MET_tagSB_data");
    // h_C_data = (TH1F*)fSingleLept->Get("MET_antitagSR_data"); h_D_data = (TH1F*)fSingleLept->Get("MET_antitagSB_data");


    // h_A_sum = (TH1F*)fSingleLept->Get("MET_doubletagSR_sum"); h_B_sum = (TH1F*)fSingleLept->Get("MET_doubletagSB_sum");
    // h_A1_sum = (TH1F*)fSingleLept->Get("MET_tagSR_sum"); h_B1_sum = (TH1F*)fSingleLept->Get("MET_tagSB_sum");
    // h_C_sum = (TH1F*)fSingleLept->Get("MET_antitagSR_sum"); h_D_sum = (TH1F*)fSingleLept->Get("MET_antitagSB_sum");

    h_A5_sum = (TH1F*)fSingleLept->Get("MET5_doubletagSR_sum"); h_B5_sum = (TH1F*)fSingleLept->Get("MET5_doubletagSB_sum");
    h_A15_sum = (TH1F*)fSingleLept->Get("MET5_tagSR_sum"); h_B15_sum = (TH1F*)fSingleLept->Get("MET5_tagSB_sum");
    h_C5_sum = (TH1F*)fSingleLept->Get("MET5_antitagSR_sum"); h_D5_sum = (TH1F*)fSingleLept->Get("MET5_antitagSB_sum");

    h_A_SnglT = (TH1F*)fSingleLept->Get("MET_doubletagSR_SnglT"); h_B_SnglT = (TH1F*)fSingleLept->Get("MET_doubletagSB_SnglT");
    h_A1_SnglT = (TH1F*)fSingleLept->Get("MET_tagSR_SnglT"); h_B1_SnglT = (TH1F*)fSingleLept->Get("MET_tagSB_SnglT");
    h_C_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSR_SnglT"); h_D_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSB_SnglT");

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

    h_COpt1_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_SnglT"); h_DOpt1_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_SnglT");
    h_COpt1_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_TT"); h_DOpt1_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_TT");
    h_COpt1_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_WJets"); h_DOpt1_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_WJets");
    h_COpt2_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSROpt2_SnglT"); h_DOpt2_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt2_SnglT");
    h_COpt2_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt2_TT"); h_DOpt2_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt2_TT");
    h_COpt2_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt2_WJets"); h_DOpt2_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt2_WJets");
    h_COpt3_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSROpt3_SnglT"); h_DOpt3_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt3_SnglT");
    h_COpt3_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt3_TT"); h_DOpt3_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt3_TT");
    h_COpt3_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt3_WJets"); h_DOpt3_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt3_WJets");
    h_COpt4_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSROpt4_SnglT"); h_DOpt4_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt4_SnglT");
    h_COpt4_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt4_TT"); h_DOpt4_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt4_TT");
    h_COpt4_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt4_WJets"); h_DOpt4_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt4_WJets");
    h_COpt5_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSROpt5_SnglT"); h_DOpt5_SnglT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt5_SnglT");
    h_COpt5_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt5_TT"); h_DOpt5_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt5_TT");
    h_COpt5_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt5_WJets"); h_DOpt5_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt5_WJets");

    h_C5Opt1_sum = (TH1F*)fSingleLept->Get("MET5_antitagSROpt1_sum"); h_D5Opt1_sum = (TH1F*)fSingleLept->Get("MET5_antitagSBOpt1_sum");
    h_C5Opt2_sum = (TH1F*)fSingleLept->Get("MET5_antitagSROpt2_sum"); h_D5Opt2_sum = (TH1F*)fSingleLept->Get("MET5_antitagSBOpt2_sum");
    h_C5Opt3_sum = (TH1F*)fSingleLept->Get("MET5_antitagSROpt3_sum"); h_D5Opt3_sum = (TH1F*)fSingleLept->Get("MET5_antitagSBOpt3_sum");
    h_C5Opt4_sum = (TH1F*)fSingleLept->Get("MET5_antitagSROpt4_sum"); h_D5Opt4_sum = (TH1F*)fSingleLept->Get("MET5_antitagSBOpt4_sum");

    // h_A_sum = (TH1F*)h_A_TT->Clone("h_A_sum");  h_B_sum = (TH1F*)h_B_TT->Clone("h_B_sum");
    // h_A1_sum = (TH1F*)h_A1_TT->Clone("h_A1_sum");  h_B1_sum = (TH1F*)h_B1_TT->Clone("h_B1_sum");
    // h_C_sum = (TH1F*)h_C_TT->Clone("h_C_sum");  h_D_sum = (TH1F*)h_D_TT->Clone("h_D_sum");
    // h_A_sum->Add(h_A_WJets); h_B_sum->Add(h_B_WJets);
    // h_A1_sum->Add(h_A1_WJets); h_B1_sum->Add(h_B1_WJets);
    // h_C_sum->Add(h_C_WJets); h_D_sum->Add(h_D_WJets);
    //Adding uncertainties to 0-event bins

    //TH1F make0EventUncSum_1l(vector<TH1F*> dem_histos); //in this order: SnglT,TT, WJets
    h_A_sum = make0EventUncSum_1l({h_A_SnglT,h_A_TT,h_A_WJets});
    h_B_sum = make0EventUncSum_1l({h_B_SnglT,h_B_TT,h_B_WJets});
    h_A1_sum = make0EventUncSum_1l({h_A1_SnglT,h_A1_TT,h_A1_WJets});
    h_B1_sum = make0EventUncSum_1l({h_B1_SnglT,h_B1_TT,h_B1_WJets});
    h_C_sum = make0EventUncSum_1l({h_C_SnglT,h_C_TT,h_C_WJets});
    h_D_sum = make0EventUncSum_1l({h_D_SnglT,h_D_TT,h_D_WJets,});


    h_A_sum->Write("ABCDRegions/h_A_sum");
    h_A1_sum->Write("ABCDRegions/h_A1_sum");
    h_B_sum->Write("ABCDRegions/h_B_sum");
    h_B1_sum->Write("ABCDRegions/h_B1_sum");
    h_C_sum->Write("ABCDRegions/h_C_sum");
    h_D_sum->Write("ABCDRegions/h_D_sum");


    h_COpt1_sum = make0EventUncSum_1l({h_COpt1_SnglT,h_COpt1_TT,h_COpt1_WJets});
    h_DOpt1_sum = make0EventUncSum_1l({h_DOpt1_SnglT,h_DOpt1_TT,h_DOpt1_WJets});
    h_COpt2_sum = make0EventUncSum_1l({h_COpt2_SnglT,h_COpt2_TT,h_COpt2_WJets});
    h_DOpt2_sum = make0EventUncSum_1l({h_DOpt2_SnglT,h_DOpt2_TT,h_DOpt2_WJets});
    h_COpt3_sum = make0EventUncSum_1l({h_COpt3_SnglT,h_COpt3_TT,h_COpt3_WJets});
    h_DOpt3_sum = make0EventUncSum_1l({h_DOpt3_SnglT,h_DOpt3_TT,h_DOpt3_WJets});
    h_COpt4_sum = make0EventUncSum_1l({h_COpt4_SnglT,h_COpt4_TT,h_COpt4_WJets});
    h_DOpt4_sum = make0EventUncSum_1l({h_DOpt4_SnglT,h_DOpt4_TT,h_DOpt4_WJets});
    h_COpt5_sum = make0EventUncSum_1l({h_COpt5_SnglT,h_COpt5_TT,h_COpt5_WJets});
    h_DOpt5_sum = make0EventUncSum_1l({h_DOpt5_SnglT,h_DOpt5_TT,h_DOpt5_WJets});


    h_COpt1_sum->Write("ABCDRegions/h_C_M1_sum");
    h_COpt2_sum->Write("ABCDRegions/h_C_M0_sum");
    h_COpt3_sum->Write("ABCDRegions/h_C_T0_sum");
    h_COpt4_sum->Write("ABCDRegions/h_C_L1M0_sum");
    h_COpt5_sum->Write("ABCDRegions/h_C_M1T0_sum");
    h_DOpt1_sum->Write("ABCDRegions/h_D_M1_sum");
    h_DOpt2_sum->Write("ABCDRegions/h_D_M0_sum");
    h_DOpt3_sum->Write("ABCDRegions/h_D_T0_sum");
    h_DOpt4_sum->Write("ABCDRegions/h_D_L1M0_sum");
    h_DOpt5_sum->Write("ABCDRegions/h_D_M1T0_sum");


    h_J1M_doubletagSR_sum = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSR_sum"); h_J2M_doubletagSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSR_sum");
    h_J1M_doubletagSB_sum = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSB_sum"); h_J2M_doubletagSB_sum = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSB_sum");
    h_J1M_antitagSR_sum = (TH1F*)fSingleLept->Get("J1pt_M_antitagSR_sum"); h_J2M_antitagSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_antitagSR_sum");
    h_J1M_antitagSB_sum = (TH1F*)fSingleLept->Get("J1pt_M_antitagSB_sum"); h_J2M_antitagSB_sum = (TH1F*)fSingleLept->Get("J2pt_M_antitagSB_sum");
    h_J1M_tagSR_sum = (TH1F*)fSingleLept->Get("J1pt_M_tagSR_sum"); h_J2M_tagSR_sum = (TH1F*)fSingleLept->Get("J2pt_M_tagSR_sum");
    h_J1M_tagSB_sum = (TH1F*)fSingleLept->Get("J1pt_M_tagSB_sum"); h_J2M_tagSB_sum = (TH1F*)fSingleLept->Get("J2pt_M_tagSB_sum");

    h_J1J2M_doubletagSR_sum = (TH1F*)h_J1M_doubletagSR_sum->Clone("J1J2M_doubletagSR_sum"); h_J1J2M_doubletagSR_sum->Add(h_J2M_doubletagSR_sum);
    h_J1J2M_doubletagSB_sum = (TH1F*)h_J1M_doubletagSB_sum->Clone("J1J2M_doubletagSB_sum"); h_J1J2M_doubletagSB_sum->Add(h_J2M_doubletagSB_sum);
    h_J1J2M_tagSR_sum = (TH1F*)h_J1M_tagSR_sum->Clone("J1J2M_tagSR_sum"); h_J1J2M_tagSR_sum->Add(h_J2M_tagSR_sum);
    h_J1J2M_tagSB_sum = (TH1F*)h_J1M_tagSB_sum->Clone("J1J2M_tagSB_sum"); h_J1J2M_tagSB_sum->Add(h_J2M_tagSB_sum);
    h_J1J2M_antitagSR_sum = (TH1F*)h_J1M_antitagSR_sum->Clone("J1J2M_antitagSR_sum"); h_J1J2M_antitagSR_sum->Add(h_J2M_antitagSR_sum);
    h_J1J2M_antitagSB_sum = (TH1F*)h_J1M_antitagSB_sum->Clone("J1J2M_antitagSB_sum"); h_J1J2M_antitagSB_sum->Add(h_J2M_antitagSB_sum);

    h_J1M_doubletagSR_SnglT = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSR_SnglT"); h_J2M_doubletagSR_SnglT = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSR_SnglT");
    h_J1M_doubletagSB_SnglT = (TH1F*)fSingleLept->Get("J1pt_M_doubletagSB_SnglT"); h_J2M_doubletagSB_SnglT = (TH1F*)fSingleLept->Get("J2pt_M_doubletagSB_SnglT");
    h_J1M_antitagSR_SnglT = (TH1F*)fSingleLept->Get("J1pt_M_antitagSR_SnglT"); h_J2M_antitagSR_SnglT = (TH1F*)fSingleLept->Get("J2pt_M_antitagSR_SnglT");
    h_J1M_antitagSB_SnglT = (TH1F*)fSingleLept->Get("J1pt_M_antitagSB_SnglT"); h_J2M_antitagSB_SnglT = (TH1F*)fSingleLept->Get("J2pt_M_antitagSB_SnglT");
    h_J1M_tagSR_SnglT = (TH1F*)fSingleLept->Get("J1pt_M_tagSR_SnglT"); h_J2M_tagSR_SnglT = (TH1F*)fSingleLept->Get("J2pt_M_tagSR_SnglT");
    h_J1M_tagSB_SnglT = (TH1F*)fSingleLept->Get("J1pt_M_tagSB_SnglT"); h_J2M_tagSB_SnglT = (TH1F*)fSingleLept->Get("J2pt_M_tagSB_SnglT");

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

    h_J1Pt_doubletagSR_SnglT = (TH1F*)fSingleLept->Get("J1pt_Pt_doubletagSR_SnglT"); h_J2Pt_doubletagSR_SnglT = (TH1F*)fSingleLept->Get("J2pt_Pt_doubletagSR_SnglT");
    h_J1Pt_doubletagSB_SnglT = (TH1F*)fSingleLept->Get("J1pt_Pt_doubletagSB_SnglT"); h_J2Pt_doubletagSB_SnglT = (TH1F*)fSingleLept->Get("J2pt_Pt_doubletagSB_SnglT");
    h_J1Pt_antitagSR_SnglT = (TH1F*)fSingleLept->Get("J1pt_Pt_antitagSR_SnglT"); h_J2Pt_antitagSR_SnglT = (TH1F*)fSingleLept->Get("J2pt_Pt_antitagSR_SnglT");
    h_J1Pt_antitagSB_SnglT = (TH1F*)fSingleLept->Get("J1pt_Pt_antitagSB_SnglT"); h_J2Pt_antitagSB_SnglT = (TH1F*)fSingleLept->Get("J2pt_Pt_antitagSB_SnglT");
    h_J1Pt_tagSR_SnglT = (TH1F*)fSingleLept->Get("J1pt_Pt_tagSR_SnglT"); h_J2Pt_tagSR_SnglT = (TH1F*)fSingleLept->Get("J2pt_Pt_tagSR_SnglT");
    h_J1Pt_tagSB_SnglT = (TH1F*)fSingleLept->Get("J1pt_Pt_tagSB_SnglT"); h_J2Pt_tagSB_SnglT = (TH1F*)fSingleLept->Get("J2pt_Pt_tagSB_SnglT");


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
    h_J2M_mjBins_doubletagSR_SnglT = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSR_SnglT"); h_J2M_mjBins_doubletagSB_SnglT = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSB_SnglT");
    h_J2M_mjBins_antitagSR_SnglT = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSR_SnglT"); h_J2M_mjBins_antitagSB_SnglT = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSB_SnglT");

    // h_J2M_mjBins_doubletagSR_TT_Di = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_TT_Di"); h_J2M_mjBins_doubletagSB_TT_Di = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_TT_Di");
    // h_J2M_mjBins_antitagSR_TT_Di = (TH1F*)f->Get("J2_M_jetBins_antitagSR_TT_Di"); h_J2M_mjBins_antitagSB_TT_Di = (TH1F*)f->Get("J2_M_jetBins_antitagSB_TT_Di");
    // h_J2M_mjBins_doubletagSR_TT_SL = (TH1F*)f->Get("J2_M_jetBins_doubletagSR_TT_SL"); h_J2M_mjBins_doubletagSB_TT_SL = (TH1F*)f->Get("J2_M_jetBins_doubletagSB_TT_SL");
    // h_J2M_mjBins_antitagSR_TT_SL = (TH1F*)f->Get("J2_M_jetBins_antitagSR_TT_SL"); h_J2M_mjBins_antitagSB_TT_SL = (TH1F*)f->Get("J2_M_jetBins_antitagSB_TT_SL");

    h_J2M_mjBins_doubletagSR_WJets = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSR_WJets"); h_J2M_mjBins_doubletagSB_WJets = (TH1F*)fSingleLept->Get("J2_M_jetBins_doubletagSB_WJets");
    h_J2M_mjBins_antitagSR_WJets = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSR_WJets"); h_J2M_mjBins_antitagSB_WJets = (TH1F*)fSingleLept->Get("J2_M_jetBins_antitagSB_WJets");
    //
    // h_DOpt1_sum = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_sum"); h_COpt1_sum = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_sum");
    // h_DOpt2_sum = (TH1F*)fSingleLept->Get("MET_antitagSBOpt2_sum"); h_COpt2_sum = (TH1F*)fSingleLept->Get("MET_antitagSROpt2_sum");
    // h_DOpt3_sum = (TH1F*)fSingleLept->Get("MET_antitagSBOpt3_sum"); h_COpt3_sum = (TH1F*)fSingleLept->Get("MET_antitagSROpt3_sum");
    // h_DOpt4_sum = (TH1F*)fSingleLept->Get("MET_antitagSBOpt4_sum"); h_COpt4_sum = (TH1F*)fSingleLept->Get("MET_antitagSROpt4_sum");

    // h_DOpt1_TT = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_TT"); h_DOpt1_WJets = (TH1F*)fSingleLept->Get("MET_antitagSBOpt1_WJets");
    // h_COpt1_TT = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_TT"); h_COpt1_WJets = (TH1F*)fSingleLept->Get("MET_antitagSROpt1_WJets");
  }

  vector<TH1F*> histos_ABCD_data = {h_A_data, h_B_data, h_C_data, h_D_data};
  vector<TH1F*> histos_A1B1CD_data = {h_A1_data, h_B1_data, h_C_data, h_D_data};

  vector<TH1F*> histos_allRegions_sum = {h_A_sum, h_A1_sum, h_B_sum, h_B1_sum, h_C_sum, h_D_sum};

  vector<TH1F*> histosPhoton_ABCD_sum = {hP_A_sum, hP_B_sum, hP_C_sum, hP_D_sum};
  vector<TH1F*> histosPhoton_A1B1CD_sum = {hP_A1_sum, hP_B1_sum, hP_C_sum, hP_D_sum};


  vector<TH1F*> histos_ABCD_sum = {h_A_sum, h_B_sum, h_C_sum, h_D_sum};
  vector<TH1F*> histos_ABCDlowPU_sum = {h_AlowPU_sum, h_BlowPU_sum, h_ClowPU_sum, h_DlowPU_sum};
  vector<TH1F*> histos_ABCDhighPU_sum = {h_AhighPU_sum, h_BhighPU_sum, h_ChighPU_sum, h_DhighPU_sum};
  vector<TH1F*> histos_ABCD_QCD = {h_A_QCD, h_B_QCD, h_C_QCD, h_D_QCD};
  vector<TH1F*> histos_ABCD_GJets = {h_A_GJets, h_B_GJets, h_C_GJets, h_D_GJets};
  vector<TH1F*> histos_ABCD_SnglT = {h_A_SnglT, h_B_SnglT, h_C_SnglT, h_D_SnglT};
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
  vector<TH1F*> histos_A1B1CDlowPU_sum = {h_A1lowPU_sum, h_B1lowPU_sum, h_ClowPU_sum, h_DlowPU_sum};
  vector<TH1F*> histos_A1B1CDhighPU_sum = {h_A1highPU_sum, h_B1highPU_sum, h_ChighPU_sum, h_DhighPU_sum};
  vector<TH1F*> histos_A1B1CD_QCD = {h_A1_QCD, h_B1_QCD, h_C_QCD, h_D_QCD};
  vector<TH1F*> histos_A1B1CD_GJets = {h_A1_GJets, h_B1_GJets, h_C_GJets, h_D_GJets};
  vector<TH1F*> histos_A1B1CD_SnglT = {h_A1_SnglT, h_B1_SnglT, h_C_SnglT, h_D_SnglT};
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

  vector<TH1F*> histos_RPF_DM = {h_avgM_A_sum,h_avgM_B_sum,h_avgM_C_sum,h_avgM_D_sum};
  vector<TH1F*> histos_RPFsingle_DM = {h_avgM_A1_sum,h_avgM_B1_sum,h_avgM_C_sum,h_avgM_D_sum};

  vector<TH1F*> histos_Rpf_J1_sum = {h_J1M_doubletagSR_sum, h_J1M_doubletagSB_sum, h_J1M_antitagSR_sum, h_J1M_antitagSB_sum};
  vector<TH1F*> histos_Rpf_J2_sum = {h_J2M_doubletagSR_sum, h_J2M_doubletagSB_sum, h_J2M_antitagSR_sum, h_J2M_antitagSB_sum};
  vector<TH1F*> histos_Rpf_J2_mJBins_sum = {h_J2M_mjBins_doubletagSR_sum, h_J2M_mjBins_doubletagSB_sum, h_J2M_mjBins_antitagSR_sum, h_J2M_mjBins_antitagSB_sum};
  vector<TH1F*> histos_Rpfsingle_J1_sum = {h_J1M_tagSR_sum, h_J1M_tagSB_sum, h_J1M_antitagSR_sum, h_J1M_antitagSB_sum};
  vector<TH1F*> histos_Rpfsingle_J2_sum = {h_J2M_tagSR_sum, h_J2M_tagSB_sum, h_J2M_antitagSR_sum, h_J2M_antitagSB_sum};

  vector<TH1F*> histos_Rpf_J1J2_sum = {h_J1J2M_doubletagSR_sum, h_J1J2M_doubletagSB_sum, h_J1J2M_antitagSR_sum, h_J1J2M_antitagSB_sum};
  vector<TH1F*> histos_Rpfsingle_J1J2_sum = {h_J1J2M_tagSR_sum, h_J1J2M_tagSB_sum, h_J1J2M_antitagSR_sum, h_J1J2M_antitagSB_sum};


  vector<TH1F*> histos_Rpf_J1_QCD = {h_J1M_doubletagSR_QCD, h_J1M_doubletagSB_QCD, h_J1M_antitagSR_QCD, h_J1M_antitagSB_QCD};
  vector<TH1F*> histos_Rpf_J2_QCD = {h_J2M_doubletagSR_QCD, h_J2M_doubletagSB_QCD, h_J2M_antitagSR_QCD, h_J2M_antitagSB_QCD};
  vector<TH1F*> histos_Rpf_J2_mJBins_QCD = {h_J2M_mjBins_doubletagSR_QCD, h_J2M_mjBins_doubletagSB_QCD, h_J2M_mjBins_antitagSR_QCD, h_J2M_mjBins_antitagSB_QCD};
  vector<TH1F*> histos_Rpfsingle_J1_QCD = {h_J1M_tagSR_QCD, h_J1M_tagSB_QCD, h_J1M_antitagSR_QCD, h_J1M_antitagSB_QCD};
  vector<TH1F*> histos_Rpfsingle_J2_QCD = {h_J2M_tagSR_QCD, h_J2M_tagSB_QCD, h_J2M_antitagSR_QCD, h_J2M_antitagSB_QCD};

  vector<TH1F*> histos_Rpf_J1_GJets = {h_J1M_doubletagSR_GJets, h_J1M_doubletagSB_GJets, h_J1M_antitagSR_GJets, h_J1M_antitagSB_GJets};
  vector<TH1F*> histos_Rpf_J2_GJets = {h_J2M_doubletagSR_GJets, h_J2M_doubletagSB_GJets, h_J2M_antitagSR_GJets, h_J2M_antitagSB_GJets};
  vector<TH1F*> histos_Rpf_J2_mJBins_GJets = {h_J2M_mjBins_doubletagSR_GJets, h_J2M_mjBins_doubletagSB_GJets, h_J2M_mjBins_antitagSR_GJets, h_J2M_mjBins_antitagSB_GJets};
  vector<TH1F*> histos_Rpfsingle_J1_GJets = {h_J1M_tagSR_GJets, h_J1M_tagSB_GJets, h_J1M_antitagSR_GJets, h_J1M_antitagSB_GJets};
  vector<TH1F*> histos_Rpfsingle_J2_GJets = {h_J2M_tagSR_GJets, h_J2M_tagSB_GJets, h_J2M_antitagSR_GJets, h_J2M_antitagSB_GJets};

  vector<TH1F*> histos_Rpf_J1_SnglT = {h_J1M_doubletagSR_SnglT, h_J1M_doubletagSB_SnglT, h_J1M_antitagSR_SnglT, h_J1M_antitagSB_SnglT};
  vector<TH1F*> histos_Rpf_J2_SnglT = {h_J2M_doubletagSR_SnglT, h_J2M_doubletagSB_SnglT, h_J2M_antitagSR_SnglT, h_J2M_antitagSB_SnglT};
  vector<TH1F*> histos_Rpf_J2_mJBins_SnglT = {h_J2M_mjBins_doubletagSR_SnglT, h_J2M_mjBins_doubletagSB_SnglT, h_J2M_mjBins_antitagSR_SnglT, h_J2M_mjBins_antitagSB_SnglT};
  vector<TH1F*> histos_Rpfsingle_J1_SnglT = {h_J1M_tagSR_SnglT, h_J1M_tagSB_SnglT, h_J1M_antitagSR_SnglT, h_J1M_antitagSB_SnglT};
  vector<TH1F*> histos_Rpfsingle_J2_SnglT = {h_J2M_tagSR_SnglT, h_J2M_tagSB_SnglT, h_J2M_antitagSR_SnglT, h_J2M_antitagSB_SnglT};

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
  vector<TH1F*> histos_j1mass_SnglT    = {h_J1M_doubletagSR_SnglT,h_J1M_doubletagSB_SnglT,h_J1M_tagSR_SnglT,h_J1M_tagSB_SnglT, h_J1M_antitagSR_SnglT, h_J1M_antitagSB_SnglT};
  vector<TH1F*> histos_j1mass_TT    = {h_J1M_doubletagSR_TT,h_J1M_doubletagSB_TT,h_J1M_tagSR_TT,h_J1M_tagSB_TT, h_J1M_antitagSR_TT, h_J1M_antitagSB_TT};
  vector<TH1F*> histos_j1mass_TT_Di    = {h_J1M_doubletagSR_TT_Di,h_J1M_doubletagSB_TT_Di,h_J1M_tagSR_TT_Di,h_J1M_tagSB_TT_Di, h_J1M_antitagSR_TT_Di, h_J1M_antitagSB_TT_Di};
  vector<TH1F*> histos_j1mass_TT_SL    = {h_J1M_doubletagSR_TT_SL,h_J1M_doubletagSB_TT_SL,h_J1M_tagSR_TT_SL,h_J1M_tagSB_TT_SL, h_J1M_antitagSR_TT_SL, h_J1M_antitagSB_TT_SL};
  vector<TH1F*> histos_j1mass_QCD   = {h_J1M_doubletagSR_QCD,h_J1M_doubletagSB_QCD,h_J1M_tagSR_QCD,h_J1M_tagSB_QCD, h_J1M_antitagSR_QCD, h_J1M_antitagSB_QCD};
  vector<TH1F*> histos_j1mass_T5HH1600 = {h_J1M_doubletagSR_T5HH1600,h_J1M_doubletagSB_T5HH1600,h_J1M_tagSR_T5HH1600,h_J1M_tagSB_T5HH1600, h_J1M_antitagSR_T5HH1600, h_J1M_antitagSB_T5HH1600};
  vector<TH1F*> histos_j1mass_T5HH2000 = {h_J1M_doubletagSR_T5HH2000,h_J1M_doubletagSB_T5HH2000,h_J1M_tagSR_T5HH2000,h_J1M_tagSB_T5HH2000, h_J1M_antitagSR_T5HH2000, h_J1M_antitagSB_T5HH2000};
  vector<TH1F*> histos_j1mass_T5HH2200 = {h_J1M_doubletagSR_T5HH2200,h_J1M_doubletagSB_T5HH2200,h_J1M_tagSR_T5HH2200,h_J1M_tagSB_T5HH2200, h_J1M_antitagSR_T5HH2200, h_J1M_antitagSB_T5HH2200};

  vector<TH1F*> histos_j2mass_ZJets = {h_J2M_doubletagSR_ZJets,h_J2M_doubletagSB_ZJets,h_J2M_tagSR_ZJets,h_J2M_tagSB_ZJets, h_J2M_antitagSR_ZJets, h_J2M_antitagSB_ZJets};
  vector<TH1F*> histos_j2mass_WJets = {h_J2M_doubletagSR_WJets,h_J2M_doubletagSB_WJets,h_J2M_tagSR_WJets,h_J2M_tagSB_WJets, h_J2M_antitagSR_WJets, h_J2M_antitagSB_WJets};
  vector<TH1F*> histos_j2mass_SnglT    = {h_J2M_doubletagSR_SnglT,h_J2M_doubletagSB_SnglT,h_J2M_tagSR_SnglT,h_J2M_tagSB_SnglT, h_J2M_antitagSR_SnglT, h_J2M_antitagSB_SnglT};
  vector<TH1F*> histos_j2mass_TT    = {h_J2M_doubletagSR_TT,h_J2M_doubletagSB_TT,h_J2M_tagSR_TT,h_J2M_tagSB_TT, h_J2M_antitagSR_TT, h_J2M_antitagSB_TT};
  vector<TH1F*> histos_j2mass_TT_Di    = {h_J2M_doubletagSR_TT_Di,h_J2M_doubletagSB_TT_Di,h_J2M_tagSR_TT_Di,h_J2M_tagSB_TT_Di, h_J2M_antitagSR_TT_Di, h_J2M_antitagSB_TT_Di};
  vector<TH1F*> histos_j2mass_TT_SL    = {h_J2M_doubletagSR_TT_SL,h_J2M_doubletagSB_TT_SL,h_J2M_tagSR_TT_SL,h_J2M_tagSB_TT_SL, h_J2M_antitagSR_TT_SL, h_J2M_antitagSB_TT_SL};
  vector<TH1F*> histos_j2mass_QCD   = {h_J2M_doubletagSR_QCD,h_J2M_doubletagSB_QCD,h_J2M_tagSR_QCD,h_J2M_tagSB_QCD, h_J2M_antitagSR_QCD, h_J2M_antitagSB_QCD};
  vector<TH1F*> histos_j2mass_T5HH1600 = {h_J2M_doubletagSR_T5HH1600,h_J2M_doubletagSB_T5HH1600,h_J2M_tagSR_T5HH1600,h_J2M_tagSB_T5HH1600, h_J2M_antitagSR_T5HH1600, h_J2M_antitagSB_T5HH1600};
  vector<TH1F*> histos_j2mass_T5HH2000 = {h_J2M_doubletagSR_T5HH2000,h_J2M_doubletagSB_T5HH2000,h_J2M_tagSR_T5HH2000,h_J2M_tagSB_T5HH2000, h_J2M_antitagSR_T5HH2000, h_J2M_antitagSB_T5HH2000};
  vector<TH1F*> histos_j2mass_T5HH2200 = {h_J2M_doubletagSR_T5HH2200,h_J2M_doubletagSB_T5HH2200,h_J2M_tagSR_T5HH2200,h_J2M_tagSB_T5HH2200, h_J2M_antitagSR_T5HH2200, h_J2M_antitagSB_T5HH2200};

  vector<TH1F*> histos_j1pt_ZJets = {h_J1Pt_doubletagSR_ZJets,h_J1Pt_doubletagSB_ZJets,h_J1Pt_tagSR_ZJets,h_J1Pt_tagSB_ZJets, h_J1Pt_antitagSR_ZJets, h_J1Pt_antitagSB_ZJets};
  vector<TH1F*> histos_j1pt_WJets = {h_J1Pt_doubletagSR_WJets,h_J1Pt_doubletagSB_WJets,h_J1Pt_tagSR_WJets,h_J1Pt_tagSB_WJets, h_J1Pt_antitagSR_WJets, h_J1Pt_antitagSB_WJets};
  vector<TH1F*> histos_j1pt_SnglT    = {h_J1Pt_doubletagSR_SnglT,h_J1Pt_doubletagSB_SnglT,h_J1Pt_tagSR_SnglT,h_J1Pt_tagSB_SnglT, h_J1Pt_antitagSR_SnglT, h_J1Pt_antitagSB_SnglT};
  vector<TH1F*> histos_j1pt_TT    = {h_J1Pt_doubletagSR_TT,h_J1Pt_doubletagSB_TT,h_J1Pt_tagSR_TT,h_J1Pt_tagSB_TT, h_J1Pt_antitagSR_TT, h_J1Pt_antitagSB_TT};
  vector<TH1F*> histos_j1pt_TT_Di    = {h_J1Pt_doubletagSR_TT_Di,h_J1Pt_doubletagSB_TT_Di,h_J1Pt_tagSR_TT_Di,h_J1Pt_tagSB_TT_Di, h_J1Pt_antitagSR_TT_Di, h_J1Pt_antitagSB_TT_Di};
  vector<TH1F*> histos_j1pt_TT_SL    = {h_J1Pt_doubletagSR_TT_SL,h_J1Pt_doubletagSB_TT_SL,h_J1Pt_tagSR_TT_SL,h_J1Pt_tagSB_TT_SL, h_J1Pt_antitagSR_TT_SL, h_J1Pt_antitagSB_TT_SL};
  vector<TH1F*> histos_j1pt_QCD   = {h_J1Pt_doubletagSR_QCD,h_J1Pt_doubletagSB_QCD,h_J1Pt_tagSR_QCD,h_J1Pt_tagSB_QCD, h_J1Pt_antitagSR_QCD, h_J1Pt_antitagSB_QCD};
  vector<TH1F*> histos_j1pt_T5HH1600 = {h_J1Pt_doubletagSR_T5HH1600,h_J1Pt_doubletagSB_T5HH1600,h_J1Pt_tagSR_T5HH1600,h_J1Pt_tagSB_T5HH1600, h_J1Pt_antitagSR_T5HH1600, h_J1Pt_antitagSB_T5HH1600};
  vector<TH1F*> histos_j1pt_T5HH2000 = {h_J1Pt_doubletagSR_T5HH2000,h_J1Pt_doubletagSB_T5HH2000,h_J1Pt_tagSR_T5HH2000,h_J1Pt_tagSB_T5HH2000, h_J1Pt_antitagSR_T5HH2000, h_J1Pt_antitagSB_T5HH2000};
  vector<TH1F*> histos_j1pt_T5HH2200 = {h_J1Pt_doubletagSR_T5HH2200,h_J1Pt_doubletagSB_T5HH2200,h_J1Pt_tagSR_T5HH2200,h_J1Pt_tagSB_T5HH2200, h_J1Pt_antitagSR_T5HH2200, h_J1Pt_antitagSB_T5HH2200};

  vector<TH1F*> histos_j2pt_ZJets = {h_J2Pt_doubletagSR_ZJets,h_J2Pt_doubletagSB_ZJets,h_J2Pt_tagSR_ZJets,h_J2Pt_tagSB_ZJets, h_J2Pt_antitagSR_ZJets, h_J2Pt_antitagSB_ZJets};
  vector<TH1F*> histos_j2pt_WJets = {h_J2Pt_doubletagSR_WJets,h_J2Pt_doubletagSB_WJets,h_J2Pt_tagSR_WJets,h_J2Pt_tagSB_WJets, h_J2Pt_antitagSR_WJets, h_J2Pt_antitagSB_WJets};
  vector<TH1F*> histos_j2pt_SnglT    = {h_J2Pt_doubletagSR_SnglT,h_J2Pt_doubletagSB_SnglT,h_J2Pt_tagSR_SnglT,h_J2Pt_tagSB_SnglT, h_J2Pt_antitagSR_SnglT, h_J2Pt_antitagSB_SnglT};
  vector<TH1F*> histos_j2pt_TT    = {h_J2Pt_doubletagSR_TT,h_J2Pt_doubletagSB_TT,h_J2Pt_tagSR_TT,h_J2Pt_tagSB_TT, h_J2Pt_antitagSR_TT, h_J2Pt_antitagSB_TT};
  vector<TH1F*> histos_j2pt_TT_Di    = {h_J2Pt_doubletagSR_TT_Di,h_J2Pt_doubletagSB_TT_Di,h_J2Pt_tagSR_TT_Di,h_J2Pt_tagSB_TT_Di, h_J2Pt_antitagSR_TT_Di, h_J2Pt_antitagSB_TT_Di};
  vector<TH1F*> histos_j2pt_TT_SL    = {h_J2Pt_doubletagSR_TT_SL,h_J2Pt_doubletagSB_TT_SL,h_J2Pt_tagSR_TT_SL,h_J2Pt_tagSB_TT_SL, h_J2Pt_antitagSR_TT_SL, h_J2Pt_antitagSB_TT_SL};
  vector<TH1F*> histos_j2pt_QCD   = {h_J2Pt_doubletagSR_QCD,h_J2Pt_doubletagSB_QCD,h_J2Pt_tagSR_QCD,h_J2Pt_tagSB_QCD, h_J2Pt_antitagSR_QCD, h_J2Pt_antitagSB_QCD};
  vector<TH1F*> histos_j2pt_T5HH1600 = {h_J2Pt_doubletagSR_T5HH1600,h_J2Pt_doubletagSB_T5HH1600,h_J2Pt_tagSR_T5HH1600,h_J2Pt_tagSB_T5HH1600, h_J2Pt_antitagSR_T5HH1600, h_J2Pt_antitagSB_T5HH1600};
  vector<TH1F*> histos_j2pt_T5HH2000 = {h_J2Pt_doubletagSR_T5HH2000,h_J2Pt_doubletagSB_T5HH2000,h_J2Pt_tagSR_T5HH2000,h_J2Pt_tagSB_T5HH2000, h_J2Pt_antitagSR_T5HH2000, h_J2Pt_antitagSB_T5HH2000};
  vector<TH1F*> histos_j2pt_T5HH2200 = {h_J2Pt_doubletagSR_T5HH2200,h_J2Pt_doubletagSB_T5HH2200,h_J2Pt_tagSR_T5HH2200,h_J2Pt_tagSB_T5HH2200, h_J2Pt_antitagSR_T5HH2200, h_J2Pt_antitagSB_T5HH2200};


  // vector<TH1F*> h_METShape_bkgSum_1H = {h_A1_sum,h_B1_sum,h_C_sum,h_COpt4_sum,h_D_sum,h_DOpt4_sum};
  // vector<TH1F*> h_METShape_bkgSum_1H = {h_A1_sum,h_B1_sum,h_C_sum,h_COpt2_sum,h_D_sum,h_DOpt2_sum};
  // vector<TH1F*> h_METShape_bkgSum_2H = {h_A_sum,h_B_sum,h_C_sum,h_COpt1_sum,h_D_sum,h_DOpt1_sum};
  vector<TH1F*> h_METShape_bkgSum_1H = {h_A1_sum,h_B1_sum,h_C_sum,h_COpt3_sum,h_D_sum,h_DOpt3_sum};
  vector<TH1F*> h_METShape_bkgSum_2H = {h_A_sum,h_B_sum,h_C_sum,h_COpt3_sum,h_D_sum,h_DOpt3_sum};

  vector<TH1F*> h_METShapePhoton_bkgSum_1H = {hP_A1_sum,hP_B1_sum,hP_C_sum,hP_COpt3_sum,hP_D_sum,hP_DOpt3_sum};
  vector<TH1F*> h_METShapePhoton_bkgSum_2H = {hP_A_sum,hP_B_sum,hP_C_sum,hP_COpt3_sum,hP_D_sum,hP_DOpt3_sum};


  // vector<TH1F*> h_METShape5_bkgSum_2H = {h_A5_sum,h_B5_sum,h_C5_sum,h_C5Opt1_sum,h_D5_sum,h_D5Opt1_sum};
  // vector<TH1F*> h_METShape5_bkgSum_1H = {h_A15_sum,h_B15_sum,h_C5_sum,h_C5Opt2_sum,h_D5_sum,h_D5Opt2_sum};
  // vector<TH1F*> h_METShape5_bkgSum_1H = {h_A15_sum,h_B15_sum,h_C5_sum,h_C5Opt4_sum,h_D5_sum,h_D5Opt4_sum};
  vector<TH1F*> h_METShape5_bkgSum_2H = {h_A5_sum,h_B5_sum,h_C5_sum,h_C5Opt3_sum,h_D5_sum,h_D5Opt3_sum};
  vector<TH1F*> h_METShape5_bkgSum_1H = {h_A15_sum,h_B15_sum,h_C5_sum,h_C5Opt3_sum,h_D5_sum,h_D5Opt3_sum};


  vector<TH1F*> h_METShape_SnglT_1H = {h_A1_SnglT,h_B1_SnglT,h_C_SnglT,h_COpt3_SnglT,h_D_SnglT,h_DOpt3_SnglT};
  vector<TH1F*> h_METShape_SnglT_2H = {h_A_SnglT,h_B_SnglT,h_C_SnglT,h_COpt3_SnglT,h_D_SnglT,h_DOpt3_SnglT};

  vector<TH1F*> h_METShape_TT_1H = {h_A1_TT,h_B1_TT,h_C_TT,h_COpt3_TT,h_D_TT,h_DOpt3_TT};
  vector<TH1F*> h_METShape_TT_2H = {h_A_TT,h_B_TT,h_C_TT,h_COpt3_TT,h_D_TT,h_DOpt3_TT};

  vector<TH1F*> h_METShape_ZJets_1H = {h_A1_ZJets,h_B1_ZJets,h_C_ZJets,h_COpt3_ZJets,h_D_ZJets,h_DOpt3_ZJets};
  vector<TH1F*> h_METShape_ZJets_2H = {h_A_ZJets,h_B_ZJets,h_C_ZJets,h_COpt3_ZJets,h_D_ZJets,h_DOpt3_ZJets};

  vector<TH1F*> h_METShape_WJets_1H = {h_A1_WJets,h_B1_WJets,h_C_WJets,h_COpt3_WJets,h_D_WJets,h_DOpt3_WJets};
  vector<TH1F*> h_METShape_WJets_2H = {h_A_WJets,h_B_WJets,h_C_WJets,h_COpt3_WJets,h_D_WJets,h_DOpt3_WJets};

  vector<TH1F*> h_METShape_GJets_1H = {h_A1_GJets,h_B1_GJets,h_C_GJets,h_COpt3_GJets,h_D_GJets,h_DOpt3_GJets};
  vector<TH1F*> h_METShape_GJets_2H = {h_A_GJets,h_B_GJets,h_C_GJets,h_COpt3_GJets,h_D_GJets,h_DOpt3_GJets};

  vector<TH1F*> h_METShape_QCD_1H = {h_A1_QCD,h_B1_QCD,h_C_QCD,h_COpt3_QCD,h_D_QCD,h_DOpt3_QCD};
  vector<TH1F*> h_METShape_QCD_2H = {h_A_QCD,h_B_QCD,h_C_QCD,h_COpt3_QCD,h_D_QCD,h_DOpt3_QCD};

  if (runMassCorrelations){
    std::cout<<"Running mass correlations..."<<std::endl;
    massCorrelations(histos_allRegions_sum,"BkgSum");
  }

  if (runPies){
    std::cout<<"Running pies..."<<std::endl;

    if (whichRegion=="signal"){
      pieChart({h_A_QCD,h_B_QCD}, {h_A_WJets,h_B_WJets}, {h_A_ZJets,h_B_ZJets}, {h_A_TT,h_B_TT},{h_A_SnglT,h_B_SnglT}, "2H", "all");
      pieChart({h_A_QCD}, {h_A_WJets}, {h_A_ZJets}, {h_A_TT},{h_A_SnglT}, "2HSR", "all");
      pieChart({h_B_QCD}, {h_B_WJets}, {h_B_ZJets}, {h_B_TT},{h_B_SnglT}, "2HSB", "all");
      pieChart({h_A1_QCD,h_B1_QCD}, {h_A1_WJets,h_B1_WJets}, {h_A1_ZJets,h_B1_ZJets}, {h_A1_TT,h_B1_TT}, {h_A1_SnglT,h_B1_SnglT}, "1H", "all");
      pieChart({h_A1_QCD}, {h_A1_WJets}, {h_A1_ZJets}, {h_A1_TT},{h_A1_SnglT}, "1HSR", "all");
      pieChart({h_B1_QCD}, {h_B1_WJets}, {h_B1_ZJets}, {h_B1_TT},{h_B1_SnglT}, "1HSB", "all");
      pieChart({h_C_QCD,h_D_QCD}, {h_C_WJets,h_D_WJets}, {h_C_ZJets,h_D_ZJets}, {h_C_TT,h_D_TT}, {h_C_SnglT,h_D_SnglT}, "0H", "all");
      pieChart({h_C_QCD}, {h_C_WJets}, {h_C_ZJets}, {h_C_TT}, {h_C_SnglT}, "0HSR", "all");
      pieChart({h_D_QCD}, {h_D_WJets}, {h_D_ZJets}, {h_D_TT}, {h_D_SnglT}, "0HSB", "all");
      pieChart({h_COpt1_QCD,h_DOpt1_QCD}, {h_COpt1_WJets,h_DOpt1_WJets}, {h_COpt1_ZJets,h_DOpt1_ZJets}, {h_COpt1_TT,h_DOpt1_TT}, {h_COpt1_SnglT,h_DOpt1_SnglT}, "0H,BTagsM>1", "all");
      pieChart({h_COpt2_QCD,h_DOpt2_QCD}, {h_COpt2_WJets,h_DOpt2_WJets}, {h_COpt2_ZJets,h_DOpt2_ZJets}, {h_COpt2_TT,h_DOpt2_TT}, {h_COpt2_SnglT,h_DOpt2_SnglT}, "0H,BTagsM>0", "all");
      pieChart({h_COpt3_QCD,h_DOpt3_QCD}, {h_COpt3_WJets,h_DOpt3_WJets}, {h_COpt3_ZJets,h_DOpt3_ZJets}, {h_COpt3_TT,h_DOpt3_TT}, {h_COpt3_SnglT,h_DOpt3_SnglT}, "0H,BTagsT>0", "all");
      pieChart({h_COpt4_QCD,h_DOpt4_QCD}, {h_COpt4_WJets,h_DOpt4_WJets}, {h_COpt4_ZJets,h_DOpt4_ZJets}, {h_COpt4_TT,h_DOpt4_TT}, {h_COpt4_SnglT,h_DOpt4_SnglT}, "0H,BTagsL>1 and BTagsM>0", "all");

      pieChart({h_A_QCD,h_B_QCD}, {h_A_WJets,h_B_WJets}, {h_A_ZJets,h_B_ZJets}, {h_A_TT,h_B_TT},{h_A_SnglT,h_B_SnglT}, "2H", "bin1");
      pieChart({h_A_QCD}, {h_A_WJets}, {h_A_ZJets}, {h_A_TT},{h_A_SnglT}, "2HSR", "bin1");
      pieChart({h_B_QCD}, {h_B_WJets}, {h_B_ZJets}, {h_B_TT},{h_B_SnglT}, "2HSB", "bin1");
      pieChart({h_A1_QCD,h_B1_QCD}, {h_A1_WJets,h_B1_WJets}, {h_A1_ZJets,h_B1_ZJets}, {h_A1_TT,h_B1_TT}, {h_A1_SnglT,h_B1_SnglT}, "1H", "bin1");
      pieChart({h_A1_QCD}, {h_A1_WJets}, {h_A1_ZJets}, {h_A1_TT},{h_A1_SnglT}, "1HSR", "bin1");
      pieChart({h_B1_QCD}, {h_B1_WJets}, {h_B1_ZJets}, {h_B1_TT},{h_B1_SnglT}, "1HSB", "bin1");
      pieChart({h_C_QCD,h_D_QCD}, {h_C_WJets,h_D_WJets}, {h_C_ZJets,h_D_ZJets}, {h_C_TT,h_D_TT}, {h_C_SnglT,h_D_SnglT}, "0H", "bin1");
      pieChart({h_C_QCD}, {h_C_WJets}, {h_C_ZJets}, {h_C_TT}, {h_C_SnglT}, "0HSR", "bin1");
      pieChart({h_D_QCD}, {h_D_WJets}, {h_D_ZJets}, {h_D_TT}, {h_D_SnglT}, "0HSB", "bin1");
      pieChart({h_COpt1_QCD,h_DOpt1_QCD}, {h_COpt1_WJets,h_DOpt1_WJets}, {h_COpt1_ZJets,h_DOpt1_ZJets}, {h_COpt1_TT,h_DOpt1_TT}, {h_COpt1_SnglT,h_DOpt1_SnglT}, "0H,BTagsM>1", "bin1");
      pieChart({h_COpt2_QCD,h_DOpt2_QCD}, {h_COpt2_WJets,h_DOpt2_WJets}, {h_COpt2_ZJets,h_DOpt2_ZJets}, {h_COpt2_TT,h_DOpt2_TT}, {h_COpt2_SnglT,h_DOpt2_SnglT}, "0H,BTagsM>0", "bin1");
      pieChart({h_COpt3_QCD,h_DOpt3_QCD}, {h_COpt3_WJets,h_DOpt3_WJets}, {h_COpt3_ZJets,h_DOpt3_ZJets}, {h_COpt3_TT,h_DOpt3_TT}, {h_COpt3_SnglT,h_DOpt3_SnglT}, "0H,BTagsT>0", "bin1");
      pieChart({h_COpt4_QCD,h_DOpt4_QCD}, {h_COpt4_WJets,h_DOpt4_WJets}, {h_COpt4_ZJets,h_DOpt4_ZJets}, {h_COpt4_TT,h_DOpt4_TT}, {h_COpt4_SnglT,h_DOpt4_SnglT}, "0H,BTagsL>1 and BTagsM>0", "bin1");

      pieChart({h_A_QCD,h_B_QCD}, {h_A_WJets,h_B_WJets}, {h_A_ZJets,h_B_ZJets}, {h_A_TT,h_B_TT},{h_A_SnglT,h_B_SnglT}, "2H", "bin2");
      pieChart({h_A_QCD}, {h_A_WJets}, {h_A_ZJets}, {h_A_TT},{h_A_SnglT}, "2HSR", "bin2");
      pieChart({h_B_QCD}, {h_B_WJets}, {h_B_ZJets}, {h_B_TT},{h_B_SnglT}, "2HSB", "bin2");
      pieChart({h_A1_QCD,h_B1_QCD}, {h_A1_WJets,h_B1_WJets}, {h_A1_ZJets,h_B1_ZJets}, {h_A1_TT,h_B1_TT}, {h_A1_SnglT,h_B1_SnglT}, "1H", "bin2");
      pieChart({h_A1_QCD}, {h_A1_WJets}, {h_A1_ZJets}, {h_A1_TT},{h_A1_SnglT}, "1HSR", "bin2");
      pieChart({h_B1_QCD}, {h_B1_WJets}, {h_B1_ZJets}, {h_B1_TT},{h_B1_SnglT}, "1HSB", "bin2");
      pieChart({h_C_QCD,h_D_QCD}, {h_C_WJets,h_D_WJets}, {h_C_ZJets,h_D_ZJets}, {h_C_TT,h_D_TT}, {h_C_SnglT,h_D_SnglT}, "0H", "bin2");
      pieChart({h_C_QCD}, {h_C_WJets}, {h_C_ZJets}, {h_C_TT}, {h_C_SnglT}, "0HSR", "bin2");
      pieChart({h_D_QCD}, {h_D_WJets}, {h_D_ZJets}, {h_D_TT}, {h_D_SnglT}, "0HSB", "bin2");
      pieChart({h_COpt1_QCD,h_DOpt1_QCD}, {h_COpt1_WJets,h_DOpt1_WJets}, {h_COpt1_ZJets,h_DOpt1_ZJets}, {h_COpt1_TT,h_DOpt1_TT}, {h_COpt1_SnglT,h_DOpt1_SnglT}, "0H,BTagsM>1", "bin2");
      pieChart({h_COpt2_QCD,h_DOpt2_QCD}, {h_COpt2_WJets,h_DOpt2_WJets}, {h_COpt2_ZJets,h_DOpt2_ZJets}, {h_COpt2_TT,h_DOpt2_TT}, {h_COpt2_SnglT,h_DOpt2_SnglT}, "0H,BTagsM>0", "bin2");
      pieChart({h_COpt3_QCD,h_DOpt3_QCD}, {h_COpt3_WJets,h_DOpt3_WJets}, {h_COpt3_ZJets,h_DOpt3_ZJets}, {h_COpt3_TT,h_DOpt3_TT}, {h_COpt3_SnglT,h_DOpt3_SnglT}, "0H,BTagsT>0", "bin2");
      pieChart({h_COpt4_QCD,h_DOpt4_QCD}, {h_COpt4_WJets,h_DOpt4_WJets}, {h_COpt4_ZJets,h_DOpt4_ZJets}, {h_COpt4_TT,h_DOpt4_TT}, {h_COpt4_SnglT,h_DOpt4_SnglT}, "0H,BTagsL>1 and BTagsM>0", "bin2");

      pieChart({h_A_QCD,h_B_QCD}, {h_A_WJets,h_B_WJets}, {h_A_ZJets,h_B_ZJets}, {h_A_TT,h_B_TT},{h_A_SnglT,h_B_SnglT}, "2H", "bin3");
      pieChart({h_A_QCD}, {h_A_WJets}, {h_A_ZJets}, {h_A_TT},{h_A_SnglT}, "2HSR", "bin3");
      pieChart({h_B_QCD}, {h_B_WJets}, {h_B_ZJets}, {h_B_TT},{h_B_SnglT}, "2HSB", "bin3");
      pieChart({h_A1_QCD,h_B1_QCD}, {h_A1_WJets,h_B1_WJets}, {h_A1_ZJets,h_B1_ZJets}, {h_A1_TT,h_B1_TT}, {h_A1_SnglT,h_B1_SnglT}, "1H", "bin3");
      pieChart({h_A1_QCD}, {h_A1_WJets}, {h_A1_ZJets}, {h_A1_TT},{h_A1_SnglT}, "1HSR", "bin3");
      pieChart({h_B1_QCD}, {h_B1_WJets}, {h_B1_ZJets}, {h_B1_TT},{h_B1_SnglT}, "1HSB", "bin3");
      pieChart({h_C_QCD,h_D_QCD}, {h_C_WJets,h_D_WJets}, {h_C_ZJets,h_D_ZJets}, {h_C_TT,h_D_TT}, {h_C_SnglT,h_D_SnglT}, "0H", "bin3");
      pieChart({h_C_QCD}, {h_C_WJets}, {h_C_ZJets}, {h_C_TT}, {h_C_SnglT}, "0HSR", "bin3");
      pieChart({h_D_QCD}, {h_D_WJets}, {h_D_ZJets}, {h_D_TT}, {h_D_SnglT}, "0HSB", "bin3");
      pieChart({h_COpt1_QCD,h_DOpt1_QCD}, {h_COpt1_WJets,h_DOpt1_WJets}, {h_COpt1_ZJets,h_DOpt1_ZJets}, {h_COpt1_TT,h_DOpt1_TT}, {h_COpt1_SnglT,h_DOpt1_SnglT}, "0H,BTagsM>1", "bin3");
      pieChart({h_COpt2_QCD,h_DOpt2_QCD}, {h_COpt2_WJets,h_DOpt2_WJets}, {h_COpt2_ZJets,h_DOpt2_ZJets}, {h_COpt2_TT,h_DOpt2_TT}, {h_COpt2_SnglT,h_DOpt2_SnglT}, "0H,BTagsM>0", "bin3");
      pieChart({h_COpt3_QCD,h_DOpt3_QCD}, {h_COpt3_WJets,h_DOpt3_WJets}, {h_COpt3_ZJets,h_DOpt3_ZJets}, {h_COpt3_TT,h_DOpt3_TT}, {h_COpt3_SnglT,h_DOpt3_SnglT}, "0H,BTagsT>0", "bin3");
      pieChart({h_COpt4_QCD,h_DOpt4_QCD}, {h_COpt4_WJets,h_DOpt4_WJets}, {h_COpt4_ZJets,h_DOpt4_ZJets}, {h_COpt4_TT,h_DOpt4_TT}, {h_COpt4_SnglT,h_DOpt4_SnglT}, "0H,BTagsL>1 and BTagsM>0", "bin3");

      // pieChart({h_COpt5_QCD,h_DOpt5_QCD}, {h_COpt5_WJets,h_DOpt5_WJets}, {h_COpt5_ZJets,h_DOpt5_ZJets}, {h_COpt5_TT,h_DOpt5_TT}, {h_COpt5_SnglT,h_DOpt5_SnglT}, "0H,BTagsT>0 and BTagsM>1");

    }
    else if (whichRegion=="singleLept"){
      pieChart1l({h_A_WJets,h_B_WJets}, {h_A_TT,h_B_TT}, {h_A_SnglT,h_B_SnglT}, "2H", "all");
      pieChart1l({h_A_WJets}, {h_A_TT}, {h_A_SnglT}, "2HSR", "all");
      pieChart1l({h_B_WJets},{h_B_TT}, {h_B_SnglT}, "2HSB", "all");
      pieChart1l({h_A1_WJets,h_B1_WJets}, {h_A1_TT,h_B1_TT}, {h_A1_SnglT,h_B1_SnglT}, "1H", "all");
      pieChart1l({h_A1_WJets}, {h_A1_TT}, {h_A1_SnglT}, "1HSR", "all");
      pieChart1l({h_B1_WJets}, {h_B1_TT}, {h_B1_SnglT}, "1HSB", "all");
      pieChart1l({h_C_WJets,h_D_WJets}, {h_C_TT,h_D_TT}, {h_C_SnglT,h_C_SnglT}, "0H", "all");
      pieChart1l({h_C_WJets}, {h_C_TT}, {h_C_SnglT}, "0HSR", "all");
      pieChart1l({h_D_WJets}, {h_D_TT}, {h_D_SnglT}, "0HSB", "all");
      pieChart1l({h_COpt1_WJets,h_DOpt1_WJets}, {h_COpt1_TT,h_DOpt1_TT}, {h_COpt1_SnglT,h_DOpt1_SnglT}, "0H,BTagsM>1", "all");
      pieChart1l({h_COpt3_WJets,h_DOpt3_WJets}, {h_COpt3_TT,h_DOpt3_TT}, {h_COpt3_SnglT,h_DOpt3_SnglT}, "0H,BTagsT>0", "all");
      pieChart1l({h_COpt4_WJets,h_DOpt4_WJets}, {h_COpt4_TT,h_DOpt4_TT}, {h_COpt4_SnglT,h_DOpt4_SnglT}, "0H,BTagsL>1 and BTagsM>0", "all");
      // pieChart1l({h_COpt2_WJets,h_DOpt2_WJets}, {h_COpt2_TT,h_DOpt2_TT}, {h_COpt2_SnglT,h_DOpt2_SnglT}, "0H,BTagsM>0");

      pieChart1l({h_A_WJets,h_B_WJets}, {h_A_TT,h_B_TT}, {h_A_SnglT,h_B_SnglT}, "2H", "bin1");
      pieChart1l({h_A_WJets}, {h_A_TT}, {h_A_SnglT}, "2HSR", "bin1");
      pieChart1l({h_B_WJets},{h_B_TT}, {h_B_SnglT}, "2HSB", "bin1");
      pieChart1l({h_A1_WJets,h_B1_WJets}, {h_A1_TT,h_B1_TT}, {h_A1_SnglT,h_B1_SnglT}, "1H", "bin1");
      pieChart1l({h_A1_WJets}, {h_A1_TT}, {h_A1_SnglT}, "1HSR", "bin1");
      pieChart1l({h_B1_WJets}, {h_B1_TT}, {h_B1_SnglT}, "1HSB", "bin1");
      pieChart1l({h_C_WJets,h_D_WJets}, {h_C_TT,h_D_TT}, {h_C_SnglT,h_C_SnglT}, "0H", "bin1");
      pieChart1l({h_C_WJets}, {h_C_TT}, {h_C_SnglT}, "0HSR", "bin1");
      pieChart1l({h_D_WJets}, {h_D_TT}, {h_D_SnglT}, "0HSB", "bin1");
      pieChart1l({h_COpt1_WJets,h_DOpt1_WJets}, {h_COpt1_TT,h_DOpt1_TT}, {h_COpt1_SnglT,h_DOpt1_SnglT}, "0H,BTagsM>1", "bin1");
      pieChart1l({h_COpt3_WJets,h_DOpt3_WJets}, {h_COpt3_TT,h_DOpt3_TT}, {h_COpt3_SnglT,h_DOpt3_SnglT}, "0H,BTagsT>0", "bin1");
      pieChart1l({h_COpt4_WJets,h_DOpt4_WJets}, {h_COpt4_TT,h_DOpt4_TT}, {h_COpt4_SnglT,h_DOpt4_SnglT}, "0H,BTagsL>1 and BTagsM>0", "bin1");
      // pieChart1l({h_COpt2_WJets,h_DOpt2_WJets}, {h_COpt2_TT,h_DOpt2_TT}, {h_COpt2_SnglT,h_DOpt2_SnglT}, "0H,BTagsM>0");

      pieChart1l({h_A_WJets,h_B_WJets}, {h_A_TT,h_B_TT}, {h_A_SnglT,h_B_SnglT}, "2H", "bin2");
      pieChart1l({h_A_WJets}, {h_A_TT}, {h_A_SnglT}, "2HSR", "bin2");
      pieChart1l({h_B_WJets},{h_B_TT}, {h_B_SnglT}, "2HSB", "bin2");
      pieChart1l({h_A1_WJets,h_B1_WJets}, {h_A1_TT,h_B1_TT}, {h_A1_SnglT,h_B1_SnglT}, "1H", "bin2");
      pieChart1l({h_A1_WJets}, {h_A1_TT}, {h_A1_SnglT}, "1HSR", "bin2");
      pieChart1l({h_B1_WJets}, {h_B1_TT}, {h_B1_SnglT}, "1HSB", "bin2");
      pieChart1l({h_C_WJets,h_D_WJets}, {h_C_TT,h_D_TT}, {h_C_SnglT,h_C_SnglT}, "0H", "bin2");
      pieChart1l({h_C_WJets}, {h_C_TT}, {h_C_SnglT}, "0HSR", "bin2");
      pieChart1l({h_D_WJets}, {h_D_TT}, {h_D_SnglT}, "0HSB", "bin2");
      pieChart1l({h_COpt1_WJets,h_DOpt1_WJets}, {h_COpt1_TT,h_DOpt1_TT}, {h_COpt1_SnglT,h_DOpt1_SnglT}, "0H,BTagsM>1", "bin2");
      pieChart1l({h_COpt3_WJets,h_DOpt3_WJets}, {h_COpt3_TT,h_DOpt3_TT}, {h_COpt3_SnglT,h_DOpt3_SnglT}, "0H,BTagsT>0", "bin2");
      pieChart1l({h_COpt4_WJets,h_DOpt4_WJets}, {h_COpt4_TT,h_DOpt4_TT}, {h_COpt4_SnglT,h_DOpt4_SnglT}, "0H,BTagsL>1 and BTagsM>0", "bin2");
      // pieChart1l({h_COpt2_WJets,h_DOpt2_WJets}, {h_COpt2_TT,h_DOpt2_TT}, {h_COpt2_SnglT,h_DOpt2_SnglT}, "0H,BTagsM>0");

      pieChart1l({h_A_WJets,h_B_WJets}, {h_A_TT,h_B_TT}, {h_A_SnglT,h_B_SnglT}, "2H", "bin3");
      pieChart1l({h_A_WJets}, {h_A_TT}, {h_A_SnglT}, "2HSR", "bin3");
      pieChart1l({h_B_WJets},{h_B_TT}, {h_B_SnglT}, "2HSB", "bin3");
      pieChart1l({h_A1_WJets,h_B1_WJets}, {h_A1_TT,h_B1_TT}, {h_A1_SnglT,h_B1_SnglT}, "1H", "bin3");
      pieChart1l({h_A1_WJets}, {h_A1_TT}, {h_A1_SnglT}, "1HSR", "bin3");
      pieChart1l({h_B1_WJets}, {h_B1_TT}, {h_B1_SnglT}, "1HSB", "bin3");
      pieChart1l({h_C_WJets,h_D_WJets}, {h_C_TT,h_D_TT}, {h_C_SnglT,h_C_SnglT}, "0H", "bin3");
      pieChart1l({h_C_WJets}, {h_C_TT}, {h_C_SnglT}, "0HSR", "bin3");
      pieChart1l({h_D_WJets}, {h_D_TT}, {h_D_SnglT}, "0HSB", "bin3");
      pieChart1l({h_COpt1_WJets,h_DOpt1_WJets}, {h_COpt1_TT,h_DOpt1_TT}, {h_COpt1_SnglT,h_DOpt1_SnglT}, "0H,BTagsM>1", "bin3");
      pieChart1l({h_COpt3_WJets,h_DOpt3_WJets}, {h_COpt3_TT,h_DOpt3_TT}, {h_COpt3_SnglT,h_DOpt3_SnglT}, "0H,BTagsT>0", "bin3");
      pieChart1l({h_COpt4_WJets,h_DOpt4_WJets}, {h_COpt4_TT,h_DOpt4_TT}, {h_COpt4_SnglT,h_DOpt4_SnglT}, "0H,BTagsL>1 and BTagsM>0", "bin3");
      // pieChart1l({h_COpt2_WJets,h_DOpt2_WJets}, {h_COpt2_TT,h_DOpt2_TT}, {h_COpt2_SnglT,h_DOpt2_SnglT}, "0H,BTagsM>0");
    }

    else if (whichRegion=="photon"){ //GJets, QCD
      pieChartPhoton({h_A_GJets,h_B_GJets}, {h_A_QCD,h_B_QCD}, "2H", "all");
      pieChartPhoton({h_A_GJets}, {h_A_QCD}, "2HSR", "all");
      pieChartPhoton({h_B_GJets},{h_B_QCD},"2HSB", "all");
      pieChartPhoton({h_A1_GJets,h_B1_GJets}, {h_A1_QCD,h_B1_QCD}, "1H", "all");
      pieChartPhoton({h_A1_GJets}, {h_A1_QCD}, "1HSR", "all");
      pieChartPhoton({h_B1_GJets}, {h_B1_QCD}, "1HSB", "all");
      pieChartPhoton({h_C_GJets,h_D_GJets}, {h_C_QCD,h_D_QCD}, "0H", "all");
      pieChartPhoton({h_C_GJets}, {h_C_QCD}, "0HSR", "all");
      pieChartPhoton({h_D_GJets}, {h_D_QCD}, "0HSB", "all");
      pieChartPhoton({h_COpt1_GJets,h_DOpt1_GJets}, {h_COpt1_QCD,h_DOpt1_QCD}, "0H,BTagsM>1", "all");
      pieChartPhoton({h_COpt3_GJets,h_DOpt3_GJets}, {h_COpt3_QCD,h_DOpt3_QCD}, "0H,BTagsT>0", "all");
      pieChartPhoton({h_COpt4_GJets,h_DOpt4_GJets}, {h_COpt4_QCD,h_DOpt4_QCD}, "0H,BTagsL>1 and BTagsM>0", "all");

      pieChartPhoton({h_A_GJets,h_B_GJets}, {h_A_QCD,h_B_QCD}, "2H", "bin1");
      pieChartPhoton({h_A_GJets}, {h_A_QCD}, "2HSR", "bin1");
      pieChartPhoton({h_B_GJets},{h_B_QCD}, "2HSB", "bin1");
      pieChartPhoton({h_A1_GJets,h_B1_GJets}, {h_A1_QCD,h_B1_QCD}, "1H", "bin1");
      pieChartPhoton({h_A1_GJets}, {h_A1_QCD}, "1HSR", "bin1");
      pieChartPhoton({h_B1_GJets}, {h_B1_QCD}, "1HSB", "bin1");
      pieChartPhoton({h_C_GJets,h_D_GJets}, {h_C_QCD,h_D_QCD}, "0H", "bin1");
      pieChartPhoton({h_C_GJets}, {h_C_QCD}, "0HSR", "bin1");
      pieChartPhoton({h_D_GJets}, {h_D_QCD}, "0HSB", "bin1");
      pieChartPhoton({h_COpt1_GJets,h_DOpt1_GJets}, {h_COpt1_QCD,h_DOpt1_QCD}, "0H,BTagsM>1", "bin1");
      pieChartPhoton({h_COpt3_GJets,h_DOpt3_GJets}, {h_COpt3_QCD,h_DOpt3_QCD}, "0H,BTagsT>0", "bin1");
      pieChartPhoton({h_COpt4_GJets,h_DOpt4_GJets}, {h_COpt4_QCD,h_DOpt4_QCD}, "0H,BTagsL>1 and BTagsM>0", "bin1");

      pieChartPhoton({h_A_GJets,h_B_GJets}, {h_A_QCD,h_B_QCD}, "2H", "bin2");
      pieChartPhoton({h_A_GJets}, {h_A_QCD}, "2HSR", "bin2");
      pieChartPhoton({h_B_GJets},{h_B_QCD}, "2HSB", "bin2");
      pieChartPhoton({h_A1_GJets,h_B1_GJets}, {h_A1_QCD,h_B1_QCD}, "1H", "bin2");
      pieChartPhoton({h_A1_GJets}, {h_A1_QCD}, "1HSR", "bin2");
      pieChartPhoton({h_B1_GJets}, {h_B1_QCD}, "1HSB", "bin2");
      pieChartPhoton({h_C_GJets,h_D_GJets}, {h_C_QCD,h_D_QCD}, "0H", "bin2");
      pieChartPhoton({h_C_GJets}, {h_C_QCD}, "0HSR", "bin2");
      pieChartPhoton({h_D_GJets}, {h_D_QCD}, "0HSB", "bin2");
      pieChartPhoton({h_COpt1_GJets,h_DOpt1_GJets}, {h_COpt1_QCD,h_DOpt1_QCD}, "0H,BTagsM>1", "bin2");
      pieChartPhoton({h_COpt3_GJets,h_DOpt3_GJets}, {h_COpt3_QCD,h_DOpt3_QCD}, "0H,BTagsT>0", "bin2");
      pieChartPhoton({h_COpt4_GJets,h_DOpt4_GJets}, {h_COpt4_QCD,h_DOpt4_QCD}, "0H,BTagsL>1 and BTagsM>0", "bin2");

      pieChartPhoton({h_A_GJets,h_B_GJets}, {h_A_QCD,h_B_QCD}, "2H", "bin3");
      pieChartPhoton({h_A_GJets}, {h_A_QCD},  "2HSR", "bin3");
      pieChartPhoton({h_B_GJets},{h_B_QCD},  "2HSB", "bin3");
      pieChartPhoton({h_A1_GJets,h_B1_GJets}, {h_A1_QCD,h_B1_QCD}, "1H", "bin3");
      pieChartPhoton({h_A1_GJets}, {h_A1_QCD}, "1HSR", "bin3");
      pieChartPhoton({h_B1_GJets}, {h_B1_QCD}, "1HSB", "bin3");
      pieChartPhoton({h_C_GJets,h_D_GJets}, {h_C_QCD,h_D_QCD}, "0H", "bin3");
      pieChartPhoton({h_C_GJets}, {h_C_QCD}, "0HSR", "bin3");
      pieChartPhoton({h_D_GJets}, {h_D_QCD}, "0HSB", "bin3");
      pieChartPhoton({h_COpt1_GJets,h_DOpt1_GJets}, {h_COpt1_QCD,h_DOpt1_QCD}, "0H,BTagsM>1", "bin3");
      pieChartPhoton({h_COpt3_GJets,h_DOpt3_GJets}, {h_COpt3_QCD,h_DOpt3_QCD}, "0H,BTagsT>0", "bin3");
      pieChartPhoton({h_COpt4_GJets,h_DOpt4_GJets}, {h_COpt4_QCD,h_DOpt4_QCD}, "0H,BTagsL>1 and BTagsM>0", "bin3");
    }
  }


  if (runFullBkg) {
    std::cout<<"Running full bkg closure..."<<std::endl;

    makeFullBkgClosure(histos_ABCD_sum, "BkgSum", "Double",year);
    makeFullBkgClosure(histos_A1B1CD_sum, "BkgSum", "Single",year);
    if (runPure) {
      makePureBkgClosure(histos_ABCD_sum, "BkgSum", "Double",year);
      makePureBkgClosure(histos_A1B1CD_sum, "BkgSum", "Single",year);
    }
  }

  if (runABCDPlots){
    std::cout<<"Running ABCD plots..."<<std::endl;

    if (whichRegion=="signal"){
      makeABCDPlot(histos_ABCD_sum, "BkgSum", "Double",year);
      makeABCDPlot(histos_A1B1CD_sum, "BkgSum", "Single",year);

      makeABCDPlot(histos_ABCDlowPU_sum, "BkgSumlowPU", "Double",year);
      makeABCDPlot(histos_A1B1CDlowPU_sum, "BkgSumlowPU", "Single",year);

      makeABCDPlot(histos_ABCDhighPU_sum, "BkgSumhighPU", "Double",year);
      makeABCDPlot(histos_A1B1CDhighPU_sum, "BkgSumhighPU", "Single",year);



      makeABCDPlot(histos_ABCD_QCD, "QCD", "Double",year);
      makeABCDPlot(histos_ABCD_SnglT, "SnglT", "Double",year);
      makeABCDPlot(histos_ABCD_TT, "TT", "Double",year);
      makeABCDPlot(histos_ABCD_WJets, "WJets", "Double",year);
      makeABCDPlot(histos_ABCD_ZJets, "ZJets", "Double",year);

      makeABCDPlot(histos_A1B1CD_QCD, "QCD", "Single",year);
      makeABCDPlot(histos_A1B1CD_SnglT, "SnglT", "Single",year);
      makeABCDPlot(histos_A1B1CD_TT, "TT", "Single",year);
      makeABCDPlot(histos_A1B1CD_WJets, "WJets", "Single",year);
      makeABCDPlot(histos_A1B1CD_ZJets, "ZJets", "Single",year);
    }

    else if (whichRegion=="singleLept"){
      makeABCDPlot(histos_ABCD_sum, "BkgSum", "Double",year);
      makeABCDPlot(histos_A1B1CD_sum, "BkgSum", "Single",year);

      // makeABCDPlot(histos_ABCD_data, "Data", "Double",year);
      // makeABCDPlot(histos_A1B1CD_data, "Data", "Single",year);

      makeABCDPlot(histos_ABCD_SnglT, "SnglT", "Double",year);
      makeABCDPlot(histos_ABCD_TT, "TT", "Double",year);
      makeABCDPlot(histos_ABCD_WJets, "WJets", "Double",year);

      makeABCDPlot(histos_A1B1CD_SnglT, "SnglT", "Single",year);
      makeABCDPlot(histos_A1B1CD_TT, "TT", "Single",year);
      makeABCDPlot(histos_A1B1CD_WJets, "WJets", "Single",year);
    }
    else if (whichRegion=="photon"){
      makeABCDPlot(histosPhoton_ABCD_sum, "BkgSum", "Double",year);
      makeABCDPlot(histosPhoton_A1B1CD_sum, "BkgSum", "Single",year);

      // makeABCDPlot(histos_ABCD_sum, "BkgSum", "Double",year);
      // makeABCDPlot(histos_A1B1CD_sum, "BkgSum", "Single",year);
      // makeABCDPlot(histos_ABCD_QCD, "QCD", "Double",year);
      // makeABCDPlot(histos_ABCD_GJets, "GJets", "Double",year);
      // makeABCDPlot(histos_A1B1CD_QCD, "QCD", "Single",year);
      // makeABCDPlot(histos_A1B1CD_GJets, "GJets", "Single",year);
    }
  }

  if (runRPFPlots){
    std::cout<<"Running RPF plots..."<<std::endl;

    if (whichRegion=="signal"){
      //DM
      if (runDM) {
        makeRpfPlot(histos_RPF_DM, "BkgSum_DM", "Avg", "Double",year);
        makeRpfPlot(histos_RPFsingle_DM, "BkgSum_DM", "Avg", "Single",year);
      }


      //J1+J2
      makeRpfPlot(histos_Rpf_J1J2_sum, "BkgSum", "Both", "Double",year);
      makeRpfPlot(histos_Rpfsingle_J1J2_sum, "BkgSum", "Both", "Single",year);


      //J1 Double tag region
      makeRpfPlot(histos_Rpf_J1_sum, "BkgSum", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_QCD, "QCD", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_SnglT, "SnglT", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_TT, "TT", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_WJets, "WJets", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_ZJets, "ZJets", "J1", "Double",year);

      //J1 Single tag region
      makeRpfPlot(histos_Rpfsingle_J1_sum, "BkgSum", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_QCD, "QCD", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_SnglT, "SnglT", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_TT, "TT", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_WJets, "WJets", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_ZJets, "ZJets", "J1", "Single",year);

      //J2 double tag region
      makeRpfPlot(histos_Rpf_J2_sum, "BkgSum", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_QCD, "QCD", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_SnglT, "SnglT", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_TT, "TT", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_WJets, "WJets", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_ZJets, "ZJets", "J2", "Double",year);

      //J2 double tag region, with three mass bins
      makeRpfPlot(histos_Rpf_J2_mJBins_sum, "BkgSum_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_QCD, "QCD_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_SnglT, "SnglT_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_TT, "TT_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_WJets, "WJets_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_ZJets, "ZJets_mJ", "J2", "Double",year);


      //J2 single tag region
      makeRpfPlot(histos_Rpfsingle_J2_sum, "BkgSum", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_QCD, "QCD", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_SnglT, "SnglT", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_TT, "TT", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_WJets, "WJets", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_ZJets, "ZJets", "J2", "Single",year);
    }
    if (whichRegion=="singleLept"){
      //DM
      if (runDM) {
        makeRpfPlot(histos_RPF_DM, "BkgSum_DM", "Avg", "Double",year);
        makeRpfPlot(histos_RPFsingle_DM, "BkgSum_DM", "Avg", "Single",year);
      }

      //J1+J2
      makeRpfPlot(histos_Rpf_J1J2_sum, "BkgSum", "Both", "Double",year);
      makeRpfPlot(histos_Rpfsingle_J1J2_sum, "BkgSum", "Both", "Single",year);


      //J1 Double tag region
      makeRpfPlot(histos_Rpf_J1_sum, "BkgSum", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_SnglT, "SnglT", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_TT, "TT", "J1", "Double",year);
      makeRpfPlot(histos_Rpf_J1_WJets, "WJets", "J1", "Double",year);

      //J1 Single tag region
      makeRpfPlot(histos_Rpfsingle_J1_sum, "BkgSum", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_SnglT, "SnglT", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_TT, "TT", "J1", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J1_WJets, "WJets", "J1", "Single",year);

      //J2 double tag region
      makeRpfPlot(histos_Rpf_J2_sum, "BkgSum", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_SnglT, "SnglT", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_TT, "TT", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_WJets, "WJets", "J2", "Double",year);

      //J2 single tag region
      makeRpfPlot(histos_Rpfsingle_J2_sum, "BkgSum", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_SnglT, "SnglT", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_TT, "TT", "J2", "Single",year);
      makeRpfPlot(histos_Rpfsingle_J2_WJets, "WJets", "J2", "Single",year);

      //J2 double tag region, with three mass bins
      makeRpfPlot(histos_Rpf_J2_mJBins_sum, "BkgSum_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_SnglT, "SnglT_mJ", "J2", "Double",year);
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

      //J2 double tag region, with three mass bins
      makeRpfPlot(histos_Rpf_J2_mJBins_sum, "BkgSum_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_QCD, "QCD_mJ", "J2", "Double",year);
      makeRpfPlot(histos_Rpf_J2_mJBins_GJets, "GJets_mJ", "J2", "Double",year);
    }
  }


  if (runStacks && whichRegion=="signal") {
    std::cout<<"Running stacks..."<<std::endl;

    makeStackPlot(h_baseline_j1mass_QCD,h_baseline_j1mass_TT,h_baseline_j1mass_WJets,h_baseline_j1mass_ZJets,h_baseline_j1mass_SnglT,h_baseline_j1mass_T5HH1600,h_baseline_j1mass_T5HH2000,h_baseline_j1mass_T5HH2200,h_baseline_j1mass_TChiHH600,h_baseline_j1mass_TChiHH800,h_baseline_j1mass_TChiHH1000,"leadmass");
    makeStackPlot(h_baseline_j2mass_QCD,h_baseline_j2mass_TT,h_baseline_j2mass_WJets,h_baseline_j2mass_ZJets,h_baseline_j2mass_SnglT,h_baseline_j2mass_T5HH1600,h_baseline_j2mass_T5HH2000,h_baseline_j2mass_T5HH2200,h_baseline_j2mass_TChiHH600,h_baseline_j2mass_TChiHH800,h_baseline_j2mass_TChiHH1000,"subleadmass");
    makeStackPlot(h_baseline_j1pt_QCD,h_baseline_j1pt_TT,h_baseline_j1pt_WJets,h_baseline_j1pt_ZJets,h_baseline_j1pt_SnglT,h_baseline_j1pt_T5HH1600,h_baseline_j1pt_T5HH2000,h_baseline_j1pt_T5HH2200,h_baseline_j1pt_TChiHH600,h_baseline_j1pt_TChiHH800,h_baseline_j1pt_TChiHH1000,"leadpt");
    makeStackPlot(h_baseline_j2pt_QCD,h_baseline_j2pt_TT,h_baseline_j2pt_WJets,h_baseline_j2pt_ZJets,h_baseline_j2pt_SnglT,h_baseline_j2pt_T5HH1600,h_baseline_j2pt_T5HH2000,h_baseline_j2pt_T5HH2200,h_baseline_j2pt_TChiHH600,h_baseline_j2pt_TChiHH800,h_baseline_j2pt_TChiHH1000,"subleadpt");

    // makeMETStack(h_baseline_QCD, h_baseline_TT,h_baseline_WJets,h_baseline_ZJets, h_baseline_SnglT, h_baseline_T5HH1600,h_baseline_T5HH2000,h_baseline_T5HH2200,h_baseline_TChiHH600,h_baseline_TChiHH800,h_baseline_TChiHH1000,"baseline");
    // makeMETStack(h_2H_QCD, h_2H_TT,h_2H_WJets,h_2H_ZJets, h_2H_SnglT, h_2H_T5HH1600,h_2H_T5HH2000,h_2H_T5HH2200,h_2H_TChiHH600,h_2H_TChiHH800,h_2H_TChiHH1000,"2H");
    // makeMETStack(h_1H_QCD, h_1H_TT,h_1H_WJets,h_1H_ZJets, h_1H_SnglT, h_1H_T5HH1600,h_1H_T5HH2000,h_1H_T5HH2200,h_1H_TChiHH600,h_1H_TChiHH800,h_1H_TChiHH1000,"1H");
    // makeMETStack(h_0H_QCD, h_0H_TT,h_0H_WJets,h_0H_ZJets, h_0H_SnglT, h_0H_T5HH1600,h_0H_T5HH2000,h_0H_T5HH2200,h_0H_TChiHH600,h_0H_TChiHH800,h_0H_TChiHH1000,"0H");
    // makeMETStack(h_0Hb_QCD, h_0Hb_TT,h_0Hb_WJets,h_0Hb_ZJets, h_0Hb_SnglT, h_0Hb_T5HH1600,h_0Hb_T5HH2000,h_0Hb_T5HH2200,h_0Hb_TChiHH600,h_0Hb_TChiHH800,h_0Hb_TChiHH1000,"0Hb");
  }

  if (runStacks && whichRegion=="singleLept") {
    std::cout<<"Running stacks..."<<std::endl;

    makeSingleLeptStackPlot(histos_j1mass_TT,histos_j1mass_WJets,histos_j1mass_SnglT,"leadmass");
    makeSingleLeptStackPlot(histos_j2mass_TT,histos_j2mass_WJets,histos_j2mass_SnglT,"subleadmass");
    makeSingleLeptStackPlot(histos_j1pt_TT,histos_j1pt_WJets,histos_j1pt_SnglT,"leadpt");
    makeSingleLeptStackPlot(histos_j2pt_TT,histos_j2pt_WJets,histos_j2pt_SnglT,"subleadpt");
  }



  if (runMETNorm){
    std::cout<<"Running MET shape comparison..."<<std::endl;

    if (whichRegion=="signal"){
      makeMETNorm(h_METShape_bkgSum_2H,"BkgSum","Double");
      makeMETNorm(h_METShape_TT_2H,"TT","Double");
      makeMETNorm(h_METShape_WJets_2H,"WJets","Double");
      makeMETNorm(h_METShape_ZJets_2H,"ZJets","Double");
      makeMETNorm(h_METShape_SnglT_2H,"QCD","Double");
      makeMETNorm(h_METShape_QCD_2H,"QCD","Double");

      makeMETNorm(h_METShape_bkgSum_1H, "BkgSum","Single");
      makeMETNorm(h_METShape_TT_1H, "TT","Single");
      makeMETNorm(h_METShape_WJets_1H, "WJets","Single");
      makeMETNorm(h_METShape_ZJets_1H, "ZJets","Single");
      makeMETNorm(h_METShape_SnglT_1H, "QCD","Single");
      makeMETNorm(h_METShape_QCD_1H, "QCD","Single");

      //Just the SR and the cuts region
      makeMETNormCompare(h_METShape5_bkgSum_2H,"BkgSum5","Double");

      makeMETNormCompare(h_METShape_bkgSum_2H,"BkgSum","Double");
      makeMETNormCompare(h_METShape_TT_2H,"TT","Double");
      makeMETNormCompare(h_METShape_WJets_2H,"WJets","Double");
      makeMETNormCompare(h_METShape_ZJets_2H,"ZJets","Double");
      makeMETNormCompare(h_METShape_SnglT_2H,"QCD","Double");
      makeMETNormCompare(h_METShape_QCD_2H,"QCD","Double");

      makeMETNormCompare(h_METShape5_bkgSum_1H, "BkgSum5","Single");

      makeMETNormCompare(h_METShape_bkgSum_1H, "BkgSum","Single");
      makeMETNormCompare(h_METShape_TT_1H, "TT","Single");
      makeMETNormCompare(h_METShape_WJets_1H, "WJets","Single");
      makeMETNormCompare(h_METShape_ZJets_1H, "ZJets","Single");
      makeMETNormCompare(h_METShape_SnglT_1H, "SnglT","Single");
      makeMETNormCompare(h_METShape_QCD_1H, "QCD","Single");
    }
    else if (whichRegion=="singleLept"){
      makeMETNorm(h_METShape_bkgSum_2H,"BkgSum","Double");
      makeMETNorm(h_METShape_SnglT_2H,"SnglT","Double");
      makeMETNorm(h_METShape_TT_2H,"TT","Double");
      makeMETNorm(h_METShape_WJets_2H,"WJets","Double");

      makeMETNorm(h_METShape_bkgSum_1H, "BkgSum","Single");
      makeMETNorm(h_METShape_SnglT_1H, "SnglT","Single");
      makeMETNorm(h_METShape_TT_1H, "TT","Single");
      makeMETNorm(h_METShape_WJets_1H, "WJets","Single");

      //Just the SR and the cuts region
      makeMETNormCompare(h_METShape5_bkgSum_2H,"BkgSum5","Double");

      makeMETNormCompare(h_METShape_bkgSum_2H,"BkgSum","Double");
      makeMETNormCompare(h_METShape_SnglT_2H,"SnglT","Double");
      makeMETNormCompare(h_METShape_TT_2H,"TT","Double");
      makeMETNormCompare(h_METShape_WJets_2H,"WJets","Double");

      makeMETNormCompare(h_METShape5_bkgSum_1H, "BkgSum5","Single");

      makeMETNormCompare(h_METShape_bkgSum_1H, "BkgSum","Single");
      makeMETNormCompare(h_METShape_SnglT_1H, "SnglT","Single");
      makeMETNormCompare(h_METShape_TT_1H, "TT","Single");
      makeMETNormCompare(h_METShape_WJets_1H, "WJets","Single");
    }
    else if (whichRegion=="photon"){
      // makeMETNorm(h_METShape_bkgSum_2H,"BkgSum","Double");
      // makeMETNorm(h_METShape_QCD_2H,"QCD","Double");
      // makeMETNorm(h_METShape_GJets_2H,"GJets","Double");
      //
      // makeMETNorm(h_METShape_bkgSum_1H, "BkgSum","Single");
      // makeMETNorm(h_METShape_QCD_1H, "QCD","Single");
      // makeMETNorm(h_METShape_GJets_1H, "GJets","Single");

      //Just the SR and the cuts region
      // makeMETNormCompare(h_METShape5_bkgSum_2H,"BkgSum5","Double");

      makeMETNormCompare(h_METShapePhoton_bkgSum_2H,"BkgSum","Double");
      // makeMETNormCompare(h_METShape_QCD_2H,"QCD","Double");
      // makeMETNormCompare(h_METShape_GJets_2H,"GJets","Double");

      // makeMETNormCompare(h_METShape5_bkgSum_1H, "BkgSum5","Single");

      makeMETNormCompare(h_METShapePhoton_bkgSum_1H, "BkgSum","Single");
      // makeMETNormCompare(h_METShape_QCD_1H, "QCD","Single");
      // makeMETNormCompare(h_METShape_GJets_1H, "GJets","Single");
    }
  }

  if (runTableOfYields){
    std::cout<<"Running table of yields..."<<std::endl;

    if (whichRegion=="signal"){
      if (runPure) {
        tableOfYieldsPure(histos_ABCD_sum, "BkgSum", "Double");
        tableOfYieldsPure(histos_A1B1CD_sum, "BkgSum", "Single");
      }
      tableOfYields(histos_ABCD_sum, "BkgSum", "Double");
      tableOfYields(histos_A1B1CD_sum, "BkgSum", "Single");

      tableOfYields(histos_ABCDlowPU_sum, "BkgSumlowPU", "Double");
      tableOfYields(histos_A1B1CDlowPU_sum, "BkgSumlowPU", "Single");

      tableOfYields(histos_ABCDhighPU_sum, "BkgSumhighPU", "Double");
      tableOfYields(histos_A1B1CDhighPU_sum, "BkgSumhighPU", "Single");

      tableOfYields(histos_ABCD_QCD, "QCD", "Double");
      tableOfYields(histos_ABCD_SnglT, "SnglT", "Double");
      tableOfYields(histos_ABCD_TT, "TT", "Double");
      tableOfYields(histos_ABCD_WJets, "WJets", "Double");
      tableOfYields(histos_ABCD_ZJets, "ZJets", "Double");

      tableOfYields(histos_A1B1CD_QCD, "QCD", "Single");
      tableOfYields(histos_A1B1CD_SnglT, "SnglT", "Single");
      tableOfYields(histos_A1B1CD_TT, "TT", "Single");
      tableOfYields(histos_A1B1CD_WJets, "WJets", "Single");
      tableOfYields(histos_A1B1CD_ZJets, "ZJets", "Single");
    }

    else if (whichRegion=="singleLept"){
      if (runPure) {
        tableOfYieldsPure(histos_ABCD_sum, "BkgSum", "Double");
        tableOfYieldsPure(histos_A1B1CD_sum, "BkgSum", "Single");
      }

      tableOfYields(histos_ABCD_sum, "BkgSum", "Double");
      tableOfYields(histos_A1B1CD_sum, "BkgSum", "Single");
      tableOfYields(histos_ABCD_SnglT, "SnglT", "Double");
      tableOfYields(histos_ABCD_TT, "TT", "Double");
      tableOfYields(histos_ABCD_WJets, "WJets", "Double");
      tableOfYields(histos_A1B1CD_SnglT, "SnglT", "Single");
      tableOfYields(histos_A1B1CD_TT, "TT", "Single");
      tableOfYields(histos_A1B1CD_WJets, "WJets", "Single");
    }
    else if (whichRegion=="photon"){
      tableOfYields(histos_ABCD_sum, "BkgSum", "Double");
      tableOfYields(histos_A1B1CD_sum, "BkgSum", "Single");

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

  DrawOverflow(histo_A); DrawOverflow(histo_B); DrawOverflow(histo_C); DrawOverflow(histo_D);

  bool isData = false;
  TString histoName = dem_histos[0]->GetName();
  if (histoName.Contains("data")) isData=true;


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
  if (whichRegion!="photon") graph->GetXaxis()->SetRangeUser(300,1400);
  else graph->GetXaxis()->SetRangeUser(100,1400);
  graph->GetYaxis()->SetRangeUser(-1.0,5.0);
  if (tagType=="Single") graph->GetYaxis()->SetRangeUser(0.0,2.0);
  graph->Draw("APE");
  can_h->Modified(); can_h->Update();
  if (whichRegion!="photon") graph->GetXaxis()->SetRangeUser(300,1400);
  else graph->GetXaxis()->SetRangeUser(100,1400);
  can_h->Modified(); can_h->Update();

  if (whichRegion!="photon") {
    TLine *line = new TLine(300,1.0,1400,1.0);
    line->SetLineColor(kRed); line->SetLineStyle(2);
    line->Draw("same");
  }
  else {
    TLine *line = new TLine(100,1.0,1400,1.0);
    line->SetLineColor(kRed); line->SetLineStyle(2);
    line->Draw("same");
  }


  TString savename = "ABCD_"+tagType+"_"+bkgType+"_"+year;
  cdABCD->cd();
  can_h->Write(savename);
  if (savePDFs) can_h->SaveAs("ANPlots/"+savename+".pdf","PDF");
  delete can_h;
}

void makeFullBkgClosure(vector<TH1F*> dem_histos, TString bkgType = "", TString tagType = "", TString year = "") {
  TH1F *histo_A = (TH1F*)dem_histos[0]->Clone("histo_A"); TH1F *histo_B = (TH1F*)dem_histos[1]->Clone("histo_B");
  TH1F *histo_C = (TH1F*)dem_histos[2]->Clone("histo_C"); TH1F *histo_D = (TH1F*)dem_histos[3]->Clone("histo_D");

  DrawOverflow(histo_A); DrawOverflow(histo_B); DrawOverflow(histo_C); DrawOverflow(histo_D);

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

  //V18
  if (whichRegion=="signal"){
    if (tagType=="Double"){
      kappa = 1.08; kappa_error = 0.20;
      bkgFrac1 = 0.851;   bkgFrac1_error = 0.016;
      bkgFrac2 = 0.119;  bkgFrac2_error = 0.005;
      bkgFrac3 = 0.031; bkgFrac3_error = 0.003;
    }
    else if (tagType=="Single"){
      kappa = 1.03; kappa_error = 0.04;
      bkgFrac1 = 0.851;   bkgFrac1_error = 0.016;
      bkgFrac2 = 0.119;  bkgFrac2_error = 0.005;
      bkgFrac3 = 0.031; bkgFrac3_error = 0.003;
    }
  }

  if (whichRegion=="singleLept"){
    if (tagType=="Double"){
      kappa = 0.79; kappa_error = 0.13;
      bkgFrac1 = 0.871;  bkgFrac1_error = 0.016;
      bkgFrac2 = 0.109; bkgFrac2_error = 0.007;
      bkgFrac3 = 0.020; bkgFrac3_error = 0.002;
    }
    else if (tagType=="Single"){
      kappa = 0.94; kappa_error = 0.05;
      bkgFrac1 = 0.871;  bkgFrac1_error = 0.016;
      bkgFrac2 = 0.109; bkgFrac2_error = 0.007;
      bkgFrac3 = 0.020; bkgFrac3_error = 0.002;
    }
  }


  //V17
  // if (whichRegion=="signal"){
  //   if (tagType=="Double"){
  //     kappa = 1.1; kappa_error = 0.18;
  //     bkgFrac1 = 0.857;   bkgFrac1_error = 0.023;
  //     bkgFrac2 = 0.118;  bkgFrac2_error = 0.0079;
  //     bkgFrac3 = 0.0251; bkgFrac3_error = 0.0033;
  //   }
  //   else if (tagType=="Single"){
  //     kappa = 0.95; kappa_error = 0.05;
  //     bkgFrac1 = 0.821;   bkgFrac1_error = 0.0111;
  //     bkgFrac2 = 0.136;  bkgFrac2_error = 0.0040;
  //     bkgFrac3 = 0.044; bkgFrac3_error = 0.0020;
  //   }
  // }
  //
  // if (whichRegion=="singleLept"){
  //   if (tagType=="Double"){
  //     kappa = 0.66; kappa_error = 0.137;
  //     bkgFrac1 = 0.923;  bkgFrac1_error = 0.0149;
  //     bkgFrac2 = 0.067; bkgFrac2_error = 0.0038;
  //     bkgFrac3 = 0.009; bkgFrac3_error = 0.0013;
  //   }
  //   else if (tagType=="Single"){
  //     kappa = 0.91; kappa_error = 0.0424;
  //     bkgFrac1 = 0.911;  bkgFrac1_error = 0.0089;
  //     bkgFrac2 = 0.077; bkgFrac2_error = 0.0024;
  //     bkgFrac3 = 0.012; bkgFrac3_error = 0.0008;
  //   }
  // }

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
  if (whichRegion!="photon") graph->GetXaxis()->SetRangeUser(300,1400);
  else graph->GetXaxis()->SetRangeUser(100,1400);
  can_h->Modified(); can_h->Update();

  TLine *line = new TLine(300,1.0,1400,1.0);
  line->SetLineColor(kRed); line->SetLineStyle(2);
  line->Draw("same");

  cdClose->cd();
  can_h->Write("FullBkg_"+tagType);
  if (savePDFs) can_h->SaveAs("ANPlots/FullBkg_"+tagType+".pdf","PDF");
  delete can_h;
}

void makePureBkgClosure(vector<TH1F*> dem_histos, TString bkgType = "", TString tagType = "", TString year = "") {
  TH1F *histo_A = (TH1F*)dem_histos[0]->Clone("histo_A"); TH1F *histo_B = (TH1F*)dem_histos[1]->Clone("histo_B");
  TH1F *histo_C = (TH1F*)dem_histos[2]->Clone("histo_C"); TH1F *histo_D = (TH1F*)dem_histos[3]->Clone("histo_D");

  DrawOverflow(histo_A); DrawOverflow(histo_B); DrawOverflow(histo_C); DrawOverflow(histo_D);

  TH1F *histo_pred = (TH1F*)histo_B->Clone("histo_pred");
  histo_pred->Multiply(histo_C);histo_pred->Divide(histo_D);

  DrawOverflow(histo_pred);

  TH1F *h_finalBkg = (TH1F*)histo_pred->Clone("h_finalBkg");

  float kappa1; float kappa_err1;
  float kappa2; float kappa_err2;
  float kappa3; float kappa_err3;


  //V18
  if (whichRegion=="signal"){
    if (tagType=="Double"){
      kappa1 = 1.17; kappa_err1 = 0.24;
      kappa2 = 0.46; kappa_err2 = 0.19;
      kappa3 = 0.84; kappa_err3 = 0.59;
    }
    else if (tagType=="Single"){
      kappa1 = 1.04; kappa_err1 = 0.04;
      kappa2 = 0.98; kappa_err2 = 0.10;
      kappa3 = 0.86; kappa_err3 = 0.14;
    }
  }

  if (whichRegion=="singleLept"){
    if (tagType=="Double"){
      kappa1 = 0.80; kappa_err1 = 0.11;
      kappa2 = 0.69; kappa_err2 = 0.28;
      kappa3 = 0.80; kappa_err3 = 0.85;
    }
    else if (tagType=="Single"){
      kappa1 = 0.95; kappa_err1 = 0.05;
      kappa2 = 0.89; kappa_err2 = 0.13;
      kappa3 = 0.76; kappa_err3 = 0.28;
    }
  }

  float bin1Content = histo_pred->GetBinContent(1)*kappa1;
  float bin2Content = histo_pred->GetBinContent(2)*kappa2;
  float bin3Content = histo_pred->GetBinContent(3)*kappa3;


  h_finalBkg->SetBinContent(1,bin1Content);
  h_finalBkg->SetBinContent(2,bin2Content);
  h_finalBkg->SetBinContent(3,bin3Content);

  float bin1Error = bin1Content*sqrt( TMath::Power(histo_pred->GetBinError(1)/histo_pred->GetBinContent(1),2) + TMath::Power(kappa_err1/kappa1,2) );
  float bin2Error = bin2Content*sqrt( TMath::Power(histo_pred->GetBinError(2)/histo_pred->GetBinContent(2),2) + TMath::Power(kappa_err2/kappa2,2) );
  float bin3Error = bin3Content*sqrt( TMath::Power(histo_pred->GetBinError(3)/histo_pred->GetBinContent(3),2) + TMath::Power(kappa_err3/kappa3,2) );


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
    legend->AddEntry(h_finalBkg,"Pred: B2*C/D*kappa","f");
  }
  else if (tagType=="Single") {
    legend->AddEntry(histo_A,"MC: A1","lp");
    legend->AddEntry(h_finalBkg,"Pred: B1*C/D*kappa","f");
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
  if (whichRegion!="photon") graph->GetXaxis()->SetRangeUser(300,1400);
  else graph->GetXaxis()->SetRangeUser(100,1400);
  graph->GetYaxis()->SetRangeUser(-1.0,5.0);
  if (tagType=="Single") graph->GetYaxis()->SetRangeUser(0.0,2.0);
  graph->Draw("APE");
  can_h->Modified(); can_h->Update();
  if (whichRegion!="photon") graph->GetXaxis()->SetRangeUser(300,1400);
  else graph->GetXaxis()->SetRangeUser(100,1400);
  can_h->Modified(); can_h->Update();

  TLine *line = new TLine(300,1.0,1400,1.0);
  line->SetLineColor(kRed); line->SetLineStyle(2);
  line->Draw("same");

  cdClose->cd();
  can_h->Write("PureABCDBkg_"+tagType);
  delete can_h;
}

void makeRpfPlot(vector<TH1F*> dem_histos, TString bkgType = "", TString jetType = "", TString tagType = "", TString year = "") {
  TH1F * num   = (TH1F*)dem_histos[1]->Clone("num"+tagType+"_"+jetType+"_"+bkgType); //this should be 2H SB
  TH1F * denom = (TH1F*)dem_histos[3]->Clone("denom"+tagType+"_"+jetType+"_"+bkgType);  //this should be 0H SB
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

  // graph->SetMinimum(-0.1); graph->GetXaxis()->SetRangeUser(50,250);
  graph->SetMinimum(-0.1); graph->GetXaxis()->SetRangeUser(60,260);
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
    legend->AddEntry(num, "1H SR+SB (numerator)", "PL");
    legend->AddEntry(denom, "0H SR+SB (denominator)", "PL");
    graph->SetTitle(";Soft-drop mass [GeV];R(1H/0H)");
  }
  else if (tagType == "Double"){
    legend->AddEntry(num, "2H SR+SB (numerator)", "PL");
    legend->AddEntry(denom, "0H SR+SB (denominator)", "PL");
    graph->SetTitle(";Soft-drop mass [GeV];R(2H/0H)");
  }

  legend->Draw();

  TLine *line2 = new TLine(95.0,0.0,95.0,1.0);
  line2->SetLineColor(kBlack);line2->SetLineStyle(2);
  line2->Draw("same");

  TLine *line3 = new TLine(145.0,0.0,145.0,1.0);
  line3->SetLineColor(kBlack);line3->SetLineStyle(2);
  line3->Draw("same");

  pad2->cd();
  graph->Draw("APE");
  graph->GetYaxis()->SetRangeUser(-0.1,0.5);
  if (bkgType.Contains("_DM")) graph->GetXaxis()->SetTitle("Avg. soft-drop mass [GeV]");
  can_h->Modified();
  can_h->Update();

  TF1*f0=new TF1("f0","pol0", 60,160);
  graph->Fit("f0","Q","R",60,160);
  double p0 = f0->GetParameter(0); //this does, indeed, work
  double error = f0->GetParError(0);
  graph->GetFunction("f0")->SetBit(TF1::kNotDraw);

  // TLine *line = new TLine(50.0,p0,250.0,p0);
  TLine *line = new TLine(60.0,p0,260.0,p0);
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
  line2->Draw("same");
  line3->Draw("same");

  cdRPFFinal->cd();
  can_h->Write(thisType);
  if (savePDFs) can_h->SaveAs("ANPlots/"+thisType+".pdf", "PDF");

  cdRPF->cd();
  num->Write("Num_"+thisType);
  denom->Write("Denom_"+thisType);
  graph->Write("RPFRatio_"+thisType);
  delete can_h;
}


void makeStackPlot(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets, TH1F* h_SnglT, TH1F* h_T5HH1600,TH1F* h_T5HH2000,TH1F* h_T5HH2200,
TH1F* h_TChiHH600,TH1F* h_TChiHH800,TH1F* h_TChiHH1000,TString which = ""){
  TH1F *h_QCD_stack = (TH1F*)h_QCD->Clone("QCDAll");
  // for (int i=1;i<h_QCD.size()) h_QCD_stack->Add(h_QCD[i]);
  h_QCD_stack->SetFillColor(kGray); h_QCD_stack->SetMarkerStyle(21); h_QCD_stack->SetMarkerColor(kGray);

  TH1F *h_TT_stack = (TH1F*)h_TT->Clone("TTAll");
  // for (int i=1;i<h_TT.size()) h_TT_stack->Add(h_TT[i]);
  h_TT_stack->SetFillColor(kCyan);h_TT_stack->SetMarkerStyle(21); h_TT_stack->SetMarkerColor(kCyan);

  TH1F *h_WJets_stack = (TH1F*)h_WJets->Clone("WJetsAll");
  // for (int i=1;i<h_WJets.size()) h_WJets_stack->Add(h_WJets[i]);
  h_WJets_stack->SetFillColor(kBlue);h_WJets_stack->SetMarkerStyle(21); h_WJets_stack->SetMarkerColor(kBlue);

  TH1F *h_ZJets_stack = (TH1F*)h_ZJets->Clone("WZetsAll");
  // for (int i=1;i<h_ZJets.size()) h_ZJets_stack->Add(h_ZJets[i]);
  h_ZJets_stack->SetFillColor(kGreen+2); h_ZJets_stack->SetMarkerStyle(21); h_ZJets_stack->SetMarkerColor(kGreen+2);

  TH1F *h_SnglT_stack = (TH1F*)h_SnglT->Clone("SnglTAll");
  // for (int i=1;i<h_SnglT.size()) h_T_stack->Add(h_SnglT[i]);
  h_SnglT_stack->SetFillColor(kCyan+2); h_SnglT_stack->SetMarkerStyle(21); h_SnglT_stack->SetMarkerColor(kCyan+2);

  TH1F *h_T5HH1600_stack = (TH1F*)h_T5HH1600->Clone("T5HH1600");
  // for (int i=1;i<h_T5HH1600.size()) h_T5HH1600_stack->Add(h_T5HH1600[i]);
  h_T5HH1600_stack->SetMarkerStyle(1); h_T5HH1600_stack->SetMarkerColor(kRed);
  h_T5HH1600_stack->SetLineColor(kRed); h_T5HH1600_stack->SetLineStyle(1); h_T5HH1600_stack->SetLineWidth(2);

  TH1F *h_T5HH2000_stack = (TH1F*)h_T5HH2000->Clone("T5HH2000");
  // for (int i=1;i<h_T5HH2000.size()) h_T5HH2000_stack->Add(h_T5HH2000[i]);
  h_T5HH2000_stack->SetMarkerStyle(1); h_T5HH2000_stack->SetMarkerColor(kBlack);
  h_T5HH2000_stack->SetLineColor(kBlack); h_T5HH2000_stack->SetLineStyle(1); h_T5HH2000_stack->SetLineWidth(2);

  TH1F *h_T5HH2200_stack = (TH1F*)h_T5HH2200->Clone("T5HH2200");
  // for (int i=1;i<h_T5HH2200.size()) h_T5HH2200_stack->Add(h_T5HH2200[i]);
  h_T5HH2200_stack->SetMarkerStyle(1); h_T5HH2200_stack->SetMarkerColor(kMagenta);
  h_T5HH2200_stack->SetLineColor(kMagenta); h_T5HH2200_stack->SetLineStyle(1); h_T5HH2200_stack->SetLineWidth(2);

  TH1F *h_TChiHH600_stack = (TH1F*)h_TChiHH600->Clone("TChiHH600");
  // for (int i=1;i<h_TChiHH600.size()) h_TChiHH600_stack->Add(h_TChiHH600[i]);
  h_TChiHH600_stack->SetMarkerStyle(1); h_TChiHH600_stack->SetMarkerColor(kRed);
  h_TChiHH600_stack->SetLineColor(kRed); h_TChiHH600_stack->SetLineStyle(2); h_TChiHH600_stack->SetLineWidth(2);

  TH1F *h_TChiHH800_stack = (TH1F*)h_TChiHH800->Clone("TChiHH800");
  // for (int i=1;i<h_TChiHH800.size()) h_TChiHH800_stack->Add(h_TChiHH800[i]);
  h_TChiHH800_stack->SetMarkerStyle(1); h_TChiHH800_stack->SetMarkerColor(kBlack);
  h_TChiHH800_stack->SetLineColor(kBlack); h_TChiHH800_stack->SetLineStyle(2); h_TChiHH800_stack->SetLineWidth(2);

  TH1F *h_TChiHH1000_stack = (TH1F*)h_TChiHH1000->Clone("TChiHH1000");
  // for (int i=1;i<h_TChiHH1000.size()) h_TChiHH1000_stack->Add(h_TChiHH1000[i]);
  h_TChiHH1000_stack->SetMarkerStyle(1); h_TChiHH1000_stack->SetMarkerColor(kMagenta);
  h_TChiHH1000_stack->SetLineColor(kMagenta); h_TChiHH1000_stack->SetLineStyle(2); h_TChiHH1000_stack->SetLineWidth(2);

  THStack * doubleBStack = new THStack("hs","");
  doubleBStack->Add(h_ZJets_stack);
  doubleBStack->Add(h_WJets_stack);
  doubleBStack->Add(h_TT_stack);
  doubleBStack->Add(h_SnglT_stack);
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

  if (which == "DoubleB") doubleBStack->GetXaxis()->SetTitle("Deep bb-tag");
  else if (which == "leadmass") doubleBStack->GetXaxis()->SetTitle("Lead jet soft drop mass [GeV]");
  else if (which == "subleadmass") doubleBStack->GetXaxis()->SetTitle("Sublead jet soft drop mass [GeV]");
  else if (which == "leadpt") doubleBStack->GetXaxis()->SetTitle("Lead jet p_{T} [GeV]");
  else if (which == "subleadpt") doubleBStack->GetXaxis()->SetTitle("Sublead jet p_{T} [GeV]");

  h_T5HH1600_stack->Draw("same l"); h_T5HH2000_stack->Draw("same l"); h_T5HH2200_stack->Draw("same l");
  h_TChiHH600_stack->Draw("same l"); h_TChiHH800_stack->Draw("same l"); h_TChiHH1000_stack->Draw("same l");

  TLegend* legend = new TLegend(0.35,0.8,0.9,0.9) ;
  legend->AddEntry(h_QCD_stack, "QCD", "f");
  legend->AddEntry(h_TT_stack, "TT", "f");
  legend->AddEntry(h_SnglT_stack, "SnglT", "f");
  legend->AddEntry(h_WJets_stack, "WJets", "f");
  legend->AddEntry(h_ZJets_stack, "ZJets", "f");
  legend->AddEntry(h_T5HH1600_stack, "T5HH1600", "l");
  legend->AddEntry(h_T5HH2000_stack, "T5HH2000", "l");
  legend->AddEntry(h_T5HH2200_stack, "T5HH2200", "l");
  legend->AddEntry(h_TChiHH600_stack, "TChiHH600", "l");
  legend->AddEntry(h_TChiHH800_stack, "TChiHH800", "l");
  legend->AddEntry(h_TChiHH1000_stack, "TChiHH1000", "l");
  legend->SetNColumns(4);
  legend->SetBorderSize(0);
  legend->Draw("same");
  can_h->Update();can_h->Modified();

  TString savename = which+"_Stack";
  cdOther->cd();
  can_h->Write(savename);
  if (savePDFs) can_h->SaveAs("ANPlots/"+savename+".pdf","PDF");
  delete can_h;
}

void makeSingleLeptStackPlot(vector<TH1F*> h_TT,vector<TH1F*> h_WJets,vector<TH1F*> h_SnglT, TString which = ""){
  TH1F * h_TT_stack = (TH1F*)h_TT[0]->Clone("TTAll");
  for (int i=1;i<h_TT.size();i++) h_TT_stack->Add(h_TT[i]);
  h_TT_stack->SetFillColor(kCyan); h_TT_stack->SetMarkerStyle(21); h_TT_stack->SetMarkerColor(kCyan);

  TH1F *h_WJets_stack = (TH1F*)h_WJets[0]->Clone("WJetsAll");
  for (int i=1;i<h_WJets.size();i++) h_WJets_stack->Add(h_WJets[i]);
  h_WJets_stack->SetFillColor(kBlue); h_WJets_stack->SetMarkerStyle(21); h_WJets_stack->SetMarkerColor(kBlue);

  TH1F * h_SnglT_stack = (TH1F*)h_SnglT[0]->Clone("SnglTAll");
  for (int i=1;i<h_SnglT.size();i++) h_SnglT_stack->Add(h_SnglT[i]);
  h_SnglT_stack->SetFillColor(kCyan+2); h_SnglT_stack->SetMarkerStyle(21); h_SnglT_stack->SetMarkerColor(kCyan+2);


  THStack * doubleBStack = new THStack("hs","");
  doubleBStack->Add(h_WJets_stack);
  doubleBStack->Add(h_TT_stack);
  doubleBStack->Add(h_SnglT_stack);

  TString graphName = which+"_Stack";
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  styleCanvas(can_h);
  can_h->SetLogy();

  doubleBStack->Draw("hist");
  doubleBStack->GetYaxis()->SetTitle("Events");
  doubleBStack->SetMinimum(0.1); doubleBStack->SetMaximum(300);

  if (which == "DoubleB") doubleBStack->GetXaxis()->SetTitle("Deep bb-tag");
  else if (which == "tau32") doubleBStack->GetXaxis()->SetTitle("Tau32");
  else if (which == "leadmass") doubleBStack->GetXaxis()->SetTitle("Lead jet soft drop mass [GeV]");
  else if (which == "subleadmass") doubleBStack->GetXaxis()->SetTitle("Sublead jet soft drop mass [GeV]");
  else if (which == "leadpt") doubleBStack->GetXaxis()->SetTitle("Lead jet p_{T} [GeV]");
  else if (which == "subleadpt") doubleBStack->GetXaxis()->SetTitle("Sublead jet p_{T} [GeV]");

  TLegend* legend = new TLegend(0.35,0.8,0.9,0.9) ;
  legend->AddEntry(h_TT_stack, "TT", "f");
  legend->AddEntry(h_WJets_stack, "WJets", "f");
  legend->AddEntry(h_SnglT_stack, "SnglT", "f");
  legend->SetNColumns(4);
  legend->SetBorderSize(0);
  legend->Draw("same");
  can_h->Update(); can_h->Modified();

  TString savename = which+"_SingleLept_Stack";
  cdOther->cd();
  can_h->Write(savename);
  if (savePDFs) can_h->SaveAs("ANPlots/"+savename+".pdf","PDF");
  delete can_h;
}

void makeMETStack(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets,TH1F* h_SnglT, TH1F* h_T5HH1600,TH1F* h_T5HH2000,TH1F* h_T5HH2200,TH1F* h_TChiHH600,TH1F* h_TChiHH800, TH1F* h_TChiHH1000,TString tagType = ""){
  h_QCD->SetFillColor(kGray); h_QCD->SetMarkerStyle(21); h_QCD->SetMarkerColor(kGray);
  h_TT->SetFillColor(kCyan); h_TT->SetMarkerStyle(21); h_TT->SetMarkerColor(kCyan);
  h_SnglT->SetFillColor(kCyan+2); h_SnglT->SetMarkerStyle(21); h_SnglT->SetMarkerColor(kCyan+2);
  h_WJets->SetFillColor(kBlue);h_WJets->SetMarkerStyle(21); h_WJets->SetMarkerColor(kBlue);
  h_ZJets->SetFillColor(kGreen+2);h_ZJets->SetMarkerStyle(21); h_ZJets->SetMarkerColor(kGreen+2);

  h_T5HH1600->SetMarkerStyle(1); h_T5HH1600->SetMarkerColor(kRed); h_T5HH1600->SetLineColor(kRed);
  h_T5HH2000->SetMarkerStyle(1); h_T5HH2000->SetMarkerColor(kBlack); h_T5HH2000->SetLineColor(kBlack);
  h_T5HH2200->SetMarkerStyle(1); h_T5HH2200->SetMarkerColor(kMagenta); h_T5HH2200->SetLineColor(kMagenta);

  h_TChiHH600->SetMarkerStyle(1); h_TChiHH600->SetMarkerColor(kRed); h_TChiHH600->SetLineColor(kRed); h_TChiHH600->SetLineStyle(2);
  h_TChiHH800->SetMarkerStyle(1); h_TChiHH800->SetMarkerColor(kBlack); h_TChiHH800->SetLineColor(kBlack); h_TChiHH800->SetLineStyle(2);
  h_TChiHH1000->SetMarkerStyle(1); h_TChiHH1000->SetMarkerColor(kMagenta); h_TChiHH1000->SetLineColor(kMagenta); h_TChiHH1000->SetLineStyle(2);


  THStack * METStack = new THStack("hs","");
  METStack->Add(h_ZJets);
  METStack->Add(h_WJets);
  METStack->Add(h_TT);
  METStack->Add(h_SnglT);
  METStack->Add(h_QCD);


  TString graphName = "METStack_"+tagType;
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  styleCanvas(can_h);
  can_h->SetLogy();

  METStack->Draw("hist");
  METStack->SetTitle(";MET [GeV];Events");
  METStack->SetMinimum(0.1); METStack->SetMaximum(300);

  h_T5HH1600->Draw("same l"); h_T5HH2000->Draw("same l"); h_T5HH2200->Draw("same l");
  h_TChiHH600->Draw("same l"); h_TChiHH800->Draw("same l"); h_TChiHH1000->Draw("same l");

  TLegend* legend = new TLegend(0.35,0.8,0.9,0.9) ;
  legend->AddEntry(h_QCD, "QCD", "f");
  legend->AddEntry(h_TT, "TT", "f");
  legend->AddEntry(h_SnglT, "SnglT", "f");
  legend->AddEntry(h_WJets, "WJets", "f");
  legend->AddEntry(h_ZJets, "ZJets", "f");
  legend->AddEntry(h_T5HH1600, "T5HH1600", "l");
  legend->AddEntry(h_T5HH2000, "T5HH2000", "l");
  legend->AddEntry(h_T5HH2200, "T5HH2200", "l");
  legend->AddEntry(h_TChiHH600, "TChiHH600", "l");
  legend->AddEntry(h_TChiHH800, "TChiHH800", "l");
  legend->AddEntry(h_TChiHH1000, "TChiHH1000", "l");

  legend->SetNColumns(4);
  legend->SetBorderSize(0);

  legend->Draw("same");
  can_h->Update(); can_h->Modified();

  TString savename = "MET"+tagType+"_Stack";
  cdOther->cd();
  can_h->Write(savename);
  if (savePDFs) can_h->SaveAs("ANPlots/"+savename+".pdf","PDF");
  delete can_h;
}

void makeMETNorm(vector<TH1F*> dem_histos, TString bkgType = "", TString tagType = "") {
  //dem_histos will be [A,B,C,COpt,D,DOpt]
  TH1F * h_A   = (TH1F*)dem_histos[0]->Clone("h_A"); TH1F * h_B   = (TH1F*)dem_histos[1]->Clone("h_B");
  TH1F * h_C = (TH1F*)dem_histos[2]->Clone("h_C"); TH1F * h_COpt = (TH1F*)dem_histos[3]->Clone("h_COpt");
  TH1F * h_D = (TH1F*)dem_histos[4]->Clone("h_D"); TH1F * h_DOpt = (TH1F*)dem_histos[5]->Clone("h_DOpt");

  DrawOverflow(h_A); DrawOverflow(h_B); DrawOverflow(h_C); DrawOverflow(h_D); DrawOverflow(h_COpt); DrawOverflow(h_DOpt);

  TH1F *h_0H = (TH1F*)h_C->Clone("h_0H"); h_0H->Add(h_D);
  TH1F *h_0HOpt = (TH1F*)h_COpt->Clone("h_0HOpt"); h_0HOpt->Add(h_DOpt);

  h_A->SetMarkerStyle(20); h_A->SetMarkerSize(0.85); h_A->SetLineColor(kRed); h_A->SetMarkerColor(kRed);
  h_B->SetMarkerStyle(20); h_B->SetMarkerSize(0.85); h_B->SetLineColor(kMagenta); h_B->SetMarkerColor(kMagenta);
  h_C->SetMarkerStyle(20); h_C->SetMarkerSize(0.85); h_C->SetLineColor(kBlue); h_C->SetMarkerColor(kBlue);
  h_D->SetMarkerStyle(20); h_D->SetMarkerSize(0.85); h_D->SetLineColor(kGreen+2); h_D->SetMarkerColor(kGreen+2);
  h_0HOpt->SetMarkerStyle(20); h_0HOpt->SetMarkerSize(0.85); h_0HOpt->SetLineColor(kBlack); h_0HOpt->SetMarkerColor(kBlack);


  TString plotName;
  TString thisRegion;
  TString thisCut;

  if (whichRegion=="signal") thisRegion = "0-lepton";
  else if (whichRegion=="singleLept") thisRegion = "1-lepton";

  // if (tagType=="Single") thisCut = "BTagsL>1&&BTagsM>0";
  // else if (tagType=="Double") thisCut = "BTagsM>1";
  thisCut = "BTagsT>0";

  plotName = thisRegion+": "+bkgType+";MET [GeV];Normalized to area";
  h_A->SetTitle(plotName);
  h_0HOpt->SetTitle(plotName);

  TString graphName = bkgType+"_"+tagType+"_METShape";
  TGraphAsymmErrors * graph = new TGraphAsymmErrors(h_A, h_0HOpt, "pois");
  styleGraph(graph);
  graph->SetMinimum(-0.1);
  if (whichRegion!="photon") graph->GetXaxis()->SetRangeUser(300,1400);
  else graph->GetXaxis()->SetRangeUser(100,1400);
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
  h_A->GetXaxis()->SetLabelSize(0);
  h_A->GetYaxis()->SetTitleOffset(0.8);

  //try normalizing to 1
  h_A->Scale(1/h_A->Integral());
  h_B->Scale(1/h_B->Integral());
  h_C->Scale(1/h_C->Integral());
  h_D->Scale(1/h_D->Integral());
  h_0HOpt->Scale(1/h_0HOpt->Integral());

  h_A->Draw();
  h_B->Draw("same");
  h_C->Draw("same");
  h_D->Draw("same");
  h_0HOpt->Draw("same");

  TLegend* legend = new TLegend(0.5,0.65,0.95,0.88,NULL,"brNDC") ;
  legend->SetBorderSize(0);
  TString legendEntryName = tagType+" SR";
  if (tagType=="Single") {
    legend->AddEntry(h_A, "1H SR", "PL");
    legend->AddEntry(h_B, "1H SB", "PL");
  }
  else if (tagType=="Double") {
    legend->AddEntry(h_A, "2H SR", "PL");
    legend->AddEntry(h_B, "2H SB", "PL");
  }
  legend->AddEntry(h_C, "0H SR", "PL");
  legend->AddEntry(h_D, "0H SB", "PL");
  legend->AddEntry(h_0HOpt, "0H SB+SR,"+thisCut, "PL");

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
  can_h->Write(graphName);
  delete can_h;
}

void makeMETNormCompare(vector<TH1F*> dem_histos, TString bkgType = "", TString tagType = "") {
  //dem_histos will be [SR,0HSB,0HSBCuts,0HSR,0HSRCuts]
  TH1F * h_SR   = (TH1F*)dem_histos[0]->Clone("h_SR");
  TH1F * h_SB   = (TH1F*)dem_histos[1]->Clone("h_SB");
  TH1F * h_0HSR_noCuts = (TH1F*)dem_histos[2]->Clone("h_0HSR_noCuts");
  TH1F * h_0HSR_withCuts = (TH1F*)dem_histos[3]->Clone("h_0HSR_withCuts");
  TH1F * h_0HSB_noCuts = (TH1F*)dem_histos[4]->Clone("h_0HSB_noCuts");
  TH1F * h_0HSB_withCuts = (TH1F*)dem_histos[5]->Clone("h_0HSB_withCuts");

  TH1F * h_SRNoScale   = (TH1F*)dem_histos[0]->Clone("h_SRNoScale");
  //Adding SR and SB together
  h_0HSB_noCuts->Add(h_0HSR_noCuts);
  h_0HSB_withCuts->Add(h_0HSR_withCuts);

  DrawOverflow(h_SR); DrawOverflow(h_SB); DrawOverflow(h_0HSB_noCuts); DrawOverflow(h_0HSB_withCuts);

  h_SR->SetMarkerStyle(20); h_SR->SetMarkerSize(0.85); h_SR->SetLineColor(kRed); h_SR->SetMarkerColor(kRed);
  h_0HSB_withCuts->SetMarkerStyle(20); h_0HSB_withCuts->SetMarkerSize(0.85); h_0HSB_withCuts->SetLineColor(kBlack); h_0HSB_withCuts->SetMarkerColor(kBlack);
  h_0HSB_noCuts->SetMarkerStyle(20); h_0HSB_noCuts->SetMarkerSize(0.85); h_0HSB_noCuts->SetLineColor(kBlue); h_0HSB_noCuts->SetMarkerColor(kBlue);

  // if (tagType=="Single") h_SR->SetTitle("BTagsL>1&&BTagsM>0;MET [GeV];Normalized to area");
  // else h_SR->SetTitle("BTagsM>1;MET [GeV];Normalized to area");

  // h_SR->SetTitle("BTagsT>0;MET [GeV];Normalized to area");
  // if (whichRegion == "singleLept") h_SR->SetTitle("1-l, BTagsT>0;MET [GeV];Normalized to area");
  h_SR->SetTitle(";MET [GeV];Normalized to area");
  h_SR->GetYaxis()->SetTitleOffset(0.8);
  h_SR->GetXaxis()->SetTitleSize(0);

  // if (whichRegion == "singleLept" && tagType=="Single") h_SR->SetTitle("1-l, BTagsL>1&&BTagsM>0;MET [GeV];Normalized to area");
  // else if (whichRegion == "singleLept" && tagType=="Double") h_SR->SetTitle("1-l, BTagsM>1;MET [GeV];Normalized to area");

  TString graphName = bkgType+"_"+tagType+"_METShape";
  TGraphAsymmErrors * graph = new TGraphAsymmErrors(h_SR, h_0HSB_withCuts, "pois");
  styleGraph(graph);
  graph->SetMinimum(-0.1);
  if (whichRegion!="photon") graph->GetXaxis()->SetRangeUser(300,1400);
  else graph->GetXaxis()->SetRangeUser(100,1400);
  graph->SetTitle(";MET [GeV];2HSR / 0HCuts");
  if (tagType=="single") graph->SetTitle(";MET [GeV];1HSR / 0HCuts");

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
  if (showNoCutLine) legend->AddEntry(h_0HSB_noCuts, "0H SB+SR, no cuts", "PL");
  // if (tagType=="Single") legend->AddEntry(h_0HSB_withCuts, "0H SB+SR,BTagsL>1&&BTagsM>0", "PL");
  // else if (tagType=="Double") legend->AddEntry(h_0HSB_withCuts, "0H SB+SR,BTagsM>1", "PL");
  legend->AddEntry(h_0HSB_withCuts, "0H SB+SR,BTagsT>0", "PL");
  legend->Draw();

  pad2->cd();
  graph->Draw("APE");
  graph->GetYaxis()->SetTitleOffset(0.25); graph->GetYaxis()->SetTitleSize(0.13);
  graph->GetXaxis()->SetTitleOffset(0.9); graph->GetXaxis()->SetTitleSize(0.17);
  graph->GetXaxis()->SetLabelSize(0.14); graph->GetYaxis()->SetLabelSize(0.10);
  graph->GetYaxis()->SetNdivisions(505);
  graph->SetMinimum(0.0);
  graph->Draw("APE");


  double p0; double error;
  if (whichRegion!="photon"){
    TF1*f0=new TF1("f0","pol0", 300,1400);
    graph->Fit("f0","Q","R",300,1400);
    p0 = f0->GetParameter(0); //this does, indeed, work
    error = f0->GetParError(0);
    graph->GetFunction("f0")->SetBit(TF1::kNotDraw);
  }
  else {
    TF1*f0=new TF1("f0","pol0", 100,1400);
    graph->Fit("f0","Q","R",100,1400);
    p0 = f0->GetParameter(0); //this does, indeed, work
    error = f0->GetParError(0);
    graph->GetFunction("f0")->SetBit(TF1::kNotDraw);
  }

  TLine *line = new TLine(300.0,p0,1400.0,p0);
  line->SetLineColor(kRed);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  std::string lineConst = to_string(p0); std::string lineErr = to_string(error);
  TString constString = "Const: "+lineConst+" #pm "+ lineErr;
  TLatex *t = new TLatex(0.2,0.85,constString);
  t->SetNDC(); t->SetTextFont(52);
  t->SetTextSize(0.13);


  //Change for MET shape
  if (whichRegion=="singleLept" && bkgType=="BkgSum"){
    const int n = 3;
    double x[n]   = {400.0,600.0,1050.0};
    double y[n]   = {0.0,0.0,0.0};
    double exl[n] = {100.0,100.0,350.0};
    double exh[n] = {100.0,100.0,350.0};
    double eyl[n] = {0.0,0.0,0.0};
    double eyh[n] = {0.0,0.0,0.0};

    for (int i=0;i<n;i++){
      y[i] = p0;
      double yieldsPercent = sqrt(h_SRNoScale->GetBinContent(i+1))/h_SRNoScale->GetBinContent(i+1);
      double err = p0-p0*(1.0-yieldsPercent);
      double lowErr = err; if (lowErr>p0) lowErr = p0;
      eyl[i] = lowErr;
      eyh[i] = err;
    }
    TGraphAsymmErrors * graph_errorBands = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

    graph_errorBands->SetFillColor(kBlue); graph_errorBands->SetFillStyle(3445);
    graph_errorBands->SetMarkerSize(0); graph_errorBands->SetLineWidth(1);
    graph_errorBands->SetLineColor(kBlue);
    graph_errorBands->Draw("E2");
  }
  else if (whichRegion=="photon" && bkgType=="BkgSum"){
    const int n = 4;
    double x[n]   = {200.0,400.0,600.0,1050.0};
    double y[n]   = {0.0,0.0,0.0,0.0};
    double exl[n] = {100.0,100.0,100.0,350.0};
    double exh[n] = {100.0,100.0,100.0,350.0};
    double eyl[n] = {0.0,0.0,0.0,0.0};
    double eyh[n] = {0.0,0.0,0.0,0.0};

    for (int i=0;i<n;i++){
      y[i] = p0;
      double yieldsPercent = sqrt(h_SRNoScale->GetBinContent(i+1))/h_SRNoScale->GetBinContent(i+1);
      double err = p0-p0*(1.0-yieldsPercent);
      double lowErr = err; if (lowErr>p0) lowErr = p0;
      eyl[i] = lowErr;
      eyh[i] = err;
    }
    TGraphAsymmErrors * graph_errorBands = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

    graph_errorBands->SetFillColor(kBlue); graph_errorBands->SetFillStyle(3445);
    graph_errorBands->SetMarkerSize(0); graph_errorBands->SetLineWidth(1);
    graph_errorBands->SetLineColor(kBlue);
    graph_errorBands->Draw("E2");
  }
  else {
    line->Draw("same");
    t->Draw();
  }

  cdOther->cd();
  can_h->Write(graphName+"_Comp");
  TString thisRegion = "0l"; if (whichRegion=="singleLept")thisRegion="1l";
  TString thisSR = "2H"; if (tagType =="single")thisSR="1H";
  if (savePDFs && bkgType!="BkgSum") can_h->SaveAs("ANPlots/METShape_"+thisRegion+"_"+thisSR+"_"+bkgType+".pdf","PDF");
  else if (savePDFs) can_h->SaveAs("ANPlots/METShape_"+thisRegion+"_"+thisSR+".pdf","PDF");
  delete can_h;

  if (bkgType =="BkgSum") tableOfMETNorm(dem_histos, bkgType, tagType);
}

void tableOfYields(vector<TH1F*> dem_histos, TString bkgType, TString tagType){  //A, B, C, D
  TH1F * histo_A = (TH1F*)dem_histos[0]->Clone("histo_A");
  TH1F * histo_B = (TH1F*)dem_histos[1]->Clone("histo_B");
  TH1F * histo_C = (TH1F*)dem_histos[2]->Clone("histo_C");
  TH1F * histo_D = (TH1F*)dem_histos[3]->Clone("histo_D");

  DrawOverflow(histo_A); DrawOverflow(histo_B); DrawOverflow(histo_C); DrawOverflow(histo_D);

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
  if (tagType=="Double") yieldsFileName = "Yields/2HYields_"+whichRegion+"_"+bkgType+"TightCuts.txt";
  else if (tagType=="Single") yieldsFileName = "Yields/1HYields_"+whichRegion+"_"+bkgType+"TightCuts.txt";

  ofstream yields;
  yields.open(yieldsFileName);

  yields<<"% "<< bkgType <<endl;
  yields <<"\\documentclass[11pt, oneside]{article}\n\n";
  yields<<"\\begin{document}\n"<<"\\begin{table}\n";
  yields<<"\\begin{tabular}{ |c|c| }\n"<<"\\hline\n";
  yields<<"Region &MC Yields \\\\  \\hline \\hline\n";
  if (tagType=="Double") {
    yields<<"2H SR & "<< a_int <<" $\\pm$ "<< a_error <<"\\\\ \n";
    yields<<"2H SB & "<< b_int <<" $\\pm$ "<< b_error <<"\\\\ \n";
  }
  else if (tagType=="Single") {
    yields<<"1H SR & "<< a_int <<" $\\pm$ "<< a_error <<"\\\\ \n";
    yields<<"1H SB & "<< b_int <<" $\\pm$ "<< b_error <<"\\\\ \n";
  }
  yields<<"0H SR & "<< c_int <<" $\\pm$ "<< c_error <<"\\\\ \n";
  yields<<"0H SB & "<< d_int <<" $\\pm$ "<< d_error <<"\\\\ \n";
  yields<<"Prediction & "<< bkgNorm <<" $\\pm$ "<< bkgNorm_error <<"\\\\ \\hline\n";
  yields<<"\\hline\n";
  yields<<"Kappa & "<< kappa << " $\\pm$ "<< kappa_error <<"\\\\ \\hline\n";
  yields<<"\\end{tabular}\n"<<"\\end{table}\n"<<"\\end{document}\n";
  yields.close();
}

void tableOfYieldsPure(vector<TH1F*> dem_histos, TString bkgType, TString tagType){  //A, B, C, D
  //Need to report yields and kappa per MET bin
  TH1F * histo_A = (TH1F*)dem_histos[0]->Clone("histo_A");
  TH1F * histo_B = (TH1F*)dem_histos[1]->Clone("histo_B");
  TH1F * histo_C = (TH1F*)dem_histos[2]->Clone("histo_C");
  TH1F * histo_D = (TH1F*)dem_histos[3]->Clone("histo_D");

  DrawOverflow(histo_A); DrawOverflow(histo_B); DrawOverflow(histo_C); DrawOverflow(histo_D);

  TH1F *histo_pred = (TH1F*)histo_B->Clone("histo_pred");
  histo_pred->Multiply(histo_C);histo_pred->Divide(histo_D);
  DrawOverflow(histo_pred); //Make sure overflow bin is merged with last MET bin


  double a_int; double a_error; a_int = histo_A->IntegralAndError(0,5,a_error,"");
  double b_int; double b_error; b_int = histo_B->IntegralAndError(0,5,b_error,"");
  double c_int; double c_error; c_int = histo_C->IntegralAndError(0,5,c_error,"");
  double d_int; double d_error; d_int = histo_D->IntegralAndError(0,5,d_error,"");
  double pred_int; double pred_error; pred_int = histo_pred->IntegralAndError(0,5,pred_error,"");


  TH1F *h_kappas = (TH1F*)histo_A->Clone("h_kappas");
  h_kappas->Divide(histo_pred);
  DrawOverflow(h_kappas); //Make sure overflow bin is merged with last MET bin

  // double kappa1 = h_kappas->GetBinContent(1); double kappa2 = h_kappas->GetBinContent(2); double kappa3 = h_kappas->GetBinContent(3);
  // double kappa_err1 = h_kappas->GetBinError(1); double kappa_err2 = h_kappas->GetBinError(2); double kappa_err3 = h_kappas->GetBinError(3);

  TString yieldsFileName;
  if (tagType=="Double") yieldsFileName = "Yields/2HPureYields_"+whichRegion+"_"+bkgType+"TightCuts.txt";
  else if (tagType=="Single") yieldsFileName = "Yields/1HPureYields_"+whichRegion+"_"+bkgType+"TightCuts.txt";

  ofstream yields;
  yields.open(yieldsFileName);

  yields<<"% Pure ABCD, per MET bin"<<endl;
  yields<<"% Background: "<< bkgType <<endl;
  yields<<"% Region: "<< whichRegion <<endl;

  yields <<"\\documentclass[12pt]{article}\n\n";
  yields<<"\\usepackage{numprint}\n"<<"\\begin{document}\n"<<"\\begin{table}\n";
  yields<<"\\npdecimalsign{.}\n"<<"\\nprounddigits{2}\n"<<"\\centering\n";
  yields<<"\\begin{tabular}{ c|n{4}{2}n{2}{2}|n{2}{2}n{2}{2}|n{2}{2}n{2}{2} }\n"<<"\\hline\n";
  yields<<"MET region & \\multicolumn{2}{c|}{Region A} & \\multicolumn{2}{c|}{Prediction} & \\multicolumn{2}{c|}{Kappa} \\\\ \n";
  yields<<"\\hline \\hline\n";
  if (tagType=="Double") yields<<"\\multicolumn{7}{c}{2H SR} \\\\ \n";
  else if (tagType=="Single") yields<<"\\multicolumn{7}{c}{1H SR} \\\\ \n";
  yields<<"\\hline\n";
  yields<<"All MET & "<<  a_int<<"& \\pm "<< a_error << "& " << pred_int<<" & \\pm "<<pred_error<<"&  & \\\\ \n";
  yields<<"$[300,500]$ GeV & "<<  histo_A->GetBinContent(1)<<"& \\pm "<< histo_A->GetBinError(1)<< "& " << histo_pred->GetBinContent(1)<<" & \\pm "<<histo_pred->GetBinError(1)<< "& " << h_kappas->GetBinContent(1)<<" & \\pm "<<h_kappas->GetBinError(1)<<" \\\\ \n";
  yields<<"$[500,700]$ GeV & "<<  histo_A->GetBinContent(2)<<"& \\pm "<< histo_A->GetBinError(2)<< "& " << histo_pred->GetBinContent(2)<<" & \\pm "<<histo_pred->GetBinError(2)<< "& " << h_kappas->GetBinContent(2)<<" & \\pm "<<h_kappas->GetBinError(2)<<" \\\\ \n";
  yields<<"$[700+]$ GeV & "   <<  histo_A->GetBinContent(3)<<"& \\pm "<< histo_A->GetBinError(3)<< "& " << histo_pred->GetBinContent(3)<<" & \\pm "<<histo_pred->GetBinError(3)<< "& " << h_kappas->GetBinContent(3)<<" & \\pm "<<h_kappas->GetBinError(3)<<" \\\\ \n";
  yields<<"\\hline\n";
  yields<<"\\multicolumn{7}{c}{Regions B, C, and D} \\\\ \n";
  yields<<"\\hline\n";

  yields<<"Region B & "<<  b_int<<"& \\pm "<< b_error << "&  &  &  & \\\\ \n";
  yields<<"$[300,500]$ GeV & "<<  histo_B->GetBinContent(1)<<"& \\pm "<< histo_B->GetBinError(1)<< "&  &  &  &  \\\\ \n";
  yields<<"$[500,700]$ GeV & "<<  histo_B->GetBinContent(2)<<"& \\pm "<< histo_B->GetBinError(2)<< "&  &  &  &  \\\\ \n";
  yields<<"$[700+]$ GeV & "<<  histo_B->GetBinContent(3)<<"& \\pm "<< histo_B->GetBinError(3)<< "&  &  &  &  \\\\ \n";
  yields<<"\\hline\n";

  yields<<"Region C & "<<  c_int<<"& \\pm "<< c_error << "&  &  &  & \\\\ \n";
  yields<<"$[300,500]$ GeV & "<<  histo_C->GetBinContent(1)<<"& \\pm "<< histo_C->GetBinError(1)<< "&  &  &  &  \\\\ \n";
  yields<<"$[500,700]$ GeV & "<<  histo_C->GetBinContent(2)<<"& \\pm "<< histo_C->GetBinError(2)<< "&  &  &  &  \\\\ \n";
  yields<<"$[700+]$ GeV & "<<  histo_C->GetBinContent(3)<<"& \\pm "<< histo_C->GetBinError(3)<< "&  &  &  &  \\\\ \n";
  yields<<"\\hline\n";

  yields<<"Region D & "<<  d_int<<"& \\pm "<< d_error << "&  &  &  & \\\\ \n";
  yields<<"$[300,500]$ GeV & "<<  histo_D->GetBinContent(1)<<"& \\pm "<< histo_D->GetBinError(1)<< "&  &  &  &  \\\\ \n";
  yields<<"$[500,700]$ GeV & "<<  histo_D->GetBinContent(2)<<"& \\pm "<< histo_D->GetBinError(2)<< "&  &  &  &  \\\\ \n";
  yields<<"$[700+]$ GeV & "<<  histo_D->GetBinContent(3)<<"& \\pm "<< histo_D->GetBinError(3)<< "&  &  &  &  \\\\ \n";
  yields<<"\\hline\n";
  yields<<"\\end{tabular}\n"<<"\\end{table}\n"<<"\\end{document}\n";
  yields.close();
}

void tableOfMETNorm(vector<TH1F*> dem_histos, TString bkgType, TString tagType){  //A, D, DOpt, C, COpt
  TH1F * histo_A = (TH1F*)dem_histos[0]->Clone("histo_A");
  TH1F * histo_B = (TH1F*)dem_histos[1]->Clone("histo_B");
  TH1F * histo_CNoCuts = (TH1F*)dem_histos[2]->Clone("histo_CNoCuts");
  TH1F * histo_CWithCuts = (TH1F*)dem_histos[3]->Clone("histo_CWithCuts");
  TH1F * histo_DNoCuts = (TH1F*)dem_histos[4]->Clone("histo_DNoCuts");
  TH1F * histo_DWithCuts = (TH1F*)dem_histos[5]->Clone("histo_DWithCuts");

  DrawOverflow(histo_A); DrawOverflow(histo_B); DrawOverflow(histo_CNoCuts); DrawOverflow(histo_CWithCuts); DrawOverflow(histo_DNoCuts); DrawOverflow(histo_DWithCuts);

  TH1F *h_2H_Norm = (TH1F*) histo_A->Clone("h_2H_Norm");
  h_2H_Norm->Scale(1/h_2H_Norm->Integral(1,4));
  TH1F *h_0H_Norm = (TH1F*) histo_DWithCuts->Clone("h_0H_Norm");
  h_0H_Norm->Add(histo_CWithCuts);
  h_0H_Norm->Scale(1/h_0H_Norm->Integral(1,4));

  histo_DWithCuts->Add(histo_CWithCuts);

  TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo_A, histo_DWithCuts, "pois");
  double p0; double error;
  if (whichRegion=="photon") {
    TF1*f0=new TF1("f0","pol0", 100,1400);
    graph->Fit("f0","Q","R",100,1400);
    p0 = f0->GetParameter(0); //this does, indeed, work
    error = f0->GetParError(0);
  }
  else {
    TF1*f0=new TF1("f0","pol0", 300,1400);
    graph->Fit("f0","Q","R",300,1400);
    p0 = f0->GetParameter(0); //this does, indeed, work
    error = f0->GetParError(0);
  }


  double bin1Obs = histo_A->GetBinContent(1)/histo_DWithCuts->GetBinContent(1);
  double bin2Obs = histo_A->GetBinContent(2)/histo_DWithCuts->GetBinContent(2);
  double bin3Obs = histo_A->GetBinContent(3)/histo_DWithCuts->GetBinContent(3);
  double bin4Obs; double bin4Dev;

  double bin1Dev = 1.0 + (bin1Obs-p0)/p0;
  double bin2Dev = 1.0 + (bin2Obs-p0)/p0;
  double bin3Dev = 1.0 + (bin3Obs-p0)/p0;
  if (whichRegion=="photon") {
    bin4Obs = histo_A->GetBinContent(4)/histo_DWithCuts->GetBinContent(4);
    bin4Dev = 1.0 + (bin4Obs-p0)/p0;
  }


  TString yieldsFileName;
  if (tagType=="Double") yieldsFileName = "Yields/2HMETShape_"+whichRegion+"_"+bkgType+"TightCuts.txt";
  else if (tagType=="Single") yieldsFileName = "Yields/1HMETShape_"+whichRegion+"_"+bkgType+"TightCuts.txt";

  ofstream yields;
  yields.open(yieldsFileName);

  yields<<"% Background: "<< bkgType <<endl;
  yields<<"% Region: "<< whichRegion <<endl;
  yields<<"% Best fit: "<< p0 <<" +/- "<<error <<endl;

  yields <<"\\documentclass[12pt]{article}\n\n";
  yields<<"\\usepackage{numprint}\n"<<"\\begin{document}\n"<<"\\begin{table}\n";
  yields<<"\\npdecimalsign{.}\n"<<"\\nprounddigits{3}\n"<<"\\centering\n";
  if (whichRegion=="signal"){
    yields<<"\\begin{tabular}{ |c|n{4}{3}n{2}{3}|n{3}{3}n{2}{3}| }\n"<<"\\hline\n";
    yields<<"MET region & \\multicolumn{2}{c|}{MC Yields} & \\multicolumn{2}{c|}{Bkg Frac.}  \\\\ \n";
    yields<<"\\hline \\hline\n";
    if (tagType=="Double") yields<<"\\multicolumn{5}{|c|}{2H SR} \\\\ \n";
    else if (tagType=="Single") yields<<"\\multicolumn{5}{|c|}{1H SR} \\\\ \n";
    yields<<"\\hline\n";
    yields<<"$[300,500]$ GeV & "<<  histo_A->GetBinContent(1)<<"& \\pm "<< histo_A->GetBinError(1)<< "& " << h_2H_Norm->GetBinContent(1)<<" & \\pm "<<h_2H_Norm->GetBinError(1)<<" \\\\ \n";
    yields<<"$[500,700]$ GeV & "<<  histo_A->GetBinContent(2)<<"& \\pm "<< histo_A->GetBinError(2)<< "& " << h_2H_Norm->GetBinContent(2)<<" & \\pm "<<h_2H_Norm->GetBinError(2)<<" \\\\ \n";
    yields<<"$[700+]$ GeV & "<<  histo_A->GetBinContent(3)<<"& \\pm "<< histo_A->GetBinError(3)<< "& " << h_2H_Norm->GetBinContent(3)<<" & \\pm "<<h_2H_Norm->GetBinError(3)<<" \\\\ \n";
    yields<<"\\hline\n";
    yields<<"\\multicolumn{5}{|c|}{0H, $B_{T}>$0} \\\\ \n";
    // if (tagType=="Double") yields<<"\\multicolumn{5}{c}{0H, $B_{M}>$1} \\\\ \n";
    // else if (tagType=="Single") yields<<"\\multicolumn{5}{c}{0H, $B_{L}>$1 and $B_{M}>$0} \\\\ \n";
    yields<<"\\hline\n";
    yields<<"$[300,500]$ GeV & "<<  histo_DWithCuts->GetBinContent(1)<<" & \\pm "<< histo_DWithCuts->GetBinError(1)<< "& " << h_0H_Norm->GetBinContent(1)<<" & \\pm "<<h_0H_Norm->GetBinError(1)<<"\\\\ \n";
    yields<<"$[500,700]$ GeV & "<<  histo_DWithCuts->GetBinContent(2)<<" & \\pm "<< histo_DWithCuts->GetBinError(2)<< "& " << h_0H_Norm->GetBinContent(2)<<" & \\pm "<<h_0H_Norm->GetBinError(2)<<" \\\\ \n";
    yields<<"$[700+]$ GeV & "<<  histo_DWithCuts->GetBinContent(3)<<" & \\pm "<< histo_DWithCuts->GetBinError(3)<< "& " << h_0H_Norm->GetBinContent(3)<<" & \\pm "<<h_0H_Norm->GetBinError(3)<<" \\\\ \n";
  }
  else if (whichRegion=="singleLept"){
    yields<<"\\begin{tabular}{ |c|n{4}{3}n{2}{3}|n{3}{3}n{2}{3}|n{4}{3}| }\n"<<"\\hline\n";
    yields<<"MET region & \\multicolumn{2}{c|}{MC Yields} & \\multicolumn{2}{c|}{Bkg Frac.} & \\multicolumn{1}{c|}{Deviation}  \\\\ \n";
    yields<<"\\hline \\hline\n";
    if (tagType=="Double") yields<<"\\multicolumn{6}{|c|}{2H SR} \\\\ \n";
    else if (tagType=="Single") yields<<"\\multicolumn{6}{|c|}{1H SR} \\\\ \n";
    yields<<"\\hline\n";
    yields<<"$[300,500]$ GeV & "<<  histo_A->GetBinContent(1)<<"& \\pm "<< histo_A->GetBinError(1)<< "& " << h_2H_Norm->GetBinContent(1)<<" & \\pm "<<h_2H_Norm->GetBinError(1)<<" & "<< " \\\\ \n";
    yields<<"$[500,700]$ GeV & "<<  histo_A->GetBinContent(2)<<"& \\pm "<< histo_A->GetBinError(2)<< "& " << h_2H_Norm->GetBinContent(2)<<" & \\pm "<<h_2H_Norm->GetBinError(2)<<" & "<< " \\\\ \n";
    yields<<"$[700+]$ GeV & "<<  histo_A->GetBinContent(3)<<"& \\pm "<< histo_A->GetBinError(3)<< "& " << h_2H_Norm->GetBinContent(3)<<" & \\pm "<<h_2H_Norm->GetBinError(3)<<" & "<< " \\\\ \n";
    yields<<"\\hline\n";
    yields<<"\\multicolumn{6}{|c|}{0H, $B_{T}>$0} \\\\ \n";
    // if (tagType=="Double") yields<<"\\multicolumn{6}{c}{0H, $B_{M}>$1} \\\\ \n";
    // else if (tagType=="Single") yields<<"\\multicolumn{6}{c}{0H, $B_{M}>$0} \\\\ \n";
    yields<<"\\hline\n";
    yields<<"$[300,500]$ GeV & "<<  histo_DWithCuts->GetBinContent(1)<<" & \\pm "<< histo_DWithCuts->GetBinError(1)<< "& " << h_0H_Norm->GetBinContent(1)<<" & \\pm "<<h_0H_Norm->GetBinError(1)<<" & "<< bin1Dev<<" \\\\ \n";
    yields<<"$[500,700]$ GeV & "<<  histo_DWithCuts->GetBinContent(2)<<" & \\pm "<< histo_DWithCuts->GetBinError(2)<< "& " << h_0H_Norm->GetBinContent(2)<<" & \\pm "<<h_0H_Norm->GetBinError(2)<<" & "<< bin2Dev<<" \\\\ \n";
    yields<<"$[700+]$ GeV & "<<  histo_DWithCuts->GetBinContent(3)<<" & \\pm "<< histo_DWithCuts->GetBinError(3)<< "& " << h_0H_Norm->GetBinContent(3)<<" & \\pm "<<h_0H_Norm->GetBinError(3)<<" & "<< bin3Dev<< " \\\\ \n";
  }
  else if (whichRegion=="photon"){
    yields<<"\\begin{tabular}{ |c|n{4}{3}n{2}{3}|n{3}{3}n{2}{3}|n{4}{3}| }\n"<<"\\hline\n";
    yields<<"MET region & \\multicolumn{2}{c|}{MC Yields} & \\multicolumn{2}{c|}{Bkg Frac.} & \\multicolumn{1}{c|}{Deviation}  \\\\ \n";
    yields<<"\\hline \\hline\n";
    if (tagType=="Double") yields<<"\\multicolumn{6}{|c|}{2H SR} \\\\ \n";
    else if (tagType=="Single") yields<<"\\multicolumn{6}{|c|}{1H SR} \\\\ \n";
    yields<<"\\hline\n";
    yields<<"$[100,300]$ GeV & "<<  histo_A->GetBinContent(1)<<"& \\pm "<< histo_A->GetBinError(1)<< "& " << h_2H_Norm->GetBinContent(1)<<" & \\pm "<<h_2H_Norm->GetBinError(1)<<" & "<< " \\\\ \n";
    yields<<"$[300,500]$ GeV & "<<  histo_A->GetBinContent(2)<<"& \\pm "<< histo_A->GetBinError(2)<< "& " << h_2H_Norm->GetBinContent(2)<<" & \\pm "<<h_2H_Norm->GetBinError(2)<<" & "<< " \\\\ \n";
    yields<<"$[500,700]$ GeV & "<<  histo_A->GetBinContent(3)<<"& \\pm "<< histo_A->GetBinError(3)<< "& " << h_2H_Norm->GetBinContent(3)<<" & \\pm "<<h_2H_Norm->GetBinError(3)<<" & "<< " \\\\ \n";
    yields<<"$[700+]$ GeV & "<<  histo_A->GetBinContent(4)<<"& \\pm "<< histo_A->GetBinError(4)<< "& " << h_2H_Norm->GetBinContent(4)<<" & \\pm "<<h_2H_Norm->GetBinError(4)<<" & "<< " \\\\ \n";
    yields<<"\\hline\n";
    yields<<"\\multicolumn{6}{|c|}{0H, $B_{T}>$0} \\\\ \n";
    // if (tagType=="Double") yields<<"\\multicolumn{6}{c}{0H, $B_{M}>$1} \\\\ \n";
    // else if (tagType=="Single") yields<<"\\multicolumn{6}{c}{0H, $B_{M}>$0} \\\\ \n";
    yields<<"\\hline\n";
    yields<<"$[100,300]$ GeV & "<<  histo_DWithCuts->GetBinContent(1)<<" & \\pm "<< histo_DWithCuts->GetBinError(1)<< "& " << h_0H_Norm->GetBinContent(1)<<" & \\pm "<<h_0H_Norm->GetBinError(1)<<" & "<< bin1Dev<<" \\\\ \n";
    yields<<"$[300,500]$ GeV & "<<  histo_DWithCuts->GetBinContent(2)<<" & \\pm "<< histo_DWithCuts->GetBinError(2)<< "& " << h_0H_Norm->GetBinContent(2)<<" & \\pm "<<h_0H_Norm->GetBinError(2)<<" & "<< bin2Dev<<" \\\\ \n";
    yields<<"$[500,700]$ GeV & "<<  histo_DWithCuts->GetBinContent(3)<<" & \\pm "<< histo_DWithCuts->GetBinError(3)<< "& " << h_0H_Norm->GetBinContent(3)<<" & \\pm "<<h_0H_Norm->GetBinError(3)<<" & "<< bin3Dev<<" \\\\ \n";
    yields<<"$[700+]$ GeV & "<<  histo_DWithCuts->GetBinContent(4)<<" & \\pm "<< histo_DWithCuts->GetBinError(4)<< "& " << h_0H_Norm->GetBinContent(4)<<" & \\pm "<<h_0H_Norm->GetBinError(4)<<" & "<< bin4Dev<< " \\\\ \n";
  }
  else std::cout<<"Something wrong with region in table of MET norms!\n";

  yields<<"\\hline\n";
  yields<<"\\end{tabular}\n"<<"\\npnoround\n"<<"\\end{table}\n"<<"\\end{document}\n";
  yields.close();
}


void pieChart(vector<TH1F*> h_QCD, vector<TH1F*> h_WJets, vector<TH1F*> h_ZJets, vector<TH1F*> h_TTJets,vector<TH1F*> h_SnglT, TString regionLabel, TString bin){
  TH1F * h_Q = (TH1F*)h_QCD[0]->Clone("h_Q"); for (int i=1;i<h_QCD.size();i++) h_Q->Add(h_QCD[i]);
  TH1F * h_W = (TH1F*)h_WJets[0]->Clone("h_W"); for (int i=1;i<h_WJets.size();i++) h_W->Add(h_WJets[i]);
  TH1F * h_Z = (TH1F*)h_ZJets[0]->Clone("h_Z"); for (int i=1;i<h_ZJets.size();i++) h_Z->Add(h_ZJets[i]);
  TH1F * h_TT = (TH1F*)h_TTJets[0]->Clone("h_TT");  for (int i=1;i<h_TTJets.size();i++) h_TT->Add(h_TTJets[i]);
  TH1F * h_T = (TH1F*)h_SnglT[0]->Clone("h_T");  for (int i=1;i<h_SnglT.size();i++) h_T->Add(h_SnglT[i]);

  DrawOverflow(h_Q); DrawOverflow(h_W); DrawOverflow(h_Z); DrawOverflow(h_TT); DrawOverflow(h_T);

  double Q_yields;double W_yields;double Z_yields;double TT_yields;double T_yields;

  if (bin=="all"){
    Q_yields = h_Q->Integral();
    W_yields = h_W->Integral();
    Z_yields = h_Z->Integral();
    TT_yields = h_TT->Integral();
    T_yields = h_T->Integral();
  }

  if (bin=="bin1"){
    Q_yields = h_Q->GetBinContent(1);
    W_yields = h_W->GetBinContent(1);
    Z_yields = h_Z->GetBinContent(1);
    TT_yields = h_TT->GetBinContent(1);
    T_yields = h_T->GetBinContent(1);
  }

  if (bin=="bin2"){
    Q_yields = h_Q->GetBinContent(2);
    W_yields = h_W->GetBinContent(2);
    Z_yields = h_Z->GetBinContent(2);
    TT_yields = h_TT->GetBinContent(2);
    T_yields = h_T->GetBinContent(2);
  }

  if (bin=="bin3"){
    Q_yields = h_Q->GetBinContent(3)+h_Q->GetBinContent(4);
    W_yields = h_W->GetBinContent(3)+h_W->GetBinContent(4);
    Z_yields = h_Z->GetBinContent(3)+h_Z->GetBinContent(4);
    TT_yields = h_TT->GetBinContent(3)+h_TT->GetBinContent(4);
    T_yields = h_T->GetBinContent(3)+h_T->GetBinContent(4);
  }

  TString yieldsFileName;
  yieldsFileName = whichRegion+", "+regionLabel;
  TCanvas * can_h = new TCanvas(yieldsFileName,yieldsFileName, 50, 50, 1200, 1200);

  double vals[] = {Q_yields,W_yields,Z_yields,TT_yields,T_yields};
  int colors[] = {kGray,kAzure+1,kOrange-2,kGreen+3,kSpring-5};
  const char *labels[] = {"QCD","WJets","ZJets","TT","SnglT"};
  int nvals = sizeof(vals)/sizeof(vals[0]);
  TString binDef;
  if (bin=="all") binDef = "allMET";
  else if (bin=="bin1")binDef = "MET300to500";
  else if (bin=="bin2")binDef = "MET500to700";
  else if (bin=="bin3")binDef = "MET700up";
  TString pieName = yieldsFileName + ": " + binDef;
  TPie *pie1 = new TPie("pie1", pieName,nvals,vals,colors, labels);
  pie1->SetAngularOffset(355.);
  // pie1->SetTitle("Over all MET");
  pie1->SetAngle3D(45.);
  pie1->SetHeight(0.04);
  pie1->SetCircle(0.5,0.5,.25);
  pie1->SetLabelsOffset(.03);
  pie1->SetLabelFormat("%perc");
  pie1->Draw("3d");

  TString savename = whichRegion+"_"+regionLabel+"_"+binDef;
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
  if (savePDFs && whichRegion=="signal") can_h->SaveAs("ANPlots/pie_0l_"+regionLabel+"_"+binDef+".pdf","PDF");

  delete pie1;
  delete can_h;
}

void pieChart1l(vector<TH1F*> h_WJets, vector<TH1F*> h_TTJets, vector<TH1F*> h_SnglT, TString regionLabel, TString bin){
  TH1F * h_W = (TH1F*)h_WJets[0]->Clone("h_W"); for (int i=1;i<h_WJets.size();i++) h_W->Add(h_WJets[i]);
  TH1F * h_TT = (TH1F*)h_TTJets[0]->Clone("h_TT");  for (int i=1;i<h_TTJets.size();i++) h_TT->Add(h_TTJets[i]);
  TH1F * h_T = (TH1F*)h_SnglT[0]->Clone("h_T");  for (int i=1;i<h_SnglT.size();i++) h_T->Add(h_SnglT[i]);

  DrawOverflow(h_W); DrawOverflow(h_TT); DrawOverflow(h_T);
  double W_yields;double TT_yields;double T_yields;

  if (bin=="all"){
    W_yields = h_W->Integral();
    TT_yields = h_TT->Integral();
    T_yields = h_T->Integral();
  }

  if (bin=="bin1"){
    W_yields = h_W->GetBinContent(1);
    TT_yields = h_TT->GetBinContent(1);
    T_yields = h_T->GetBinContent(1);
  }

  if (bin=="bin2"){
    W_yields = h_W->GetBinContent(2);
    TT_yields = h_TT->GetBinContent(2);
    T_yields = h_T->GetBinContent(2);
  }

  if (bin=="bin3"){
    W_yields = h_W->GetBinContent(3)+h_W->GetBinContent(4);
    TT_yields = h_TT->GetBinContent(3)+h_TT->GetBinContent(4);
    T_yields = h_T->GetBinContent(3)+h_T->GetBinContent(4);
  }

  TString yieldsFileName;
  yieldsFileName = whichRegion+", "+regionLabel;
  TCanvas * can_h = new TCanvas(yieldsFileName,yieldsFileName, 50, 50, 1200, 1200);

  double vals[] = {W_yields,TT_yields,T_yields};
  int colors[] = {kAzure+1,kGreen+3,kSpring-5};
  const char *labels[] = {"WJets","TT","SnglT"};
  int nvals = sizeof(vals)/sizeof(vals[0]);
  TString binDef;
  if (bin=="all") binDef = "allMET";
  else if (bin=="bin1")binDef = "MET300to500";
  else if (bin=="bin2")binDef = "MET500to700";
  else if (bin=="bin3")binDef = "MET700up";
  TString pieName = yieldsFileName + ": " + binDef;
  TPie *pie1 = new TPie("pie1", pieName,nvals,vals,colors, labels);
  // pie1->SetTitle("MET [700+] GeV");
  pie1->SetAngularOffset(355.);
  pie1->SetAngle3D(45.);
  pie1->SetHeight(0.04);
  pie1->SetCircle(0.5,0.5,.25);
  pie1->SetLabelsOffset(.03);
  pie1->SetLabelFormat("%perc");
  pie1->Draw("3d");

  TString savename = whichRegion+"_"+regionLabel+"_"+binDef;
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
  if (savePDFs && whichRegion=="singleLept") can_h->SaveAs("ANPlots/pie_1l_"+regionLabel+"_"+binDef+".pdf","PDF");


  delete pie1;delete can_h;
}

void pieChartPhoton(vector<TH1F*> h_GJets, vector<TH1F*> h_QCD,  TString regionLabel, TString bin){
  TH1F * h_G = (TH1F*)h_GJets[0]->Clone("h_G"); for (int i=1;i<h_GJets.size();i++) h_G->Add(h_GJets[i]);
  TH1F * h_Q = (TH1F*)h_QCD[0]->Clone("h_Q");  for (int i=1;i<h_QCD.size();i++) h_Q->Add(h_QCD[i]);

  DrawOverflow(h_G); DrawOverflow(h_Q);
  double G_yields;double Q_yields;

  if (bin=="all"){
    G_yields = h_G->Integral();
    Q_yields = h_Q->Integral();
  }

  if (bin=="bin1"){
    G_yields = h_G->GetBinContent(1);
    Q_yields = h_Q->GetBinContent(1);
  }

  if (bin=="bin2"){
    G_yields = h_G->GetBinContent(2);
    Q_yields = h_Q->GetBinContent(2);
  }

  if (bin=="bin3"){
    G_yields = h_G->GetBinContent(3)+h_G->GetBinContent(4);
    Q_yields = h_Q->GetBinContent(3)+h_Q->GetBinContent(4);
  }

  TString yieldsFileName;
  yieldsFileName = whichRegion+", "+regionLabel;
  TCanvas * can_h = new TCanvas(yieldsFileName,yieldsFileName, 50, 50, 1200, 1200);

  double vals[] = {G_yields,Q_yields};
  int colors[] = {kAzure+1,kGreen+3};
  const char *labels[] = {"GJets","QCD"};
  int nvals = sizeof(vals)/sizeof(vals[0]);
  TString binDef;
  if (bin=="all") binDef = "allMET";
  else if (bin=="bin1")binDef = "MET300to500";
  else if (bin=="bin2")binDef = "MET500to700";
  else if (bin=="bin3")binDef = "MET700up";
  TString pieName = yieldsFileName + ": " + binDef;
  TPie *pie1 = new TPie("pie1", pieName,nvals,vals,colors, labels);
  // pie1->SetTitle("MET [700+] GeV");
  pie1->SetAngularOffset(355.);
  pie1->SetAngle3D(45.);
  pie1->SetHeight(0.04);
  pie1->SetCircle(0.5,0.5,.25);
  pie1->SetLabelsOffset(.03);
  pie1->SetLabelFormat("%perc");
  pie1->Draw("3d");

  TString savename = whichRegion+"_"+regionLabel+"_"+binDef;
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
  if (savePDFs && whichRegion=="photon") can_h->SaveAs("ANPlots/pie_photon_"+regionLabel+"_"+binDef+".pdf","PDF");

  delete pie1;delete can_h;
}

TH1F* make0EventUncSum(vector<TH1F*> dem_histos) {
  //Histograms should be in this order: QCD, TT, WJets,ZJets, SnglT
  // TH1F * h_sum   = (TH1F*)dem_histos[0]->Clone("h_sum");
  TH1F * h_QCD   = (TH1F*)dem_histos[0]->Clone("h_QCD");
  TH1F * h_TT    = (TH1F*)dem_histos[1]->Clone("h_TT");
  TH1F * h_WJets = (TH1F*)dem_histos[2]->Clone("h_WJets");
  TH1F * h_ZJets = (TH1F*)dem_histos[3]->Clone("h_ZJets");
  TH1F * h_SnglT = (TH1F*)dem_histos[4]->Clone("h_SnglT");

  TH1F *h_sum = (TH1F*)h_QCD->Clone("h_sum");
  h_sum->Add(h_TT);
  h_sum->Add(h_WJets);
  h_sum->Add(h_ZJets);
  h_sum->Add(h_SnglT);

  DrawOverflow(h_sum);

  return h_sum;
}

TH1F *make0EventUncSum_1l(vector<TH1F*> dem_histos) {
  //Histograms should be in this order: TT, WJets, SnglT
  TH1F * h_TT    = (TH1F*)dem_histos[0]->Clone("h_TT");
  TH1F * h_WJets = (TH1F*)dem_histos[1]->Clone("h_WJets");
  TH1F * h_SnglT = (TH1F*)dem_histos[2]->Clone("h_SnglT");

  TH1F *h_sum = (TH1F*)h_TT->Clone("h_sum");
  h_sum->Add(h_WJets);
  h_sum->Add(h_SnglT);
  DrawOverflow(h_sum);
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

void DrawOverflow(TH1F* &h){
  int lastBin = h->GetNbinsX();
  double lastBinError = sqrt(h->GetBinError(lastBin)*h->GetBinError(lastBin) + h->GetBinError(lastBin+1)*h->GetBinError(lastBin+1));
  double lastBinContent = h->GetBinContent(lastBin)+h->GetBinContent(lastBin+1);
  h->SetBinContent(lastBin, lastBinContent);
  h->SetBinError(lastBin, lastBinError);
  h->SetBinContent(lastBin+1, 0.0);
  h->SetBinError(lastBin+1, 0.0);
}

void massCorrelations(vector<TH1F*> dem_histos, TString bkgType = ""){
  TH1F *histo_A = (TH1F*)dem_histos[0]->Clone("histo_A"); TH1F *histo_A1 = (TH1F*)dem_histos[1]->Clone("histo_A1");
  TH1F *histo_B = (TH1F*)dem_histos[2]->Clone("histo_B"); TH1F *histo_B1 = (TH1F*)dem_histos[3]->Clone("histo_B1");
  TH1F *histo_C = (TH1F*)dem_histos[4]->Clone("histo_C"); TH1F *histo_D = (TH1F*)dem_histos[5]->Clone("histo_D");

  DrawOverflow(histo_A); DrawOverflow(histo_A1); DrawOverflow(histo_B); DrawOverflow(histo_B1); DrawOverflow(histo_C); DrawOverflow(histo_D);

  //calculate ratios (SIG/SB) and errors
  double a_int; double a_error; a_int = histo_A->IntegralAndError(0,5,a_error,"");
  double a1_int; double a1_error; a1_int = histo_A1->IntegralAndError(0,5,a1_error,"");
  double b_int; double b_error; b_int = histo_B->IntegralAndError(0,5,b_error,"");
  double b1_int; double b1_error; b1_int = histo_B1->IntegralAndError(0,5,b1_error,"");
  double c_int; double c_error; c_int = histo_C->IntegralAndError(0,5,c_error,"");
  double d_int; double d_error; d_int = histo_D->IntegralAndError(0,5,d_error,"");


  float ratio_0H = c_int / d_int;
  float err_0H = ratio_0H * sqrt((c_error/c_int)*(c_error/c_int) + (d_error/d_int)*(d_error/d_int));
  float ratio_1H = a1_int / b1_int;
  float err_1H = ratio_1H * sqrt((a1_error/a1_int)*(a1_error/a1_int) + (b1_error/b1_int)*(b1_error/b1_int));
  float ratio_2H = a_int / b_int;
  float err_2H = ratio_2H * sqrt((a_error/a_int)*(a_error/a_int) + (b_error/b_int)*(b_error/b_int));

  //adding data statistics error bars
  float err_0H_data = ratio_0H * sqrt((sqrt(c_int)/c_int)*(sqrt(c_int)/c_int) + (sqrt(d_int)/d_int)*(sqrt(d_int)/d_int));
  float err_1H_data = ratio_1H * sqrt((sqrt(a1_int)/a1_int)*(sqrt(a1_int)/a1_int) + (sqrt(b1_int)/b1_int)*(sqrt(b1_int)/b1_int));
  float err_2H_data = ratio_2H * sqrt((sqrt(a_int)/a_int)*(sqrt(a_int)/a_int) + (sqrt(b_int)/b_int)*(sqrt(b_int)/b_int));

  TString name = "massCorr"+bkgType;
  TString title = whichRegion+": Mass to bb-tag correlation check;;SIG/SB";
  TH1F * massCorrPlot = new TH1F(name, title,3,0,3);

  massCorrPlot->SetBinContent(1,ratio_0H); massCorrPlot->SetBinContent(2,ratio_1H); massCorrPlot->SetBinContent(3,ratio_2H);
  massCorrPlot->SetBinError(1,err_0H); massCorrPlot->SetBinError(2,err_1H); massCorrPlot->SetBinError(3,err_2H);


  const int n = 3;
  double x[n]   = {0.5,1.5,2.5}; double y[n]   = {0.0,0.0,0.0};
  double exl[n] = {0.5,0.5,0.5}; double exh[n] = {0.5,0.5,0.5};
  double eyl[n] = {0.0,0.0,0.0}; double eyh[n] = {0.0,0.0,0.0};
  y[0] = ratio_0H; y[1] = ratio_1H; y[2] = ratio_2H;
  eyl[0] = err_0H_data; eyl[1] = err_1H_data; eyl[2] = err_2H_data;
  eyh[0] = err_0H_data; eyh[1] = err_1H_data; eyh[2] = err_2H_data;


  TGraphAsymmErrors * graph_errorBands = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

  // graph_errorBands->SetTitle("TGraphAsymmErrors Example");
  graph_errorBands->SetFillColor(kBlue); graph_errorBands->SetFillStyle(3445);
  graph_errorBands->SetMarkerSize(0); graph_errorBands->SetLineWidth(1);
  graph_errorBands->SetLineColor(kBlue);

  massCorrPlot->GetXaxis()->SetBinLabel(1,"0H");
  massCorrPlot->GetXaxis()->SetBinLabel(2,"1H");
  massCorrPlot->GetXaxis()->SetBinLabel(3,"2H");
  massCorrPlot->GetYaxis()->SetTitleSize(0.045);
  massCorrPlot->GetYaxis()->SetLabelSize(0.03);
  // massCorrPlot->GetYaxis()->SetTitleOffset(1.2);
  massCorrPlot->GetXaxis()->SetLabelSize(0.07);
  massCorrPlot->SetStats(0);

  TF1*f0=new TF1("f0","pol0");
  massCorrPlot->Fit("f0","Q");
  double p0 = f0->GetParameter(0); double error = f0->GetParError(0);
  // graph->GetFunction("f0")->SetBit(TF1::kNotDraw);

  // TLine *line = new TLine(60.0,p0,260.0,p0);
  // line->SetLineColor(kRed);
  // line->SetLineWidth(3);
  // line->SetLineStyle(2);
  // line->Draw("same");
  std::string lineConst = to_string(p0); std::string lineErr = to_string(error);
  TString constString = "p0 = "+lineConst+" #pm "+ lineErr;
  TLatex *t = new TLatex(0.2,0.85,constString);
  t->SetNDC(); //t->SetTextFont(52);
  t->SetTextSize(0.05);


  TString massCorrFileName;
  massCorrFileName = "massCorr, "+bkgType;
  TCanvas * can_h = new TCanvas(massCorrFileName,massCorrFileName, 50, 50, 1200, 1200);
  massCorrPlot->SetMarkerStyle(20); massCorrPlot->SetMarkerColor(kBlack); massCorrPlot->SetLineColor(kBlack);
  massCorrPlot->Draw();
  graph_errorBands->Draw("E2");
  massCorrPlot->Draw("same");
  t->Draw("same");
  can_h->Write("massCorr");
  if (savePDFs && whichRegion=="signal") can_h->SaveAs("ANPlots/SIGSBRatio_0l.pdf","PDF");
  else if (savePDFs && whichRegion=="singleLept") can_h->SaveAs("ANPlots/SIGSBRatio_1l.pdf","PDF");
  else if (savePDFs && whichRegion=="photon") can_h->SaveAs("ANPlots/SIGSBRatio_photon.pdf","PDF");
}
