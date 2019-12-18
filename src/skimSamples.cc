// ROOT/custom libraries
#include "TChain.h"
#include "RA2bTree.cc"
#include "TString.h"

// STL libraries
#include <iostream>
#include <vector>

//static const TString BASE_DIR="/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV17/";
static const TString BASE_DIR="root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/Run2ProductionV18a/";
// static const TString BASE_DIR="/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV16/tree_signal";
// static const TString BASE_DIR="/mnt/hadoop/store/user/aperloff/SusyRA2Analysis2015/Skims/Run2ProductionV12/tree_signal";


class skimSamples {

  public :

  TChain *WJets,*ZJets,*QCD,*SnglT,*TT,*GJets,*GJets0p4,*Other,*DY,*TTinc;
  //TChain *T5HH750, *T5HH1000, *T5HH1100,*T5HH1200,*T5HH1300,*T5HH1400,*T5HH1500,*T5HH1600,*T5HH1700,*T5HH1800,*T5HH1900,*T5HH2000,*T5HH2100;
  //TChain *T5qqqqZH750, *T5qqqqZH1000, *T5qqqqZH1100, *T5qqqqZH1200, *T5qqqqZH1300, *T5qqqqZH1400, *T5qqqqZH1500, *T5qqqqZH1600, *T5qqqqZH1700, *T5qqqqZH1800, *T5qqqqZH1900, *T5qqqqZH2000, *T5qqqqZH2100, *T5qqqqZH2200;
  TChain *T5qqqqZH1300, *T5qqqqZH1700, *T5qqqqZH2100;

//TChain*TChiHH225,*TChiHH400,*TChiHH700;
  TChain *TChiHH127, *TChiHH150, *TChiHH175,*TChiHH200,*TChiHH225, *TChiHH250, *TChiHH275, *TChiHH300, *TChiHH325,*TChiHH350, *TChiHH375, *TChiHH400, *TChiHH425, *TChiHH450,*TChiHH475,*TChiHH500, *TChiHH525, *TChiHH550,*TChiHH575, *TChiHH600, *TChiHH625, *TChiHH650, *TChiHH675, *TChiHH700,*TChiHH725, *TChiHH750,*TChiHH775, *TChiHH800, *TChiHH825, *TChiHH850, *TChiHH875, *TChiHH900, *TChiHH925, *TChiHH950, *TChiHH975, *TChiHH1000;
  TChain*TChiHH1025, *TChiHH1050,*TChiHH1075,*TChiHH1100,*TChiHH1125, *TChiHH1150,*TChiHH1175,*TChiHH1200,*TChiHH1225, *TChiHH1250,*TChiHH1275,*TChiHH1300,*TChiHH1325, *TChiHH1350,*TChiHH1375,*TChiHH1400,*TChiHH1425, *TChiHH1450,*TChiHH1475,*TChiHH1500;
  TChain *data;
  std::vector<RA2bTree*> ntuples,signalNtuples;
  RA2bTree* dataNtuple;
  std::vector<TString> sampleName, signalSampleName;
  std::vector<TString> dataSampleName;
  std::vector<int> fillColor, lineColor, sigLineColor;

  enum region {kSignal,kPhoton,kSLm,kSLe,kLowDphi, kNumRegions};
  TString regionNames[kNumRegions]={"signal","photon","SLm","SLe","kLowDphi"};

  TString skimType;

  skimSamples(region r=kSignal, TString Year="MC2016"){
    skimType="";

    if( r == kSignal ){
        skimType=BASE_DIR+"tree_signal/";
    }
    if( r == kPhoton ){
        skimType=BASE_DIR+"tree_GJet_CleanVars/";
    }
    if( r == kSLm ){
        skimType=BASE_DIR+"tree_SLm/";
    }
    if( r == kSLe ){
        skimType=BASE_DIR+"tree_SLe/";
    }
    if(r==kLowDphi){
        skimType=BASE_DIR+"tree_LDP/";
    }

    ///////////////////////////////////////////////////////////////////////
    // - - - - - - - - - - BACKGROUND INPUTS - - - - - - - - - - - - - - //
    ///////////////////////////////////////////////////////////////////////

    //This version of the ZJets 2016 has a bugged AK8 collection. Can rather use 2017 or 2018
    std::vector<TString> ZJetsFileNames;

    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-100to200_"+Year+".root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-200to400_"+Year+".root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-400to600_"+Year+".root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-600to800_"+Year+".root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-800to1200_"+Year+".root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-1200to2500_"+Year+".root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-2500toInf_"+Year+".root");

    ZJets = new TChain("tree");
    for( int i = 0 ; i < ZJetsFileNames.size() ; i++ ){
      ZJets->Add(skimType+"/"+ZJetsFileNames[i]);
    }
    if( r == kSignal || r == kLowDphi ){
     // ntuples.push_back(new RA2bTree(ZJets));//Need this one
     // sampleName.push_back("ZJets");
      fillColor.push_back(kGreen+1);
      lineColor.push_back(1);
    }

    std::vector<TString> WJetsFileNames;
    WJetsFileNames.push_back("tree_WJetsToLNu_HT-100to200_"+Year+".root");
    WJetsFileNames.push_back("tree_WJetsToLNu_HT-200to400_"+Year+".root");
    WJetsFileNames.push_back("tree_WJetsToLNu_HT-400to600_"+Year+".root");
    WJetsFileNames.push_back("tree_WJetsToLNu_HT-600to800_"+Year+".root");
    WJetsFileNames.push_back("tree_WJetsToLNu_HT-800to1200_"+Year+".root");
    WJetsFileNames.push_back("tree_WJetsToLNu_HT-1200to2500_"+Year+".root");
    WJetsFileNames.push_back("tree_WJetsToLNu_HT-2500toInf_"+Year+".root");

    WJets = new TChain("tree");
    for( int i = 0 ; i < WJetsFileNames.size() ; i++ ){
      WJets->Add(skimType+"/"+WJetsFileNames[i]);
    }
    if( r == kSignal || r == kSLm || r == kSLe || r == kLowDphi ){
  //    ntuples.push_back(new RA2bTree(WJets));//Need this one
   //   sampleName.push_back("WJets");
      fillColor.push_back(kBlue);
      lineColor.push_back(1);
    }


    std::vector<TString> TTincFileNames;
    TTincFileNames.push_back("tree_TTJets.root");
    TTinc = new TChain("tree");
    for( int i = 0 ; i < TTincFileNames.size() ; i++ ){
      TTinc->Add(skimType+"/"+TTincFileNames[i]);
    }

    std::vector<TString> TTFileNames;

    //Use only these for boosted ////////////////////////////////
    TTFileNames.push_back("tree_TTJets_HT-600to800_"+Year+".root");
    TTFileNames.push_back("tree_TTJets_HT-800to1200_"+Year+".root");
    TTFileNames.push_back("tree_TTJets_HT-1200to2500_"+Year+".root");
    TTFileNames.push_back("tree_TTJets_HT-2500toInf_"+Year+".root");
    TTFileNames.push_back("tree_TTJets_SingleLeptFromT_"+Year+".root");
    TTFileNames.push_back("tree_TTJets_SingleLeptFromTbar_"+Year+".root");
    TTFileNames.push_back("tree_TTJets_DiLept_"+Year+".root");

    ////////////////////////////////////////////////////////////////

    // TTFileNames.push_back("tree_TTTT.root");
    // TTFileNames.push_back("tree_TTWJetsToLNu.root");
    // TTFileNames.push_back("tree_TTWJetsToQQ.root");
    // TTFileNames.push_back("tree_TTZToLLNuNu.root");
    // TTFileNames.push_back("tree_TTZToQQ.root");
    // TTFileNames.push_back("tree_TTGJets.root"); //doesn't exist for 2016
/*
    std::vector<TString> SnglTFileNames;
   SnglTFileNames.push_back("tree_ST_s-channel.root");
    SnglTFileNames.push_back("tree_ST_t-channel_antitop.root");
    SnglTFileNames.push_back("tree_ST_t-channel_top.root");
    SnglTFileNames.push_back("tree_ST_tW_antitop.root");
    SnglTFileNames.push_back("tree_ST_tW_top.root");
    SnglT = new TChain("tree");
    for( int i = 0 ; i < SnglTFileNames.size() ; i++ ) {
      SnglT->Add(skimType+"/"+SnglTFileNames[i]);
    }
    if( r == kSignal || r == kSLm || r == kSLe ){
      ntuples.push_back(new RA2bTree(SnglT));
      sampleName.push_back("SnglT");
      fillColor.push_back(kOrange);
      lineColor.push_back(1);
    }
*/

    TT = new TChain("tree");
    for( int i = 0 ; i < TTFileNames.size() ; i++ ){
      TT->Add(skimType+"/"+TTFileNames[i]);
    }
    if( r == kSignal || r == kSLm || r == kSLe || r == kLowDphi ){
      //ntuples.push_back(new RA2bTree(TT));//Need this one
      //sampleName.push_back("TT");
      fillColor.push_back(kCyan);
      lineColor.push_back(kCyan);
    }


       std::vector<TString> GJets0p4FileNames;
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-100to200_");
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-200to400_");
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-400to600_");
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-600toInf_");
        GJets0p4 = new TChain("tree");
        for( int i = 0 ; i < GJets0p4FileNames.size() ; i++ ){
            GJets0p4->Add(skimType+"/"+GJets0p4FileNames[i]+Year+".root");
        }
        if( r == kPhoton ){
            ntuples.push_back(new RA2bTree(GJets0p4));
            sampleName.push_back("GJets");
            fillColor.push_back(kGreen);
            lineColor.push_back(1);
        }
    std::vector<TString> QCDFileNames;
    QCDFileNames.push_back("tree_QCD_HT-200to300_"+Year+".root");
    QCDFileNames.push_back("tree_QCD_HT-300to500_"+Year+".root");
    QCDFileNames.push_back("tree_QCD_HT-500to700_"+Year+".root");
    QCDFileNames.push_back("tree_QCD_HT-700to1000_"+Year+".root");
    QCDFileNames.push_back("tree_QCD_HT-1000to1500_"+Year+".root");
    QCDFileNames.push_back("tree_QCD_HT-1500to2000_"+Year+".root");
    QCDFileNames.push_back("tree_QCD_HT-2000toInf_"+Year+".root");


    QCD = new TChain("tree");
    for( int i = 0 ; i < QCDFileNames.size() ; i++ ){
      QCD->Add(skimType+"/"+QCDFileNames[i]);
    }
    if( r == kSignal || r == kPhoton || r == kLowDphi ){
     //ntuples.push_back(new RA2bTree(QCD));//Need this one
     // sampleName.push_back("QCD");
      fillColor.push_back(kGray);
      lineColor.push_back(1);
    }
      
    TChiHH127 = new TChain("tree");
    TChiHH150 = new TChain("tree");
    TChiHH175 = new TChain("tree");
    TChiHH200 = new TChain("tree");
    TChiHH225 = new TChain("tree");
    TChiHH250 = new TChain("tree");
    TChiHH275 = new TChain("tree");
    TChiHH300 = new TChain("tree");
    TChiHH325 = new TChain("tree");
    TChiHH350 = new TChain("tree");
    TChiHH375 = new TChain("tree");
    TChiHH400 = new TChain("tree");
    TChiHH425 = new TChain("tree");
    TChiHH450 = new TChain("tree");
    TChiHH475 = new TChain("tree");
    TChiHH500 = new TChain("tree");
    TChiHH525 = new TChain("tree");
    TChiHH550 = new TChain("tree");
    TChiHH575 = new TChain("tree");
    TChiHH600 = new TChain("tree");
    TChiHH625 = new TChain("tree");
    TChiHH650 = new TChain("tree");
    TChiHH675 = new TChain("tree");
    TChiHH700 = new TChain("tree");
    TChiHH725 = new TChain("tree");
    TChiHH750 = new TChain("tree");
    TChiHH775 = new TChain("tree");
    TChiHH800 = new TChain("tree");
    TChiHH825 = new TChain("tree");
    TChiHH850 = new TChain("tree");
    TChiHH875 = new TChain("tree");
    TChiHH900 = new TChain("tree");
    TChiHH925 = new TChain("tree");
    TChiHH950 = new TChain("tree");
    TChiHH975 = new TChain("tree");
    TChiHH1000 = new TChain("tree");
    TChiHH1025 = new TChain("tree");
    TChiHH1050 = new TChain("tree");
    TChiHH1075 = new TChain("tree");
    TChiHH1100 = new TChain("tree");
    TChiHH1125 = new TChain("tree");
    TChiHH1150 = new TChain("tree");
    TChiHH1175 = new TChain("tree");
    TChiHH1200 = new TChain("tree");
    TChiHH1225 = new TChain("tree");
    TChiHH1250 = new TChain("tree");
    TChiHH1275 = new TChain("tree");
    TChiHH1200 = new TChain("tree");
    TChiHH1225 = new TChain("tree");
    TChiHH1250 = new TChain("tree");
    TChiHH1275 = new TChain("tree");
    TChiHH1300 = new TChain("tree");
    TChiHH1325 = new TChain("tree");
    TChiHH1350 = new TChain("tree");
    TChiHH1375 = new TChain("tree");
    TChiHH1400 = new TChain("tree");
    TChiHH1425 = new TChain("tree");
    TChiHH1450 = new TChain("tree");
    TChiHH1475 = new TChain("tree");
    TChiHH1500 = new TChain("tree");
    if( r == kSignal ){
    //  T5qqqqZH1300->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoHiggsEvents/tree_T5qqqqZH_1300_1_"+Year+"_fast.root");
     // T5qqqqZH1700->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoHiggsEvents/tree_T5qqqqZH_1700_1_"+Year+"_fast.root");
     // T5qqqqZH2100->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/TwoHiggsEvents/tree_T5qqqqZH_2100_1_"+Year+"_fast.root");
      TChiHH127->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_127_1_MC2016_fast.root");
      TChiHH150->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_150_1_MC2016_fast.root");
      TChiHH175->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_175_1_MC2016_fast.root");
      TChiHH200->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_200_1_MC2016_fast.root");
      TChiHH225->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_225_1_MC2016_fast.root");
      TChiHH250->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_250_1_MC2016_fast.root");
      TChiHH275->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_275_1_MC2016_fast.root");
      TChiHH300->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_300_1_MC2016_fast.root");
      TChiHH325->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_325_1_MC2016_fast.root");
      TChiHH350->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_350_1_MC2016_fast.root");
      TChiHH375->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_375_1_MC2016_fast.root");
      TChiHH400->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_400_1_MC2016_fast.root");
      TChiHH425->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_425_1_MC2016_fast.root");
      TChiHH450->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_450_1_MC2016_fast.root");
      TChiHH475->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_475_1_MC2016_fast.root");
      TChiHH500->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_500_1_MC2016_fast.root");
      TChiHH525->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_525_1_MC2016_fast.root");
      TChiHH550->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_550_1_MC2016_fast.root");
      TChiHH575->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_575_1_MC2016_fast.root");
      TChiHH600->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_600_1_MC2016_fast.root");
      TChiHH625->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_625_1_MC2016_fast.root");
      TChiHH650->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_650_1_MC2016_fast.root");
      TChiHH675->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_675_1_MC2016_fast.root");
      TChiHH700->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_700_1_MC2016_fast.root");
      TChiHH725->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_725_1_MC2016_fast.root");
      TChiHH750->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_750_1_MC2016_fast.root");
      TChiHH775->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_775_1_MC2016_fast.root");
      TChiHH800->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_800_1_MC2016_fast.root");
      TChiHH825->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_825_1_MC2016_fast.root");
      TChiHH850->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_850_1_MC2016_fast.root");
      TChiHH875->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_875_1_MC2016_fast.root");
      TChiHH900->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_900_1_MC2016_fast.root");
      TChiHH925->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_925_1_MC2016_fast.root");
      TChiHH950->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_950_1_MC2016_fast.root");
      TChiHH975->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_975_1_MC2016_fast.root");
      TChiHH1000->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1000_1_MC2016_fast.root");

      TChiHH1025->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1025_1_MC2016_fast.root");
      TChiHH1050->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1050_1_MC2016_fast.root");
      TChiHH1075->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1075_1_MC2016_fast.root");
      TChiHH1100->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1100_1_MC2016_fast.root");
      TChiHH1125->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1125_1_MC2016_fast.root");
      TChiHH1150->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1150_1_MC2016_fast.root");
      TChiHH1175->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1175_1_MC2016_fast.root");
      TChiHH1200->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1200_1_MC2016_fast.root");
      TChiHH1225->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1225_1_MC2016_fast.root");
      TChiHH1250->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1250_1_MC2016_fast.root");
      TChiHH1275->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1275_1_MC2016_fast.root");
      TChiHH1300->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1300_1_MC2016_fast.root");
      TChiHH1325->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1325_1_MC2016_fast.root");
      TChiHH1350->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1350_1_MC2016_fast.root");
      TChiHH1375->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1375_1_MC2016_fast.root");
      TChiHH1400->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1400_1_MC2016_fast.root");
      TChiHH1425->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1425_1_MC2016_fast.root");
      TChiHH1450->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1450_1_MC2016_fast.root");
      TChiHH1475->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1475_1_MC2016_fast.root");
      TChiHH1500->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_1500_1_MC2016_fast.root");
 //     TChiHH700->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/tree_TChiHH_HToBB_HToBB_700_1_MC2016_fast.root");
      //ntuples.push_back(new RA2bTree(T5qqqqZH1300));
      //ntuples.push_back(new RA2bTree(T5qqqqZH1700));
      //ntuples.push_back(new RA2bTree(T5qqqqZH2100));
      ntuples.push_back(new RA2bTree(TChiHH127));
      ntuples.push_back(new RA2bTree(TChiHH150));
      ntuples.push_back(new RA2bTree(TChiHH175));
      ntuples.push_back(new RA2bTree(TChiHH200));
      ntuples.push_back(new RA2bTree(TChiHH225));
      ntuples.push_back(new RA2bTree(TChiHH250));
      ntuples.push_back(new RA2bTree(TChiHH275));
      ntuples.push_back(new RA2bTree(TChiHH300));
      ntuples.push_back(new RA2bTree(TChiHH325));
      ntuples.push_back(new RA2bTree(TChiHH350));
      ntuples.push_back(new RA2bTree(TChiHH375));
      ntuples.push_back(new RA2bTree(TChiHH400));
      ntuples.push_back(new RA2bTree(TChiHH425));
      ntuples.push_back(new RA2bTree(TChiHH450));
      ntuples.push_back(new RA2bTree(TChiHH475));
      ntuples.push_back(new RA2bTree(TChiHH500));
      ntuples.push_back(new RA2bTree(TChiHH525));
      ntuples.push_back(new RA2bTree(TChiHH550));
      ntuples.push_back(new RA2bTree(TChiHH575));
      ntuples.push_back(new RA2bTree(TChiHH600));
      ntuples.push_back(new RA2bTree(TChiHH625));
      ntuples.push_back(new RA2bTree(TChiHH650));
      ntuples.push_back(new RA2bTree(TChiHH675));
      ntuples.push_back(new RA2bTree(TChiHH700));
      ntuples.push_back(new RA2bTree(TChiHH725));
      ntuples.push_back(new RA2bTree(TChiHH750));
      ntuples.push_back(new RA2bTree(TChiHH775));
      ntuples.push_back(new RA2bTree(TChiHH800));
      ntuples.push_back(new RA2bTree(TChiHH825));
      ntuples.push_back(new RA2bTree(TChiHH850));
      ntuples.push_back(new RA2bTree(TChiHH875));
      ntuples.push_back(new RA2bTree(TChiHH900));
      ntuples.push_back(new RA2bTree(TChiHH925));
      ntuples.push_back(new RA2bTree(TChiHH950));
      ntuples.push_back(new RA2bTree(TChiHH975));
      ntuples.push_back(new RA2bTree(TChiHH1000));
      ntuples.push_back(new RA2bTree(TChiHH1025));
      ntuples.push_back(new RA2bTree(TChiHH1050));
      ntuples.push_back(new RA2bTree(TChiHH1075));
      ntuples.push_back(new RA2bTree(TChiHH1100));
      ntuples.push_back(new RA2bTree(TChiHH1125));
      ntuples.push_back(new RA2bTree(TChiHH1150));
      ntuples.push_back(new RA2bTree(TChiHH1175));
      ntuples.push_back(new RA2bTree(TChiHH1200));
      ntuples.push_back(new RA2bTree(TChiHH1225));
      ntuples.push_back(new RA2bTree(TChiHH1250));
      ntuples.push_back(new RA2bTree(TChiHH1275));
      ntuples.push_back(new RA2bTree(TChiHH1300));
      ntuples.push_back(new RA2bTree(TChiHH1325));
      ntuples.push_back(new RA2bTree(TChiHH1350));
      ntuples.push_back(new RA2bTree(TChiHH1375));
      ntuples.push_back(new RA2bTree(TChiHH1400));
      ntuples.push_back(new RA2bTree(TChiHH1425));
      ntuples.push_back(new RA2bTree(TChiHH1450));
      ntuples.push_back(new RA2bTree(TChiHH1475));
      ntuples.push_back(new RA2bTree(TChiHH1500));
      sampleName.push_back("TChiHH127");
     for(unsigned int i=0; i<55; ++i)sampleName.push_back(TString::Format("TChiHH%d", 150+i*25));
       for(unsigned int i=0; i<sampleName.size(); ++i)sigLineColor.push_back(kRed);
    }
  }; //end Skim samples

  RA2bTree* findNtuple(TString name){
    for( int iSam = 0 ; iSam < sampleName.size() ; iSam++ ){
      if( sampleName[iSam] == name )
      return ntuples[iSam] ;
    }
    for( int iSam = 0 ; iSam < signalSampleName.size() ; iSam++ ){
      if( signalSampleName[iSam] == name )
      return signalNtuples[iSam] ;
    }
    return NULL;
  };

};
