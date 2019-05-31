// ROOT/custom libraries
#include "TChain.h"
#include "RA2bTree.cc"
#include "TString.h"

// STL libraries
#include <iostream>
#include <vector>

static const TString BASE_DIR="/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV17/";
// static const TString BASE_DIR="/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV16/tree_signal";
// static const TString BASE_DIR="/mnt/hadoop/store/user/aperloff/SusyRA2Analysis2015/Skims/Run2ProductionV12/tree_signal";


class skimSamples {

  public :

  TChain *WJets,*ZJets,*QCD,*SnglT,*TT,*GJets,*GJets0p4,*Other,*DY,*TTinc;
  //TChain *T5HH750, *T5HH1000, *T5HH1100,*T5HH1200,*T5HH1300,*T5HH1400,*T5HH1500,*T5HH1600,*T5HH1700,*T5HH1800,*T5HH1900,*T5HH2000,*T5HH2100;
  //TChain *T5qqqqZH750, *T5qqqqZH1000, *T5qqqqZH1100, *T5qqqqZH1200, *T5qqqqZH1300, *T5qqqqZH1400, *T5qqqqZH1500, *T5qqqqZH1600, *T5qqqqZH1700, *T5qqqqZH1800, *T5qqqqZH1900, *T5qqqqZH2000, *T5qqqqZH2100, *T5qqqqZH2200;
  TChain *T5qqqqZH1300, *T5qqqqZH1700, *T5qqqqZH2100;
  //TChain *TChiHH127, *TChiHH150, *TChiHH175,*TChiHH200,*TChiHH225, *TChiHH250, *TChiHH275, *TChiHH300, *TChiHH325,*TChiHH350, *TChiHH375, *TChiHH400, *TChiHH425, *TChiHH450,*TChiHH475,*TChiHH500, *TChiHH525, *TChiHH550,*TChiHH575, *TChiHH600, *TChiHH625, *TChiHH650, *TChiHH675, *TChiHH700,*TChiHH725, *TChiHH750,*TChiHH775, *TChiHH800, *TChiHH825, *TChiHH850, *TChiHH875, *TChiHH900, *TChiHH925, *TChiHH950, *TChiHH975, *TChiHH1000;
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
        skimType=BASE_DIR+"tree_GJet/";
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
      ntuples.push_back(new RA2bTree(ZJets));//Need this one
      sampleName.push_back("ZJets");
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
      ntuples.push_back(new RA2bTree(WJets));//Need this one
      sampleName.push_back("WJets");
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



    TT = new TChain("tree");
    for( int i = 0 ; i < TTFileNames.size() ; i++ ){
      TT->Add(skimType+"/"+TTFileNames[i]);
    }
    if( r == kSignal || r == kSLm || r == kSLe || r == kLowDphi ){
      ntuples.push_back(new RA2bTree(TT));//Need this one
      sampleName.push_back("TT");
      fillColor.push_back(kCyan);
      lineColor.push_back(kCyan);
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
      ntuples.push_back(new RA2bTree(QCD));//Need this one
      sampleName.push_back("QCD");
      fillColor.push_back(kGray);
      lineColor.push_back(1);
    }

    /*
    T5qqqqZH1300 = new TChain("tree");
    T5qqqqZH1700 = new TChain("tree");
    T5qqqqZH2100 = new TChain("tree");
    if( r == kSignal ){
      T5qqqqZH1300->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/BooostedZSkims/tree_signal/tree_T5qqqqZH_1300_"+Year+".root");
      T5qqqqZH1700->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/BooostedZSkims/tree_signal/tree_T5qqqqZH_1700_"+Year+".root");
      T5qqqqZH2100->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/BooostedZSkims/tree_signal/tree_T5qqqqZH_2100_"+Year+".root");
      // T5qqqqZH2100->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/T5ZH/RA2ProductionV16/tree_T5qqqqZH_2100_"+Year+".root");

      ntuples.push_back(new RA2bTree(T5qqqqZH1300));
      ntuples.push_back(new RA2bTree(T5qqqqZH1700));
      ntuples.push_back(new RA2bTree(T5qqqqZH2100));

      sampleName.push_back("T5qqqqZH1300");
      sampleName.push_back("T5qqqqZH1700");
      sampleName.push_back("T5qqqqZH2100");
      for(unsigned int i=0; i<sampleName.size(); ++i)sigLineColor.push_back(kRed);
    }
    */

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
