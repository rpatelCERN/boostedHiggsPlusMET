// ROOT/custom libraries
#include "TChain.h"
#include "RA2bTree.cc"
#include "TString.h"

// STL libraries
#include <iostream>
#include <vector>

bool onCUBl = true;
// static const TString BASE_DIR="/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV16/tree_signal";
static const TString BASE_DIR="/mnt/hadoop/store/user/aperloff/SusyRA2Analysis2015/Skims/Run2ProductionV12/tree_signal";


class skimSamples {

  public :

  TChain *WJets,*ZJets,*QCD,*SnglT,*TT,*GJets,*GJets0p4,*Other,*DY,*TTinc;
  //TChain *T5HH750, *T5HH1000, *T5HH1100,*T5HH1200,*T5HH1300,*T5HH1400,*T5HH1500,*T5HH1600,*T5HH1700,*T5HH1800,*T5HH1900,*T5HH2000,*T5HH2100;
  //TChain *T5qqqqZH750, *T5qqqqZH1000, *T5qqqqZH1100, *T5qqqqZH1200, *T5qqqqZH1300, *T5qqqqZH1400, *T5qqqqZH1500, *T5qqqqZH1600, *T5qqqqZH1700, *T5qqqqZH1800, *T5qqqqZH1900, *T5qqqqZH2000, *T5qqqqZH2100, *T5qqqqZH2200;
  TChain *T5qqqqZH1700, *T5qqqqZH2100;
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

  skimSamples(region r=kSignal){
    skimType="";
    if( r == kSignal ){
      skimType=BASE_DIR;
    }

    ///////////////////////////////////////////////////////////////////////
    // - - - - - - - - - - BACKGROUND INPUTS - - - - - - - - - - - - - - //
    ///////////////////////////////////////////////////////////////////////
    /*
    //This version of the ZJets 2016 has a bugged AK8 collection. Can rather use 2017 or 2018
    std::vector<TString> ZJetsFileNames;
    if (onCUBl){
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-100to200.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-200to400.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-400to600.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-600to800.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-800to1200.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-1200to2500.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-2500toInf.root");
    }
    else {
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-100to200_MC2017.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-200to400_MC2017.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-400to600_MC2017.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-600to800_MC2017.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-800to1200_MC2017.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-1200to2500_MC2017.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-2500toInf_MC2017.root");
    }



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
    if (onCUBl){
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-100to200.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-200to400.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-400to600.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-600to800.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-800to1200.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-1200to2500.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-2500toInf.root");
    }
    else {
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-100to200_MC2017.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-200to400_MC2017.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-400to600_MC2017.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-600to800_MC2017.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-800to1200_MC2017.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-1200to2500_MC2017.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-2500toInf_MC2017.root");
    }


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
    if (onCUBl){
      TTFileNames.push_back("tree_TTJets_HT-600to800.root");
      TTFileNames.push_back("tree_TTJets_HT-800to1200.root");
      TTFileNames.push_back("tree_TTJets_HT-1200to2500.root");
      TTFileNames.push_back("tree_TTJets_HT-2500toInf.root");
      TTFileNames.push_back("tree_TTJets_SingleLeptFromT.root");
      TTFileNames.push_back("tree_TTJets_SingleLeptFromTbar.root");
      TTFileNames.push_back("tree_TTJets_DiLept.root");
      // TTFileNames.push_back("tree_TTTT.root");
      // TTFileNames.push_back("tree_TTWJetsToLNu.root");
      // TTFileNames.push_back("tree_TTWJetsToQQ.root");
      // TTFileNames.push_back("tree_TTZToLLNuNu.root");
      // TTFileNames.push_back("tree_TTZToQQ.root");
      // TTFileNames.push_back("tree_TTGJets.root");
    }
    else {
      TTFileNames.push_back("tree_TTJets_HT-600to800_MC2017.root");
      TTFileNames.push_back("tree_TTJets_HT-800to1200_MC2017.root");
      TTFileNames.push_back("tree_TTJets_HT-1200to2500_MC2017.root");
      TTFileNames.push_back("tree_TTJets_HT-2500toInf_MC2017.root");
      //TTFileNames.push_back("tree_TTJets_MC2017.root"); //for resolved
      TTFileNames.push_back("tree_TTJets_SingleLeptFromT_MC2017.root");
      TTFileNames.push_back("tree_TTJets_SingleLeptFromTbar_MC2017.root");
      TTFileNames.push_back("tree_TTJets_DiLept_MC2017.root");
      // TTFileNames.push_back("tree_TTTT_MC2017.root");
      // TTFileNames.push_back("tree_TTWJetsToLNu_MC2017.root");
      // TTFileNames.push_back("tree_TTWJetsToQQ_MC2017.root");
      // TTFileNames.push_back("tree_TTZToLLNuNu_MC2017.root");
      // TTFileNames.push_back("tree_TTZToQQ_MC2017.root");
      // TTFileNames.push_back("tree_TTGJets_MC2017.root");
    }


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
    */

    std::vector<TString> QCDFileNames;
    if (onCUBl){
      QCDFileNames.push_back("tree_QCD_HT-200to300.root");
      QCDFileNames.push_back("tree_QCD_HT-300to500.root");
      QCDFileNames.push_back("tree_QCD_HT-500to700.root");
      QCDFileNames.push_back("tree_QCD_HT-700to1000.root");
      QCDFileNames.push_back("tree_QCD_HT-1000to1500.root");
      QCDFileNames.push_back("tree_QCD_HT-1500to2000.root");
      QCDFileNames.push_back("tree_QCD_HT-2000toInf.root");
    }
    else {
      QCDFileNames.push_back("tree_QCD_HT-200to300_MC2017.root");
      QCDFileNames.push_back("tree_QCD_HT-300to500_MC2017.root");
      QCDFileNames.push_back("tree_QCD_HT-500to700_MC2017.root");
      QCDFileNames.push_back("tree_QCD_HT-700to1000_MC2017.root");
      QCDFileNames.push_back("tree_QCD_HT-1000to1500_MC2017.root");
      QCDFileNames.push_back("tree_QCD_HT-1500to2000_MC2017.root");
      QCDFileNames.push_back("tree_QCD_HT-2000toInf_MC2017.root");
    }

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
    if( r == kSignal && !onCUBl){
      T5qqqqZH1700 = new TChain("tree");
      T5qqqqZH2100 = new TChain("tree");
      T5qqqqZH1700->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/T5ZH/RA2ProductionV16/tree_T5qqqqZH_1700_MC2016.root");
      T5qqqqZH2100->Add("root://cmseos.fnal.gov//store/user/rgp230/SUSY/T5ZH/RA2ProductionV16/tree_T5qqqqZH_2100_MC2016.root");

      ntuples.push_back(new RA2bTree(T5qqqqZH1700));
      ntuples.push_back(new RA2bTree(T5qqqqZH2100));

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
