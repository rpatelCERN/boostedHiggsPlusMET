// ROOT/custom libraries
#include "TChain.h"
#include "RA2bTree.cc"
#include "TString.h"

// STL libraries
#include <iostream>
#include <vector>

// static const TString BASE_DIR="root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV12/";
static const TString BASE_DIR="/mnt/hadoop/store/user/aperloff/SusyRA2Analysis2015/Skims/Run2ProductionV12/";
static const TString emily_Dir_Bkg="/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/background/";
static const TString emily_Dir_TChiSignal="/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/";

class skimSamples{

  public :

  TChain *WJets,*ZJets,*QCD,*SnglT,*TT,*GJets,*GJets0p4,*Other,*DY,*TTinc;
  TChain *T5HH750, *T5HH1000, *T5HH1100,*T5HH1200,*T5HH1300,*T5HH1400,*T5HH1500,*T5HH1600,*T5HH1700,*T5HH1800,*T5HH1900,*T5HH2000,*T5HH2100,*T5HH2200;
  TChain *TChiHH127, *TChiHH150, *TChiHH175,*TChiHH200,*TChiHH225, *TChiHH250, *TChiHH275, *TChiHH300, *TChiHH325,*TChiHH350, *TChiHH375, *TChiHH400, *TChiHH425, *TChiHH450,*TChiHH475,*TChiHH500, *TChiHH525, *TChiHH550,*TChiHH575, *TChiHH600, *TChiHH625, *TChiHH650, *TChiHH675, *TChiHH700,*TChiHH725, *TChiHH750,*TChiHH775, *TChiHH800, *TChiHH825, *TChiHH850, *TChiHH875, *TChiHH900, *TChiHH925, *TChiHH950, *TChiHH975, *TChiHH1000;
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
      // skimType=BASE_DIR+"tree_signal/";
      skimType=emily_Dir_Bkg;
    }
    if( r == kPhoton ){
      // skimType="root://cmseos.fnal.gov//store/user/fojensen/boostedSkims_19062017/Run2ProductionV12/tree_GJet/";
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

    //bools that determine which processes are run, as a subset of signal
    bool run_TT = true;
    bool run_QCD = true;
    bool run_WJets = true;
    bool run_ZJets = true;
    bool run_AllTChiHH = true;
    bool run_SomeTChiHH = false;
    bool run_T5HH = true;


    ///////////////////////////////////////////////////////////////////////
    // - - - - - - - - - - BACKGROUND INPUTS - - - - - - - - - - - - - - //
    ///////////////////////////////////////////////////////////////////////

    /*
    std::vector<TString> OtherFileNames;
    OtherFileNames.push_back("tree_WWTo1L1Nu2Q.root");
    OtherFileNames.push_back("tree_WWTo2L2Nu.root");
    OtherFileNames.push_back("tree_WWZ.root");
    OtherFileNames.push_back("tree_WZTo1L1Nu2Q.root");
    OtherFileNames.push_back("tree_WZTo1L3Nu.root");
    OtherFileNames.push_back("tree_WZZ.root");
    OtherFileNames.push_back("tree_ZZTo2L2Q.root");
    OtherFileNames.push_back("tree_ZZTo2Q2Nu.root");
    OtherFileNames.push_back("tree_ZZZ.root");
    OtherFileNames.push_back("tree_TTTT.root");
    OtherFileNames.push_back("tree_TTWJetsToLNu.root");
    OtherFileNames.push_back("tree_TTWJetsToQQ.root");
    OtherFileNames.push_back("tree_TTGJets.root");
    OtherFileNames.push_back("tree_TTZToLLNuNu.root");
    OtherFileNames.push_back("tree_TTZToQQ.root");
    Other = new TChain("tree");
    for( int i = 0 ; i < OtherFileNames.size() ; i++ ){
      Other->Add(skimType+"/"+OtherFileNames[i]);
    }
    if( r == kSignal || r == kSLm || r == kSLe || r == kLowDphi || r == kPhoton ){
      ntuples.push_back(new RA2bTree(Other));
      sampleName.push_back("Other");
      fillColor.push_back(kRed+1);
      lineColor.push_back(1);
    }
    */

    if (run_ZJets){
      std::vector<TString> ZJetsFileNames;
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-100to200.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-200to400.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-400to600.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-600to800.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-800to1200.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-1200to2500.root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-2500toInf.root");
      ZJets = new TChain("tree");
      for( int i = 0 ; i < ZJetsFileNames.size() ; i++ ){
        ZJets->Add(skimType+"/"+ZJetsFileNames[i]);
      }
      if( r == kSignal || r == kLowDphi ){
        ntuples.push_back(new RA2bTree(ZJets));
        sampleName.push_back("ZJets");
        fillColor.push_back(kGreen+1);
        lineColor.push_back(1);
      }
    }

    if (run_WJets){
      std::vector<TString> WJetsFileNames;
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-100to200.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-1200to2500.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-200to400.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-2500toInf.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-400to600.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-600to800.root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-800to1200.root");
      WJets = new TChain("tree");
      for( int i = 0 ; i < WJetsFileNames.size() ; i++ ){
        WJets->Add(skimType+"/"+WJetsFileNames[i]);
      }
      if( r == kSignal || r == kSLm || r == kSLe || r == kLowDphi ){
        ntuples.push_back(new RA2bTree(WJets));
        sampleName.push_back("WJets");
        fillColor.push_back(kBlue);
        lineColor.push_back(1);
      }
    }
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

    std::vector<TString> TTincFileNames;
    TTincFileNames.push_back("tree_TTJets.root");
    TTinc = new TChain("tree");
    for( int i = 0 ; i < TTincFileNames.size() ; i++ ){
      TTinc->Add(skimType+"/"+TTincFileNames[i]);
    }
    */
    // if( r == kSignal || r == kSLm || r == kSLe || r == kLowDphi ){
    //   ntuples.push_back(new RA2bTree(TTinc));
    //   sampleName.push_back("TT");
    //   fillColor.push_back(kCyan);
    //   lineColor.push_back(kCyan);
    // }

    if (run_TT){
      std::vector<TString> TTFileNames;
      TTFileNames.push_back("tree_TTJets_HT-600to800.root");
      TTFileNames.push_back("tree_TTJets_HT-800to1200.root");
      TTFileNames.push_back("tree_TTJets_HT-1200to2500.root");
      TTFileNames.push_back("tree_TTJets_HT-2500toInf.root");
      TTFileNames.push_back("tree_TTJets_SingleLeptFromT.root");
      TTFileNames.push_back("tree_TTJets_SingleLeptFromTbar.root");
      TTFileNames.push_back("tree_TTJets_DiLept.root");

      // TTFileNames.push_back("tree_TTJets.root");
      // TTFileNames.push_back("tree_TTJets_DiLept_genMET-150.root");
      // TTFileNames.push_back("tree_TTJets_SingleLeptFromT_genMET-150.root");
      // TTFileNames.push_back("tree_TTJets_SingleLeptFromTbar_genMET-150.root");


      TT = new TChain("tree");
      for( int i = 0 ; i < TTFileNames.size() ; i++ ){
        TT->Add(skimType+"/"+TTFileNames[i]);
      }
      if( r == kSignal || r == kSLm || r == kSLe || r == kLowDphi ){
        ntuples.push_back(new RA2bTree(TT));
        sampleName.push_back("TT");
        fillColor.push_back(kCyan);
        lineColor.push_back(kCyan);
      }
    }
    /*
    std::vector<TString> DYFileNames;
    DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-100to200.root");
    DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-200to400.root");
    DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-400to600.root");
    DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-600toInf.root");
    DY = new TChain("tree");
    for( int i = 0 ; i < DYFileNames.size() ; i++ ){
      DY->Add(skimType+"/"+DYFileNames[i]);
      //DY->Add(skimTypeLDP+"/"+DYFileNames[i]);
    }
    //ntuples.push_back(new RA2bTree(DY));
    //sampleName.push_back("DY");
    //fillColor.push_back(kGreen);
    //lineColor.push_back(1);

    std::vector<TString> GJets0p4FileNames;
    GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-100to200.root");
    GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-200to400.root");
    GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-400to600.root");
    GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-600toInf.root");
    GJets0p4 = new TChain("tree");
    for( int i = 0 ; i < GJets0p4FileNames.size() ; i++ ){
      GJets0p4->Add(skimType+"/"+GJets0p4FileNames[i]);
    }
    if( r == kPhoton ){
      ntuples.push_back(new RA2bTree(GJets0p4));
      sampleName.push_back("GJets");
      fillColor.push_back(kGreen);
      lineColor.push_back(1);
    }

    std::vector<TString> GJetsFileNames;
    GJetsFileNames.push_back("tree_GJets_HT-100to200.root");
    GJetsFileNames.push_back("tree_GJets_HT-200to400.root");
    GJetsFileNames.push_back("tree_GJets_HT-400to600.root");
    GJetsFileNames.push_back("tree_GJets_HT-600toInf.root");
    GJets = new TChain("tree");
    for( int i = 0 ; i < GJetsFileNames.size() ; i++ ){
      GJets->Add(skimType+"/"+GJetsFileNames[i]);
    }

    // if( r == kPhoton ){
    //   ntuples.push_back(new RA2bTree(GJets));
    //   sampleName.push_back("GJets");
    //   fillColor.push_back(kGreen);
    //   lineColor.push_back(1);
    // }
    */

    if (run_QCD){
      std::vector<TString> QCDFileNames;
      QCDFileNames.push_back("tree_QCD_HT-200to300.root");
      QCDFileNames.push_back("tree_QCD_HT-300to500.root");
      QCDFileNames.push_back("tree_QCD_HT-500to700.root");
      QCDFileNames.push_back("tree_QCD_HT-700to1000.root");
      QCDFileNames.push_back("tree_QCD_HT-1000to1500.root");
      QCDFileNames.push_back("tree_QCD_HT-1500to2000.root");
      QCDFileNames.push_back("tree_QCD_HT-2000toInf.root");
      QCD = new TChain("tree");
      for( int i = 0 ; i < QCDFileNames.size() ; i++ ){
        QCD->Add(skimType+"/"+QCDFileNames[i]);
      }
      if( r == kSignal || r == kPhoton || r == kLowDphi ){
        ntuples.push_back(new RA2bTree(QCD));
        sampleName.push_back("QCD");
        fillColor.push_back(kGray);
        lineColor.push_back(1);
      }
    }

    /*
    ////////////////////////////////////////////////////////////
    // - - - - - - - - - - - DATA INPUTS - - - - - - - - - -  //
    ////////////////////////////////////////////////////////////

    std::vector<TString> HTMHTFileNames;
    HTMHTFileNames.push_back("tree_HTMHT_re2016B.root");
    HTMHTFileNames.push_back("tree_HTMHT_re2016C.root");
    HTMHTFileNames.push_back("tree_HTMHT_re2016D.root");
    HTMHTFileNames.push_back("tree_HTMHT_re2016E.root");
    HTMHTFileNames.push_back("tree_HTMHT_re2016F.root");
    HTMHTFileNames.push_back("tree_HTMHT_re2016G.root");
    HTMHTFileNames.push_back("tree_HTMHT_re2016H2.root");
    HTMHTFileNames.push_back("tree_HTMHT_re2016H3.root");
    if( r == kSignal || r == kLowDphi ){
      data = new TChain("tree");
      for( int i = 0 ; i < HTMHTFileNames.size() ; i++ ){
        data->Add(skimType+"/"+HTMHTFileNames[i]);
      }
      dataNtuple = new RA2bTree(data);
      ntuples.push_back(dataNtuple);
      sampleName.push_back("data");
      fillColor.push_back(kWhite);
      lineColor.push_back(1);
    }
    */
    /*
    std::vector<TString> SingleElectronNames;
    SingleElectronNames.push_back("tree_SingleElectron_re2016C.root");
    SingleElectronNames.push_back("tree_SingleElectron_re2016D.root");
    SingleElectronNames.push_back("tree_SingleElectron_re2016E.root");
    SingleElectronNames.push_back("tree_SingleElectron_re2016F.root");
    SingleElectronNames.push_back("tree_SingleElectron_re2016G.root");
    SingleElectronNames.push_back("tree_SingleElectron_re2016H2.root");
    SingleElectronNames.push_back("tree_SingleElectron_re2016H3.root");
    if( r == kSLe ){
      data = new TChain("tree");
      for( int i = 0 ; i < SingleElectronNames.size() ; i++ ){
        data->Add(skimType+"/"+SingleElectronNames[i]);
      }
      dataNtuple = new RA2bTree(data);
    }

    std::vector<TString> SingleMuonNames;
    SingleMuonNames.push_back("tree_SingleMuon_re2016B.root");
    SingleMuonNames.push_back("tree_SingleMuon_re2016C.root");
    SingleMuonNames.push_back("tree_SingleMuon_re2016D.root");
    SingleMuonNames.push_back("tree_SingleMuon_re2016E.root");
    SingleMuonNames.push_back("tree_SingleMuon_re2016F.root");
    SingleMuonNames.push_back("tree_SingleMuon_re2016G.root");
    SingleMuonNames.push_back("tree_SingleMuon_re2016H2.root");
    SingleMuonNames.push_back("tree_SingleMuon_re2016H3.root");
    if( r == kSLm ){
      data = new TChain("tree");
      for( int i = 0 ; i < SingleMuonNames.size() ; i++ ){
        data->Add(skimType+"/"+SingleMuonNames[i]);
      }
      dataNtuple = new RA2bTree(data);
    }

    std::vector<TString> SinglePhotonFileNames;
    SinglePhotonFileNames.push_back("tree_SinglePhoton_re2016B.root");
    SinglePhotonFileNames.push_back("tree_SinglePhoton_re2016C.root");
    SinglePhotonFileNames.push_back("tree_SinglePhoton_re2016D.root");
    SinglePhotonFileNames.push_back("tree_SinglePhoton_re2016E.root");
    SinglePhotonFileNames.push_back("tree_SinglePhoton_re2016F.root");
    SinglePhotonFileNames.push_back("tree_SinglePhoton_re2016G.root");
    SinglePhotonFileNames.push_back("tree_SinglePhoton_re2016H2.root");
    SinglePhotonFileNames.push_back("tree_SinglePhoton_re2016H3.root");
    if( r == kPhoton ){
      data = new TChain("tree");
      for( int i = 0 ; i < SinglePhotonFileNames.size() ; i++ ){
        data->Add(skimType+"/"+SinglePhotonFileNames[i]);
      }
      dataNtuple = new RA2bTree(data);
    }
    */

    if (run_T5HH){
      T5HH750 = new TChain("tree");
      T5HH1000 = new TChain("tree");
      T5HH1100 = new TChain("tree");
      T5HH1200 = new TChain("tree");
      T5HH1300 = new TChain("tree");
      T5HH1400 = new TChain("tree");
      T5HH1500 = new TChain("tree");
      T5HH1600 = new TChain("tree");
      T5HH1700 = new TChain("tree");
      T5HH1800 = new TChain("tree");
      T5HH1900 = new TChain("tree");
      T5HH2000 = new TChain("tree");
      T5HH2100 = new TChain("tree");
      T5HH2200 = new TChain("tree");

      T5HH750->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_750.root");
      T5HH1000->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_1000.root");
      T5HH1100->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_1100.root");
      T5HH1200->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_1200.root");
      T5HH1300->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_1300.root");
      T5HH1400->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_1400.root");
      T5HH1500->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_1500.root");
      T5HH1600->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_1600.root");
      T5HH1700->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_1700.root");
      T5HH1800->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_1800.root");
      T5HH1900->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_1900.root");
      T5HH2000->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_2000.root");
      T5HH2100->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_2100.root");
      T5HH2200->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/StrongProduction/tree_T5qqqqZH_2200.root");

      ntuples.push_back(new RA2bTree(T5HH750));
      ntuples.push_back(new RA2bTree(T5HH1000));
      ntuples.push_back(new RA2bTree(T5HH1100));
      ntuples.push_back(new RA2bTree(T5HH1200));
      ntuples.push_back(new RA2bTree(T5HH1300));
      ntuples.push_back(new RA2bTree(T5HH1400));
      ntuples.push_back(new RA2bTree(T5HH1500));
      ntuples.push_back(new RA2bTree(T5HH1600));
      ntuples.push_back(new RA2bTree(T5HH1700));
      ntuples.push_back(new RA2bTree(T5HH1800));
      ntuples.push_back(new RA2bTree(T5HH1900));
      ntuples.push_back(new RA2bTree(T5HH2000));
      ntuples.push_back(new RA2bTree(T5HH2100));
      ntuples.push_back(new RA2bTree(T5HH2200));

      sampleName.push_back("T5HH750");
      sampleName.push_back("T5HH1000");
      sampleName.push_back("T5HH1100");
      sampleName.push_back("T5HH1200");
      sampleName.push_back("T5HH1300");
      sampleName.push_back("T5HH1400");
      sampleName.push_back("T5HH1500");
      sampleName.push_back("T5HH1600");
      sampleName.push_back("T5HH1700");
      sampleName.push_back("T5HH1800");
      sampleName.push_back("T5HH1900");
      sampleName.push_back("T5HH2000");
      sampleName.push_back("T5HH2100");
      sampleName.push_back("T5HH2200");

      for (unsigned int i=0; i<sampleName.size(); ++i) {
        sigLineColor.push_back(kRed);
        fillColor.push_back(kRed);
      }
    } //if run_T5HH

    if (run_AllTChiHH || run_SomeTChiHH){
      //TChiHH samples
      // TChain *TChiHH400, *TChiHH425, *TChiHH450,*TChiHH475,*TChiHH500, *TChiHH525, *TChiHH550,*TChiHH575, *TChiHH600, *TChiHH625, *TChiHH650, *TChiHH675, *TChiHH700,*TChiHH725, *TChiHH750,*TChiHH775, *TChiHH800, *TChiHH825, *TChiHH850, *TChiHH875, *TChiHH900, *TChiHH925, *TChiHH950, *TChiHH975, *TChiHH1000;
      TChiHH127 = new TChain("tree");TChiHH150 = new TChain("tree");
      TChiHH175 = new TChain("tree");TChiHH200 = new TChain("tree");
      TChiHH225 = new TChain("tree");TChiHH250 = new TChain("tree");
      TChiHH275 = new TChain("tree");TChiHH300 = new TChain("tree");
      TChiHH325 = new TChain("tree");TChiHH350 = new TChain("tree");
      TChiHH375 = new TChain("tree");TChiHH400 = new TChain("tree");
      TChiHH425 = new TChain("tree");TChiHH450 = new TChain("tree");
      TChiHH475 = new TChain("tree");TChiHH500 = new TChain("tree");
      TChiHH525 = new TChain("tree");TChiHH550 = new TChain("tree");
      TChiHH575 = new TChain("tree");TChiHH600 = new TChain("tree");
      TChiHH625 = new TChain("tree");TChiHH650 = new TChain("tree");
      TChiHH675 = new TChain("tree");TChiHH700 = new TChain("tree");
      TChiHH725 = new TChain("tree");TChiHH750 = new TChain("tree");
      TChiHH775 = new TChain("tree");TChiHH800 = new TChain("tree");
      TChiHH825 = new TChain("tree");TChiHH850 = new TChain("tree");
      TChiHH875 = new TChain("tree");TChiHH900 = new TChain("tree");
      TChiHH925 = new TChain("tree");TChiHH950 = new TChain("tree");
      TChiHH975 = new TChain("tree");TChiHH1000 = new TChain("tree");


      TChiHH127->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_127_1_fast.root");
      TChiHH150->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_150_1_fast.root");
      TChiHH175->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_175_1_fast.root");
      TChiHH200->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_200_1_fast.root");
      TChiHH225->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_225_1_fast.root");
      TChiHH250->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_250_1_fast.root");
      TChiHH275->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_275_1_fast.root");
      TChiHH300->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_300_1_fast.root");
      TChiHH325->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_325_1_fast.root");
      TChiHH350->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_350_1_fast.root");
      TChiHH375->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_375_1_fast.root");

      TChiHH400->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_400_1_fast.root");
      TChiHH425->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_425_1_fast.root");
      TChiHH450->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_450_1_fast.root");
      TChiHH475->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_475_1_fast.root");
      TChiHH500->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_500_1_fast.root");
      TChiHH525->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_525_1_fast.root");
      TChiHH550->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_550_1_fast.root");
      TChiHH575->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_575_1_fast.root");
      TChiHH600->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_600_1_fast.root");
      TChiHH625->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_625_1_fast.root");
      TChiHH650->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_650_1_fast.root");
      TChiHH675->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_675_1_fast.root");
      TChiHH700->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_700_1_fast.root");
      TChiHH725->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_725_1_fast.root");
      TChiHH750->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_750_1_fast.root");
      TChiHH775->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_775_1_fast.root");
      TChiHH800->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_800_1_fast.root");
      TChiHH825->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_825_1_fast.root");
      TChiHH850->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_850_1_fast.root");
      TChiHH875->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_875_1_fast.root");
      TChiHH900->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_900_1_fast.root");
      TChiHH925->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_925_1_fast.root");
      TChiHH950->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_950_1_fast.root");
      TChiHH975->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_975_1_fast.root");
      TChiHH1000->Add("/nfs/data37/cms/emacdonald/WorkingArea/CombinedHiggs/samples/TChiHH/tree_TChiHH_HToBB_HToBB_1000_1_fast.root");

      if (run_AllTChiHH){
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
      }

      ntuples.push_back(new RA2bTree(TChiHH400));
      ntuples.push_back(new RA2bTree(TChiHH700));
      ntuples.push_back(new RA2bTree(TChiHH1000));

      if (run_AllTChiHH){
        sampleName.push_back("TChiHH127");
        sampleName.push_back("TChiHH150");
        sampleName.push_back("TChiHH175");
        sampleName.push_back("TChiHH200");
        sampleName.push_back("TChiHH225");
        sampleName.push_back("TChiHH250");
        sampleName.push_back("TChiHH275");
        sampleName.push_back("TChiHH300");
        sampleName.push_back("TChiHH325");
        sampleName.push_back("TChiHH350");
        sampleName.push_back("TChiHH375");
        sampleName.push_back("TChiHH425");
        sampleName.push_back("TChiHH450");
        sampleName.push_back("TChiHH475");
        sampleName.push_back("TChiHH500");
        sampleName.push_back("TChiHH525");
        sampleName.push_back("TChiHH550");
        sampleName.push_back("TChiHH575");
        sampleName.push_back("TChiHH600");
        sampleName.push_back("TChiHH625");
        sampleName.push_back("TChiHH650");
        sampleName.push_back("TChiHH675");
        sampleName.push_back("TChiHH725");
        sampleName.push_back("TChiHH750");
        sampleName.push_back("TChiHH775");
        sampleName.push_back("TChiHH800");
        sampleName.push_back("TChiHH825");
        sampleName.push_back("TChiHH850");
        sampleName.push_back("TChiHH875");
        sampleName.push_back("TChiHH900");
        sampleName.push_back("TChiHH925");
        sampleName.push_back("TChiHH950");
        sampleName.push_back("TChiHH975");
      }

      sampleName.push_back("TChiHH400");
      sampleName.push_back("TChiHH700");
      sampleName.push_back("TChiHH1000");

      for (unsigned int i=0; i<sampleName.size(); ++i) {
        sigLineColor.push_back(kRed);
        fillColor.push_back(kRed);
      }
    }

  };


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
