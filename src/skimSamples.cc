// ROOT/custom libraries
#include "TChain.h"
#include "RA2bTree.cc"
#include "TString.h"

#include <iostream>
#include <vector>

static const TString BASE_DIR="/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV17/";

static const TString V18Signal_DIR="/eos/uscms/store/user/rgp230/SUSY/TChiHHV17/Skims/tree_signal/";
static const TString V18T5HH_DIR = "TwoHiggsEvents/";

class skimSamples {
  public :
  TChain *WJets,*ZJets,*QCD,*SnglT,*TT,*TT_Di,*TT_SingleLept,*GJets,*GJets0p4,*Other,*DY,*TTinc;
  TChain *T5HH750, *T5HH1000, *T5HH1100,*T5HH1200,*T5HH1300,*T5HH1400,*T5HH1500,*T5HH1600,*T5HH1700,*T5HH1800,*T5HH1900,*T5HH2000,*T5HH2100,*T5HH2200,*T5HH2300,*T5HH2400,*T5HH2500;
  TChain *TChiHH127, *TChiHH150, *TChiHH175,*TChiHH200,*TChiHH225, *TChiHH250, *TChiHH275, *TChiHH300, *TChiHH325,*TChiHH350, *TChiHH375, *TChiHH400, *TChiHH425, *TChiHH450,*TChiHH475,*TChiHH500, *TChiHH525, *TChiHH550,*TChiHH575, *TChiHH600, *TChiHH625, *TChiHH650, *TChiHH675, *TChiHH700,*TChiHH725, *TChiHH750,*TChiHH775, *TChiHH800, *TChiHH825, *TChiHH850, *TChiHH875, *TChiHH900, *TChiHH925, *TChiHH950, *TChiHH975, *TChiHH1000;
  TChain *TChiHH1025, *TChiHH1050, *TChiHH1075, *TChiHH1100, *TChiHH1125,*TChiHH1150, *TChiHH1175, *TChiHH1200, *TChiHH1225, *TChiHH1250, *TChiHH1275, *TChiHH1300, *TChiHH1325, *TChiHH1350, *TChiHH1375, *TChiHH1400, *TChiHH1425, *TChiHH1450, *TChiHH1475, *TChiHH1500;
  TChain*TChiHH;
  TChain *data;
  TChain *data2017;
  TChain *data2018;
  std::vector<RA2bTree*> ntuples,signalNtuples;
  RA2bTree* dataNtuple;  RA2bTree* dataNtuple2017; RA2bTree* dataNtuple2018;
  std::vector<TString> sampleName, signalSampleName;
  std::vector<TString> dataSampleName;
  std::vector<int> fillColor, lineColor, sigLineColor;
  int NLSP_, LSP_;

  enum region {kSignal,kSignalOnly,kSLm,kSLe,kPhoton,kLowDphi, kNumRegions};
  TString regionNames[kNumRegions]={"signal", "signalOnly","SLm", "SLe", "photon", "kLowDphi"};
  TString skimType;

  skimSamples(region r=kSignal, TString Year="MC2016", int NLSP=-1, int LSP=-1) {
    skimType="";
    if (r == kSignal || r == kSignalOnly) {
      skimType=BASE_DIR+"tree_signal";
    }
    if (r == kSignalOnly) {
      NLSP_ = NLSP;
      LSP_ = LSP;
    }
    if (r == kSLm) {
      skimType=BASE_DIR+"tree_SLm";
    }
    if (r == kSLe) {
      skimType=BASE_DIR+"tree_SLe";
    }
    if (r == kPhoton) {
      skimType=BASE_DIR+"tree_GJet_CleanVars";
    }
    if (r==kLowDphi) {
      skimType=BASE_DIR+"tree_LDP";
    }

    //bools that determine which processes are run, as a subset of signal
    bool run_singleT = false;
    bool run_TT = true;
    bool run_TT_Di = false;
    bool run_TT_SL = false;
    bool run_QCD = true;
    bool run_GJets = true;
    bool run_WJets = true;
    bool run_ZJets = true;
    bool run_AllTChiHH = false;
    bool run_SomeTChiHH = false;
    bool run_TChiHH2D = false;
    bool run_AllT5HH = false;
    bool run_SomeT5HH = false;
    bool runData = false;

    if (r == kSignalOnly) {
      run_TT = false;
      run_QCD = false;
      run_WJets = false;
      run_ZJets = false;
      run_AllTChiHH = false;
      run_SomeTChiHH = false;
      run_TChiHH2D = true;
      run_AllT5HH = false;
      std::cout<<"Check Signal Only "<<std::endl;
    }

    ///////////////////////////////////////////////////////////////////////
    // - - - - - - - - - - BACKGROUND INPUTS - - - - - - - - - - - - - - //
    ///////////////////////////////////////////////////////////////////////
    if (run_ZJets) {
      std::vector<TString> ZJetsFileNames;
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-100to200_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-200to400_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-400to600_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-600to800_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-800to1200_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-1200to2500_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-2500toInf_"+Year+".root");
      ZJets = new TChain("tree");
      for (int i = 0 ; i < ZJetsFileNames.size() ; i++) {
        ZJets->Add(skimType+"/"+ZJetsFileNames[i]);
      }
      if (r == kSignal || r == kLowDphi) {
        ntuples.push_back(new RA2bTree(ZJets));
        sampleName.push_back("ZJets");
        fillColor.push_back(kGreen+1);
        lineColor.push_back(1);
      }
    }

    if (run_WJets) {
      std::vector<TString> WJetsFileNames;
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-100to200_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-200to400_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-400to600_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-600to800_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-800to1200_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-1200to2500_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-2500toInf_"+Year+".root");
      WJets = new TChain("tree");
      for (int i = 0 ; i < WJetsFileNames.size() ; i++) {
        WJets->Add(skimType+"/"+WJetsFileNames[i]);
      }
      if (r == kSignal || r == kSLm || r == kSLe || r == kLowDphi) {
        ntuples.push_back(new RA2bTree(WJets));
        sampleName.push_back("WJets");
        fillColor.push_back(kBlue);
        lineColor.push_back(1);
      }
    }

    if (run_TT) {
      std::vector<TString> TTFileNames;
      // TTFileNames.push_back("tree_TTJets_HT-600to800_"+Year+".root");
      // TTFileNames.push_back("tree_TTJets_HT-800to1200_"+Year+".root");
      // TTFileNames.push_back("tree_TTJets_HT-1200to2500_"+Year+".root");
      // TTFileNames.push_back("tree_TTJets_HT-2500toInf_"+Year+".root");

      TTFileNames.push_back("tree_TTJets_DiLept_"+Year+".root");
      TTFileNames.push_back("tree_TTJets_SingleLeptFromT_"+Year+".root");
      TTFileNames.push_back("tree_TTJets_SingleLeptFromTbar_"+Year+".root");
      if (Year.Contains("MC2018")) {
        TTFileNames.push_back("tree_TTJets_DiLept_genMET-80_"+Year+".root");
        TTFileNames.push_back("tree_TTJets_SingleLeptFromT_genMET-80_"+Year+".root");
        TTFileNames.push_back("tree_TTJets_SingleLeptFromTbar_genMET-80_"+Year+".root");
      }
      else {
        TTFileNames.push_back("tree_TTJets_DiLept_genMET-150_"+Year+".root");
        TTFileNames.push_back("tree_TTJets_SingleLeptFromT_genMET-150_"+Year+".root");
        TTFileNames.push_back("tree_TTJets_SingleLeptFromTbar_genMET-150_"+Year+".root");
      }

      TT = new TChain("tree");
      for (int i = 0 ; i < TTFileNames.size() ; i++) {
        TT->Add(skimType+"/"+TTFileNames[i]);
      }
      if (r == kSignal || r == kSLm || r == kSLe || r == kLowDphi) {
        ntuples.push_back(new RA2bTree(TT));
        sampleName.push_back("TT");
        fillColor.push_back(kCyan);
        lineColor.push_back(kCyan);
      }
    }

    if (run_TT_Di) {
      std::vector<TString> TTFileNames;
      TTFileNames.push_back("tree_TTJets_DiLept_"+Year+".root");
      if (Year.Contains("MC2018")) {
        TTFileNames.push_back("tree_TTJets_DiLept_genMET-80_"+Year+".root");
      }
      else {
        TTFileNames.push_back("tree_TTJets_DiLept_genMET-150_"+Year+".root");
      }

      TT_Di = new TChain("tree");
      for (int i = 0 ; i < TTFileNames.size() ; i++) {
        TT_Di->Add(skimType+"/"+TTFileNames[i]);
      }
      if (r == kSignal || r == kSLm || r == kSLe || r == kLowDphi) {
        ntuples.push_back(new RA2bTree(TT_Di));
        sampleName.push_back("TT_Di");
        fillColor.push_back(kCyan);
        lineColor.push_back(kCyan);
      }
    }

    if (run_TT_SL) {
      //TT SingleLept
      std::vector<TString> TT_SingleLept_FileNames;
      TT_SingleLept_FileNames.push_back("tree_TTJets_SingleLeptFromT_"+Year+".root");
      TT_SingleLept_FileNames.push_back("tree_TTJets_SingleLeptFromTbar_"+Year+".root");
      if (Year.Contains("MC2018")) {
        TT_SingleLept_FileNames.push_back("tree_TTJets_SingleLeptFromT_genMET-80_"+Year+".root");
        TT_SingleLept_FileNames.push_back("tree_TTJets_SingleLeptFromTbar_genMET-80_"+Year+".root");
      }
      else {
        TT_SingleLept_FileNames.push_back("tree_TTJets_SingleLeptFromT_genMET-150_"+Year+".root");
        TT_SingleLept_FileNames.push_back("tree_TTJets_SingleLeptFromTbar_genMET-150_"+Year+".root");
      }
      TT_SingleLept = new TChain("tree");
      for (int i = 0 ; i < TT_SingleLept_FileNames.size() ; i++) {
        TT_SingleLept->Add(skimType+"/"+TT_SingleLept_FileNames[i]);
      }
      if (r == kSignal || r == kSLm || r == kSLe || r == kLowDphi) {
        ntuples.push_back(new RA2bTree(TT_SingleLept));
        sampleName.push_back("TT_SL");
        fillColor.push_back(kCyan);
        lineColor.push_back(kCyan);
      }
    }

    if (run_QCD) {
      std::vector<TString> QCDFileNames;
      QCDFileNames.push_back("tree_QCD_HT-200to300_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-300to500_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-500to700_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-700to1000_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-1000to1500_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-1500to2000_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-2000toInf_"+Year+".root");
      QCD = new TChain("tree");
      for (int i = 0 ; i < QCDFileNames.size() ; i++) {
        QCD->Add(skimType+"/"+QCDFileNames[i]);
      }
      if (r == kSignal || r == kPhoton || r == kLowDphi) {
        ntuples.push_back(new RA2bTree(QCD));
        sampleName.push_back("QCD");
        fillColor.push_back(kGray);
        lineColor.push_back(1);
      }
    }

    if (run_singleT) {
      std::vector<TString> SnglTFileNames;
      SnglTFileNames.push_back("tree_ST_s-channel_"+Year+".root");
      SnglTFileNames.push_back("tree_ST_t-channel_antitop_"+Year+".root");
      SnglTFileNames.push_back("tree_ST_t-channel_top_"+Year+".root");
      SnglTFileNames.push_back("tree_ST_tW_antitop_"+Year+".root");
      SnglTFileNames.push_back("tree_ST_tW_top_"+Year+".root");
      SnglT = new TChain("tree");
      for (int i = 0 ; i < SnglTFileNames.size() ; i++) {
        SnglT->Add(skimType+"/"+SnglTFileNames[i]);
      }
      if (r == kSignal || r == kSLm || r == kSLe) {
        ntuples.push_back(new RA2bTree(SnglT));
        sampleName.push_back("SnglT");
        fillColor.push_back(kOrange);
        lineColor.push_back(1);
      }
    }

    if (run_GJets) {
      std::vector<TString> GJets0p4FileNames;
      GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-100to200_"+Year+".root");
      GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-200to400_"+Year+".root");
      GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-400to600_"+Year+".root");
      GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-600toInf_"+Year+".root");
      GJets0p4 = new TChain("tree");
      for (int i = 0 ; i < GJets0p4FileNames.size() ; i++) {
        GJets0p4->Add(skimType+"/"+GJets0p4FileNames[i]);
      }
      if (r == kPhoton) {
        ntuples.push_back(new RA2bTree(GJets0p4));
        sampleName.push_back("GJets");
        fillColor.push_back(kGreen);
        lineColor.push_back(1);
      }
    }


    ////////////////////////////////////////////////////////////
    // - - - - - - - - - - - DATA INPUTS - - - - - - - - - -  //
    ////////////////////////////////////////////////////////////
    /*
    std::vector<TString> METFileNames;
    METFileNames.push_back("tree_MET_2016B.root");
    METFileNames.push_back("tree_MET_2016C.root");
    METFileNames.push_back("tree_MET_2016D.root");
    METFileNames.push_back("tree_MET_2016E.root");
    METFileNames.push_back("tree_MET_2016F.root");
    METFileNames.push_back("tree_MET_2016G.root");
    METFileNames.push_back("tree_MET_2016H.root");
    // METFileNames.push_back("tree_HTMHT_re2016H3.root");
    if (r == kSignal || r == kLowDphi) {
      data = new TChain("tree");
      for (int i = 0 ; i < METFileNames.size() ; i++) {
        data->Add(skimType+"/"+METFileNames[i]);
      }
      dataNtuple = new RA2bTree(data);
      //ntuples.push_back(dataNtuple);
      //sampleName.push_back("data");
      fillColor.push_back(kWhite);
      lineColor.push_back(1);
    }

    METFileNames.resize(0);
    METFileNames.push_back("tree_MET_2017B.root");
    METFileNames.push_back("tree_MET_2017C.root");
    METFileNames.push_back("tree_MET_2017D.root");
    METFileNames.push_back("tree_MET_2017E.root");
    METFileNames.push_back("tree_MET_2017F.root");
    // METFileNames.push_back("tree_HTMHT_re2016H3.root");
    if (r == kSignal || r == kLowDphi) {
      data2017 = new TChain("tree");
      for (int i = 0 ; i < METFileNames.size() ; i++) {
        data2017->Add(skimType+"/"+METFileNames[i]);
      }
      dataNtuple = new RA2bTree(data2017);
      //ntuples.push_back(dataNtuple);
      //sampleName.push_back("data2017");
      fillColor.push_back(kWhite);
      lineColor.push_back(1);
    }

    METFileNames.resize(0);
    METFileNames.push_back("tree_MET_2018A.root");
    METFileNames.push_back("tree_MET_2018B.root");
    METFileNames.push_back("tree_MET_2018C.root");
    METFileNames.push_back("tree_MET_2018D.root");
    // METFileNames.push_back("tree_HTMHT_re2016H3.root");
    if (r == kSignal || r == kLowDphi) {
      data2018 = new TChain("tree");
      for (int i = 0 ; i < METFileNames.size() ; i++) {
        data2018->Add(skimType+"/"+METFileNames[i]);
      }
      dataNtuple = new RA2bTree(data2018);
      //ntuples.push_back(dataNtuple);
      //sampleName.push_back("data2018");
      fillColor.push_back(kWhite);
      lineColor.push_back(1);
    }
    */
    std::vector<TString> SingleElectronNames;
    std::vector<TString> SingleMuonNames;
    std::vector<TString> SinglePhotonFileNames;

    if (Year.Contains("2016") && runData){
      SingleElectronNames.push_back("tree_SingleElectron_2016B.root");
      SingleElectronNames.push_back("tree_SingleElectron_2016C.root");
      SingleElectronNames.push_back("tree_SingleElectron_2016D.root");
      SingleElectronNames.push_back("tree_SingleElectron_2016E.root");
      SingleElectronNames.push_back("tree_SingleElectron_2016F.root");
      SingleElectronNames.push_back("tree_SingleElectron_2016G.root");
      SingleElectronNames.push_back("tree_SingleElectron_2016H.root");
      if (r == kSLe) {
        data = new TChain("tree");
        for (int i = 0 ; i < SingleElectronNames.size() ; i++) {
          data->Add(skimType+"/"+SingleElectronNames[i]);
        }
        dataNtuple = new RA2bTree(data);
        ntuples.push_back(dataNtuple);
        sampleName.push_back("data");
        fillColor.push_back(kBlack);
        lineColor.push_back(1);
      }

      SingleMuonNames.push_back("tree_SingleMuon_2016B.root");
      SingleMuonNames.push_back("tree_SingleMuon_2016C.root");
      SingleMuonNames.push_back("tree_SingleMuon_2016D.root");
      SingleMuonNames.push_back("tree_SingleMuon_2016E.root");
      SingleMuonNames.push_back("tree_SingleMuon_2016F.root");
      SingleMuonNames.push_back("tree_SingleMuon_2016G.root");
      SingleMuonNames.push_back("tree_SingleMuon_2016H.root");
      if (r == kSLm) {
        data = new TChain("tree");
        for (int i = 0 ; i < SingleMuonNames.size() ; i++) {
          data->Add(skimType+"/"+SingleMuonNames[i]);
        }
        dataNtuple = new RA2bTree(data);
        ntuples.push_back(dataNtuple);
        sampleName.push_back("data");
        fillColor.push_back(kBlack);
        lineColor.push_back(1);
      }

      SinglePhotonFileNames.push_back("tree_SinglePhoton_2016B.root");
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2016C.root");
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2016D.root");
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2016E.root");
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2016F.root");
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2016G.root");
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2016H.root");
      if (r == kPhoton) {
        data = new TChain("tree");
        for (int i = 0 ; i < SinglePhotonFileNames.size() ; i++) {
          data->Add(skimType+"/"+SinglePhotonFileNames[i]);
        }
        dataNtuple = new RA2bTree(data);
        ntuples.push_back(dataNtuple);
        sampleName.push_back("data");
        fillColor.push_back(kBlack);
        lineColor.push_back(1);
      }

    }

    if (Year.Contains("2017") && runData){
      // SingleElectronNames.resize(0);
      SingleElectronNames.push_back("tree_SingleElectron_2017B.root");
      SingleElectronNames.push_back("tree_SingleElectron_2017C.root");
      SingleElectronNames.push_back("tree_SingleElectron_2017D.root");
      SingleElectronNames.push_back("tree_SingleElectron_2017E.root");
      SingleElectronNames.push_back("tree_SingleElectron_2017F.root");
      if (r == kSLe) {
        data = new TChain("tree");
        for (int i = 0 ; i < SingleElectronNames.size() ; i++) {
          data->Add(skimType+"/"+SingleElectronNames[i]);
        }
        dataNtuple = new RA2bTree(data);
        ntuples.push_back(dataNtuple);
        sampleName.push_back("data");
        fillColor.push_back(kBlack);
        lineColor.push_back(1);
      }

      // SingleMuonNames.resize(0);
      SingleMuonNames.push_back("tree_SingleMuon_2017B.root");
      SingleMuonNames.push_back("tree_SingleMuon_2017C.root");
      SingleMuonNames.push_back("tree_SingleMuon_2017D.root");
      SingleMuonNames.push_back("tree_SingleMuon_2017E.root");
      SingleMuonNames.push_back("tree_SingleMuon_2017F.root");
      if (r == kSLm) {
        data = new TChain("tree");
        for (int i = 0 ; i < SingleMuonNames.size() ; i++) {
          data->Add(skimType+"/"+SingleMuonNames[i]);
        }
        dataNtuple = new RA2bTree(data);
        ntuples.push_back(dataNtuple);
        sampleName.push_back("data");
        fillColor.push_back(kBlack);
        lineColor.push_back(1);
      }

      // SinglePhotonFileNames.resize(0);
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2017B.root");
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2017C.root");
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2017D.root");
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2017E.root");
      SinglePhotonFileNames.push_back("tree_SinglePhoton_2017F.root");
      if (r == kPhoton) {
        data = new TChain("tree");
        for (int i = 0 ; i < SinglePhotonFileNames.size() ; i++) {
          data->Add(skimType+"/"+SinglePhotonFileNames[i]);
        }
        dataNtuple = new RA2bTree(data);
        ntuples.push_back(dataNtuple);
        sampleName.push_back("data");
        fillColor.push_back(kBlack);
        lineColor.push_back(1);
      }
    }


    if (Year.Contains("2018") && runData){
      // SingleElectronNames.resize(0);
      SingleElectronNames.push_back("tree_EGamma_2018A.root");
      SingleElectronNames.push_back("tree_EGamma_2018B.root");
      SingleElectronNames.push_back("tree_EGamma_2018C.root");
      SingleElectronNames.push_back("tree_EGamma_2018D.root");
      if ( r == kSLe) {
        data = new TChain("tree");
        for (int i = 0 ; i < SingleElectronNames.size() ; i++) {
          data->Add(skimType+"/"+SingleElectronNames[i]);
        }
        dataNtuple = new RA2bTree(data);
        ntuples.push_back(dataNtuple);
        sampleName.push_back("data");
        fillColor.push_back(kBlack);
        lineColor.push_back(1);
      }

      // SingleMuonNames.resize(0);
      SingleMuonNames.push_back("tree_SingleMuon_2018A.root");
      SingleMuonNames.push_back("tree_SingleMuon_2018B.root");
      SingleMuonNames.push_back("tree_SingleMuon_2018C.root");
      SingleMuonNames.push_back("tree_SingleMuon_2018D.root");
      if (r == kSLm) {
        data = new TChain("tree");
        for (int i = 0 ; i < SingleMuonNames.size() ; i++) {
          data->Add(skimType+"/"+SingleMuonNames[i]);
        }
        dataNtuple = new RA2bTree(data);
        ntuples.push_back(dataNtuple);
        sampleName.push_back("data");
        fillColor.push_back(kBlack);
        lineColor.push_back(1);
      }

      // SinglePhotonFileNames.resize(0);
      SinglePhotonFileNames.push_back("tree_EGamma_2018A.root");
      SinglePhotonFileNames.push_back("tree_EGamma_2018B.root");
      SinglePhotonFileNames.push_back("tree_EGamma_2018C.root");
      SinglePhotonFileNames.push_back("tree_EGamma_2018D.root");
      if (r == kPhoton) {
        data = new TChain("tree");
        for (int i = 0 ; i < SinglePhotonFileNames.size() ; i++) {
          data->Add(skimType+"/"+SinglePhotonFileNames[i]);
        }
        dataNtuple = new RA2bTree(data);
        ntuples.push_back(dataNtuple);
        sampleName.push_back("data");
        fillColor.push_back(kBlack);
        lineColor.push_back(1);
      }
    }


    /////////////////////////////////////////////////////////////////////
    // - - - - - - - - - - - - - - Signal  - - - - - - - - - - - - - - //
    /////////////////////////////////////////////////////////////////////

    if (run_AllT5HH || run_SomeT5HH) {

      if (run_AllT5HH) {
        // T5HH750 = new TChain("tree");
        T5HH1000 = new TChain("tree");
        T5HH1100 = new TChain("tree");
        T5HH1200 = new TChain("tree");
        T5HH1400 = new TChain("tree");
        T5HH1500 = new TChain("tree");
        T5HH1600 = new TChain("tree");
        T5HH1800 = new TChain("tree");
        T5HH1900 = new TChain("tree");
        T5HH2000 = new TChain("tree");
        T5HH2200 = new TChain("tree");
        T5HH2300 = new TChain("tree");
        T5HH2400 = new TChain("tree");
        T5HH2500 = new TChain("tree");
      }
      T5HH1300 = new TChain("tree");
      T5HH1700 = new TChain("tree");
      T5HH2100 = new TChain("tree");

      if (run_AllT5HH) {
        // T5HH750->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_750_1_MC2016_fast.root");
        T5HH1000->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_1000_1_MC2016_fast.root");
        T5HH1100->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_1100_1_MC2016_fast.root");
        T5HH1200->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_1200_1_MC2016_fast.root");
        T5HH1400->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_1400_1_MC2016_fast.root");
        T5HH1500->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_1500_1_MC2016_fast.root");
        T5HH1600->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_1600_1_MC2016_fast.root");
        T5HH1800->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_1800_1_MC2016_fast.root");
        T5HH1900->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_1900_1_MC2016_fast.root");
        T5HH2000->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_2000_1_MC2016_fast.root");
        T5HH2200->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_2200_1_MC2016_fast.root");
        T5HH2300->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_2300_1_MC2016_fast.root");
        T5HH2400->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_2400_1_MC2016_fast.root");
        T5HH2500->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_2500_1_MC2016_fast.root");
      }
      T5HH1300->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_1300_1_MC2016_fast.root");
      T5HH1700->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_1700_1_MC2016_fast.root");
      T5HH2100->Add(V18Signal_DIR+V18T5HH_DIR+"tree_T5qqqqZH_2100_1_MC2016_fast.root");

      if (run_AllT5HH) {
        // ntuples.push_back(new RA2bTree(T5HH750));
        signalNtuples.push_back(new RA2bTree(T5HH1000));
        signalNtuples.push_back(new RA2bTree(T5HH1100));
        signalNtuples.push_back(new RA2bTree(T5HH1200));
        signalNtuples.push_back(new RA2bTree(T5HH1400));
        signalNtuples.push_back(new RA2bTree(T5HH1500));
        signalNtuples.push_back(new RA2bTree(T5HH1600));
        signalNtuples.push_back(new RA2bTree(T5HH1800));
        signalNtuples.push_back(new RA2bTree(T5HH1900));
        signalNtuples.push_back(new RA2bTree(T5HH2000));
        signalNtuples.push_back(new RA2bTree(T5HH2200));
        signalNtuples.push_back(new RA2bTree(T5HH2300));
        signalNtuples.push_back(new RA2bTree(T5HH2400));
        signalNtuples.push_back(new RA2bTree(T5HH2500));
      }
      signalNtuples.push_back(new RA2bTree(T5HH1300));
      signalNtuples.push_back(new RA2bTree(T5HH1700));
      signalNtuples.push_back(new RA2bTree(T5HH2100));

      if (run_AllT5HH) {
        // sampleName.push_back("T5HH750");
        signalSampleName.push_back("T5HH1000");
        signalSampleName.push_back("T5HH1100");
        signalSampleName.push_back("T5HH1200");
        signalSampleName.push_back("T5HH1400");
        signalSampleName.push_back("T5HH1500");
        signalSampleName.push_back("T5HH1600");
        signalSampleName.push_back("T5HH1800");
        signalSampleName.push_back("T5HH1900");
        signalSampleName.push_back("T5HH2000");
        signalSampleName.push_back("T5HH2200");
        signalSampleName.push_back("T5HH2300");
        signalSampleName.push_back("T5HH2400");
        signalSampleName.push_back("T5HH2500");
      }
      signalSampleName.push_back("T5HH1300");
      signalSampleName.push_back("T5HH1700");
      signalSampleName.push_back("T5HH2100");


      for (unsigned int i=0; i<signalSampleName.size(); ++i) {
        sigLineColor.push_back(kRed);
        fillColor.push_back(kRed);
      }
    } //if run_AllT5HH

    if (run_TChiHH2D) {
      TChiHH = new TChain("tree");
      TChiHH->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_"+TString::Format("%d_%d_",NLSP_,LSP_)+Year+"_fast.root");
      signalNtuples.push_back(new RA2bTree(TChiHH));
      signalSampleName.push_back(TString::Format("TChiHH%d_LSP%d", NLSP_, LSP_));
      sigLineColor.push_back(kRed);
      fillColor.push_back(kRed);
    }

    if (run_AllTChiHH || run_SomeTChiHH) {
      // TChiHH127 = new TChain("tree");
      TChiHH150 = new TChain("tree");
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
      TChiHH1025 = new TChain("tree");TChiHH1050 = new TChain("tree");
      TChiHH1075 = new TChain("tree");TChiHH1100 = new TChain("tree");
      TChiHH1125 = new TChain("tree");TChiHH1150 = new TChain("tree");
      TChiHH1175 = new TChain("tree");TChiHH1200 = new TChain("tree");
      TChiHH1225 = new TChain("tree");TChiHH1250 = new TChain("tree");
      TChiHH1275 = new TChain("tree");TChiHH1300 = new TChain("tree");
      TChiHH1325 = new TChain("tree");TChiHH1350 = new TChain("tree");
      TChiHH1375 = new TChain("tree");TChiHH1400 = new TChain("tree");
      TChiHH1425 = new TChain("tree");TChiHH1450 = new TChain("tree");
      TChiHH1475 = new TChain("tree");TChiHH1500 = new TChain("tree");


      // TChiHH127->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_127_1_MC2016_fast.root");
      TChiHH150->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_150_1_MC2016_fast.root");
      TChiHH175->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_175_1_MC2016_fast.root");
      TChiHH200->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_200_1_MC2016_fast.root");
      TChiHH225->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_225_1_MC2016_fast.root");
      TChiHH250->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_250_1_MC2016_fast.root");
      TChiHH275->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_275_1_MC2016_fast.root");
      TChiHH300->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_300_1_MC2016_fast.root");
      TChiHH325->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_325_1_MC2016_fast.root");
      TChiHH350->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_350_1_MC2016_fast.root");
      TChiHH375->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_375_1_MC2016_fast.root");
      TChiHH400->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_400_1_MC2016_fast.root");
      TChiHH425->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_425_1_MC2016_fast.root");
      TChiHH450->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_450_1_MC2016_fast.root");
      TChiHH475->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_475_1_MC2016_fast.root");
      TChiHH500->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_500_1_MC2016_fast.root");
      TChiHH525->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_525_1_MC2016_fast.root");
      TChiHH550->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_550_1_MC2016_fast.root");
      TChiHH575->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_575_1_MC2016_fast.root");
      TChiHH600->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_600_1_MC2016_fast.root");
      TChiHH625->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_625_1_MC2016_fast.root");
      TChiHH650->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_650_1_MC2016_fast.root");
      TChiHH675->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_675_1_MC2016_fast.root");
      TChiHH700->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_700_1_MC2016_fast.root");
      TChiHH725->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_725_1_MC2016_fast.root");
      TChiHH750->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_750_1_MC2016_fast.root");
      TChiHH775->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_775_1_MC2016_fast.root");
      TChiHH800->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_800_1_MC2016_fast.root");
      TChiHH825->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_825_1_MC2016_fast.root");
      TChiHH850->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_850_1_MC2016_fast.root");
      TChiHH875->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_875_1_MC2016_fast.root");
      TChiHH900->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_900_1_MC2016_fast.root");
      TChiHH925->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_925_1_MC2016_fast.root");
      TChiHH950->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_950_1_MC2016_fast.root");
      TChiHH975->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_975_1_MC2016_fast.root");
      TChiHH1000->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1000_1_MC2016_fast.root");
      TChiHH1025->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1025_1_MC2016_fast.root");
      TChiHH1050->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1050_1_MC2016_fast.root");
      TChiHH1075->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1075_1_MC2016_fast.root");
      TChiHH1100->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1100_1_MC2016_fast.root");
      TChiHH1125->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1125_1_MC2016_fast.root");
      TChiHH1150->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1150_1_MC2016_fast.root");
      TChiHH1175->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1175_1_MC2016_fast.root");
      TChiHH1200->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1200_1_MC2016_fast.root");
      TChiHH1225->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1225_1_MC2016_fast.root");
      TChiHH1250->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1250_1_MC2016_fast.root");
      TChiHH1275->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1275_1_MC2016_fast.root");
      TChiHH1300->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1300_1_MC2016_fast.root");
      TChiHH1325->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1325_1_MC2016_fast.root");
      TChiHH1350->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1350_1_MC2016_fast.root");
      TChiHH1375->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1375_1_MC2016_fast.root");
      TChiHH1400->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1400_1_MC2016_fast.root");
      TChiHH1425->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1425_1_MC2016_fast.root");
      TChiHH1450->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1450_1_MC2016_fast.root");
      TChiHH1475->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1475_1_MC2016_fast.root");
      TChiHH1500->Add(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_1500_1_MC2016_fast.root");

      if (run_SomeTChiHH){
        signalNtuples.push_back(new RA2bTree(TChiHH400));
        signalNtuples.push_back(new RA2bTree(TChiHH700));
        signalNtuples.push_back(new RA2bTree(TChiHH900));
        signalNtuples.push_back(new RA2bTree(TChiHH1200));
        signalSampleName.push_back("TChiHH400");
        signalSampleName.push_back("TChiHH700");
        signalSampleName.push_back("TChiHH900");
        signalSampleName.push_back("TChiHH1200");
      }


      if (run_AllTChiHH) {
        // ntuples.push_back(new RA2bTree(TChiHH127));
        signalNtuples.push_back(new RA2bTree(TChiHH150));
        signalNtuples.push_back(new RA2bTree(TChiHH175));
        signalNtuples.push_back(new RA2bTree(TChiHH200));
        signalNtuples.push_back(new RA2bTree(TChiHH225));
        signalNtuples.push_back(new RA2bTree(TChiHH250));
        signalNtuples.push_back(new RA2bTree(TChiHH275));
        signalNtuples.push_back(new RA2bTree(TChiHH300));
        signalNtuples.push_back(new RA2bTree(TChiHH325));
        signalNtuples.push_back(new RA2bTree(TChiHH350));
        signalNtuples.push_back(new RA2bTree(TChiHH375));
        signalNtuples.push_back(new RA2bTree(TChiHH400));
        signalNtuples.push_back(new RA2bTree(TChiHH425));
        signalNtuples.push_back(new RA2bTree(TChiHH450));
        signalNtuples.push_back(new RA2bTree(TChiHH475));
        signalNtuples.push_back(new RA2bTree(TChiHH500));
        signalNtuples.push_back(new RA2bTree(TChiHH525));
        signalNtuples.push_back(new RA2bTree(TChiHH550));
        signalNtuples.push_back(new RA2bTree(TChiHH575));
        signalNtuples.push_back(new RA2bTree(TChiHH600));
        signalNtuples.push_back(new RA2bTree(TChiHH625));
        signalNtuples.push_back(new RA2bTree(TChiHH650));
        signalNtuples.push_back(new RA2bTree(TChiHH675));
        signalNtuples.push_back(new RA2bTree(TChiHH700));
        signalNtuples.push_back(new RA2bTree(TChiHH725));
        signalNtuples.push_back(new RA2bTree(TChiHH750));
        signalNtuples.push_back(new RA2bTree(TChiHH775));
        signalNtuples.push_back(new RA2bTree(TChiHH800));
        signalNtuples.push_back(new RA2bTree(TChiHH825));
        signalNtuples.push_back(new RA2bTree(TChiHH850));
        signalNtuples.push_back(new RA2bTree(TChiHH875));
        signalNtuples.push_back(new RA2bTree(TChiHH900));
        signalNtuples.push_back(new RA2bTree(TChiHH925));
        signalNtuples.push_back(new RA2bTree(TChiHH950));
        signalNtuples.push_back(new RA2bTree(TChiHH975));
        signalNtuples.push_back(new RA2bTree(TChiHH1000));
        signalNtuples.push_back(new RA2bTree(TChiHH1025));
        signalNtuples.push_back(new RA2bTree(TChiHH1050));
        signalNtuples.push_back(new RA2bTree(TChiHH1075));
        signalNtuples.push_back(new RA2bTree(TChiHH1100));
        signalNtuples.push_back(new RA2bTree(TChiHH1125));
        signalNtuples.push_back(new RA2bTree(TChiHH1150));
        signalNtuples.push_back(new RA2bTree(TChiHH1175));
        signalNtuples.push_back(new RA2bTree(TChiHH1200));
        signalNtuples.push_back(new RA2bTree(TChiHH1225));
        signalNtuples.push_back(new RA2bTree(TChiHH1250));
        signalNtuples.push_back(new RA2bTree(TChiHH1275));
        signalNtuples.push_back(new RA2bTree(TChiHH1300));
        signalNtuples.push_back(new RA2bTree(TChiHH1325));
        signalNtuples.push_back(new RA2bTree(TChiHH1350));
        signalNtuples.push_back(new RA2bTree(TChiHH1375));
        signalNtuples.push_back(new RA2bTree(TChiHH1400));
        signalNtuples.push_back(new RA2bTree(TChiHH1425));
        signalNtuples.push_back(new RA2bTree(TChiHH1450));
        signalNtuples.push_back(new RA2bTree(TChiHH1475));
        signalNtuples.push_back(new RA2bTree(TChiHH1500));

        // signalSampleName.push_back("TChiHH127");
        signalSampleName.push_back("TChiHH150");
        signalSampleName.push_back("TChiHH175");
        signalSampleName.push_back("TChiHH200");
        signalSampleName.push_back("TChiHH225");
        signalSampleName.push_back("TChiHH250");
        signalSampleName.push_back("TChiHH275");
        signalSampleName.push_back("TChiHH300");
        signalSampleName.push_back("TChiHH325");
        signalSampleName.push_back("TChiHH350");
        signalSampleName.push_back("TChiHH375");
        signalSampleName.push_back("TChiHH400");
        signalSampleName.push_back("TChiHH425");
        signalSampleName.push_back("TChiHH450");
        signalSampleName.push_back("TChiHH475");
        signalSampleName.push_back("TChiHH500");
        signalSampleName.push_back("TChiHH525");
        signalSampleName.push_back("TChiHH550");
        signalSampleName.push_back("TChiHH575");
        signalSampleName.push_back("TChiHH600");
        signalSampleName.push_back("TChiHH625");
        signalSampleName.push_back("TChiHH650");
        signalSampleName.push_back("TChiHH675");
        signalSampleName.push_back("TChiHH700");
        signalSampleName.push_back("TChiHH725");
        signalSampleName.push_back("TChiHH750");
        signalSampleName.push_back("TChiHH775");
        signalSampleName.push_back("TChiHH800");
        signalSampleName.push_back("TChiHH825");
        signalSampleName.push_back("TChiHH850");
        signalSampleName.push_back("TChiHH875");
        signalSampleName.push_back("TChiHH900");
        signalSampleName.push_back("TChiHH925");
        signalSampleName.push_back("TChiHH950");
        signalSampleName.push_back("TChiHH975");
        signalSampleName.push_back("TChiHH1000");
        signalSampleName.push_back("TChiHH1025");
        signalSampleName.push_back("TChiHH1050");
        signalSampleName.push_back("TChiHH1075");
        signalSampleName.push_back("TChiHH1100");
        signalSampleName.push_back("TChiHH1125");
        signalSampleName.push_back("TChiHH1150");
        signalSampleName.push_back("TChiHH1175");
        signalSampleName.push_back("TChiHH1200");
        signalSampleName.push_back("TChiHH1225");
        signalSampleName.push_back("TChiHH1250");
        signalSampleName.push_back("TChiHH1275");
        signalSampleName.push_back("TChiHH1300");
        signalSampleName.push_back("TChiHH1325");
        signalSampleName.push_back("TChiHH1350");
        signalSampleName.push_back("TChiHH1375");
        signalSampleName.push_back("TChiHH1400");
        signalSampleName.push_back("TChiHH1425");
        signalSampleName.push_back("TChiHH1450");
        signalSampleName.push_back("TChiHH1475");
        signalSampleName.push_back("TChiHH1500");
      }


      for (unsigned int i=0; i<signalSampleName.size(); ++i) {
        sigLineColor.push_back(kRed);
        fillColor.push_back(kRed);
      }
    }
  }; //end Skim samples

  RA2bTree* findNtuple(TString name) {
    for (int iSam = 0 ; iSam < sampleName.size() ; iSam++) {
      if (sampleName[iSam] == name)
      return ntuples[iSam] ;
    }
    for (int iSam = 0 ; iSam < signalSampleName.size() ; iSam++) {
      if (signalSampleName[iSam] == name)
      return signalNtuples[iSam] ;
    }
    return NULL;
  };
};
