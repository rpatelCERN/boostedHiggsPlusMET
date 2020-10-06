// ROOT/custom libraries
#include "TChain.h"
#include "RA2bTree.cc"
#include "TString.h"

#include <cstdlib>

#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <limits>

#include <getopt.h>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TDirectory.h"

#include <iostream>
#include <vector>

std::vector<std::string> split(const std::string& s, char delimiter) {
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter)) { tokens.push_back(token);}
   return tokens;
}

static const TString BASE_DIR="/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV18/";
static const TString V18Signal_DIR="/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV18/scan/tree_signal/";

class skimSamples {
  public :
  TChain *WJets,*ZJets,*QCD,*SnglT,*TT,*TT_Di,*TT_SingleLept,*GJets,*GJets0p4,*Other,*DY,*TTinc;
  TChain *T5HH, *TChiHH, *TChiHH1D;
  TChain *data; TChain *data2017; TChain *data2018;

  std::vector<RA2bTree*> ntuples,signalNtuples;
  RA2bTree* dataNtuple;  RA2bTree* dataNtuple2017; RA2bTree* dataNtuple2018;
  std::vector<TString> sampleName, signalSampleName;
  std::vector<TString> dataSampleName;

  std::vector<int> fillColor, lineColor, sigLineColor;

  enum region {kSignal,kSignalOnly,kSLm,kSLe,kPhoton,kLowDphi, kNumRegions};
  TString regionNames[kNumRegions]={"signal", "signalOnly","SLm", "SLe", "photon", "kLowDphi"};
  TString skimType;
  int LSP=0;

  skimSamples(region r=kSignal, TString Year="MC2016") {
    skimType="";
    if (r == kSignal || r == kSignalOnly) {
      skimType=BASE_DIR+"tree_signal";
    }
    if (r == kSignalOnly) {  }
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
    bool run_singleT = true;
    bool run_TT = true;
    bool run_TT_Di = false;
    bool run_TT_SL = false;
    bool run_QCD = true;
    bool run_GJets = true;
    bool run_WJets = true;
    bool run_ZJets = true;
    bool run_TChiHH1D = false;
    bool run_TChiHH2D = false;
    bool run_T5HH = false;
    bool runData = false;

    if (r == kSignalOnly) {
      run_TT = false;
      run_QCD = false;
      run_WJets = false;
      run_ZJets = false;
      run_TChiHH1D = true;
      run_TChiHH2D = false;
      run_T5HH = true;
      if (LSP>0) {
        std::cout<<"Check Signal Only (Higgsino 2D)"<<std::endl;
        run_TChiHH1D = false;
        run_TChiHH2D = true;
        run_T5HH = false;
      }
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

    if (Year.Contains("2016") && runData) {
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
    } // end 2016 data

    if (Year.Contains("2017") && runData) {
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
    } // end 2017 data


    if (Year.Contains("2018") && runData) {
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
    } //end 2018 data


    /////////////////////////////////////////////////////////////////////
    // - - - - - - - - - - - - - - Signal  - - - - - - - - - - - - - - //
    /////////////////////////////////////////////////////////////////////

    if (run_T5HH) { //goes from 1000-2500 GeV, all 3 years available
      // for (int i=1000;i<1400;i+=100) { //tree_T5qqqqZH-mGluino-1000to2500_1000_1_MC2016.root
      for (int i=1000;i<2600;i+=100) { //tree_T5qqqqZH-mGluino-1000to2500_1000_1_MC2016.root
        T5HH = new TChain("tree");
        TString fileName = V18Signal_DIR+"tree_T5qqqqZH-mGluino-1000to2500_"+TString::Format("%d",i)+"_1_"+Year+".root";
        T5HH->Add(fileName);
        ntuples.push_back(new RA2bTree(T5HH));
        sampleName.push_back("T5qqqqZH_"+TString::Format("%d",i)+"_1");
        fillColor.push_back(kRed);
        lineColor.push_back(kRed);
      }
    }

    if (run_TChiHH2D) {
      ifstream file("higgsino2DFileNames.txt");
      string line;
      TString fileName;
      while(std::getline(file, line)) {
        // std::cout<<"Line: "<<line<<std::endl;
        std::vector<std::string> x = split(line, '_');
        int hino_mass = std::stoi(x[5]);
        int LSP_mass = std::stoi(x[6]);
        // std::cout<<"mNLSP: "<<hino_mass<<", mLSP: "<<LSP_mass<<std::endl;

        TChiHH = new TChain("tree");
        fileName = V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_2D_"+TString::Format("%d_%d_",hino_mass,LSP_mass)+Year+"_fast.root";
        TChiHH->Add(fileName);
        ntuples.push_back(new RA2bTree(TChiHH));
        sampleName.push_back(TString::Format("TChiHH%d_LSP%d", hino_mass, LSP_mass));
        sigLineColor.push_back(kRed);
        fillColor.push_back(kRed);
        if (hino_mass==1200 && LSP_mass==950) break;
      }
    }

    if (run_TChiHH1D) { // from 150-1500 GeV, every 25 GeV, all three years
      // for (int i=1000;i<1400;i+=100) {
      for (int i=150;i<1525;i+=25) {
        TChiHH1D = new TChain("tree");
        TString fileName = V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_"+TString::Format("%d",i)+"_1_"+Year+"_fast.root";
        TChiHH1D->Add(fileName);
        ntuples.push_back(new RA2bTree(TChiHH1D));
        sampleName.push_back("TChiHH_"+TString::Format("%d",i)+"_1");
        fillColor.push_back(kRed);
        lineColor.push_back(kRed);
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
