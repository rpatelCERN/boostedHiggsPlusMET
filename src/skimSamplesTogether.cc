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

class skimSamplesTogether {
  public :
  TChain *WJets,*ZJets,*QCD,*SnglT,*TT,*TT_Di,*TT_SingleLept,*GJets,*GJets0p4,*Other,*DY,*TTinc;
  TChain *totalBkg;
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

  skimSamplesTogether(region r=kSignal, TString Year="MC2016") {
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


    ///////////////////////////////////////////////////////////////////////
    // - - - - - - - - - - BACKGROUND INPUTS - - - - - - - - - - - - - - //
    ///////////////////////////////////////////////////////////////////////
    std::vector<TString> ZJetsFileNames;
    std::vector<TString> WJetsFileNames;
    std::vector<TString> TTFileNames;
    std::vector<TString> SnglTFileNames;
    std::vector<TString> QCDFileNames;
    std::vector<TString> GJets0p4FileNames;


    if (r == kSignal) {
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-100to200_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-200to400_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-400to600_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-600to800_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-800to1200_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-1200to2500_"+Year+".root");
      ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-2500toInf_"+Year+".root");

      WJetsFileNames.push_back("tree_WJetsToLNu_HT-100to200_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-200to400_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-400to600_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-600to800_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-800to1200_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-1200to2500_"+Year+".root");
      WJetsFileNames.push_back("tree_WJetsToLNu_HT-2500toInf_"+Year+".root");

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

      SnglTFileNames.push_back("tree_ST_s-channel_"+Year+".root");
      SnglTFileNames.push_back("tree_ST_t-channel_antitop_"+Year+".root");
      SnglTFileNames.push_back("tree_ST_t-channel_top_"+Year+".root");
      SnglTFileNames.push_back("tree_ST_tW_antitop_"+Year+".root");
      SnglTFileNames.push_back("tree_ST_tW_top_"+Year+".root");
    }
    if (r == kSignal || r == kPhoton) {
      QCDFileNames.push_back("tree_QCD_HT-200to300_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-300to500_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-500to700_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-700to1000_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-1000to1500_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-1500to2000_"+Year+".root");
      QCDFileNames.push_back("tree_QCD_HT-2000toInf_"+Year+".root");
    }

    if (r == kPhoton) {
      GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-100to200_"+Year+".root");
      GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-200to400_"+Year+".root");
      GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-400to600_"+Year+".root");
      GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-600toInf_"+Year+".root");
    }

    totalBkg = new TChain("tree");
    if (r == kSignal) {
      for (int i = 0 ; i < ZJetsFileNames.size() ; i++) {
        totalBkg->Add(skimType+"/"+ZJetsFileNames[i]);
      }
      for (int i = 0 ; i < WJetsFileNames.size() ; i++) {
        totalBkg->Add(skimType+"/"+WJetsFileNames[i]);
      }
      for (int i = 0 ; i < TTFileNames.size() ; i++) {
        totalBkg->Add(skimType+"/"+TTFileNames[i]);
      }
      for (int i = 0 ; i < SnglTFileNames.size() ; i++) {
        totalBkg->Add(skimType+"/"+SnglTFileNames[i]);
      }
    }
    if (r == kSignal || r == kPhoton) {
      for (int i = 0 ; i < QCDFileNames.size() ; i++) {
        totalBkg->Add(skimType+"/"+QCDFileNames[i]);
      }
    }
    if (r == kPhoton) {
      for (int i = 0 ; i < GJets0p4FileNames.size() ; i++) {
        totalBkg->Add(skimType+"/"+GJets0p4FileNames[i]);
      }
    }
    if (r != kSignalOnly) {
      ntuples.push_back(new RA2bTree(totalBkg));
      sampleName.push_back("TotalBkg");
      fillColor.push_back(kGreen+1);
      lineColor.push_back(1);
    }


    /////////////////////////////////////////////////////////////////////
    // - - - - - - - - - - - - - - Signal  - - - - - - - - - - - - - - //
    /////////////////////////////////////////////////////////////////////

    if (r == kSignal || r==kSignalOnly) {
      for (int i=1000;i<2600;i+=100) { //tree_T5qqqqZH-mGluino-1000to2500_1000_1_MC2016.root
        T5HH = new TChain("tree");
        TString fileName = V18Signal_DIR+"tree_T5qqqqZH-mGluino-1000to2500_"+TString::Format("%d",i)+"_1_"+Year+".root";
        T5HH->Add(fileName);
        ntuples.push_back(new RA2bTree(T5HH));
        sampleName.push_back("T5HH_"+TString::Format("%d",i));
        fillColor.push_back(kRed);
        lineColor.push_back(kRed);
      }

      for (int i=150;i<1525;i+=25) {
        TChiHH1D = new TChain("tree");
        TString fileName = V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_"+TString::Format("%d",i)+"_1_"+Year+"_fast.root";
        TChiHH1D->Add(fileName);
        ntuples.push_back(new RA2bTree(TChiHH1D));
        sampleName.push_back("TChiHH_"+TString::Format("%d",i));
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
