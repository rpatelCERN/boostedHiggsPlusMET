#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "THStack.h"
#include "TPad.h"
#include "TLorentzVector.h"

#include <vector>
#include <map>
#include <iostream>
#include <assert.h>

#include "plotterUtils.cc"
#include "skimSamples.cc"
#include "definitions.cc"
#include "RA2bTree.cc"
#include "ALPHABET.h"
// #include "TriggerCorrector.h"
#include "apply_trigeffs.cpp"

using namespace std;
using namespace alphabet;

int main(int argc, char** argv) {
  std::cout<<"In ALPHABET: start"<<std::endl;

  int region(0);
  int runVeto(0);
  bool runData = false;

  // gROOT->ProcessLine(".L tdrstyle.C");
  // gROOT->ProcessLine("setTDRStyle()");

  region = atoi(argv[1]);
  TString Year(argv[2]);
  runVeto = atoi(argv[3]);
  cout <<"In ALPHABET: region: "<<region<<", veto?: "<<runVeto<<endl;

  // TriggerCorrector trigcorror; TriggerCorrector trigcorrorHT; TriggerCorrector trigcorrorFakeMHT;
	// TriggerCorrector trigcorrorPhoBarrel; TriggerCorrector trigcorrorPhoEndcap;

	// if (Year=="MC2016") {
    // trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2016");
    // trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2016");
    // trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2016");
    // trigcorrorPhoBarrel.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonBarrelLoose_hists_Run2016_JetHT");
    // trigcorrorPhoEndcap.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonEndcapLoose_hists_Run2016_JetHT");
 	// }
	// if (Year=="MC2017") {
    // trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2017");
    // trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2017");
    // trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2017");
    // trigcorrorPhoBarrel.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonBarrelLoose_hists_Run2017_JetHT");
    // trigcorrorPhoEndcap.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonEndcapLoose_hists_Run2017_JetHT");
  // }
	// if (Year=="MC2018") {
    // trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2018");
    // trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2018");
    // trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2018");
    // trigcorrorPhoBarrel.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonBarrelLoose_hists_Run2018_JetHT");
    // trigcorrorPhoEndcap.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonEndcapLoose_hists_Run2018_JetHT");
  // }

  skimSamples* skims_;
  if (region == 0 ) skims_ = new skimSamples(skimSamples::kSignal, Year);
  else if (region == 1 ) skims_ = new skimSamples(skimSamples::kSignalOnly, Year);
  else if (region == 2 ) skims_ = new skimSamples(skimSamples::kSLm, Year);
  else if (region == 3 ) skims_ = new skimSamples(skimSamples::kSLe, Year);
  else if (region == 4 ) skims_ = new skimSamples(skimSamples::kPhoton, Year);
  else if (region == 5 ) skims_ = new skimSamples(skimSamples::kLowDphi, Year);
  else assert(1);

  std::cout<<"In ALPHABET: after skimSamples"<<std::endl;

  typedef bool(*cuts)(RA2bTree*);
  vector<cuts> baselineCuts;


  if (region == 0 || region == 1) baselineCuts.push_back(*boostedBaselineCut<RA2bTree>); //boostedBaselineCut_loose
  else if (region == 2) baselineCuts.push_back(*singleMuBaselineCut<RA2bTree>);
  else if (region == 3) baselineCuts.push_back(*singleEleBaselineCut<RA2bTree>);
  else if (region == 4) baselineCuts.push_back(*photonBaselineCut<RA2bTree>); //photonBaselineCut_loose (usual cuts)
  else if (region == 5) baselineCuts.push_back(*lowDphiBaselineCut<RA2bTree>);
  else assert(1);

  skimSamples skims = *skims_;

  TFile* outputFile;
  TString regionName;
  TString cutName="";
  TString METcut = "_MET300";
  // TString METcut = "_PUWeights";

  if (region == 0) regionName="";
  if (region == 1) {
    regionName = "_1Dsignal";
    // regionName = "_2Dsignal";
  }
  if (region == 2) regionName="_singleMu";
  if (region == 3) regionName="_singleEle";
  if (region == 4) regionName="_photon";
  if (region == 5) regionName="_lowDphi";
  if (runVeto) cutName = "_resVeto";
  if (runData) METcut = "_data";

  outputFile = new TFile("ALPHABET"+Year+"_V18"+regionName+METcut+cutName+"_2BoostedH_test.root","RECREATE");

  typedef plot<RA2bTree> plot;
  vector<vector<plot> > plots;

  double mJbins[4]={baselineMassLow,HmassWindowLow,HmassWindowHigh,baselineMassHigh};

  double METbinsPhoton[5]={100.,300.,500.,700.,1400.};
  double METbinsLow[5]={250.,300.,500.,700.,1400.};
  double METbinsExtra[6]={300.,400.,500.,600.,700.,1400.}; //these are for looking at MET shape cut options
  double METbins[4]={300.,500.,700.,1400.};

  TString METbinsPhoton_string[5]={"100","300","500","700","1400"};
  TString METbinsLow_string[5]={"250","300","500","700","1400"};
  TString METbins_string[4]={"300","500","700","1400"};

  for (int i = 0; i < numMETbins ; i++) {
    TString tag="_";
    tag+=METbins_string[i];
    vector<plot> plotsTemp;
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_doubletagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_doubletagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_tagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_tagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSROpt1"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSBOpt1"+tag,"m_{J} [GeV]",3,mJbins));
    plots.push_back(plotsTemp);
  }

  //vector<plot> tempPlots;
  plot METPhoton_Plot(*fillMET<RA2bTree>,"METPhoton","MET [GeV]",4,METbinsPhoton);
  plot MET_Plot(*fillMET<RA2bTree>,"MET","MET [GeV]",3,METbins);
  plot HT_Plot(*fillHT<RA2bTree>,"HT","HT [GeV]",100,0.,5000.);
  plot METall_Plot(*fillMET<RA2bTree>,"METall","MET [GeV]",100,0.,3000.);
  plot MET5_Plot(*fillMET<RA2bTree>,"MET5","MET [GeV]",5,METbinsExtra);

  plot J1pt_Ptplot(*fillLeadingJetPt<RA2bTree>,"J1pt_Pt","p_{T,J} [GeV]",50,300.,1300.);
  plot J1pt_Mplot(*fillLeadingJetMass<RA2bTree>,"J1pt_M","m_{J} [GeV]",80,60.,260.);
  plot J1_M_jetBins(*fillLeadingJetMass<RA2bTree>,"J1_M_jetBins","m_{J} [GeV]",3,mJbins);
  plot J2pt_Ptplot(*fillSubLeadingJetPt<RA2bTree>,"J2pt_Pt","p_{T,J} [GeV]",50,300.,1300.);
  plot J2pt_Mplot(*fillSubLeadingJetMass<RA2bTree>,"J2pt_M","m_{J} [GeV]",80,60.,260.);
  plot J2_M_jetBins(*fillSubLeadingJetMass<RA2bTree>,"J2_M_jetBins","m_{J} [GeV]",3,mJbins);
  plot J1_deepbbtag(*fillLeadingdeepBBtag<RA2bTree>,"LeadDeepBBTag","Lead jet deep bb-tag",50,0.0,1.0);
  plot J2_deepbbtag(*fillSubLeadingdeepBBtag<RA2bTree>,"SubLeadDeepBBTag","Sub-lead jet deep bb-tag",50,0.0,1.0);
  plot nVert(*fillNVtx<RA2bTree>,"NVert","NVtx",100,0.0,100.0);

  vector<plot> baselinePlots;
  baselinePlots.push_back(plot(MET_Plot));
  baselinePlots.push_back(plot(METall_Plot));
  baselinePlots.push_back(plot(METPhoton_Plot));
  baselinePlots.push_back(plot(MET5_Plot));
  baselinePlots.push_back(plot(HT_Plot));
  baselinePlots.push_back(plot(J1pt_Ptplot));
  baselinePlots.push_back(plot(J2pt_Ptplot));
  baselinePlots.push_back(plot(J1pt_Mplot));
  baselinePlots.push_back(plot(J2pt_Mplot));
  baselinePlots.push_back(plot(J1_M_jetBins));
  baselinePlots.push_back(plot(J2_M_jetBins));
  baselinePlots.push_back(plot(J1_deepbbtag));
  baselinePlots.push_back(plot(J2_deepbbtag));
  baselinePlots.push_back(plot(nVert));

  vector<plot> doubletagSRPlots;
  doubletagSRPlots.push_back(plot(MET_Plot));
  doubletagSRPlots.push_back(plot(METall_Plot));
  doubletagSRPlots.push_back(plot(METPhoton_Plot));
  doubletagSRPlots.push_back(plot(MET5_Plot));
  doubletagSRPlots.push_back(plot(HT_Plot));
  doubletagSRPlots.push_back(plot(J1pt_Ptplot));
  doubletagSRPlots.push_back(plot(J2pt_Ptplot));
  doubletagSRPlots.push_back(plot(J1pt_Mplot));
  doubletagSRPlots.push_back(plot(J2pt_Mplot));
  doubletagSRPlots.push_back(plot(J1_M_jetBins));
  doubletagSRPlots.push_back(plot(J2_M_jetBins));
  doubletagSRPlots.push_back(plot(J1_deepbbtag));
  doubletagSRPlots.push_back(plot(J2_deepbbtag));
  doubletagSRPlots.push_back(plot(nVert));

  vector<plot> doubletagSBPlots;
  doubletagSBPlots.push_back(plot(MET_Plot));
  doubletagSBPlots.push_back(plot(METall_Plot));
  doubletagSBPlots.push_back(plot(METPhoton_Plot));
  doubletagSBPlots.push_back(plot(MET5_Plot));
  doubletagSBPlots.push_back(plot(HT_Plot));
  doubletagSBPlots.push_back(plot(J1pt_Ptplot));
  doubletagSBPlots.push_back(plot(J2pt_Ptplot));
  doubletagSBPlots.push_back(plot(J1pt_Mplot));
  doubletagSBPlots.push_back(plot(J2pt_Mplot));
  doubletagSBPlots.push_back(plot(J1_M_jetBins));
  doubletagSBPlots.push_back(plot(J2_M_jetBins));
  doubletagSBPlots.push_back(plot(J1_deepbbtag));
  doubletagSBPlots.push_back(plot(J2_deepbbtag));
  doubletagSBPlots.push_back(plot(nVert));

  vector<plot> tagSRPlots;
  tagSRPlots.push_back(plot(MET_Plot));
  tagSRPlots.push_back(plot(METall_Plot));
  tagSRPlots.push_back(plot(METPhoton_Plot));
  tagSRPlots.push_back(plot(MET5_Plot));
  tagSRPlots.push_back(plot(HT_Plot));
  tagSRPlots.push_back(plot(J1pt_Ptplot));
  tagSRPlots.push_back(plot(J2pt_Ptplot));
  tagSRPlots.push_back(plot(J1pt_Mplot));
  tagSRPlots.push_back(plot(J2pt_Mplot));
  tagSRPlots.push_back(plot(J1_M_jetBins));
  tagSRPlots.push_back(plot(J2_M_jetBins));
  tagSRPlots.push_back(plot(J1_deepbbtag));
  tagSRPlots.push_back(plot(J2_deepbbtag));
  tagSRPlots.push_back(plot(nVert));

  vector<plot> tagSBPlots;
  tagSBPlots.push_back(plot(MET_Plot));
  tagSBPlots.push_back(plot(METall_Plot));
  tagSBPlots.push_back(plot(METPhoton_Plot));
  tagSBPlots.push_back(plot(MET5_Plot));
  tagSBPlots.push_back(plot(HT_Plot));
  tagSBPlots.push_back(plot(J1pt_Ptplot));
  tagSBPlots.push_back(plot(J2pt_Ptplot));
  tagSBPlots.push_back(plot(J1pt_Mplot));
  tagSBPlots.push_back(plot(J2pt_Mplot));
  tagSBPlots.push_back(plot(J1_M_jetBins));
  tagSBPlots.push_back(plot(J2_M_jetBins));
  tagSBPlots.push_back(plot(J1_deepbbtag));
  tagSBPlots.push_back(plot(J2_deepbbtag));
  tagSBPlots.push_back(plot(nVert));

  vector<plot> antitagSRPlots;
  antitagSRPlots.push_back(plot(MET_Plot));
  antitagSRPlots.push_back(plot(METall_Plot));
  antitagSRPlots.push_back(plot(METPhoton_Plot));
  antitagSRPlots.push_back(plot(MET5_Plot));
  antitagSRPlots.push_back(plot(HT_Plot));
  antitagSRPlots.push_back(plot(J1pt_Ptplot));
  antitagSRPlots.push_back(plot(J2pt_Ptplot));
  antitagSRPlots.push_back(plot(J1pt_Mplot));
  antitagSRPlots.push_back(plot(J2pt_Mplot));
  antitagSRPlots.push_back(plot(J1_M_jetBins));
  antitagSRPlots.push_back(plot(J2_M_jetBins));
  antitagSRPlots.push_back(plot(J1_deepbbtag));
  antitagSRPlots.push_back(plot(J2_deepbbtag));
  antitagSRPlots.push_back(plot(nVert));

  vector<plot> antitagSBPlots;
  antitagSBPlots.push_back(plot(MET_Plot));
  antitagSBPlots.push_back(plot(METall_Plot));
  antitagSBPlots.push_back(plot(METPhoton_Plot));
  antitagSBPlots.push_back(plot(MET5_Plot));
  antitagSBPlots.push_back(plot(HT_Plot));
  antitagSBPlots.push_back(plot(J1pt_Ptplot));
  antitagSBPlots.push_back(plot(J2pt_Ptplot));
  antitagSBPlots.push_back(plot(J1pt_Mplot));
  antitagSBPlots.push_back(plot(J2pt_Mplot));
  antitagSBPlots.push_back(plot(J1_M_jetBins));
  antitagSBPlots.push_back(plot(J2_M_jetBins));
  antitagSBPlots.push_back(plot(J1_deepbbtag));
  antitagSBPlots.push_back(plot(J2_deepbbtag));
  antitagSBPlots.push_back(plot(nVert));

  vector<plot> antitagSROpt1Plots; //this is BTagsT>0 now, so used for the MET shape
  antitagSROpt1Plots.push_back(plot(MET_Plot));
  antitagSROpt1Plots.push_back(plot(METall_Plot));
  antitagSROpt1Plots.push_back(plot(METPhoton_Plot));
  antitagSROpt1Plots.push_back(plot(MET5_Plot));
  antitagSROpt1Plots.push_back(plot(HT_Plot));
  antitagSROpt1Plots.push_back(plot(J1pt_Ptplot));
  antitagSROpt1Plots.push_back(plot(J2pt_Ptplot));
  antitagSROpt1Plots.push_back(plot(J1pt_Mplot));
  antitagSROpt1Plots.push_back(plot(J2pt_Mplot));
  antitagSROpt1Plots.push_back(plot(J1_M_jetBins));
  antitagSROpt1Plots.push_back(plot(J2_M_jetBins));
  antitagSROpt1Plots.push_back(plot(J1_deepbbtag));
  antitagSROpt1Plots.push_back(plot(J2_deepbbtag));
  antitagSROpt1Plots.push_back(plot(nVert));

  vector<plot> antitagSBOpt1Plots;
  antitagSBOpt1Plots.push_back(plot(MET_Plot));
  antitagSBOpt1Plots.push_back(plot(METall_Plot));
  antitagSBOpt1Plots.push_back(plot(METPhoton_Plot));
  antitagSBOpt1Plots.push_back(plot(MET5_Plot));
  antitagSBOpt1Plots.push_back(plot(HT_Plot));
  antitagSBOpt1Plots.push_back(plot(J1pt_Ptplot));
  antitagSBOpt1Plots.push_back(plot(J2pt_Ptplot));
  antitagSBOpt1Plots.push_back(plot(J1pt_Mplot));
  antitagSBOpt1Plots.push_back(plot(J2pt_Mplot));
  antitagSBOpt1Plots.push_back(plot(J1_M_jetBins));
  antitagSBOpt1Plots.push_back(plot(J2_M_jetBins));
  antitagSBOpt1Plots.push_back(plot(J1_deepbbtag));
  antitagSBOpt1Plots.push_back(plot(J2_deepbbtag));
  antitagSBOpt1Plots.push_back(plot(nVert));


  //Region counters
  int count_0HSR = 0; int count_0HSB = 0;
  int count_1HSR = 0; int count_1HSB = 0;
  int count_2HSR = 0; int count_2HSB = 0;

  // background MC samples - 0 lepton regions
  for (int iSample = 0; iSample < skims.ntuples.size(); iSample++) {
    RA2bTree* ntuple = skims.ntuples[iSample];

    // for (int iBin = 0; iBin < numMETbins ; iBin++) {
    //   for (int iPlot = 0; iPlot < plots[iBin].size(); iPlot++) {
    //     plots[iBin][iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
    //     plots[iBin][iPlot].setFillColor(ntuple,skims.fillColor[iSample]);
    //   }
    // }
    for (int i = 0; i < baselinePlots.size(); i++) {
      baselinePlots[i].addNtuple(ntuple,"baseline_"+skims.sampleName[iSample]);
      baselinePlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < doubletagSRPlots.size(); i++) {
      doubletagSRPlots[i].addNtuple(ntuple,"doubletagSR_"+skims.sampleName[iSample]);
      doubletagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < doubletagSBPlots.size(); i++) {
      doubletagSBPlots[i].addNtuple(ntuple,"doubletagSB_"+skims.sampleName[iSample]);
      doubletagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < tagSRPlots.size(); i++) {
      tagSRPlots[i].addNtuple(ntuple,"tagSR_"+skims.sampleName[iSample]);
      tagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < tagSBPlots.size(); i++) {
      tagSBPlots[i].addNtuple(ntuple,"tagSB_"+skims.sampleName[iSample]);
      tagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < antitagSRPlots.size(); i++) {
      antitagSRPlots[i].addNtuple(ntuple,"antitagSR_"+skims.sampleName[iSample]);
      antitagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < antitagSBPlots.size(); i++) {
      antitagSBPlots[i].addNtuple(ntuple,"antitagSB_"+skims.sampleName[iSample]);
      antitagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < antitagSROpt1Plots.size(); i++) {
      antitagSROpt1Plots[i].addNtuple(ntuple,"antitagSROpt1_"+skims.sampleName[iSample]);
      antitagSROpt1Plots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < antitagSBOpt1Plots.size(); i++) {
      antitagSBOpt1Plots[i].addNtuple(ntuple,"antitagSBOpt1_"+skims.sampleName[iSample]);
      antitagSBOpt1Plots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }

    int numEvents = ntuple->fChain->GetEntries();
    ntupleBranchStatus<RA2bTree>(ntuple); //these are set in definitions
    int bin = -1;
    double weight = 1.0;
    float trigWeight = 1.0;
    bool passBaseline;
    double jetMass1,jetMass2;

    int TotalEvents = 0;
    TString filename = ntuple->fChain->GetFile()->GetName();
    double this_lumi = 35862.824;
    if ( filename.Contains("2016") ) this_lumi = 35922.0;
    else if ( filename.Contains("2017") ) this_lumi = 41529.0;
    else if ( filename.Contains("2018") ) this_lumi = 59740.0;
    // else if ( filename.Contains("2017") ) this_lumi = 101269.0;

    // TH1* h_njetsisr;
    if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH") ) {
      TFile *fin = new TFile(filename,"READ");
      TH1F *nEventsHisto = (TH1F*)fin->Get("nEventProc");
      TotalEvents = nEventsHisto->GetBinContent(1);
    }


    for (int iEvt = 0; iEvt < numEvents ; iEvt++) {
      // if (iEvt>200000) break;
      ntuple->GetEntry(iEvt);
      if (iEvt % 100000 == 0) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << numEvents << endl;

      passBaseline=true;
      for (auto baselineCut : baselineCuts) {
        if (!passBaseline ) continue;
        passBaseline&=baselineCut(ntuple);
      }
      if (!passBaseline) continue;

      // if (( filename.Contains("SingleLept") || filename.Contains("DiLept") ) && ntuple->madHT>600. )continue;
      if (filename.Contains("MC2018")) {
        if (( filename.Contains("SingleLeptFromT_MC") || filename.Contains("SingleLeptFromTbar_MC")  || filename.Contains("DiLept_MC") ) && ntuple->GenMET>80.) continue;
      }
      else {
        if (( filename.Contains("SingleLeptFromT_MC") || filename.Contains("SingleLeptFromTbar_MC")  || filename.Contains("DiLept_MC") ) && ntuple->GenMET>150.) continue;
      }

      if (filename.Contains("T5qqqqZH") && getNumGenHiggses(ntuple)!=2) continue;

      if (region==0 || region == 1) {
        // trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc)*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
        // if (skims.sampleName[iSample]=="QCD") trigWeight=trigcorrorFakeMHT.GetCorrection(ntuple->MHT,trigunc);
        if (Year == "MC2016") trigWeight = Higfuncs::get_0l_trigeff2016(ntuple->MET,ntuple->HT);
        else if (Year == "MC2017") trigWeight = Higfuncs::get_0l_trigeff2017(ntuple->MET,ntuple->HT);
        else if (Year == "MC2018") trigWeight = Higfuncs::get_0l_trigeff2018(ntuple->MET,ntuple->HT);
        if (skims.sampleName[iSample]=="QCD")  {
          if (Year == "MC2016") trigWeight = Higfuncs::get_0l_fakemet_trigeff2016(ntuple->MET,ntuple->HT);
          else if (Year == "MC2017") trigWeight = Higfuncs::get_0l_fakemet_trigeff2017(ntuple->MET,ntuple->HT);
          else if (Year == "MC2018") trigWeight = Higfuncs::get_0l_fakemet_trigeff2018(ntuple->MET,ntuple->HT);
        }
      }

      else if (region == 2) {
        // trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc)*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
        if (Year == "MC2016") trigWeight = Higfuncs::get_1mu_trigeff2016(ntuple->MET,ntuple->HT);
        else if (Year == "MC2017") trigWeight = Higfuncs::get_1mu_trigeff2017(ntuple->MET,ntuple->HT);
        else if (Year == "MC2018") trigWeight = Higfuncs::get_1mu_trigeff2018(ntuple->MET,ntuple->HT);

      }
      else if (region == 3) {
        // trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc)*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
        if (Year == "MC2016") trigWeight = Higfuncs::get_1el_trigeff2016(ntuple->MET,ntuple->HT);
        else if (Year == "MC2017") trigWeight = Higfuncs::get_1el_trigeff2017(ntuple->MET,ntuple->HT);
        else if (Year == "MC2018") trigWeight = Higfuncs::get_1el_trigeff2018(ntuple->MET,ntuple->HT);
      }
      else if (region == 4 ) { //photon trigger weights?
         // trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc)*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      }
      else if (region == 5) {
        // trigWeight=trigcorrorFakeMHT.GetCorrection(ntuple->MHT,trigunc);
        if (Year == "MC2016") trigWeight = Higfuncs::get_0l_fakemet_trigeff2016(ntuple->MET,ntuple->HT);
        else if (Year == "MC2017") trigWeight = Higfuncs::get_0l_fakemet_trigeff2017(ntuple->MET,ntuple->HT);
        else if (Year == "MC2018") trigWeight = Higfuncs::get_0l_fakemet_trigeff2018(ntuple->MET,ntuple->HT);
      }
      else {
        trigWeight=1.0;
      }


      double isrweight = 1.0;
      if (region == 1 && filename.Contains("T5qqqqZH")) isrweight = SignalISRCorrection(ntuple); //this is for strong SUSY
      else if (region == 1 && filename.Contains("TChiHH_HToBB")) isrweight = SignalEWKISRCorrection(ntuple); //this is for EWK SUSY
      else if (filename.Contains("TTJets") && Year=="MC2016")    isrweight = SignalISRCorrection(ntuple); //this is for tt, 2016 MC fullSIM only
      else isrweight = 1.0;

      double puweight = 1.0;
      // puweight = ntuple->puWeight;

      weight = ntuple->Weight*this_lumi*trigWeight*isrweight*puweight;

      if ( filename.Contains("TChiHH_HToBB") ) weight = weight*0.5823329*0.5823329/TotalEvents;
      else if ( filename.Contains("T5qqqqZH") )  weight = weight/TotalEvents*4.0; //times 4 because only using 1/4 of events

      //Toggle whether or not we veto resolved events
      if (runVeto && resolvedBaselineCut(ntuple)) continue;

      // bin = -1;
      // for (int iBin = 0; iBin < numMETbins ; iBin++) {
      //   if (ntuple->MET > lowestMET ) {
      //     if (ntuple->MET > numMETbins*(binWidth-1)+lowestMET) bin = numMETbins-1;
      //     else bin = int((ntuple->MET-lowestMET)/binWidth);
      //   }
      // }

      // if (bin<0) continue;


      for (int i = 0; i < baselinePlots.size(); i++) {baselinePlots[i].fill(ntuple,weight);}

      // first check for 2H
      if (doubleTaggingLooseCut(ntuple)) {

        if (doubletagSRCut(ntuple)) {
          count_2HSR++;
          // plots[bin][0].fill(ntuple,weight);
          for (int i = 0; i < doubletagSRPlots.size(); i++) {doubletagSRPlots[i].fill(ntuple,weight);}
        }
        else if (doubletagSBCut(ntuple)) {
          count_2HSB++;
          // plots[bin][1].fill(ntuple,weight);
          for (int i = 0; i < doubletagSBPlots.size(); i++) {doubletagSBPlots[i].fill(ntuple,weight);}
        }
      } //end 2H

      //then check for 1H
      if (singleHiggsTagLooseCut(ntuple)) {
        if (tagSRCut( ntuple )) {
          count_1HSR++;
          // plots[bin][2].fill(ntuple,weight);
          for (int i = 0; i < tagSRPlots.size(); i++) tagSRPlots[i].fill(ntuple,weight);
        }
        else if (tagSBCut( ntuple )) {
          count_1HSB++;
          // plots[bin][3].fill(ntuple,weight);
          for (int i = 0; i < tagSBPlots.size(); i++) tagSBPlots[i].fill(ntuple,weight);
        }
      } //end 1H

      //then check for 0H
      if (antiTaggingLooseCut(ntuple)) {
        if (antitagSRCut( ntuple )) {
          count_0HSR++;
          // plots[bin][4].fill(ntuple,weight);
          for (int i = 0; i < antitagSRPlots.size(); i++) antitagSRPlots[i].fill(ntuple,weight);
          if (antitagSRCut_opt1( ntuple )) for (int i = 0; i < antitagSROpt1Plots.size(); i++) antitagSROpt1Plots[i].fill(ntuple,weight);
        }
        else {
          if (antitagSBCut( ntuple )) {
            count_0HSB++;
            // plots[bin][5].fill(ntuple,weight);
            for (int i = 0; i < antitagSBPlots.size(); i++) antitagSBPlots[i].fill(ntuple,weight);
            if (antitagSBCut_opt1( ntuple )) for (int i = 0; i < antitagSBOpt1Plots.size(); i++) antitagSBOpt1Plots[i].fill(ntuple,weight);

          }
        }
      } //end 0H region
    } // end event loop
    outputFile->cd();
  } // end sample loop


  // Begin data
  if (region!=0 && region!=1 && runData) {
    RA2bTree* ntuple = skims.dataNtuple;
    // for (int iBin = 0; iBin < numMETbins ; iBin++) {
    //   for (int iPlot = 0; iPlot < plots[iBin].size(); iPlot++) {
    //     plots[iBin][iPlot].addDataNtuple(ntuple,"data");
    //   }
    // }
    for (int i = 0; i < baselinePlots.size(); i++) {
      baselinePlots[i].addDataNtuple(ntuple,"baseline_data");
    }
    for (int i = 0; i < doubletagSRPlots.size(); i++) {
      doubletagSRPlots[i].addDataNtuple(ntuple,"doubletagSR_data");
    }
    for (int i = 0; i < doubletagSBPlots.size(); i++) {
      doubletagSBPlots[i].addDataNtuple(ntuple,"doubletagSB_data");
    }
    for (int i = 0; i < tagSRPlots.size(); i++) {
      tagSRPlots[i].addDataNtuple(ntuple,"tagSR_data");
    }
    for (int i = 0; i < tagSBPlots.size(); i++) {
      tagSBPlots[i].addDataNtuple(ntuple,"tagSB_data");
    }
    for (int i = 0; i < antitagSRPlots.size(); i++) {
      antitagSRPlots[i].addDataNtuple(ntuple,"antitagSR_data");
    }
    for (int i = 0; i < antitagSBPlots.size(); i++) {
      antitagSBPlots[i].addDataNtuple(ntuple,"antitagSB_data");
    }
    for (int i = 0; i < antitagSROpt1Plots.size(); i++) {
      antitagSROpt1Plots[i].addDataNtuple(ntuple,"antitagSROpt1_data");
    }
    for (int i = 0; i < antitagSBOpt1Plots.size(); i++) {
      antitagSBOpt1Plots[i].addDataNtuple(ntuple,"antitagSBOpt1_data");
    }

    int numEvents = ntuple->fChain->GetEntries();
    TString filename = ntuple->fChain->GetFile()->GetName();

    ntupleBranchStatus<RA2bTree>(ntuple);

    bool passBaseline;
    double jetMass1,jetMass2;
    for (int iEvt = 0; iEvt < numEvents; iEvt++) {
      ntuple->GetEntry(iEvt);
      if (iEvt % 100000 == 0 ) cout << "data: " << iEvt << "/" << numEvents << endl;


      passBaseline=true;
      for (auto baselineCut : baselineCuts) {
        if (!passBaseline ) continue;
        passBaseline&=baselineCut(ntuple);
      }
      if (!passBaseline) continue;

      // double weight = 1.0;
      // double this_lumi = 35862.824;
      // if ( filename.Contains("2016") ) this_lumi = 35922.0;
      // else if ( filename.Contains("2017") ) this_lumi = 41529.0;
      // else if ( filename.Contains("2018") ) this_lumi = 59740.0;
      // // else if ( filename.Contains("2017") ) this_lumi = 101269.0;
      //
      // weight = ntuple->Weight*this_lumi;

      //Toggle whether or not we veto resolved events
      if (runVeto && resolvedBaselineCut(ntuple)) continue;

      // int bin = -1;
      // for (int iBin = 0; iBin < numMETbins; iBin++) {
      //   if (ntuple->MET > lowestMET) {
      //     if (ntuple->MET > lowestMET+binWidth*(numMETbins-1) )
      //       bin = numMETbins-1;
      //     else
      //       bin = int((ntuple->MET-lowestMET)/binWidth);
      //   }
      // }
      // if (bin<0) continue;

      for (int i = 0; i < baselinePlots.size(); i++) {baselinePlots[i].fillData(ntuple);}

      // first check for 2H
      if (doubleTaggingLooseCut(ntuple) && region!=1) {
        if (doubletagSRCut(ntuple) && region!=0) {
          // plots[bin][0].fillData(ntuple);
          for (int i = 0; i < doubletagSRPlots.size(); i++) {doubletagSRPlots[i].fillData(ntuple);}
        }
        else if (doubletagSBCut(ntuple)) {
          for (int i = 0; i < doubletagSBPlots.size(); i++) {doubletagSBPlots[i].fillData(ntuple);}
        }
      } //end 2H

      //then check for 1H
      if (singleHiggsTagLooseCut(ntuple) && region!=1) {
        if (tagSRCut( ntuple ) && region!=0) {
          for (int i = 0; i < tagSRPlots.size(); i++) tagSRPlots[i].fillData(ntuple);
        }
        else if (tagSBCut( ntuple )) {
          for (int i = 0; i < tagSBPlots.size(); i++) tagSBPlots[i].fillData(ntuple);
        }
      } //end 1H

      //then check for 0H
      if (antiTaggingLooseCut(ntuple) && region!=1) {
        if (antitagSRCut( ntuple )) {
          for (int i = 0; i < antitagSRPlots.size(); i++) antitagSRPlots[i].fillData(ntuple);
          if (antitagSRCut_opt1( ntuple )) for (int i = 0; i < antitagSROpt1Plots.size(); i++) antitagSROpt1Plots[i].fillData(ntuple);
        }
        else {
          if (antitagSBCut( ntuple )) {
            for (int i = 0; i < antitagSBPlots.size(); i++) antitagSBPlots[i].fillData(ntuple);
            if (antitagSBCut_opt1( ntuple )) for (int i = 0; i < antitagSBOpt1Plots.size(); i++) antitagSBOpt1Plots[i].fillData(ntuple);
          }
        }
      } //end 0H region
    } // end event loop
  } // end not signal, and ifData


  bool sumBkgs = true;
  if (region == 1) sumBkgs = false;
  if (runData) sumBkgs = false;


  //this just saves a TON of plots per MET bin... maybe these are useful?
  // for (int iBin = 0; iBin < numMETbins; iBin++) {
  //   for (int iPlot = 0; iPlot < plots[iBin].size(); iPlot++) {
  //     outputFile->cd();
  //     plots[iBin][iPlot].buildSum();
  //     plots[iBin][iPlot].Write();
  //     if (sumBkgs) plots[iBin][iPlot].sum->Write();
  //   }
  // }

  for (int i = 0; i < baselinePlots.size(); i++) {
    outputFile->cd();
    baselinePlots[i].Write();
    if (sumBkgs) {
      baselinePlots[i].buildSum("baseline");
      baselinePlots[i].sum->Write();
    }
  }

  for (int i = 0; i < doubletagSRPlots.size(); i++) {
    outputFile->cd();
    doubletagSRPlots[i].Write();
    if (sumBkgs) {
      doubletagSRPlots[i].buildSum("doubletagSR");
      doubletagSRPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < doubletagSBPlots.size(); i++) {
    outputFile->cd();
    doubletagSBPlots[i].Write();
    if (sumBkgs) {
      doubletagSBPlots[i].buildSum("doubletagSB");
      doubletagSBPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < tagSRPlots.size(); i++) {
    outputFile->cd();
    tagSRPlots[i].Write();
    if (sumBkgs) {
      tagSRPlots[i].buildSum("tagSR");
      tagSRPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < tagSBPlots.size(); i++) {
    outputFile->cd();
    tagSBPlots[i].Write();
    if (sumBkgs) {
      tagSBPlots[i].buildSum("tagSB");
      tagSBPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < antitagSRPlots.size(); i++) {
    outputFile->cd();
    antitagSRPlots[i].Write();
    if (sumBkgs) {
      antitagSRPlots[i].buildSum("antitagSR");
      antitagSRPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < antitagSBPlots.size(); i++) {
    outputFile->cd();
    antitagSBPlots[i].Write();
    if (sumBkgs) {
      antitagSBPlots[i].buildSum("antitagSB");
      antitagSBPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < antitagSROpt1Plots.size(); i++) {
    outputFile->cd();
    antitagSROpt1Plots[i].Write();
    if (sumBkgs) {
      antitagSROpt1Plots[i].buildSum("antitagSROpt1");
      antitagSROpt1Plots[i].sum->Write();
    }
  }
  for (int i = 0; i < antitagSBOpt1Plots.size(); i++) {
    outputFile->cd();
    antitagSBOpt1Plots[i].Write();
    if (sumBkgs) {
      antitagSBOpt1Plots[i].buildSum("antitagSBOpt1");
      antitagSBOpt1Plots[i].sum->Write();
    }
  }
  outputFile->Close();

  std::cout<<"# unscaled events: (w/ PU weights) for "<< Year <<": \n";
  std::cout<<"2HSR: "<< count_2HSR <<", 2HSB: "<< count_2HSB <<", 1HSR: "<< count_1HSR <<", 1HSB: "<< count_1HSB <<", 0HSR: "<< count_0HSR <<", 0HSB: "<< count_0HSB<<std::endl;
}
