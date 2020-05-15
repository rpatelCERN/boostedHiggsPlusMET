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
#include "TriggerCorrector.h"

using namespace std;
using namespace alphabet;

int main(int argc, char** argv) {

  int region(0);
  bool looseCuts(false);
  int MAX_EVENTS(99999999);
  int runVeto(0);
  bool runData = false;

  region = atoi(argv[1]);

  gROOT->ProcessLine(".L tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  TString Year(argv[2]);
  runVeto = atoi(argv[3]);
  cout <<"region: "<<region<<", veto?: "<<runVeto<<endl;


  TriggerCorrector trigcorror; TriggerCorrector trigcorrorHT; TriggerCorrector trigcorrorFakeMHT;
  trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2016");
  trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2016");
  trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2016");


  skimSamples* skims_;
  if (region == 0 ) skims_ = new skimSamples(skimSamples::kSignal, Year);
  else if (region == 1 ) {
    int NLSP=atoi(argv[4]);
    int LSP=atoi(argv[5]);
    std::cout<<" Parent "<<NLSP <<" LSP "<<LSP<<std::endl;
    skims_ = new skimSamples(skimSamples::kSignalOnly, Year,NLSP,LSP);
  }
  else if (region == 2 ) skims_ = new skimSamples(skimSamples::kSLm, Year);
  else if (region == 3 ) skims_ = new skimSamples(skimSamples::kSLe, Year);
  else if (region == 4 ) skims_ = new skimSamples(skimSamples::kPhoton, Year);
  else if (region == 5 ) skims_ = new skimSamples(skimSamples::kLowDphi, Year);
  else assert(1);

  typedef bool(*cuts)(RA2bTree*);
  vector<cuts> baselineCuts;


  if (region == 0 || region == 1) baselineCuts.push_back(*boostedBaselineCut<RA2bTree>);
  else if (region == 2) baselineCuts.push_back(*singleMuBaselineCut<RA2bTree>);
  else if (region == 3) baselineCuts.push_back(*singleEleBaselineCut<RA2bTree>);
  else if (region == 4) baselineCuts.push_back(*photonBaselineCut<RA2bTree>);
  else if (region == 5) baselineCuts.push_back(*lowDphiBaselineCut<RA2bTree>);
  else assert(1);

  skimSamples skims = *skims_;

  typedef plot<RA2bTree> plot;
  vector<vector<plot> > plots;

  double mJbins[4]={50.,85.,135.,250.};

  double METbinsPhoton[5]={100.,300.,500.,700.,1400.};
  double METbinsLow[5]={250.,300.,500.,700.,1400.};
  double METbins[4]={300.,500.,700.,1400.};

  TString METbinsPhoton_string[5]={"100","300","500","700","1400"};
  TString METbinsLow_string[5]={"250","300","500","700","1400"};
  TString METbins_string[4]={"300","500","700","1400"};

  // for (int i = 0; i < 3 ; i++) {
  for (int i = 0; i < numMETbins_V12 ; i++) {
    TString tag="_";
    tag+=METbins_string[i];
    // tag+=lowestMET+i*binWidth;
    vector<plot> plotsTemp;
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_doubletagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_doubletagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_tagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_tagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSROpt1"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSROpt2"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSBOpt1"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSBOpt2"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_notagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_notagSB"+tag,"m_{J} [GeV]",3,mJbins));

    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_doubletagSRLead"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSRLead"+tag,"m_{J} [GeV]",3,mJbins));

    plots.push_back(plotsTemp);
  }

  //vector<plot> tempPlots;
  // plot MET_Plot(*fillMET<RA2bTree>,"MET250","m_{J} [GeV]",4,METbinsLow);
  plot MET_Plot(*fillMET<RA2bTree>,"MET","m_{J} [GeV]",3,METbins);

  plot J1pt_Ptplot(*fillLeadingJetPt<RA2bTree>,"J1pt_Pt","p_{T,J} [GeV]",50,300.,1300.);
  plot J1pt_Mplot(*fillLeadingJetMass<RA2bTree>,"J1pt_M","m_{J} [GeV]",80,50.,250.);
  plot J1_M_jetBins(*fillLeadingJetMass<RA2bTree>,"J1_M_jetBins","m_{J} [GeV]",3,mJbins);

  plot J2pt_Ptplot(*fillSubLeadingJetPt<RA2bTree>,"J2pt_Pt","p_{T,J} [GeV]",50,300.,1300.);
  plot J2pt_Mplot(*fillSubLeadingJetMass<RA2bTree>,"J2pt_M","m_{J} [GeV]",80,50.,250.);
  plot J2_M_jetBins(*fillSubLeadingJetMass<RA2bTree>,"J2_M_jetBins","m_{J} [GeV]",3,mJbins);

  plot J1_deepbbtag(*fillLeadingdeepBBtag<RA2bTree>,"LeadDeepBBTag","Lead jet (unsupp) deep bb-tag",50,0.0,1.0);
  plot J2_deepbbtag(*fillSubLeadingdeepBBtag<RA2bTree>,"SubLeadDeepBBTag","Sub-lead jet (unsupp) deep bb-tag",50,0.0,1.0);

  // TH2D * MET_v_photonPT;
  //plots where lead jet is in SR and sublead is anywhere
  //Need both for 2H and 0H

  vector<plot> doubletagSRPlots;
  doubletagSRPlots.push_back(plot(MET_Plot));
  doubletagSRPlots.push_back(plot(J1pt_Ptplot));
  doubletagSRPlots.push_back(plot(J2pt_Ptplot));
  doubletagSRPlots.push_back(plot(J1pt_Mplot));
  doubletagSRPlots.push_back(plot(J2pt_Mplot));
  doubletagSRPlots.push_back(plot(J1_M_jetBins));
  doubletagSRPlots.push_back(plot(J2_M_jetBins));
  doubletagSRPlots.push_back(plot(J1_deepbbtag));
  doubletagSRPlots.push_back(plot(J2_deepbbtag));


  vector<plot> doubletagSBPlots;
  doubletagSBPlots.push_back(plot(MET_Plot));
  doubletagSBPlots.push_back(plot(J1pt_Ptplot));
  doubletagSBPlots.push_back(plot(J2pt_Ptplot));
  doubletagSBPlots.push_back(plot(J1pt_Mplot));
  doubletagSBPlots.push_back(plot(J2pt_Mplot));
  doubletagSBPlots.push_back(plot(J1_M_jetBins));
  doubletagSBPlots.push_back(plot(J2_M_jetBins));
  doubletagSBPlots.push_back(plot(J1_deepbbtag));
  doubletagSBPlots.push_back(plot(J2_deepbbtag));


  vector<plot> tagSRPlots;
  tagSRPlots.push_back(plot(MET_Plot));
  tagSRPlots.push_back(plot(J1pt_Ptplot));
  tagSRPlots.push_back(plot(J2pt_Ptplot));
  tagSRPlots.push_back(plot(J1pt_Mplot));
  tagSRPlots.push_back(plot(J2pt_Mplot));
  tagSRPlots.push_back(plot(J1_M_jetBins));
  tagSRPlots.push_back(plot(J2_M_jetBins));
  tagSRPlots.push_back(plot(J1_deepbbtag));
  tagSRPlots.push_back(plot(J2_deepbbtag));


  vector<plot> tagSBPlots;
  tagSBPlots.push_back(plot(MET_Plot));
  tagSBPlots.push_back(plot(J1pt_Ptplot));
  tagSBPlots.push_back(plot(J2pt_Ptplot));
  tagSBPlots.push_back(plot(J1pt_Mplot));
  tagSBPlots.push_back(plot(J2pt_Mplot));
  tagSBPlots.push_back(plot(J1_M_jetBins));
  tagSBPlots.push_back(plot(J2_M_jetBins));
  tagSBPlots.push_back(plot(J1_deepbbtag));
  tagSBPlots.push_back(plot(J2_deepbbtag));


  vector<plot> antitagSRPlots;
  antitagSRPlots.push_back(plot(MET_Plot));
  antitagSRPlots.push_back(plot(J1pt_Ptplot));
  antitagSRPlots.push_back(plot(J2pt_Ptplot));
  antitagSRPlots.push_back(plot(J1pt_Mplot));
  antitagSRPlots.push_back(plot(J2pt_Mplot));
  antitagSRPlots.push_back(plot(J1_M_jetBins));
  antitagSRPlots.push_back(plot(J2_M_jetBins));
  antitagSRPlots.push_back(plot(J1_deepbbtag));
  antitagSRPlots.push_back(plot(J2_deepbbtag));


  vector<plot> antitagSBPlots;
  antitagSBPlots.push_back(plot(MET_Plot));
  antitagSBPlots.push_back(plot(J1pt_Ptplot));
  antitagSBPlots.push_back(plot(J2pt_Ptplot));
  antitagSBPlots.push_back(plot(J1pt_Mplot));
  antitagSBPlots.push_back(plot(J2pt_Mplot));
  antitagSBPlots.push_back(plot(J1_M_jetBins));
  antitagSBPlots.push_back(plot(J2_M_jetBins));
  antitagSBPlots.push_back(plot(J1_deepbbtag));
  antitagSBPlots.push_back(plot(J2_deepbbtag));

  vector<plot> antitagSROpt1Plots;
  antitagSROpt1Plots.push_back(plot(MET_Plot));
  antitagSROpt1Plots.push_back(plot(J1pt_Ptplot));
  antitagSROpt1Plots.push_back(plot(J2pt_Ptplot));
  antitagSROpt1Plots.push_back(plot(J1pt_Mplot));
  antitagSROpt1Plots.push_back(plot(J2pt_Mplot));
  antitagSROpt1Plots.push_back(plot(J1_M_jetBins));
  antitagSROpt1Plots.push_back(plot(J2_M_jetBins));
  antitagSROpt1Plots.push_back(plot(J1_deepbbtag));
  antitagSROpt1Plots.push_back(plot(J2_deepbbtag));

  vector<plot> antitagSBOpt1Plots;
  antitagSBOpt1Plots.push_back(plot(MET_Plot));
  antitagSBOpt1Plots.push_back(plot(J1pt_Ptplot));
  antitagSBOpt1Plots.push_back(plot(J2pt_Ptplot));
  antitagSBOpt1Plots.push_back(plot(J1pt_Mplot));
  antitagSBOpt1Plots.push_back(plot(J2pt_Mplot));
  antitagSBOpt1Plots.push_back(plot(J1_M_jetBins));
  antitagSBOpt1Plots.push_back(plot(J2_M_jetBins));
  antitagSBOpt1Plots.push_back(plot(J1_deepbbtag));
  antitagSBOpt1Plots.push_back(plot(J2_deepbbtag));

  vector<plot> antitagSROpt2Plots;
  antitagSROpt2Plots.push_back(plot(MET_Plot));
  antitagSROpt2Plots.push_back(plot(J1pt_Ptplot));
  antitagSROpt2Plots.push_back(plot(J2pt_Ptplot));
  antitagSROpt2Plots.push_back(plot(J1pt_Mplot));
  antitagSROpt2Plots.push_back(plot(J2pt_Mplot));
  antitagSROpt2Plots.push_back(plot(J1_M_jetBins));
  antitagSROpt2Plots.push_back(plot(J2_M_jetBins));
  antitagSROpt2Plots.push_back(plot(J1_deepbbtag));
  antitagSROpt2Plots.push_back(plot(J2_deepbbtag));

  vector<plot> antitagSBOpt2Plots;
  antitagSBOpt2Plots.push_back(plot(MET_Plot));
  antitagSBOpt2Plots.push_back(plot(J1pt_Ptplot));
  antitagSBOpt2Plots.push_back(plot(J2pt_Ptplot));
  antitagSBOpt2Plots.push_back(plot(J1pt_Mplot));
  antitagSBOpt2Plots.push_back(plot(J2pt_Mplot));
  antitagSBOpt2Plots.push_back(plot(J1_M_jetBins));
  antitagSBOpt2Plots.push_back(plot(J2_M_jetBins));
  antitagSBOpt2Plots.push_back(plot(J1_deepbbtag));
  antitagSBOpt2Plots.push_back(plot(J2_deepbbtag));

  vector<plot> notagSRPlots;
  notagSRPlots.push_back(plot(MET_Plot));
  notagSRPlots.push_back(plot(J1pt_Ptplot));
  notagSRPlots.push_back(plot(J2pt_Ptplot));
  notagSRPlots.push_back(plot(J1pt_Mplot));
  notagSRPlots.push_back(plot(J2pt_Mplot));
  notagSRPlots.push_back(plot(J1_M_jetBins));
  notagSRPlots.push_back(plot(J2_M_jetBins));
  notagSRPlots.push_back(plot(J1_deepbbtag));
  notagSRPlots.push_back(plot(J2_deepbbtag));


  vector<plot> notagSBPlots;
  notagSBPlots.push_back(plot(MET_Plot));
  notagSBPlots.push_back(plot(J1pt_Ptplot));
  notagSBPlots.push_back(plot(J2pt_Ptplot));
  notagSBPlots.push_back(plot(J1pt_Mplot));
  notagSBPlots.push_back(plot(J2pt_Mplot));
  notagSBPlots.push_back(plot(J1_M_jetBins));
  notagSBPlots.push_back(plot(J2_M_jetBins));
  notagSBPlots.push_back(plot(J1_deepbbtag));
  notagSBPlots.push_back(plot(J2_deepbbtag));


  vector<plot> doubletagSRLeadPlots;
  doubletagSRLeadPlots.push_back(plot(MET_Plot));
  doubletagSRLeadPlots.push_back(plot(J1pt_Ptplot));
  doubletagSRLeadPlots.push_back(plot(J2pt_Ptplot));
  doubletagSRLeadPlots.push_back(plot(J1pt_Mplot));
  doubletagSRLeadPlots.push_back(plot(J2pt_Mplot));
  doubletagSRLeadPlots.push_back(plot(J1_M_jetBins));
  doubletagSRLeadPlots.push_back(plot(J2_M_jetBins));
  doubletagSRLeadPlots.push_back(plot(J1_deepbbtag));
  doubletagSRLeadPlots.push_back(plot(J2_deepbbtag));

  vector<plot> antitagSRLeadPlots;
  antitagSRLeadPlots.push_back(plot(MET_Plot));
  antitagSRLeadPlots.push_back(plot(J1pt_Ptplot));
  antitagSRLeadPlots.push_back(plot(J2pt_Ptplot));
  antitagSRLeadPlots.push_back(plot(J1pt_Mplot));
  antitagSRLeadPlots.push_back(plot(J2pt_Mplot));
  antitagSRLeadPlots.push_back(plot(J1_M_jetBins));
  antitagSRLeadPlots.push_back(plot(J2_M_jetBins));
  antitagSRLeadPlots.push_back(plot(J1_deepbbtag));
  antitagSRLeadPlots.push_back(plot(J2_deepbbtag));

  vector<plot> doubletagSBLeadPlots;
  doubletagSBLeadPlots.push_back(plot(MET_Plot));
  doubletagSBLeadPlots.push_back(plot(J1pt_Ptplot));
  doubletagSBLeadPlots.push_back(plot(J2pt_Ptplot));
  doubletagSBLeadPlots.push_back(plot(J1pt_Mplot));
  doubletagSBLeadPlots.push_back(plot(J2pt_Mplot));
  doubletagSBLeadPlots.push_back(plot(J1_M_jetBins));
  doubletagSBLeadPlots.push_back(plot(J2_M_jetBins));
  doubletagSBLeadPlots.push_back(plot(J1_deepbbtag));
  doubletagSBLeadPlots.push_back(plot(J2_deepbbtag));

  vector<plot> antitagSBLeadPlots;
  antitagSBLeadPlots.push_back(plot(MET_Plot));
  antitagSBLeadPlots.push_back(plot(J1pt_Ptplot));
  antitagSBLeadPlots.push_back(plot(J2pt_Ptplot));
  antitagSBLeadPlots.push_back(plot(J1pt_Mplot));
  antitagSBLeadPlots.push_back(plot(J2pt_Mplot));
  antitagSBLeadPlots.push_back(plot(J1_M_jetBins));
  antitagSBLeadPlots.push_back(plot(J2_M_jetBins));
  antitagSBLeadPlots.push_back(plot(J1_deepbbtag));
  antitagSBLeadPlots.push_back(plot(J2_deepbbtag));

  TFile* outputFile;
  TString regionName;
  TString cutName="";
  TString METcut="";

  if (looseCuts ) cutName="_looseCuts";
  if (region == 0 ) regionName="";
  if (region == 1 ){
    // regionName="_signalOnly";
    int NLSP=atoi(argv[3]);
    int LSP=atoi(argv[4]);
    regionName=TString::Format("_TChiHH%d_LSP%d", atoi(argv[3]), atoi(argv[4]));
  }
  if (region == 2 ) regionName="_singleMu";
  if (region == 3 ) regionName="_singleEle";
  if (region == 4 ) regionName="_photon";
  if (region == 5 ) regionName="_lowDphi";

  if (runVeto) cutName = "_resVeto";



  // METcut = "_MET250";
  METcut = "_MET300";
  if (runData) METcut = "_data";
  // regionName = "_MET250";
  // regionName = "_TTgenMET150";


  // outputFile = new TFile("ALPHABET_ZJets16_withWeights.root","RECREATE");
  outputFile = new TFile("ALPHABET"+Year+"_V17"+regionName+METcut+cutName+"_2BoostedH.root","RECREATE");
  // outputFile = new TFile("ALPHABETBoost_V17_signalOnly.root","RECREATE");
  // outputFile = new TFile("ALPHABETBoost_MC2016_V17_oldBBTag.root","RECREATE");

  // background MC samples - 0 lepton regions
  // for (int iSample = 0; iSample < 0; iSample++) {
  for (int iSample = 0; iSample < skims.ntuples.size(); iSample++) {
    RA2bTree* ntuple = skims.ntuples[iSample];

    // MET_v_photonPT = new TH2D ("MET_v_photonPT_"+skims.sampleName[iSample],";MET; photon pT",40,100,900,40,100,900);


    for (int iBin = 0; iBin < numMETbins ; iBin++) {
      for (int iPlot = 0; iPlot < plots[iBin].size(); iPlot++) {
        plots[iBin][iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
        plots[iBin][iPlot].setFillColor(ntuple,skims.fillColor[iSample]);
      }
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
    for (int i = 0; i < antitagSROpt2Plots.size(); i++) {
      antitagSROpt2Plots[i].addNtuple(ntuple,"antitagSROpt2_"+skims.sampleName[iSample]);
      antitagSROpt2Plots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < antitagSBOpt2Plots.size(); i++) {
      antitagSBOpt2Plots[i].addNtuple(ntuple,"antitagSBOpt2_"+skims.sampleName[iSample]);
      antitagSBOpt2Plots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < notagSRPlots.size(); i++) {
      notagSRPlots[i].addNtuple(ntuple,"notagSR_"+skims.sampleName[iSample]);
      notagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < notagSBPlots.size(); i++) {
      notagSBPlots[i].addNtuple(ntuple,"notagSB_"+skims.sampleName[iSample]);
      notagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < doubletagSRLeadPlots.size(); i++) {
      doubletagSRLeadPlots[i].addNtuple(ntuple,"doubletagSRLead_"+skims.sampleName[iSample]);
      doubletagSRLeadPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < antitagSRLeadPlots.size(); i++) {
      antitagSRLeadPlots[i].addNtuple(ntuple,"antitagSRLead_"+skims.sampleName[iSample]);
      antitagSRLeadPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < doubletagSBLeadPlots.size(); i++) {
      doubletagSBLeadPlots[i].addNtuple(ntuple,"doubletagSBLead_"+skims.sampleName[iSample]);
      doubletagSBLeadPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0; i < antitagSBLeadPlots.size(); i++) {
      antitagSBLeadPlots[i].addNtuple(ntuple,"antitagSBLead_"+skims.sampleName[iSample]);
      antitagSBLeadPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }

    int numEvents = ntuple->fChain->GetEntries();
    ntupleBranchStatus<RA2bTree>(ntuple); //these are set in definitions
    int bin = -1;
    double weight=0.;
    float trigWeight=1.0;
    float trigunc=0;
    bool passBaseline;
    double jetMass1,jetMass2;

    int TotalEvents = 0;
    TString filename;
    filename = ntuple->fChain->GetFile()->GetName();
    double this_lumi = 35862.824;
    if ( filename.Contains("2016") ) this_lumi = 35922.0;
    else if ( filename.Contains("2017") ) this_lumi = 41529.0;
    else if ( filename.Contains("2018") ) this_lumi = 59740.0;
    // else if ( filename.Contains("2017") ) this_lumi = 101269.0;

    if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH") ) {
      this_lumi = 137191.0;
      TFile *fin = new TFile(filename,"READ");
      TH1F *nEventsHisto = (TH1F*)fin->Get("nEventProc");
      TotalEvents = nEventsHisto->GetBinContent(1);
    }

    // for (int iEvt = 0; iEvt < 0 ; iEvt++) {
    for (int iEvt = 0; iEvt < numEvents ; iEvt++) {
      ntuple->GetEntry(iEvt);
      if (iEvt % 100000 == 0 ) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << min(MAX_EVENTS,numEvents) << endl;
      if (region==0 || region == 1) {
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
        if(skims.sampleName[iSample]=="QCD") trigWeight=trigcorrorFakeMHT.GetCorrection(ntuple->MHT,trigunc);
      }
      else if (region == 2 ) {
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      }
      else if (region == 3 ) {
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      }
      // else if (region == 4 ) { //check this
      //   trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
      //   trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      // }
      else if (region == 5 ) {
        if(skims.sampleName[iSample]!="QCD") {
          trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
          trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
        }
        else {
          trigWeight=1.0;
        }
      }

      passBaseline=true;
      for (auto baselineCut : baselineCuts ) {
        if (! passBaseline ) continue;
        passBaseline&=baselineCut(ntuple);
      }

      if (! passBaseline ) continue;
      if (ntuple->MET<300.) continue; //MET cut

      // if (( filename.Contains("SingleLept") || filename.Contains("DiLept") ) && ntuple->madHT>600. )continue;
      if (filename.Contains("MC2018")) {
        if (( filename.Contains("SingleLeptFromT_MC") || filename.Contains("SingleLeptFromTbar_MC")  || filename.Contains("DiLept_MC") ) && ntuple->GenMET>80.) continue;
      }
      else {
        if (( filename.Contains("SingleLeptFromT_MC") || filename.Contains("SingleLeptFromTbar_MC")  || filename.Contains("DiLept_MC") ) && ntuple->GenMET>150.) continue;
      }

      // weight = ntuple->Weight*this_lumi;
      weight = ntuple->Weight*this_lumi*trigWeight;

      // MET_v_photonPT->Fill(ntuple->MET,ntuple->Photons->at(0).Pt(),weight);
      if ( filename.Contains("TChiHH_HToBB") ) weight = weight*0.5823329*0.5823329/TotalEvents;
      else if ( filename.Contains("T5qqqqZH") )  weight = weight/TotalEvents*4.0; //times 4 because only using 1/4 of events

      //if (skims.sampleName[iSample] == "TT" ) {
      //    weight *= ISRweights(ntuple);
      //}

      //Toggle whether or not we veto resolved events
      if (runVeto){
        if (resolvedBaselineCut(ntuple) ) continue;
      }

      bin = -1;
      // //our newest MET bins would be 250-300, 300-500, 500-700, and 700+
      // float thisMET_forBin = ntuple->MET;
      // // if (thisMET_forBin>100.&&thisMET_forBin<=300.) bin=0;
      // if (thisMET_forBin>250.&&thisMET_forBin<=300.) bin=0;
      // else if (thisMET_forBin>300.&&thisMET_forBin<=500.) bin=1;
      // else if (thisMET_forBin>500.&&thisMET_forBin<=700.) bin=2;
      // else if (thisMET_forBin>700.) bin=3;

      for (int iBin = 0; iBin < numMETbins_V12 ; iBin++) {
        if (ntuple->MET > lowestMET_V12 ) {
          if (ntuple->MET > numMETbins_V12*(binWidth_V12-1)+lowestMET_V12) bin = numMETbins_V12-1;
          else bin = int((ntuple->MET-lowestMET_V12)/binWidth_V12);
        }
      }

      if (bin < 0 ) continue;

      // No-tag plots, used for MET shape
      if (doubleMassCut(ntuple) ) { //this requires both lead and sublead jet to be in the mass window
        plots[bin][10].fill(ntuple,weight);
        for (int i = 0; i < notagSRPlots.size(); i++) {
          notagSRPlots[i].fill(ntuple,weight);
        }
      }
      else { //one of the jets failed the SR mass cut
        if (SBMassCut(ntuple) ) { //either J1 or J2 are in SB (or both)
          plots[bin][11].fill(ntuple,weight);
          for (int i = 0; i < notagSBPlots.size(); i++) {
            notagSBPlots[i].fill(ntuple,weight);
          }
        }
      }


      // Require lead jet to be in the mass SR, but sublead can be anywhere after baseline cuts
      // Then divide by 2H and 0H
      if (isMassSR(ntuple,0)) {
        if (doubleTaggingLooseCut(ntuple)) {
          plots[bin][12].fill(ntuple,weight);
          for (int i = 0; i < doubletagSRLeadPlots.size(); i++) {doubletagSRLeadPlots[i].fill(ntuple,weight);}
        }
        else if (antiTaggingLooseCut(ntuple)) {
          plots[bin][13].fill(ntuple,weight);
          for (int i = 0; i < antitagSRLeadPlots.size(); i++) {antitagSRLeadPlots[i].fill(ntuple,weight);}
        }
      }

      if (isMassSR(ntuple,1)) {
        if (doubleTaggingLooseCut(ntuple)) {
          for (int i = 0; i < doubletagSBLeadPlots.size(); i++) {doubletagSBLeadPlots[i].fill(ntuple,weight);}
        }
        else if (antiTaggingLooseCut(ntuple)) {
          for (int i = 0; i < antitagSBLeadPlots.size(); i++) {antitagSBLeadPlots[i].fill(ntuple,weight);}
        }
      }

      // first check for 2H
      if (doubleTaggingLooseCut(ntuple)) {
        if (doubletagSRCut(ntuple)) {
          plots[bin][0].fill(ntuple,weight);
          for (int i = 0; i < doubletagSRPlots.size(); i++) {doubletagSRPlots[i].fill(ntuple,weight);}
        }
        else if (doubletagSBCut(ntuple)) {
          plots[bin][1].fill(ntuple,weight);
          for (int i = 0; i < doubletagSBPlots.size(); i++) {doubletagSBPlots[i].fill(ntuple,weight);}
        }
      } //end 2H

      //then check for 1H
      if (singleHiggsTagLooseCut(ntuple)) {
        if (tagSRCut( ntuple )) {
          plots[bin][2].fill(ntuple,weight);
          for (int i = 0; i < tagSRPlots.size(); i++) tagSRPlots[i].fill(ntuple,weight);
        }
        else if (tagSBCut( ntuple )) {
          plots[bin][3].fill(ntuple,weight);
          for (int i = 0; i < tagSBPlots.size(); i++) tagSBPlots[i].fill(ntuple,weight);
        }
      } //end 1H

      //then check for 0H
      if (antiTaggingLooseCut(ntuple)) {
        if (antitagSRCut( ntuple )) {
          plots[bin][4].fill(ntuple,weight);
          for (int i = 0; i < antitagSRPlots.size(); i++) antitagSRPlots[i].fill(ntuple,weight);
          if (antitagSRCut_opt1( ntuple )) {
            plots[bin][6].fill(ntuple,weight);
            for (int i = 0; i < antitagSROpt1Plots.size(); i++) antitagSROpt1Plots[i].fill(ntuple,weight);
          }
          if (antitagSRCut_opt2( ntuple )) {
            plots[bin][7].fill(ntuple,weight);
            for (int i = 0; i < antitagSROpt2Plots.size(); i++) antitagSROpt2Plots[i].fill(ntuple,weight);
          }
        }
        else {
          if (antitagSBCut( ntuple )) {
            plots[bin][5].fill(ntuple,weight);
            for (int i = 0; i < antitagSBPlots.size(); i++) antitagSBPlots[i].fill(ntuple,weight);
            if (antitagSBCut_opt1( ntuple )) {
              plots[bin][8].fill(ntuple,weight);
              for (int i = 0; i < antitagSBOpt1Plots.size(); i++) antitagSBOpt1Plots[i].fill(ntuple,weight);
            }
            if (antitagSBCut_opt2( ntuple )) {
              plots[bin][9].fill(ntuple,weight);
              for (int i = 0; i < antitagSBOpt2Plots.size(); i++) antitagSBOpt2Plots[i].fill(ntuple,weight);
            }
          }
        }
      } //end 0H region

      // if (doubletagSRCut_V12(ntuple)) {for (int i = 0; i < doubletagSRPlots.size(); i++) {doubletagSRPlots[i].fill(ntuple,weight);} }
      // else if (doubletagSBCut_V12(ntuple)) {for (int i = 0; i < doubletagSBPlots.size(); i++) {doubletagSBPlots[i].fill(ntuple,weight);} }
      // else if (tagSRCut_V12(ntuple)) {for (int i = 0; i < tagSRPlots.size(); i++) {tagSRPlots[i].fill(ntuple,weight);} }
      // else if (tagSBCut_V12(ntuple)) {for (int i = 0; i < tagSBPlots.size(); i++) {tagSBPlots[i].fill(ntuple,weight);} }
      // else if (antitagSRCut_V12(ntuple)) {for (int i = 0; i < antitagSRPlots.size(); i++) {antitagSRPlots[i].fill(ntuple,weight);} }
      // else if (antitagSBCut_V12(ntuple)) {for (int i = 0; i < antitagSBPlots.size(); i++) {antitagSBPlots[i].fill(ntuple,weight);} }

    }// end event loop
    outputFile->cd();
  }// end sample loop


  // Begin data
  if (region!=0 && region!=1 && runData) {
    RA2bTree* ntuple = skims.dataNtuple;
    for (int iBin = 0; iBin < numMETbins_V12 ; iBin++) {
      for (int iPlot = 0; iPlot < plots[iBin].size(); iPlot++) {
        plots[iBin][iPlot].addDataNtuple(ntuple,"data");
      }
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
    ntupleBranchStatus<RA2bTree>(ntuple);

    bool passBaseline;
    double jetMass1,jetMass2;
    for (int iEvt = 0; iEvt < numEvents; iEvt++) {
      ntuple->GetEntry(iEvt);
      if (iEvt % 100000 == 0 ) cout << "data: " << iEvt << "/" << numEvents << endl;

      passBaseline=true;
      for (auto baselineCut : baselineCuts ) {
        if (! passBaseline ) continue;
        passBaseline&=baselineCut(ntuple);
      }
      if (! passBaseline ) continue;
      if (ntuple->MET<300.) continue; //MET cut


      // if (region == 0 || region == 1) {
      //   if (!signalTriggerCut(ntuple) ) continue;
      // }
      // else if (region == 2) {
      //   if (!singleMuTriggerCut(ntuple) ) continue;
      // }
      // else if (region == 3) {
      //   if (!singleEleTriggerCut(ntuple) ) continue;
      // }
      // else if (region == 5 ) {
      //   if (!lowDphiTriggerCut(ntuple) ) continue;
      // }

      int bin = -1;
      for (int iBin = 0; iBin < numMETbins_V12 ; iBin++) {
        if (ntuple->MET > lowestMET_V12 ) {
          if (ntuple->MET > lowestMET_V12+binWidth_V12*(numMETbins_V12-1) )
          bin = numMETbins_V12-1;
          else
          bin = int((ntuple->MET-lowestMET_V12)/binWidth_V12);
        }
      }
      if (bin < 0 ) continue;


      // first check for 2H
      if (doubleTaggingLooseCut(ntuple) && region!=1) {
        if (doubletagSRCut(ntuple) && region!=0) {
          // plots[bin][0].fillData(ntuple);
          for (int i = 0; i < doubletagSRPlots.size(); i++) {doubletagSRPlots[i].fillData(ntuple);}
        }
        else if (doubletagSBCut(ntuple)) {
          // plots[bin][1].fillData(ntuple);
          for (int i = 0; i < doubletagSBPlots.size(); i++) {doubletagSBPlots[i].fillData(ntuple);}
        }
      } //end 2H

      //then check for 1H
      if (singleHiggsTagLooseCut(ntuple) && region!=1) {
        if (tagSRCut( ntuple ) && region!=0) {
          // plots[bin][2].fillData(ntuple);
          for (int i = 0; i < tagSRPlots.size(); i++) tagSRPlots[i].fillData(ntuple);
        }
        else if (tagSBCut( ntuple )) {
          // plots[bin][3].fillData(ntuple);
          for (int i = 0; i < tagSBPlots.size(); i++) tagSBPlots[i].fillData(ntuple);
        }
      } //end 1H

      //then check for 0H
      if (antiTaggingLooseCut(ntuple) && region!=1) {
        if (antitagSRCut( ntuple )) {
          // plots[bin][4].fillData(ntuple);
          for (int i = 0; i < antitagSRPlots.size(); i++) antitagSRPlots[i].fillData(ntuple);
          if (antitagSRCut_opt1( ntuple )) {
            // plots[bin][6].fillData(ntuple);
            for (int i = 0; i < antitagSROpt1Plots.size(); i++) antitagSROpt1Plots[i].fillData(ntuple);
          }
          if (antitagSRCut_opt2( ntuple )) {
            // plots[bin][7].fillData(ntuple);
            for (int i = 0; i < antitagSROpt2Plots.size(); i++) antitagSROpt2Plots[i].fillData(ntuple);
          }
        }
        else {
          if (antitagSBCut( ntuple )) {
            // plots[bin][5].fillData(ntuple);
            for (int i = 0; i < antitagSBPlots.size(); i++) antitagSBPlots[i].fillData(ntuple);
            if (antitagSBCut_opt1( ntuple )) {
              // plots[bin][8].fillData(ntuple);
              for (int i = 0; i < antitagSBOpt1Plots.size(); i++) antitagSBOpt1Plots[i].fillData(ntuple);
            }
            if (antitagSBCut_opt2( ntuple )) {
              // plots[bin][9].fillData(ntuple);
              for (int i = 0; i < antitagSBOpt2Plots.size(); i++) antitagSBOpt2Plots[i].fillData(ntuple);
            }
          }
        }
      } //end 0H region

    }// end event loop
    // end Data
  } // end if region!=0 or 1


  bool sumBkgs = true;
  if (region==1) sumBkgs = false;
  if (runData) sumBkgs = false;


  //this just saves a TON of plots per MET bin... maybe these are useful?
  // for (int iBin = 0; iBin < numMETbins_V12; iBin++) {
  //   for (int iPlot = 0; iPlot < plots[iBin].size(); iPlot++) {
  //     outputFile->cd();
  //     plots[iBin][iPlot].buildSum();
  //     plots[iBin][iPlot].Write();
  //     if (sumBkgs) plots[iBin][iPlot].sum->Write();
  //   }
  // }
  // MET_v_photonPT->Write();

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
  for (int i = 0; i < antitagSROpt2Plots.size(); i++) {
    outputFile->cd();
    antitagSROpt2Plots[i].Write();
    if (sumBkgs) {
      antitagSROpt2Plots[i].buildSum("antitagSROpt2");
      antitagSROpt2Plots[i].sum->Write();
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
  for (int i = 0; i < antitagSBOpt2Plots.size(); i++) {
    outputFile->cd();
    antitagSBOpt2Plots[i].Write();
    if (sumBkgs) {
      antitagSBOpt2Plots[i].buildSum("antitagSBOpt2");
      antitagSBOpt2Plots[i].sum->Write();
    }
  }

  for (int i = 0; i < notagSRPlots.size(); i++) {
    outputFile->cd();
    notagSRPlots[i].Write();
    if (sumBkgs) {
      notagSRPlots[i].buildSum("notagSR");
      notagSRPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < notagSBPlots.size(); i++) {
    outputFile->cd();
    notagSBPlots[i].Write();
    if (sumBkgs) {
      notagSBPlots[i].buildSum("notagSB");
      notagSBPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < doubletagSRLeadPlots.size(); i++) {
    outputFile->cd();
    doubletagSRLeadPlots[i].Write();
    if (sumBkgs) {
      doubletagSRLeadPlots[i].buildSum("doubletagSRLead");
      doubletagSRLeadPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < antitagSRLeadPlots.size(); i++) {
    outputFile->cd();
    antitagSRLeadPlots[i].Write();
    if (sumBkgs) {
      antitagSRLeadPlots[i].buildSum("antitagSRLead");
      antitagSRLeadPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < doubletagSBLeadPlots.size(); i++) {
    outputFile->cd();
    doubletagSBLeadPlots[i].Write();
    if (sumBkgs) {
      doubletagSBLeadPlots[i].buildSum("doubletagSBLead");
      doubletagSBLeadPlots[i].sum->Write();
    }
  }
  for (int i = 0; i < antitagSBLeadPlots.size(); i++) {
    outputFile->cd();
    antitagSBLeadPlots[i].Write();
    if (sumBkgs) {
      antitagSBLeadPlots[i].buildSum("antitagSBLead");
      antitagSBLeadPlots[i].sum->Write();
    }
  }

  outputFile->Close();
}
