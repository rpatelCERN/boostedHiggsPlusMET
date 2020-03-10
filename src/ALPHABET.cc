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

int main(int argc, char** argv){

  int region(0);
  bool looseCuts(false);
  int MAX_EVENTS(99999999);

  region = atoi(argv[1]);

  gROOT->ProcessLine(".L tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  TString Year(argv[2]);
  TriggerCorrector trigcorror;
  TriggerCorrector trigcorrorHT;
  TriggerCorrector trigcorrorFakeMHT;
  trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2016");
  trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2016");
  trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2016");


  skimSamples* skims_;
  if (region == 0 ) skims_ = new skimSamples(skimSamples::kSignal, Year);
  else if (region == 1 ) skims_ = new skimSamples(skimSamples::kSLm, Year);
  else if (region == 2 ) skims_ = new skimSamples(skimSamples::kSLe, Year);
  else if (region == 3 ) skims_ = new skimSamples(skimSamples::kLowDphi, Year);
  else assert(1);

  typedef bool(*cuts)(RA2bTree*);
  vector<cuts> baselineCuts;

  if (looseCuts ){
    baselineCuts.push_back(*FiltersCut<RA2bTree>);
    if (region == 3 ){
      baselineCuts.push_back(*lowDPhiCuts<RA2bTree>);
    }
    else {
      baselineCuts.push_back(*DeltaPhiCuts<RA2bTree>);
    }
    if (region == 1 ){
      baselineCuts.push_back(*singleMuCut<RA2bTree>);
    }
    if (region == 2 ){
      baselineCuts.push_back(*singleEleCut<RA2bTree>);
    }
    baselineCuts.push_back(*METHTlooseCut<RA2bTree>);
    baselineCuts.push_back(*AK8MultCut<RA2bTree>);
  }
  else {
    if (region == 0 ){
      baselineCuts.push_back(*boostedBaselineCut<RA2bTree>);
    }
    else if (region == 1){
      baselineCuts.push_back(*singleMuBaselineCut<RA2bTree>);
    }
    else if (region == 2){
      baselineCuts.push_back(*singleEleBaselineCut<RA2bTree>);
    }
    else if (region == 3){
      baselineCuts.push_back(*lowDphiBaselineCut<RA2bTree>);
    }
    else assert(1);
  }

  skimSamples skims = *skims_;

  typedef plot<RA2bTree> plot;

  double mJbins[4]={50.,85.,135.,250.};
  double METbins[5]={250.,300.,500.,700.,1400.};
  // double METbins[4]={300.,500.,700.,1400.};
  vector<vector<plot> > plots;

  TString METbins_string[5]={"250","300","500","700","1400"};

  for (int i = 0 ; i < numMETbins ; i++){
    TString tag="_";
    tag+=METbins_string[i];
    // tag+=lowestMET+i*binWidth;
    vector<plot> plotsTemp;
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_tagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_tagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_antitagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_doubletagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_doubletagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_notagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_notagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_doubletagSB2"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_anttagSRWeight"+tag,"m_{J} [GeV]",3,mJbins));
    plots.push_back(plotsTemp);
  }

  //vector<plot> tempPlots;
  plot MET_Plot(*fillMET<RA2bTree>,"MET","m_{J} [GeV]",4,METbins);
  // plot MET_Plot(*fillMET<RA2bTree>,"MET","m_{J} [GeV]",3,METbins);
  plot J1pt_Ptplot(*fillLeadingJetPt<RA2bTree>,"J1pt_Pt","p_{T,J} [GeV]",50,300.,1300.);
  plot J2pt_Ptplot(*fillSubLeadingJetPt<RA2bTree>,"J2pt_Pt","p_{T,J} [GeV]",50,300.,1300.);
  plot J1pt_Mplot(*fillLeadingJetMass<RA2bTree>,"J1pt_M","m_{J} [GeV]",80,50.,250.);
  plot J2pt_Mplot(*fillSubLeadingJetMass<RA2bTree>,"J2pt_M","m_{J} [GeV]",80,50.,250.);
  plot ClosestMass(*fillClosestJetMass<RA2bTree>,"ClosestMass","m_{J} [GeV]",3,mJbins);
  plot FarthestMass(*fillFarthestJetMass<RA2bTree>,"FarthestMass","m_{J} [GeV]",3,mJbins);
  plot METRatio(*fillMETRatio<RA2bTree>,"METRatio","MET/MHT",50,0,700);

  vector<plot> doubletagSRPlots;
  doubletagSRPlots.push_back(plot(MET_Plot));
  doubletagSRPlots.push_back(plot(J1pt_Ptplot));
  doubletagSRPlots.push_back(plot(J2pt_Ptplot));
  doubletagSRPlots.push_back(plot(J1pt_Mplot));
  doubletagSRPlots.push_back(plot(J2pt_Mplot));
  doubletagSRPlots.push_back(plot(ClosestMass));
  doubletagSRPlots.push_back(plot(FarthestMass));
  doubletagSRPlots.push_back(plot(METRatio));

  vector<plot> doubletagSBPlots;
  doubletagSBPlots.push_back(plot(MET_Plot));
  doubletagSBPlots.push_back(plot(J1pt_Ptplot));
  doubletagSBPlots.push_back(plot(J2pt_Ptplot));
  doubletagSBPlots.push_back(plot(J1pt_Mplot));
  doubletagSBPlots.push_back(plot(J2pt_Mplot));
  doubletagSBPlots.push_back(plot(ClosestMass));
  doubletagSBPlots.push_back(plot(FarthestMass));
  doubletagSBPlots.push_back(plot(METRatio));

  vector<plot> doubletagSB2Plots;
  doubletagSB2Plots.push_back(plot(MET_Plot));
  doubletagSB2Plots.push_back(plot(J1pt_Ptplot));
  doubletagSB2Plots.push_back(plot(J2pt_Ptplot));
  doubletagSB2Plots.push_back(plot(J1pt_Mplot));
  doubletagSB2Plots.push_back(plot(J2pt_Mplot));
  doubletagSB2Plots.push_back(plot(ClosestMass));
  doubletagSB2Plots.push_back(plot(FarthestMass));
  doubletagSB2Plots.push_back(plot(METRatio));

  vector<plot> tagSRPlots;
  tagSRPlots.push_back(plot(MET_Plot));
  tagSRPlots.push_back(plot(J1pt_Ptplot));
  tagSRPlots.push_back(plot(J2pt_Ptplot));
  tagSRPlots.push_back(plot(J1pt_Mplot));
  tagSRPlots.push_back(plot(J2pt_Mplot));
  tagSRPlots.push_back(plot(ClosestMass));
  tagSRPlots.push_back(plot(FarthestMass));
  tagSRPlots.push_back(plot(METRatio));

  vector<plot> tagSBPlots;
  tagSBPlots.push_back(plot(MET_Plot));
  tagSBPlots.push_back(plot(J1pt_Ptplot));
  tagSBPlots.push_back(plot(J2pt_Ptplot));
  tagSBPlots.push_back(plot(J1pt_Mplot));
  tagSBPlots.push_back(plot(J2pt_Mplot));
  tagSBPlots.push_back(plot(ClosestMass));
  tagSBPlots.push_back(plot(FarthestMass));
  tagSBPlots.push_back(plot(METRatio));

  vector<plot> antitagSRPlots;
  antitagSRPlots.push_back(plot(MET_Plot));
  antitagSRPlots.push_back(plot(J1pt_Ptplot));
  antitagSRPlots.push_back(plot(J2pt_Ptplot));
  antitagSRPlots.push_back(plot(J1pt_Mplot));
  antitagSRPlots.push_back(plot(J2pt_Mplot));
  antitagSRPlots.push_back(plot(ClosestMass));
  antitagSRPlots.push_back(plot(FarthestMass));
  antitagSRPlots.push_back(plot(METRatio));

  vector<plot> antitagSRWeightPlots;
  antitagSRWeightPlots.push_back(plot(MET_Plot));
  antitagSRWeightPlots.push_back(plot(J1pt_Ptplot));
  antitagSRWeightPlots.push_back(plot(J2pt_Ptplot));
  antitagSRWeightPlots.push_back(plot(J1pt_Mplot));
  antitagSRWeightPlots.push_back(plot(J2pt_Mplot));
  antitagSRWeightPlots.push_back(plot(ClosestMass));
  antitagSRWeightPlots.push_back(plot(FarthestMass));
  antitagSRWeightPlots.push_back(plot(METRatio));

  vector<plot> antitagSBPlots;
  antitagSBPlots.push_back(plot(MET_Plot));
  antitagSBPlots.push_back(plot(J1pt_Ptplot));
  antitagSBPlots.push_back(plot(J2pt_Ptplot));
  antitagSBPlots.push_back(plot(J1pt_Mplot));
  antitagSBPlots.push_back(plot(J2pt_Mplot));
  antitagSBPlots.push_back(plot(ClosestMass));
  antitagSBPlots.push_back(plot(FarthestMass));
  antitagSBPlots.push_back(plot(METRatio));


  vector<plot> notagSRPlots;
  notagSRPlots.push_back(plot(MET_Plot));
  notagSRPlots.push_back(plot(J1pt_Ptplot));
  notagSRPlots.push_back(plot(J2pt_Ptplot));
  notagSRPlots.push_back(plot(J1pt_Mplot));
  notagSRPlots.push_back(plot(J2pt_Mplot));
  notagSRPlots.push_back(plot(ClosestMass));
  notagSRPlots.push_back(plot(FarthestMass));
  notagSRPlots.push_back(plot(METRatio));

  vector<plot> notagSBPlots;
  notagSBPlots.push_back(plot(MET_Plot));
  notagSBPlots.push_back(plot(J1pt_Ptplot));
  notagSBPlots.push_back(plot(J2pt_Ptplot));
  notagSBPlots.push_back(plot(J1pt_Mplot));
  notagSBPlots.push_back(plot(J2pt_Mplot));
  notagSBPlots.push_back(plot(ClosestMass));
  notagSBPlots.push_back(plot(FarthestMass));
  notagSBPlots.push_back(plot(METRatio));


  TFile* outputFile;
  TString regionName;
  TString cutName="";
  if (looseCuts )
  cutName="_looseCuts";
  if (region == 0 )
  regionName="";
  if (region == 1 )
  regionName="_singleMu";
  if (region == 2 )
  regionName="_singleEle";
  if (region == 3 )
  regionName="_lowDphi";

  // regionName = "_resVeto";
  outputFile = new TFile("ALPHABETBoost"+Year+"_V17"+regionName+".root","RECREATE");
  // outputFile = new TFile("ALPHABETBoost_V17_signalOnly.root","RECREATE");
  // outputFile = new TFile("ALPHABETBoost_MC2016_V17_oldBBTag.root","RECREATE");

  // background MC samples - 0 lepton regions
  for (int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++){
    RA2bTree* ntuple = skims.ntuples[iSample];

    for (int iBin = 0 ; iBin < numMETbins ; iBin++){
      for (int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
        plots[iBin][iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
        plots[iBin][iPlot].setFillColor(ntuple,skims.fillColor[iSample]);
      }
    }
    for (int i = 0 ; i < doubletagSRPlots.size() ; i++){
      doubletagSRPlots[i].addNtuple(ntuple,"doubletagSR_"+skims.sampleName[iSample]);
      doubletagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < doubletagSBPlots.size() ; i++){
      doubletagSBPlots[i].addNtuple(ntuple,"doubletagSB_"+skims.sampleName[iSample]);
      doubletagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < doubletagSB2Plots.size() ; i++){
      doubletagSB2Plots[i].addNtuple(ntuple,"doubletagSB2_"+skims.sampleName[iSample]);
      doubletagSB2Plots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < tagSRPlots.size() ; i++){
      tagSRPlots[i].addNtuple(ntuple,"tagSR_"+skims.sampleName[iSample]);
      tagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < tagSBPlots.size() ; i++){
      tagSBPlots[i].addNtuple(ntuple,"tagSB_"+skims.sampleName[iSample]);
      tagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < antitagSRPlots.size() ; i++){
      antitagSRPlots[i].addNtuple(ntuple,"antitagSR_"+skims.sampleName[iSample]);
      antitagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < antitagSRWeightPlots.size() ; i++){
      antitagSRWeightPlots[i].addNtuple(ntuple,"antitagSRWeight_"+skims.sampleName[iSample]);
      antitagSRWeightPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < antitagSBPlots.size() ; i++){
      antitagSBPlots[i].addNtuple(ntuple,"antitagSB_"+skims.sampleName[iSample]);
      antitagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < notagSRPlots.size() ; i++){
      notagSRPlots[i].addNtuple(ntuple,"notagSR_"+skims.sampleName[iSample]);
      notagSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < notagSBPlots.size() ; i++){
      notagSBPlots[i].addNtuple(ntuple,"notagSB_"+skims.sampleName[iSample]);
      notagSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }

    int numEvents = ntuple->fChain->GetEntries();
    ntupleBranchStatus<RA2bTree>(ntuple);
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
    if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH") ){
      TFile *fin = new TFile(filename,"READ");
      TH1F *nEventsHisto = (TH1F*)fin->Get("nEventProc");
      TotalEvents = nEventsHisto->GetBinContent(1);
    }

    for (int iEvt = 0 ; iEvt < numEvents ; iEvt++){
      ntuple->GetEntry(iEvt);
      if (iEvt % 10000 == 0 ) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << min(MAX_EVENTS,numEvents) << endl;

      if (region==0){
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
        if(skims.sampleName[iSample]=="QCD") trigWeight=trigcorrorFakeMHT.GetCorrection(ntuple->MHT,trigunc);
      }
      else if (region == 1 ){
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      }
      else if (region == 2 ){
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      }
      else if (region == 3 ){
        if(skims.sampleName[iSample]!="QCD"){
          trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
          trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
        }
        else{
          trigWeight=1.0;
        }
      }

      passBaseline=true;
      for (auto baselineCut : baselineCuts ){
        if (! passBaseline ) continue;
        passBaseline&=baselineCut(ntuple);
      }
      if (! passBaseline ) continue;


      if (( filename.Contains("SingleLept") || filename.Contains("DiLept") ) && ntuple->madHT>600. )continue;

      // weight = ntuple->Weight*this_lumi;
      weight = ntuple->Weight*this_lumi*trigWeight;
      if ( filename.Contains("TChiHH_HToBB") ) weight = weight*0.5823329*0.5823329/TotalEvents;
      else if ( filename.Contains("T5qqqqZH") )  weight = weight/TotalEvents;

      //if (skims.sampleName[iSample] == "TT" ){
      //    weight *= ISRweights(ntuple);
      //}

      //Toggle whether or not we veto resolved events
      // if ( resolvedBaselineCut(ntuple) ) continue;

      bin = -1;
      //let's see...
      //our newest MET bins would be 250-300, 300-500, 500-700, and 700+
      float thisMET_forBin = ntuple->MET;
      if (thisMET_forBin>250.&&thisMET_forBin<=300.) bin=0;
      else if (thisMET_forBin>300.&&thisMET_forBin<=500.) bin=1;
      else if (thisMET_forBin>500.&&thisMET_forBin<=700.) bin=2;
      else if (thisMET_forBin>700.) bin=3;


      // for (int iBin = 0 ; iBin < numMETbins ; iBin++){
      //   if (ntuple->MET > lowestMET ){
      //     if (ntuple->MET > numMETbins*(binWidth-1)+lowestMET )
      //     bin = numMETbins-1;
      //     else
      //     bin = int((ntuple->MET-lowestMET)/binWidth);
      //   }
      // }

      if (bin < 0 ) continue;

      //insert notag plots here
      if (doubleMassCut(ntuple) ){ //this requires both lead and sublead jet to be in the mass window
        plots[bin][6].fill(ntuple,weight);
        for (int i = 0 ; i < notagSRPlots.size() ; i++){
          notagSRPlots[i].fill (ntuple,weight);
        }
      }

      //presumably, this will only fill if the other doesn't,
      //the only cuts on here are baseline cuts, meaning, 2 AK8 jets with mass 50-250, and pT>300
      //So at least one of the jets is not in the mass window
      else if (!doubleMassCut( ntuple ) ){
        plots[bin][7].fill(ntuple,weight);
        for (int i = 0 ; i < notagSBPlots.size() ; i++)
        notagSBPlots[i].fill (ntuple,weight);
      }

      //For fit of J2 mass shape
      if (doubletagSB2Cut( ntuple ) ){
        plots[bin][8].fill(ntuple,weight);
        for (int i = 0 ; i < doubletagSB2Plots.size() ; i++)
        doubletagSB2Plots[i].fill (ntuple,weight);
      }

      if (doubletagSRCut(ntuple) ){
        plots[bin][4].fill(ntuple,weight);
        for (int i = 0 ; i < doubletagSRPlots.size() ; i++){
          doubletagSRPlots[i].fill (ntuple,weight);
        }
      }
      else if (doubletagSBCut( ntuple ) ){
        plots[bin][5].fill(ntuple,weight);
        for (int i = 0 ; i < doubletagSBPlots.size() ; i++)
        doubletagSBPlots[i].fill (ntuple,weight);
      }
      else if (tagSRCut( ntuple ) ){
        plots[bin][0].fill(ntuple,weight);
        for (int i = 0 ; i < tagSRPlots.size() ; i++)
        tagSRPlots[i].fill (ntuple,weight);
      }
      else if (tagSBCut( ntuple ) ){
        plots[bin][1].fill(ntuple,weight);
        for (int i = 0 ; i < tagSBPlots.size() ; i++)
        tagSBPlots[i].fill (ntuple,weight);
      }
      else if (antitagSRCut( ntuple ) ){
        float LeadJetWeightRPF = antitagSRWeight(ntuple);
        plots[bin][2].fill(ntuple,weight);
        for (int i = 0 ; i < antitagSRPlots.size() ; i++)
        antitagSRPlots[i].fill (ntuple,weight);

        //Also fill with weight from RPF fit for lead jet mass
        plots[bin][9].fill(ntuple,weight*LeadJetWeightRPF);
        for (int i = 0 ; i < antitagSRWeightPlots.size() ; i++)
        antitagSRWeightPlots[i].fill (ntuple,weight*LeadJetWeightRPF);
      }
      else if (antitagSBCut( ntuple ) ){
        plots[bin][3].fill(ntuple,weight);
        for (int i = 0 ; i < antitagSBPlots.size() ; i++)
        antitagSBPlots[i].fill (ntuple,weight);
      } // end if-else-if block for tagging regions
    }// end event loop
    outputFile->cd();
  }// end sample loop


  // // Begin data
  // RA2bTree* ntuple = skims.dataNtuple;
  //
  // for (int iBin = 0 ; iBin < numMETbins ; iBin++){
  //   for (int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
  //     plots[iBin][iPlot].addDataNtuple(ntuple,"data");
  //   }
  // }
  // for (int i = 0 ; i < doubletagSRPlots.size() ; i++){
  //   doubletagSRPlots[i].addDataNtuple(ntuple,"doubletagSR_data");
  // }
  // for (int i = 0 ; i < doubletagSBPlots.size() ; i++){
  //   doubletagSBPlots[i].addDataNtuple(ntuple,"doubletagSB_data");
  // }
  // for (int i = 0 ; i < tagSRPlots.size() ; i++){
  //   tagSRPlots[i].addDataNtuple(ntuple,"tagSR_data");
  // }
  // for (int i = 0 ; i < tagSBPlots.size() ; i++){
  //   tagSBPlots[i].addDataNtuple(ntuple,"tagSB_data");
  // }
  // for (int i = 0 ; i < antitagSRPlots.size() ; i++){
  //   antitagSRPlots[i].addDataNtuple(ntuple,"antitagSR_data");
  // }
  // for (int i = 0 ; i < antitagSBPlots.size() ; i++){
  //   antitagSBPlots[i].addDataNtuple(ntuple,"antitagSB_data");
  // }
  //
  // int numEvents = ntuple->fChain->GetEntries();
  // ntupleBranchStatus<RA2bTree>(ntuple);
  // bool passBaseline;
  // double jetMass1,jetMass2;
  // for (int iEvt = 0 ; iEvt < min(MAX_EVENTS,numEvents) ; iEvt++){
  //   ntuple->GetEntry(iEvt);
  //   if (iEvt % 100000 == 0 ) cout << "data: " << iEvt << "/" << min(MAX_EVENTS,numEvents) << endl;
  //
  //   passBaseline=true;
  //   for (auto baselineCut : baselineCuts ){
  //     passBaseline&=baselineCut(ntuple);
  //   }
  //   if (! passBaseline ) continue;
  //
  //   if (region == 0 ){
  //     if (!signalTriggerCut(ntuple) ) continue;
  //   }else if (region == 1){
  //     if (!singleMuTriggerCut(ntuple) ) continue;
  //   }else if (region == 2){
  //     if (!singleEleTriggerCut(ntuple) ) continue;
  //   }else if (region == 3 ){
  //     if (!lowDphiTriggerCut(ntuple) ) continue;
  //   }
  //
  //   int bin = -1;
  //   for (int iBin = 0 ; iBin < numMETbins ; iBin++){
  //     if (ntuple->MET > lowestMET ){
  //       if (ntuple->MET > lowestMET+binWidth*(numMETbins-1) )
  //       bin = numMETbins-1;
  //       else
  //       bin = int((ntuple->MET-lowestMET)/binWidth);
  //     }
  //   }
  //   if (bin < 0 ) continue;
  //
  //   if (doubletagSRCut( ntuple ) ){
  //     if (region != 0 ){
  //       plots[bin][4].fillData(ntuple);
  //       for (int i = 0 ; i < doubletagSRPlots.size() ; i++)
  //       doubletagSRPlots[i].fillData(ntuple);
  //     }
  //   }else if (doubletagSBCut( ntuple ) ){
  //     plots[bin][5].fillData(ntuple);
  //     for (int i = 0 ; i < doubletagSBPlots.size() ; i++)
  //     doubletagSBPlots[i].fillData(ntuple);
  //   }else if (tagSRCut( ntuple ) ){
  //     if (region != 0 ){
  //       plots[bin][0].fillData(ntuple);
  //       for (int i = 0 ; i < tagSRPlots.size() ; i++)
  //       tagSRPlots[i].fillData(ntuple);
  //     }
  //   }else if (tagSBCut( ntuple ) ){
  //     plots[bin][1].fillData(ntuple);
  //     for (int i = 0 ; i < tagSBPlots.size() ; i++)
  //     tagSBPlots[i].fillData(ntuple);
  //   }else if (antitagSRCut( ntuple ) ){
  //     plots[bin][2].fillData(ntuple);
  //     for (int i = 0 ; i < antitagSRPlots.size() ; i++)
  //     antitagSRPlots[i].fillData(ntuple);
  //   }else if (antitagSBCut( ntuple ) ){
  //     plots[bin][3].fillData(ntuple);
  //     for (int i = 0 ; i < antitagSBPlots.size() ; i++)
  //     antitagSBPlots[i].fillData(ntuple);
  //   }// end if-else-if block for tagging regions
  // }// end event loop
  // // end Data


  bool sumBkgs = true;

  for (int iBin = 0 ; iBin < numMETbins; iBin++){
    for (int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
      outputFile->cd();
      plots[iBin][iPlot].buildSum();
      plots[iBin][iPlot].Write();
      if (sumBkgs) plots[iBin][iPlot].sum->Write();
    }
  }

  for (int i = 0 ; i < doubletagSRPlots.size() ; i++){
    outputFile->cd();
    doubletagSRPlots[i].buildSum("doubletagSR");
    doubletagSRPlots[i].Write();
    if (sumBkgs) doubletagSRPlots[i].sum->Write();
  }

  for (int i = 0 ; i < doubletagSBPlots.size() ; i++){
    outputFile->cd();
    doubletagSBPlots[i].buildSum("doubletagSB");
    doubletagSBPlots[i].Write();
    if (sumBkgs) doubletagSBPlots[i].sum->Write();
  }
  for (int i = 0 ; i < doubletagSB2Plots.size() ; i++){
    outputFile->cd();
    doubletagSB2Plots[i].buildSum("doubletagSB2");
    doubletagSB2Plots[i].Write();
    if (sumBkgs) doubletagSB2Plots[i].sum->Write();
  }
  for (int i = 0 ; i < tagSRPlots.size() ; i++){
    outputFile->cd();
    tagSRPlots[i].buildSum("tagSR");
    tagSRPlots[i].Write();
    if (sumBkgs) tagSRPlots[i].sum->Write();
  }
  for (int i = 0 ; i < tagSBPlots.size() ; i++){
    outputFile->cd();
    tagSBPlots[i].buildSum("tagSB");
    tagSBPlots[i].Write();
    if (sumBkgs) tagSBPlots[i].sum->Write();
  }
  for (int i = 0 ; i < antitagSRPlots.size() ; i++){
    outputFile->cd();
    antitagSRPlots[i].buildSum("antitagSR");
    antitagSRPlots[i].Write();
    if (sumBkgs) antitagSRPlots[i].sum->Write();
  }
  for (int i = 0 ; i < antitagSRWeightPlots.size() ; i++){
    outputFile->cd();
    antitagSRWeightPlots[i].buildSum("antitagSRWeight");
    antitagSRWeightPlots[i].Write();
    if (sumBkgs) antitagSRWeightPlots[i].sum->Write();
  }
  for (int i = 0 ; i < antitagSBPlots.size() ; i++){
    outputFile->cd();
    antitagSBPlots[i].buildSum("antitagSB");
    antitagSBPlots[i].Write();
    if (sumBkgs) antitagSBPlots[i].sum->Write();
  }
  for (int i = 0 ; i < notagSRPlots.size() ; i++){
    outputFile->cd();
    notagSRPlots[i].buildSum("notagSR");
    notagSRPlots[i].Write();
    if (sumBkgs) notagSRPlots[i].sum->Write();
  }
  for (int i = 0 ; i < notagSBPlots.size() ; i++){
    outputFile->cd();
    notagSBPlots[i].buildSum("notagSB");
    notagSBPlots[i].Write();
    if (sumBkgs) notagSBPlots[i].sum->Write();
  }
  outputFile->Close();
}
