#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "THStack.h"
#include "TPad.h"

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
  if (region == 0) skims_ = new skimSamples(skimSamples::kSignal, Year);
  else if (region == 1) skims_ = new skimSamples(skimSamples::kSLm, Year);
  else if (region == 2) skims_ = new skimSamples(skimSamples::kSLe, Year);
  else if (region == 3) skims_ = new skimSamples(skimSamples::kLowDphi, Year);
  else assert(1);

  typedef bool(*cuts)(RA2bTree*);
  vector<cuts> baselineCuts;


  if (region == 0) {
    baselineCuts.push_back(*resolvedBaselineCut<RA2bTree>);
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


  skimSamples skims = *skims_;

  typedef plot<RA2bTree> plot;

  double mJbins[4]={50.,85.,135.,250.};
  // double METbins[4]={300.,500.,700.,1000.};
  double METbins[5]={150.,200.,300.,500.,700.};
  vector<vector<plot> > plots;

  for (int i = 0 ; i < numMETbins ; i++) {
    TString tag="_";
    tag+=lowestMET+i*binWidth;
    vector<plot> plotsTemp;
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_tagSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_tagSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_twobSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_twobSB"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_fourbSR"+tag,"m_{J} [GeV]",3,mJbins));
    plotsTemp.push_back(plot(*fillLeadingJetMass<RA2bTree>,"mJ_fourbSB"+tag,"m_{J} [GeV]",3,mJbins));
    plots.push_back(plotsTemp);
  }

  //vector<plot> tempPlots;
  plot MET_Plot(*fillMET<RA2bTree>,"MET","m_{J} [GeV]",4,METbins);
  plot J1pt_Ptplot(*fillLeadingJetPt<RA2bTree>,"J1pt_Pt","p_{T,J} [GeV]",50,300.,1300.);
  plot J2pt_Ptplot(*fillSubLeadingJetPt<RA2bTree>,"J2pt_Pt","p_{T,J} [GeV]",50,300.,1300.);
  plot J1pt_Mplot(*fillLeadingJetMass<RA2bTree>,"J1pt_M","m_{J} [GeV]",50,50.,250.);
  plot J2pt_Mplot(*fillSubLeadingJetMass<RA2bTree>,"J2pt_M","m_{J} [GeV]",50,50.,250.);
  plot ClosestMass(*fillClosestJetMass<RA2bTree>,"ClosestMass","m_{J} [GeV]",3,mJbins);
  plot FarthestMass(*fillFarthestJetMass<RA2bTree>,"FarthestMass","m_{J} [GeV]",3,mJbins);

  vector<plot> fourbSRPlots;
  fourbSRPlots.push_back(plot(MET_Plot));
  fourbSRPlots.push_back(plot(J1pt_Ptplot));
  fourbSRPlots.push_back(plot(J2pt_Ptplot));
  fourbSRPlots.push_back(plot(J1pt_Mplot));
  fourbSRPlots.push_back(plot(J2pt_Mplot));
  fourbSRPlots.push_back(plot(ClosestMass));
  fourbSRPlots.push_back(plot(FarthestMass));

  vector<plot> fourbSBPlots;
  fourbSBPlots.push_back(plot(MET_Plot));
  fourbSBPlots.push_back(plot(J1pt_Ptplot));
  fourbSBPlots.push_back(plot(J2pt_Ptplot));
  fourbSBPlots.push_back(plot(J1pt_Mplot));
  fourbSBPlots.push_back(plot(J2pt_Mplot));
  fourbSBPlots.push_back(plot(ClosestMass));
  fourbSBPlots.push_back(plot(FarthestMass));

  vector<plot> threebSRPlots;
  threebSRPlots.push_back(plot(MET_Plot));
  threebSRPlots.push_back(plot(J1pt_Ptplot));
  threebSRPlots.push_back(plot(J2pt_Ptplot));
  threebSRPlots.push_back(plot(J1pt_Mplot));
  threebSRPlots.push_back(plot(J2pt_Mplot));
  threebSRPlots.push_back(plot(ClosestMass));
  threebSRPlots.push_back(plot(FarthestMass));

  vector<plot> threebSBPlots;
  threebSBPlots.push_back(plot(MET_Plot));
  threebSBPlots.push_back(plot(J1pt_Ptplot));
  threebSBPlots.push_back(plot(J2pt_Ptplot));
  threebSBPlots.push_back(plot(J1pt_Mplot));
  threebSBPlots.push_back(plot(J2pt_Mplot));
  threebSBPlots.push_back(plot(ClosestMass));
  threebSBPlots.push_back(plot(FarthestMass));

  vector<plot> twobSRPlots;
  twobSRPlots.push_back(plot(MET_Plot));
  twobSRPlots.push_back(plot(J1pt_Ptplot));
  twobSRPlots.push_back(plot(J2pt_Ptplot));
  twobSRPlots.push_back(plot(J1pt_Mplot));
  twobSRPlots.push_back(plot(J2pt_Mplot));
  twobSRPlots.push_back(plot(ClosestMass));
  twobSRPlots.push_back(plot(FarthestMass));

  vector<plot> twobSBPlots;
  twobSBPlots.push_back(plot(MET_Plot));
  twobSBPlots.push_back(plot(J1pt_Ptplot));
  twobSBPlots.push_back(plot(J2pt_Ptplot));
  twobSBPlots.push_back(plot(J1pt_Mplot));
  twobSBPlots.push_back(plot(J2pt_Mplot));
  twobSBPlots.push_back(plot(ClosestMass));
  twobSBPlots.push_back(plot(FarthestMass));

  // background MC samples - 0 lepton regions
  for (int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++){
    RA2bTree* ntuple = skims.ntuples[iSample];

    for (int iBin = 0 ; iBin < numMETbins ; iBin++){
      for (int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
        plots[iBin][iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
        plots[iBin][iPlot].setFillColor(ntuple,skims.fillColor[iSample]);
      }
    }
    for (int i = 0 ; i < fourbSRPlots.size() ; i++) {
      fourbSRPlots[i].addNtuple(ntuple,"fourbSR_"+skims.sampleName[iSample]);
      fourbSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < fourbSBPlots.size() ; i++) {
      fourbSBPlots[i].addNtuple(ntuple,"fourbSB_"+skims.sampleName[iSample]);
      fourbSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }

    for (int i = 0 ; i < threebSRPlots.size() ; i++) {
      threebSRPlots[i].addNtuple(ntuple,"threebSR_"+skims.sampleName[iSample]);
      threebSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < threebSBPlots.size() ; i++) {
      threebSBPlots[i].addNtuple(ntuple,"threebSB_"+skims.sampleName[iSample]);
      threebSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < twobSRPlots.size() ; i++) {
      twobSRPlots[i].addNtuple(ntuple,"twobSR_"+skims.sampleName[iSample]);
      twobSRPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
    }
    for (int i = 0 ; i < twobSBPlots.size() ; i++) {
      twobSBPlots[i].addNtuple(ntuple,"twobSB_"+skims.sampleName[iSample]);
      twobSBPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
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
    if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH")) {
      TFile *fin = new TFile(filename,"READ");
      TH1F *nEventsHisto = (TH1F*)fin->Get("nEventProc");
      TotalEvents = nEventsHisto->GetBinContent(1);
    }

    for (int iEvt = 0 ; iEvt < numEvents ; iEvt++) {
      ntuple->GetEntry(iEvt);
      if (iEvt % 10000 == 0) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << min(MAX_EVENTS,numEvents) << endl;

      if (region==0) {
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
        if(skims.sampleName[iSample]=="QCD") trigWeight=trigcorrorFakeMHT.GetCorrection(ntuple->MHT,trigunc);
      }
      else if (region == 1) {
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      }
      else if (region == 2) {
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      }
      else if (region == 3) {
        if(skims.sampleName[iSample]!="QCD") {
          trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
          trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
        }
        else{
          trigWeight=1.0;
        }
      }

      passBaseline=true;
      for (auto baselineCut : baselineCuts) {
        if (! passBaseline) continue;
        passBaseline&=baselineCut(ntuple);
      }
      if (! passBaseline) continue;

      if ((filename.Contains("SingleLept") || filename.Contains("DiLept")) && ntuple->madHT>600.) continue;
      bin = -1;

      weight = ntuple->Weight*this_lumi*trigWeight;
      if (filename.Contains("TChiHH_HToBB")) weight = weight*0.5823329*0.5823329/TotalEvents;
      else if (filename.Contains("T5qqqqZH")) weight = weight/TotalEvents;

      //if (skims.sampleName[iSample] == "TT") {
      //    weight *= ISRweights(ntuple);
      //}

      for (int iBin = 0 ; iBin < numMETbins ; iBin++) {
        if (ntuple->MET > lowestMET) {
          if (ntuple->MET > numMETbins*(binWidth-1)+lowestMET) bin = numMETbins-1;
          else bin = int((ntuple->MET-lowestMET)/binWidth);
        }
      }

      if (bin < 0) continue;

      if (fourbSRCut(ntuple)) {
        plots[bin][4].fill(ntuple,weight);
        for (int i = 0 ; i < fourbSRPlots.size() ; i++) {
          fourbSRPlots[i].fill (ntuple,weight);
        }
      }
      else if (fourbSBCut(ntuple)) {
        plots[bin][5].fill(ntuple,weight);
        for (int i = 0 ; i < fourbSBPlots.size() ; i++)
        fourbSBPlots[i].fill (ntuple,weight);
      }

      else if (threebSRCut(ntuple)) {
        plots[bin][0].fill(ntuple,weight);
        for (int i = 0 ; i < threebSRPlots.size() ; i++)
        threebSRPlots[i].fill (ntuple,weight);
      }
      else if (threebSBCut(ntuple)) {
        plots[bin][1].fill(ntuple,weight);
        for (int i = 0 ; i < threebSBPlots.size() ; i++)
        threebSBPlots[i].fill (ntuple,weight);
      }
      else if (twobSRCut(ntuple)) {
        plots[bin][2].fill(ntuple,weight);
        for (int i = 0 ; i < twobSRPlots.size() ; i++)
        twobSRPlots[i].fill (ntuple,weight);
      }
      else if (twobSBCut(ntuple)) {
        plots[bin][3].fill(ntuple,weight);
        for (int i = 0 ; i < twobSBPlots.size() ; i++)
        twobSBPlots[i].fill (ntuple,weight);
      } // end if-else-if block for tagging regions

    }// end event loop
  }// end sample loop


  // data
  // RA2bTree* ntuple = skims.dataNtuple;
  //
  // for (int iBin = 0 ; iBin < numMETbins ; iBin++){
  //   for (int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
  //     plots[iBin][iPlot].addDataNtuple(ntuple,"data");
  //   }
  // }
  // for (int i = 0 ; i < fourbSRPlots.size() ; i++) {
  //   fourbSRPlots[i].addDataNtuple(ntuple,"fourbSR_data");
  // }
  // for (int i = 0 ; i < fourbSBPlots.size() ; i++) {
  //   fourbSBPlots[i].addDataNtuple(ntuple,"fourbSB_data");
  // }
  // for (int i = 0 ; i < threebSRPlots.size() ; i++) {
  //   threebSRPlots[i].addDataNtuple(ntuple,"threebSR_data");
  // }
  // for (int i = 0 ; i < threebSBPlots.size() ; i++) {
  //   threebSBPlots[i].addDataNtuple(ntuple,"threebSB_data");
  // }
  // for (int i = 0 ; i < twobSRPlots.size() ; i++) {
  //   twobSRPlots[i].addDataNtuple(ntuple,"twobSR_data");
  // }
  // for (int i = 0 ; i < twobSBPlots.size() ; i++) {
  //   twobSBPlots[i].addDataNtuple(ntuple,"twobSB_data");
  // }
  //
  // int numEvents = ntuple->fChain->GetEntries();
  // ntupleBranchStatus<RA2bTree>(ntuple);
  // bool passBaseline;
  // double jetMass1,jetMass2;
  // for (int iEvt = 0 ; iEvt < min(MAX_EVENTS,numEvents) ; iEvt++) {
  //   ntuple->GetEntry(iEvt);
  //   if (iEvt % 100000 == 0) cout << "data: " << iEvt << "/" << min(MAX_EVENTS,numEvents) << endl;
  //
  //   passBaseline=true;
  //   for (auto baselineCut : baselineCuts) {
  //     passBaseline&=baselineCut(ntuple);
  //   }
  //   if (! passBaseline) continue;
  //
  //   if (region == 0) {
  //     if (!signalTriggerCut(ntuple)) continue;
  //   }else if (region == 1){
  //     if (!singleMuTriggerCut(ntuple)) continue;
  //   }else if (region == 2){
  //     if (!singleEleTriggerCut(ntuple)) continue;
  //   }else if (region == 3) {
  //     if (!lowDphiTriggerCut(ntuple)) continue;
  //   }
  //
  //   int bin = -1;
  //   for (int iBin = 0 ; iBin < numMETbins ; iBin++) {
  //     if (ntuple->MET > lowestMET) {
  //       if (ntuple->MET > lowestMET+binWidth*(numMETbins-1))
  //       bin = numMETbins-1;
  //       else
  //       bin = int((ntuple->MET-lowestMET)/binWidth);
  //     }
  //   }
  //   if (bin < 0) continue;
  //
  //   if (fourbSRCut(ntuple)) {
  //     if (region != 0) {
  //       plots[bin][4].fillData(ntuple);
  //       for (int i = 0 ; i < fourbSRPlots.size() ; i++)
  //       fourbSRPlots[i].fillData(ntuple);
  //     }
  //   }else if (fourbSBCut(ntuple)) {
  //     plots[bin][5].fillData(ntuple);
  //     for (int i = 0 ; i < fourbSBPlots.size() ; i++)
  //     fourbSBPlots[i].fillData(ntuple);
  //   }else if (threebSRCut(ntuple)) {
  //     if (region != 0) {
  //       plots[bin][0].fillData(ntuple);
  //       for (int i = 0 ; i < threebSRPlots.size() ; i++)
  //       threebSRPlots[i].fillData(ntuple);
  //     }
  //   }else if (threebSBCut(ntuple)) {
  //     plots[bin][1].fillData(ntuple);
  //     for (int i = 0 ; i < threebSBPlots.size() ; i++)
  //     threebSBPlots[i].fillData(ntuple);
  //   }else if (twobSRCut(ntuple)) {
  //     plots[bin][2].fillData(ntuple);
  //     for (int i = 0 ; i < twobSRPlots.size() ; i++)
  //     twobSRPlots[i].fillData(ntuple);
  //   }else if (twobSBCut(ntuple)) {
  //     plots[bin][3].fillData(ntuple);
  //     for (int i = 0 ; i < twobSBPlots.size() ; i++)
  //     twobSBPlots[i].fillData(ntuple);
  //   }// end if-else-if block for tagging regions
  // }// end event loop


  TFile* outputFile;
  TString regionName;
  TString cutName="";
  if (looseCuts) cutName="_looseCuts";
  if (region == 0) regionName="";
  if (region == 1) regionName="_singleMu";
  if (region == 2) regionName="_singleEle";
  if (region == 3) regionName="_lowDphi";

  outputFile = new TFile("ALPHABETRes"+Year+"_V18"+regionName+".root","RECREATE");

  bool sumBkgs = false;

  for (int iBin = 0 ; iBin < numMETbins; iBin++){
    for (int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
      outputFile->cd();
      plots[iBin][iPlot].buildSum();
      plots[iBin][iPlot].Write();
      if (sumBkgs) plots[iBin][iPlot].sum->Write();
    }
  }

  for (int i = 0 ; i < fourbSRPlots.size() ; i++) {
    outputFile->cd();
    fourbSRPlots[i].buildSum("fourbSR");
    fourbSRPlots[i].Write();
    if (sumBkgs) fourbSRPlots[i].sum->Write();
  }
  for (int i = 0 ; i < fourbSBPlots.size() ; i++) {
    outputFile->cd();
    fourbSBPlots[i].buildSum("fourbSB");
    fourbSBPlots[i].Write();
    if (sumBkgs) fourbSBPlots[i].sum->Write();
  }
  for (int i = 0 ; i < threebSRPlots.size() ; i++) {
    outputFile->cd();
    threebSRPlots[i].buildSum("threebSR");
    threebSRPlots[i].Write();
    if (sumBkgs) threebSRPlots[i].sum->Write();
  }
  for (int i = 0 ; i < threebSBPlots.size() ; i++) {
    outputFile->cd();
    threebSBPlots[i].buildSum("threebSB");
    threebSBPlots[i].Write();
    if (sumBkgs) threebSBPlots[i].sum->Write();
  }
  for (int i = 0 ; i < twobSRPlots.size() ; i++) {
    outputFile->cd();
    twobSRPlots[i].buildSum("twobSR");
    twobSRPlots[i].Write();
    if (sumBkgs) twobSRPlots[i].sum->Write();
  }
  for (int i = 0 ; i < twobSBPlots.size() ; i++) {
    outputFile->cd();
    twobSBPlots[i].buildSum("twobSB");
    twobSBPlots[i].Write();
    if (sumBkgs) twobSBPlots[i].sum->Write();
  }

  outputFile->Close();
}
