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
  if (region == 0) skims_ = new skimSamples(skimSamples::kSignal, Year);
  else if (region == 1) skims_ = new skimSamples(skimSamples::kSignalOnly, Year);
  else if (region == 2) skims_ = new skimSamples(skimSamples::kSLm, Year);
  else if (region == 3) skims_ = new skimSamples(skimSamples::kSLe, Year);
  else if (region == 4) skims_ = new skimSamples(skimSamples::kLowDphi, Year);
  else assert(1);

  typedef bool(*cuts)(RA2bTree*);
  vector<cuts> baselineCuts;


  if (region == 0 || region == 1 ) {
    baselineCuts.push_back(*resolvedBaselineCut<RA2bTree>);
  }
  else if (region == 2){
    baselineCuts.push_back(*singleMuBaselineCut<RA2bTree>);
  }
  else if (region == 3){
    baselineCuts.push_back(*singleEleBaselineCut<RA2bTree>);
  }
  else if (region == 4){
    baselineCuts.push_back(*lowDphiBaselineCut<RA2bTree>);
  }
  else assert(1);


  skimSamples skims = *skims_;

  typedef plot<RA2bTree> plot;

  double mJbins[4]={50.,85.,135.,250.};
  // double METbins[4]={300.,500.,700.,1000.};
  double METbins[5]={150.,200.,300.,400.,600.};


  //vector<plot> tempPlots;
  plot MET_Plot(*fillMET<RA2bTree>,"MET","m_{J} [GeV]",4,METbins);


  vector<plot> fourbSRHighPlots;
  fourbSRHighPlots.push_back(plot(MET_Plot));

  vector<plot> fourbSRLowPlots;
  fourbSRLowPlots.push_back(plot(MET_Plot));


  vector<plot> fourbSBHighPlots;
  fourbSBHighPlots.push_back(plot(MET_Plot));

  vector<plot> fourbSBLowPlots;
  fourbSBLowPlots.push_back(plot(MET_Plot));


  vector<plot> threebSRHighPlots;
  threebSRHighPlots.push_back(plot(MET_Plot));

  vector<plot> threebSRLowPlots;
  threebSRLowPlots.push_back(plot(MET_Plot));


  vector<plot> threebSBHighPlots;
  threebSBHighPlots.push_back(plot(MET_Plot));

  vector<plot> threebSBLowPlots;
  threebSBLowPlots.push_back(plot(MET_Plot));


  vector<plot> twobSRHighPlots;
  twobSRHighPlots.push_back(plot(MET_Plot));

  vector<plot> twobSRLowPlots;
  twobSRLowPlots.push_back(plot(MET_Plot));


  vector<plot> twobSBHighPlots;
  twobSBHighPlots.push_back(plot(MET_Plot));

  vector<plot> twobSBLowPlots;
  twobSBLowPlots.push_back(plot(MET_Plot));


  TFile* outputFile;
  TString regionName;
  TString cutName="";
  if (looseCuts) cutName="_looseCuts";
  if (region == 0) regionName="";
  if (region == 1) regionName="_signalOnly";
  if (region == 2) regionName="_singleMu";
  if (region == 3) regionName="_singleEle";
  if (region == 4) regionName="_lowDphi";

  outputFile = new TFile("ALPHABETRes"+Year+"_V17"+regionName+".root","RECREATE");


  // background MC samples - 0 lepton regions
  for (int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++){
    RA2bTree* ntuple = skims.ntuples[iSample];

    // for (int iBin = 0 ; iBin < numMETbins ; iBin++){
    //   for (int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
    //     plots[iBin][iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
    //     plots[iBin][iPlot].setFillColor(ntuple,skims.fillColor[iSample]);
    //   }
    // }
    for (int i = 0 ; i < fourbSRHighPlots.size() ; i++) {
      fourbSRHighPlots[i].addNtuple(ntuple,"fourbSRHigh_"+skims.sampleName[iSample]);
      fourbSRHighPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
      fourbSRLowPlots[i].addNtuple(ntuple,"fourbSRLow_"+skims.sampleName[iSample]);
      fourbSRLowPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);

      fourbSBHighPlots[i].addNtuple(ntuple,"fourbSBHigh_"+skims.sampleName[iSample]);
      fourbSBHighPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
      fourbSBLowPlots[i].addNtuple(ntuple,"fourbSBLow_"+skims.sampleName[iSample]);
      fourbSBLowPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);

      threebSRHighPlots[i].addNtuple(ntuple,"threebSRHigh_"+skims.sampleName[iSample]);
      threebSRHighPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
      threebSRLowPlots[i].addNtuple(ntuple,"threebSRLow_"+skims.sampleName[iSample]);
      threebSRLowPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);

      threebSBHighPlots[i].addNtuple(ntuple,"threebSBHigh_"+skims.sampleName[iSample]);
      threebSBHighPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
      threebSBLowPlots[i].addNtuple(ntuple,"threebSBLow_"+skims.sampleName[iSample]);
      threebSBLowPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);

      twobSRHighPlots[i].addNtuple(ntuple,"twobSRHigh_"+skims.sampleName[iSample]);
      twobSRHighPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
      twobSRLowPlots[i].addNtuple(ntuple,"twobSRLow_"+skims.sampleName[iSample]);
      twobSRLowPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);

      twobSBHighPlots[i].addNtuple(ntuple,"twobSBHigh_"+skims.sampleName[iSample]);
      twobSBHighPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
      twobSBLowPlots[i].addNtuple(ntuple,"twobSBLow_"+skims.sampleName[iSample]);
      twobSBLowPlots[i].setFillColor(ntuple,skims.fillColor[iSample]);
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
    else if ( filename.Contains("2017") ) this_lumi = 101269.0;

    // else if ( filename.Contains("2017") ) this_lumi = 41529.0;
    // else if ( filename.Contains("2018") ) this_lumi = 59740.0;

    if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH")) {
      this_lumi = 137191.0;
      TFile *fin = new TFile(filename,"READ");
      TH1F *nEventsHisto = (TH1F*)fin->Get("nEventProc");
      TotalEvents = nEventsHisto->GetBinContent(1);
    }

    for (int iEvt = 0 ; iEvt < numEvents ; iEvt++) {
      ntuple->GetEntry(iEvt);
      if (iEvt % 10000 == 0) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << min(MAX_EVENTS,numEvents) << endl;

      if (region==0 || region==1) {
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
        if(skims.sampleName[iSample]=="QCD") trigWeight=trigcorrorFakeMHT.GetCorrection(ntuple->MHT,trigunc);
      }
      else if (region == 2) {
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      }
      else if (region == 3) {
        trigWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
        trigWeight=trigWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      }
      else if (region == 4) {
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

      // if ((filename.Contains("SingleLept") || filename.Contains("DiLept")) && ntuple->madHT>600.) continue;
      if (( filename.Contains("SingleLeptFromT_MC") || filename.Contains("SingleLeptFromTbar_MC")  || filename.Contains("DiLept_MC") ) && ntuple->GenMET>150. )continue;


      weight = ntuple->Weight*this_lumi*trigWeight;
      if (filename.Contains("TChiHH_HToBB")) weight = weight*0.5823329*0.5823329/TotalEvents;
      else if (filename.Contains("T5qqqqZH")) weight = weight/TotalEvents*4;

      //if (skims.sampleName[iSample] == "TT") {
      //    weight *= ISRweights(ntuple);
      //}
      bin = -1;

      // for (int iBin = 0 ; iBin < numMETbins ; iBin++) {
      //   if (ntuple->MET > lowestMET) {
      //     if (ntuple->MET > numMETbins*(binWidth-1)+lowestMET) bin = numMETbins-1;
      //     else bin = int((ntuple->MET-lowestMET)/binWidth);
      //   }
      // }

      //our res MET bins are: 150-200, 200-300, 300-400, 400+
      float thisMET_forBin = ntuple->MET;
      if (thisMET_forBin>150.&&thisMET_forBin<=200.) bin=0;
      else if (thisMET_forBin>200.&&thisMET_forBin<=300.) bin=1;
      else if (thisMET_forBin>300.&&thisMET_forBin<=400.) bin=2;
      else if (thisMET_forBin>400.) bin=3;

      if (bin < 0) continue;

      if (fourbSRCut(ntuple)) {
        // plots[bin][4].fill(ntuple,weight);
        if (deltaRLow(ntuple)) {
          for (int i = 0 ; i < fourbSRLowPlots.size() ; i++) {
            fourbSRLowPlots[i].fill (ntuple,weight);
          }
        }
        else if (deltaRHigh(ntuple)) {
          for (int i = 0 ; i < fourbSRHighPlots.size() ; i++) {
            fourbSRHighPlots[i].fill (ntuple,weight);
          }
        }
      }

      else if (fourbSBCut(ntuple)) {
        // plots[bin][5].fill(ntuple,weight);

        if (deltaRLow(ntuple)) {
          for (int i = 0 ; i < fourbSBLowPlots.size() ; i++) {
            fourbSBLowPlots[i].fill (ntuple,weight);
          }
        }
        else if (deltaRHigh(ntuple)) {
          for (int i = 0 ; i < fourbSBHighPlots.size() ; i++) {
            fourbSBHighPlots[i].fill (ntuple,weight);
          }
        }
      }

      else if (threebSRCut(ntuple)) {
        // plots[bin][0].fill(ntuple,weight);

        if (deltaRLow(ntuple)) {
          for (int i = 0 ; i < threebSRLowPlots.size() ; i++) {
            threebSRLowPlots[i].fill (ntuple,weight);
          }
        }
        else if (deltaRHigh(ntuple)) {
          for (int i = 0 ; i < threebSRHighPlots.size() ; i++) {
            threebSRHighPlots[i].fill (ntuple,weight);
          }
        }
      }


      else if (threebSBCut(ntuple)) {
        // plots[bin][1].fill(ntuple,weight);
        if (deltaRLow(ntuple)) {
          for (int i = 0 ; i < threebSBLowPlots.size() ; i++) {
            threebSBLowPlots[i].fill (ntuple,weight);
          }
        }
        else if (deltaRHigh(ntuple)) {
          for (int i = 0 ; i < threebSBHighPlots.size() ; i++) {
            threebSBHighPlots[i].fill (ntuple,weight);
          }
        }
      }



      else if (twobSRCut(ntuple)) {
        // plots[bin][2].fill(ntuple,weight);
        if (deltaRLow(ntuple)) {
          for (int i = 0 ; i < twobSRLowPlots.size() ; i++) {
            twobSRLowPlots[i].fill (ntuple,weight);
          }
        }
        else if (deltaRHigh(ntuple)) {
          for (int i = 0 ; i < twobSRHighPlots.size() ; i++) {
            twobSRHighPlots[i].fill (ntuple,weight);
          }
        }
      }


      else if (twobSBCut(ntuple)) {
        // plots[bin][3].fill(ntuple,weight);
        if (deltaRLow(ntuple)) {
          for (int i = 0 ; i < twobSBLowPlots.size() ; i++) {
            twobSBLowPlots[i].fill (ntuple,weight);
          }
        }
        else if (deltaRHigh(ntuple)) {
          for (int i = 0 ; i < twobSBHighPlots.size() ; i++) {
            twobSBHighPlots[i].fill (ntuple,weight);
          }
        }
      }


    }// end event loop
    outputFile->cd();
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



  bool sumBkgs = true;
  if (region == 1 ) sumBkgs=false;

  // for (int iBin = 0 ; iBin < numMETbins; iBin++){
  //   for (int iPlot = 0 ; iPlot < plots[iBin].size() ; iPlot++){
  //     outputFile->cd();
  //     plots[iBin][iPlot].buildSum();
  //     plots[iBin][iPlot].Write();
  //     if (sumBkgs) plots[iBin][iPlot].sum->Write();
  //   }
  // }

  for (int i = 0 ; i < fourbSRLowPlots.size() ; i++) {
    outputFile->cd();
    fourbSRLowPlots[i].buildSum("fourbSRLow");
    fourbSRLowPlots[i].Write();
    if (sumBkgs) fourbSRLowPlots[i].sum->Write();

    fourbSRHighPlots[i].buildSum("fourbSRHigh");
    fourbSRHighPlots[i].Write();
    if (sumBkgs) fourbSRHighPlots[i].sum->Write();
  }

  for (int i = 0 ; i < fourbSBLowPlots.size() ; i++) {
    outputFile->cd();
    fourbSBLowPlots[i].buildSum("fourbSBLow");
    fourbSBLowPlots[i].Write();
    if (sumBkgs) fourbSBLowPlots[i].sum->Write();

    fourbSBHighPlots[i].buildSum("fourbSBHigh");
    fourbSBHighPlots[i].Write();
    if (sumBkgs) fourbSBHighPlots[i].sum->Write();
  }

  for (int i = 0 ; i < threebSRLowPlots.size() ; i++) {
    outputFile->cd();
    threebSRLowPlots[i].buildSum("threebSRLow");
    threebSRLowPlots[i].Write();
    if (sumBkgs) threebSRLowPlots[i].sum->Write();

    threebSRHighPlots[i].buildSum("threebSRHigh");
    threebSRHighPlots[i].Write();
    if (sumBkgs) threebSRHighPlots[i].sum->Write();
  }

  for (int i = 0 ; i < threebSBLowPlots.size() ; i++) {
    outputFile->cd();
    threebSBLowPlots[i].buildSum("threebSBLow");
    threebSBLowPlots[i].Write();
    if (sumBkgs) threebSBLowPlots[i].sum->Write();

    threebSBHighPlots[i].buildSum("threebSBHigh");
    threebSBHighPlots[i].Write();
    if (sumBkgs) threebSBHighPlots[i].sum->Write();
  }

  for (int i = 0 ; i < twobSRLowPlots.size() ; i++) {
    outputFile->cd();
    twobSRLowPlots[i].buildSum("twobSRLow");
    twobSRLowPlots[i].Write();
    if (sumBkgs) twobSRLowPlots[i].sum->Write();

    twobSRHighPlots[i].buildSum("twobSRHigh");
    twobSRHighPlots[i].Write();
    if (sumBkgs) twobSRHighPlots[i].sum->Write();
  }

  for (int i = 0 ; i < twobSBLowPlots.size() ; i++) {
    outputFile->cd();
    twobSBLowPlots[i].buildSum("twobSBLow");
    twobSBLowPlots[i].Write();
    if (sumBkgs) twobSBLowPlots[i].sum->Write();

    twobSBHighPlots[i].buildSum("twobSBHigh");
    twobSBHighPlots[i].Write();
    if (sumBkgs) twobSBHighPlots[i].sum->Write();
  }

  outputFile->Close();
}
