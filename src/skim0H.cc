#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "THStack.h"
#include "TString.h"
#include "TPad.h"
#include "TLorentzVector.h"

#include <vector>
#include <map>
#include <iostream>
#include <assert.h>

// #include "plotterUtils.cc"
#include "skimSamplesTogether.cc"
#include "definitions.cc"
#include "RA2bTree.cc"
#include "TriggerCorrector.h"

using namespace std;

int main(int argc, char** argv) {
  int region(0);
  region = atoi(argv[1]);
  TString Year(argv[2]);


  // reg = static_cast<skimSamples::region>(reg_);
  skimSamplesTogether* skims_;
  if (region == 0 ) skims_ = new skimSamplesTogether(skimSamplesTogether::kSignal, Year);
  else if (region == 1) {
    skims_ = new skimSamplesTogether(skimSamplesTogether::kSignalOnly, Year);
  }
  else if (region == 2 ) skims_ = new skimSamplesTogether(skimSamplesTogether::kSLm, Year);
  else if (region == 3 ) skims_ = new skimSamplesTogether(skimSamplesTogether::kSLe, Year);
  else if (region == 4 ) skims_ = new skimSamplesTogether(skimSamplesTogether::kPhoton, Year);
  else if (region == 5 ) skims_ = new skimSamplesTogether(skimSamplesTogether::kLowDphi, Year);
  else assert(1);

  TriggerCorrector trigcorror; TriggerCorrector trigcorrorHT; TriggerCorrector trigcorrorFakeMHT;
	TriggerCorrector trigcorrorPhoBarrel; TriggerCorrector trigcorrorPhoEndcap;

	if (Year=="MC2016") {
    trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2016");
    trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2016");
    trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2016");
    trigcorrorPhoBarrel.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonBarrelLoose_hists_Run2016_JetHT");
    trigcorrorPhoEndcap.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonEndcapLoose_hists_Run2016_JetHT");
 	}
	if (Year=="MC2017") {
    trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2017");
    trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2017");
    trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2017");
    trigcorrorPhoBarrel.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonBarrelLoose_hists_Run2017_JetHT");
    trigcorrorPhoEndcap.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonEndcapLoose_hists_Run2017_JetHT");
  }
	if (Year=="MC2018") {
    trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2018");
    trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2018");
    trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2018");
    trigcorrorPhoBarrel.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonBarrelLoose_hists_Run2018_JetHT");
    trigcorrorPhoEndcap.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root", "teff_SinglePhotonEndcapLoose_hists_Run2018_JetHT");
  }

  skimSamplesTogether skims = *skims_;
  TString regionName;

  if (region == 0 ) regionName = "signal";
  else if (region == 1 ) regionName = regionName = "signalOnly";
  else if (region == 2 ) regionName = regionName = "singleMu";
  else if (region == 3 ) regionName = regionName = "singleEle";
  else if (region == 4 ) regionName = regionName = "photon";
  else if (region == 5 ) regionName = regionName = "lowDPhi";

  // regionName="signal";

  TFile* outputFile = new TFile("testing_"+Year+"_"+regionName+".root","RECREATE");

  // background MC samples - 0 lepton regions
  for ( int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++) {

    cout << skims.sampleName[iSample] << endl;
		RA2bTree* ntuple = skims.ntuples[iSample];
    int numEvents = ntuple->fChain->GetEntries();
		ntupleBranchStatus<RA2bTree>(ntuple);

    TTree* newtree = new TTree(skims.sampleName[iSample],skims.sampleName[iSample]);
    int nAK4=0; //these have pT>30 and ets<2.4 cuts
    int nAK8=0; //these have pT>170 and ets<2.4 cuts
    int nISR=0;
    int nGenZs = 0; int nGenHs = 0;
    int BTagsL=0; int BTagsM=0; int BTagsT=0;
    float MET = 0; float HT = 0; float MHT = 0;
    float Weight = 0; float lumi = 0;
    float J1_mass = 0; float J2_mass = 0;
    float J1_pt = 0; float J2_pt = 0;
    float J1AK4_pt = 0; float J2AK4_pt = 0;
    float J1_deepbbtag = 0; float J2_deepbbtag = 0;
    int isSR = 0;
    int isSB = 0;
    int isSB2 = 0;
    int isSBBoth = 0;
    int isDoubleTag = 0;
    int isAntiTag = 0; //this should allow me to do that one test later
    float photonPt=0;

    TBranch *b_nAK4, *b_nAK8;
    TBranch *b_nISR;
    TBranch *b_nGenZs, *b_nGenHs;
    TBranch *b_BTagsL, *b_BTagsM, *b_BTagsT;
    TBranch *b_MET, *b_HT, *b_MHT;
    TBranch *b_Weight, *b_lumi;
    TBranch *b_J1_mass, *b_J2_mass, *b_J1_pt, *b_J2_pt;
    TBranch *b_J1AK4_pt, *b_J2AK4_pt;
    TBranch *b_J1_deepbbtag, *b_J2_deepbbtag;
    TBranch *b_isSR, *b_isSB, *b_isSB2, *b_isSBBoth;
    TBranch *b_isDoubleTag, *b_isAntiTag;
    TBranch *b_photonPt;

    newtree->Branch("nAK4", &nAK4, "nAK4/I");
    newtree->Branch("nAK8", &nAK8, "nAK8/I");
    newtree->Branch("nISR", &nISR, "nISR/I");
    newtree->Branch("MET",&MET,"MET/F");
    newtree->Branch("HT",&HT,"HT/F");
    newtree->Branch("MHT",&MHT,"MHT/F");
    newtree->Branch("Weight", &Weight, "Weight/F");
    newtree->Branch("Lumi", &lumi, "lumi/F");
    newtree->Branch("nGenZs", &nGenZs, "nGenZs/I");
    newtree->Branch("nGenHs", &nGenHs, "nGenHs/I");
    newtree->Branch("BTagsL", &BTagsL, "BTagsL/I");
    newtree->Branch("BTagsM", &BTagsM, "BTagsM/I");
    newtree->Branch("BTagsT", &BTagsT, "BTagsT/I");
    newtree->Branch("J1_mass", &J1_mass, "J1_mass/F");
    newtree->Branch("J2_mass", &J2_mass, "J2_mass/F");
    newtree->Branch("J1_pt", &J1_pt, "J1_pt/F");
    newtree->Branch("J2_pt", &J2_pt, "J2_pt/F");
    newtree->Branch("J1AK4_pt", &J1AK4_pt, "J1AK4_pt/F");
    newtree->Branch("J2AK4_pt", &J2AK4_pt, "J2AK4_pt/F");
    newtree->Branch("J1_deepbbtag", &J1_deepbbtag, "J1_deepbbtag/F");
    newtree->Branch("J2_deepbbtag", &J2_deepbbtag, "J2_deepbbtag/F");
    newtree->Branch("isSR", &isSR, "isSR/I");
    newtree->Branch("isSB", &isSB, "isSB/I");
    newtree->Branch("isSB2", &isSB2, "isSB2/I");
    newtree->Branch("isSBBoth", &isSBBoth, "isSBBoth/I");
    newtree->Branch("isDoubleTag", &isDoubleTag, "isDoubleTag/I");
    newtree->Branch("isAntiTag", &isAntiTag, "isAntiTag/I");
    if (region==4) newtree->Branch("photonPt", &photonPt, "photonPt/F");

    newtree->SetBranchAddress("nAK4",&nAK4, &b_nAK4);
    newtree->SetBranchAddress("nAK8",&nAK8, &b_nAK8);
    newtree->SetBranchAddress("nISR",&nISR, &b_nISR);
    newtree->SetBranchAddress("MET", &MET, &b_MET);
    newtree->SetBranchAddress("HT", &HT, &b_HT);
    newtree->SetBranchAddress("MHT", &MHT, &b_MHT);
    newtree->SetBranchAddress("Weight", &Weight, &b_Weight);
    newtree->SetBranchAddress("Lumi", &lumi, &b_lumi);
    newtree->SetBranchAddress("nGenZs",&nGenZs, &b_nGenZs);
    newtree->SetBranchAddress("nGenHs",&nGenHs, &b_nGenHs);
    newtree->SetBranchAddress("BTagsL", &BTagsL, &b_BTagsL);
    newtree->SetBranchAddress("BTagsM", &BTagsM, &b_BTagsM);
    newtree->SetBranchAddress("BTagsT", &BTagsT, &b_BTagsT);
    newtree->SetBranchAddress("J1_mass", &J1_mass, &b_J1_mass);
    newtree->SetBranchAddress("J2_mass", &J2_mass, &b_J2_mass);
    newtree->SetBranchAddress("J1_pt", &J1_pt, &b_J1_pt);
    newtree->SetBranchAddress("J2_pt", &J2_pt, &b_J2_pt);
    newtree->SetBranchAddress("J1AK4_pt", &J1AK4_pt, &b_J1AK4_pt);
    newtree->SetBranchAddress("J2AK4_pt", &J2AK4_pt, &b_J2AK4_pt);
    newtree->SetBranchAddress("J1_deepbbtag", &J1_deepbbtag, &b_J1_deepbbtag);
    newtree->SetBranchAddress("J2_deepbbtag", &J2_deepbbtag, &b_J2_deepbbtag);
    newtree->SetBranchAddress("isSR", &isSR, &b_isSR);
    newtree->SetBranchAddress("isSB", &isSB, &b_isSB);
    newtree->SetBranchAddress("isSB2", &isSB2, &b_isSB2);
    newtree->SetBranchAddress("isSBBoth", &isSBBoth, &b_isSBBoth);
    newtree->SetBranchAddress("isDoubleTag", &isDoubleTag, &b_isDoubleTag);
    newtree->SetBranchAddress("isAntiTag", &isAntiTag, &b_isAntiTag);
    if (region==4) newtree->SetBranchAddress("photonPt", &photonPt, &b_photonPt);

    int TotalEvents = 0;
    TString filename = ntuple->fChain->GetFile()->GetName();
    float triggerWeight=1.0; float trigunc=0;

    if ( filename.Contains("2016") ) lumi = 35922.0;
    // else if ( filename.Contains("2017") ) lumi = 101269.0; //2017+2018
    else if ( filename.Contains("2017") ) lumi = 41529.0;
    else if ( filename.Contains("2018") ) lumi = 59740.0;

    if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH") ) {
      TFile *fin = new TFile(filename,"READ");
      TH1F *nEventsHisto = (TH1F*)fin->Get("nEventProc");
      TotalEvents = nEventsHisto->GetBinContent(1);
      std::cout<<"Total events: "<<TotalEvents<<std::endl;
    }

    for ( int iEvt = 0 ; iEvt < numEvents ; iEvt++ ) {
      ntuple->GetEntry(iEvt);
      if ( iEvt % 100000 == 0 ) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << numEvents << endl;

      if ((region == 0 ||region==1) && !boostedBaselineCut(ntuple) ) continue;
      else if (region == 2 && !singleMuBaselineCut(ntuple) ) continue;
      else if (region == 3 && !singleEleBaselineCut(ntuple) ) continue;
      else if (region == 4 && !photonBaselineCut(ntuple) ) continue;

      // if ( ( filename.Contains("SingleLept") || filename.Contains("DiLept") ) && ntuple->madHT>600. ) continue;
      // if (( filename.Contains("SingleLeptFromT_MC") || filename.Contains("SingleLeptFromTbar_MC")  || filename.Contains("DiLept_MC") ) && ntuple->GenMET>150. )continue;
      if (filename.Contains("MC2018")) {
        if (( filename.Contains("SingleLeptFromT_MC") || filename.Contains("SingleLeptFromTbar_MC")  || filename.Contains("DiLept_MC") ) && ntuple->GenMET>80.) continue;
      }
      else {
        if (( filename.Contains("SingleLeptFromT_MC") || filename.Contains("SingleLeptFromTbar_MC")  || filename.Contains("DiLept_MC") ) && ntuple->GenMET>150.) continue;
      }

      triggerWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
      triggerWeight=triggerWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
      if (filename.Contains("QCD")) triggerWeight=trigcorrorFakeMHT.GetCorrection(ntuple->MHT,trigunc);


      nAK4=0; nAK8=0;
      BTagsL=0; BTagsM=0; BTagsT=0;
      nGenZs = 0; nGenHs = 0;
      MET = 0; HT = 0; MHT = 0;
      Weight = 0;
      J1_mass = 0; J2_mass = 0;
      J1_pt = 0; J2_pt = 0;
      J1AK4_pt = 0; J2AK4_pt = 0;
      J1_deepbbtag = 0; J2_deepbbtag = 0;
      isSR = 0; isSB = 0; isSB2 = 0; isSBBoth = 0;
      isDoubleTag = 0; isAntiTag = 0;
      photonPt = 0;

      nGenHs = getNumGenHiggses(ntuple);
      if (filename.Contains("T5qqqqZH") && nGenHs!=2) continue;

      MET=ntuple->MET; HT=ntuple->HT; MHT=ntuple->MHT;
      Weight=ntuple->Weight*triggerWeight;
      if ( filename.Contains("TChiHH_HToBB") ) Weight = Weight*0.5823329*0.5823329/TotalEvents;
      else if ( filename.Contains("T5qqqqZH") )  Weight = Weight/TotalEvents*4.0; //times 4 because only using 1/4 of events

      nGenZs = getNumGenZs(ntuple);
      if (region==4) photonPt = ntuple->Photons->at(0).Pt();

      nAK4 = 0; nAK8 = 0;
      nISR = ntuple->NJetsISR;

      vector<int> thisNBs = numDeepBs(ntuple);
      BTagsL=thisNBs.at(0); BTagsM=thisNBs.at(1); BTagsT=thisNBs.at(2);
      nAK4 = numJets(ntuple);

      for (int i=0; i<ntuple->JetsAK8->size();i++) {
        if (ntuple->JetsAK8->at(i).Pt()<170.0 || abs(ntuple->JetsAK8->at(i).Eta())>2.4 || ntuple->JetsAK8_softDropMass->at(i)<1.0) continue;
        nAK8++;
      }

      J1_mass = ntuple->JetsAK8_softDropMass->at(0);
      J2_mass = ntuple->JetsAK8_softDropMass->at(1);
      J1_pt = ntuple->JetsAK8->at(0).Pt();
      J2_pt = ntuple->JetsAK8->at(1).Pt();
      J1AK4_pt = ntuple->Jets->at(0).Pt();
      J2AK4_pt = ntuple->Jets->at(1).Pt();
      J1_deepbbtag = ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(0);
      J2_deepbbtag = ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(1);


      if (doubleMassCut(ntuple)) isSR = 1;
      else {
        if (SBMassCut(ntuple)) isSB = 1;
        if (SB2MassCut(ntuple)) isSB2 = 1;
        if (SBBothMassCut(ntuple)) isSBBoth = 1;
      }

      if (doubletagSRCut(ntuple) || doubletagSBCut(ntuple)) isDoubleTag = 1;
      else if (antitagSRCut(ntuple) || antitagSBCut(ntuple)) isAntiTag = 1;

      newtree->Fill();
    } //end event loop

    outputFile->cd();
    newtree->Write(skims.sampleName[iSample]);
    delete newtree;
  } //end sample loop

  outputFile->Close();
  return 0;
}
