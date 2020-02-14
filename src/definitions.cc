#include "TLorentzVector.h"
#include "TRandom3.h"
// constants
// ==============================================
double bbtagCut_unsupp = 0.70; //this is for deep bb tag - 2016 unsupported, matching signal efficiency to real deep (Med is 0.89)
double deepBBTagCut = 0.7; //(Med is 0.86 or 0.89)
double bbtagCut = 0.3; //this is for nominal MVA bb tag  (Med is 0.6)
TFile* puWeightFile = new TFile("../data/PileupHistograms_0121_69p2mb_pm4p6.root");
TH1F* puWeightHist = (TH1F*) puWeightFile->Get("pu_weights_down");
// - - - - - - weights for WJets, GJets, - - - - - - - -
// - - - - - - and ZJets NLO Pt distribution - - - - - -
TFile* NLOWeightFile = new TFile("../data/kfactors.root");
TH1F* GJets_NLO = (TH1F*) NLOWeightFile->Get("GJets_1j_NLO/nominal_G");
TH1F* GJets_LO = (TH1F*) NLOWeightFile->Get("GJets_LO/inv_pt_G");
TH1F* WJets_NLO = (TH1F*) NLOWeightFile->Get("WJets_012j_NLO/nominal");
TH1F* WJets_LO = (TH1F*) NLOWeightFile->Get("WJets_LO/inv_pt");
TH1F* ZJets_NLO = (TH1F*) NLOWeightFile->Get("ZJets_01j_NLO/nominal");
TH1F* ZJets_LO = (TH1F*) NLOWeightFile->Get("ZJets_LO/inv_pt");
// ==============================================

double CalcdPhi(double phi1 , double phi2) {
  double dPhi = phi1-phi2;
  if (dPhi < -TMath::Pi())
  dPhi += 2*TMath::Pi();
  if (dPhi > TMath::Pi())
  dPhi -= 2*TMath::Pi();
  return fabs(dPhi);
}

template<typename ntupleType>void ntupleBranchStatus(ntupleType* ntuple) {
  ntuple->fChain->SetBranchStatus("*",0);
  ntuple->fChain->SetBranchStatus("Muons",1);
  ntuple->fChain->SetBranchStatus("Electrons",1);
  ntuple->fChain->SetBranchStatus("NMuons",1);
  ntuple->fChain->SetBranchStatus("NElectrons",1);
  ntuple->fChain->SetBranchStatus("isoElectronTracks",1);
  ntuple->fChain->SetBranchStatus("isoMuonTracks",1);
  ntuple->fChain->SetBranchStatus("isoPionTracks",1);
  ntuple->fChain->SetBranchStatus("Photon*",1);
  ntuple->fChain->SetBranchStatus("DeltaPhi*",1);

  ntuple->fChain->SetBranchStatus("MHT",1);
  ntuple->fChain->SetBranchStatus("HT",1);
  ntuple->fChain->SetBranchStatus("NJets",1);
  ntuple->fChain->SetBranchStatus("BTags",1);
  ntuple->fChain->SetBranchStatus("MET",1);
  ntuple->fChain->SetBranchStatus("MHT",1);
  ntuple->fChain->SetBranchStatus("METPhi",1);

  ntuple->fChain->SetBranchStatus("JetsAK8*",1);
  ntuple->fChain->SetBranchStatus("Jets*",1);
  ntuple->fChain->SetBranchStatus("Weight",1);
  ntuple->fChain->SetBranchStatus("TrueNumInteractions",1);
  ntuple->fChain->SetBranchStatus("TriggerPass",1);
  ntuple->fChain->SetBranchStatus("CaloMET",1);
  ntuple->fChain->SetBranchStatus("NVtx",1);
  ntuple->fChain->SetBranchStatus("NumInteractions",1);
  ntuple->fChain->SetBranchStatus("nAllVertices",1);
  ntuple->fChain->SetBranchStatus("JetID*",1);
  ntuple->fChain->SetBranchStatus("madHT",1);
  ntuple->fChain->SetBranchStatus("NJetsISR",1);
  ntuple->fChain->SetBranchStatus("madMinDeltaRStatus",1);
  ntuple->fChain->SetBranchStatus("madMinPhotonDeltaR",1);
  ntuple->fChain->SetBranchStatus("GenParticles*",1);

  ntuple->fChain->SetBranchStatus("HBHENoiseFilter",1);
  ntuple->fChain->SetBranchStatus("eeBadScFilter",1);
  ntuple->fChain->SetBranchStatus("NVtx",1);
  ntuple->fChain->SetBranchStatus("BadPFMuonFilter",1);
  ntuple->fChain->SetBranchStatus("globalSuperTightHalo2016Filter",1);
  ntuple->fChain->SetBranchStatus("LowNeutralJetFilter",1);
  ntuple->fChain->SetBranchStatus("HTRatioDPhiTightFilter",1);
  ntuple->fChain->SetBranchStatus("FakeJetFilter",1);
}

/***************************************************************/
/* - - - - - - - - - - - - gen-level cuts - - - - - - - - - -  */
/***************************************************************/
template<typename ntupleType> bool genWmatched(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size() == 0) return false;

  for (int i=0 ; i < ntuple->GenParticles->size(); i++) {
    if (abs(ntuple->GenParticles_PdgId->at(i)) == 24 && ntuple->JetsAK8->at(0).DeltaR(ntuple->GenParticles->at(i))<0.4) return true;
  }
  return false;
}

template<typename ntupleType> bool genTmatched(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size() == 0) return false;

  for (int i=0 ; i < ntuple->GenParticles->size(); i++) {
    if (abs(ntuple->GenParticles_PdgId->at(i)) == 6 && ntuple->JetsAK8->at(0).DeltaR(ntuple->GenParticles->at(i))<0.4) return true;
  }
  return false;
}

template<typename ntupleType> int getNumGenHiggses(ntupleType* ntuple) {
  int numHiggses=0;
  for (int i=0 ; i < ntuple->GenParticles->size(); i++) {
    if (ntuple->GenParticles_PdgId->at(i) == 25 &&
    ntuple->GenParticles_ParentId->at(i) == 1000023 &&
    ntuple->GenParticles_Status->at(i) == 22)
    numHiggses++;
  }
  return numHiggses;
}

template<typename ntupleType> int getNumGenZs(ntupleType* ntuple) {
  int numZs=0;
  for (int i=0 ; i < ntuple->GenParticles->size(); i++) {
    if (ntuple->GenParticles_PdgId->at(i) == 23 &&
    ntuple->GenParticles_ParentId->at(i) == 1000023 &&
    ntuple->GenParticles_Status->at(i) == 22)
    numZs++;
  }
  return numZs;
}

template<typename ntupleType> bool genLevelHHcut(ntupleType* ntuple) {
  int numHiggses=getNumGenHiggses(ntuple),numZs=getNumGenZs(ntuple);
  if (numHiggses==2 && numZs==0) return true;
  else return false;
}

template<typename ntupleType> bool genLevelZHcut(ntupleType* ntuple) {
  int numHiggses=getNumGenHiggses(ntuple),numZs=getNumGenZs(ntuple);
  if (numHiggses==1 && numZs==1) return true;
  else return false;
}

template<typename ntupleType> bool genLevelZZcut(ntupleType* ntuple) {
  int numHiggses=getNumGenHiggses(ntuple),numZs=getNumGenZs(ntuple);
  if (numHiggses==0 && numZs==2) return true;
  else return false;
}

/***************************************************************/
/* - - - - - - - - - - - - custom weights - - - - - - - - - -  */
/***************************************************************/
template<typename ntupleType> double GJetsNLOWeights(ntupleType* ntuple) {
  if (ntuple->Photons->size() == 0) return 0.;
  double photon_pt = -999.;//ntuple->Photons->at(0).Pt();
  int photonIndex=-1;
  for (unsigned int p = 0 ; p < ntuple->GenParticles->size(); p++) {
    if (abs(ntuple->GenParticles_PdgId->at(p)) == 22) {
      if (photonIndex < 0)
      photonIndex = p;
      else if (ntuple->GenParticles->at(p).Pt() > ntuple->GenParticles->at(photonIndex).Pt())
      photonIndex = p;
    }
  }
  photon_pt = ntuple->GenParticles->at(photonIndex).Pt();

  if (photon_pt>150.) {
    double LO = GJets_LO->GetBinContent(GJets_LO->FindBin(photon_pt));
    double NLO = GJets_NLO->GetBinContent(GJets_NLO->FindBin(photon_pt));
    return (LO==0?0.:NLO/LO);
  }
  else return GJets_NLO->GetBinContent(1)/GJets_LO->GetBinContent(1);
}

template<typename ntupleType> double WJetsNLOWeights(ntupleType* ntuple) {
  double Wpt=-999.;
  for (unsigned int p = 0 ; p < ntuple->GenParticles->size(); p++) {
    if (abs(ntuple->GenParticles_PdgId->at(p)) == 24)
    Wpt = ntuple->GenParticles->at(p).Pt();
  }
  if (Wpt>150.) {
    double LO = WJets_LO->GetBinContent(WJets_LO->FindBin(Wpt));
    double NLO = WJets_NLO->GetBinContent(WJets_NLO->FindBin(Wpt));
    return (LO==0?0.:NLO/LO/1.21);
  }
  else return WJets_NLO->GetBinContent(1)/WJets_LO->GetBinContent(1)/1.21;
}

template<typename ntupleType> double ZJetsNLOWeights(ntupleType* ntuple) {
  double Zpt=-999.;
  for (unsigned int p = 0 ; p < ntuple->GenParticles->size(); p++) {
    if (abs(ntuple->GenParticles_PdgId->at(p)) == 23)
    Zpt = ntuple->GenParticles->at(p).Pt();
  }
  if (Zpt>150.) {
    double LO = ZJets_LO->GetBinContent(ZJets_LO->FindBin(Zpt));
    double NLO = ZJets_NLO->GetBinContent(ZJets_NLO->FindBin(Zpt));
    return (LO==0?0.:NLO/LO/1.23);
  }
  else return ZJets_NLO->GetBinContent(1)/ZJets_LO->GetBinContent(1)/1.23;
}

template<typename ntupleType> double singleMuonTrigWeights(ntupleType* ntuple) {
  if (ntuple->Muons->size() == 0) return 0.;
  else if (ntuple->HT > 300.) {
    if (ntuple->HT < 500.) {
      if (ntuple->Muons->at(0).Pt() > 25.) {
        if (ntuple->Muons->at(0).Pt() < 30.) {return .787;}
        else if (ntuple->Muons->at(0).Pt() < 50.) {return .843;}
        else {return .908;}
      }
      else {return 0.;}
    }
    else if (ntuple->HT > 500.) {return .949;}
    else {return 0.;}
  }
}

template<typename ntupleType> double singleElectronTrigWeights(ntupleType* ntuple) {
  if (ntuple->Electrons->size() == 0) return 0.;
  else if (ntuple->HT > 450.) {
    if (ntuple->Electrons->at(0).Pt() > 25.) {
      if (ntuple->Electrons->at(0).Pt() < 30.) {return 0.794;}
      else if (ntuple->Electrons->at(0).Pt() < 40.) {return 0.826;}
      else if (ntuple->Electrons->at(0).Pt() < 50.) {return 0.872;}
      else if (ntuple->Electrons->at(0).Pt() < 75.) {return 0.884;}
      else if (ntuple->Electrons->at(0).Pt() < 100.) {return 0.913;}
      else {return 0.947;}
    }
    else {return 0.;}
  }
  else if (ntuple->HT > 300.) {
    if (ntuple->Electrons->at(0).Pt() > 25.) {
      if (ntuple->Electrons->at(0).Pt() < 30.) {return 0.572;}
      else if (ntuple->Electrons->at(0).Pt() < 40.) {return 0.775;}
      else if (ntuple->Electrons->at(0).Pt() < 50.) {return 0.858;}
      else if (ntuple->Electrons->at(0).Pt() < 75.) {return 0.861;}
      else if (ntuple->Electrons->at(0).Pt() < 100.) {return 0.932;}
      else {return 1.;}
    }
    else {return 0.;}
  }
  else {return 0.;}
}

template<typename ntupleType> double lowDphiTrigWeights(ntupleType* ntuple) {
  if (ntuple->MET>100.) {
    if (ntuple->MET<200.) {return 0.500;}
    else if (ntuple->MET<300.) {return 0.712;}
    else if (ntuple->MET<400.) {return 0.806;}
    else if (ntuple->MET<500.) {return 0.874;}
    else if (ntuple->MET<700.) {return 0.866;}
    else {return 0.766;}
  }
  else return 0.;
}

template<typename ntupleType> double customPUweights(ntupleType* ntuple) {
  int nVtx = ntuple->TrueNumInteractions;
  return puWeightHist->GetBinContent(puWeightHist->GetXaxis()->FindBin(min(ntuple->TrueNumInteractions,puWeightHist->GetBinLowEdge(puWeightHist->GetNbinsX()+1))));
}

enum ISRweightType {kNom,kUp,kDn};
template<typename ntupleType> double ISRweights(ntupleType* ntuple, ISRweightType wType = kNom) {

  double wanted_w_isr=1.;
  double wanted_sys_isr[2] = {1.,1.};

  TString sample = ntuple->fChain->GetFile()->GetName();
  // these are taken from here:
  // https://github.com/manuelfs/babymaker/blob/3a57e1bace6c52832fe40e401cf37bc6b50923c3/bmaker/genfiles/src/change_weights.cxx#L156-L175
  // via Manuel Franco Sevilla
  if (sample.Contains("TTJets_HT-600to800")) {
    wanted_w_isr = 0.7838;
    wanted_sys_isr[0] = 0.8965;
    wanted_sys_isr[1] = 0.6604;
  }
  else if (sample.Contains("TTJets_HT-800to1200")) {
    wanted_w_isr = 0.7600;
    wanted_sys_isr[0] = 0.8851;
    wanted_sys_isr[1] = 0.6230;
  }
  else if (sample.Contains("TTJets_HT-1200to2500")) {
    wanted_w_isr = 0.7365;
    wanted_sys_isr[0] = 0.8739;
    wanted_sys_isr[1] = 0.5861;
  }
  else if (sample.Contains("TTJets_HT-2500toInf")) {
    wanted_w_isr = 0.7254;
    wanted_sys_isr[0] = 0.8686;
    wanted_sys_isr[1] = 0.5687;
  }
  else { //  if (sample.Contains("TTJets_SingleLept") or sample.Contains("TTJets_DiLept")) {
    // these numbers should really only be applied to the inclusive sample
    wanted_w_isr = 1.071;
    wanted_sys_isr[0] = 1.071;
    wanted_sys_isr[1] = 1.071;
  }

  double D;
  if (wType == kNom) D = wanted_w_isr;
  else D = wanted_sys_isr[wType-1];

  double w[6] = {0.920,0.821,0.715,0.662,0.561,0.511};
  if (ntuple->NJetsISR == 0) return D;
  else if (ntuple->NJetsISR >= 6) return w[5]*D;
  else return w[ntuple->NJetsISR]*D;
}

//////////////////////
// Lepton functions //
//////////////////////

template<typename ntupleType> double computeMuonMT(ntupleType* ntuple) {
  if (ntuple->Muons->size() == 0) return -9999.;
  double lepPt = ntuple->Muons->at(0).Pt();
  double lepPhi = ntuple->Muons->at(0).Phi();
  double MET = ntuple->MET;
  double METPhi = ntuple->METPhi;
  return sqrt(2*lepPt*MET * (1 - cos(METPhi-lepPhi)));
}

template<typename ntupleType> double computeElectronMT(ntupleType* ntuple) {
  if (ntuple->Electrons->size() == 0) return -9999.;
  double lepPt = ntuple->Electrons->at(0).Pt();
  double lepPhi = ntuple->Electrons->at(0).Phi();
  double MET = ntuple->MET;
  double METPhi = ntuple->METPhi;
  return sqrt(2*lepPt*MET * (1 - cos(METPhi-lepPhi)));
}


///////////////////////////
// RECO Muon definitions //
///////////////////////////
template<typename ntupleType> int numMuons(ntupleType* ntuple) {
  return ntuple->Muons->size();
}

template<typename ntupleType> double muonLeadJetdR(ntupleType* ntuple) {
  if (ntuple->Muons->size() == 1) {
    return ntuple->JetsAK8->size()>=1?ntuple->Muons->at(0).DeltaR(ntuple->JetsAK8->at(0)):-99.;
  }
  else return -99.;
}

template<typename ntupleType> double muonSubleadJetdR(ntupleType* ntuple) {
  if (ntuple->Muons->size() == 1) {
    return ntuple->JetsAK8->size()>=2?ntuple->Muons->at(0).DeltaR(ntuple->JetsAK8->at(1)):-99.;
  }
  else return -99.;
}

template<typename ntupleType> double leadJetMuondR_mass(ntupleType* ntuple) {
  if (ntuple->Muons->size() == 1) {
    if (ntuple->Muons->at(0).DeltaR(ntuple->JetsAK8->at(0)) < ntuple->Muons->at(0).DeltaR(ntuple->JetsAK8->at(1)))
    return ntuple->JetsAK8_softDropMass->at(1);
    else return ntuple->JetsAK8_softDropMass->at(0);
  }
  else return -99.;
}

template<typename ntupleType> double subleadJetMuondR_mass(ntupleType* ntuple) {
  if (ntuple->Muons->size() == 1) {
    if (ntuple->Muons->at(0).DeltaR(ntuple->JetsAK8->at(0)) < ntuple->Muons->at(0).DeltaR(ntuple->JetsAK8->at(1)))
    return ntuple->JetsAK8_softDropMass->at(0);
    else return ntuple->JetsAK8_softDropMass->at(1);
  }
  else return -99.;
}

template<typename ntupleType> double leadJetMuondR_bbdisc(ntupleType* ntuple) {
  if (ntuple->Muons->size() == 1) {
    if (ntuple->Muons->at(0).DeltaR(ntuple->JetsAK8->at(0)) < ntuple->Muons->at(0).DeltaR(ntuple->JetsAK8->at(1)))
    return ntuple->JetsAK8_doubleBDiscriminator->at(1);
    else return ntuple->JetsAK8_doubleBDiscriminator->at(0);
  }
  else return -99.;
}

template<typename ntupleType> double subleadJetMuondR_bbdisc(ntupleType* ntuple) {
  if (ntuple->Muons->size() == 1) {
    if (ntuple->Muons->at(0).DeltaR(ntuple->JetsAK8->at(0)) < ntuple->Muons->at(0).DeltaR(ntuple->JetsAK8->at(1)))
    return ntuple->JetsAK8_doubleBDiscriminator->at(0);
    else return ntuple->JetsAK8_doubleBDiscriminator->at(1);
  }
  else return -99.;
}

template<typename ntupleType> double fillNumMuons(ntupleType* ntuple) {
  return double(ntuple->Muons->size());
}

///////////////////////////////////////////////
// - - - - - - - RECO ELECTRONS - - - - - -  //
///////////////////////////////////////////////
template<typename ntupleType> double electronLeadJetdR(ntupleType* ntuple) {
  if (ntuple->Electrons->size() == 1) {
    return ntuple->Electrons->at(0).DeltaR(ntuple->JetsAK8->at(0));
  }
  else return -99.;
}

template<typename ntupleType> double electronSubleadJetdR(ntupleType* ntuple) {
  if (ntuple->Electrons->size() == 1) {
    return ntuple->Electrons->at(0).DeltaR(ntuple->JetsAK8->at(1));
  }
  else return -99.;
}

template<typename ntupleType> int numElectrons(ntupleType* ntuple) {
  return ntuple->Electrons->size();
}

template<typename ntupleType> double fillLepPt(ntupleType* ntuple) {
  if (ntuple->Electrons->size() > 0 && ntuple->Muons->size() == 0){
    return ntuple->Electrons->at(0).Pt();
  }
  else if (ntuple->Electrons->size() ==0 && ntuple->Muons->size() > 0){
    return ntuple->Muons->at(0).Pt();
  }
  else return -9999.;
}

template<typename ntupleType> double fillLepEta(ntupleType* ntuple) {
  if (ntuple->Electrons->size() > 0 && ntuple->Muons->size() == 0){
    return ntuple->Electrons->at(0).Eta();
  }
  else if (ntuple->Electrons->size() ==0 && ntuple->Muons->size() > 0){
    return ntuple->Muons->at(0).Eta();
  }
  else return -9999.;
}

template<typename ntupleType> double fillLepActivity(ntupleType* ntuple) {
  if (ntuple->Electrons->size() > 0 && ntuple->Muons->size() == 0){
    return ntuple->Electrons_MT2Activity->at(0);
  }
  else if (ntuple->Electrons->size() ==0 && ntuple->Muons->size() > 0){
    return ntuple->Muons_MT2Activity->at(0);
  }
  else return -9999.;
}

////////////////////////////////////////////////////////////
// - - - - - - - - EVENT LEVEL VARIABLES - - - - - - - -  //
////////////////////////////////////////////////////////////
template<typename ntupleType> double fillJetPt1(ntupleType* ntuple) {
  if (ntuple->Jets->size() >= 1) return ntuple->Jets->at(0).Pt();
  else return -999.;
}

template<typename ntupleType> double fillJetPt2(ntupleType* ntuple) {
  if (ntuple->Jets->size() >= 2) return ntuple->Jets->at(1).Pt();
  else return -999.;
}

template<typename ntupleType> double fillJetPt3(ntupleType* ntuple) {
  if (ntuple->Jets->size() >= 3) return ntuple->Jets->at(2).Pt();
  else return -999.;
}

template<typename ntupleType> double fillJetPt4(ntupleType* ntuple) {
  if (ntuple->Jets->size() >= 4) return ntuple->Jets->at(3).Pt();
  else return -999.;
}

template<typename ntupleType> double fillNVtx(ntupleType* ntuple) {
  return ntuple->NVtx;
}

template<typename ntupleType> double fillnAllVertices(ntupleType* ntuple) {
  return ntuple->nAllVertices;
}

template<typename ntupleType> double fillNumInteractions(ntupleType* ntuple) {
  return ntuple->NumInteractions;
}

template<typename ntupleType> double fillMadHT(ntupleType* ntuple) {
  return ntuple->madHT;
}

template<typename ntupleType> double fillDeltaPhi1(ntupleType* ntuple) {
  return ntuple->DeltaPhi1;
}

template<typename ntupleType> double fillDeltaPhi2(ntupleType* ntuple) {
  return ntuple->DeltaPhi2;
}

template<typename ntupleType> double fillDeltaPhi3(ntupleType* ntuple) {
  return ntuple->DeltaPhi3;
}

template<typename ntupleType> double fillDeltaPhi4(ntupleType* ntuple) {
  return ntuple->DeltaPhi4;
}


template<typename ntupleType> double fillHT(ntupleType* ntuple) {
  return ntuple->HT;
}

template<typename ntupleType> double fillMHT(ntupleType* ntuple) {
  return ntuple->MHT;
}

template<typename ntupleType> double fillMET(ntupleType* ntuple) {
  return ntuple->MET;
}

template<typename ntupleType> double fillMETRatio(ntupleType* ntuple) {
  return ntuple->MET/ntuple->MHT;
}

template<typename ntupleType> double fillOne(ntupleType* ntuple) {
  return 1.;
}

template<typename ntupleType> double fillNJets(ntupleType* ntuple) {
  return ntuple->NJets;
}

template<typename ntupleType> double fillBTags(ntupleType* ntuple) {
  return ntuple->BTags;
}


template<typename ntupleType> int numJets(ntupleType* ntuple) { //Returns the number jets within pT and eta
  int NJets = 0;
  for (unsigned int j=0; j<ntuple->Jets->size();++j) {
    if (ntuple->Jets->at(j).Pt()<30 || fabs(ntuple->Jets->at(j).Eta())>2.4) continue;
    NJets++;
  }
  return NJets;
}

//Using deepCSV, returns the number of b's based on loose, medium, and tight WPs
template<typename ntupleType> std::vector<int> numDeepBs(ntupleType* ntuple) {
  int NJets = numJets(ntuple);
  if (NJets<4 || NJets>5) return {0,0,0};
  double CSVBtagLoose = 0.2217;
  double CSVBtagMed   = 0.6321;
  double CSVBtagTight = 0.8953;

  int BTagsL = 0; int BTagsM = 0; int BTagsT = 0;
  for (unsigned int j=0; j<ntuple->Jets->size();++j) {
    if (ntuple->Jets->at(j).Pt()<30 || fabs(ntuple->Jets->at(j).Eta())>2.4) continue;
    float this_CSV_value = ntuple->Jets_bJetTagDeepCSVprobb->at(j)+ntuple->Jets_bJetTagDeepCSVprobbb->at(j);
    if (this_CSV_value > CSVBtagTight) BTagsT++;
    if (this_CSV_value > CSVBtagMed) BTagsM++;
    if (this_CSV_value > CSVBtagLoose) BTagsL++;
  }
  std::vector<int> thisNBs = {BTagsL, BTagsM, BTagsT};
  return thisNBs;
}


template<typename ntupleType> std::vector<int> resHCandidates(ntupleType* ntuple) { //returns the indices for the 4 AK4 jets that make up H candidates (1,2) and (3,4)
  vector<int> thisNBs = numDeepBs(ntuple); int numTight = thisNBs.at(2);
  int NJets = numJets(ntuple);
  if (NJets<4 || NJets>5 || numTight<2) return {-11,-11,-11,-11};
  double HighestValuesTest[] = {-11.0,-11.0,-11.0,-11.0};
  int JetIndices[] = {-1,-1,-1,-1};
  for (unsigned int j=0; j<ntuple->Jets->size();++j) {
    if (ntuple->Jets->at(j).Pt()<30 || fabs(ntuple->Jets->at(j).Eta())>2.4) continue;

    double *current_min = min_element(HighestValuesTest,HighestValuesTest+4);
    int min_pos = distance(HighestValuesTest,min_element(HighestValuesTest,HighestValuesTest+4));

    double CSVofJet = (ntuple->Jets_bJetTagDeepCSVprobb->at(j)+ntuple->Jets_bJetTagDeepCSVprobbb->at(j));
    if (CSVofJet > *current_min) {
      HighestValuesTest[min_pos] = CSVofJet;
      JetIndices[min_pos] = j;
    }
  } //end loop to find jets with 4 highest CSV

  TLorentzVector JetComboTest1a = ntuple->Jets->at(JetIndices[0])+ntuple->Jets->at(JetIndices[1]);
  TLorentzVector JetComboTest1b = ntuple->Jets->at(JetIndices[2])+ntuple->Jets->at(JetIndices[3]);
  TLorentzVector JetComboTest2a = ntuple->Jets->at(JetIndices[0])+ntuple->Jets->at(JetIndices[2]);
  TLorentzVector JetComboTest2b = ntuple->Jets->at(JetIndices[1])+ntuple->Jets->at(JetIndices[3]);
  TLorentzVector JetComboTest3a = ntuple->Jets->at(JetIndices[0])+ntuple->Jets->at(JetIndices[3]);
  TLorentzVector JetComboTest3b = ntuple->Jets->at(JetIndices[1])+ntuple->Jets->at(JetIndices[2]);

  double MassDiff1 = abs(JetComboTest1a.M()-JetComboTest1b.M());
  double MassDiff2 = abs(JetComboTest2a.M()-JetComboTest2b.M());
  double MassDiff3 = abs(JetComboTest3a.M()-JetComboTest3b.M());

  std::vector<int> theChosen; //ordered by combo, so 1st+2nd jet, and 3rd+4th jet, which has the smallest mass difference
  if (MassDiff1<MassDiff2 && MassDiff1<MassDiff3) theChosen = {JetIndices[0],JetIndices[1],JetIndices[2],JetIndices[3]};
  else if (MassDiff2<MassDiff1 && MassDiff2<MassDiff3) theChosen = {JetIndices[0],JetIndices[2],JetIndices[1],JetIndices[3]};
  else if (MassDiff3<MassDiff1 && MassDiff3<MassDiff2) theChosen = {JetIndices[0],JetIndices[3],JetIndices[1],JetIndices[2]};
  return theChosen;
}

template<typename ntupleType> float resMassDiff(ntupleType* ntuple) {
  vector<int> theseJets = resHCandidates(ntuple);
  if (theseJets.at(0) == -11) return -999.;
  TLorentzVector Higgs1 = ntuple->Jets->at(theseJets.at(0))+ntuple->Jets->at(theseJets.at(1));
  TLorentzVector Higgs2 = ntuple->Jets->at(theseJets.at(2))+ntuple->Jets->at(theseJets.at(3));
  float massDiff = fabs(Higgs1.M()-Higgs2.M());
  return massDiff;
}

template<typename ntupleType> float resAvgMass(ntupleType* ntuple) {
  vector<int> theseJets = resHCandidates(ntuple);
  if (theseJets.at(0) == -11) return -999.;
  TLorentzVector Higgs1 = ntuple->Jets->at(theseJets.at(0))+ntuple->Jets->at(theseJets.at(1));
  TLorentzVector Higgs2 = ntuple->Jets->at(theseJets.at(2))+ntuple->Jets->at(theseJets.at(3));
  float avgMass = (Higgs1.M()+Higgs2.M())/2 ;
  return avgMass;
}

template<typename ntupleType> float resDeltaRMax(ntupleType* ntuple) {
  vector<int> theseJets = resHCandidates(ntuple);
  if (theseJets.at(0) == -11) return -999.;
  float deltaEta1 = ntuple->Jets->at(theseJets.at(0)).Eta() - ntuple->Jets->at(theseJets.at(1)).Eta();
  float deltaPhi1 = CalcdPhi(ntuple->Jets->at(theseJets.at(0)).Phi() , ntuple->Jets->at(theseJets.at(1)).Phi());
  float deltaEta2 = ntuple->Jets->at(theseJets.at(2)).Eta() - ntuple->Jets->at(theseJets.at(3)).Eta();
  float deltaPhi2 = CalcdPhi(ntuple->Jets->at(theseJets.at(2)).Phi() , ntuple->Jets->at(theseJets.at(3)).Phi());

  float deltaR1 = sqrt((deltaEta1*deltaEta1)+(deltaPhi1*deltaPhi1));
  float deltaR2 = sqrt((deltaEta2*deltaEta2)+(deltaPhi2*deltaPhi2));
  float deltaR_max = max(deltaR1, deltaR2);
  return deltaR_max;
}

template<typename ntupleType> double fillDeltaRMax(ntupleType* ntuple) {
  return resDeltaRMax(ntuple);
}

////////////////////////////////
// HIGHEST PT JET PROPERTIES  //
////////////////////////////////
template<typename ntupleType> double fillLeadingJetMinDRB(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size() == 0) return -999.;
  else {
    double minDRB = 999.;
    double DRB = 999.;
    for (int i = 0 ; i < ntuple->Jets->size(); i++) {
      DRB = 999.;
      if (abs(ntuple->Jets_partonFlavor->at(i)) == 5) {
        DRB = ntuple->JetsAK8->at(0).DeltaR(ntuple->Jets->at(i));
      }
      if (DRB < minDRB) minDRB = DRB;
    }
    return minDRB;
  }
}

template<typename ntupleType> double fillLeadingJetMass(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return -99999.;
  return ntuple->JetsAK8_softDropMass->at(0);
}

template<typename ntupleType> double fillLeadingJetFlavor(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return -99999.;
  if (ntuple->JetsAK8_NumBhadrons->at(0)==2) return 21.;
  else if (ntuple->JetsAK8_NumBhadrons->at(0)==1) return 5.;
  else return 1.;
}


template<typename ntupleType> double fillLeadingNbHadrons(ntupleType* ntuple) {
  return ntuple->JetsAK8->size()>=1?ntuple->JetsAK8_NumBhadrons->at(0):-999.;
}

template<typename ntupleType> double fillLeadingJetPt(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  return ntuple->JetsAK8->at(0).Pt();
}

template<typename ntupleType> double fillLeadingBBtag(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  return ntuple->JetsAK8_doubleBDiscriminator->at(0);
}

template<typename ntupleType> double fillLeadingdeepBBtag(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  // return ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(0);
  return ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(0); //change for V18

}

template<typename ntupleType> double fillTau32(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  return ntuple->JetsAK8_NsubjettinessTau3->at(0)/ntuple->JetsAK8_NsubjettinessTau2->at(0);
}

template<typename ntupleType> double fillLeadingTau21(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  return ntuple->JetsAK8_NsubjettinessTau2->at(0)/ntuple->JetsAK8_NsubjettinessTau1->at(0);
}

//////////////////////////////////////////
// SECOND HIGHEST PT AK8 JET PROPERTIES //
//////////////////////////////////////////
template<typename ntupleType> double fillSubLeadingJetMass(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<=1) return-99999.;
  return ntuple->JetsAK8_softDropMass->at(1);
}

template<typename ntupleType> double fillSubLeadingJetFlavor(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<=1) return-99999.;
  if (ntuple->JetsAK8_NumBhadrons->at(1)==2) return 21.;
  else if (ntuple->JetsAK8_NumBhadrons->at(1)==1) return 5.;
  else return 1.;
}

template<typename ntupleType> double fillSubLeadingNbHadrons(ntupleType* ntuple) {
  return ntuple->JetsAK8->size()>=2?ntuple->JetsAK8_NumBhadrons->at(1):-999.;
}

template<typename ntupleType> double fillSubLeadingJetPt(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<=1) return -99999.;
  return ntuple->JetsAK8->at(1).Pt();
}

template<typename ntupleType> double fillSubLeadingBBtag(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<=1) return-99999.;
  return ntuple->JetsAK8_doubleBDiscriminator->at(1);
}

template<typename ntupleType> double fillSubLeadingdeepBBtag(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<=1) return-99999.;
  // return ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(1);
  return ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(1); //change for V18
}

template<typename ntupleType> double fillSubLeadingTau21(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<=1) return-99999.;
  return ntuple->JetsAK8_NsubjettinessTau2->at(1)/ntuple->JetsAK8_NsubjettinessTau1->at(1);
}

////////////////////////////////////////
// HIGHEST BBtag AK8 JET PROPERTIES  ///
////////////////////////////////////////
template<typename ntupleType> double fillLeadingBBtagJetMass(ntupleType* ntuple) { //change for V18
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_doubleBDiscriminator->at(0) < ntuple->JetsAK8_doubleBDiscriminator->at(1));
    return ntuple->JetsAK8_softDropMass->at(index);
  }
}

template<typename ntupleType> double fillLeadingBBtagJetFlavor(ntupleType* ntuple) {  //change for V18
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_doubleBDiscriminator->at(0) < ntuple->JetsAK8_doubleBDiscriminator->at(1));
    if (ntuple->JetsAK8_NumBhadrons->at(index)==2) return 21.;
    else if (ntuple->JetsAK8_NumBhadrons->at(index)==1) return 5.;
    else return 1.;
  }
}

template<typename ntupleType> double fillLeadingBBtagJetPt(ntupleType* ntuple) {  //change for V18
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_doubleBDiscriminator->at(0) < ntuple->JetsAK8_doubleBDiscriminator->at(1));
    return ntuple->JetsAK8->at(index).Pt();
  }
}

template<typename ntupleType> double fillLeadingBBtagJetBBtag(ntupleType* ntuple) {  //change for V18
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_doubleBDiscriminator->at(0) < ntuple->JetsAK8_doubleBDiscriminator->at(1));
    return ntuple->JetsAK8_doubleBDiscriminator->at(index);
  }
}

template<typename ntupleType> double fillLeadingBBtagJetTau21(ntupleType* ntuple) {  //change for V18
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_doubleBDiscriminator->at(0) < ntuple->JetsAK8_doubleBDiscriminator->at(1));
    return ntuple->JetsAK8_NsubjettinessTau2->at(index)/ntuple->JetsAK8_NsubjettinessTau1->at(index);
  }
}

/////////////////////////////////////////////
// SECOND HIGHEST BBtag AK8 JET PROPERTIES //
/////////////////////////////////////////////
template<typename ntupleType> double fillSubLeadingBBtagJetMass(ntupleType* ntuple) {  //change for V18
  int index = int(ntuple->JetsAK8_doubleBDiscriminator->at(0) > ntuple->JetsAK8_doubleBDiscriminator->at(1));
  return ntuple->JetsAK8_softDropMass->at(index);
}

template<typename ntupleType> double fillSubLeadingBBtagJetFlavor(ntupleType* ntuple) {  //change for V18
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_doubleBDiscriminator->at(0) > ntuple->JetsAK8_doubleBDiscriminator->at(1));
    if (ntuple->JetsAK8_NumBhadrons->at(index)==2) return 21.;
    else if (ntuple->JetsAK8_NumBhadrons->at(index)==1) return 5.;
    else return 1.;
  }
}

template<typename ntupleType> double fillSubLeadingBBtagJetPt(ntupleType* ntuple) {  //change for V18
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_doubleBDiscriminator->at(0) > ntuple->JetsAK8_doubleBDiscriminator->at(1));
    return ntuple->JetsAK8->at(index).Pt();
  }
}

template<typename ntupleType> double fillSubLeadingBBtagJetBBtag(ntupleType* ntuple) {  //change for V18
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_doubleBDiscriminator->at(0) > ntuple->JetsAK8_doubleBDiscriminator->at(1));
    return ntuple->JetsAK8_doubleBDiscriminator->at(index);
  }
}

template<typename ntupleType> double fillSubLeadingBBtagJetTau21(ntupleType* ntuple) {  //change for V18
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_doubleBDiscriminator->at(0) > ntuple->JetsAK8_doubleBDiscriminator->at(1));
    return ntuple->JetsAK8_NsubjettinessTau2->at(index)/ntuple->JetsAK8_NsubjettinessTau1->at(index);
  }
}

///////////////////////////////////////
// HIGHEST BBtag AK8 JET PROPERTIES ///
///////////////////////////////////////
template<typename ntupleType> double fillLeadingMassJetMass(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_softDropMass->at(0) < ntuple->JetsAK8_softDropMass->at(1));
    return ntuple->JetsAK8_softDropMass->at(index);
  }
}

template<typename ntupleType> double fillLeadingMassJetFlavor(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_softDropMass->at(0) < ntuple->JetsAK8_softDropMass->at(1));
    if (ntuple->JetsAK8_NumBhadrons->at(index)==2) return 21.;
    else if (ntuple->JetsAK8_NumBhadrons->at(index)==1) return 5.;
    else return 1.;
  }
}

template<typename ntupleType> double fillLeadingMassJetPt(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_softDropMass->at(0) < ntuple->JetsAK8_softDropMass->at(1));
    return ntuple->JetsAK8->at(index).Pt();
  }
}
template<typename ntupleType> double fillLeadingMassJetBBtag(ntupleType* ntuple) { //change for V18
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_softDropMass->at(0) < ntuple->JetsAK8_softDropMass->at(1));
    return ntuple->JetsAK8_doubleBDiscriminator->at(index);
  }
}

template<typename ntupleType> double fillLeadingMassJetTau21(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_softDropMass->at(0) < ntuple->JetsAK8_softDropMass->at(1));
    return ntuple->JetsAK8_NsubjettinessTau2->at(index)/ntuple->JetsAK8_NsubjettinessTau1->at(index);
  }
}

////////////////////////////////////////////
// SECOND HIGHEST Mass AK8 JET PROPERTIES //
////////////////////////////////////////////
template<typename ntupleType> double fillSubLeadingMassJetMass(ntupleType* ntuple) {
  int index = int(ntuple->JetsAK8_softDropMass->at(0) > ntuple->JetsAK8_softDropMass->at(1));
  return ntuple->JetsAK8_softDropMass->at(index);
}

template<typename ntupleType> double fillSubLeadingMassJetFlavor(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_softDropMass->at(0) > ntuple->JetsAK8_softDropMass->at(1));
    if (ntuple->JetsAK8_NumBhadrons->at(index)==2) return 21.;
    else if (ntuple->JetsAK8_NumBhadrons->at(index)==1) return 5.;
    else return 1.;
  }
}

template<typename ntupleType> double fillSubLeadingMassJetPt(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_softDropMass->at(0) > ntuple->JetsAK8_softDropMass->at(1));
    return ntuple->JetsAK8->at(index).Pt();
  }
}

template<typename ntupleType> double fillSubLeadingMassJetBBtag(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_softDropMass->at(0) > ntuple->JetsAK8_softDropMass->at(1));
    return ntuple->JetsAK8_doubleBDiscriminator->at(index);
  }
}

template<typename ntupleType> double fillSubLeadingMassJetTau21(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()==0) return-99999.;
  int index;
  if (ntuple->JetsAK8->size()==1) index = 0;
  else {
    index = int(ntuple->JetsAK8_softDropMass->at(0) > ntuple->JetsAK8_softDropMass->at(1));
    return ntuple->JetsAK8_NsubjettinessTau2->at(index)/ntuple->JetsAK8_NsubjettinessTau1->at(index);
  }
}

////////////////////////////////////////////
// SECOND HIGHEST Mass AK8 JET PROPERTIES //
////////////////////////////////////////////
template<typename ntupleType> double fillClosestJetMass(ntupleType* ntuple) {
  if (ntuple->JetsAK8_softDropMass->size() == 0)
  return -999.;
  else if (ntuple->JetsAK8_softDropMass->size() == 1)
  return ntuple->JetsAK8_softDropMass->at(0);
  else {
    double J1diff,J2diff;
    J1diff = ntuple->JetsAK8_softDropMass->at(0)-110.;
    J2diff = ntuple->JetsAK8_softDropMass->at(1)-110.;
    return fabs(J1diff)>fabs(J2diff) ?  ntuple->JetsAK8_softDropMass->at(1) :  ntuple->JetsAK8_softDropMass->at(0);
  }
}

template<typename ntupleType> double fillFarthestJetMass(ntupleType* ntuple) {
  if (ntuple->JetsAK8_softDropMass->size() == 0)
  return -999.;
  else if (ntuple->JetsAK8_softDropMass->size() == 1)
  return ntuple->JetsAK8_softDropMass->at(0);
  else {
    double J1diff,J2diff;
    J1diff = ntuple->JetsAK8_softDropMass->at(0)-110.;
    J2diff = ntuple->JetsAK8_softDropMass->at(1)-110.;
    return fabs(J1diff)<fabs(J2diff) ?  ntuple->JetsAK8_softDropMass->at(1) :  ntuple->JetsAK8_softDropMass->at(0);
  }
}

/////////////////
// OTHER STUFF //
/////////////////
template<typename ntupleType> double fillAnalysisBins(ntupleType* ntuple) {
  double MET = ntuple->MET;
  double HT = ntuple->HT;

  if (MET > 300. && MET < 600.) {
    if (HT > 300. && HT < 1000.) {return 1.;}
    else if (HT > 1000. && HT < 2000.) {return 2.;}
    else if (HT > 2000.) {return 3.;}
    else return -1.;
  }
  else if (MET > 600. && MET < 1000.) {
    if (HT > 600. && HT < 1000.) {return 4.;}
    else if (HT > 1000. && HT < 2000.) {return 5.;}
    else if (HT > 2000.) {return 6.;}
    else return -1.;
  }
  else if (MET > 1000.) {
    if (HT > 1000. && HT < 2000.) {return 7.;}
    else if (HT > 2000.) {return 8.;}
    else return -1.;
  }
  else return -1.;
}

template<typename ntupleType> double fillRA2b10Bins(ntupleType* ntuple) {
  double MET = ntuple->met_pt;
  double HT = fillHT(ntuple);

  if (MET > 300. && MET < 350.) {
    if (HT > 300. && HT < 500.) {return 1.;}
    else if (HT > 500. && HT < 1000.) {return 2.;}
    else if (HT > 1000.) {return 3.;}
    else return -1.;
  }
  else if (MET > 350. && MET < 500.) {
    if (HT > 350. && HT < 500.) {return 4.;}
    else if (HT > 500. && HT < 1000.) {return 5.;}
    else if (HT > 1000.) {return 6.;}
    else return -1.;
  }
  else if (MET > 500. && MET < 750.) {
    if (HT > 500. && HT < 1000.) {return 7.;}
    else if (HT > 1000.) {return 8.;}
    else return -1.;
  }
  else if (MET > 750.) {
    if (HT > 750. && HT < 1500.) {return 9.;}
    else if (HT > 1500.) {return 10.;}
    else return -1.;
  }
  else return -1.;
}

template<typename ntupleType> double fillRA2b160Bins(ntupleType* ntuple) {
  int BTags = int(ntuple->BTags);
  int NJets = int(ntuple->NJets);

  if (NJets >= 3 && NJets <=4) {
    if (BTags == 0) return fillRA2b10Bins(ntuple);
    else if (BTags == 1) return 10.+fillRA2b10Bins(ntuple);
    else if (BTags == 2) return 20.+fillRA2b10Bins(ntuple);
    else if (BTags >= 3) return 30.+fillRA2b10Bins(ntuple);
  }
  else if (NJets >= 5 && NJets <= 6) {
    if (BTags == 0) return 40.+fillRA2b10Bins(ntuple);
    else if (BTags == 1) return 50.+fillRA2b10Bins(ntuple);
    else if (BTags == 2) return 60.+fillRA2b10Bins(ntuple);
    else if (BTags >= 3) return 70.+fillRA2b10Bins(ntuple);
  }
  else if (NJets >= 7 && NJets <= 8) {
    if (BTags == 0) return 80.+fillRA2b10Bins(ntuple);
    else if (BTags == 1) return 90.+fillRA2b10Bins(ntuple);
    else if (BTags == 2) return 100.+fillRA2b10Bins(ntuple);
    else if (BTags >= 3) return 110.+fillRA2b10Bins(ntuple);
  }
  else if (NJets >= 9) {
    if (BTags == 0) return 120.+fillRA2b10Bins(ntuple);
    else if (BTags == 1) return 130.+fillRA2b10Bins(ntuple);
    else if (BTags == 2) return 140.+fillRA2b10Bins(ntuple);
    else if (BTags >= 3) return 150.+fillRA2b10Bins(ntuple);
  }
  else return -1.;
}

template<typename ntupleType> bool ptBinCut(double pt , int ithBin) {
  if (ithBin > 5) return false;
  double ptCut[6] = {300.,400.,500.,700.,1000.,999999.};
  return pt>ptCut[ithBin] && pt<ptCut[ithBin+1];
}

template<typename ntupleType> bool RA2bBaselineCut(ntupleType* ntuple) {
  double DeltaPhi1 = ntuple->DeltaPhi1;
  double DeltaPhi2 = ntuple->DeltaPhi2;
  double DeltaPhi3 = ntuple->DeltaPhi3;
  double DeltaPhi4 = ntuple->DeltaPhi4;
  double HT = ntuple->HT;
  double MET = ntuple->MET;
  int NJets = ntuple->NJets;

  return (NJets == 3 && MET > 300. && HT > 300. && DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3) || (NJets > 3 && MET > 300. && HT > 300. && DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3);
}

template<typename ntupleType> bool FiltersCut(ntupleType* ntuple) {
  TString filename = ntuple->fChain->GetFile()->GetName();
  if (filename.Contains("TChiHH_HToBB")) {
    return (
      ntuple->NVtx>0 &&
      ntuple->MET/ntuple->CaloMET < 5. &&
      ntuple->LowNeutralJetFilter==1 &&
      ntuple->HTRatioDPhiTightFilter==1 &&
      ntuple->FakeJetFilter==1
   );
  }
  else {
    return (
      ntuple->HBHENoiseFilter==1 &&
      ntuple->eeBadScFilter==1 &&
      ntuple->NVtx>0 &&
      ntuple->MET/ntuple->CaloMET < 5. &&
      ntuple->BadPFMuonFilter == 1 &&
      ntuple->globalSuperTightHalo2016Filter==1 &&
      ntuple->LowNeutralJetFilter==1 &&
      ntuple->HTRatioDPhiTightFilter==1 &&
      ntuple->FakeJetFilter==1
   );
  }
}

template<typename ntupleType> bool AK8MultCut(ntupleType* ntuple) {
  return ntuple->JetsAK8->size()>1 ;
}

template<typename ntupleType> bool BTagsCut(ntupleType* ntuple) {
  return ntuple->BTags>0 ;
}

template<typename ntupleType> bool DeltaPhi1Cut(ntupleType* ntuple) {
  return ntuple->DeltaPhi1>0.5;
}

template<typename ntupleType> bool DeltaPhi2Cut(ntupleType* ntuple) {
  return ntuple->DeltaPhi2>0.5;
}

template<typename ntupleType> bool DeltaPhi3Cut(ntupleType* ntuple) {
  return ntuple->DeltaPhi3>0.3;
}

template<typename ntupleType> bool DeltaPhi4Cut(ntupleType* ntuple) {
  return ntuple->DeltaPhi4>0.3;
}

template<typename ntupleType> bool DeltaPhiCuts(ntupleType* ntuple) {
  return (
    DeltaPhi1Cut(ntuple) &&
    DeltaPhi2Cut(ntuple) &&
    DeltaPhi3Cut(ntuple) &&
    DeltaPhi4Cut(ntuple)
 );
}

template<typename ntupleType> bool lowDPhiCuts(ntupleType* ntuple) {
  return !DeltaPhiCuts(ntuple);
}

template<typename ntupleType> bool METHTlooseCut(ntupleType* ntuple) {
  return (ntuple->MET > 100. && ntuple->HT > 300.);
}

template<typename ntupleType> bool METHTCut(ntupleType* ntuple) {
  return (ntuple->MET > 300. && ntuple->HT > 600.);
}

template<typename ntupleType> bool AK8JetPtCut(ntupleType* ntuple) {
  return (
    ntuple->JetsAK8->size() >= 2 &&
    ntuple->JetsAK8->at(0).Pt() > 300. &&
    ntuple->JetsAK8->at(1).Pt() > 300.
 );
}

template<typename ntupleType> bool AK8JetLooseMassCut(ntupleType* ntuple) {
  return (
    ntuple->JetsAK8_softDropMass->at(0) > 50. &&
    ntuple->JetsAK8_softDropMass->at(0) < 250. &&
    ntuple->JetsAK8_softDropMass->at(1) > 50. &&
    ntuple->JetsAK8_softDropMass->at(1) < 250.
 );
}


template<typename ntupleType> bool baselineCut(ntupleType* ntuple) {
  int numAK4Jets_cuts = numJets(ntuple);
  int numTightBs = numDeepBs(ntuple).at(2);
  bool isMaybeRes = (numTightBs>=2 && numAK4Jets_cuts>=4 && numAK4Jets_cuts<6);
  bool isMaybeBoost = false;
  if (ntuple->JetsAK8->size()>=2) isMaybeBoost = (ntuple->JetsAK8->at(0).Pt()>300.0 && ntuple->JetsAK8->at(1).Pt()>300.0);

  return (
    ntuple->MET > 150. &&
    (isMaybeRes || isMaybeBoost) &&
    DeltaPhiCuts(ntuple) &&
    ntuple->NMuons==0 && ntuple->NElectrons==0 &&
    ntuple->isoElectronTracks==0 && ntuple->isoMuonTracks==0 && ntuple->isoPionTracks==0 &&
    // FiltersCut(ntuple) &&
    ntuple->JetID == 1
 );
}

template<typename ntupleType> bool cutflowLeptVeto(ntupleType* ntuple) {
  return (ntuple->NMuons==0 && ntuple->NElectrons==0);
}

template<typename ntupleType> bool cutflowDPhiCut(ntupleType* ntuple) {
  return (DeltaPhiCuts(ntuple) && ntuple->isoElectronTracks+ntuple->isoMuonTracks +ntuple->isoPionTracks==0);
}

template<typename ntupleType> bool cutflowBoostBase(ntupleType* ntuple) {
  return (
    ntuple->MET > 300. && ntuple->HT > 300. &&
    // FiltersCut(ntuple) &&
    DeltaPhiCuts(ntuple) &&
    ntuple->NMuons==0 && ntuple->NElectrons==0 &&
    ntuple->isoElectronTracks+ntuple->isoMuonTracks +ntuple->isoPionTracks==0
 );
}

template<typename ntupleType> bool cutflowBoostHT(ntupleType* ntuple) {
  return (ntuple->HT > 600.);
}

template<typename ntupleType> bool cutflowBoost2AK8(ntupleType* ntuple) {
  return (ntuple->JetsAK8->size()>1);
}

template<typename ntupleType> bool cutflowBoostPt(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<2) return false;
  else return (ntuple->JetsAK8->at(0).Pt() > 300. && ntuple->JetsAK8->at(1).Pt() > 300.);
}

template<typename ntupleType> bool cutflowBoostLooseMass(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<2) return false;
  else return (
    ntuple->JetsAK8_softDropMass->at(0) > 50. &&
    ntuple->JetsAK8_softDropMass->at(0) < 250. &&
    ntuple->JetsAK8_softDropMass->at(1) > 50. &&
    ntuple->JetsAK8_softDropMass->at(1) < 250.
 );
}

template<typename ntupleType> bool cutflowBoostBBTag(ntupleType* ntuple) { //change for V18
  TString filename = ntuple->fChain->GetFile()->GetName();
  if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH")) { //V18 signal has newest deep bbtag, so use that in cuts
    float this_bbtag = 0.70;
    if (ntuple->JetsAK8->size()<2 || ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->size()<2) return false;
    else return (ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(0) > this_bbtag && ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(1) > this_bbtag);
  }
  else { //for background, use unsupported version of deep bbtag
    float this_bbtag = 0.70;
    if (ntuple->JetsAK8->size()<2 || ntuple->JetsAK8_deepDoubleBDiscriminatorH->size()<2) return false;
    else return (ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(0) > this_bbtag && ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(1) > this_bbtag);
  }
}

template<typename ntupleType> bool cutflowBoostTightMass(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<2) return false;
  else return (
    ntuple->JetsAK8_softDropMass->at(0) > 85. &&
    ntuple->JetsAK8_softDropMass->at(0) < 135. &&
    ntuple->JetsAK8_softDropMass->at(1) > 85. &&
    ntuple->JetsAK8_softDropMass->at(1) < 135.
 );
}

template<typename ntupleType> bool boostedBaselineCut(ntupleType* ntuple) {
  return (
    baselineCut(ntuple) &&
    ntuple->MET > 300. &&
    ntuple->HT > 600. &&
    ntuple->JetsAK8->size() >= 2 &&
    ntuple->JetsAK8->at(0).Pt() > 300. &&
    ntuple->JetsAK8_softDropMass->at(0) > 50. &&
    ntuple->JetsAK8_softDropMass->at(0) < 250. &&
    ntuple->JetsAK8->at(1).Pt() > 300. &&
    ntuple->JetsAK8_softDropMass->at(1) > 50. &&
    ntuple->JetsAK8_softDropMass->at(1) < 250.
 );
}

template<typename ntupleType> bool resolvedBaselineCut(ntupleType* ntuple) {
  float avgMass = resAvgMass(ntuple);
  float massDiff = resMassDiff(ntuple);
  float deltaR = resDeltaRMax(ntuple);

  return (
    baselineCut(ntuple) &&
    avgMass>=0.0 && avgMass<=250.0 &&
    massDiff>=0.0 && massDiff<=40.0 &&
    deltaR>=0.0 && deltaR<=2.2
 );
}

template<typename ntupleType> bool singleMuCut(ntupleType* ntuple) {
  if (ntuple->Muons->size() != 1 || ntuple->Electrons->size() != 0) return false;
  double MT = computeMuonMT(ntuple);
  return (ntuple->Muons->at(0).Pt()>25. && MT < 100.);
}

template<typename ntupleType> bool singleMuBaselineCut(ntupleType* ntuple) {

  return (
    singleMuCut(ntuple) &&
    ntuple->MET > 100. &&
    ntuple->HT > 600. &&
    ntuple->JetsAK8->size() >= 2 &&
    ntuple->JetsAK8->at(0).Pt() > 300. &&
    ntuple->JetsAK8_softDropMass->at(0) > 50. &&
    ntuple->JetsAK8_softDropMass->at(0) < 250. &&
    ntuple->JetsAK8->at(1).Pt() > 300. &&
    ntuple->JetsAK8_softDropMass->at(1) > 50. &&
    ntuple->JetsAK8_softDropMass->at(1) < 250.&&
    DeltaPhiCuts(ntuple) &&
    FiltersCut(ntuple) &&
    ntuple->JetID == 1
 );
}

template<typename ntupleType> bool singleEleCut(ntupleType* ntuple) {
  if (ntuple->Muons->size() != 0 || ntuple->Electrons->size() != 1) return false;
  double MT = computeElectronMT(ntuple);
  return (ntuple->Electrons->at(0).Pt()>25. && MT < 100.);
}

template<typename ntupleType> bool singleEleBaselineCut(ntupleType* ntuple) {
  return (
    singleEleCut(ntuple) &&
    ntuple->MET > 100. &&
    ntuple->HT > 600. &&
    ntuple->JetsAK8->size() >= 2 &&
    ntuple->JetsAK8->at(0).Pt() > 300. &&
    ntuple->JetsAK8_softDropMass->at(0) > 50. &&
    ntuple->JetsAK8_softDropMass->at(0) < 250. &&
    ntuple->JetsAK8->at(1).Pt() > 300. &&
    ntuple->JetsAK8_softDropMass->at(1) > 50. &&
    ntuple->JetsAK8_softDropMass->at(1) < 250.&&
    DeltaPhiCuts(ntuple) &&
    FiltersCut(ntuple) &&
    ntuple->JetID == 1
 );
}

template<typename ntupleType> bool lowDphiBaselineCut(ntupleType* ntuple) {
  return (
    ntuple->MET > 300. &&
    ntuple->HT > 600. &&
    ntuple->JetsAK8->size() >= 2 &&
    ntuple->JetsAK8->at(0).Pt() > 300. &&
    ntuple->JetsAK8_softDropMass->at(0) > 50. &&
    ntuple->JetsAK8_softDropMass->at(0) < 250. &&
    ntuple->JetsAK8->at(1).Pt() > 300. &&
    ntuple->JetsAK8_softDropMass->at(1) > 50. &&
    ntuple->JetsAK8_softDropMass->at(1) < 250.&&
    ! DeltaPhiCuts(ntuple) &&
    FiltersCut(ntuple) &&
    ntuple->JetID == 1
 );
}

template<typename ntupleType> bool photonBaselineCut(ntupleType* ntuple) {
  return (
    ntuple->Photons->size()==1 &&
    ntuple->Photons->at(0).Pt() > 100. &&
    ntuple->Photons_fullID->size() == 1 &&
    ntuple->Photons_fullID->at(0) == 1 &&
    ntuple->MET > 100. &&
    ntuple->HT > 400. &&
    ntuple->JetsAK8->size()>=2 &&
    ntuple->JetsAK8->at(0).Pt() > 300. &&
    ntuple->JetsAK8_softDropMass->at(0) > 50. &&
    ntuple->JetsAK8_softDropMass->at(0) < 250. &&
    ntuple->JetsAK8->at(1).Pt() > 300. &&
    ntuple->JetsAK8_softDropMass->at(1) > 50. &&
    ntuple->JetsAK8_softDropMass->at(1) < 250.&&
    ntuple->DeltaPhi1>0.5 &&
    ntuple->DeltaPhi2>0.5 &&
    ntuple->DeltaPhi3>0.3 &&
    ntuple->DeltaPhi4>0.3 &&
    ntuple->isoElectronTracks==0 &&
    ntuple->isoMuonTracks == 0 &&
    ntuple->isoPionTracks == 0 &&
    ntuple->Electrons->size() == 0 &&
    ntuple->Muons->size() == 0 &&
    FiltersCut(ntuple) &&
    ntuple->JetID == 1
 );
}

template<typename ntupleType> bool photonBaselineCut_loose(ntupleType* ntuple) {
  return (
    ntuple->Photons->size()==1 &&
    ntuple->Photons->at(0).Pt() > 100. &&
    ntuple->Photons_fullID->size() == 1 &&
    ntuple->Photons_fullID->at(0) == 1 &&
    ntuple->MET > 100. &&
    ntuple->HT > 400. &&
    ntuple->JetsAK8->size()>=2 &&
    ntuple->DeltaPhi1>0.5 &&
    ntuple->DeltaPhi2>0.5 &&
    ntuple->DeltaPhi3>0.3 &&
    ntuple->DeltaPhi4>0.3 &&
    ntuple->isoElectronTracks==0 &&
    ntuple->isoMuonTracks == 0 &&
    ntuple->isoPionTracks == 0 &&
    ntuple->Electrons->size() == 0 &&
    ntuple->Muons->size() == 0 &&
    FiltersCut(ntuple) &&
    ntuple->JetID == 1
 );
}

template<typename ntupleType> bool singleHiggsTagLooseCut(ntupleType* ntuple) { //change for V18
  float this_bbtag = 0.70;
  return ((ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(0) > this_bbtag)
  && (ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(1) < this_bbtag)) ||
  ((ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(0) < this_bbtag)
  && (ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(1) > this_bbtag));
}

template<typename ntupleType> bool antiTaggingLooseCut(ntupleType* ntuple) { //change for V18
  float this_bbtag = 0.70;
  return (
    (ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(0) < this_bbtag) &&
    (ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(1) < this_bbtag)
 );
}

template<typename ntupleType> bool doubleTaggingLooseCut(ntupleType* ntuple) { //change for V18
  float this_bbtag = 0.70;
  return (
    ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(0) > this_bbtag &&
    ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(1) > this_bbtag
 );
}

template<typename ntupleType> bool doubleMassCut(ntupleType* ntuple) {
  return (
    ntuple->JetsAK8_softDropMass->at(0) > 85. &&
    ntuple->JetsAK8_softDropMass->at(0) < 135. &&
    ntuple->JetsAK8_softDropMass->at(1) > 85. &&
    ntuple->JetsAK8_softDropMass->at(1) < 135.
 );
}

template<typename ntupleType> bool singleHiggsTagCut(ntupleType* ntuple) {  //change for V18
  float this_bbtag = 0.70;
  return (
    (ntuple->JetsAK8_softDropMass->at(0) > 85. &&
    ntuple->JetsAK8_softDropMass->at(0) < 135. &&
    ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(0) > this_bbtag) ||
    (ntuple->JetsAK8_softDropMass->at(1) > 85. &&
    ntuple->JetsAK8_softDropMass->at(1) < 135. &&
    ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(1) > this_bbtag)
 );
}

template<typename ntupleType> bool doubleHiggsTagCut(ntupleType* ntuple) {  //change for V18
  float this_bbtag = 0.70;
  return (
    ntuple->JetsAK8_softDropMass->at(0) > 85. &&
    ntuple->JetsAK8_softDropMass->at(0) < 135. &&
    ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(0) > this_bbtag &&
    ntuple->JetsAK8_softDropMass->at(1) > 85. &&
    ntuple->JetsAK8_softDropMass->at(1) < 135. &&
    ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(1) > this_bbtag
 );
}

template<typename ntupleType> bool tagSR(ntupleType* ntuple, int i) {  //change for V18
  if (ntuple->JetsAK8->size()<2) return false;
  TString filename = ntuple->fChain->GetFile()->GetName();
  if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH")) { //Use supported deep bbtag
    if (ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->size() <= i ||
    ntuple->JetsAK8_softDropMass->size() <= i) return false;
    float this_bbtag = 0.70;
    return (
      ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(i) > this_bbtag &&
      ntuple->JetsAK8_softDropMass->at(i) > 85. &&
      ntuple->JetsAK8_softDropMass->at(i) < 135.
   );
  }
  else { //use unsupported
    if (ntuple->JetsAK8_deepDoubleBDiscriminatorH->size() <= i ||
    ntuple->JetsAK8_softDropMass->size() <= i) return false;
    float this_bbtag = 0.70;
    return (
      ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(i) > this_bbtag &&
      ntuple->JetsAK8_softDropMass->at(i) > 85. &&
      ntuple->JetsAK8_softDropMass->at(i) < 135.
   );
  }
}

template<typename ntupleType> bool tagSB(ntupleType* ntuple, int i) {  //change for V18
  if (ntuple->JetsAK8->size()<2) return false;
  TString filename = ntuple->fChain->GetFile()->GetName();
  if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH")) { //Use supported deep bbtag
    if (ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->size() <= i ||
    ntuple->JetsAK8_softDropMass->size() <= i) return false;
    float this_bbtag = 0.70;
    return (
      ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(i) > this_bbtag &&
      ((ntuple->JetsAK8_softDropMass->at(i) < 85. &&
      ntuple->JetsAK8_softDropMass->at(i) > 50.) ||
      (ntuple->JetsAK8_softDropMass->at(i) > 135. &&
      ntuple->JetsAK8_softDropMass->at(i) < 250.))
   );
  }
  else { //use unsupported bbtag
    if (ntuple->JetsAK8_deepDoubleBDiscriminatorH->size() <= i ||
    ntuple->JetsAK8_softDropMass->size() <= i) return false;
    float this_bbtag = 0.70;
    return (
      ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(i) > this_bbtag &&
      ((ntuple->JetsAK8_softDropMass->at(i) < 85. &&
      ntuple->JetsAK8_softDropMass->at(i) > 50.) ||
      (ntuple->JetsAK8_softDropMass->at(i) > 135. &&
      ntuple->JetsAK8_softDropMass->at(i) < 250.))
   );
  }
}

template<typename ntupleType> bool antitagSR(ntupleType* ntuple, int i) {  //change for V18
  if (ntuple->JetsAK8->size()<2) return false;
  TString filename = ntuple->fChain->GetFile()->GetName();
  if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH")) { //Use supported deep bbtag
    if (ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->size() <= i ||
    ntuple->JetsAK8_softDropMass->size() <= i) return false;
    float this_bbtag = 0.70;
    return (
      ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(i) < this_bbtag &&
      (ntuple->JetsAK8_softDropMass->at(i) > 85. &&
      ntuple->JetsAK8_softDropMass->at(i) < 135.)
   );
  }
  else { //use unsupported deep bbtag
    if (ntuple->JetsAK8_deepDoubleBDiscriminatorH->size() <= i ||
    ntuple->JetsAK8_softDropMass->size() <= i) return false;
    float this_bbtag = 0.70;
    return (
      ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(i) < this_bbtag &&
      (ntuple->JetsAK8_softDropMass->at(i) > 85. &&
      ntuple->JetsAK8_softDropMass->at(i) < 135.)
   );
  }
}

template<typename ntupleType> bool antitagSB(ntupleType* ntuple, int i) {  //change for V18
  TString filename = ntuple->fChain->GetFile()->GetName();
  if (filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH")) { //Use supported deep bbtag
    if (ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->size() <= i ||
    ntuple->JetsAK8_softDropMass->size() <= i) return false;
    float this_bbtag = 0.70;
    return (
      ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(i) < this_bbtag &&
      ((ntuple->JetsAK8_softDropMass->at(i) < 85. &&
      ntuple->JetsAK8_softDropMass->at(i) > 50.) ||
      (ntuple->JetsAK8_softDropMass->at(i) > 135. &&
      ntuple->JetsAK8_softDropMass->at(i) < 250.))
   );
  }
  else { //use unsupported deep bbtag
    if (ntuple->JetsAK8_deepDoubleBDiscriminatorH->size() <= i ||
    ntuple->JetsAK8_softDropMass->size() <= i) return false;
    float this_bbtag = 0.70;
    return (
      ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(i) < this_bbtag &&
      ((ntuple->JetsAK8_softDropMass->at(i) < 85. &&
      ntuple->JetsAK8_softDropMass->at(i) > 50.) ||
      (ntuple->JetsAK8_softDropMass->at(i) > 135. &&
      ntuple->JetsAK8_softDropMass->at(i) < 250.))
   );
  }

}

template<typename ntupleType> bool antitagSRCut(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<2) return false;
  return (antitagSR(ntuple,0) && antitagSR(ntuple,1));
}

template<typename ntupleType> bool newantitagSRCut(ntupleType* ntuple) {
  return (newantitagSR(ntuple,0) && newantitagSR(ntuple,1));
}

template<typename ntupleType> bool antitagSBCut(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<2) return false;
  return ((antitagSB(ntuple,0) && antitagSB(ntuple,1)) ||
  (antitagSB(ntuple,0) && antitagSR(ntuple,1)) ||
  (antitagSR(ntuple,0) && antitagSB(ntuple,1)));
}

template<typename ntupleType> bool tagSRCut(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<2) return false;
  return ((tagSR(ntuple,0) && antitagSR(ntuple,1)) ||
  (antitagSR(ntuple,0) && tagSR(ntuple,1)));
}

template<typename ntupleType> bool tagSBCut(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<2) return false;
  return (
    (tagSB(ntuple,0) && antitagSB(ntuple,1)) ||
    (tagSR(ntuple,0) && antitagSB(ntuple,1)) ||
    (tagSB(ntuple,0) && antitagSR(ntuple,1)) ||
    (antitagSB(ntuple,0) && tagSB(ntuple,1)) ||
    (antitagSR(ntuple,0) && tagSB(ntuple,1)) ||
    (antitagSB(ntuple,0) && tagSR(ntuple,1))
 );
}

template<typename ntupleType> bool doubletagSRCut(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<2) return false;
  return (tagSR(ntuple,0) && tagSR(ntuple,1));
}

template<typename ntupleType> bool doubletagSBCut(ntupleType* ntuple) {
  if (ntuple->JetsAK8->size()<2) return false;
  return (
    (tagSB(ntuple,0) && tagSB(ntuple,1)) ||
    (tagSB(ntuple,0) && tagSR(ntuple,1)) ||
    (tagSR(ntuple,0) && tagSB(ntuple,1))
 );
}

//Resolved ABCD region cuts
template<typename ntupleType> bool fourbSRCut(ntupleType* ntuple) {
  float avgMass = resAvgMass(ntuple);
  float massDiff = resMassDiff(ntuple);
  float deltaR = resDeltaRMax(ntuple);
  vector<int> thisNBs = numDeepBs(ntuple);
  int numLoose = thisNBs.at(0);  int numMed = thisNBs.at(1);  int numTight = thisNBs.at(2);

  return (
    baselineCut(ntuple) &&
    (avgMass>=100.0 && avgMass<=140.0) &&
    massDiff>=0.0 && massDiff<=40.0 &&
    deltaR>=0.0 && deltaR<=2.2 &&
    numLoose>=4 && numMed>=3 && numTight>=2
 );
}

template<typename ntupleType> bool fourbSBCut(ntupleType* ntuple) {
  float avgMass = resAvgMass(ntuple);
  float massDiff = resMassDiff(ntuple);
  float deltaR = resDeltaRMax(ntuple);
  vector<int> thisNBs = numDeepBs(ntuple);
  int numLoose = thisNBs.at(0);  int numMed = thisNBs.at(1);  int numTight = thisNBs.at(2);

  return (
    baselineCut(ntuple) &&
    ((avgMass>=0.0 && avgMass<=100.0) || (avgMass>=140.0 && avgMass<=250.0)) &&
    massDiff>=0.0 && massDiff<=40.0 &&
    deltaR>=0.0 && deltaR<=2.2 &&
    numLoose>=4 && numMed>=3 && numTight>=2
 );
}

template<typename ntupleType> bool threebSRCut(ntupleType* ntuple) {
  float avgMass = resAvgMass(ntuple);
  float massDiff = resMassDiff(ntuple);
  float deltaR = resDeltaRMax(ntuple);
  vector<int> thisNBs = numDeepBs(ntuple);
  int numLoose = thisNBs.at(0);  int numMed = thisNBs.at(1);  int numTight = thisNBs.at(2);

  return (
    baselineCut(ntuple) &&
    (avgMass>=100.0 && avgMass<=140.0) &&
    massDiff>=0.0 && massDiff<=40.0 &&
    deltaR>=0.0 && deltaR<=2.2 &&
    numMed==3 && numTight>=2
 );
}

template<typename ntupleType> bool threebSBCut(ntupleType* ntuple) {
  float avgMass = resAvgMass(ntuple);
  float massDiff = resMassDiff(ntuple);
  float deltaR = resDeltaRMax(ntuple);
  vector<int> thisNBs = numDeepBs(ntuple);
  int numLoose = thisNBs.at(0); int numMed = thisNBs.at(1); int numTight = thisNBs.at(2);

  return (
    baselineCut(ntuple) &&
    ((avgMass>=0.0 && avgMass<=100.0) || (avgMass>=140.0 && avgMass<=250.0)) &&
    massDiff>=0.0 && massDiff<=40.0 &&
    deltaR>=0.0 && deltaR<=2.2 &&
    numMed==3 && numTight>=2
 );
}

template<typename ntupleType> bool twobSRCut(ntupleType* ntuple) {
  float avgMass = resAvgMass(ntuple);
  float massDiff = resMassDiff(ntuple);
  float deltaR = resDeltaRMax(ntuple);
  vector<int> thisNBs = numDeepBs(ntuple);
  int numLoose = thisNBs.at(0); int numMed = thisNBs.at(1); int numTight = thisNBs.at(2);

  return (
    baselineCut(ntuple) &&
    (avgMass>=100.0 && avgMass<=140.0) &&
    massDiff>=0.0 && massDiff<=40.0 &&
    deltaR>=0.0 && deltaR<=2.2 &&
    numMed==2 && numTight>=2
 );
}

template<typename ntupleType> bool twobSBCut(ntupleType* ntuple) {
  float avgMass = resAvgMass(ntuple);
  float massDiff = resMassDiff(ntuple);
  float deltaR = resDeltaRMax(ntuple);
  vector<int> thisNBs = numDeepBs(ntuple);
  int numLoose = thisNBs.at(0);  int numMed = thisNBs.at(1);  int numTight = thisNBs.at(2);

  return (
    baselineCut(ntuple) &&
    ((avgMass>=0.0 && avgMass<=100.0) || (avgMass>=140.0 && avgMass<=250.0)) &&
    massDiff>=0.0 && massDiff<=40.0 &&
    deltaR>=0.0 && deltaR<=2.2 &&
    numMed==2 && numTight>=2
 );
}

////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - photon specializations - - - - - - - - - - - - //
////////////////////////////////////////////////////////////////////////
template<typename ntupleType> double fillPhotonPt(ntupleType* ntuple) {
  if (ntuple->Photons->size() == 0) return -999.;
  else return ntuple->Photons->at(0).Pt();
}

/////////////////////////////////////////////////
// - - - - - - - Trigger Cuts - - - - - - - -  //
/////////////////////////////////////////////////
template<typename ntupleType> bool signalTriggerCut(ntupleType* ntuple) {
  return ntuple->TriggerPass->at(42) == 1 || ntuple->TriggerPass->at(43) == 1 || ntuple->TriggerPass->at(44) == 1 || ntuple->TriggerPass->at(45) == 1 ;
}

template<typename ntupleType> bool singleMuTriggerCut(ntupleType* ntuple) {
  return (ntuple->TriggerPass->at(20)==1 || ntuple->TriggerPass->at(21)==1 || ntuple->TriggerPass->at(22)==1 || ntuple->TriggerPass->at(23)==1 || ntuple->TriggerPass->at(24)==1 || ntuple->TriggerPass->at(28)==1 || ntuple->TriggerPass->at(29)==1);
}

template<typename ntupleType> bool singleEleTriggerCut(ntupleType* ntuple) {
  return ntuple->TriggerPass->at(6) == 1 || ntuple->TriggerPass->at(7) == 1 || ntuple->TriggerPass->at(8) == 1 || ntuple->TriggerPass->at(9) == 1 || ntuple->TriggerPass->at(10) == 1 || ntuple->TriggerPass->at(11) == 1 || ntuple->TriggerPass->at(12) == 1 || ntuple->TriggerPass->at(13) == 1 || ntuple->TriggerPass->at(14) == 1 ;
}

template<typename ntupleType> bool lowDphiTriggerCut(ntupleType* ntuple) {
  return ntuple->TriggerPass->at(42) == 1 || ntuple->TriggerPass->at(43) == 1 || ntuple->TriggerPass->at(44) == 1 || ntuple->TriggerPass->at(45) == 1 ;
}

template<typename ntupleType> bool photonTriggerCut(ntupleType* ntuple) {
  return ntuple->TriggerPass->at(53) == 1 || ntuple->TriggerPass->at(54) == 1 ; // || ntuple->TriggerPass->at(51) || ntuple->TriggerPass->at(52);
}

template<typename ntupleType> int getClosestGenHiggses(ntupleType* ntuple, double jeteta, double jetphi) {
  float dRMin=999999.;
  for (int i=0 ; i < ntuple->GenParticles->size(); i++) {
    if (ntuple->GenParticles_PdgId->at(i) == 25 &&
    ntuple->GenParticles_ParentId->at(i) == 1000023 &&
    ntuple->GenParticles_Status->at(i) == 22) {
      float dR=sqrt((jeteta-ntuple->GenParticles->at(i).Eta())*(jeteta-ntuple->GenParticles->at(i).Eta()) +(jetphi-ntuple->GenParticles->at(i).Phi())*(jetphi-ntuple->GenParticles->at(i).Phi()));
      if (dRMin>dR)dRMin=dR;
    }
  }
  return dRMin;
}

template<typename ntupleType> int getClosestGenZ(ntupleType* ntuple, double jeteta, double jetphi) {
  float dRMin=999999.;
  for (int i=0 ; i < ntuple->GenParticles->size(); i++) {
    if (ntuple->GenParticles_PdgId->at(i) == 23 &&
    ntuple->GenParticles_ParentId->at(i) == 1000023 &&
    ntuple->GenParticles_Status->at(i) == 22) {
      float dR=sqrt((jeteta-ntuple->GenParticles->at(i).Eta())*(jeteta-ntuple->GenParticles->at(i).Eta()) +(jetphi-ntuple->GenParticles->at(i).Phi())*(jetphi-ntuple->GenParticles->at(i).Phi()));
      if (dRMin>dR)dRMin=dR;
    }
  }
  return dRMin;
}

template<typename ntupleType>double ResolutionSmear(ntupleType* ntuple, int j,unsigned int seed, bool SFUp=false) {
  TRandom3 rand(seed);
  double sigmaJMR=0;
  if (ntuple->JetsAK8_NumBhadrons->at(j)!=2) return ntuple->JetsAK8_softDropMass->at(j);
  if (ntuple->JetsAK8->at(j).Pt()>300. && ntuple->JetsAK8->at(j).Pt()<=600.) sigmaJMR=12.55;
  if (ntuple->JetsAK8->at(j).Pt()>600. && ntuple->JetsAK8->at(j).Pt()<=800.) sigmaJMR=10.24;
  if (ntuple->JetsAK8->at(j).Pt()>800. && ntuple->JetsAK8->at(j).Pt()<=1000.) sigmaJMR= 9.85;
  if (ntuple->JetsAK8->at(j).Pt()>1000.) sigmaJMR=9.44;
  sigmaJMR = sigmaJMR/110.;
  double sigmaJMRSF = 1.23;
  if (SFUp) sigmaJMRSF = sigmaJMRSF-0.18;
  double dRHiggs = getClosestGenHiggses(ntuple, ntuple->JetsAK8->at(j).Eta(), ntuple->JetsAK8->at(j).Phi());
  double dRZ = getClosestGenZ(ntuple, ntuple->JetsAK8->at(j).Eta(), ntuple->JetsAK8->at(j).Phi());
  if (dRHiggs > dRZ) {
    if (ntuple->JetsAK8->at(j).Pt()>300. && ntuple->JetsAK8->at(j).Pt()<=600.) sigmaJMR = 8.27;
    if (ntuple->JetsAK8->at(j).Pt()>600. && ntuple->JetsAK8->at(j).Pt()<=800.) sigmaJMR = 7.13;
    if (ntuple->JetsAK8->at(j).Pt()>800. && ntuple->JetsAK8->at(j).Pt()<=1000.) sigmaJMR = 6.83;
    if (ntuple->JetsAK8->at(j).Pt()>1000.) sigmaJMR = 6.90;
    sigmaJMR = sigmaJMR/83.;
  }
  double gausSmear = rand.Gaus(0, sigmaJMR)*sqrt((sigmaJMRSF*sigmaJMRSF -1));
  double smearmass = (gausSmear+1.)*ntuple->JetsAK8_softDropMass->at(j);
  return smearmass;
}
template<typename ntupleType>double SignalISRCorrection(ntupleType* ntuple) {
  float ISRWeights[7] = {1.0, 0.920, 0.821, 0.715, 0.662, 0.561,0.511};
  if (ntuple->NJetsISR==0) return ISRWeights[0];
  if (ntuple->NJetsISR==1) return ISRWeights[1];
  if (ntuple->NJetsISR==2) return ISRWeights[2];
  if (ntuple->NJetsISR==3) return ISRWeights[3];
  if (ntuple->NJetsISR==4) return ISRWeights[4];
  if (ntuple->NJetsISR==5) return ISRWeights[5];
  if (ntuple->NJetsISR>=6) return ISRWeights[6];
}
template<typename ntupleType> double doubleBSF(ntupleType* ntuple,int j) {
  double doubleBSF=1.0;
  if (ntuple->JetsAK8->at(j).Pt()>300. && ntuple->JetsAK8->at(j).Pt()<=350)doubleBSF=0.96;
  if (ntuple->JetsAK8->at(j).Pt()>350. && ntuple->JetsAK8->at(j).Pt()<=430) doubleBSF=1.00;
  if (ntuple->JetsAK8->at(j).Pt()>430.) doubleBSF=1.01;
  return doubleBSF;
}
template<typename ntupleType> double doubleBSFUp(ntupleType* ntuple,int j) {
  double doubleBSF=1.0;
  if (ntuple->JetsAK8->at(j).Pt()>300. && ntuple->JetsAK8->at(j).Pt()<=350) doubleBSF=0.96+0.03;
  if (ntuple->JetsAK8->at(j).Pt()>350. && ntuple->JetsAK8->at(j).Pt()<=430) doubleBSF=1.00+0.024;
  if (ntuple->JetsAK8->at(j).Pt()>430. && ntuple->JetsAK8->at(j).Pt()<840) doubleBSF=1.01+0.02;
  if (ntuple->JetsAK8->at(j).Pt()>840) doubleBSF=1.01+0.04;
  return doubleBSF;
}
template<typename ntupleType> double doubleBSFDn(ntupleType* ntuple,int j) {
  double doubleBSF=1.0;
  if (ntuple->JetsAK8->at(j).Pt()>300. && ntuple->JetsAK8->at(j).Pt()<=350) doubleBSF=0.96-0.02;
  if (ntuple->JetsAK8->at(j).Pt()>350. && ntuple->JetsAK8->at(j).Pt()<=430) doubleBSF=1.00-0.043;
  if (ntuple->JetsAK8->at(j).Pt()>430. && ntuple->JetsAK8->at(j).Pt()<840) doubleBSF=1.01-0.04;
  if (ntuple->JetsAK8->at(j).Pt()>840) doubleBSF=1.01-0.08;
  return doubleBSF;
}
