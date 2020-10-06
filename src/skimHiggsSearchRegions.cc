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

#include "plotterUtils.cc"
#include "skimSamples.cc"
#include "definitions.cc"
#include "RA2bTree.cc"
#include "TriggerCorrector.h"


using namespace std;

int main(int argc, char** argv) {
  skimSamples::region reg;
  int reg_(0);
  bool looseCuts(false);
  reg = static_cast<skimSamples::region>(reg_);

  TriggerCorrector trigcorror;
  TriggerCorrector trigcorrorHT;
  TriggerCorrector trigcorrorFakeMHT;
  trigcorror.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSingleEl_effRun2016");
  trigcorrorFakeMHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsMHTFromSinglePho_effRun2016");
  trigcorrorHT.SetEff("../data/triggersRa2bRun2_v2_withTEffs.root","hPassMhtMet6packVsHTFromSingleEl_effRun2016");


  skimSamples* skims_ = new skimSamples(reg);
	typedef bool(*cuts)(RA2bTree*);
  vector<cuts> baselineCuts;

	baselineCuts.push_back(*baselineCut<RA2bTree>);

  skimSamples skims = *skims_;
  TString regionName;
  TString cutName="baseline";
  regionName="signal";
  // TFile* outputFile = new TFile("output_V17.root","RECREATE");
  TFile* outputFile = new TFile("output_V18.root","RECREATE");

  // background MC samples - 0 lepton regions
  for ( int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++) {

    cout << skims.sampleName[iSample] << endl;
		RA2bTree* ntuple = skims.ntuples[iSample];
    int numEvents = ntuple->fChain->GetEntries();
		ntupleBranchStatus<RA2bTree>(ntuple);

    TTree* newtree = new TTree("newtree","newtree");
    int boostedBaseline = 0; int resolvedBaseline = 0;
    int nAK4=0; //these have pT>30 and ets<2.4 cuts
    int nBoostedH=0; int nDiAK8=0; int nResolvedH=0;
    int nGenZs = 0; int nGenHs = 0;
    int BTagsL=0; int BTagsM=0; int BTagsT=0;
    float MET = 0; float HT = 0; float MHT = 0;
    float Weight = 0; float trigWeight = 0;
    float H1_pt_res = 0; float H2_pt_res = 0;
    float H1_pt_boost = 0; float H2_pt_boost = 0;

    float boost_avgM = 0; float boost_deltaM = 0;
    float boost_J1mass = 0; float boost_J2mass = 0;
    float boost_J1bbtag = 0; float boost_J2bbtag = 0;
    float res_avgM=0; float res_deltaM=0; float res_deltaRmax=0;
    float res_J1btag=0; float res_J2btag=0;
    float res_J3btag=0; float res_J4btag=0;

    // defining SR and CR for both boosted and resolved cases
    int boost_doubletagSR = 0; int boost_doubletagSB = 0;
    int boost_tagSR = 0; int boost_tagSB = 0;
    int boost_antitagSR = 0; int boost_antitagSB = 0;
    int res_4bSR = 0; int res_4bSB = 0;
    int res_3bSR = 0; int res_3bSB = 0;
    int res_2bSR = 0; int res_2bSB = 0;

    // vector<TLorentzVector> CleanJets;
    // vector<float>  CleanJets_bDiscriminatorCSV;

    TBranch *b_nAK4;
    TBranch *b_nDiAK8, *b_nBoostedH, *b_nResolvedH;
    TBranch *b_nGenZs, *b_nGenHs;
    TBranch *b_BTagsL, *b_BTagsM, *b_BTagsT;
    TBranch *b_H1_pt_res, *b_H2_pt_res;
    TBranch *b_H1_pt_boost, *b_H2_pt_boost;
    TBranch *b_MET, *b_HT, *b_MHT;
    TBranch *b_res_J1btag, *b_res_J2btag, *b_res_J3btag, *b_res_J4btag;
    TBranch *b_res_avgM, *b_res_deltaM, *b_res_deltaRmax;
    TBranch *b_Weight, *b_trigWeight ;
    TBranch *b_dPhi, *b_leptVeto;
    TBranch *b_boostedBaseline, *b_resolvedBaseline;

    TBranch *b_boost_J1mass, *b_boost_J2mass, *b_boost_J1bbtag, *b_boost_J2bbtag;
    TBranch *b_boost_doubletagSR, *b_boost_doubletagSB, *b_boost_tagSR, *b_boost_tagSB;
    TBranch *b_boost_antitagSR, *b_boost_antitagSB;
    TBranch *b_res_4bSR, *b_res_4bSB, *b_res_3bSR, *b_res_3bSB, *b_res_2bSR, *b_res_2bSB;
    TBranch *b_boost_avgM, *b_boost_deltaM;

    newtree->Branch("boostedBaseline", &boostedBaseline, "boostedBaseline/I");
    newtree->Branch("resolvedBaseline", &resolvedBaseline, "resolvedBaseline/I");
    newtree->Branch("nAK4", &nAK4, "nAK4/I");
    newtree->Branch("nDiAK8", &nDiAK8, "nDiAK8/I");
    newtree->Branch("nBoostedH", &nBoostedH, "nBoostedH/I");
    newtree->Branch("nResolvedH", &nResolvedH, "nResolvedH/I");
    newtree->Branch("MET",&MET,"MET/F");
    newtree->Branch("HT",&HT,"HT/F");
    newtree->Branch("MHT",&MHT,"MHT/F");
    newtree->Branch("Weight", &Weight, "Weight/F");
    newtree->Branch("trigWeight", &trigWeight, "trigWeight/F");
    newtree->Branch("nGenZs", &nGenZs, "nGenZs/I");
    newtree->Branch("nGenHs", &nGenHs, "nGenHs/I");
    newtree->Branch("BTagsL", &BTagsL, "BTagsL/I");
    newtree->Branch("BTagsM", &BTagsM, "BTagsM/I");
    newtree->Branch("BTagsT", &BTagsT, "BTagsT/I");
    newtree->Branch("H1_pt_res", &H1_pt_res, "H1_pt_res/F");
    newtree->Branch("H1_pt_boost", &H1_pt_boost, "H1_pt_boost/F");
    newtree->Branch("H2_pt_res", &H2_pt_res, "H2_pt_res/F");
    newtree->Branch("H2_pt_boost", &H2_pt_boost, "H2_pt_boost/F");

    newtree->Branch("res_J1btag", &res_J1btag, "res_J1btag/F");
    newtree->Branch("res_J2btag", &res_J2btag, "res_J2btag/F");
    newtree->Branch("res_J3btag", &res_J3btag, "res_J3btag/F");
    newtree->Branch("res_J4btag", &res_J4btag, "res_J4btag/F");
    newtree->Branch("res_avgM", &res_avgM, "res_avgM/F");
    newtree->Branch("res_deltaM", &res_deltaM, "res_deltaM/F");
    newtree->Branch("res_deltaRmax", &res_deltaRmax, "res_deltaRmax/F");

    newtree->Branch("boost_J1mass", &boost_J1mass, "boost_J1mass/F");
    newtree->Branch("boost_J2mass", &boost_J2mass, "boost_J2mass/F");
    newtree->Branch("boost_J1bbtag", &boost_J1mass, "boost_J1bbtag/F");
    newtree->Branch("boost_J2bbtag", &boost_J1mass, "boost_J2bbtag/F");
    newtree->Branch("boost_avgM", &boost_avgM, "boost_avgM/F");
    newtree->Branch("boost_deltaM", &boost_deltaM, "boost_deltaM/F");

    newtree->Branch("boost_doubletagSR", &boost_doubletagSR, "boost_doubletagSR/I");
    newtree->Branch("boost_doubletagSB", &boost_doubletagSB, "boost_doubletagSB/I");
    newtree->Branch("boost_tagSR", &boost_tagSR, "boost_tagSR/I");
    newtree->Branch("boost_tagSB", &boost_tagSB, "boost_tagSB/I");
    newtree->Branch("boost_antitagSR", &boost_antitagSR, "boost_antitagSR/I");
    newtree->Branch("boost_antitagSB", &boost_antitagSB, "boost_antitagSB/I");
    newtree->Branch("res_4bSR", &res_4bSR, "res_4bSR/I");
    newtree->Branch("res_4bSB", &res_4bSB, "res_4bSB/I");
    newtree->Branch("res_3bSR", &res_3bSR, "res_3bSR/I");
    newtree->Branch("res_3bSB", &res_3bSB, "res_3bSB/I");
    newtree->Branch("res_2bSR", &res_2bSR, "res_2bSR/I");
    newtree->Branch("res_2bSB", &res_2bSB, "res_2bSB/I");

    newtree->SetBranchAddress("boostedBaseline", &boostedBaseline, &b_boostedBaseline);
    newtree->SetBranchAddress("resolvedBaseline", &resolvedBaseline, &b_resolvedBaseline);
    newtree->SetBranchAddress("nAK4",&nAK4, &b_nAK4);
    newtree->SetBranchAddress("nDiAK8",&nDiAK8, &b_nDiAK8);
    newtree->SetBranchAddress("nBoostedH",&nBoostedH, &b_nBoostedH);
    newtree->SetBranchAddress("nResolvedH",&nResolvedH, &b_nResolvedH);
    newtree->SetBranchAddress("MET", &MET, &b_MET);
    newtree->SetBranchAddress("HT", &HT, &b_HT);
    newtree->SetBranchAddress("MHT", &MHT, &b_MHT);
    newtree->SetBranchAddress("Weight", &Weight, &b_Weight);
    newtree->SetBranchAddress("trigWeight", &trigWeight, &b_trigWeight);
    newtree->SetBranchAddress("nGenZs",&nGenZs, &b_nGenZs);
    newtree->SetBranchAddress("nGenHs",&nGenHs, &b_nGenHs);
    newtree->SetBranchAddress("BTagsL", &BTagsL, &b_BTagsL);
    newtree->SetBranchAddress("BTagsM", &BTagsM, &b_BTagsM);
    newtree->SetBranchAddress("BTagsT", &BTagsT, &b_BTagsT);
    newtree->SetBranchAddress("H1_pt_res", &H1_pt_res, &b_H1_pt_res);
    newtree->SetBranchAddress("H2_pt_res", &H2_pt_res, &b_H2_pt_res);
    newtree->SetBranchAddress("H1_pt_boost", &H1_pt_boost, &b_H1_pt_boost);
    newtree->SetBranchAddress("H2_pt_boost", &H2_pt_boost, &b_H2_pt_boost);

    newtree->SetBranchAddress("res_J1btag", &res_J1btag, &b_res_J1btag);
    newtree->SetBranchAddress("res_J2btag", &res_J2btag, &b_res_J2btag);
    newtree->SetBranchAddress("res_J3btag", &res_J3btag, &b_res_J3btag);
    newtree->SetBranchAddress("res_J4btag", &res_J4btag, &b_res_J4btag);
    newtree->SetBranchAddress("res_avgM", &res_avgM, &b_res_avgM);
    newtree->SetBranchAddress("res_deltaM", &res_deltaM, &b_res_deltaM);
    newtree->SetBranchAddress("res_deltaRmax", &res_deltaRmax, &b_res_deltaRmax);

    newtree->SetBranchAddress("boost_J1mass", &boost_J1mass, &b_boost_J1mass);
    newtree->SetBranchAddress("boost_J2mass", &boost_J2mass, &b_boost_J2mass);
    newtree->SetBranchAddress("boost_J1bbtag", &boost_J1bbtag, &b_boost_J1bbtag);
    newtree->SetBranchAddress("boost_J2bbtag", &boost_J2bbtag, &b_boost_J2bbtag);
    newtree->SetBranchAddress("boost_avgM", &boost_avgM, &b_boost_avgM);
    newtree->SetBranchAddress("boost_deltaM", &boost_deltaM, &b_boost_deltaM);

    newtree->SetBranchAddress("boost_doubletagSR", &boost_doubletagSR, &b_boost_doubletagSR);
    newtree->SetBranchAddress("boost_doubletagSB", &boost_doubletagSB, &b_boost_doubletagSB);
    newtree->SetBranchAddress("boost_tagSR", &boost_tagSR, &b_boost_tagSR);
    newtree->SetBranchAddress("boost_tagSB", &boost_tagSB, &b_boost_tagSB);
    newtree->SetBranchAddress("boost_antitagSR", &boost_antitagSR, &b_boost_antitagSR);
    newtree->SetBranchAddress("boost_antitagSB", &boost_antitagSB, &b_boost_antitagSB);
    newtree->SetBranchAddress("res_4bSR", &res_4bSR, &b_res_4bSR);
    newtree->SetBranchAddress("res_4bSB", &res_4bSB, &b_res_4bSB);
    newtree->SetBranchAddress("res_3bSR", &res_3bSR, &b_res_3bSR);
    newtree->SetBranchAddress("res_3bSB", &res_3bSB, &b_res_3bSB);
    newtree->SetBranchAddress("res_2bSR", &res_2bSR, &b_res_2bSR);
    newtree->SetBranchAddress("res_2bSB", &res_2bSB, &b_res_2bSB);


    float minAK8_JetPt = 300;
    float bbtagCut = 0.70;
		bool runRealDeepBB = true;
    TString filename = ntuple->fChain->GetFile()->GetName();

    //Apply trigger weights
    float triggerWeight=1.0;
    float trigunc=0;
    triggerWeight=trigcorror.GetCorrection(ntuple->MHT,trigunc);
    triggerWeight=triggerWeight*trigcorrorHT.GetCorrection(ntuple->HT,trigunc);
    if(skims.sampleName[iSample]=="QCD") triggerWeight=trigcorrorFakeMHT.GetCorrection(ntuple->MHT,trigunc);
    trigWeight=triggerWeight;

		//If using supported deep bb-tag
		// if ( filename.Contains("TChiHH_HToBB") || filename.Contains("T5qqqqZH")) runRealDeepBB = true;

    for ( int iEvt = 0 ; iEvt < numEvents ; iEvt++ ) {
      ntuple->GetEntry(iEvt);
      if ( iEvt % 100000 == 0 ) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << numEvents << endl;

      //basic basleline cut (rather 2 tight BTags and 4-5 jets or 2+ AK8 jets with pT>300)
      if (!baselineCut(ntuple)) continue;

      // if ( ( filename.Contains("SingleLept") || filename.Contains("DiLept") ) && ntuple->madHT>600. ) continue;

      MET=ntuple->MET;
      HT=ntuple->HT;
      MHT=ntuple->MHT;
      Weight=ntuple->Weight;
      nGenZs = getNumGenZs(ntuple);
      nGenHs = getNumGenHiggses(ntuple);
      std::vector<int> theseBs = numDeepBs(ntuple);
      BTagsL = theseBs.at(0);
      BTagsM = theseBs.at(1);
      BTagsT = theseBs.at(2);

	    bbtagCut=0.70;

      nAK4 = 0;
      nBoostedH= 0;
      nDiAK8 = 0;
      nResolvedH= 0;
      boostedBaseline = 0;
      resolvedBaseline = 0;

      boost_doubletagSR = 0;
      boost_doubletagSB = 0;
      boost_tagSR = 0;
      boost_tagSB = 0;
      boost_antitagSR = 0;
      boost_antitagSB = 0;
      res_4bSR = 0;
      res_4bSB = 0;
      res_3bSR = 0;
      res_3bSB = 0;
      res_2bSR = 0;
      res_2bSB = 0;

      if (boostedBaselineCut(ntuple)) {
        boostedBaseline = 1;
        if (doubletagSRCut(ntuple) ) boost_doubletagSR = 1;
        else if (tagSRCut(ntuple) ) boost_tagSR = 1;
        else if (doubletagSBCut(ntuple) ) boost_doubletagSB = 1;
        else if (tagSBCut(ntuple) ) boost_tagSB = 1;
        else if (antitagSRCut(ntuple) ) boost_antitagSR = 1;
        else if (antitagSBCut(ntuple) ) boost_antitagSB = 1;
      }
      if (resolvedBaselineCut(ntuple)) {
        resolvedBaseline = 1;
        if (fourbSRCut(ntuple)) res_4bSR = 1;
        else if (fourbSBCut(ntuple)) res_4bSB = 1;
        else if (threebSRCut(ntuple)) res_3bSR = 1;
        else if (threebSBCut(ntuple)) res_3bSB = 1;
        else if (twobSRCut(ntuple)) res_2bSR = 1;
        else if (twobSBCut(ntuple)) res_2bSB = 1;
      }


      // CleanJets.clear();
      // CleanJets_bDiscriminatorCSV.clear();

      // Stand-alone boosted case
      if ( boostedBaselineCut(ntuple) ) {
        if (MET<300) continue;
        if (ntuple->JetsAK8_softDropMass->at(0)<1.0 || ntuple->JetsAK8_softDropMass->at(1)<1.0) continue;

        double J1_pT = ntuple->JetsAK8->at(0).Pt();double J2_pT = ntuple->JetsAK8->at(1).Pt();
        double J1_phi = ntuple->JetsAK8->at(0).Phi();double J2_phi = ntuple->JetsAK8->at(1).Phi();
        double J1_eta = ntuple->JetsAK8->at(0).Eta();double J2_eta = ntuple->JetsAK8->at(1).Eta();
        double J1_softDropMass = ntuple->JetsAK8_softDropMass->at(0);double J2_softDropMass = ntuple->JetsAK8_softDropMass->at(1);

				double J1_doubleB = ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(0);
				double J2_doubleB = ntuple->JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb->at(1);


        boost_J1mass = J1_softDropMass; boost_J2mass = J2_softDropMass;
        boost_J1bbtag = J1_doubleB; boost_J2bbtag = J2_doubleB;


        bool J1_isBoost = (J1_pT>300.0 && J1_softDropMass>=95.0 && J1_softDropMass<145.0 && J1_doubleB>bbtagCut);
        bool J2_isBoost = (J2_pT>300.0 && J2_softDropMass>=95.0 && J2_softDropMass<145.0 && J2_doubleB>bbtagCut);

        bool J1_isBoost_antiB = (J1_pT>300.0 && J1_softDropMass>=95.0 && J1_softDropMass<145.0 && J1_doubleB<=bbtagCut);
        bool J2_isBoost_antiB = (J2_pT>300.0 && J2_softDropMass>=95.0 && J2_softDropMass<145.0 && J2_doubleB<=bbtagCut);

        if (J1_isBoost && J2_isBoost) {
          nBoostedH=2;
          H1_pt_boost = J1_pT;
          H2_pt_boost = J2_pT;
          boost_avgM = (J1_softDropMass+J2_softDropMass)/2;
          boost_deltaM = fabs(J1_softDropMass-J2_softDropMass);
        }
        else if (J1_isBoost && J2_isBoost_antiB) {
          nBoostedH=1; nDiAK8=1;
          H1_pt_boost = J1_pT;
          H2_pt_boost = J2_pT;
        }
        else if (J2_isBoost && J1_isBoost_antiB) {
          nBoostedH=1; nDiAK8=1;
          H1_pt_boost = J1_pT;
          H2_pt_boost = J2_pT;
        }
      }


      // //Stand-alone resolved case (using deepCSV)
		  // double CSVBtagLoose = 0.2217;
		  // double CSVBtagMed   = 0.6321;
		  // double CSVBtagTight = 0.8953;
      // BTagsL= 0;
      // BTagsM= 0;
      // BTagsT= 0;
      //
      // //Save NBTags for different points, save clean jet collection, determine overlap with AK8 jets in the case of 1 boosted
      // for (unsigned int j=0; j<ntuple->Jets->size();++j) {
      //   if (ntuple->Jets->at(j).Pt()<30 || fabs(ntuple->Jets->at(j).Eta())>2.4) continue;
      //   ++nAK4;
			// 	float CSV = ntuple->Jets_bJetTagDeepCSVprobb->at(j)+ntuple->Jets_bJetTagDeepCSVprobbb->at(j);
      //   if (CSV > CSVBtagLoose) BTagsL++;
      //   if (CSV > CSVBtagMed)   BTagsM++;
      //   if (CSV > CSVBtagTight) BTagsT++;
      //   CleanJets.push_back(ntuple->Jets->at(j));
      //   CleanJets_bDiscriminatorCSV.push_back(CSV);
      // } //end loop over AK4 jets
      //
      // if (CleanJets.size()>=4 && CleanJets.size()<6) {
      //   //Find the four jets with the highest b-discriminator
      //   float HighestValuesTest[] = {-11,-11,-11,-11};
      //   int IndicesTest[] = {-1,-1,-1,-1};
      //   for (unsigned int j=0; j<CleanJets.size();++j) {
      //     float *current_min = min_element(HighestValuesTest,HighestValuesTest+4);
      //     int min_pos = distance(HighestValuesTest,min_element(HighestValuesTest,HighestValuesTest+4));
      //     float *CSVofJet = &(CleanJets_bDiscriminatorCSV.at(j));
      //     if (*CSVofJet > *current_min) {
      //       HighestValuesTest[min_pos] = *CSVofJet;
      //       IndicesTest[min_pos] = j;
      //     }
      //   } //end loop to find jets with 4 highest CSV
      //
      //   //Save CSV values for problem-solving things
      //   std::sort(HighestValuesTest, HighestValuesTest+4, std::greater<double>());
      //   res_J1btag = HighestValuesTest[0];
      //   res_J2btag = HighestValuesTest[1];
      //   res_J3btag = HighestValuesTest[2];
      //   res_J4btag = HighestValuesTest[3];
      //
      //   //loop through candidate combinations
      //   TLorentzVector JetComboTest1a = CleanJets[IndicesTest[0]]+CleanJets[IndicesTest[1]];
      //   TLorentzVector JetComboTest1b = CleanJets[IndicesTest[2]]+CleanJets[IndicesTest[3]];
      //   TLorentzVector JetComboTest2a = CleanJets[IndicesTest[0]]+CleanJets[IndicesTest[2]];
      //   TLorentzVector JetComboTest2b = CleanJets[IndicesTest[1]]+CleanJets[IndicesTest[3]];
      //   TLorentzVector JetComboTest3a = CleanJets[IndicesTest[0]]+CleanJets[IndicesTest[3]];
      //   TLorentzVector JetComboTest3b = CleanJets[IndicesTest[1]]+CleanJets[IndicesTest[2]];
      //
      //   double MassDiff1 = abs(JetComboTest1a.M()-JetComboTest1b.M());
      //   double MassDiff2 = abs(JetComboTest2a.M()-JetComboTest2b.M());
      //   double MassDiff3 = abs(JetComboTest3a.M()-JetComboTest3b.M());
      //   double MassAvg1 = (JetComboTest1a.M()+JetComboTest1b.M())/2;
      //   double MassAvg2 = (JetComboTest2a.M()+JetComboTest2b.M())/2;
      //   double MassAvg3 = (JetComboTest3a.M()+JetComboTest3b.M())/2;
      //   double deltaR1, deltaR2, deltaEta1, deltaEta2, deltaPhi1, deltaPhi2 = -1;
      //   double thisMassDiff, thisMassAvg;
      //   double resCanMass1, resCanMass2;
      //
      //   if (MassDiff1<MassDiff2 && MassDiff1<MassDiff3){
      //     deltaEta1 = (CleanJets[IndicesTest[0]].Eta()-CleanJets[IndicesTest[1]].Eta());
      //     deltaPhi1 = CalcdPhi(CleanJets[IndicesTest[0]].Phi(),CleanJets[IndicesTest[1]].Phi());
      //     deltaEta2 = (CleanJets[IndicesTest[2]].Eta()-CleanJets[IndicesTest[3]].Eta());
      //     deltaPhi2 = CalcdPhi(CleanJets[IndicesTest[2]].Phi(),CleanJets[IndicesTest[3]].Phi());
      //     thisMassDiff = MassDiff1; thisMassAvg = MassAvg1;
      //     resCanMass1 = JetComboTest1a.M(); resCanMass2 = JetComboTest1b.M();
      //     H1_pt_res = JetComboTest1a.Pt(); H2_pt_res = JetComboTest1b.Pt();
      //   }
      //   else if (MassDiff2<MassDiff1 && MassDiff2<MassDiff3){
      //     deltaEta1 = (CleanJets[IndicesTest[0]].Eta()-CleanJets[IndicesTest[2]].Eta());
      //     deltaPhi1 = CalcdPhi(CleanJets[IndicesTest[0]].Phi(),CleanJets[IndicesTest[2]].Phi());
      //     deltaEta2 = (CleanJets[IndicesTest[1]].Eta()-CleanJets[IndicesTest[3]].Eta());
      //     deltaPhi2 = CalcdPhi(CleanJets[IndicesTest[1]].Phi(),CleanJets[IndicesTest[3]].Phi());
      //     thisMassDiff = MassDiff2; thisMassAvg = MassAvg2;
      //     resCanMass1 = JetComboTest2a.M(); resCanMass2 = JetComboTest2b.M();
      //     H1_pt_res = JetComboTest2a.Pt(); H2_pt_res = JetComboTest2b.Pt();
      //   }
      //   else if (MassDiff3<MassDiff1 && MassDiff3<MassDiff2){
      //     deltaEta1 = (CleanJets[IndicesTest[0]].Eta()-CleanJets[IndicesTest[3]].Eta());
      //     deltaPhi1 = CalcdPhi(CleanJets[IndicesTest[0]].Phi(),CleanJets[IndicesTest[3]].Phi());
      //     deltaEta2 = (CleanJets[IndicesTest[1]].Eta()-CleanJets[IndicesTest[2]].Eta());
      //     deltaPhi2 = CalcdPhi(CleanJets[IndicesTest[1]].Phi(),CleanJets[IndicesTest[2]].Phi());
      //     thisMassDiff = MassDiff3; thisMassAvg = MassAvg3;
      //     resCanMass1 = JetComboTest3a.M(); resCanMass2 = JetComboTest3b.M();
      //     H1_pt_res = JetComboTest3a.Pt(); H2_pt_res = JetComboTest3b.Pt();
      //   }
      //
      //   deltaR1 = sqrt((deltaEta1*deltaEta1)+(deltaPhi1*deltaPhi1));
      //   deltaR2 = sqrt((deltaEta2*deltaEta2)+(deltaPhi2*deltaPhi2));
      //
      //   //only cuts here are basline, 4-5 AK4 jets (w cuts), and 2 tight b's
      //   res_deltaRmax = max(deltaR1, deltaR2);
      //   res_avgM = thisMassAvg; res_deltaM = thisMassDiff;
      // } // end Resolved baseline

      newtree->Fill();
      nBoostedH=0;
      nDiAK8=0;
      nAK4=0;
      nResolvedH=0;
      BTagsL=0; BTagsM=0; BTagsT=0;

      boostedBaseline = 0;
      resolvedBaseline = 0;
      boost_doubletagSR = 0; boost_doubletagSB = 0;
      boost_tagSR = 0; boost_tagSB = 0;
      boost_antitagSR = 0; boost_antitagSB = 0;
      res_4bSR = 0; res_4bSB = 0;
      res_3bSR = 0; res_3bSB = 0;
      res_2bSR = 0; res_2bSB = 0;

      // CleanJets.clear();
      // CleanJets_bDiscriminatorCSV.clear();
    } //end event loop

    outputFile->cd();
    newtree->Write(skims.sampleName[iSample]);
    delete newtree;
  } //end sample loop

  outputFile->Close();
  return 0;
}
