#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "THStack.h"
#include "TString.h"
#include "TPad.h"

#include <vector>
#include <map>
#include <iostream>
#include <assert.h>

#include "plotterUtils.cc"
#include "skimSamples.cc"
#include "definitions.cc"
#include "RA2bTree.cc"
#include "TriggerEfficiencySextet.cc"

using namespace std;
double FindHiggs( vector<TLorentzVector> Jets, std::vector<unsigned int >MedTags);
double BuildBCombinations( vector<TLorentzVector> Jets, std::vector<unsigned int >LooseTags, std::vector<unsigned int >MedTags,std::vector<unsigned int >TightTags);

int main(int argc, char** argv) {

  skimSamples::region reg;
  int reg_(0);
  bool looseCuts(false);
  reg = static_cast<skimSamples::region>(reg_);

  // gROOT->ProcessLine(".L tdrstyle.C");
  // gROOT->ProcessLine("setTDRStyle()");

  skimSamples* skims_ = new skimSamples(reg);
  typedef bool(*cuts)(RA2bTree*);
  vector<cuts> baselineCuts;

  baselineCuts.push_back(*baselineCut<RA2bTree>);

  skimSamples skims = *skims_;
  TString regionName;
  TString cutName="baseline";
  regionName="signal";
  TFile* outputFile = new TFile("output.root","RECREATE");

  // background MC samples - 0 lepton regions
  for ( int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++) {
    cout << skims.sampleName[iSample] << endl;
    RA2bTree* ntuple = skims.ntuples[iSample];
    int numEvents = ntuple->fChain->GetEntries();
    ntupleBranchStatus<RA2bTree>(ntuple);

    TTree* newtree = new TTree("newtree","newtree");
    int nAK8_200=0; int nAK8_300=0; int nAK8_boost=0;
    int nAK4=0; int nAK4_cuts=0;
    int nBoostedH=0; int nDiAK8=0;
    int nResolvedH=0; int nBJets=0;
    int nResolvedH_Iso=0; int nBJets_Iso=0;
    int nGenZs = 0; int nGenHs = 0;
    int BTagsL=0; int BTagsM=0; int BTagsT=0;
    float MET = 0; float HT = 0;
    float CSV_leading=0; float CSV_subleading=0;
    float CSV_3leading=0; float CSV_4leading=0;
    float res_avgM=0; float res_deltaM=0; float res_deltaRmax=0;
    float com_mBoosted=0; float com_mResolved=0; float com_deltaR=0; float com_avgM=0;
    float Weight = 0;


    vector<TLorentzVector> CleanJets;
    vector<float>  CleanJets_bDiscriminatorCSV;
    vector<TLorentzVector> BJets_Iso;
    vector<float>  BJetsIso_bDiscriminatorCSV;

    // vector<TLorentzVector> AK4BJets;
    // double HiggsCan_AvgMass=0;
    // double deltaR_max=0;
    // double HiggsCan_MassDiff=0;
    // vector<TLorentzVector> *JetsAK8;
    // vector<TLorentzVector> *Jets;
    // vector<TLorentzVector> CleanJets;
    // vector<double>  CleanJets_bDiscriminatorCSV;
    // vector<double>  *Jets_bDiscriminatorCSV;
    // double Weight = 0;
    // vector<unsigned int> BJets_Iso;
    // vector<double> boostedHiggsCan_eta;
    // vector<double> boostedHiggsCan_phi;
    // vector<double> *resolvedHCan_mass;
    // vector<double> *boostedHCan_mass;
    // // vector<double> *boostedHCan_pt;
    // vector<double> *AK8_doubleBDiscriminator;

    TBranch *b_nAK8_200, *b_nAK8_300, *b_nAK8_boost;
    TBranch *b_nAK4, *b_nAK4_cuts;
    TBranch *b_nDiAK8, *b_nBoostedH, *b_nResolvedH;
    TBranch *b_nBJets, *b_nBJets_Iso, *b_nResolvedH_Iso;
    TBranch *b_nGenZs, *b_nGenHs;
    TBranch *b_BTagsL, *b_BTagsM, *b_BTagsT;
    TBranch *b_MET, *b_HT;
    TBranch *b_CSV_leading, *b_CSV_subleading, *b_CSV_3leading, *b_CSV_4leading;
    TBranch *b_res_avgM, *b_res_deltaM, *b_res_deltaRmax;
    TBranch *b_com_mBoosted, *b_com_mResolved, *b_com_deltaR, *b_com_avgM;
    TBranch *b_Weight;
    // TBranch *b_deltaR_max;
    // TBranch *b_BTags, *b_Jets_bDiscriminatorCSV, *b_Weight;
    // // TBranch *b_boostedHCan_pt, *b_boostedHCan_mass, *b_resolvedHCan_mass;
    // TBranch *b_boostedHCan_mass, *b_resolvedHCan_mass;
    // TBranch *b_HiggsCan_MassDiff, *b_HiggsCan_AvgMass;
    // TBranch *b_deltaMAK8;
    // TBranch *b_AK8_doubleBDiscriminator;

    newtree->Branch("nAK8_200", &nAK8_200, "nAK8_200/I");
    newtree->Branch("nAK8_300", &nAK8_300, "nAK8_300/I");
    newtree->Branch("nAK8_boost", &nAK8_boost, "nAK8_boost/I");
    newtree->Branch("nAK4", &nAK4, "nAK4/I");
    newtree->Branch("nAK4_cuts", &nAK4_cuts, "nAK4_cuts/I");
    newtree->Branch("nDiAK8", &nDiAK8, "nDiAK8/I");
    newtree->Branch("nBoostedH", &nBoostedH, "nBoostedH/I");
    newtree->Branch("nResolvedH", &nResolvedH, "nResolvedH/I");
    newtree->Branch("nBJets", &nBJets, "nBJets/I");
    newtree->Branch("nBJets_Iso", &nBJets_Iso, "nBJets_Iso/I");
    newtree->Branch("nResolvedH_Iso", &nResolvedH_Iso, "nResolvedH_Iso/I");
    newtree->Branch("nGenZs", &nGenZs, "nGenZs/I");
    newtree->Branch("nGenHs", &nGenHs, "nGenHs/I");
    newtree->Branch("BTagsL", &BTagsL, "BTagsL/I");
    newtree->Branch("BTagsM", &BTagsM, "BTagsM/I");
    newtree->Branch("BTagsT", &BTagsT, "BTagsT/I");
    newtree->Branch("MET",&MET,"MET/F");
    newtree->Branch("HT",&HT,"HT/F");
    newtree->Branch("CSV_leading", &CSV_leading, "CSV_leading/F");
    newtree->Branch("CSV_subleading", &CSV_subleading, "CSV_subleading/F");
    newtree->Branch("CSV_3leading", &CSV_3leading, "CSV_3leading/F");
    newtree->Branch("CSV_4leading", &CSV_4leading, "CSV_4leading/F");
    newtree->Branch("res_avgM", &res_avgM, "res_avgM/F");
    newtree->Branch("res_deltaM", &res_deltaM, "res_deltaM/F");
    newtree->Branch("res_deltaRmax", &res_deltaRmax, "res_deltaRmax/F");
    newtree->Branch("com_mBoosted", &com_mBoosted, "com_mBoosted/F");
    newtree->Branch("com_mResolved", &com_mResolved, "com_mResolved/F");
    newtree->Branch("com_deltaR", &com_deltaR, "com_deltaR/F");
    newtree->Branch("com_avgM", &com_avgM, "com_avgM/F");
    newtree->Branch("Weight", &Weight, "Weight/F");

    // newtree->Branch("BTags", &BTags, "BTags/I");
    // newtree->Branch("deltaR_max", &deltaR_max, "deltaR_max/D");
    // newtree->Branch("Jets_bDiscriminatorCSV",&Jets_bDiscriminatorCSV);
    // newtree->Branch("AK8_doubleBDiscriminator", &AK8_doubleBDiscriminator);
    // newtree->Branch("boostedHCan_mass", &boostedHCan_mass);
    // // newtree->Branch("boostedHCan_pt", &boostedHCan_pt);
    // newtree->Branch("resolvedHCan_mass", &resolvedHCan_mass);
    // newtree->Branch("HiggsCan_AvgMass", &HiggsCan_AvgMass, "HiggsCan_AvgMass/D");
    // newtree->Branch("HiggsCan_MassDiff", &HiggsCan_MassDiff, "HiggsCan_MassDiff/D");
    // newtree->Branch("deltaMAK8", &deltaMAK8, "deltaMAK8/D");


    newtree->SetBranchAddress("nAK8_200",&nAK8_200, &b_nAK8_200);
    newtree->SetBranchAddress("nAK8_300",&nAK8_300, &b_nAK8_300);
    newtree->SetBranchAddress("nAK8_boost",&nAK8_boost, &b_nAK8_boost);
    newtree->SetBranchAddress("nAK4",&nAK4, &b_nAK4);
    newtree->SetBranchAddress("nAK4_cuts",&nAK4_cuts, &b_nAK4_cuts);
    newtree->SetBranchAddress("nDiAK8",&nDiAK8, &b_nDiAK8);
    newtree->SetBranchAddress("nBoostedH",&nBoostedH, &b_nBoostedH);
    newtree->SetBranchAddress("nResolvedH",&nResolvedH, &b_nResolvedH);
    newtree->SetBranchAddress("nBJets", &nBJets, &b_nBJets);
    newtree->SetBranchAddress("nBJets_Iso", &nBJets_Iso, &b_nBJets_Iso);
    newtree->SetBranchAddress("nResolvedH_Iso",&nResolvedH_Iso, &b_nResolvedH_Iso);//added
    newtree->SetBranchAddress("nGenZs",&nGenZs, &b_nGenZs);
    newtree->SetBranchAddress("nGenHs",&nGenHs, &b_nGenHs);
    newtree->SetBranchAddress("BTagsL", &BTagsL, &b_BTagsL);
    newtree->SetBranchAddress("BTagsM", &BTagsM, &b_BTagsM);
    newtree->SetBranchAddress("BTagsT", &BTagsT, &b_BTagsT);
    newtree->SetBranchAddress("MET", &MET, &b_MET);
    newtree->SetBranchAddress("HT", &HT, &b_HT);
    newtree->SetBranchAddress("Weight", &Weight, &b_Weight);
    newtree->SetBranchAddress("CSV_leading", &CSV_leading, &b_CSV_leading);
    newtree->SetBranchAddress("CSV_subleading", &CSV_subleading, &b_CSV_subleading);
    newtree->SetBranchAddress("CSV_3leading", &CSV_3leading, &b_CSV_3leading);
    newtree->SetBranchAddress("CSV_4leading", &CSV_4leading, &b_CSV_4leading);

    newtree->SetBranchAddress("res_avgM", &res_avgM, &b_res_avgM);
    newtree->SetBranchAddress("res_deltaM", &res_deltaM, &b_res_deltaM);
    newtree->SetBranchAddress("res_deltaRmax", &res_deltaRmax, &b_res_deltaRmax);
    newtree->SetBranchAddress("com_mBoosted", &com_mBoosted, &b_com_mBoosted);
    newtree->SetBranchAddress("com_mResolved", &com_mResolved, &b_com_mResolved);
    newtree->SetBranchAddress("com_deltaR", &com_deltaR, &b_com_deltaR);
    newtree->SetBranchAddress("com_avgM", &com_avgM, &b_com_avgM);

    // // newtree->SetBranchAddress("BTags", &BTags, &b_BTags);
    // newtree->SetBranchAddress("Jets_bDiscriminatorCSV",&Jets_bDiscriminatorCSV, &b_Jets_bDiscriminatorCSV);
    // newtree->SetBranchAddress("deltaR_max", &deltaR_max, &b_deltaR_max);
    // newtree->SetBranchAddress("resolvedHCan_mass", &resolvedHCan_mass, &b_resolvedHCan_mass);
    // newtree->SetBranchAddress("boostedHCan_mass", &boostedHCan_mass, &b_boostedHCan_mass);
    // // newtree->SetBranchAddress("boostedHCan_pt", &boostedHCan_pt, &b_boostedHCan_pt);
    // newtree->SetBranchAddress("HiggsCan_AvgMass", &HiggsCan_AvgMass, &b_HiggsCan_AvgMass);
    // newtree->SetBranchAddress("HiggsCan_MassDiff", &HiggsCan_MassDiff, &b_HiggsCan_MassDiff);
    // newtree->SetBranchAddress("deltaMAK8", &deltaMAK8, &b_deltaMAK8);
    // newtree->SetBranchAddress("AK8_doubleBDiscriminator", &AK8_doubleBDiscriminator, &b_AK8_doubleBDiscriminator);


    float minAK8_JetPt=300;
    float bbtagCut=0.3;
    TString filename = ntuple->fChain->GetFile()->GetName();
    for ( int iEvt = 0 ; iEvt < numEvents ; iEvt++ ) {
      ntuple->GetEntry(iEvt);
      if ( iEvt % 100000 == 0 ) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << numEvents << endl;

      if (!baselineCut(ntuple)) continue; //at least 2 b jets or 1 AK8 jet, MET>150, some filter cuts
      if ( ( filename.Contains("SingleLept") || filename.Contains("DiLept") ) && ntuple->madHT>600. ) continue;
      MET=ntuple->MET;
      HT=ntuple->HT;
      Weight=ntuple->Weight;
      nGenZs = getNumGenZs(ntuple);
      nGenHs = getNumGenHiggses(ntuple);


      nAK8_200= 0;
      nAK8_300= 0;
      nAK8_boost = 0;
      nAK4= 0;
      nAK4_cuts = 0;
      nBoostedH= 0;
      nDiAK8 = 0;
      nResolvedH= 0;
      nBJets = 0;
      nResolvedH_Iso= 0;
      nBJets_Iso = 0;
      double com_boost_eta, com_boost_phi;


      BJets_Iso.clear();
      CleanJets.clear();
      CleanJets_bDiscriminatorCSV.clear();
      BJetsIso_bDiscriminatorCSV.clear();



      if (MET<300) bbtagCut=0.6;
      // Counting nAK8 jets
      for (unsigned int fj=0; fj<ntuple->JetsAK8->size();++fj) {
        double pT = ntuple->JetsAK8->at(fj).Pt();
        if (pT<200) continue;
        nAK8_200++;
        if (pT<minAK8_JetPt) continue;
        nAK8_300++;
        if (ntuple->JetsAK8_doubleBDiscriminator->at(fj)<bbtagCut) continue;
        nAK8_boost++;
      }

      // Stand-alone boosted case
      if ( boostedBaselineCut(ntuple) ) {
        double J1_pT = ntuple->JetsAK8->at(0).Pt();double J2_pT = ntuple->JetsAK8->at(1).Pt();
        double J1_phi = ntuple->JetsAK8->at(0).Phi();double J2_phi = ntuple->JetsAK8->at(1).Phi();
        double J1_eta = ntuple->JetsAK8->at(0).Eta();double J2_eta = ntuple->JetsAK8->at(1).Eta();
        double J1_prunedMass = ntuple->JetsAK8_prunedMass->at(0);double J2_prunedMass = ntuple->JetsAK8_prunedMass->at(1);
        double J1_doubleB = ntuple->JetsAK8_doubleBDiscriminator->at(0);double J2_doubleB = ntuple->JetsAK8_doubleBDiscriminator->at(1);

        bool J1_isBoost = (J1_pT>300.0 && J1_prunedMass>=85.0 && J1_prunedMass<135.0 && J1_doubleB>bbtagCut);
        bool J2_isBoost = (J2_pT>300.0 && J2_prunedMass>=85.0 && J2_prunedMass<135.0 && J2_doubleB>bbtagCut);

        bool J1_isBoost_antiB = (J1_pT>300.0 && J1_prunedMass>=85.0 && J1_prunedMass<135.0 && J1_doubleB<=bbtagCut);
        bool J2_isBoost_antiB = (J2_pT>300.0 && J2_prunedMass>=85.0 && J2_prunedMass<135.0 && J2_doubleB<=bbtagCut);

        if (J1_isBoost && J2_isBoost) {
          nBoostedH=2;
        }
        else if (J1_isBoost && J2_isBoost_antiB) {
          nBoostedH=1; nDiAK8=1;
        }
        else if (J2_isBoost && J1_isBoost_antiB) {
          nBoostedH=1; nDiAK8=1;
        }
        // else if (J1_isBoost){
        //   nBoostedH=1;
        //   com_mBoosted = J1_prunedMass;
        //   com_boost_eta = J1_eta;
        //   com_boost_phi = J1_phi;
        // }
        // else if (J2_isBoost){
        //   nBoostedH=1;
        //   com_mBoosted = J2_prunedMass;
        //   com_boost_eta = J2_eta;
        //   com_boost_phi = J2_phi;
        // }
      }


      //Stand-alone semi-resolved case, looking for a single boosted candidate (semi-lead fails the pT, or there's only one AK8)
      //baseline has two cases: 2+ AK8 jets, and only 1
      //2+ AK8 looks at lead and sublead jets - only one can pass the pT cut. If only one passes, check on it's mass and double-b cut
      //Only 1 AK8: check if passes pT, mass, and double-b cuts
      if ( semiResBaselineCut(ntuple) ) {

        bool J1_isBoost = false; bool J2_isBoost = false;

        double J1_pT = ntuple->JetsAK8->at(0).Pt();
        double J1_phi = ntuple->JetsAK8->at(0).Phi();
        double J1_eta = ntuple->JetsAK8->at(0).Eta();
        double J1_prunedMass = ntuple->JetsAK8_prunedMass->at(0);
        double J1_doubleB = ntuple->JetsAK8_doubleBDiscriminator->at(0);

        J1_isBoost = (J1_pT>300.0 && J1_prunedMass>=85.0 && J1_prunedMass<135.0 && J1_doubleB>bbtagCut);

        double J2_pT; double J2_phi; double J2_eta; double J2_prunedMass; double J2_doubleB;

        if (ntuple->JetsAK8->size() > 1 && !J1_isBoost) {
          double J2_pT = ntuple->JetsAK8->at(1).Pt();
          double J2_phi = ntuple->JetsAK8->at(1).Phi();
          double J2_eta = ntuple->JetsAK8->at(1).Eta();
          double J2_prunedMass = ntuple->JetsAK8_prunedMass->at(1);
          double J2_doubleB = ntuple->JetsAK8_doubleBDiscriminator->at(1);
          J2_isBoost = (J2_pT>300.0 && J2_prunedMass>=85.0 && J2_prunedMass<135.0 && J2_doubleB>bbtagCut);
        }

        if (J1_isBoost){
          nBoostedH=1;
          com_mBoosted = J1_prunedMass;
          com_boost_eta = J1_eta;
          com_boost_phi = J1_phi;
        }
        else if (J2_isBoost){
          nBoostedH=1;
          com_mBoosted = J2_prunedMass;
          com_boost_eta = J2_eta;
          com_boost_phi = J2_phi;
        }
      }



      //Stand-alone resolved case
      //based on Moriond Recommendations for combined CSV: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Recommendation_for_combining_b_a
      double CSVBtagLoose = 0.5426;
      double CSVBtagMed   = 0.8484;
      double CSVBtagTight = 0.9535;
      BTagsL= 0;
      BTagsM= 0;
      BTagsT= 0;

      //Save NBTags for different points, save clean jet collection, determine overlap with AK8 jets in the case of 1 boosted
      for (unsigned int j=0; j<ntuple->Jets->size();++j) {
        ++nAK4;
        if (ntuple->Jets->at(j).Pt()<30 || fabs(ntuple->Jets->at(j).Eta())>2.4) continue;
        ++nAK4_cuts;
        float CSV = ntuple->Jets_bDiscriminatorCSV->at(j);
        if (CSV > CSVBtagLoose) BTagsL++;
        if (CSV > CSVBtagMed) {
          BTagsM++;
          nBJets++;
        }
        if (CSV > CSVBtagTight) BTagsT++;
        CleanJets.push_back(ntuple->Jets->at(j));
        CleanJets_bDiscriminatorCSV.push_back(ntuple->Jets_bDiscriminatorCSV->at(j));
        if (ntuple->Jets_bDiscriminatorCSV->at(j)<CSVBtagMed) continue;
        // AK4BJets.push_back(ntuple->Jets->at(j));

        //Check overlap between boosted-AK8 and AK4 jet collection for the semi-resolved case only
        if (nBoostedH==1 && nDiAK8!=1) {
          bool Overlap=false;

          float deltaEta = ntuple->Jets->at(j).Eta() - com_boost_eta;
          float deltaPhi = CalcdPhi(ntuple->Jets->at(j).Phi(), com_boost_phi);
          float dR = sqrt((deltaEta*deltaEta)+(deltaPhi*deltaPhi));
          if (dR<0.8) {
            Overlap = true; continue;
          }
          if (!Overlap) {
            BJets_Iso.push_back(ntuple->Jets->at(j));
            nBJets_Iso++; //Num of b-tagged AK4 jets that are isolated from the AK8 jets, with pT and eta cuts
            BJetsIso_bDiscriminatorCSV.push_back(ntuple->Jets_bDiscriminatorCSV->at(j));
          }
        } //For combined case
      } //end loop over AK4 jets


      //Look over all possible resolved Higgs candidates, without requiring any isolation (completely ignore presence of AK8 jets)
      //So this is the Resolved Stand-alone
      if (BTagsT>1 && CleanJets.size()>=4 && CleanJets.size()<6) {
        //Find the four jets with the highest b-discriminator
        float HighestValuesTest[] = {-11,-11,-11,-11};
        int IndicesTest[] = {-1,-1,-1,-1};
        for (unsigned int j=0; j<CleanJets.size();++j) {
          float *current_min = min_element(HighestValuesTest,HighestValuesTest+4);
          int min_pos = distance(HighestValuesTest,min_element(HighestValuesTest,HighestValuesTest+4));
          float *CSVofJet = &(CleanJets_bDiscriminatorCSV.at(j));
          if (*CSVofJet > *current_min) {
            HighestValuesTest[min_pos] = *CSVofJet;
            IndicesTest[min_pos] = j;
          }
        } //end loop to find jets with 4 highest CSV

        //Save CSV values for problem-solving things
        std::sort(HighestValuesTest, HighestValuesTest+4, std::greater<double>());
        CSV_leading = HighestValuesTest[0];
        CSV_subleading = HighestValuesTest[1];
        CSV_3leading = HighestValuesTest[2];
        CSV_4leading = HighestValuesTest[3];

        //loop through candidate combinations
        TLorentzVector JetComboTest1a = CleanJets[IndicesTest[0]]+CleanJets[IndicesTest[1]];
        TLorentzVector JetComboTest1b = CleanJets[IndicesTest[2]]+CleanJets[IndicesTest[3]];
        TLorentzVector JetComboTest2a = CleanJets[IndicesTest[0]]+CleanJets[IndicesTest[2]];
        TLorentzVector JetComboTest2b = CleanJets[IndicesTest[1]]+CleanJets[IndicesTest[3]];
        TLorentzVector JetComboTest3a = CleanJets[IndicesTest[0]]+CleanJets[IndicesTest[3]];
        TLorentzVector JetComboTest3b = CleanJets[IndicesTest[1]]+CleanJets[IndicesTest[2]];

        double MassDiff1 = abs(JetComboTest1a.M()-JetComboTest1b.M());
        double MassDiff2 = abs(JetComboTest2a.M()-JetComboTest2b.M());
        double MassDiff3 = abs(JetComboTest3a.M()-JetComboTest3b.M());
        double MassAvg1 = (JetComboTest1a.M()+JetComboTest1b.M())/2;
        double MassAvg2 = (JetComboTest2a.M()+JetComboTest2b.M())/2;
        double MassAvg3 = (JetComboTest3a.M()+JetComboTest3b.M())/2;
        double deltaR1, deltaR2, deltaEta1, deltaEta2, deltaPhi1, deltaPhi2 = -1;
        double thisMassDiff, thisMassAvg;
        double resCanMass1, resCanMass2;

        if (MassDiff1<MassDiff2 && MassDiff1<MassDiff3){
          deltaEta1 = (CleanJets[IndicesTest[0]].Eta()-CleanJets[IndicesTest[1]].Eta());
          deltaPhi1 = CalcdPhi(CleanJets[IndicesTest[0]].Phi(),CleanJets[IndicesTest[1]].Phi());
          deltaEta2 = (CleanJets[IndicesTest[2]].Eta()-CleanJets[IndicesTest[3]].Eta());
          deltaPhi2 = CalcdPhi(CleanJets[IndicesTest[2]].Phi(),CleanJets[IndicesTest[3]].Phi());
          thisMassDiff = MassDiff1; thisMassAvg = MassAvg1;
          resCanMass1 = JetComboTest1a.M(); resCanMass2 = JetComboTest1b.M();
        }
        else if (MassDiff2<MassDiff1 && MassDiff2<MassDiff3){
          deltaEta1 = (CleanJets[IndicesTest[0]].Eta()-CleanJets[IndicesTest[2]].Eta());
          deltaPhi1 = CalcdPhi(CleanJets[IndicesTest[0]].Phi(),CleanJets[IndicesTest[2]].Phi());
          deltaEta2 = (CleanJets[IndicesTest[1]].Eta()-CleanJets[IndicesTest[3]].Eta());
          deltaPhi2 = CalcdPhi(CleanJets[IndicesTest[1]].Phi(),CleanJets[IndicesTest[3]].Phi());
          thisMassDiff = MassDiff2; thisMassAvg = MassAvg2;
          resCanMass1 = JetComboTest2a.M(); resCanMass2 = JetComboTest2b.M();
        }
        else if (MassDiff3<MassDiff1 && MassDiff3<MassDiff2){
          deltaEta1 = (CleanJets[IndicesTest[0]].Eta()-CleanJets[IndicesTest[3]].Eta());
          deltaPhi1 = CalcdPhi(CleanJets[IndicesTest[0]].Phi(),CleanJets[IndicesTest[3]].Phi());
          deltaEta2 = (CleanJets[IndicesTest[1]].Eta()-CleanJets[IndicesTest[2]].Eta());
          deltaPhi2 = CalcdPhi(CleanJets[IndicesTest[1]].Phi(),CleanJets[IndicesTest[2]].Phi());
          thisMassDiff = MassDiff3; thisMassAvg = MassAvg3;
          resCanMass1 = JetComboTest3a.M(); resCanMass2 = JetComboTest3b.M();
        }

        deltaR1 = sqrt((deltaEta1*deltaEta1)+(deltaPhi1*deltaPhi1));
        deltaR2 = sqrt((deltaEta2*deltaEta2)+(deltaPhi2*deltaPhi2));

        //only cuts here are basline, 4-5 AK4 jets (w cuts), and 2 tight b's
        res_deltaRmax = max(deltaR1, deltaR2);
        res_avgM = thisMassAvg; res_deltaM = thisMassDiff;

        if (thisMassAvg>100 && thisMassAvg<=140 && abs(thisMassDiff)<40) {
          nResolvedH=2;
        }
        else if ((resCanMass1>100 && resCanMass1<140) || (resCanMass2>100 && resCanMass2<140) ) nResolvedH=1;
      } // end Resolved baseline




        // Combined case - When there is one boosted, first look for a resolved candidate
        if (nBoostedH==1 && nDiAK8!=1 && BTagsM>=2 && BJets_Iso.size()>0) {
          //Find the four jets with the highest b-discriminator
          float HighestCSVs[] = {-11,-11,-11,-11};
          int JetIndex[] = {-1,-1,-1,-1};

          //Loop over the b-jets isolated from the boosted AK8
          for (unsigned int j=0; j<BJets_Iso.size();++j) {
            float *current_min = min_element(HighestCSVs,HighestCSVs+4);
            int min_pos = distance(HighestCSVs,min_element(HighestCSVs,HighestCSVs+4));
            float *CSVofJet = &(BJetsIso_bDiscriminatorCSV.at(j));
            if (*CSVofJet > *current_min) {
              HighestCSVs[min_pos] = *CSVofJet;
              JetIndex[min_pos] = j;
            }
          } //end loop to find jets with 4 highest CSV

          std::sort(HighestCSVs, HighestCSVs+4, std::greater<double>());
          float CSV0 = HighestCSVs[0]; float CSV1 = HighestCSVs[1];
          float CSV2 = HighestCSVs[2]; float CSV3 = HighestCSVs[3];

          vector<float> comboMass = {-11,-11,-11,-11,-11,-11};
          vector<float> comboDeltaR = {-11,-11,-11,-11,-11,-11};
          float deltaEta, deltaPhi, deltaR = -1;

          //loop through all 6 candidate combinations
          TLorentzVector JetCombo1, JetCombo2, JetCombo3, JetCombo4, JetCombo5, JetCombo6;
          if (CSV0>0 && CSV1>0) {
            JetCombo1 = BJets_Iso[JetIndex[0]]+BJets_Iso[JetIndex[1]];
            if (JetCombo1.M() > 100 && JetCombo1.M() <140) {
              comboMass.at(0)=JetCombo1.M();
              deltaEta = (BJets_Iso[JetIndex[0]].Eta()-BJets_Iso[JetIndex[1]].Eta());
              deltaPhi = CalcdPhi(BJets_Iso[JetIndex[0]].Phi(),BJets_Iso[JetIndex[1]].Phi());
              deltaR = sqrt((deltaEta*deltaEta)+(deltaPhi*deltaPhi));
              comboDeltaR.at(0)=deltaR;
            }
          }
          if (CSV2>0 && CSV3>0) {
            JetCombo2 = BJets_Iso[JetIndex[2]]+BJets_Iso[JetIndex[3]];
            if (JetCombo2.M() > 100 && JetCombo2.M() <140) {
              comboMass.at(1)=JetCombo2.M();
              deltaEta = (BJets_Iso[JetIndex[2]].Eta()-BJets_Iso[JetIndex[3]].Eta());
              deltaPhi = CalcdPhi(BJets_Iso[JetIndex[2]].Phi(),BJets_Iso[JetIndex[3]].Phi());
              deltaR = sqrt((deltaEta*deltaEta)+(deltaPhi*deltaPhi));
              comboDeltaR.at(1)=deltaR;
            }
          }
          if (CSV0>0 && CSV2>0) {
            JetCombo3 = BJets_Iso[JetIndex[0]]+BJets_Iso[JetIndex[2]];
            if (JetCombo3.M() > 100 && JetCombo3.M() <140) {
              comboMass.at(2)=JetCombo3.M();
              deltaEta = (BJets_Iso[JetIndex[0]].Eta()-BJets_Iso[JetIndex[2]].Eta());
              deltaPhi = CalcdPhi(BJets_Iso[JetIndex[0]].Phi(),BJets_Iso[JetIndex[2]].Phi());
              deltaR = sqrt((deltaEta*deltaEta)+(deltaPhi*deltaPhi));
              comboDeltaR.at(2)=deltaR;
            }
          }
          if (CSV1>0 && CSV3>0) {
            JetCombo4 = BJets_Iso[JetIndex[1]]+BJets_Iso[JetIndex[3]];
            if (JetCombo4.M() > 100 && JetCombo4.M() <140) {
              comboMass.at(3)=JetCombo4.M();
              deltaEta = (BJets_Iso[JetIndex[1]].Eta()-BJets_Iso[JetIndex[3]].Eta());
              deltaPhi = CalcdPhi(BJets_Iso[JetIndex[1]].Phi(),BJets_Iso[JetIndex[3]].Phi());
              deltaR = sqrt((deltaEta*deltaEta)+(deltaPhi*deltaPhi));
              comboDeltaR.at(3)=deltaR;
            }
          }
          if (CSV0>0 && CSV3>0) {
            JetCombo5 = BJets_Iso[JetIndex[0]]+BJets_Iso[JetIndex[3]];
            if (JetCombo5.M() > 100 && JetCombo5.M() <140) {
              comboMass.at(4)=JetCombo5.M();
              deltaEta = (BJets_Iso[JetIndex[0]].Eta()-BJets_Iso[JetIndex[3]].Eta());
              deltaPhi = CalcdPhi(BJets_Iso[JetIndex[0]].Phi(),BJets_Iso[JetIndex[3]].Phi());
              deltaR = sqrt((deltaEta*deltaEta)+(deltaPhi*deltaPhi));
              comboDeltaR.at(4)=deltaR;
            }
          }
          if (CSV1>0 && CSV2>0) {
            JetCombo6 = BJets_Iso[JetIndex[1]]+BJets_Iso[JetIndex[2]];
            if (JetCombo6.M() > 100 && JetCombo6.M() <140) {
              comboMass.at(5)=JetCombo6.M();
              deltaEta = (BJets_Iso[JetIndex[1]].Eta()-BJets_Iso[JetIndex[2]].Eta());
              deltaPhi = CalcdPhi(BJets_Iso[JetIndex[1]].Phi(),BJets_Iso[JetIndex[2]].Phi());
              deltaR = sqrt((deltaEta*deltaEta)+(deltaPhi*deltaPhi));
              comboDeltaR.at(5)=deltaR;
            }
          }

          // Determining which will be the resolved candidate saved, based on which has the smallest deltaM to the boosted
          float tmp_deltaM = 999;
          int this_combo = -1;
          for (int i=0; i<comboMass.size();i++) {
            float this_mass = comboMass.at(i);
            if (this_mass<0) continue;

            float this_deltaM = abs(this_mass - com_mBoosted);
            if (this_deltaM<tmp_deltaM) {
              tmp_deltaM = this_deltaM;
              this_combo = i;
            }
          } //end loop over combo

          if (this_combo == -1) { //if a resolved candidate cannot be found, we're just grabbing the b-jet with the highest CSV
            nResolvedH_Iso = 0; //this is implied, but keeping it here for notes-sake
            com_mResolved = BJets_Iso[JetIndex[0]].M();
            deltaEta = (BJets_Iso[JetIndex[0]].Eta() - com_boost_eta);
            deltaPhi = CalcdPhi(BJets_Iso[JetIndex[0]].Phi(), com_boost_phi);
            com_deltaR = sqrt((deltaEta*deltaEta)+(deltaPhi*deltaPhi));
          }
          else {
            nResolvedH_Iso = 1;
            com_mResolved = comboMass.at(this_combo);
            com_deltaR = comboDeltaR.at(this_combo);
            com_avgM = abs( comboMass.at(this_combo) - com_mBoosted );
          }
        } //end condition for combined case



      newtree->Fill();
      nBoostedH=0;
      nDiAK8=0;
      nAK8_200=0;
      nAK8_300=0;
      nAK8_boost=0;
      nAK4=0;
      nAK4_cuts=0;
      nBJets=0;
      nResolvedH=0;
      nResolvedH_Iso=0;
      nBJets_Iso=0;
      BTagsL=0;
      BTagsM=0;
      BTagsT=0;
      BJets_Iso.clear();
      CleanJets.clear();
      CleanJets_bDiscriminatorCSV.clear();
      BJetsIso_bDiscriminatorCSV.clear();


    } //end event loop

    outputFile->cd();
    newtree->Write(skims.sampleName[iSample]);
    delete newtree;

  } //end sample loop

  outputFile->Close();
  return 0;
}

double FindHiggs( vector<TLorentzVector> Jets, std::vector<unsigned int >MedTags) {
  TLorentzVector tempLead;
  TLorentzVector tempSubLead;
  TLorentzVector ResolvedCandidate;
  vector<double> MaybeHiggs_mass;
  vector<double> MaybeHiggs_AbsMass;
  double DasMass = -999.;
  double HiggsMass = 125.;
  for (unsigned int m=0; m<MedTags.size(); ++m) {
    unsigned int jetindex=MedTags[m];
    tempLead=Jets.at(jetindex);
    for (unsigned int n=(m+1); n<MedTags.size(); ++n) {
      unsigned int subjetindex=MedTags[n];
      tempSubLead=Jets.at(jetindex);
      ResolvedCandidate=tempLead+tempSubLead;
      if (abs(ResolvedCandidate.M()-HiggsMass)>60) continue;
      //std::cout<<"Maybe Higgs mass: "<<ResolvedCandidate.M()<<std::endl;
      MaybeHiggs_mass.push_back(ResolvedCandidate.M());
      MaybeHiggs_AbsMass.push_back(abs(ResolvedCandidate.M()-HiggsMass));
    } //end loop over sub-jet
  } // end loop over "lead" jet

  if (MaybeHiggs_AbsMass.size()>0) {
    //int size = (int)MaybeHiggs_AbsMass.size();
    auto min = min_element(MaybeHiggs_AbsMass.begin(),MaybeHiggs_AbsMass.end());
    int min_pos = distance(MaybeHiggs_AbsMass.begin(),min);
    //std::cout<<"Jet pairing with mass closest to Higgs has mass of: "<< MaybeHiggs_mass[min_pos]<<std::endl;
    DasMass = MaybeHiggs_mass[min_pos];
  }
  else{
    //finally what if there are no tags? Just take two AK4 Jets and build the mass
    if (Jets.size()>1) {
      ResolvedCandidate=Jets.at(0)+Jets.at(1);
      DasMass = ResolvedCandidate.M();
    }

  }
  if (MedTags.size()>=2	|| Jets.size()>1)return DasMass;
  else return -999.;
  //	}
}

double BuildBCombinations( vector<TLorentzVector> Jets, std::vector<unsigned int >LooseTags, std::vector<unsigned int >MedTags,std::vector<unsigned int >TightTags) {
  //	if (nBoostedCandidates==1) {
  TLorentzVector tempLead;
  TLorentzVector tempSubLead;
  //Highest Purity b-tags
  for (unsigned int t=0; t<TightTags.size(); ++t) {
    unsigned int jetindex=TightTags.at(t);
    if (t==0) tempLead=Jets.at(jetindex);
    if (t==1) tempSubLead=Jets.at(jetindex); //Tight-Tight bb candidates
    if (t>1)break;
  }
  if (TightTags.size()<2) {
    for (unsigned int m=0; m<MedTags.size(); ++m) {
      unsigned int jetindex=MedTags[m];

      if (m==0 && MedTags.size()==0) tempLead=Jets.at(jetindex);
      if (m==0 && MedTags.size()==1) tempSubLead=Jets.at(jetindex); //Tight-Medium bb candidates
      if (m==1 && MedTags.size()==0) tempSubLead=Jets.at(jetindex); //Medium-Medium bb candidates
      if (m>1)break;
    }
  }
  if (TightTags.size()<1 && MedTags.size()<1) { //don't have at least 1 Medium or 1 Tight Tag
  for (unsigned int l=0; l<LooseTags.size(); ++l) {
    unsigned int jetindex=LooseTags[l];
    if (l>1)break;
    if (l==1 && (MedTags.size()==1 || TightTags.size()==1 )) tempSubLead=Jets.at(jetindex); //Medium-Loose bb candidate
    if (TightTags.size()>0) continue;
    if (MedTags.size()>0) continue;
    if (MedTags.size()==0 && l==0) tempLead=Jets.at(jetindex);
    if (l==1) tempSubLead=Jets.at(jetindex); //Loose-Loose bb candidates
  }
}
TLorentzVector ResolvedCandidate;
if (TightTags.size()+MedTags.size()>2) {//have at least 2 b-tags
  ResolvedCandidate=tempLead+tempSubLead;
}
else{
  //finally what if there are no tags? Just take two AK4 Jets and build the mass
  if (Jets.size()>1) {
    ResolvedCandidate=Jets.at(0)+Jets.at(1);
  }

}
//	if (TightTags.size()+MedTags.size()+LooseTags.size()>2   || Jets.size()>1)std::cout<<"Mass of Candidate "<<ResolvedCandidate.M()<<std::endl;
if (TightTags.size()+MedTags.size()+LooseTags.size()>=2	|| Jets.size()>1)return ResolvedCandidate.M();
else return -999.;
//	}
}
