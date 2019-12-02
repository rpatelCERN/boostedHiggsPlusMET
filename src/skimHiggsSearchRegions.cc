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


int main(int argc, char** argv) {
	skimSamples::region reg;
	int reg_(0);
	bool looseCuts(false);
	reg = static_cast<skimSamples::region>(reg_);

	//gROOT->ProcessLine(".L tdrstyle.C");
	//gROOT->ProcessLine("setTDRStyle()");

	skimSamples* skims_ = new skimSamples(reg);
	typedef bool(*cuts)(RA2bTree*);
	vector<cuts> baselineCuts;

	baselineCuts.push_back(*baselineCut<RA2bTree>);

	skimSamples skims = *skims_;
	TString regionName;
	TString cutName="baseline";
	regionName="signal";

	TFile* outputFile = new TFile("BoostedCutflow.root","RECREATE");

	// background MC samples - 0 lepton regions
	for ( int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++) {
		cout << skims.sampleName[iSample] << endl;
		RA2bTree* ntuple = skims.ntuples[iSample];
		int numEvents = ntuple->fChain->GetEntries();
		ntupleBranchStatus<RA2bTree>(ntuple);

		TTree* newtree = new TTree("newtree","newtree");
		int boostBool1 = 0; int boostBool2 = 0; int boostBool3 = 0;
		int boostBool4 = 0; int boostBool5 = 0; int boostBool6 = 0;
		int boostBool7 = 0;

		int resBool1 = 0; int resBool2 = 0; int resBool3 = 0;
		int resBool4 = 0; int resBool5 = 0; int resBool6 = 0;
		int resBool7 = 0; int resBool8 = 0; int resBool9 = 0;

		float MET_ntuple=0;
		float HT_ntuple = 0;
		int BTags_ntuple = 0;
		double Weight_ntuple = 0;
		int nGenHs = 0;
		int nGenZs = 0;
		vector<TLorentzVector> JetsBCan;
		vector<double>  JetsBCan_bDiscriminatorCSV;
		double CSV_leading = 0; double CSV_subleading = 0;
		double CSV_3leading = 0; double CSV_4leading = 0;
		int BTagsL = 0; int BTagsM = 0; int BTagsT = 0;
		double deltaR_max = 0;
		

		TBranch *b_MET, *b_HT;
		TBranch *b_BTags, *b_Weight;
		TBranch *b_boostbool1, *b_boostbool2, *b_boostbool3, *b_boostbool4, *b_boostbool5, *b_boostbool6, *b_boostbool7;
		TBranch *b_resBool1, *b_resBool2, *b_resBool3, *b_resBool4, *b_resBool5, *b_resBool6, *b_resBool7, *b_resBool8, *b_resBool9;
		TBranch *b_nGenZs, *b_nGenHs;


		newtree->Branch("boostBool1", &boostBool1, "boostBool1/I");
		newtree->Branch("boostBool2", &boostBool2, "boostBool2/I");
		newtree->Branch("boostBool3", &boostBool3, "boostBool3/I");
		newtree->Branch("boostBool4", &boostBool4, "boostBool4/I");
		newtree->Branch("boostBool5", &boostBool5, "boostBool5/I");
		newtree->Branch("boostBool6", &boostBool6, "boostBool6/I");
		newtree->Branch("boostBool7", &boostBool7, "boostBool7/I");

		newtree->Branch("resBool1", &resBool1, "resBool1/I");
		newtree->Branch("resBool2", &resBool2, "resBool2/I");
		newtree->Branch("resBool3", &resBool3, "resBool3/I");
		newtree->Branch("resBool4", &resBool4, "resBool4/I");
		newtree->Branch("resBool5", &resBool5, "resBool5/I");
		newtree->Branch("resBool6", &resBool6, "resBool6/I");
		newtree->Branch("resBool7", &resBool7, "resBool7/I");
		newtree->Branch("resBool8", &resBool8, "resBool8/I");
		newtree->Branch("resBool9", &resBool9, "resBool9/I");

		newtree->Branch("BTags", &BTags_ntuple, "BTags_ntuple/I");
		newtree->Branch("HT",&HT_ntuple,"HT_ntuple/F");
		newtree->Branch("MET",&MET_ntuple,"MET_ntuple/F");
		newtree->Branch("Weight", &Weight_ntuple, "Weight_ntuple/D");
		newtree->Branch("nGenZs", &nGenZs, "nGenZs/I");
		newtree->Branch("nGenHs", &nGenHs, "nGenHs/I");


		newtree->SetBranchAddress("boostBool1",&boostBool1, &b_boostbool1);
		newtree->SetBranchAddress("boostBool2",&boostBool2, &b_boostbool2);
		newtree->SetBranchAddress("boostBool3",&boostBool3, &b_boostbool3);
		newtree->SetBranchAddress("boostBool4",&boostBool4, &b_boostbool4);
		newtree->SetBranchAddress("boostBool5",&boostBool5, &b_boostbool5);
		newtree->SetBranchAddress("boostBool6",&boostBool6, &b_boostbool6);
		newtree->SetBranchAddress("boostBool7",&boostBool7, &b_boostbool7);

		newtree->SetBranchAddress("resBool1",&resBool1, &b_resBool1);
		newtree->SetBranchAddress("resBool2",&resBool2, &b_resBool2);
		newtree->SetBranchAddress("resBool3",&resBool3, &b_resBool3);
		newtree->SetBranchAddress("resBool4",&resBool4, &b_resBool4);
		newtree->SetBranchAddress("resBool5",&resBool5, &b_resBool5);
		newtree->SetBranchAddress("resBool6",&resBool6, &b_resBool6);
		newtree->SetBranchAddress("resBool7",&resBool7, &b_resBool7);
		newtree->SetBranchAddress("resBool8",&resBool8, &b_resBool8);
		newtree->SetBranchAddress("resBool9",&resBool9, &b_resBool9);

		newtree->SetBranchAddress("BTags", &BTags_ntuple, &b_BTags);
		newtree->SetBranchAddress("Weight", &Weight_ntuple, &b_Weight);
		newtree->SetBranchAddress("MET", &MET_ntuple, &b_MET);
		newtree->SetBranchAddress("HT", &HT_ntuple, &b_HT);
		newtree->SetBranchAddress("nGenZs",&nGenZs, &b_nGenZs);
		newtree->SetBranchAddress("nGenHs",&nGenHs, &b_nGenHs);



		float bbtagCut=0.7;
		TString filename = ntuple->fChain->GetFile()->GetName();
		//for ( int iEvt = 0 ; iEvt < 1000 ; iEvt++ ) {
		for ( int iEvt = 0 ; iEvt < numEvents ; iEvt++ ) {
			//if (iEvt > 100) break;
			ntuple->GetEntry(iEvt);
			if ( iEvt % 100000 == 0 ) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << numEvents << endl;


			//if( !signalTriggerCut(ntuple) ) continue; //this kills 32% of all events


			MET_ntuple=ntuple->MET;
			HT_ntuple=ntuple->HT;
			BTags_ntuple=ntuple->BTags;
			//Weight_ntuple=ntuple->Weight;
			nGenZs = getNumGenZs(ntuple);
			nGenHs = getNumGenHiggses(ntuple);

			//Adding trigger efficiency and PU weights
			double weight=0.;
			float trigWeight=1.0;
			std::vector<double> EfficiencyCenterUpDown = Eff_MetMhtSextetReal_CenterUpDown(ntuple->HT, ntuple->MHT, ntuple->NJets);
			trigWeight=EfficiencyCenterUpDown[0];
			Weight_ntuple = ntuple->Weight;//*trigWeight*customPUweights(ntuple); //do I want this or no??
			//Weight_ntuple = ntuple->Weight*trigWeight;

			JetsBCan.clear();
			JetsBCan_bDiscriminatorCSV.clear();

			//if (!baselineCut(ntuple)) continue;
			if (baselineCut(ntuple)) {
				boostBool1 = 1;
				if ( ( filename.Contains("SingleLept") || filename.Contains("DiLept") ) && ntuple->madHT>600. ) continue;



				//if (MET_ntuple>200) {
					if (MET_ntuple>300. && HT_ntuple>600.) {
					boostBool2 = 1;

					// int AK8_size = ntuple->JetsAK8->size();
					//if (MET_ntuple>250) {
					if (ntuple->JetsAK8->size()>=2) {
						boostBool3 = 1;


						if (ntuple->JetsAK8->at(0).Pt()>300. && ntuple->JetsAK8->at(1).Pt()>300. ) {
							boostBool4 = 1;

							if (ntuple->JetsAK8_softDropMass->at(0)>50. &&  ntuple->JetsAK8_softDropMass->at(0)<250. && ntuple->JetsAK8_softDropMass->at(1)>50. &&  ntuple->JetsAK8_softDropMass->at(1)<250.) {
								boostBool5 = 1;

								if (ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(0)>bbtagCut && ntuple->JetsAK8_deepDoubleBDiscriminatorH->at(1)>bbtagCut) {
									boostBool6 = 1;

									if (ntuple->JetsAK8_softDropMass->at(0)>85. &&  ntuple->JetsAK8_softDropMass->at(0)<135. && ntuple->JetsAK8_softDropMass->at(1)>85. &&  ntuple->JetsAK8_softDropMass->at(1)<135.) {
										boostBool7 = 1;
									} //end boostbool 7
								} // end boostbool 6
							} // end boostbool5
						} //end boostbool 4

					} //end boostBool 3
				} //end boostBool2
			} //end boostBool1

			//cutflow for resolved
			if (resolvedBaselineCut(ntuple)) {
			// resBool1 = 1;
			//if ( ( filename.Contains("SingleLept") || filename.Contains("DiLept") ) && ntuple->HT>600. ) continue; //check HT or madHT


			// //based on Moriond Recommendations for combined CSV: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Recommendation_for_combining_b_a
			// double CSVBtagLoose = 0.5426;
			// double CSVBtagMed   = 0.8484;
			// double CSVBtagTight = 0.9535;

			//based on SUSY Recommendations for deep CSV: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy#Supported_Algorithms_and_Operati
			double CSVBtagLoose = 0;
			double CSVBtagMed   = 0;
			double CSVBtagTight = 0;

			if (filename.Contains("MC2016") ) {
			CSVBtagLoose = 0.2217; //2016
			CSVBtagMed   = 0.6321; //2016
			CSVBtagTight = 0.8953; //2016
		}

		else if ( filename.Contains("MC2017") ) {
		CSVBtagLoose = 0.1522; //2017
		CSVBtagMed   = 0.4941; //2017
		CSVBtagTight = 0.8001; //2017
	}

	else {
	std::cout<<"File name does not contain MC2016 nor MC2017!"<<std::endl;
}

BTagsL = 0;
BTagsM = 0;
BTagsT = 0;
for(unsigned int j=0; j<ntuple->Jets->size();++j){
if(ntuple->Jets->at(j).Pt()<30 || fabs(ntuple->Jets->at(j).Eta())>2.4) continue;
JetsBCan.push_back(ntuple->Jets->at(j));

// JetsBCan_bDiscriminatorCSV.push_back(ntuple->Jets_bDiscriminatorCSV->at(j));
// float this_CSV_value = ntuple->Jets_bDiscriminatorCSV->at(j);
float this_CSV_value = ntuple->Jets_bJetTagDeepCSVprobb->at(j)+ntuple->Jets_bJetTagDeepCSVprobbb->at(j) ;
JetsBCan_bDiscriminatorCSV.push_back(this_CSV_value);

if (this_CSV_value > CSVBtagLoose) BTagsL++;
if (this_CSV_value > CSVBtagMed)   BTagsM++;
if (this_CSV_value > CSVBtagTight) BTagsT++;
}

if ( (JetsBCan.size()==4 || JetsBCan.size()==5) && BTagsT>=2) {
resBool1 = 1;

//BTags requirement and nJets>=4 already in baseline
double HighestValuesTest[] = {-11,-11,-11,-11};
int IndicesTest[] = {-1,-1,-1,-1};
for(unsigned int j=0; j<JetsBCan.size();++j){
double *current_min = min_element(HighestValuesTest,HighestValuesTest+4);
int min_pos = distance(HighestValuesTest,min_element(HighestValuesTest,HighestValuesTest+4));
double *CSVofJet = &(JetsBCan_bDiscriminatorCSV.at(j));
if (*CSVofJet > *current_min){
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

//cout<<"Test: "<<CSV_leading<<", "<<CSV_subleading<<", "<<CSV_3leading<<", "<<CSV_4leading<<endl;
//loop through candidate combinations
TLorentzVector JetComboTest1a = JetsBCan[IndicesTest[0]]+JetsBCan[IndicesTest[1]];
TLorentzVector JetComboTest1b = JetsBCan[IndicesTest[2]]+JetsBCan[IndicesTest[3]];
TLorentzVector JetComboTest2a = JetsBCan[IndicesTest[0]]+JetsBCan[IndicesTest[2]];
TLorentzVector JetComboTest2b = JetsBCan[IndicesTest[1]]+JetsBCan[IndicesTest[3]];
TLorentzVector JetComboTest3a = JetsBCan[IndicesTest[0]]+JetsBCan[IndicesTest[3]];
TLorentzVector JetComboTest3b = JetsBCan[IndicesTest[1]]+JetsBCan[IndicesTest[2]];

double MassDiff1 = abs(JetComboTest1a.M()-JetComboTest1b.M());
double MassDiff2 = abs(JetComboTest2a.M()-JetComboTest2b.M());
double MassDiff3 = abs(JetComboTest3a.M()-JetComboTest3b.M());
double MassAvg1 = (JetComboTest1a.M()+JetComboTest1b.M())/2;
double MassAvg2 = (JetComboTest2a.M()+JetComboTest2b.M())/2;
double MassAvg3 = (JetComboTest3a.M()+JetComboTest3b.M())/2;
double deltaR1 = -1;
double deltaR2 = -1;
//Find smallest mass difference of the 3 candidates
if (MassDiff1<MassDiff2 && MassDiff1<MassDiff3){
double deltaEta1 = (JetsBCan[IndicesTest[0]].Eta()-JetsBCan[IndicesTest[1]].Eta());
double deltaPhi1 = CalcdPhi(JetsBCan[IndicesTest[0]].Phi(),JetsBCan[IndicesTest[1]].Phi());
double deltaEta2 = (JetsBCan[IndicesTest[2]].Eta()-JetsBCan[IndicesTest[3]].Eta());
double deltaPhi2 = CalcdPhi(JetsBCan[IndicesTest[2]].Phi(),JetsBCan[IndicesTest[3]].Phi());
deltaR1 = sqrt((deltaEta1*deltaEta1)+(deltaPhi1*deltaPhi1));
deltaR2 = sqrt((deltaEta2*deltaEta2)+(deltaPhi2*deltaPhi2));
deltaR_max = max(deltaR1, deltaR2);
if (abs(MassDiff1)<40) { resBool2 = 1;
if (deltaR_max<2.2) { resBool3 = 1;
if (MassAvg1>100 && MassAvg1<=140){ resBool4 = 1;
if ( (BTagsT>=2 && BTagsM==3 && BTagsL==3) || (BTagsT>=2 && BTagsM>=3 && BTagsL>=4) ) { resBool5 = 1;
if (BTagsT>=2 && BTagsM>=3 && BTagsL>=4){ resBool6 = 1;

if (MET_ntuple>200){ resBool7 = 1;
if (MET_ntuple>300){ resBool8 = 1;
if (MET_ntuple>450){ resBool9 = 1;
}//end resbool9
}	//end resbool8
}//end resbool7
} //end resbool6
} //end resBool5
} //end res bool 4
} //end resbool3
}//end resbool2
}
else if (MassDiff2<MassDiff1 && MassDiff2<MassDiff3){
double deltaEta1 = (JetsBCan[IndicesTest[0]].Eta()-JetsBCan[IndicesTest[2]].Eta());
double deltaPhi1 = CalcdPhi(JetsBCan[IndicesTest[0]].Phi(),JetsBCan[IndicesTest[2]].Phi());
double deltaEta2 = (JetsBCan[IndicesTest[1]].Eta()-JetsBCan[IndicesTest[3]].Eta());
double deltaPhi2 = CalcdPhi(JetsBCan[IndicesTest[1]].Phi(),JetsBCan[IndicesTest[3]].Phi());
deltaR1 = sqrt((deltaEta1*deltaEta1)+(deltaPhi1*deltaPhi1));
deltaR2 = sqrt((deltaEta2*deltaEta2)+(deltaPhi2*deltaPhi2));
deltaR_max = max(deltaR1, deltaR2);
if (abs(MassDiff2)<40) { resBool2 = 1;
if (deltaR_max<2.2) { resBool3 = 1;
if (MassAvg2>100 && MassAvg2<=140){ resBool4 = 1;
if ( (BTagsT>=2 && BTagsM==3 && BTagsL==3) || (BTagsT>=2 && BTagsM>=3 && BTagsL>=4) ) { resBool5 = 1;
if (BTagsT>=2 && BTagsM>=3 && BTagsL>=4){ resBool6 = 1;
if (MET_ntuple>200){ resBool7 = 1;
if (MET_ntuple>300){ resBool8 = 1;
if (MET_ntuple>450){ resBool9 = 1;
}//end resbool9
}	//end resbool8
}//end resbool7
} //end resbool6
} //end resBool5
} //end res bool 4
} //end resbool3
}//end resbool2
}
else if (MassDiff3<MassDiff1 && MassDiff3<MassDiff2){
double deltaEta1 = (JetsBCan[IndicesTest[0]].Eta()-JetsBCan[IndicesTest[3]].Eta());
double deltaPhi1 = CalcdPhi(JetsBCan[IndicesTest[0]].Phi(),JetsBCan[IndicesTest[3]].Phi());
double deltaEta2 = (JetsBCan[IndicesTest[1]].Eta()-JetsBCan[IndicesTest[2]].Eta());
double deltaPhi2 = CalcdPhi(JetsBCan[IndicesTest[1]].Phi(),JetsBCan[IndicesTest[2]].Phi());
deltaR1 = sqrt((deltaEta1*deltaEta1)+(deltaPhi1*deltaPhi1));
deltaR2 = sqrt((deltaEta2*deltaEta2)+(deltaPhi2*deltaPhi2));
deltaR_max = max(deltaR1, deltaR2);
if (abs(MassDiff3)<40) { resBool2 = 1;
if (deltaR_max<2.2) { resBool3 = 1;
if (MassAvg3>100 && MassAvg3<=140){ resBool4 = 1;
if ( (BTagsT>=2 && BTagsM==3 && BTagsL==3) || (BTagsT>=2 && BTagsM>=3 && BTagsL>=4) ) { resBool5 = 1;
if (BTagsT>=2 && BTagsM>=3 && BTagsL>=4){ resBool6 = 1;
if (MET_ntuple>200){ resBool7 = 1;
if (MET_ntuple>300){ resBool8 = 1;
if (MET_ntuple>450){ resBool9 = 1;
}//end resbool9
}	//end resbool8
}//end resbool7
} //end resbool6
} //end resBool5
} //end res bool 4
} //end resbool3
}//end resbool2
}
}  //end resBool1, nJets 4-6

} //end resBool1, baseline selection

newtree->Fill();
boostBool1 = 0;
boostBool2 = 0;
boostBool3 = 0;
boostBool4 = 0;
boostBool5 = 0;
boostBool6 = 0;
boostBool7 = 0;

resBool1 = 0;
resBool2 = 0;
resBool3 = 0;
resBool4 = 0;
resBool5 = 0;
resBool6 = 0;
resBool7 = 0;
resBool8 = 0;
resBool9 = 0;
} //end event loop


outputFile->cd();
newtree->Write(skims.sampleName[iSample]);
delete newtree;
} //end sample loop
outputFile->Close();

return 0;
}
