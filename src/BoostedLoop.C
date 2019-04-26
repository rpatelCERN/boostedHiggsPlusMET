#define BoostedLoop_cxx
#include "BoostedLoop.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>


vector<float> BoostedLoop::Loop(float EventWeight, string process = "", bool save_files = true){
  vector<float> Cutflow;
  //string histname = process + "_MET";
  string histname = process;
  TH1F* hist_Boost = new TH1F(histname.c_str(),histname.c_str(),7, 0.5,7.5);
  if (fChain == 0) return Cutflow;

  TH1F* hist_bool1 = new TH1F("histbool1","histbool1",5, 0,1000);
  TH1F* hist_bool2 = new TH1F("histbool2","histbool2",5, 0,1000);
  TH1F* hist_bool3 = new TH1F("histbool3","histbool3",5, 0,1000);
  TH1F* hist_bool4 = new TH1F("histbool4","histbool4",5, 0,1000);
  TH1F* hist_bool5 = new TH1F("histbool5","histbool5",5, 0,1000);
  TH1F* hist_bool6 = new TH1F("histbool6","histbool6",5, 0,1000);
  TH1F* hist_bool7 = new TH1F("histbool7","histbool7",5, 0,1000);


  int nentries = fChain->GetEntriesFast();

  int nbytes = 0, nb = 0;
  for (int jentry=0; jentry<nentries;jentry++) {
    int ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    float EventWeightHere = Weight*EventWeight;

    if(process.find("T5qqqq")!=string::npos) {
      if (nGenHs==2 && boostBool1) hist_bool1->Fill(MET,EventWeightHere);
      if (nGenHs==2 && boostBool2) hist_bool2->Fill(MET,EventWeightHere);
      if (nGenHs==2 && boostBool3) hist_bool3->Fill(MET,EventWeightHere);
      if (nGenHs==2 && boostBool4) hist_bool4->Fill(MET,EventWeightHere);
      if (nGenHs==2 && boostBool5) hist_bool5->Fill(MET,EventWeightHere);
      if (nGenHs==2 && boostBool6) hist_bool6->Fill(MET,EventWeightHere);
      if (nGenHs==2 && boostBool7) hist_bool7->Fill(MET,EventWeightHere);
    }
    else {
      if (boostBool1) hist_bool1->Fill(MET,EventWeightHere);
      if (boostBool2) hist_bool2->Fill(MET,EventWeightHere);
      if (boostBool3) hist_bool3->Fill(MET,EventWeightHere);
      if (boostBool4) hist_bool4->Fill(MET,EventWeightHere);
      if (boostBool5) hist_bool5->Fill(MET,EventWeightHere);
      if (boostBool6) hist_bool6->Fill(MET,EventWeightHere);
      if (boostBool7) hist_bool7->Fill(MET,EventWeightHere);
    }

  } //end loop over entries

  vector<TH1F*> bool_hists_cutflow = {hist_bool1,hist_bool2,hist_bool3,hist_bool4,hist_bool5,hist_bool6,hist_bool7};
  float value = 0;

  for (int i=1; i<=7; i++){
    value =  bool_hists_cutflow[i-1]->Integral(1,bool_hists_cutflow[i-1]->GetNbinsX()+1);
    hist_Boost->SetBinContent(i,value);
    Cutflow.push_back(value);
  }

  string savename = "boostedcutflow_"+process+".root";
  if (save_files) hist_Boost->SaveAs( (savename).c_str() );
  else hist_Boost->SaveAs( (savename).c_str() );
  delete hist_bool1; delete hist_bool2; delete hist_bool3; delete hist_bool4; delete hist_bool5; delete hist_bool6; delete hist_bool7;

  return Cutflow;
} //end Loop()
