#define myHiggsSkim_cxx
#include "myHiggsSkim.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>
#include <TString.h>
#include <iostream>
#include <fstream>



// TString directory = "BoostToRes/";
TString directory = "ResToBoost/";

float bins_together[]={150,200,300,500,700,1000};
int binnum = 5;

// float bins_boost[]={0,500,700,999999};
float bins_boost[]={300,500,700,999999};
int binnum_boost = 3;

// float bins_res[]={150,200,300,450,600};
float bins_res[]={0,200,300,450,999999};
int binnum_res = 4;


vector<TH1F*> myHiggsSkim::Loop(float EventWeight, TString process = "", bool save_files = true) {
  vector<TH1F*> dem_hists;
   if (fChain == 0) return dem_hists;

   TString histname = process + "_MET";
   TH1F *hist_Boost = new TH1F(histname,histname,binnum, bins_boost);
   TH1F *hist_DiAK8 = new TH1F(histname,histname,binnum, bins_boost);
   TH1F *hist_Res4b_lowR = new TH1F(histname,histname,binnum_res, bins_res);
   TH1F *hist_Res4b_highR = new TH1F(histname,histname,binnum_res, bins_res);
   TH1F *hist_Res3b_lowR = new TH1F(histname,histname,binnum_res, bins_res);
   TH1F *hist_Res3b_highR = new TH1F(histname,histname,binnum_res, bins_res);
   dem_hists = {hist_Boost,hist_DiAK8,hist_Res4b_lowR,hist_Res4b_highR,hist_Res3b_lowR,hist_Res3b_highR};


   int nentries = fChain->GetEntriesFast();
   int nbytes = 0, nb = 0;

   for (int jentry=0; jentry<nentries;jentry++) {
     int ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     float EventWeightHere = Weight*EventWeight;

     bool is_4b = (BTagsT>=2 && BTagsM>=3 && BTagsL>=4);
     bool is_3b = (BTagsT>=2 && BTagsM==3);
     bool is_resolved = (res_deltaRmax<2.2  && nAK4>=4 && nAK4<6 && res_deltaM<40 && res_avgM>100 && res_avgM<=140) && (is_4b||is_3b);
     // bool is_common = (HT>600 && nBoostedH==2 && nResolvedH==2 && nAK4>=4 && nAK4<6 && res_deltaM<40 && res_avgM>100 && res_avgM<=140);


     if(process.Contains("T5HH")) {
       if (nGenHs!=2) continue;

       if (is_resolved && is_4b && res_deltaRmax<1.1) hist_Res4b_lowR->Fill(MET,EventWeightHere);
       else if (is_resolved && is_4b) hist_Res4b_highR->Fill(MET,EventWeightHere);
       else if (is_resolved && is_3b && res_deltaRmax<1.1) hist_Res3b_lowR->Fill(MET,EventWeightHere);
       else if (is_resolved && is_3b) hist_Res3b_highR->Fill(MET,EventWeightHere);
       else if (nBoostedH>=2 && MET>300) hist_Boost->Fill(MET,EventWeightHere);
       else if (nBoostedH==1 && nDiAK8==1 && MET>300) hist_DiAK8->Fill(MET,EventWeightHere);
     }

     else {
       if (is_resolved && is_4b && res_deltaRmax<1.1) hist_Res4b_lowR->Fill(MET,EventWeightHere);
       else if (is_resolved && is_4b) hist_Res4b_highR->Fill(MET,EventWeightHere);
       else if (is_resolved && is_3b && res_deltaRmax<1.1) hist_Res3b_lowR->Fill(MET,EventWeightHere);
       else if (is_resolved && is_3b) hist_Res3b_highR->Fill(MET,EventWeightHere);
       else if (nBoostedH>=2 && MET>300) hist_Boost->Fill(MET,EventWeightHere);
       else if (nBoostedH==1 && nDiAK8==1 && MET>300) hist_DiAK8->Fill(MET,EventWeightHere);
     }
   } //end loop over entries

   TString savename = process+".root";
   if (save_files){
     hist_Boost->SaveAs(directory+"Boost/"+savename);
     hist_DiAK8->SaveAs(directory+"DiAK8/"+savename);
     hist_Res4b_lowR->SaveAs(directory+"Res4b_lowR/"+savename);
     hist_Res4b_highR->SaveAs(directory+"Res4b_highR/"+savename);
     hist_Res3b_lowR->SaveAs(directory+"Res3b_lowR/"+savename);
     hist_Res3b_highR->SaveAs(directory+"Res3b_highR/"+savename);
   }
 return dem_hists;
 } //end Loop()

void myHiggsSkim::eff(float EventWeight, TString process = "") {
     TString histname_ResNum = process + "_ResNum";
     TString histname_BoostNum = process + "_BoostNum";
     TString histname_ResDenom = process + "_ResDenom";
     TString histname_BoostDenom = process + "_BoostDenom";
     TH1F * hist_ResNum = new TH1F(histname_ResNum,histname_ResNum,22,0,1100); //fill with lead jet pT (if H-tagged)
     TH1F * hist_BoostNum = new TH1F(histname_BoostNum,histname_BoostNum,22,0,1100);//fill with lead jet pT (if H-tagged)
     TH1F * hist_ResDenom = new TH1F(histname_ResDenom,histname_ResDenom,22,0,1100); //fill with lead jet pT (if H-tagged)
     TH1F * hist_BoostDenom = new TH1F(histname_BoostDenom,histname_BoostDenom,22,0,1100); //fill with lead jet pT (if H-tagged)

     int nentries = fChain->GetEntriesFast();
     int nbytes = 0, nb = 0;

     for (int jentry=0; jentry<nentries;jentry++) {
       int ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       float EventWeightHere = Weight*EventWeight;

       bool is_4b = (BTagsT>=2 && BTagsM>=3 && BTagsL>=4);
       bool is_3b = (BTagsT>=2 && BTagsM==3);
       bool is_resolved = (res_deltaRmax<2.2 && nResolvedH==2 && nAK4>=4 && nAK4<6 && res_deltaM<40 && res_avgM>100 && res_avgM<=140) && (is_4b||is_3b);
       bool is_boost = nBoostedH>=2 || (nBoostedH==1 && nDiAK8==1) && MET>300;
       float lead_res_pT = H1_pt_res; if (H2_pt_res>H1_pt_res) lead_res_pT=H2_pt_res;
       if (lead_res_pT>=1100) lead_res_pT=1099.0;
       if (H1_pt_boost>=1100) H1_pt_boost=1099.0;

       if(process.Contains("T5HH")) {
         if (nGenHs==2) hist_BoostDenom->Fill(H1_pt_boost,EventWeightHere);
         if (nGenHs==2) hist_ResDenom->Fill(lead_res_pT,EventWeightHere);
         if (nGenHs==2 && is_boost) hist_BoostNum->Fill(H1_pt_boost,EventWeightHere);
         if (nGenHs==2 && is_resolved) hist_ResNum->Fill(lead_res_pT,EventWeightHere);
       }
       else {
         hist_BoostDenom->Fill(H1_pt_boost,EventWeightHere);
         hist_ResDenom->Fill(lead_res_pT,EventWeightHere);
         if (is_boost) hist_BoostNum->Fill(H1_pt_boost,EventWeightHere);
         if (is_resolved) hist_ResNum->Fill(lead_res_pT,EventWeightHere);
       }
     } //end loop over entries
}

/*
TChiHH225_ResDenom->Add(TChiHH400_ResDenom);TChiHH225_ResDenom->Add(TChiHH700_ResDenom);TChiHH225_ResDenom->Add(TChiHH900_ResDenom);TChiHH225_ResDenom->Add(TChiHH925_ResDenom);
TChiHH225_ResDenom->Add(TChiHH950_ResDenom);TChiHH225_ResDenom->Add(TChiHH975_ResDenom);TChiHH225_ResDenom->Add(TChiHH1000_ResDenom);
TChiHH225_ResNum->Add(TChiHH400_ResNum);TChiHH225_ResNum->Add(TChiHH700_ResNum);TChiHH225_ResNum->Add(TChiHH900_ResNum);TChiHH225_ResNum->Add(TChiHH925_ResNum);
TChiHH225_ResNum->Add(TChiHH950_ResNum);TChiHH225_ResNum->Add(TChiHH975_ResNum);TChiHH225_ResNum->Add(TChiHH1000_ResNum);
TChiHH225_BoostDenom->Add(TChiHH400_BoostDenom);TChiHH225_BoostDenom->Add(TChiHH700_BoostDenom);TChiHH225_BoostDenom->Add(TChiHH900_BoostDenom);TChiHH225_BoostDenom->Add(TChiHH925_BoostDenom);
TChiHH225_BoostDenom->Add(TChiHH950_BoostDenom);TChiHH225_BoostDenom->Add(TChiHH975_BoostDenom);TChiHH225_BoostDenom->Add(TChiHH1000_BoostDenom);
TChiHH225_BoostNum->Add(TChiHH400_BoostNum);TChiHH225_BoostNum->Add(TChiHH700_BoostNum);TChiHH225_BoostNum->Add(TChiHH900_BoostNum);TChiHH225_BoostNum->Add(TChiHH925_BoostNum);
TChiHH225_BoostNum->Add(TChiHH950_BoostNum);TChiHH225_BoostNum->Add(TChiHH975_BoostNum);TChiHH225_BoostNum->Add(TChiHH1000_BoostNum);

TGraphAsymmErrors resEff(TChiHH225_ResNum,TChiHH225_ResDenom)
TGraphAsymmErrors boostEff(TChiHH225_BoostNum,TChiHH225_BoostDenom)
resEff.SetLineColor(kRed); boostEff.SetLineColor(kBlack);
resEff.SetTitle(";leading H p_{T} [GeV];efficiency")
resEff.SetName("resEffName")
boostEff.SetName("boostEffName")
resEff.Draw("AP")
boostEff.Draw("P same")


TLegend *leg = new TLegend(0.65,0.60,0.9,0.8);
leg->SetFillStyle(1001);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->AddEntry(resEffName,"IsResolved","l");
leg->AddEntry(boostEffName,"IsBoosted","l");
leg->Draw("same")
*/
