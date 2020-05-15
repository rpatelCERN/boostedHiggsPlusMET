#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPave.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include "myHiggsSkim.C"

using namespace std;

void runProcesses(TString path);
void MET_stack(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets, TH1F* h_TChiHH700,TH1F* h_TChiHH1000, TString which);
void runMETStacks(TString path);


void MakeTable(TString path="./") {
  runProcesses(path);
  // runMETStacks("./ResToBoost/");
}

void runProcesses(TString path="./") {
  if(!TString(path).EndsWith("/")) path+="/";
  TFile* inf = TFile::Open(path+"output_V18.root","READ");

  bool save_files = false;
  bool save_stacks = false;
  bool make_table = false;

  // TString which_process = "all";
  TString which_process = "some";


  vector<TString> processes = {"",""};
  vector<double> nEvents_signal = {0,0};
  vector<double> nEvents_T5HH = {0,0};
  float table_array[27][6];


  if (which_process=="all") {
    cout<<"Using all higgsino mass points + backgrounds"<<endl;
    processes = {"ZJets","WJets","TT","QCD", "TChiHH127", "TChiHH150", "TChiHH175","TChiHH200","TChiHH225", "TChiHH250", "TChiHH275", "TChiHH300", "TChiHH325","TChiHH350", "TChiHH375", "TChiHH400", "TChiHH425", "TChiHH450","TChiHH475","TChiHH500", "TChiHH525", "TChiHH550","TChiHH575", "TChiHH600", "TChiHH625", "TChiHH650", "TChiHH675", "TChiHH700","TChiHH725", "TChiHH750","TChiHH775", "TChiHH800", "TChiHH825", "TChiHH850", "TChiHH875", "TChiHH900", "TChiHH925", "TChiHH950", "TChiHH975", "TChiHH1000", "T5HH750", "T5HH1000", "T5HH1100", "T5HH1200", "T5HH1300", "T5HH1400", "T5HH1500", "T5HH1600", "T5HH1700", "T5HH1800", "T5HH1900", "T5HH2000", "T5HH2100", "T5HH2200"};
    // processes = {"TChiHH127", "TChiHH150", "TChiHH175","TChiHH200","TChiHH225", "TChiHH250", "TChiHH275", "TChiHH300", "TChiHH325","TChiHH350", "TChiHH375", "TChiHH400", "TChiHH425", "TChiHH450","TChiHH475","TChiHH500", "TChiHH525", "TChiHH550","TChiHH575", "TChiHH600", "TChiHH625", "TChiHH650", "TChiHH675", "TChiHH700","TChiHH725", "TChiHH750","TChiHH775", "TChiHH800", "TChiHH825", "TChiHH850", "TChiHH875", "TChiHH900", "TChiHH925", "TChiHH950", "TChiHH975", "TChiHH1000"};

    //V16
    nEvents_signal = {84454, 106406, 99099, 102285, 99208, 91650, 98101, 92229, 109068, 90086, 91969, 106716, 100323, 87694, 102881, 94578, 89209, 98000, 102377, 100218, 100554, 97914, 100716, 99027, 99564, 98288, 103443, 101643, 87285, 91713, 96388, 104283, 110373, 104522,97574, 101476};

  }
  else if (which_process=="some") {
    // processes = {"TT"};
    processes = {"TT","ZJets","WJets","QCD","TChiHH700","TChiHH1000"};
    // processes = {"TT","ZJets","WJets","QCD"};
    // processes = {"TChiHH200","TChiHH700","TChiHH1000"};
    // processes = {"TChiHH700","TChiHH1000"};
    // processes = {"TChiHH225"};
    // processes = {"TChiHH225", "TChiHH400","TChiHH700","TChiHH1000"};
    // processes = {"T5HH1300","T5HH1700","T5HH2100"};
    // processes = {"TT","ZJets","WJets","QCD","TChiHH700","TChiHH1000","T5HH1700","T5HH2100"};
     // nEvents_signal = {102285, 99027, 101476};
     // nEvents_signal = {99027, 101476}; //V12

     // nEvents_signal = {36212.0}; //V18
     // nEvents_signal = {36212,30770,34955.0, 36171.0}; //V18

     nEvents_signal = {34955.0, 36171.0}; //V18


     nEvents_T5HH = {126575.0,133459.0,121390.0}; //V18

    //For all TChiHH signals
    // processes = {"TChiHH127", "TChiHH150", "TChiHH175","TChiHH200","TChiHH225", "TChiHH250", "TChiHH275", "TChiHH300", "TChiHH325","TChiHH350", "TChiHH375", "TChiHH400", "TChiHH425", "TChiHH450","TChiHH475","TChiHH500", "TChiHH525", "TChiHH550","TChiHH575", "TChiHH600", "TChiHH625", "TChiHH650", "TChiHH675", "TChiHH700","TChiHH725", "TChiHH750","TChiHH775", "TChiHH800", "TChiHH825", "TChiHH850", "TChiHH875", "TChiHH900", "TChiHH925", "TChiHH950", "TChiHH975", "TChiHH1000"};
    // nEvents_signal = {84454, 106406, 99099, 102285, 99208, 91650, 98101, 92229, 109068, 90086, 91969, 106716, 100323, 87694, 102881, 94578, 89209, 98000, 102377, 100218, 100554, 97914, 100716, 99027, 99564, 98288, 103443, 101643, 87285, 91713, 96388, 104283, 110373, 104522,97574, 101476};

    //For some TChiHH signals
    // processes = {"TChiHH575", "TChiHH600", "TChiHH625", "TChiHH650", "TChiHH675", "TChiHH700","TChiHH725", "TChiHH750","TChiHH775", "TChiHH800", "TChiHH825", "TChiHH850", "TChiHH875", "TChiHH900", "TChiHH925", "TChiHH950", "TChiHH975", "TChiHH1000"};
    // nEvents_signal = {102377,100218, 100554, 97914, 100716, 99027, 99564, 98288, 103443, 101643, 87285, 91713, 96388, 104283, 110373, 104522,97574, 101476};

    //All T5HH
    //processes = {"T5HH750", "T5HH1000", "T5HH1100", "T5HH1200", "T5HH1300", "T5HH1400", "T5HH1500", "T5HH1600", "T5HH1700", "T5HH1800", "T5HH1900", "T5HH2000", "T5HH2100", "T5HH2200"};
  }
  else {
    cout<<"Improper process! Use some or all "<<endl;
    return;
  }

  THStack* stack_Boost = new THStack("stack_Boost","stack_Boost");
  THStack* stack_DiAK8 = new THStack("stack_DiAK8","stack_DiAK8");
  THStack* stack_Res4b_lowR = new THStack("stack_Res4b_lowR","stack_Res4b_lowR");
  THStack* stack_Res4b_highR = new THStack("stack_Res4b_highR","stack_Res4b_highR");
  THStack* stack_Res3b_lowR = new THStack("stack_Res3b_lowR","stack_Res3b_lowR");
  THStack* stack_Res3b_highR = new THStack("stack_Res3b_highR","stack_Res3b_highR");

  TLegend *legend = new TLegend(0.6,0.6,0.9,0.9);

  TString var = "MET";
  TString direct = "";
  vector<int> colors = {kOrange+1,kGreen+2,kBlue-2,kYellow+2};
  vector<TH1F*> dem_hists;
  int signal_counter = 0;
  int color_counter = 0;
  int other_counter = 0;

  for(unsigned int p = 0; p< processes.size(); p++) {
    // double this_lumi = 137000.0;
    double this_lumi = 35862.824;
    double EventWeight = this_lumi;
    if (processes[p].Contains("TChiHH")) EventWeight = this_lumi*TMath::Power(0.582,2)/nEvents_signal[signal_counter];
    else if (processes[p].Contains("T5HH")) EventWeight = this_lumi*4.0/nEvents_T5HH[other_counter];

    cout<<"Beginning "<< processes[p]<<endl;
    TTree *tree;
    inf->GetObject(processes[p],tree);
    myHiggsSkim thisSkim;
    thisSkim.Init(tree);
    // thisSkim.eff(EventWeight, processes[p]);

    if (processes[p].Contains("TChiHH")) {
      dem_hists = thisSkim.Loop(EventWeight, processes[p], save_files);
      signal_counter++;
    }
    else if (processes[p].Contains("T5HH")) {
      dem_hists = thisSkim.Loop(EventWeight, processes[p], save_files);
      other_counter++;
    }
    else {
      dem_hists = thisSkim.Loop(EventWeight, processes[p], save_files);
    }

    if ( ( !(processes[p].Contains("TChiHH") ) || !(processes[p].Contains("T5HH")) ) && save_stacks){
      legend->AddEntry(dem_hists.back(),processes[p],"f");//for stack plot
      for (unsigned int i=0; i<dem_hists.size();i++){
        dem_hists[i]->SetLineColor(colors[color_counter]);
        dem_hists[i]->SetFillColor(colors[color_counter]);
      }
      stack_Boost->Add(dem_hists[0]);
      stack_DiAK8->Add(dem_hists[1]);
      stack_Res4b_lowR->Add(dem_hists[2]);
      stack_Res4b_highR->Add(dem_hists[3]);
      stack_Res3b_lowR->Add(dem_hists[4]);
      stack_Res3b_highR->Add(dem_hists[5]);
    }

    if (make_table){

      //2 Boost
      table_array[0][p] = dem_hists[0]->Integral(1,dem_hists[0]->GetNbinsX()+1);
      for (int i=1; i<4; i++){
        float this_bin = dem_hists[0]->GetBinContent(i);
        float next_bin = dem_hists[0]->GetBinContent(i+1);
        if (this_bin<1e-20) this_bin=0;
        if (next_bin<1e-20) next_bin=0;
        table_array[i][p] = this_bin;
        if (i==3)table_array[i][p] = this_bin+next_bin;
      }

      //1 Boost, 1 AK8
      table_array[4][p] = dem_hists[1]->Integral(1,dem_hists[1]->GetNbinsX()+1);
      for (int i=1; i<4; i++){
        float this_bin = dem_hists[1]->GetBinContent(i);
        float next_bin = dem_hists[1]->GetBinContent(i+1);
        if (this_bin<1e-20) this_bin=0;
        if (next_bin<1e-20) next_bin=0;
        table_array[i+4][p] = this_bin;
        if (i==3)table_array[i+4][p] = this_bin+next_bin;
      }

      //Res 4b, low deltaR max
      table_array[8][p] = dem_hists[2]->Integral(1,dem_hists[2]->GetNbinsX()+1);
      for (int i=1; i<5; i++){
        float this_bin = dem_hists[2]->GetBinContent(i);
        float next_bin = dem_hists[2]->GetBinContent(i+1);
        if (this_bin<1e-20) this_bin=0;
        if (next_bin<1e-20) next_bin=0;
        table_array[i+8][p] = this_bin;
        if (i==4)table_array[i+8][p] = this_bin+next_bin;
      }

      //Res 4b, high deltaR max
      table_array[13][p] = dem_hists[3]->Integral(1,dem_hists[3]->GetNbinsX()+1);
      for (int i=1; i<5; i++){
        float this_bin = dem_hists[3]->GetBinContent(i);
        float next_bin = dem_hists[3]->GetBinContent(i+1);
        if (this_bin<1e-20) this_bin=0;
        if (next_bin<1e-20) next_bin=0;
        table_array[i+13][p] = this_bin;
        if (i==4)table_array[i+13][p] = this_bin+next_bin;
      }

      //Res 3b, low deltaR max
      table_array[18][p] = dem_hists[4]->Integral(1,dem_hists[4]->GetNbinsX()+1);
      for (int i=1; i<5; i++){
        float this_bin = dem_hists[4]->GetBinContent(i);
        float next_bin = dem_hists[4]->GetBinContent(i+1);
        if (this_bin<1e-20) this_bin=0;
        if (next_bin<1e-20) next_bin=0;
        table_array[i+18][p] = this_bin;
        if (i==4)table_array[i+18][p] = this_bin+next_bin;
      }

      //Res 3b, high deltaR max
      table_array[23][p] = dem_hists[5]->Integral(1,dem_hists[5]->GetNbinsX()+1);
      for (int i=1; i<5; i++){
        float this_bin = dem_hists[5]->GetBinContent(i);
        float next_bin = dem_hists[5]->GetBinContent(i+1);
        if (this_bin<1e-20) this_bin=0;
        if (next_bin<1e-20) next_bin=0;
        table_array[i+23][p] = this_bin;
        if (i==4)table_array[i+23][p] = this_bin+next_bin;
      }
    }
    delete tree;
    color_counter++;
  } //end loop through processes

  //For Rishi's MET stack plots
  if (save_stacks){
    TCanvas *c1 = new TCanvas("c1","stacked hists, boost",10,10,1000,800);
    stack_Boost->Draw("hist");
    c1->SetLogy();
    legend->Draw("same");
    c1->Modified();
    c1->SaveAs(directory+"Boost/stack_MET.root");

    TCanvas *c2 = new TCanvas("c2","stacked hists, DiAK8",10,10,1000,800);
    stack_DiAK8->Draw("hist");
    c2->SetLogy();
    legend->Draw("same");
    c2->Modified();
    c2->SaveAs(directory+"DiAK8/stack_MET.root");

    TCanvas *c3 = new TCanvas("c3","stacked hists, 4b resolved, lowDR",10,10,1000,800);
    stack_Res4b_lowR->Draw("hist");
    c3->SetLogy();
    legend->Draw("same");
    c3->Modified();
    c3->SaveAs(directory+"Res4b_lowR/stack_MET.root");

    TCanvas *c4 = new TCanvas("c4","stacked hists, 4b resolved, highDR",10,10,1000,800);
    stack_Res4b_highR->Draw("hist");
    c4->SetLogy();
    legend->Draw("same");
    c4->Modified();
    c4->SaveAs(directory+"Res4b_highR/stack_MET.root");

    TCanvas *c5 = new TCanvas("c5","stacked hists, 3b resolved, lowDR",10,10,1000,800);
    stack_Res3b_lowR->Draw("hist");
    c5->SetLogy();
    legend->Draw("same");
    c5->Modified();
    c5->SaveAs(directory+"Res3b_lowR/stack_MET.root");

    TCanvas *c6 = new TCanvas("c6","stacked hists, 3b resolved, highDR",10,10,1000,800);
    stack_Res3b_highR->Draw("hist");
    c6->SetLogy();
    legend->Draw("same");
    c6->Modified();
    c6->SaveAs(directory+"Res3b_highR/stack_MET.root");
  }


  if (make_table){   //Print table to screen by looping through the values of table_arry[row][column=process]
    ofstream myfile;
    myfile.open("FullComboYields_Test.txt");
    myfile<<endl<<endl;

    myfile<<"\\hline"<<endl;
    myfile<<"Case & TT & Zjets & WJets & QCD & TChiHH700 & TChiHH1000 \\\\ \\hline"<<endl;
    // myfile<<"Case & TT & Zjets & WJets & QCD  \\\\ \\hline"<<endl;
    // myfile<<"Case & TChiHH225 & TChiHH400 & TChiHH700 & TChiHH1000 \\\\ \\hline"<<endl;


    myfile<<"Res, 4b, lowDR";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[8][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 150-200";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[9][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 200-300";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[10][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 300-500";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[11][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: $>$500";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[12][i];
    myfile<<"\\\\"<<endl;
    myfile<<"\\hline"<<endl;

    myfile<<"Res, 4b, highDR";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[13][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 150-200";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[14][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 200-300";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[15][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 300-500";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[16][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: $>$500";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[17][i];
    myfile<<"\\\\"<<endl;
    myfile<<"\\hline"<<endl;

    myfile<<"Res, 3b, lowDR";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[18][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 150-200";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[19][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 200-300";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[20][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 300-500";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[21][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: $>$500";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[22][i];
    myfile<<"\\\\"<<endl;
    myfile<<"\\hline"<<endl;

    myfile<<"Res, 3b, highDR";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[23][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 150-200";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[24][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 200-300";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[25][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 300-500";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[26][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: $>$500";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[27][i];
    myfile<<"\\\\"<<endl;
    myfile<<"\\hline"<<endl;


    myfile<<"2 Boosted H";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[0][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 300-500";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[1][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 500-700";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[2][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: $>$700";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[3][i];
    myfile<<"\\\\"<<endl;
    myfile<<"\\hline"<<endl;

    myfile<<"1 Boosted, 1 AK8";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[4][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 300-500";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[5][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: 500-700";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[6][i];
    myfile<<"\\\\"<<endl;
    myfile<<"MET: $>$700";
    for(int i =0; i<6; i++) myfile<<" & "<<table_array[7][i];
    myfile<<"\\\\"<<endl;
    myfile<<"\\hline"<<endl;
    myfile<<endl<<endl;

    myfile.close();
  }
  inf->Close(); //is this needed?
} // MakeTable

void runMETStacks(TString path="./"){
  if(!TString(path).EndsWith("/")) path+="/";
  TFile* f_Boost_Bkg = TFile::Open(path+"Boost/Bkg.root","READ");
  TFile* f_Boost_TChiHH = TFile::Open(path+"Boost/TChiHH.root","READ");
  TFile* f_DiAK8_Bkg = TFile::Open(path+"DiAK8/Bkg.root","READ");
  TFile* f_DiAK8_TChiHH = TFile::Open(path+"DiAK8/TChiHH.root","READ");

  TFile* f_Res4bLowR_Bkg = TFile::Open(path+"Res4b_lowR/Bkg.root","READ");
  TFile* f_Res4bLowR_TChiHH = TFile::Open(path+"Res4b_lowR/TChiHH.root","READ");
  TFile* f_Res4bHighR_Bkg = TFile::Open(path+"Res4b_highR/Bkg.root","READ");
  TFile* f_Res4bHighR_TChiHH = TFile::Open(path+"Res4b_highR/TChiHH.root","READ");

  TFile* f_Res3bLowR_Bkg = TFile::Open(path+"Res3b_lowR/Bkg.root","READ");
  TFile* f_Res3bLowR_TChiHH = TFile::Open(path+"Res3b_lowR/TChiHH.root","READ");
  TFile* f_Res3bHighR_Bkg = TFile::Open(path+"Res3b_highR/Bkg.root","READ");
  TFile* f_Res3bHighR_TChiHH = TFile::Open(path+"Res3b_highR/TChiHH.root","READ");

  TH1F* h_Boost_QCD = (TH1F*)f_Boost_Bkg->Get("QCD_MET"); TH1F* h_Boost_TT = (TH1F*)f_Boost_Bkg->Get("TT_MET");
  TH1F* h_Boost_WJets = (TH1F*)f_Boost_Bkg->Get("WJets_MET"); TH1F* h_Boost_ZJets = (TH1F*)f_Boost_Bkg->Get("ZJets_MET");
  TH1F* h_DiAK8_QCD = (TH1F*)f_DiAK8_Bkg->Get("QCD_MET"); TH1F* h_DiAK8_TT = (TH1F*)f_DiAK8_Bkg->Get("TT_MET");
  TH1F* h_DiAK8_WJets = (TH1F*)f_DiAK8_Bkg->Get("WJets_MET"); TH1F* h_DiAK8_ZJets = (TH1F*)f_DiAK8_Bkg->Get("ZJets_MET");

  TH1F* h_Res4bLowR_QCD = (TH1F*)f_Res4bLowR_Bkg->Get("QCD_MET"); TH1F* h_Res4bLowR_TT = (TH1F*)f_Res4bLowR_Bkg->Get("TT_MET");
  TH1F* h_Res4bLowR_WJets = (TH1F*)f_Res4bLowR_Bkg->Get("WJets_MET"); TH1F* h_Res4bLowR_ZJets = (TH1F*)f_Res4bLowR_Bkg->Get("ZJets_MET");
  TH1F* h_Res4bHighR_QCD = (TH1F*)f_Res4bHighR_Bkg->Get("QCD_MET"); TH1F* h_Res4bHighR_TT = (TH1F*)f_Res4bHighR_Bkg->Get("TT_MET");
  TH1F* h_Res4bHighR_WJets = (TH1F*)f_Res4bHighR_Bkg->Get("WJets_MET"); TH1F* h_Res4bHighR_ZJets = (TH1F*)f_Res4bHighR_Bkg->Get("ZJets_MET");

  TH1F* h_Res3bLowR_QCD = (TH1F*)f_Res3bLowR_Bkg->Get("QCD_MET"); TH1F* h_Res3bLowR_TT = (TH1F*)f_Res3bLowR_Bkg->Get("TT_MET");
  TH1F* h_Res3bLowR_WJets = (TH1F*)f_Res3bLowR_Bkg->Get("WJets_MET"); TH1F* h_Res3bLowR_ZJets = (TH1F*)f_Res3bLowR_Bkg->Get("ZJets_MET");
  TH1F* h_Res3bHighR_QCD = (TH1F*)f_Res3bHighR_Bkg->Get("QCD_MET"); TH1F* h_Res3bHighR_TT = (TH1F*)f_Res3bHighR_Bkg->Get("TT_MET");
  TH1F* h_Res3bHighR_WJets = (TH1F*)f_Res3bHighR_Bkg->Get("WJets_MET"); TH1F* h_Res3bHighR_ZJets = (TH1F*)f_Res3bHighR_Bkg->Get("ZJets_MET");

  TH1F* h_Boost_TChiHH700 = (TH1F*)f_Boost_TChiHH->Get("TChiHH700_MET"); TH1F* h_Boost_TChiHH1000 = (TH1F*)f_Boost_TChiHH->Get("TChiHH1000_MET");
  TH1F* h_DiAK8_TChiHH700 = (TH1F*)f_DiAK8_TChiHH->Get("TChiHH700_MET"); TH1F* h_DiAK8_TChiHH1000 = (TH1F*)f_DiAK8_TChiHH->Get("TChiHH1000_MET");
  TH1F* h_Res4bLowR_TChiHH700 = (TH1F*)f_Res4bLowR_TChiHH->Get("TChiHH700_MET"); TH1F* h_Res4bLowR_TChiHH1000 = (TH1F*)f_Res4bLowR_TChiHH->Get("TChiHH1000_MET");
  TH1F* h_Res4bHighR_TChiHH700 = (TH1F*)f_Res4bHighR_TChiHH->Get("TChiHH700_MET"); TH1F* h_Res4bHighR_TChiHH1000 = (TH1F*)f_Res4bHighR_TChiHH->Get("TChiHH1000_MET");
  TH1F* h_Res3bLowR_TChiHH700 = (TH1F*)f_Res3bLowR_TChiHH->Get("TChiHH700_MET"); TH1F* h_Res3bLowR_TChiHH1000 = (TH1F*)f_Res3bLowR_TChiHH->Get("TChiHH1000_MET");
  TH1F* h_Res3bHighR_TChiHH700 = (TH1F*)f_Res3bHighR_TChiHH->Get("TChiHH700_MET"); TH1F* h_Res3bHighR_TChiHH1000 = (TH1F*)f_Res3bHighR_TChiHH->Get("TChiHH1000_MET");

  MET_stack(h_Res4bLowR_QCD, h_Res4bLowR_TT, h_Res4bLowR_WJets, h_Res4bLowR_ZJets, h_Res4bLowR_TChiHH700, h_Res4bLowR_TChiHH1000, "Res4bLowR");
  MET_stack(h_Res4bHighR_QCD, h_Res4bHighR_TT, h_Res4bHighR_WJets, h_Res4bHighR_ZJets, h_Res4bHighR_TChiHH700, h_Res4bHighR_TChiHH1000, "Res4bHighR");
  MET_stack(h_Res3bLowR_QCD, h_Res3bLowR_TT, h_Res3bLowR_WJets, h_Res3bLowR_ZJets, h_Res3bLowR_TChiHH700, h_Res3bLowR_TChiHH1000, "Res3bLowR");
  MET_stack(h_Res3bHighR_QCD, h_Res3bHighR_TT, h_Res3bHighR_WJets, h_Res3bHighR_ZJets, h_Res3bHighR_TChiHH700, h_Res3bHighR_TChiHH1000, "Res3bHighR");
  MET_stack(h_Boost_QCD, h_Boost_TT, h_Boost_WJets, h_Boost_ZJets, h_Boost_TChiHH700, h_Boost_TChiHH1000, "Boost");
  MET_stack(h_DiAK8_QCD, h_DiAK8_TT, h_DiAK8_WJets, h_DiAK8_ZJets, h_DiAK8_TChiHH700, h_DiAK8_TChiHH1000, "DiAK8");
}

void MET_stack(TH1F* h_QCD,TH1F* h_TT,TH1F* h_WJets,TH1F* h_ZJets, TH1F* h_TChiHH700,TH1F* h_TChiHH1000, TString which){
  h_QCD->SetFillColor(kGray); h_QCD->SetMarkerStyle(21); h_QCD->SetMarkerColor(kGray);
  h_TT->SetFillColor(kCyan); h_TT->SetMarkerStyle(21); h_TT->SetMarkerColor(kCyan);
  h_WJets->SetFillColor(kBlue);h_WJets->SetMarkerStyle(21); h_WJets->SetMarkerColor(kBlue);
  h_ZJets->SetFillColor(kGreen+2);h_ZJets->SetMarkerStyle(21); h_ZJets->SetMarkerColor(kGreen+2);

  h_TChiHH700->SetMarkerStyle(20); h_TChiHH700->SetMarkerColor(kRed); h_TChiHH700->SetLineColor(kRed);
  h_TChiHH700->SetLineWidth(2);
  h_TChiHH1000->SetMarkerStyle(20); h_TChiHH1000->SetMarkerColor(kBlack); h_TChiHH1000->SetLineColor(kBlack);
  h_TChiHH1000->SetLineWidth(2);

  THStack * METStack = new THStack("hs","");
  METStack->Add(h_ZJets);
  METStack->Add(h_WJets);
  METStack->Add(h_TT);
  METStack->Add(h_QCD);

  double W = 600;    double H = 600;
  double T = 0.08*H; double B = 0.12*H;
  double L = 0.08*W; double R = 0.08*W;

  TString graphName = "MET"+which;
  TCanvas * can_h = new TCanvas(graphName,graphName, 50, 50, W, H);
  can_h->SetFillColor(0); can_h->SetBorderMode(0);
  can_h->SetFrameFillStyle(0); can_h->SetFrameBorderMode(0);
  can_h->SetLeftMargin( L/W ); can_h->SetRightMargin( R/W );
  can_h->SetTopMargin( T/H ); can_h->SetBottomMargin( B/H );
  can_h->SetTickx(0); can_h->SetTicky(0);
  can_h->SetLogy();

  METStack->Draw("hist");
  METStack->SetTitle(which);
  METStack->GetYaxis()->SetTitle("Events");
  METStack->GetXaxis()->SetTitle("MET [GeV]");

  if (which=="Res4bLowR" || which=="Res3bLowR"){
    METStack->SetMaximum(30);
  }
  else if (which=="DiAK8"){
    METStack->SetMinimum(0.003);
  }

  h_TChiHH700->Draw("same");
  h_TChiHH1000->Draw("same");

  TLegend* legend = new TLegend(0.35,0.8,0.9,0.9) ;
  legend->AddEntry(h_QCD, "QCD", "f");
  legend->AddEntry(h_TT, "TT", "f");
  legend->AddEntry(h_WJets, "WJets", "f");
  legend->AddEntry(h_ZJets, "ZJets", "f");
  legend->AddEntry(h_TChiHH700, "TChiHH700", "lp");
  legend->AddEntry(h_TChiHH1000, "TChiHH1000", "lp");
  legend->SetNColumns(3);
  legend->SetBorderSize(0);

  legend->Draw("same");
  can_h->Update();
  can_h->Modified();

  TString savename = directory+"MET_"+which+".pdf";
  can_h->SaveAs(savename);
}
