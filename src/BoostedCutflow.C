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
#include "tdrstyle.C"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

//#include "INeedThis.C"
#include "BoostedLoop.C"

using namespace std;

void BoostedCutflow(string path="./") {
  setTDRStyle();
  if(!TString(path).EndsWith("/")) path+="/";
  TFile* inf = TFile::Open((path+"BoostedCutflow_V12.root").c_str(),"READ");

  bool save_files = true;
  bool make_table = true;

  float table_array[7][6] = {0};


  //vector<string> processes = {"ZJets","WJets","TT","QCD","T5qqqqZH1700","T5qqqqZH2100"};
  vector<string> processes = {"ZJets","WJets","TT","QCD"};

  string var = "MET";
  string direct = "";

  vector<float> cutflow_vector = {0};

  for(unsigned int p = 0; p< processes.size(); p++) {
    //From Rishi's repo from old times double lumi=35862.824;
    double lumi = 35862.824;
    double EventWeight_T5HH = lumi*4.0;
    double EventWeight = lumi;

    cout<<"Beginning "<< processes[p].c_str()<<endl;
    TTree *tree;
    inf->GetObject(processes[p].c_str(),tree);
    BoostedLoop ThisBoost;
    ThisBoost.Init(tree);

    if (processes[p].find("T5qqqq")!=string::npos) cutflow_vector = ThisBoost.Loop(EventWeight_T5HH, processes[p].c_str(), save_files);
    else cutflow_vector = ThisBoost.Loop(EventWeight, processes[p].c_str(), save_files);


    cout<<"Cutflow for "<<processes[p].c_str()<<": ";
    for (int i=0; i<7; i++){
      table_array[i][p] = cutflow_vector[i];
      cout<<cutflow_vector[i]<<", ";
    }
    cout<<endl;


    delete tree;
  } //end loop through processes

  float total_bkg[7] = {0};
  for(int i =0; i<=6; i++){ //loop through cutflow bools
    float total = 0;
    for(int j=0; j<=3; j++){ //loop through processes
      float this_value = table_array[i][j];
      total+=this_value;
    }

    total_bkg[i] = total;
  }


  if (make_table){
    //Print table to screen by looping through the values of table_arry[row][column=process]
    cout<<endl<<endl;
    cout<<"  _____________________________________________________________________________________________________________"<<endl;
    cout<<" |                              | Total Bkg. |  ZJets   |  WJets   |    TT    |   QCD    | T5HH1700 | T5HH2100 |"<<endl;
    cout<<" |------------------------------+------------+----------+----------+----------+----------+----------+----------|"<<endl;
    cout<<" |Baseline SUSY Hadronic Skim   |";
    cout<<setw(12)<<setprecision(6)<<total_bkg[0]<<"|";
    for(int i =0; i<=5; i++) cout<<setw(10)<<setprecision(6)<<table_array[0][i]<<"|";
    cout<<endl;
    cout<<" |------------------------------+------------+----------+----------+----------+----------+----------+----------|"<<endl;
    cout<<" |HT>600 GeV, MET>300 GeV       |";
    cout<<setw(12)<<setprecision(6)<<total_bkg[1]<<"|";
    for(int i =0; i<=5; i++) cout<<setw(10)<<setprecision(5)<<table_array[1][i]<<"|";
    cout<<endl;
    cout<<" |------------------------------+------------+----------+----------+----------+----------+----------+----------|"<<endl;
    cout<<" |2+ AK8 Jets                   |";
    cout<<setw(12)<<setprecision(6)<<total_bkg[2]<<"|";
    for(int i =0; i<=5; i++) cout<<setw(10)<<setprecision(5)<<table_array[2][i]<<"|";
    cout<<endl;
    cout<<" |------------------------------+------------+----------+----------+----------+----------+----------+----------|"<<endl;
    cout<<" |AK8 Jet pT > 300 GeV          |";
    cout<<setw(12)<<setprecision(6)<<total_bkg[3]<<"|";
    for(int i =0; i<=5; i++) cout<<setw(10)<<setprecision(5)<<table_array[3][i]<<"|";
    cout<<endl;
    cout<<" |------------------------------+------------+----------+----------+----------+----------+----------+----------|"<<endl;
    cout<<" |AK8 Jet Mass in [50,250] GeV  |";
    cout<<setw(12)<<setprecision(6)<<total_bkg[4]<<"|";
    for(int i =0; i<=5; i++) cout<<setw(10)<<setprecision(5)<<table_array[4][i]<<"|";
    cout<<endl;
    cout<<" |------------------------------+------------+----------+----------+----------+----------+----------+----------|"<<endl;
    cout<<" |AK8 Jet Double-b              |";
    cout<<setw(12)<<setprecision(6)<<total_bkg[5]<<"|";
    for(int i =0; i<=5; i++) cout<<setw(10)<<setprecision(5)<<table_array[5][i]<<"|";
    cout<<endl;
    cout<<" |------------------------------+------------+----------+----------+----------+----------+----------+----------|"<<endl;
    cout<<" |AK8 Jet Mass in [85,135] GeV  |";
    cout<<setw(12)<<setprecision(6)<<total_bkg[6]<<"|";
    for(int i =0; i<=5; i++) cout<<setw(10)<<setprecision(5)<<table_array[6][i]<<"|";
    cout<<endl;
    cout<<" ---------------------------------------------------------------------------------------------------------------"<<endl;

    cout<<endl<<endl;
  }
  inf->Close();
} // Boosted Cutflow
