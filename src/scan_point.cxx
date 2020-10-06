#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <limits>
#include <getopt.h>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TDirectory.h"

using namespace std;

std::vector<std::string> split(const std::string& s, char delimiter);
void higgsino2DCrossSection(int hig_mass, double &xsec, double &xsec_unc);

string in_dir = "/uscms_data/d3/emacdona/WorkingArea/CombinedHiggs/V18_git/CMSSW_8_1_0/src/boostedHiggsPlusMET/src/";
string out_dir = in_dir;
string model = "N1N2";

void scan_point() {
  ifstream file(in_dir+"higgsino2DFileNames.txt");
  string line;
  TString filename;
  string txtname(out_dir+"/limits_"+model+".txt");
  ofstream txtfile(txtname);

  while(std::getline(file, line)) {
    std::vector<std::string> x = split(line, '_');
    string hino_mass = x[5];
    string LSP_mass = x[6];
    std::cout<<"hino: "<<hino_mass<<", LSP: "<<LSP_mass<<std::endl;
    int hino_mass_int = std::stoi(hino_mass);
    filename = in_dir+"/higgsCombineTChiHH"+hino_mass+"_LSP"+LSP_mass+"_BothBoostedHMerge.AsymptoticLimits.mH120.root";
    double xsec, xsec_unc;
    if (model=="N1N2") higgsino2DCrossSection(hino_mass_int, xsec, xsec_unc);

    TFile limits_file(filename, "read");
    if (!limits_file.IsOpen())
      std::cout<<"Could not open limits file "<<filename<<std::endl;

    TTree *tree = static_cast<TTree*>(limits_file.Get("limit"));
    if (tree == nullptr)
      std::cout<<"Could not get limits tree"<<std::endl;

    double limit;
    tree->SetBranchAddress("limit", &limit);
    int num_entries = tree->GetEntries();
    if (num_entries != 6)
      std::cout<<"Expected 6 tree entries. Saw "<<num_entries<<std::endl;

    tree->GetEntry(0); double exp_2down = limit;
    tree->GetEntry(1); double exp_down = limit;
    tree->GetEntry(2); double exp = limit;
    tree->GetEntry(3); double exp_up = limit;
    tree->GetEntry(4); double exp_2up = limit;
    tree->GetEntry(5); double obs = limit;

    limits_file.Close();

    txtfile << setprecision(numeric_limits<double>::max_digits10)
      << ' ' << hino_mass
      << ' ' << LSP_mass
      << ' ' << xsec
      << ' ' << xsec_unc
      << ' ' << obs
      << ' ' << exp
      << ' ' << exp_up
      << ' ' << exp_down
      << ' ' << exp_2up
      << ' ' << exp_2down;
    txtfile << endl;
  } // end while loop
  txtfile.close();
} // end scan point

std::vector<std::string> split(const std::string& s, char delimiter) {
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter)) { tokens.push_back(token);}
   return tokens;
}

void higgsino2DCrossSection(int hig_mass, double &xsec, double &xsec_unc) {
  if (hig_mass ==127) { xsec = .5824*.5824*1.44725; xsec_unc = 0.0395277; return;}
  else if (hig_mass ==150) { xsec = .5824*.5824*0.71514; xsec_unc = 0.0421496; return;}
  else if (hig_mass ==175) { xsec = .5824*.5824*0.419059; xsec_unc = 0.0453279; return;}
  else if (hig_mass ==200) { xsec = .5824*.5824*0.244213; xsec_unc = 0.047925; return;}
  else if (hig_mass ==225) { xsec = .5824*.5824*0.156286; xsec_unc = 0.0502876; return;}
  else if (hig_mass ==250) { xsec = .5824*.5824*0.104252; xsec_unc = 0.0526169; return;}
  else if (hig_mass ==275) { xsec = .5824*.5824*0.0719125; xsec_unc = 0.0549666; return;}
  else if (hig_mass ==300) { xsec = .5824*.5824*0.0509994; xsec_unc = 0.0572762; return;}
  else if (hig_mass ==325) { xsec = .5824*.5824*0.0369715; xsec_unc = 0.0590317; return;}
  else if (hig_mass ==350) { xsec = .5824*.5824*0.0273286; xsec_unc = 0.0607766; return;}
  else if (hig_mass ==375) { xsec = .5824*.5824*0.0205429; xsec_unc = 0.0625031; return;}
  else if (hig_mass ==400) { xsec = .5824*.5824*0.0156691; xsec_unc = 0.0642085; return;}
  else if (hig_mass ==425) { xsec = .5824*.5824*0.0120965; xsec_unc = 0.0657801; return;}
  else if (hig_mass ==450) { xsec = .5824*.5824*0.00944017; xsec_unc = 0.0674544; return;}
  else if (hig_mass ==475) { xsec = .5824*.5824*0.00743587; xsec_unc = 0.0686033; return;}
  else if (hig_mass ==500) { xsec = .5824*.5824*0.00590757; xsec_unc = 0.0699909; return;}
  else if (hig_mass ==525) { xsec = .5824*.5824*0.00469101; xsec_unc = 0.0713704; return;}
  else if (hig_mass ==550) { xsec = .5824*.5824*0.0038167; xsec_unc = 0.0722834; return;}
  else if (hig_mass ==575) { xsec = .5824*.5824*0.003073; xsec_unc = 0.0739957; return;}
  else if (hig_mass ==600) { xsec = .5824*.5824*0.00253015; xsec_unc = 0.0754291; return;}
  else if (hig_mass ==625) { xsec = .5824*.5824*0.00206136; xsec_unc = 0.0763466; return;}
  else if (hig_mass ==650) { xsec = .5824*.5824*0.00171418; xsec_unc = 0.0775695; return;}
  else if (hig_mass ==675) { xsec = .5824*.5824*0.00140934; xsec_unc = 0.0783375; return;}
  else if (hig_mass ==700) { xsec = .5824*.5824*0.00118113; xsec_unc = 0.0796388; return;}
  else if (hig_mass ==725) { xsec = .5824*.5824*0.000979349; xsec_unc = 0.0809883; return;}
  else if (hig_mass ==750) { xsec = .5824*.5824*0.000826366; xsec_unc = 0.081879; return;}
  else if (hig_mass ==775) { xsec = .5824*.5824*0.000690208; xsec_unc = 0.0842049; return;}
  else if (hig_mass ==800) { xsec = .5824*.5824*0.000586211; xsec_unc = 0.0862527; return;}
  else if (hig_mass ==825) { xsec = .5824*.5824*0.00049277; xsec_unc = 0.0864444; return;}
  else if (hig_mass ==850) { xsec = .5824*.5824*0.000420556; xsec_unc = 0.085742; return;}
  else if (hig_mass ==875) { xsec = .5824*.5824*0.000358734; xsec_unc = 0.0889174; return;}
  else if (hig_mass ==900) { xsec = .5824*.5824*0.000305935; xsec_unc = 0.0912439; return;}
  else if (hig_mass ==925) { xsec = .5824*.5824*0.000260948; xsec_unc = 0.091372; return;}
  else if (hig_mass ==950) { xsec = .5824*.5824*0.00022285; xsec_unc = 0.0919538; return;}
  else if (hig_mass ==975) { xsec = .5824*.5824*0.000189681; xsec_unc = 0.0938108; return;}
  else if (hig_mass ==1000) { xsec = .5824*.5824*0.00016428; xsec_unc = 0.0954285; return;}
  else if (hig_mass ==1025) { xsec = .5824*.5824*0.000142206; xsec_unc = 0.0957231; return;}
  else if (hig_mass ==1050) { xsec = .5824*.5824*0.000120971; xsec_unc = 0.0968997; return;}
  else if (hig_mass ==1075) { xsec = .5824*.5824*0.000105301; xsec_unc = 0.0979041; return;}
  else if (hig_mass ==1100) { xsec = .5824*.5824*9.12469e-05; xsec_unc = 0.0964142; return;}
  else if (hig_mass ==1125) { xsec = .5824*.5824*7.9765e-05; xsec_unc = 0.099902; return;}
  else if (hig_mass ==1150) { xsec = .5824*.5824*6.78234e-05; xsec_unc = 0.101061; return;}
  else if (hig_mass ==1175) { xsec = .5824*.5824*5.9016e-05; xsec_unc = 0.102051; return;}
  else if (hig_mass ==1200) { xsec = .5824*.5824*5.16263e-05; xsec_unc = 0.102499; return;}
  else if (hig_mass ==1225) { xsec = .5824*.5824*4.5147e-05; xsec_unc = 0.10403; return;}
  else if (hig_mass ==1250) { xsec = .5824*.5824*3.88343e-05; xsec_unc = 0.105206; return;}
  else if (hig_mass ==1275) { xsec = .5824*.5824*3.41304e-05; xsec_unc = 0.10619; return;}
  else if (hig_mass ==1300) { xsec = .5824*.5824*2.99353e-05; xsec_unc = 0.10783; return;}
  else if (hig_mass ==1325) { xsec = .5824*.5824*2.63637e-05; xsec_unc = 0.108024; return;}
  else if (hig_mass ==1350) { xsec = .5824*.5824*2.26779e-05; xsec_unc = 0.109016; return;}
  else if (hig_mass ==1375) { xsec = .5824*.5824*1.99318e-05; xsec_unc = 0.109822; return;}
  else if (hig_mass ==1400) { xsec = .5824*.5824*1.75031e-05; xsec_unc = 0.111631; return;}
  else if (hig_mass ==1425) { xsec = .5824*.5824*1.53974e-05; xsec_unc = 0.111417; return;}
  else if (hig_mass ==1455) { xsec = .5824*.5824*1.3245e-05; xsec_unc = 0.112313; return;}
  else if (hig_mass ==1475) { xsec = .5824*.5824*1.16416e-05; xsec_unc = 0.113058; return;}
  else { xsec = 0; xsec_unc = 0;}
}
