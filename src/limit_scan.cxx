#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <unistd.h>
#include <getopt.h>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"

using namespace std;

int num_smooth_ = 0; // Number of times to smooth TH2D
string in_dir = "/uscms_data/d3/emacdona/WorkingArea/CombinedHiggs/V18_git/CMSSW_8_1_0/src/boostedHiggsPlusMET/src/";
string out_dir = in_dir;
string tag = "";
string model_ = "N1N2";

//pre-define all functions
// void ReadPoints(vector<double> &vmx,vector<double> &vmy,vector<double> &vxsec,vector<double> &vobs,vector<double> &vobsup,vector<double> &vobsdown, vector<double> &vexp,vector<double> &vup,vector<double> &vdown);
void ReadPoints(vector<double> &vmx,vector<double> &vmy,vector<double> &vxsec,vector<double> &vobs,vector<double> &vobsglu, vector<double> &vobsup,vector<double> &vobsdown, vector<double> &vobsup2,vector<double> &vobsdown2,vector<double> &vexp,vector<double> &vup,vector<double> &vdown,vector<double> &vup2,vector<double> &vdown2);
TH2D MakeObservedSignificancePlot(vector<double> vmx,vector<double> vmy,vector<double> vobs);
TH2D MakeExpectedSignificancePlot(vector<double> vmx,vector<double> vmy,vector<double> vobs);
// void MakeLimitPlot(vector<double> vmx,vector<double> vmy,vector<double> vlim,vector<double> vobs,vector<double> vobsup,vector<double> vobsdown,vector<double> vexp,vector<double> vup,vector<double> vdown);
void MakeLimitPlot(vector<double> vmx,vector<double> vmy,vector<double> vlim,vector<double> vobs,vector<double> vobsup,vector<double> vobsdown,vector<double> vobsup2,vector<double> vobsdown2,vector<double> vexp,vector<double> vup,vector<double> vdown,vector<double> vup2,vector<double> vdown2);
int GetNumBins(const vector<double> &pts, double width);
void GetParticleNames(string &xparticle, string &yparticle);
TLatex GetModelLabel(double x, double y);
void Style(TGraph *g, int color, int style, float width);
TGraph DrawContours(TGraph2D &g2, int color, int style, double width,int n_smooth, double val);
TGraph* joinGraphs(TGraph *graph1, TGraph *graph2);
void FixGraph(TGraph &graph);
void ReverseGraph(TGraph &graph);
void SetupColors();
void SetupSignedColors();
std::vector<std::string> split(const std::string& s, char delimiter);
void higgsinoCrossSection(int hig_mass, double &xsec, double &xsec_unc);
void gluinoCrossSection(int glu_mass, double &xsec, double &xsec_unc); //this is for gluino mass

void limit_scan() {
  // vector<double> vmx, vmy, vxsec, vobs, vobsup, vobsdown, vexp, vup, vdown;
  // ReadPoints(vmx, vmy, vxsec, vobs, vobsup, vobsdown, vexp, vup, vdown);
  vector<double> vmx, vmy, vxsec, vobs, vobsglu, vobsup, vobsdown, vobsup2, vobsdown2, vexp, vup, vdown, vup2, vdown2;
  ReadPoints(vmx, vmy, vxsec, vobs, vobsglu, vobsup, vobsdown, vobsup2, vobsdown2, vexp, vup, vdown, vup2, vdown2);

  vector<double> vlim(vxsec.size());
  for (size_t i = 0; i < vxsec.size(); ++i) {
    vlim.at(i) = vxsec.at(i) * vobs.at(i);
  }

  // TH2D hsigobs = MakeObservedSignificancePlot(vmx, vmy, vsigobs);
  // TH2D hsigexp = MakeExpectedSignificancePlot(vmx, vmy, vsigexp);

  MakeLimitPlot(vmx, vmy, vlim, vobsglu, vobsup, vobsdown, vobsup2, vobsdown2, vexp, vup, vdown, vup2, vdown2);
}

void ReadPoints(vector<double> &vmx,
                vector<double> &vmy,
                vector<double> &vxsec,
                vector<double> &vobs,
                vector<double> &vobsglu, //currently this is gluino
                vector<double> &vobsup, //currently this is gluino
                vector<double> &vobsdown, //currently this is gluino
                vector<double> &vobsup2, //currently this is gluino
                vector<double> &vobsdown2, //currently this is gluino
                vector<double> &vexp, //currently this is higgsino with bigger xsec
                vector<double> &vup, //currently this is higgsino with bigger xsec
                vector<double> &vdown, //currently this is higgsino with bigger xsec
                vector<double> &vup2, //currently this is higgsino with bigger xsec
                vector<double> &vdown2 //currently this is higgsino with bigger xsec
              ) {

  TString filename = out_dir+"limits_N1N2.txt";
  ifstream infile(filename);
  string line;

  while(getline(infile, line)) {
    istringstream iss(line);
    double pmx, pmy, pxsec, pxsecunc, pobs, pexp, pup, pdown, pup2, pdown2;
    iss >> pmx >> pmy >> pxsec >> pxsecunc >> pobs >> pexp >> pup >> pdown >> pup2 >> pdown2;
    vmx.push_back(pmx);
    vmy.push_back(pmy);
    vxsec.push_back(pxsec);

    // vobsup.push_back(pobs/(1+pxsecunc));
    // vobsdown.push_back(pobs/(1-pxsecunc));
    double glu_xsec; double glu_xsecerr;
    gluinoCrossSection(pmx, glu_xsec, glu_xsecerr);
    vobs.push_back(pobs);
    vobsglu.push_back(pexp*pxsec/glu_xsec);
    vobsup.push_back(pup*pxsec/glu_xsec);
    vobsdown.push_back(pdown*pxsec/glu_xsec);
    vobsup2.push_back(pup2*pxsec/glu_xsec);
    vobsdown2.push_back(pdown2*pxsec/glu_xsec);

    double this_xsec; double this_xsecerr;
    higgsinoCrossSection(pmx, this_xsec, this_xsecerr);
    vexp.push_back(pexp*pxsec/this_xsec);
    vup.push_back(pup*pxsec/this_xsec);
    vdown.push_back(pdown*pxsec/this_xsec);
    vup2.push_back(pup2*pxsec/this_xsec);
    vdown2.push_back(pdown2*pxsec/this_xsec);
  }
  infile.close();

  if (vmx.size() <= 2) {
    std::cout<<"Need at least 3 models to draw scan"<<std::endl;
    return;
  }
  if (vmx.size() != vmy.size()
     || vmx.size() != vxsec.size()
     || vmx.size() != vobs.size()
     || vmx.size() != vobsup.size()
     || vmx.size() != vobsdown.size()
     || vmx.size() != vexp.size()
     || vmx.size() != vup.size()
     || vmx.size() != vdown.size() ){
       std::cout<<"Error parsing text file. Model point not fully specified\n";
     }
}

TH2D MakeObservedSignificancePlot(vector<double> vmx,
                                  vector<double> vmy,
                                  vector<double> vobs) {
  SetupSignedColors();

  string xparticle, yparticle;
  GetParticleNames(xparticle, yparticle);
  string title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV];Observed Significance";
  TGraph2D g("", title.c_str(), vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));

  double the_max = 3.;
  g.SetMinimum(-the_max);
  g.SetMaximum(the_max);
  g.SetNpx((2600.-800.)/12.5);
  g.SetNpy(1600/12.5);
  g.GetHistogram()->SetTitle(title.c_str());
  g.GetHistogram()->SetTickLength(0., "Z");

  TCanvas c("can","can", 800, 900);
  c.cd();

  TLatex ltitle(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),
              "#font[62]{CMS}#scale[0.76]{#font[52]{ Supplementary}}");
  TLatex rtitle(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),
               "#scale[0.8]{137 fb^{-1} (13 TeV)}");
  ltitle.SetNDC();
  rtitle.SetNDC();
  ltitle.SetTextAlign(12);
  rtitle.SetTextAlign(32);
  TLatex model = GetModelLabel(c.GetLeftMargin()+0.03, 1.-c.GetTopMargin()-0.03);

  g.Draw("colz");

  vector<TLine> lines;
  for (double z = 0.; z < the_max; z+=0.5) {
    int style = min(5, static_cast<int>(2.*fabs(z))+1);
    double width = max(1., 3.-fabs(z));
    DrawContours(g, 1, style, width, 0, z);
    double x1 = 1.-c.GetRightMargin()+0.0047;
    double x2 = 1.-c.GetRightMargin()+0.05;
    double ybot = c.GetBottomMargin();
    double ytop = 1.-c.GetTopMargin();
    double zpos = ybot+(ytop-ybot)*(the_max+z)/(2.*the_max);
    lines.emplace_back(x1, zpos, x2, zpos);
    lines.back().SetLineColor(1);
    lines.back().SetLineStyle(style);
    lines.back().SetLineWidth(width);
    lines.back().SetNDC(true);
    if (z != 0.) {
      DrawContours(g, 1, style, width, 0, -z);
      double zneg = ybot+(ytop-ybot)*(the_max-z)/(2.*the_max);
      lines.emplace_back(x1, zneg, x2, zneg);
      lines.back().SetLineColor(1);
      lines.back().SetLineStyle(style);
      lines.back().SetLineWidth(width);
      lines.back().SetNDC(true);
    }
  }
  for (auto &l: lines) l.Draw("same");

  model.Draw("same");
  ltitle.Draw("same");
  rtitle.Draw("same");

  c.Print((model_+"_sigobs_"+tag+".pdf").c_str());
  c.Print((model_+"_sigobs_"+tag+".root").c_str());

  TH2D h = *g.GetHistogram();
  h.SetTitle("Observed Significance");
  return h;
}

TH2D MakeExpectedSignificancePlot(vector<double> vmx,
                                  vector<double> vmy,
                                  vector<double> vobs) {
  SetupColors();

  string xparticle, yparticle;
  GetParticleNames(xparticle, yparticle);
  string title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV]; Expected Significance";

  TGraph2D g("", title.c_str(), vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));

  double the_max = 0.;
  for (int i = 0; i < g.GetN(); ++i) {
    double z = g.GetZ()[i];
    if (z>the_max) the_max = z;
  }
  if (the_max > 6.) the_max = 6.;
  g.SetMinimum(0.);
  g.SetMaximum(the_max);

  g.SetNpx((2600.-800.)/12.5);
  g.SetNpy(1600/12.5);

  g.GetHistogram()->SetTitle(title.c_str());

  TCanvas c;
  c.cd();

  TLatex ltitle(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),
              "#font[62]{CMS}#scale[0.76]{#font[52]{ Supplementary}}");
  TLatex rtitle(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),
               "#scale[0.8]{137 fb^{-1} (13 TeV)}");
  ltitle.SetNDC();
  rtitle.SetNDC();
  ltitle.SetTextAlign(12);
  rtitle.SetTextAlign(32);
  TLatex model = GetModelLabel(c.GetLeftMargin()+0.03, 1.-c.GetTopMargin()-0.03);

  g.Draw("colz");

  for (double z = 0.; z < the_max; z+=1.) {
    int style = (z == 5. ? 1 : 2);
    double width = (z == 5. ? 2. : 1.);
    DrawContours(g, 1, style, width, 0, z);
  }

  model.Draw("same");
  ltitle.Draw("same");
  rtitle.Draw("same");

  c.Print((model_+"_sigexp_"+tag+".pdf").c_str());
  c.Print((model_+"_sigexp_"+tag+".root").c_str());

  TH2D h = *g.GetHistogram();
  h.SetTitle("Expected Significance");
  return h;
}

void MakeLimitPlot(vector<double> vmx,
                   vector<double> vmy,
                   vector<double> vlim,
                   vector<double> vobs,
                   vector<double> vobsup,
                   vector<double> vobsdown,
                   vector<double> vobsup2,
                   vector<double> vobsdown2,
                   vector<double> vexp,
                   vector<double> vup,
                   vector<double> vdown,
                   vector<double> vup2,
                   vector<double> vdown2
                 ) {
  SetupColors();

  string xparticle, yparticle;
  GetParticleNames(xparticle, yparticle);
  string title = ";m_{"+xparticle+"} (used as proxy for m_{#tilde{g}}) [GeV];m_{"+yparticle+"} [GeV];95% CL upper limit on cross section [pb]";

  TGraph2D glim("", title.c_str(), vlim.size(), &vmx.at(0), &vmy.at(0), &vlim.at(0));
  // TGraph2D gobs("", "Observed Limit", vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));
  // TGraph2D gobsup("", "Observed +1#sigma Limit", vobsup.size(), &vmx.at(0), &vmy.at(0), &vobsup.at(0));
  // TGraph2D gobsdown("", "Observed -1#sigma Limit", vobsdown.size(), &vmx.at(0), &vmy.at(0), &vobsdown.at(0));
  // TGraph2D gexp("", "Expected Limit", vexp.size(), &vmx.at(0), &vmy.at(0), &vexp.at(0));
  // TGraph2D gup("", "Expected +1#sigma Limit", vup.size(), &vmx.at(0), &vmy.at(0), &vup.at(0));
  // TGraph2D gdown("", "Expected -1#sigma Limit", vdown.size(), &vmx.at(0), &vmy.at(0), &vdown.at(0));


  TGraph2D gexpG("", "gexpG", vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));
  TGraph2D gexpGup("", "gexpGup", vobs.size(), &vmx.at(0), &vmy.at(0), &vobsup.at(0));
  TGraph2D gexpGdown("", "gexpGdown", vobs.size(), &vmx.at(0), &vmy.at(0), &vobsdown.at(0));
  TGraph2D gexpGup2("", "gexpGup2", vobs.size(), &vmx.at(0), &vmy.at(0), &vobsup2.at(0));
  TGraph2D gexpGdown2("", "gexpGdown2", vobs.size(), &vmx.at(0), &vmy.at(0), &vobsdown2.at(0));
  gexpG.SetMinimum(0.000001); gexpG.SetMaximum(2); gexpG.SetNpx(42.); gexpG.SetNpy(41.);
  gexpGup.SetMinimum(0.000001); gexpGup.SetMaximum(2); gexpGup.SetNpx(42.); gexpGup.SetNpy(41.);
  gexpGdown.SetMinimum(0.000001); gexpGdown.SetMaximum(2); gexpGdown.SetNpx(42.); gexpGdown.SetNpy(41.);
  gexpGup2.SetMinimum(0.000001); gexpGup2.SetMaximum(2); gexpGup2.SetNpx(42.); gexpGup2.SetNpy(41.);
  gexpGdown2.SetMinimum(0.000001); gexpGdown2.SetMaximum(2); gexpGdown2.SetNpx(42.); gexpGdown2.SetNpy(41.);


  glim.SetMinimum(0.00001); glim.SetMaximum(2);
  glim.SetNpx(42.); glim.SetNpy(41.);
  glim.SetTitle(title.c_str());
  glim.GetYaxis()->SetTitleOffset(1.3);
  glim.GetXaxis()->SetTitleOffset(1.3);


  // TLegend l(gStyle->GetPadLeftMargin()+0.05, 0.95-2.*gStyle->GetPadTopMargin(),1.-gStyle->GetPadRightMargin()-0.2, 0.95-gStyle->GetPadTopMargin());
  TLegend l(0.15,0.85,0.45,0.6);
  // l.SetNColumns(2);
  l.SetTextSize(0.04);
  l.SetBorderSize(0);

  TCanvas c("","",1200, 950);
  c.cd();
  c.SetRightMargin(0.15);

  TLatex ltitle(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),
                "#font[62]{CMS}#scale[0.76]{#font[52]{ Supplementary}}");
  TLatex rtitle(1.-c.GetRightMargin()+0.1, 1.-0.5*c.GetTopMargin(),
                "#scale[0.8]{137 fb^{-1} (13 TeV)}");
  ltitle.SetNDC();
  rtitle.SetNDC();
  ltitle.SetTextAlign(12);
  rtitle.SetTextAlign(32);

  // c.SetTopMargin(2.*c.GetTopMargin());
  c.SetLeftMargin(0.12);
  c.SetLogz();
  glim.Draw("colz");

  //My own attempt at gexp
  // TGraph2D gexp(vexp.size()); gexp.SetName("gexp");
  // TGraph2D gup(vup.size()); gup.SetName("gup");
  // TGraph2D gdown(vdown.size()); gdown.SetName("gdown");
  // TGraph2D gup2(vup2.size()); gup2.SetName("gup2");
  // TGraph2D gdown2(vdown2.size()); gdown2.SetName("gdown2");
  //
  // // TGraph2D gexpG(vobs.size()); gexpG.SetName("gexpG");
  // TGraph2D gexpGup(vobsup.size()); gexpGup.SetName("gexpGup");
  // TGraph2D gexpGdown(vobsdown.size()); gexpGdown.SetName("gexpGdown");
  // TGraph2D gexpGup2(vobsup2.size()); gexpGup2.SetName("gexpGup2");
  // TGraph2D gexpGdown2(vobsdown2.size()); gexpGdown2.SetName("gexpGdown2");

  // TH2D * SignifScan= new TH2D("SignifScan","SignifScan",41,200,1225,42,50,1100);
  //
  // for (int i=0;i<vexp.size();i++){
  //   if (vmx[i]<200.0) continue;
  //   gexp.SetPoint(gexp.GetN(), vmx[i]+4, vmy[i]+4, vexp[i]);
  //   gup.SetPoint(gup.GetN(), vmx[i]+4, vmy[i]+4, vup[i]);
  //   gdown.SetPoint(gdown.GetN(), vmx[i]+4, vmy[i]+4, vdown[i]);
  //   gup2.SetPoint(gup2.GetN(), vmx[i]+4, vmy[i]+4, vup2[i]);
  //   gdown2.SetPoint(gdown2.GetN(), vmx[i]+4, vmy[i]+4, vdown2[i]);
  // }
  // for (int i=0;i<vobs.size();i++){
  //   if (vmx[i]<200.0) continue;
  //   // gexpG.SetPoint(gexpG.GetN(), vmx[i]+4, vmy[i]+4, vobs[i]);
  //   gexpGup.SetPoint(gexpGup.GetN(), vmx[i]+4, vmy[i]+4, vobsup[i]);
  //   gexpGdown.SetPoint(gexpGdown.GetN(), vmx[i]+4, vmy[i]+4, vobsdown[i]);
  //   gexpGup2.SetPoint(gexpGup2.GetN(), vmx[i]+4, vmy[i]+4, vobsup2[i]);
  //   gexpGdown2.SetPoint(gexpGdown2.GetN(), vmx[i]+4, vmy[i]+4, vobsdown2[i]);
  //   // std::cout<<"x mass: "<< vmx[i];
  //   // std::cout<<",  y mass: "<< vmy[i]<<std::endl;
  //
  //   int thisBin = SignifScan->FindBin(vmx[i]+4.0,vmy[i]+4.0);
  //   SignifScan->SetBinContent(thisBin,vobs[i]);
  // }

  // std::cout<<"size obs: "<< vobs.size() << ", size x: "<< vmx.size()<< ", size y: "<< vmy.size()<<std::endl;

  // gexpG.Draw();
  // gexpG.SetMinimum(0.00001);gexpG.SetMaximum(300000);



  // gexpG.Draw("surf2");
  glim.Draw("colz");
  // SignifScan->SetStats(0);
  // SignifScan->GetZaxis()->SetRangeUser(1E-6,1E6);
  // SignifScan->Smooth();
  // SignifScan->Draw("colz");


  // double contours[1];
  // contours[0] = 1.0;
  // contours[1] = 0.1;
  // contours[2] = 0.01;

  // // SignifScan->DrawCopy("colz");
  // SignifScan->SetContour(1,contours);
  // SignifScan->Draw("cont3 same");
  // SignifScan->SetLineColor(kRed);

  // SignifScan->SetContour(1.0);
  // SignifScan->Draw("CONT3 SAME");
  // gexpG.Draw("colz");
  // TGraph cG = DrawContours(gexpG, 1, 1, 5, num_smooth_,0.01);


  TGraph cG = DrawContours(gexpG, 1, 1, 3, num_smooth_,1.0);
  TGraph cGup = DrawContours(gexpGup, 1, 2, 2, num_smooth_,1.0);
  TGraph cGdown = DrawContours(gexpGdown, 1, 2, 2, num_smooth_,1.0);
  TGraph cGup2 = DrawContours(gexpGup2, 1, 3, 2, num_smooth_,1.0);
  TGraph cGdown2 = DrawContours(gexpGdown2, 1, 3, 2, num_smooth_,1.0);

  // TGraph cexp = DrawContours(gexp, 2, 1, 3, num_smooth_, 1.0);
  // TGraph cup = DrawContours(gup, 2, 2, 3, num_smooth_, 1.0);
  // TGraph cdown = DrawContours(gdown, 2, 2, 3, num_smooth_, 1.0);
  // TGraph cup2 = DrawContours(gup2, 2, 3, 3, num_smooth_, 1.0);
  // TGraph cdown2 = DrawContours(gdown2, 2, 3, 3, num_smooth_, 1.0);

  // TGraph cup = DrawContours(gup, 2, 2, 5, num_smooth_,1.);
  // TGraph cdown = DrawContours(gdown, 2, 2, 5, num_smooth_,1.);
  // TGraph cexp = DrawContours(gexp, 2, 1, 5, num_smooth_, 1.0);
  // TGraph cobsup = DrawContours(gobsup, 1, 2, 5, num_smooth_,1.);
  // TGraph cobsdown = DrawContours(gobsdown, 1, 2, 5, num_smooth_,1.);
  // TGraph cobs = DrawContours(gobs, 1, 1, 5, num_smooth_, 1.);

  // l.AddEntry(&cexp, "Full higgsino xsec", "l");
  // l.AddEntry(&cup, "1-sigma band", "l");
  // l.AddEntry(&cup2, "2-sigma band", "l");

  l.AddEntry(&cG, "Gluino xsec", "l");
  l.AddEntry(&cGup, "1-sigma band", "l");
  l.AddEntry(&cGup2, "2-sigma band", "l");

  // l.AddEntry(&cexp, "Expected", "l");
  // l.AddEntry(&cobs, "Observed", "l");

  l.Draw("same");
  ltitle.Draw("same");
  rtitle.Draw("same");

  string filebase = out_dir+model_+"_limit_scan";
  if (num_smooth_>0) {
    filebase += "_smooth";
    filebase += to_string(num_smooth_);
  }

  gPad->Update();
  glim.GetZaxis()->SetTitleOffset(1.3);
  glim.GetXaxis()->SetRangeUser(0,1200);
  glim.GetYaxis()->SetRangeUser(50,1100);
  glim.GetYaxis()->SetTitleOffset(1.3);
  glim.GetXaxis()->SetTitleOffset(1.3);
  gPad->Update();
  c.Print((filebase+".pdf").c_str());
  c.Print((filebase+".root").c_str());

  TFile file((filebase+"All.root").c_str(), "recreate");
  glim.GetHistogram()->Write((model_+"ObservedExcludedXsec").c_str());
  // gexp.GetHistogram()->Write((model_+"Test").c_str());
  // gexpGlu.GetHistogram()->Write((model_+"Test2").c_str());
  // cobs.Write((model_+"ObservedLimit").c_str());
  // cobsup.Write((model_+"ObservedLimitUp").c_str());
  // cobsdown.Write((model_+"ObservedLimitDown").c_str());
  // cexp.Write((model_+"ExpectedLimit").c_str());
  // cup.Write((model_+"ExpectedLimitUp").c_str());
  // cdown.Write((model_+"ExpectedLimitDown").c_str());
  file.Close();
  cout << "\nSaved limit curves in " << filebase << ".root\n" << endl;
}

int GetNumBins(const vector<double> &pts, double width) {
  double pmin = *min_element(pts.cbegin(), pts.cend());
  double pmax = *max_element(pts.cbegin(), pts.cend());
  return max(1, min(500, static_cast<int>(ceil((pmax-pmin)/width))));
}

void GetParticleNames(string &xparticle, string &yparticle) {
  if (model_=="N1N2") {
    xparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{2}}}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  }
  else if (model_=="T1tttt") {
    xparticle = "#tilde{g}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  }
  else if (model_=="T5tttt") {
    xparticle = "#tilde{g}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  }
  else if (model_=="T2tt") {
    xparticle = "#tilde{g}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  }
  else if (model_=="T6ttWW") {
    xparticle = "sbottom";
    yparticle = "chargino";
  }
  else {
    //DBG(("Unknown model: "+model_));
    xparticle = "#tilde{g}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  }
}

TLatex GetModelLabel(double x, double y) {
  string lsp = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  string label = "";
  if (model_=="T1tttt") {
    label = "pp #rightarrow #tilde{g}#kern[0.3]{#tilde{g}}, #tilde{g} #rightarrow t#kern[0.4]{#bar{t}}#kern[0.4]{"+lsp+"}";
  }
  else if (model_=="T5tttt") {
    label = "#splitline{pp #rightarrow #tilde{g}#kern[0.3]{#tilde{g}}, #tilde{g} #rightarrow #tilde{t}_{1}t, #tilde{t}_{1} #rightarrow #bar{t}#kern[0.4]{"
      +lsp+"}}{(m#kern[0.3]{_{#lower[-0.12]{#tilde{t}_{1}}}} - m#kern[0.12]{_{"+lsp+"}} = 175 GeV)}";
  }
   else{
    //DBG(("Unknown model: "+model_));
    label = "";
  }
  TLatex l(x, y, label.c_str());
  l.SetNDC();
  l.SetTextAlign(13);
  l.SetTextSize(0.032);
  return l;
}

void Style(TGraph *g, int color, int style, float width) {
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(width);
}

TGraph DrawContours(TGraph2D &g2, int color, int style, double width,
                    int n_smooth, double val) {
  TGraph graph;

  TList *l;
  // Finding the TH2D, smoothing it, and creating a TGraph2D to get a new Delauny interpolation
  // if (n_smooth>0) {
  //   TH2D *histo2d = g2.GetHistogram();
  //   TH2D htemp("", "",
  //              100, histo2d->GetXaxis()->GetXmin(), histo2d->GetXaxis()->GetXmax(),
  //              100, histo2d->GetYaxis()->GetXmin(), histo2d->GetYaxis()->GetXmax());
  //   for (int binx=1; binx<=htemp.GetNbinsX(); ++binx) {
  //     double x = htemp.GetXaxis()->GetBinCenter(binx);
  //     for (int biny=1; biny<=htemp.GetNbinsY(); ++biny) {
  //       double y = htemp.GetYaxis()->GetBinCenter(biny);
  //       double z = g2.Interpolate(x,y);
  //       if (z!=0.) {
  //         htemp.SetBinContent(htemp.GetBin(binx, biny), z);
  //       }
  //     }
  //   }
  //
  //   for (int ind=0; ind<n_smooth; ++ind) {
  //     htemp.Smooth(1,"k5b");
  //   }
  //
  //   vector<double> vx, vy, vz;
  //   double glu_lsp = 225;
  //   for (int binx=1; binx<=htemp.GetNbinsX(); ++binx) {
  //     double x = htemp.GetXaxis()->GetBinCenter(binx);
  //     for (int biny=1; biny<=htemp.GetNbinsY(); ++biny) {
  //       double y = htemp.GetYaxis()->GetBinCenter(biny);
  //       double z = htemp.GetBinContent(htemp.GetBin(binx,biny));
  //
  //       vx.push_back(x);
  //       vy.push_back(y);
  //       int thresh = glu_lsp+30;
  //       if (model_=="T5tttt") thresh = glu_lsp+50;
  //       if (x-y>thresh) {
  //         vz.push_back(z);
  //       }else{
  //         vz.push_back(g2.Interpolate(x,y));
  //       }
  //     }
  //   }
  //
  //   TGraph2D gsmooth("gsmooth", "Cross-Section Limit", vx.size(), &vx.at(0), &vy.at(0), &vz.at(0));
  //   gsmooth.GetHistogram();
  //   l = gsmooth.GetContourList(val);
  // } else {
  //   g2.GetHistogram();
  //   l = g2.GetContourList(val);
  // }

  g2.GetHistogram();
  l = g2.GetContourList(val);
  if (l == nullptr || l->GetSize()==0) {
    std::cout<<"Contour list is empty"<<std::endl;
    return graph;
  }
  int max_points = -1;
  vector<TGraph*> g;
  int startLoop = 0;
  int endLoop = 1;
  int lSize = l->GetSize();
  if (lSize>=1){
    startLoop=1;
    endLoop=2;
  }
  for (int i = 0; i < l->GetSize(); ++i) {
    // for (int i = 0; i < l->GetSize(); ++i) {

    g.push_back(static_cast<TGraph*>(l->At(i)));
    Style(g[i], color, style, width);
    if (g[i] == nullptr) continue;
    int n_points = g[i]->GetN();
    if (n_points > max_points) {
      if (n_smooth>0) FixGraph(*(g[i]));
      TGraph* old_graph = static_cast<TGraph*>(graph.Clone());
      if (i>0) ReverseGraph(*(g[i]));
      graph = *(joinGraphs(old_graph, g[i]));
      graph.RemovePoint(graph.GetN()-1);
      max_points = n_points;
    }
    g[i]->Draw("L same");
  }
  graph.SetTitle(g2.GetTitle());
  graph.SetLineColor(color);
  graph.SetLineWidth(width);
  graph.SetLineStyle(style);
  return graph;
}

TGraph* joinGraphs(TGraph *graph1, TGraph *graph2) {
  TGraph *graph = new TGraph;
  double mglu, mlsp;
  for (int point(0); point < graph1->GetN(); point++) {
    graph1->GetPoint(point, mglu, mlsp);
    graph->SetPoint(graph->GetN(), mglu, mlsp);
  } // Points in graph1
  for (int point(0); point < graph2->GetN(); point++) {
    graph2->GetPoint(point, mglu, mlsp);
    graph->SetPoint(graph->GetN(), mglu, mlsp);
  } // Points in graph1
  graph1->GetPoint(0, mglu, mlsp);
  graph->SetPoint(graph->GetN(), mglu, mlsp);
  TString gname = graph1->GetName(); gname += graph2->GetName();
  graph->SetName(gname);

  return graph;
}

void FixGraph(TGraph &graph) {
  double glu_lsp = 225;
  if (model_=="T5tttt") glu_lsp = 265.;
  int np(graph.GetN());
  double iniglu, endglu, inilsp, endlsp;

  graph.GetPoint(0, iniglu, inilsp);
  graph.GetPoint(np-1, endglu, endlsp);

  // Reversing graph if printed towards decreasing mgluino
  if (inilsp < endlsp) {
    ReverseGraph(graph);
    endglu = iniglu;
    endlsp = inilsp;
  }

  // Adding a point so that it goes down to mLSP = 0, but not for WZ,SOS
  if (endlsp<30) {
    graph.SetPoint(graph.GetN(), endglu, 0);
    np++;
  }

  ReverseGraph(graph);
  // Adding a point at mLSP = 0, and removing points beyond the diagonal
  for (int point(0); point < np; point++) {
    double mglu, mlsp;
    graph.GetPoint(point, mglu, mlsp);
    if (mlsp > mglu-glu_lsp-5) {
      while (point <= graph.GetN() && point!=0) {
        graph.RemovePoint(graph.GetN()-1);
        np--;
      }
      break;
    }
  }
}

void ReverseGraph(TGraph &graph) {
  int np(graph.GetN());
  double mglu, mlsp;
  vector<double> mglus, mlsps;
  for (int point(np-1); point >= 0; point--) {
    graph.GetPoint(point, mglu, mlsp);
    mglus.push_back(mglu);
    mlsps.push_back(mlsp);
  }
  for (int point(0); point < np; point++)
    graph.SetPoint(point, mglus[point], mlsps[point]);
}

void SetupColors() {
  const unsigned num = 5;
  const int bands = 255;
  int colors[bands];
  double stops[num] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red[num] = {0.50, 0.50, 1.00, 1.00, 1.00};
  double green[num] = {0.50, 1.00, 1.00, 0.60, 0.50};
  double blue[num] = {1.00, 1.00, 0.50, 0.40, 0.50};
  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for (int i = 0; i < bands; ++i) colors[i] = fi+i;
  gStyle->SetNumberContours(bands);
  gStyle->SetPalette(bands, colors);
}

void SetupSignedColors() {
  const unsigned num = 3;
  const int bands = 255;
  int colors[bands];
  double stops[num] = {0.0, 0.5, 1.0};
  double red[num]   = {0.0, 1.0, 1.0};
  double green[num] = {0.0, 1.0, 0.0};
  double blue[num]  = {1.0, 1.0, 0.0};
  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for (int i = 0; i < bands; ++i) colors[i] = fi+i;
  gStyle->SetNumberContours(bands);
  gStyle->SetPalette(bands, colors);
}

std::vector<std::string> split(const std::string& s, char delimiter) {
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter)) { tokens.push_back(token);}
   return tokens;
}

void higgsinoCrossSection(int hig_mass, double &xsec, double &xsec_unc) {
  if (hig_mass ==127) { xsec = .5824*.5824*7.6022; xsec_unc = 0.0393921; return;}
  else if (hig_mass ==150) { xsec = .5824*.5824*3.83231; xsec_unc = 0.0413612; return;}
  else if (hig_mass ==175) { xsec = .5824*.5824*2.26794; xsec_unc = 0.044299; return;}
  else if (hig_mass ==200) { xsec = .5824*.5824*1.33562; xsec_unc = 0.0474362; return;}
  else if (hig_mass ==225) { xsec = .5824*.5824*0.860597; xsec_unc = 0.0504217; return;}
  else if (hig_mass ==250) { xsec = .5824*.5824*0.577314; xsec_unc = 0.0532731; return;}
  else if (hig_mass ==275) { xsec = .5824*.5824*0.400107; xsec_unc = 0.0560232; return;}
  else if (hig_mass ==300) { xsec = .5824*.5824*0.284855; xsec_unc = 0.0586867; return;}
  else if (hig_mass ==325) { xsec = .5824*.5824*0.20736; xsec_unc = 0.0613554; return;}
  else if (hig_mass ==350) { xsec = .5824*.5824*0.153841; xsec_unc = 0.0640598; return;}
  else if (hig_mass ==375) { xsec = .5824*.5824*0.116006; xsec_unc = 0.066892; return;}
  else if (hig_mass ==400) { xsec = .5824*.5824*0.0887325; xsec_unc = 0.0697517; return;}
  else if (hig_mass ==425) { xsec = .5824*.5824*0.0686963; xsec_unc = 0.0723531; return;}
  else if (hig_mass ==450) { xsec = .5824*.5824*0.0537702; xsec_unc = 0.0748325; return;}
  else if (hig_mass ==475) { xsec = .5824*.5824*0.0424699; xsec_unc = 0.0775146; return;}
  else if (hig_mass ==500) { xsec = .5824*.5824*0.0338387; xsec_unc = 0.0802572; return;}
  else if (hig_mass ==525) { xsec = .5824*.5824*0.0271867; xsec_unc = 0.0825803; return;}
  else if (hig_mass ==550) { xsec = .5824*.5824*0.0219868; xsec_unc = 0.0849278; return;}
  else if (hig_mass ==575) { xsec = .5824*.5824*0.0179062; xsec_unc = 0.087561; return;}
  else if (hig_mass ==600) { xsec = .5824*.5824*0.0146677; xsec_unc = 0.0900693; return;}
  else if (hig_mass ==625) { xsec = .5824*.5824*0.012062; xsec_unc = 0.091959; return;}
  else if (hig_mass ==650) { xsec = .5824*.5824*0.00996406; xsec_unc = 0.094065; return;}
  else if (hig_mass ==675) { xsec = .5824*.5824*0.00828246; xsec_unc = 0.0957436; return;}
  else if (hig_mass ==700) { xsec = .5824*.5824*0.00689981; xsec_unc = 0.0982894; return;}
  else if (hig_mass ==725) { xsec = .5824*.5824*0.00578355; xsec_unc = 0.0999915; return;}
  else if (hig_mass ==750) { xsec = .5824*.5824*0.0048731; xsec_unc = 0.101211; return;}
  else if (hig_mass ==775) { xsec = .5824*.5824*0.00409781; xsec_unc = 0.104646; return;}
  else if (hig_mass ==800) { xsec = .5824*.5824*0.00346143; xsec_unc = 0.107618; return;}
  else if (hig_mass ==825) { xsec = .5824*.5824*0.0029337; xsec_unc = 0.108353; return;}
  else if (hig_mass ==850) { xsec = .5824*.5824*0.0024923; xsec_unc = 0.110016; return;}
  else if (hig_mass ==875) { xsec = .5824*.5824*0.00213679; xsec_unc = 0.112636; return;}
  else if (hig_mass ==900) { xsec = .5824*.5824*0.00180616; xsec_unc = 0.1134; return;}
  else if (hig_mass ==925) { xsec = .5824*.5824*0.00155453; xsec_unc = 0.116949; return;}
  else if (hig_mass ==950) { xsec = .5824*.5824*0.00132692; xsec_unc = 0.117027; return;}
  else if (hig_mass ==975) { xsec = .5824*.5824*0.00112975; xsec_unc = 0.121244; return;}
  else if (hig_mass ==1000) { xsec = .5824*.5824*0.000968853; xsec_unc = 0.126209; return;}
  else if (hig_mass ==1025) { xsec = .5824*.5824*0.000840602; xsec_unc = 0.121654; return;}
  else if (hig_mass ==1050) { xsec = .5824*.5824*0.000731306; xsec_unc = 0.118502; return;}
  else if (hig_mass ==1075) { xsec = .5824*.5824*0.000627083; xsec_unc = 0.127723; return;}
  else if (hig_mass ==1100) { xsec = .5824*.5824*0.000538005; xsec_unc = 0.134099; return;}
  else if (hig_mass ==1125) { xsec = .5824*.5824*0.00046747; xsec_unc = 0.133755; return;}
  else if (hig_mass ==1150) { xsec = .5824*.5824*0.000405108; xsec_unc = 0.120607; return;}
  else if (hig_mass ==1175) { xsec = .5824*.5824*0.000348261; xsec_unc = 0.139744; return;}
  else if (hig_mass ==1200) { xsec = .5824*.5824*0.000299347; xsec_unc = 0.162604; return;}
  else if (hig_mass ==1225) { xsec = .5824*.5824*0.000265935; xsec_unc = 0.137575; return;}
  else if (hig_mass ==1250) { xsec = .5824*.5824*0.000240471; xsec_unc = 0.119271; return;}
  else if (hig_mass ==1275) { xsec = .5824*.5824*0.000190411; xsec_unc = 0.138061; return;}
  else if (hig_mass ==1300) { xsec = .5824*.5824*0.000160765; xsec_unc = 0.122224; return;}
  else if (hig_mass ==1325) { xsec = .5824*.5824*0.000136272; xsec_unc = 0.138533; return;}
  else if (hig_mass ==1350) { xsec = .5824*.5824*0.000111174; xsec_unc = 0.177681; return;}
  else if (hig_mass ==1375) { xsec = .5824*.5824*9.74728e-05; xsec_unc = 0.138992; return;}
  else if (hig_mass ==1400) { xsec = .5824*.5824*7.80263e-05; xsec_unc = 0.118718; return;}
  else if (hig_mass ==1425) { xsec = .5824*.5824*6.96843e-05; xsec_unc = 0.139439; return;}
  else if (hig_mass ==1450) { xsec = .5824*.5824*6.96962e-05; xsec_unc = 0.198887; return;}
  else if (hig_mass ==1475) { xsec = .5824*.5824*4.98006e-05; xsec_unc = 0.139874; return;}
  else { xsec = 0; xsec_unc = 0;}
}

void gluinoCrossSection(int glu_mass, double &xsec, double &xsec_unc){
  //Setting all xsec below 600 GeV to be the xsec of 600 GeV
  if      (glu_mass < 600.  ) { xsec = 0.113E+02; xsec_unc = 8.62 /100; return; } // we shouldn't have these points
  else if (glu_mass == 600. ) { xsec = 0.113E+02; xsec_unc = 8.62 /100; return; }
  else if (glu_mass == 605. ) { xsec = 0.108E+02; xsec_unc = 8.66 /100; return; }
  else if (glu_mass == 610. ) { xsec = 0.102E+02; xsec_unc = 8.7 /100; return; }
  else if (glu_mass == 615. ) { xsec = 0.974E+01; xsec_unc = 8.74 /100; return; }
  else if (glu_mass == 620. ) { xsec = 0.926E+01; xsec_unc = 8.78 /100; return; }
  else if (glu_mass == 625. ) { xsec = 0.881E+01; xsec_unc = 8.82 /100; return; }
  else if (glu_mass == 630. ) { xsec = 0.839E+01; xsec_unc = 8.86 /100; return; }
  else if (glu_mass == 635. ) { xsec = 0.799E+01; xsec_unc = 8.9 /100; return; }
  else if (glu_mass == 640. ) { xsec = 0.761E+01; xsec_unc = 8.94 /100; return; }
  else if (glu_mass == 645. ) { xsec = 0.725E+01; xsec_unc = 8.98 /100; return; }
  else if (glu_mass == 650. ) { xsec = 0.690E+01; xsec_unc = 9.02 /100; return; }
  else if (glu_mass == 655. ) { xsec = 0.658E+01; xsec_unc = 9.06 /100; return; }
  else if (glu_mass == 660. ) { xsec = 0.627E+01; xsec_unc = 9.1 /100; return; }
  else if (glu_mass == 665. ) { xsec = 0.598E+01; xsec_unc = 9.15 /100; return; }
  else if (glu_mass == 670. ) { xsec = 0.571E+01; xsec_unc = 9.19 /100; return; }
  else if (glu_mass == 675. ) { xsec = 0.544E+01; xsec_unc = 9.23 /100; return; }
  else if (glu_mass == 680. ) { xsec = 0.520E+01; xsec_unc = 9.27 /100; return; }
  else if (glu_mass == 685. ) { xsec = 0.496E+01; xsec_unc = 9.31 /100; return; }
  else if (glu_mass == 690. ) { xsec = 0.474E+01; xsec_unc = 9.35 /100; return; }
  else if (glu_mass == 695. ) { xsec = 0.452E+01; xsec_unc = 9.39 /100; return; }
  else if (glu_mass == 700. ) { xsec = 0.432E+01; xsec_unc = 9.43 /100; return; }
  else if (glu_mass == 705. ) { xsec = 0.413E+01; xsec_unc = 9.46 /100; return; }
  else if (glu_mass == 710. ) { xsec = 0.395E+01; xsec_unc = 9.5 /100; return; }
  else if (glu_mass == 715. ) { xsec = 0.377E+01; xsec_unc = 9.54 /100; return; }
  else if (glu_mass == 720. ) { xsec = 0.361E+01; xsec_unc = 9.58 /100; return; }
  else if (glu_mass == 725. ) { xsec = 0.345E+01; xsec_unc = 9.61 /100; return; }
  else if (glu_mass == 730. ) { xsec = 0.330E+01; xsec_unc = 9.65 /100; return; }
  else if (glu_mass == 735. ) { xsec = 0.316E+01; xsec_unc = 9.69 /100; return; }
  else if (glu_mass == 740. ) { xsec = 0.302E+01; xsec_unc = 9.72 /100; return; }
  else if (glu_mass == 745. ) { xsec = 0.289E+01; xsec_unc = 9.76 /100; return; }
  else if (glu_mass == 750. ) { xsec = 0.277E+01; xsec_unc = 9.8 /100; return; }
  else if (glu_mass == 755. ) { xsec = 0.265E+01; xsec_unc = 9.83 /100; return; }
  else if (glu_mass == 760. ) { xsec = 0.254E+01; xsec_unc = 9.87 /100; return; }
  else if (glu_mass == 765. ) { xsec = 0.243E+01; xsec_unc = 9.91 /100; return; }
  else if (glu_mass == 770. ) { xsec = 0.233E+01; xsec_unc = 9.94 /100; return; }
  else if (glu_mass == 775. ) { xsec = 0.223E+01; xsec_unc = 9.98 /100; return; }
  else if (glu_mass == 780. ) { xsec = 0.214E+01; xsec_unc = 10.01 /100; return; }
  else if (glu_mass == 785. ) { xsec = 0.205E+01; xsec_unc = 10.05 /100; return; }
  else if (glu_mass == 790. ) { xsec = 0.197E+01; xsec_unc = 10.09 /100; return; }
  else if (glu_mass == 795. ) { xsec = 0.188E+01; xsec_unc = 10.12 /100; return; }
  else if (glu_mass == 800. ) { xsec = 0.181E+01; xsec_unc = 10.16 /100; return; }
  else if (glu_mass == 805. ) { xsec = 0.173E+01; xsec_unc = 10.2 /100; return; }
  else if (glu_mass == 810. ) { xsec = 0.166E+01; xsec_unc = 10.23 /100; return; }
  else if (glu_mass == 815. ) { xsec = 0.160E+01; xsec_unc = 10.27 /100; return; }
  else if (glu_mass == 820. ) { xsec = 0.153E+01; xsec_unc = 10.31 /100; return; }
  else if (glu_mass == 825. ) { xsec = 0.147E+01; xsec_unc = 10.34 /100; return; }
  else if (glu_mass == 830. ) { xsec = 0.141E+01; xsec_unc = 10.38 /100; return; }
  else if (glu_mass == 835. ) { xsec = 0.136E+01; xsec_unc = 10.42 /100; return; }
  else if (glu_mass == 840. ) { xsec = 0.130E+01; xsec_unc = 10.45 /100; return; }
  else if (glu_mass == 845. ) { xsec = 0.125E+01; xsec_unc = 10.49 /100; return; }
  else if (glu_mass == 850. ) { xsec = 0.120E+01; xsec_unc = 10.53 /100; return; }
  else if (glu_mass == 855. ) { xsec = 0.115E+01; xsec_unc = 10.57 /100; return; }
  else if (glu_mass == 860. ) { xsec = 0.111E+01; xsec_unc = 10.6 /100; return; }
  else if (glu_mass == 865. ) { xsec = 0.107E+01; xsec_unc = 10.64 /100; return; }
  else if (glu_mass == 870. ) { xsec = 0.103E+01; xsec_unc = 10.68 /100; return; }
  else if (glu_mass == 875. ) { xsec = 0.986E+00; xsec_unc = 10.71 /100; return; }
  else if (glu_mass == 880. ) { xsec = 0.948E+00; xsec_unc = 10.75 /100; return; }
  else if (glu_mass == 885. ) { xsec = 0.912E+00; xsec_unc = 10.79 /100; return; }
  else if (glu_mass == 890. ) { xsec = 0.877E+00; xsec_unc = 10.82 /100; return; }
  else if (glu_mass == 895. ) { xsec = 0.844E+00; xsec_unc = 10.86 /100; return; }
  else if (glu_mass == 900. ) { xsec = 0.812E+00; xsec_unc = 10.89 /100; return; }
  else if (glu_mass == 905. ) { xsec = 0.781E+00; xsec_unc = 10.93 /100; return; }
  else if (glu_mass == 910. ) { xsec = 0.752E+00; xsec_unc = 10.97 /100; return; }
  else if (glu_mass == 915. ) { xsec = 0.723E+00; xsec_unc = 11.0 /100; return; }
  else if (glu_mass == 920. ) { xsec = 0.696E+00; xsec_unc = 11.04 /100; return; }
  else if (glu_mass == 925. ) { xsec = 0.670E+00; xsec_unc = 11.07 /100; return; }
  else if (glu_mass == 930. ) { xsec = 0.646E+00; xsec_unc = 11.11 /100; return; }
  else if (glu_mass == 935. ) { xsec = 0.622E+00; xsec_unc = 11.14 /100; return; }
  else if (glu_mass == 940. ) { xsec = 0.599E+00; xsec_unc = 11.18 /100; return; }
  else if (glu_mass == 945. ) { xsec = 0.577E+00; xsec_unc = 11.21 /100; return; }
  else if (glu_mass == 950. ) { xsec = 0.556E+00; xsec_unc = 11.25 /100; return; }
  else if (glu_mass == 955. ) { xsec = 0.535E+00; xsec_unc = 11.28 /100; return; }
  else if (glu_mass == 960. ) { xsec = 0.516E+00; xsec_unc = 11.32 /100; return; }
  else if (glu_mass == 965. ) { xsec = 0.497E+00; xsec_unc = 11.35 /100; return; }
  else if (glu_mass == 970. ) { xsec = 0.479E+00; xsec_unc = 11.39 /100; return; }
  else if (glu_mass == 975. ) { xsec = 0.462E+00; xsec_unc = 11.42 /100; return; }
  else if (glu_mass == 980. ) { xsec = 0.445E+00; xsec_unc = 11.46 /100; return; }
  else if (glu_mass == 985. ) { xsec = 0.430E+00; xsec_unc = 11.49 /100; return; }
  else if (glu_mass == 990. ) { xsec = 0.414E+00; xsec_unc = 11.53 /100; return; }
  else if (glu_mass == 995. ) { xsec = 0.399E+00; xsec_unc = 11.56 /100; return; }
  else if (glu_mass == 1000 ) { xsec = 0.385E+00; xsec_unc = 11.6 /100; return; }
  else if (glu_mass == 1005 ) { xsec = 0.372E+00; xsec_unc = 11.63 /100; return; }
  else if (glu_mass == 1010 ) { xsec = 0.359E+00; xsec_unc = 11.67 /100; return; }
  else if (glu_mass == 1015 ) { xsec = 0.346E+00; xsec_unc = 11.7 /100; return; }
  else if (glu_mass == 1020 ) { xsec = 0.334E+00; xsec_unc = 11.74 /100; return; }
  else if (glu_mass == 1025 ) { xsec = 0.322E+00; xsec_unc = 11.78 /100; return; }
  else if (glu_mass == 1030 ) { xsec = 0.311E+00; xsec_unc = 11.81 /100; return; }
  else if (glu_mass == 1035 ) { xsec = 0.300E+00; xsec_unc = 11.85 /100; return; }
  else if (glu_mass == 1040 ) { xsec = 0.290E+00; xsec_unc = 11.88 /100; return; }
  else if (glu_mass == 1045 ) { xsec = 0.280E+00; xsec_unc = 11.92 /100; return; }
  else if (glu_mass == 1050 ) { xsec = 0.270E+00; xsec_unc = 11.95 /100; return; }
  else if (glu_mass == 1055 ) { xsec = 0.261E+00; xsec_unc = 11.99 /100; return; }
  else if (glu_mass == 1060 ) { xsec = 0.252E+00; xsec_unc = 12.02 /100; return; }
  else if (glu_mass == 1065 ) { xsec = 0.243E+00; xsec_unc = 12.06 /100; return; }
  else if (glu_mass == 1070 ) { xsec = 0.235E+00; xsec_unc = 12.09 /100; return; }
  else if (glu_mass == 1075 ) { xsec = 0.227E+00; xsec_unc = 12.13 /100; return; }
  else if (glu_mass == 1080 ) { xsec = 0.219E+00; xsec_unc = 12.17 /100; return; }
  else if (glu_mass == 1085 ) { xsec = 0.212E+00; xsec_unc = 12.2 /100; return; }
  else if (glu_mass == 1090 ) { xsec = 0.205E+00; xsec_unc = 12.24 /100; return; }
  else if (glu_mass == 1095 ) { xsec = 0.198E+00; xsec_unc = 12.27 /100; return; }
  else if (glu_mass == 1100 ) { xsec = 0.191E+00; xsec_unc = 12.31 /100; return; }
  else if (glu_mass == 1105 ) { xsec = 0.185E+00; xsec_unc = 12.34 /100; return; }
  else if (glu_mass == 1110 ) { xsec = 0.179E+00; xsec_unc = 12.38 /100; return; }
  else if (glu_mass == 1115 ) { xsec = 0.173E+00; xsec_unc = 12.42 /100; return; }
  else if (glu_mass == 1120 ) { xsec = 0.167E+00; xsec_unc = 12.45 /100; return; }
  else if (glu_mass == 1125 ) { xsec = 0.162E+00; xsec_unc = 12.49 /100; return; }
  else if (glu_mass == 1130 ) { xsec = 0.156E+00; xsec_unc = 12.53 /100; return; }
  else if (glu_mass == 1135 ) { xsec = 0.151E+00; xsec_unc = 12.56 /100; return; }
  else if (glu_mass == 1140 ) { xsec = 0.146E+00; xsec_unc = 12.6 /100; return; }
  else if (glu_mass == 1145 ) { xsec = 0.141E+00; xsec_unc = 12.64 /100; return; }
  else if (glu_mass == 1150 ) { xsec = 0.137E+00; xsec_unc = 12.67 /100; return; }
  else if (glu_mass == 1155 ) { xsec = 0.132E+00; xsec_unc = 12.71 /100; return; }
  else if (glu_mass == 1160 ) { xsec = 0.128E+00; xsec_unc = 12.74 /100; return; }
  else if (glu_mass == 1165 ) { xsec = 0.124E+00; xsec_unc = 12.78 /100; return; }
  else if (glu_mass == 1170 ) { xsec = 0.120E+00; xsec_unc = 12.82 /100; return; }
  else if (glu_mass == 1175 ) { xsec = 0.116E+00; xsec_unc = 12.85 /100; return; }
  else if (glu_mass == 1180 ) { xsec = 0.112E+00; xsec_unc = 12.89 /100; return; }
  else if (glu_mass == 1185 ) { xsec = 0.109E+00; xsec_unc = 12.92 /100; return; }
  else if (glu_mass == 1190 ) { xsec = 0.105E+00; xsec_unc = 12.96 /100; return; }
  else if (glu_mass == 1195 ) { xsec = 0.102E+00; xsec_unc = 13.0 /100; return; }
  else if (glu_mass == 1200 ) { xsec = 0.985E-01; xsec_unc = 13.03 /100; return; }
  else if (glu_mass == 1205 ) { xsec = 0.953E-01; xsec_unc = 13.07 /100; return; }
  else if (glu_mass == 1210 ) { xsec = 0.923E-01; xsec_unc = 13.1 /100; return; }
  else if (glu_mass == 1215 ) { xsec = 0.894E-01; xsec_unc = 13.14 /100; return; }
  else if (glu_mass == 1220 ) { xsec = 0.866E-01; xsec_unc = 13.17 /100; return; }
  else if (glu_mass == 1225 ) { xsec = 0.838E-01; xsec_unc = 13.21 /100; return; }
  else if (glu_mass == 1230 ) { xsec = 0.812E-01; xsec_unc = 13.24 /100; return; }
  else if (glu_mass == 1235 ) { xsec = 0.786E-01; xsec_unc = 13.27 /100; return; }
  else if (glu_mass == 1240 ) { xsec = 0.762E-01; xsec_unc = 13.31 /100; return; }
  else if (glu_mass == 1245 ) { xsec = 0.738E-01; xsec_unc = 13.34 /100; return; }
  else if (glu_mass == 1250 ) { xsec = 0.715E-01; xsec_unc = 13.38 /100; return; }
  else if (glu_mass == 1255 ) { xsec = 0.692E-01; xsec_unc = 13.41 /100; return; }
  else if (glu_mass == 1260 ) { xsec = 0.671E-01; xsec_unc = 13.45 /100; return; }
  else if (glu_mass == 1265 ) { xsec = 0.650E-01; xsec_unc = 13.48 /100; return; }
  else if (glu_mass == 1270 ) { xsec = 0.630E-01; xsec_unc = 13.51 /100; return; }
  else if (glu_mass == 1275 ) { xsec = 0.610E-01; xsec_unc = 13.55 /100; return; }
  else if (glu_mass == 1280 ) { xsec = 0.591E-01; xsec_unc = 13.58 /100; return; }
  else if (glu_mass == 1285 ) { xsec = 0.573E-01; xsec_unc = 13.62 /100; return; }
  else if (glu_mass == 1290 ) { xsec = 0.556E-01; xsec_unc = 13.65 /100; return; }
  else if (glu_mass == 1295 ) { xsec = 0.539E-01; xsec_unc = 13.69 /100; return; }
  else if (glu_mass == 1300 ) { xsec = 0.522E-01; xsec_unc = 13.72 /100; return; }
  else if (glu_mass == 1305 ) { xsec = 0.506E-01; xsec_unc = 13.76 /100; return; }
  else if (glu_mass == 1310 ) { xsec = 0.491E-01; xsec_unc = 13.79 /100; return; }
  else if (glu_mass == 1315 ) { xsec = 0.476E-01; xsec_unc = 13.83 /100; return; }
  else if (glu_mass == 1320 ) { xsec = 0.461E-01; xsec_unc = 13.86 /100; return; }
  else if (glu_mass == 1325 ) { xsec = 0.447E-01; xsec_unc = 13.9 /100; return; }
  else if (glu_mass == 1330 ) { xsec = 0.434E-01; xsec_unc = 13.94 /100; return; }
  else if (glu_mass == 1335 ) { xsec = 0.421E-01; xsec_unc = 13.97 /100; return; }
  else if (glu_mass == 1340 ) { xsec = 0.408E-01; xsec_unc = 14.01 /100; return; }
  else if (glu_mass == 1345 ) { xsec = 0.396E-01; xsec_unc = 14.04 /100; return; }
  else if (glu_mass == 1350 ) { xsec = 0.384E-01; xsec_unc = 14.08 /100; return; }
  else if (glu_mass == 1355 ) { xsec = 0.372E-01; xsec_unc = 14.11 /100; return; }
  else if (glu_mass == 1360 ) { xsec = 0.361E-01; xsec_unc = 14.15 /100; return; }
  else if (glu_mass == 1365 ) { xsec = 0.350E-01; xsec_unc = 14.19 /100; return; }
  else if (glu_mass == 1370 ) { xsec = 0.340E-01; xsec_unc = 14.22 /100; return; }
  else if (glu_mass == 1375 ) { xsec = 0.330E-01; xsec_unc = 14.26 /100; return; }
  else if (glu_mass == 1380 ) { xsec = 0.320E-01; xsec_unc = 14.3 /100; return; }
  else if (glu_mass == 1385 ) { xsec = 0.310E-01; xsec_unc = 14.33 /100; return; }
  else if (glu_mass == 1390 ) { xsec = 0.301E-01; xsec_unc = 14.37 /100; return; }
  else if (glu_mass == 1395 ) { xsec = 0.292E-01; xsec_unc = 14.4 /100; return; }
  else if (glu_mass == 1400 ) { xsec = 0.284E-01; xsec_unc = 14.44 /100; return; }
  else if (glu_mass == 1405 ) { xsec = 0.275E-01; xsec_unc = 14.48 /100; return; }
  else if (glu_mass == 1410 ) { xsec = 0.267E-01; xsec_unc = 14.51 /100; return; }
  else if (glu_mass == 1415 ) { xsec = 0.259E-01; xsec_unc = 14.55 /100; return; }
  else if (glu_mass == 1420 ) { xsec = 0.252E-01; xsec_unc = 14.59 /100; return; }
  else if (glu_mass == 1425 ) { xsec = 0.244E-01; xsec_unc = 14.63 /100; return; }
  else if (glu_mass == 1430 ) { xsec = 0.237E-01; xsec_unc = 14.66 /100; return; }
  else if (glu_mass == 1435 ) { xsec = 0.230E-01; xsec_unc = 14.7 /100; return; }
  else if (glu_mass == 1440 ) { xsec = 0.224E-01; xsec_unc = 14.74 /100; return; }
  else if (glu_mass == 1445 ) { xsec = 0.217E-01; xsec_unc = 14.77 /100; return; }
  else if (glu_mass == 1450 ) { xsec = 0.211E-01; xsec_unc = 14.81 /100; return; }
  else if (glu_mass == 1455 ) { xsec = 0.205E-01; xsec_unc = 14.85 /100; return; }
  else if (glu_mass == 1460 ) { xsec = 0.199E-01; xsec_unc = 14.88 /100; return; }
  else if (glu_mass == 1465 ) { xsec = 0.193E-01; xsec_unc = 14.92 /100; return; }
  else if (glu_mass == 1470 ) { xsec = 0.187E-01; xsec_unc = 14.96 /100; return; }
  else if (glu_mass == 1475 ) { xsec = 0.182E-01; xsec_unc = 15.0 /100; return; }
  else if (glu_mass == 1480 ) { xsec = 0.177E-01; xsec_unc = 15.03 /100; return; }
  else if (glu_mass == 1485 ) { xsec = 0.172E-01; xsec_unc = 15.07 /100; return; }
  else if (glu_mass == 1490 ) { xsec = 0.167E-01; xsec_unc = 15.11 /100; return; }
  else if (glu_mass == 1495 ) { xsec = 0.162E-01; xsec_unc = 15.15 /100; return; }
  else if (glu_mass == 1500 ) { xsec = 0.157E-01; xsec_unc = 15.18 /100; return; }
  else if (glu_mass == 1505 ) { xsec = 0.153E-01; xsec_unc = 15.22 /100; return; }
  else if (glu_mass == 1510 ) { xsec = 0.148E-01; xsec_unc = 15.26 /100; return; }
  else if (glu_mass == 1515 ) { xsec = 0.144E-01; xsec_unc = 15.3 /100; return; }
  else if (glu_mass == 1520 ) { xsec = 0.140E-01; xsec_unc = 15.33 /100; return; }
  else if (glu_mass == 1525 ) { xsec = 0.136E-01; xsec_unc = 15.37 /100; return; }
  else if (glu_mass == 1530 ) { xsec = 0.132E-01; xsec_unc = 15.41 /100; return; }
  else if (glu_mass == 1535 ) { xsec = 0.128E-01; xsec_unc = 15.45 /100; return; }
  else if (glu_mass == 1540 ) { xsec = 0.125E-01; xsec_unc = 15.48 /100; return; }
  else if (glu_mass == 1545 ) { xsec = 0.121E-01; xsec_unc = 15.52 /100; return; }
  else if (glu_mass == 1550 ) { xsec = 0.118E-01; xsec_unc = 15.56 /100; return; }
  else if (glu_mass == 1555 ) { xsec = 0.115E-01; xsec_unc = 15.6 /100; return; }
  else if (glu_mass == 1560 ) { xsec = 0.111E-01; xsec_unc = 15.64 /100; return; }
  else if (glu_mass == 1565 ) { xsec = 0.108E-01; xsec_unc = 15.67 /100; return; }
  else if (glu_mass == 1570 ) { xsec = 0.105E-01; xsec_unc = 15.71 /100; return; }
  else if (glu_mass == 1575 ) { xsec = 0.102E-01; xsec_unc = 15.75 /100; return; }
  else if (glu_mass == 1580 ) { xsec = 0.993E-02; xsec_unc = 15.79 /100; return; }
  else if (glu_mass == 1585 ) { xsec = 0.966E-02; xsec_unc = 15.83 /100; return; }
  else if (glu_mass == 1590 ) { xsec = 0.939E-02; xsec_unc = 15.87 /100; return; }
  else if (glu_mass == 1595 ) { xsec = 0.912E-02; xsec_unc = 15.9 /100; return; }
  else if (glu_mass == 1600 ) { xsec = 0.887E-02; xsec_unc = 15.94 /100; return; }
  else if (glu_mass == 1605 ) { xsec = 0.862E-02; xsec_unc = 15.98 /100; return; }
  else if (glu_mass == 1610 ) { xsec = 0.838E-02; xsec_unc = 16.02 /100; return; }
  else if (glu_mass == 1615 ) { xsec = 0.815E-02; xsec_unc = 16.06 /100; return; }
  else if (glu_mass == 1620 ) { xsec = 0.792E-02; xsec_unc = 16.1 /100; return; }
  else if (glu_mass == 1625 ) { xsec = 0.770E-02; xsec_unc = 16.13 /100; return; }
  else if (glu_mass == 1630 ) { xsec = 0.749E-02; xsec_unc = 16.17 /100; return; }
  else if (glu_mass == 1635 ) { xsec = 0.728E-02; xsec_unc = 16.21 /100; return; }
  else if (glu_mass == 1640 ) { xsec = 0.708E-02; xsec_unc = 16.25 /100; return; }
  else if (glu_mass == 1645 ) { xsec = 0.689E-02; xsec_unc = 16.29 /100; return; }
  else if (glu_mass == 1650 ) { xsec = 0.670E-02; xsec_unc = 16.33 /100; return; }
  else if (glu_mass == 1655 ) { xsec = 0.651E-02; xsec_unc = 16.37 /100; return; }
  else if (glu_mass == 1660 ) { xsec = 0.633E-02; xsec_unc = 16.41 /100; return; }
  else if (glu_mass == 1665 ) { xsec = 0.616E-02; xsec_unc = 16.45 /100; return; }
  else if (glu_mass == 1670 ) { xsec = 0.599E-02; xsec_unc = 16.49 /100; return; }
  else if (glu_mass == 1675 ) { xsec = 0.583E-02; xsec_unc = 16.53 /100; return; }
  else if (glu_mass == 1680 ) { xsec = 0.567E-02; xsec_unc = 16.56 /100; return; }
  else if (glu_mass == 1685 ) { xsec = 0.551E-02; xsec_unc = 16.6 /100; return; }
  else if (glu_mass == 1690 ) { xsec = 0.536E-02; xsec_unc = 16.64 /100; return; }
  else if (glu_mass == 1695 ) { xsec = 0.521E-02; xsec_unc = 16.68 /100; return; }
  else if (glu_mass == 1700 ) { xsec = 0.507E-02; xsec_unc = 16.72 /100; return; }
  else if (glu_mass == 1705 ) { xsec = 0.493E-02; xsec_unc = 16.76 /100; return; }
  else if (glu_mass == 1710 ) { xsec = 0.480E-02; xsec_unc = 16.81 /100; return; }
  else if (glu_mass == 1715 ) { xsec = 0.467E-02; xsec_unc = 16.85 /100; return; }
  else if (glu_mass == 1720 ) { xsec = 0.454E-02; xsec_unc = 16.89 /100; return; }
  else if (glu_mass == 1725 ) { xsec = 0.442E-02; xsec_unc = 16.93 /100; return; }
  else if (glu_mass == 1730 ) { xsec = 0.430E-02; xsec_unc = 16.97 /100; return; }
  else if (glu_mass == 1735 ) { xsec = 0.418E-02; xsec_unc = 17.01 /100; return; }
  else if (glu_mass == 1740 ) { xsec = 0.407E-02; xsec_unc = 17.05 /100; return; }
  else if (glu_mass == 1745 ) { xsec = 0.396E-02; xsec_unc = 17.09 /100; return; }
  else if (glu_mass == 1750 ) { xsec = 0.385E-02; xsec_unc = 17.13 /100; return; }
  else if (glu_mass == 1755 ) { xsec = 0.375E-02; xsec_unc = 17.18 /100; return; }
  else if (glu_mass == 1760 ) { xsec = 0.365E-02; xsec_unc = 17.22 /100; return; }
  else if (glu_mass == 1765 ) { xsec = 0.355E-02; xsec_unc = 17.26 /100; return; }
  else if (glu_mass == 1770 ) { xsec = 0.345E-02; xsec_unc = 17.3 /100; return; }
  else if (glu_mass == 1775 ) { xsec = 0.336E-02; xsec_unc = 17.34 /100; return; }
  else if (glu_mass == 1780 ) { xsec = 0.327E-02; xsec_unc = 17.39 /100; return; }
  else if (glu_mass == 1785 ) { xsec = 0.318E-02; xsec_unc = 17.43 /100; return; }
  else if (glu_mass == 1790 ) { xsec = 0.310E-02; xsec_unc = 17.47 /100; return; }
  else if (glu_mass == 1795 ) { xsec = 0.301E-02; xsec_unc = 17.51 /100; return; }
  else if (glu_mass == 1800 ) { xsec = 0.293E-02; xsec_unc = 17.56 /100; return; }
  else if (glu_mass == 1805 ) { xsec = 0.286E-02; xsec_unc = 17.6 /100; return; }
  else if (glu_mass == 1810 ) { xsec = 0.278E-02; xsec_unc = 17.64 /100; return; }
  else if (glu_mass == 1815 ) { xsec = 0.271E-02; xsec_unc = 17.69 /100; return; }
  else if (glu_mass == 1820 ) { xsec = 0.263E-02; xsec_unc = 17.73 /100; return; }
  else if (glu_mass == 1825 ) { xsec = 0.256E-02; xsec_unc = 17.77 /100; return; }
  else if (glu_mass == 1830 ) { xsec = 0.249E-02; xsec_unc = 17.82 /100; return; }
  else if (glu_mass == 1835 ) { xsec = 0.243E-02; xsec_unc = 17.86 /100; return; }
  else if (glu_mass == 1840 ) { xsec = 0.236E-02; xsec_unc = 17.9 /100; return; }
  else if (glu_mass == 1845 ) { xsec = 0.230E-02; xsec_unc = 17.95 /100; return; }
  else if (glu_mass == 1850 ) { xsec = 0.224E-02; xsec_unc = 17.99 /100; return; }
  else if (glu_mass == 1855 ) { xsec = 0.218E-02; xsec_unc = 18.04 /100; return; }
  else if (glu_mass == 1860 ) { xsec = 0.212E-02; xsec_unc = 18.08 /100; return; }
  else if (glu_mass == 1865 ) { xsec = 0.207E-02; xsec_unc = 18.13 /100; return; }
  else if (glu_mass == 1870 ) { xsec = 0.201E-02; xsec_unc = 18.17 /100; return; }
  else if (glu_mass == 1875 ) { xsec = 0.196E-02; xsec_unc = 18.22 /100; return; }
  else if (glu_mass == 1880 ) { xsec = 0.191E-02; xsec_unc = 18.26 /100; return; }
  else if (glu_mass == 1885 ) { xsec = 0.186E-02; xsec_unc = 18.31 /100; return; }
  else if (glu_mass == 1890 ) { xsec = 0.181E-02; xsec_unc = 18.35 /100; return; }
  else if (glu_mass == 1895 ) { xsec = 0.176E-02; xsec_unc = 18.4 /100; return; }
  else if (glu_mass == 1900 ) { xsec = 0.171E-02; xsec_unc = 18.45 /100; return; }
  else if (glu_mass == 1905 ) { xsec = 0.167E-02; xsec_unc = 18.49 /100; return; }
  else if (glu_mass == 1910 ) { xsec = 0.163E-02; xsec_unc = 18.54 /100; return; }
  else if (glu_mass == 1915 ) { xsec = 0.158E-02; xsec_unc = 18.59 /100; return; }
  else if (glu_mass == 1920 ) { xsec = 0.154E-02; xsec_unc = 18.63 /100; return; }
  else if (glu_mass == 1925 ) { xsec = 0.150E-02; xsec_unc = 18.68 /100; return; }
  else if (glu_mass == 1930 ) { xsec = 0.146E-02; xsec_unc = 18.73 /100; return; }
  else if (glu_mass == 1935 ) { xsec = 0.142E-02; xsec_unc = 18.78 /100; return; }
  else if (glu_mass == 1940 ) { xsec = 0.139E-02; xsec_unc = 18.82 /100; return; }
  else if (glu_mass == 1945 ) { xsec = 0.135E-02; xsec_unc = 18.87 /100; return; }
  else if (glu_mass == 1950 ) { xsec = 0.131E-02; xsec_unc = 18.92 /100; return; }
  else if (glu_mass == 1955 ) { xsec = 0.128E-02; xsec_unc = 18.97 /100; return; }
  else if (glu_mass == 1960 ) { xsec = 0.125E-02; xsec_unc = 19.02 /100; return; }
  else if (glu_mass == 1965 ) { xsec = 0.121E-02; xsec_unc = 19.07 /100; return; }
  else if (glu_mass == 1970 ) { xsec = 0.118E-02; xsec_unc = 19.12 /100; return; }
  else if (glu_mass == 1975 ) { xsec = 0.115E-02; xsec_unc = 19.17 /100; return; }
  else if (glu_mass == 1980 ) { xsec = 0.112E-02; xsec_unc = 19.22 /100; return; }
  else if (glu_mass == 1985 ) { xsec = 0.109E-02; xsec_unc = 19.27 /100; return; }
  else if (glu_mass == 1990 ) { xsec = 0.106E-02; xsec_unc = 19.32 /100; return; }
  else if (glu_mass == 1995 ) { xsec = 0.104E-02; xsec_unc = 19.37 /100; return; }
  else if (glu_mass == 2000 ) { xsec = 0.101E-02; xsec_unc = 19.42 /100; return; }
  else if (glu_mass == 2005 ) { xsec = 0.983E-03; xsec_unc = 19.48 /100; return; }
  else if (glu_mass == 2010 ) { xsec = 0.957E-03; xsec_unc = 19.53 /100; return; }
  else if (glu_mass == 2015 ) { xsec = 0.933E-03; xsec_unc = 19.58 /100; return; }
  else if (glu_mass == 2020 ) { xsec = 0.908E-03; xsec_unc = 19.64 /100; return; }
  else if (glu_mass == 2025 ) { xsec = 0.885E-03; xsec_unc = 19.69 /100; return; }
  else if (glu_mass == 2030 ) { xsec = 0.862E-03; xsec_unc = 19.74 /100; return; }
  else if (glu_mass == 2035 ) { xsec = 0.840E-03; xsec_unc = 19.8 /100; return; }
  else if (glu_mass == 2040 ) { xsec = 0.818E-03; xsec_unc = 19.85 /100; return; }
  else if (glu_mass == 2045 ) { xsec = 0.797E-03; xsec_unc = 19.91 /100; return; }
  else if (glu_mass == 2050 ) { xsec = 0.776E-03; xsec_unc = 19.96 /100; return; }
  else if (glu_mass == 2055 ) { xsec = 0.756E-03; xsec_unc = 20.02 /100; return; }
  else if (glu_mass == 2060 ) { xsec = 0.737E-03; xsec_unc = 20.07 /100; return; }
  else if (glu_mass == 2065 ) { xsec = 0.718E-03; xsec_unc = 20.13 /100; return; }
  else if (glu_mass == 2070 ) { xsec = 0.699E-03; xsec_unc = 20.19 /100; return; }
  else if (glu_mass == 2075 ) { xsec = 0.681E-03; xsec_unc = 20.25 /100; return; }
  else if (glu_mass == 2080 ) { xsec = 0.664E-03; xsec_unc = 20.3 /100; return; }
  else if (glu_mass == 2085 ) { xsec = 0.647E-03; xsec_unc = 20.36 /100; return; }
  else if (glu_mass == 2090 ) { xsec = 0.630E-03; xsec_unc = 20.42 /100; return; }
  else if (glu_mass == 2095 ) { xsec = 0.614E-03; xsec_unc = 20.48 /100; return; }
  else if (glu_mass == 2100 ) { xsec = 0.598E-03; xsec_unc = 20.54 /100; return; }
  else if (glu_mass == 2105 ) { xsec = 0.583E-03; xsec_unc = 20.6 /100; return; }
  else if (glu_mass == 2110 ) { xsec = 0.568E-03; xsec_unc = 20.66 /100; return; }
  else if (glu_mass == 2115 ) { xsec = 0.553E-03; xsec_unc = 20.72 /100; return; }
  else if (glu_mass == 2120 ) { xsec = 0.539E-03; xsec_unc = 20.78 /100; return; }
  else if (glu_mass == 2125 ) { xsec = 0.525E-03; xsec_unc = 20.84 /100; return; }
  else if (glu_mass == 2130 ) { xsec = 0.512E-03; xsec_unc = 20.9 /100; return; }
  else if (glu_mass == 2135 ) { xsec = 0.499E-03; xsec_unc = 20.97 /100; return; }
  else if (glu_mass == 2140 ) { xsec = 0.486E-03; xsec_unc = 21.03 /100; return; }
  else if (glu_mass == 2145 ) { xsec = 0.473E-03; xsec_unc = 21.09 /100; return; }
  else if (glu_mass == 2150 ) { xsec = 0.461E-03; xsec_unc = 21.16 /100; return; }
  else if (glu_mass == 2155 ) { xsec = 0.449E-03; xsec_unc = 21.22 /100; return; }
  else if (glu_mass == 2160 ) { xsec = 0.438E-03; xsec_unc = 21.29 /100; return; }
  else if (glu_mass == 2165 ) { xsec = 0.427E-03; xsec_unc = 21.35 /100; return; }
  else if (glu_mass == 2170 ) { xsec = 0.416E-03; xsec_unc = 21.42 /100; return; }
  else if (glu_mass == 2175 ) { xsec = 0.405E-03; xsec_unc = 21.48 /100; return; }
  else if (glu_mass == 2180 ) { xsec = 0.395E-03; xsec_unc = 21.55 /100; return; }
  else if (glu_mass == 2185 ) { xsec = 0.385E-03; xsec_unc = 21.62 /100; return; }
  else if (glu_mass == 2190 ) { xsec = 0.375E-03; xsec_unc = 21.69 /100; return; }
  else if (glu_mass == 2195 ) { xsec = 0.365E-03; xsec_unc = 21.76 /100; return; }
  else if (glu_mass == 2200 ) { xsec = 0.356E-03; xsec_unc = 21.83 /100; return; }
  else if (glu_mass == 2205 ) { xsec = 0.347E-03; xsec_unc = 21.9 /100; return; }
  else if (glu_mass == 2210 ) { xsec = 0.338E-03; xsec_unc = 21.97 /100; return; }
  else if (glu_mass == 2215 ) { xsec = 0.330E-03; xsec_unc = 22.04 /100; return; }
  else if (glu_mass == 2220 ) { xsec = 0.321E-03; xsec_unc = 22.11 /100; return; }
  else if (glu_mass == 2225 ) { xsec = 0.313E-03; xsec_unc = 22.18 /100; return; }
  else if (glu_mass == 2230 ) { xsec = 0.305E-03; xsec_unc = 22.25 /100; return; }
  else if (glu_mass == 2235 ) { xsec = 0.297E-03; xsec_unc = 22.33 /100; return; }
  else if (glu_mass == 2240 ) { xsec = 0.290E-03; xsec_unc = 22.4 /100; return; }
  else if (glu_mass == 2245 ) { xsec = 0.283E-03; xsec_unc = 22.47 /100; return; }
  else if (glu_mass == 2250 ) { xsec = 0.275E-03; xsec_unc = 22.55 /100; return; }
  else if (glu_mass == 2255 ) { xsec = 0.268E-03; xsec_unc = 22.63 /100; return; }
  else if (glu_mass == 2260 ) { xsec = 0.262E-03; xsec_unc = 22.7 /100; return; }
  else if (glu_mass == 2265 ) { xsec = 0.255E-03; xsec_unc = 22.78 /100; return; }
  else if (glu_mass == 2270 ) { xsec = 0.248E-03; xsec_unc = 22.86 /100; return; }
  else if (glu_mass == 2275 ) { xsec = 0.242E-03; xsec_unc = 22.94 /100; return; }
  else if (glu_mass == 2280 ) { xsec = 0.236E-03; xsec_unc = 23.02 /100; return; }
  else if (glu_mass == 2285 ) { xsec = 0.230E-03; xsec_unc = 23.1 /100; return; }
  else if (glu_mass == 2290 ) { xsec = 0.224E-03; xsec_unc = 23.18 /100; return; }
  else if (glu_mass == 2295 ) { xsec = 0.219E-03; xsec_unc = 23.26 /100; return; }
  else if (glu_mass == 2300 ) { xsec = 0.213E-03; xsec_unc = 23.34 /100; return; }
  else if (glu_mass == 2305 ) { xsec = 0.208E-03; xsec_unc = 23.43 /100; return; }
  else if (glu_mass == 2310 ) { xsec = 0.202E-03; xsec_unc = 23.51 /100; return; }
  else if (glu_mass == 2315 ) { xsec = 0.197E-03; xsec_unc = 23.6 /100; return; }
  else if (glu_mass == 2320 ) { xsec = 0.192E-03; xsec_unc = 23.68 /100; return; }
  else if (glu_mass == 2325 ) { xsec = 0.187E-03; xsec_unc = 23.77 /100; return; }
  else if (glu_mass == 2330 ) { xsec = 0.183E-03; xsec_unc = 23.86 /100; return; }
  else if (glu_mass == 2335 ) { xsec = 0.178E-03; xsec_unc = 23.95 /100; return; }
  else if (glu_mass == 2340 ) { xsec = 0.174E-03; xsec_unc = 24.04 /100; return; }
  else if (glu_mass == 2345 ) { xsec = 0.169E-03; xsec_unc = 24.13 /100; return; }
  else if (glu_mass == 2350 ) { xsec = 0.165E-03; xsec_unc = 24.22 /100; return; }
  else if (glu_mass == 2355 ) { xsec = 0.161E-03; xsec_unc = 24.31 /100; return; }
  else if (glu_mass == 2360 ) { xsec = 0.157E-03; xsec_unc = 24.41 /100; return; }
  else if (glu_mass == 2365 ) { xsec = 0.153E-03; xsec_unc = 24.5 /100; return; }
  else if (glu_mass == 2370 ) { xsec = 0.149E-03; xsec_unc = 24.6 /100; return; }
  else if (glu_mass == 2375 ) { xsec = 0.145E-03; xsec_unc = 24.7 /100; return; }
  else if (glu_mass == 2380 ) { xsec = 0.142E-03; xsec_unc = 24.79 /100; return; }
  else if (glu_mass == 2385 ) { xsec = 0.138E-03; xsec_unc = 24.89 /100; return; }
  else if (glu_mass == 2390 ) { xsec = 0.134E-03; xsec_unc = 24.99 /100; return; }
  else if (glu_mass == 2395 ) { xsec = 0.131E-03; xsec_unc = 25.09 /100; return; }
  else if (glu_mass == 2400 ) { xsec = 0.128E-03; xsec_unc = 25.19 /100; return; }
  else if (glu_mass == 2405 ) { xsec = 0.125E-03; xsec_unc = 25.3 /100; return; }
  else if (glu_mass == 2410 ) { xsec = 0.121E-03; xsec_unc = 25.4 /100; return; }
  else if (glu_mass == 2415 ) { xsec = 0.118E-03; xsec_unc = 25.5 /100; return; }
  else if (glu_mass == 2420 ) { xsec = 0.115E-03; xsec_unc = 25.61 /100; return; }
  else if (glu_mass == 2425 ) { xsec = 0.113E-03; xsec_unc = 25.71 /100; return; }
  else if (glu_mass == 2430 ) { xsec = 0.110E-03; xsec_unc = 25.82 /100; return; }
  else if (glu_mass == 2435 ) { xsec = 0.107E-03; xsec_unc = 25.93 /100; return; }
  else if (glu_mass == 2440 ) { xsec = 0.104E-03; xsec_unc = 26.04 /100; return; }
  else if (glu_mass == 2445 ) { xsec = 0.102E-03; xsec_unc = 26.15 /100; return; }
  else if (glu_mass == 2450 ) { xsec = 0.991E-04; xsec_unc = 26.26 /100; return; }
  else if (glu_mass == 2455 ) { xsec = 0.966E-04; xsec_unc = 26.37 /100; return; }
  else if (glu_mass == 2460 ) { xsec = 0.941E-04; xsec_unc = 26.49 /100; return; }
  else if (glu_mass == 2465 ) { xsec = 0.918E-04; xsec_unc = 26.6 /100; return; }
  else if (glu_mass == 2470 ) { xsec = 0.895E-04; xsec_unc = 26.72 /100; return; }
  else if (glu_mass == 2475 ) { xsec = 0.872E-04; xsec_unc = 26.84 /100; return; }
  else if (glu_mass == 2480 ) { xsec = 0.850E-04; xsec_unc = 26.96 /100; return; }
  else if (glu_mass == 2485 ) { xsec = 0.829E-04; xsec_unc = 27.08 /100; return; }
  else if (glu_mass == 2490 ) { xsec = 0.808E-04; xsec_unc = 27.2 /100; return; }
  else if (glu_mass == 2495 ) { xsec = 0.788E-04; xsec_unc = 27.33 /100; return; }
  else if (glu_mass == 2500 ) { xsec = 0.768E-04; xsec_unc = 27.45 /100; return; }
  else if (glu_mass == 2505 ) { xsec = 0.749E-04; xsec_unc = 27.58 /100; return; }
  else if (glu_mass == 2510 ) { xsec = 0.730E-04; xsec_unc = 27.71 /100; return; }
  else if (glu_mass == 2515 ) { xsec = 0.712E-04; xsec_unc = 27.84 /100; return; }
  else if (glu_mass == 2520 ) { xsec = 0.694E-04; xsec_unc = 27.97 /100; return; }
  else if (glu_mass == 2525 ) { xsec = 0.677E-04; xsec_unc = 28.11 /100; return; }
  else if (glu_mass == 2530 ) { xsec = 0.660E-04; xsec_unc = 28.24 /100; return; }
  else if (glu_mass == 2535 ) { xsec = 0.643E-04; xsec_unc = 28.38 /100; return; }
  else if (glu_mass == 2540 ) { xsec = 0.627E-04; xsec_unc = 28.52 /100; return; }
  else if (glu_mass == 2545 ) { xsec = 0.611E-04; xsec_unc = 28.66 /100; return; }
  else if (glu_mass == 2550 ) { xsec = 0.596E-04; xsec_unc = 28.8 /100; return; }
  else if (glu_mass == 2555 ) { xsec = 0.581E-04; xsec_unc = 28.94 /100; return; }
  else if (glu_mass == 2560 ) { xsec = 0.566E-04; xsec_unc = 29.09 /100; return; }
  else if (glu_mass == 2565 ) { xsec = 0.552E-04; xsec_unc = 29.23 /100; return; }
  else if (glu_mass == 2570 ) { xsec = 0.538E-04; xsec_unc = 29.38 /100; return; }
  else if (glu_mass == 2575 ) { xsec = 0.525E-04; xsec_unc = 29.53 /100; return; }
  else if (glu_mass == 2580 ) { xsec = 0.512E-04; xsec_unc = 29.68 /100; return; }
  else if (glu_mass == 2585 ) { xsec = 0.499E-04; xsec_unc = 29.84 /100; return; }
  else if (glu_mass == 2590 ) { xsec = 0.486E-04; xsec_unc = 29.99 /100; return; }
  else if (glu_mass == 2595 ) { xsec = 0.474E-04; xsec_unc = 30.15 /100; return; }
  else if (glu_mass == 2600 ) { xsec = 0.462E-04; xsec_unc = 30.31 /100; return; }
  else if (glu_mass == 2605 ) { xsec = 0.451E-04; xsec_unc = 30.47 /100; return; }
  else if (glu_mass == 2610 ) { xsec = 0.439E-04; xsec_unc = 30.63 /100; return; }
  else if (glu_mass == 2615 ) { xsec = 0.428E-04; xsec_unc = 30.8 /100; return; }
  else if (glu_mass == 2620 ) { xsec = 0.418E-04; xsec_unc = 30.97 /100; return; }
  else if (glu_mass == 2625 ) { xsec = 0.407E-04; xsec_unc = 31.13 /100; return; }
  else if (glu_mass == 2630 ) { xsec = 0.397E-04; xsec_unc = 31.3 /100; return; }
  else if (glu_mass == 2635 ) { xsec = 0.387E-04; xsec_unc = 31.48 /100; return; }
  else if (glu_mass == 2640 ) { xsec = 0.377E-04; xsec_unc = 31.65 /100; return; }
  else if (glu_mass == 2645 ) { xsec = 0.368E-04; xsec_unc = 31.83 /100; return; }
  else if (glu_mass == 2650 ) { xsec = 0.359E-04; xsec_unc = 32.01 /100; return; }
  else if (glu_mass == 2655 ) { xsec = 0.350E-04; xsec_unc = 32.19 /100; return; }
  else if (glu_mass == 2660 ) { xsec = 0.341E-04; xsec_unc = 32.37 /100; return; }
  else if (glu_mass == 2665 ) { xsec = 0.332E-04; xsec_unc = 32.56 /100; return; }
  else if (glu_mass == 2670 ) { xsec = 0.324E-04; xsec_unc = 32.74 /100; return; }
  else if (glu_mass == 2675 ) { xsec = 0.316E-04; xsec_unc = 32.93 /100; return; }
  else if (glu_mass == 2680 ) { xsec = 0.308E-04; xsec_unc = 33.12 /100; return; }
  else if (glu_mass == 2685 ) { xsec = 0.300E-04; xsec_unc = 33.32 /100; return; }
  else if (glu_mass == 2690 ) { xsec = 0.293E-04; xsec_unc = 33.52 /100; return; }
  else if (glu_mass == 2695 ) { xsec = 0.285E-04; xsec_unc = 33.71 /100; return; }
  else if (glu_mass == 2700 ) { xsec = 0.278E-04; xsec_unc = 33.92 /100; return; }
  else if (glu_mass == 2705 ) { xsec = 0.271E-04; xsec_unc = 34.13 /100; return; }
  else if (glu_mass == 2710 ) { xsec = 0.265E-04; xsec_unc = 34.34 /100; return; }
  else if (glu_mass == 2715 ) { xsec = 0.258E-04; xsec_unc = 34.56 /100; return; }
  else if (glu_mass == 2720 ) { xsec = 0.251E-04; xsec_unc = 34.77 /100; return; }
  else if (glu_mass == 2725 ) { xsec = 0.245E-04; xsec_unc = 34.99 /100; return; }
  else if (glu_mass == 2730 ) { xsec = 0.239E-04; xsec_unc = 35.22 /100; return; }
  else if (glu_mass == 2735 ) { xsec = 0.233E-04; xsec_unc = 35.44 /100; return; }
  else if (glu_mass == 2740 ) { xsec = 0.227E-04; xsec_unc = 35.67 /100; return; }
  else if (glu_mass == 2745 ) { xsec = 0.221E-04; xsec_unc = 35.89 /100; return; }
  else if (glu_mass == 2750 ) { xsec = 0.216E-04; xsec_unc = 36.12 /100; return; }
  else if (glu_mass == 2755 ) { xsec = 0.211E-04; xsec_unc = 36.35 /100; return; }
  else if (glu_mass == 2760 ) { xsec = 0.205E-04; xsec_unc = 36.59 /100; return; }
  else if (glu_mass == 2765 ) { xsec = 0.200E-04; xsec_unc = 36.82 /100; return; }
  else if (glu_mass == 2770 ) { xsec = 0.195E-04; xsec_unc = 37.05 /100; return; }
  else if (glu_mass == 2775 ) { xsec = 0.190E-04; xsec_unc = 37.29 /100; return; }
  else if (glu_mass == 2780 ) { xsec = 0.185E-04; xsec_unc = 37.53 /100; return; }
  else if (glu_mass == 2785 ) { xsec = 0.181E-04; xsec_unc = 37.76 /100; return; }
  else if (glu_mass == 2790 ) { xsec = 0.176E-04; xsec_unc = 38.0 /100; return; }
  else if (glu_mass == 2795 ) { xsec = 0.172E-04; xsec_unc = 38.24 /100; return; }
  else if (glu_mass == 2800 ) { xsec = 0.168E-04; xsec_unc = 38.48 /100; return; }
  else if (glu_mass == 2805 ) { xsec = 0.163E-04; xsec_unc = 38.72 /100; return; }
  else if (glu_mass == 2810 ) { xsec = 0.159E-04; xsec_unc = 38.95 /100; return; }
  else if (glu_mass == 2815 ) { xsec = 0.155E-04; xsec_unc = 39.19 /100; return; }
  else if (glu_mass == 2820 ) { xsec = 0.151E-04; xsec_unc = 39.43 /100; return; }
  else if (glu_mass == 2825 ) { xsec = 0.148E-04; xsec_unc = 39.66 /100; return; }
  else if (glu_mass == 2830 ) { xsec = 0.144E-04; xsec_unc = 39.9 /100; return; }
  else if (glu_mass == 2835 ) { xsec = 0.140E-04; xsec_unc = 40.14 /100; return; }
  else if (glu_mass == 2840 ) { xsec = 0.137E-04; xsec_unc = 40.38 /100; return; }
  else if (glu_mass == 2845 ) { xsec = 0.133E-04; xsec_unc = 40.62 /100; return; }
  else if (glu_mass == 2850 ) { xsec = 0.130E-04; xsec_unc = 40.86 /100; return; }
  else if (glu_mass == 2855 ) { xsec = 0.127E-04; xsec_unc = 41.1 /100; return; }
  else if (glu_mass == 2860 ) { xsec = 0.124E-04; xsec_unc = 41.34 /100; return; }
  else if (glu_mass == 2865 ) { xsec = 0.121E-04; xsec_unc = 41.59 /100; return; }
  else if (glu_mass == 2870 ) { xsec = 0.118E-04; xsec_unc = 41.83 /100; return; }
  else if (glu_mass == 2875 ) { xsec = 0.115E-04; xsec_unc = 42.07 /100; return; }
  else if (glu_mass == 2880 ) { xsec = 0.112E-04; xsec_unc = 42.32 /100; return; }
  else if (glu_mass == 2885 ) { xsec = 0.109E-04; xsec_unc = 42.56 /100; return; }
  else if (glu_mass == 2890 ) { xsec = 0.106E-04; xsec_unc = 42.81 /100; return; }
  else if (glu_mass == 2895 ) { xsec = 0.104E-04; xsec_unc = 43.05 /100; return; }
  else if (glu_mass == 2900 ) { xsec = 0.101E-04; xsec_unc = 43.3 /100; return; }
  else if (glu_mass == 2905 ) { xsec = 0.986E-05; xsec_unc = 43.55 /100; return; }
  else if (glu_mass == 2910 ) { xsec = 0.961E-05; xsec_unc = 43.79 /100; return; }
  else if (glu_mass == 2915 ) { xsec = 0.937E-05; xsec_unc = 44.04 /100; return; }
  else if (glu_mass == 2920 ) { xsec = 0.914E-05; xsec_unc = 44.29 /100; return; }
  else if (glu_mass == 2925 ) { xsec = 0.891E-05; xsec_unc = 44.54 /100; return; }
  else if (glu_mass == 2930 ) { xsec = 0.869E-05; xsec_unc = 44.79 /100; return; }
  else if (glu_mass == 2935 ) { xsec = 0.848E-05; xsec_unc = 45.04 /100; return; }
  else if (glu_mass == 2940 ) { xsec = 0.827E-05; xsec_unc = 45.29 /100; return; }
  else if (glu_mass == 2945 ) { xsec = 0.806E-05; xsec_unc = 45.54 /100; return; }
  else if (glu_mass == 2950 ) { xsec = 0.786E-05; xsec_unc = 45.8 /100; return; }
  else if (glu_mass == 2955 ) { xsec = 0.767E-05; xsec_unc = 46.05 /100; return; }
  else if (glu_mass == 2960 ) { xsec = 0.748E-05; xsec_unc = 46.3 /100; return; }
  else if (glu_mass == 2965 ) { xsec = 0.729E-05; xsec_unc = 46.56 /100; return; }
  else if (glu_mass == 2970 ) { xsec = 0.711E-05; xsec_unc = 46.81 /100; return; }
  else if (glu_mass == 2975 ) { xsec = 0.694E-05; xsec_unc = 47.07 /100; return; }
  else if (glu_mass == 2980 ) { xsec = 0.677E-05; xsec_unc = 47.32 /100; return; }
  else if (glu_mass == 2985 ) { xsec = 0.660E-05; xsec_unc = 47.58 /100; return; }
  else if (glu_mass == 2990 ) { xsec = 0.644E-05; xsec_unc = 47.84 /100; return; }
  else if (glu_mass == 2995 ) { xsec = 0.628E-05; xsec_unc = 48.09 /100; return; }
  else if (glu_mass == 3000 ) { xsec = 0.612E-05; xsec_unc = 48.35 /100; return; }
  else {xsec = 0.; xsec_unc = 0.;}
}
