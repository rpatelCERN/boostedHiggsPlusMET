import os
from ROOT import *
import sys


V18Signal_DIR = "/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV18/scan/tree_signal/"
fNames=open("higgsino2DFileNames.txt", 'r');
fAlpha=TFile("ALPHABET_V18_2Dsignal_2BoostedH.root","r")

SignifScan=TH2F("SignifScan","SignifScan",41,200,1225,42,50,1100)
mGo=[]; mLsp=[]; limit=[];

for line in fNames:
    x = line.split('_')
    hino_mass = x[5]
    LSP_mass = x[6]
    print("For mNLSP=%s and mLSP=%s" %(hino_mass,LSP_mass))
    f1=TFile(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_2D_"+hino_mass+"_"+LSP_mass+"_MC2016_fast.root", "r")
    f2=TFile(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_2D_"+hino_mass+"_"+LSP_mass+"_MC2017_fast.root", "r")
    f3=TFile(V18Signal_DIR+"tree_TChiHH_HToBB_HToBB_2D_"+hino_mass+"_"+LSP_mass+"_MC2018_fast.root", "r")
    tree1=f1.Get("tree"); tree2=f2.Get("tree"); tree3=f3.Get("tree");
    tree1.GetEntry(0); weight=getattr(tree1,"Weight"); #this is the same for all 3 years
    procHisto1=f1.Get("nEventProc"); procHisto2=f2.Get("nEventProc"); procHisto3=f3.Get("nEventProc");
    totalEvt1_prescale = procHisto1.GetBinContent(1)
    totalEvt2_prescale = procHisto2.GetBinContent(1)
    totalEvt3_prescale = procHisto3.GetBinContent(1)
    procHisto1.Scale(weight*35922.0*0.5823329*0.5823329/totalEvt1_prescale);
    procHisto2.Scale(weight*41529.0*0.5823329*0.5823329/totalEvt2_prescale);
    procHisto3.Scale(weight*59740.0*0.5823329*0.5823329/totalEvt3_prescale);
    totalEvt1 = procHisto1.GetBinContent(1)
    totalEvt2 = procHisto2.GetBinContent(1)
    totalEvt3 = procHisto3.GetBinContent(1)
    totalEvt = totalEvt1+totalEvt2+totalEvt3

    signal=fAlpha.Get("MET_doubletagSR_TChiHH%s_LSP%s" %(hino_mass,LSP_mass));
    passCuts = float(signal.Integral(1,5))

    f1.Close(); f2.Close(); f3.Close()
    limitEff = passCuts/totalEvt

    #Add 4 to account for some of the weird masses and border edges
    thisBin = SignifScan.FindBin(Double(hino_mass)+4.0,Double(LSP_mass)+4.0)
    SignifScan.SetBinContent(thisBin,limitEff)


canv = TCanvas("","",1200,950)
canv.cd()
canv.SetRightMargin(0.15)
SignifScan.SetStats(0)
SignifScan.SetTitle('')
SignifScan.SetMinimum(0.000001);
SignifScan.SetMaximum(1.0);
xparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{2}}}";
yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
SignifScan.GetXaxis().SetTitle("m_{"+xparticle+"} [GeV]")
SignifScan.GetYaxis().SetTitle("m_{"+yparticle+"} [GeV]")
SignifScan.GetZaxis().SetTitle("Efficiency")
SignifScan.GetZaxis().SetTitleOffset(1.3);
SignifScan.GetYaxis().SetTitleOffset(1.3);
SignifScan.GetXaxis().SetTitleOffset(1.3);

ltitle = TLatex(canv.GetLeftMargin(), 1.-0.7*canv.GetTopMargin(),"#font[62]{CMS}#scale[0.76]{#font[52]{ Supplementary}}")
rtitle = TLatex(1.-canv.GetRightMargin(), 1.-0.7*canv.GetTopMargin(),"#scale[0.8]{137 fb^{-1} (13 TeV)}")
ltitle.SetNDC(); rtitle.SetNDC()
ltitle.SetTextAlign(12); rtitle.SetTextAlign(32)

canv.SetLogz()
SignifScan.Draw("colz")
ltitle.Draw("same")
rtitle.Draw("same")
canv.Print("test.pdf")
SignifScan.SaveAs("graph.root")



thisBin = SignifScan.FindBin(Double(hino_mass)+4.0,Double(LSP_mass)+4.0)
thisEff_1200_950 = SignifScan.GetBinContent(SignifScan.FindBin(1205,955))
print "Efficiency at 1200,1000: %d" %(thisEff_1200_950)

fNames.close()
fAlpha.Close()
