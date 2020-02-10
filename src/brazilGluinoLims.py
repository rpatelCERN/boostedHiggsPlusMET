#! /usr/bin/env python

import os
import glob
import math
import array
import sys
import time
import ROOT
from array import array

from ROOT import gROOT
gROOT.SetBatch(True)


import tdrstyle
tdrstyle.setTDRStyle()

## ===========================================================================================
## ===========================================================================================
## ===========================================================================================

def columnToList(fn,col):
	f = open(fn,'r');

	olist = [];
	for line in f:
		linelist = line.strip().split()
		olist.append( linelist[col] );
	return olist

def ExtractFile(iname, tag):
	f = ROOT.TFile(iname);
	t = f.Get("limit");
	lims = [];
	lims.append(tag);
	for i in range(6):
		t.GetEntry(i);
		lims.append( t.limit )
	return lims;

if __name__ == '__main__':

	idir = "./";
	results = []; results2 = []; results3 = []; results4 = [];

	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH750.AsymptoticLimits.mH750.root','750') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH1000.AsymptoticLimits.mH1000.root','1000') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH1100.AsymptoticLimits.mH1100.root','1100') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH1200.AsymptoticLimits.mH1200.root','1200') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH1300.AsymptoticLimits.mH1300.root','1300') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH1400.AsymptoticLimits.mH1400.root','1400') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH1500.AsymptoticLimits.mH1500.root','1500') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH1600.AsymptoticLimits.mH1600.root','1600') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH1700.AsymptoticLimits.mH1700.root','1700') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH1800.AsymptoticLimits.mH1800.root','1800') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH1900.AsymptoticLimits.mH1900.root','1900') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH2000.AsymptoticLimits.mH2000.root','2000') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH2100.AsymptoticLimits.mH2100.root','2100') );
	results.append( ExtractFile(idir+'/higgsCombineGluinoResolvedOnlyT5HH2200.AsymptoticLimits.mH2200.root','2200') );

	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH750.AsymptoticLimits.mH750.root','750') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1000.AsymptoticLimits.mH1000.root','1000') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1100.AsymptoticLimits.mH1100.root','1100') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1200.AsymptoticLimits.mH1200.root','1200') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1300.AsymptoticLimits.mH1300.root','1300') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1400.AsymptoticLimits.mH1400.root','1400') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1500.AsymptoticLimits.mH1500.root','1500') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1600.AsymptoticLimits.mH1600.root','1600') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1700.AsymptoticLimits.mH1700.root','1700') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1800.AsymptoticLimits.mH1800.root','1800') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1900.AsymptoticLimits.mH1900.root','1900') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH2000.AsymptoticLimits.mH2000.root','2000') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH2100.AsymptoticLimits.mH2100.root','2100') );
	# results.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH2200.AsymptoticLimits.mH2200.root','2200') );

	#higgsCombineGluinoBoostedOnlyT5HH2200
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH750.AsymptoticLimits.mH750.root','750') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1000.AsymptoticLimits.mH1000.root','1000') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1100.AsymptoticLimits.mH1100.root','1100') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1200.AsymptoticLimits.mH1200.root','1200') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1300.AsymptoticLimits.mH1300.root','1300') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1400.AsymptoticLimits.mH1400.root','1400') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1500.AsymptoticLimits.mH1500.root','1500') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1600.AsymptoticLimits.mH1600.root','1600') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1700.AsymptoticLimits.mH1700.root','1700') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1800.AsymptoticLimits.mH1800.root','1800') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH1900.AsymptoticLimits.mH1900.root','1900') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH2000.AsymptoticLimits.mH2000.root','2000') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH2100.AsymptoticLimits.mH2100.root','2100') );
	results2.append( ExtractFile(idir+'/higgsCombineGluinoBoostedOnlyT5HH2200.AsymptoticLimits.mH2200.root','2200') );

	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH750.AsymptoticLimits.mH750.root','750') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH1000.AsymptoticLimits.mH1000.root','1000') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH1100.AsymptoticLimits.mH1100.root','1100') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH1200.AsymptoticLimits.mH1200.root','1200') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH1300.AsymptoticLimits.mH1300.root','1300') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH1400.AsymptoticLimits.mH1400.root','1400') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH1500.AsymptoticLimits.mH1500.root','1500') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH1600.AsymptoticLimits.mH1600.root','1600') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH1700.AsymptoticLimits.mH1700.root','1700') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH1800.AsymptoticLimits.mH1800.root','1800') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH1900.AsymptoticLimits.mH1900.root','1900') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH2000.AsymptoticLimits.mH2000.root','2000') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH2100.AsymptoticLimits.mH2100.root','2100') );
	results3.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostOnlyT5HH2200.AsymptoticLimits.mH2200.root','2200') );

	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH750.AsymptoticLimits.mH750.root','750') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH1000.AsymptoticLimits.mH1000.root','1000') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH1100.AsymptoticLimits.mH1100.root','1100') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH1200.AsymptoticLimits.mH1200.root','1200') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH1300.AsymptoticLimits.mH1300.root','1300') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH1400.AsymptoticLimits.mH1400.root','1400') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH1500.AsymptoticLimits.mH1500.root','1500') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH1600.AsymptoticLimits.mH1600.root','1600') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH1700.AsymptoticLimits.mH1700.root','1700') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH1800.AsymptoticLimits.mH1800.root','1800') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH1900.AsymptoticLimits.mH1900.root','1900') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH2000.AsymptoticLimits.mH2000.root','2000') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH2100.AsymptoticLimits.mH2100.root','2100') );
	results4.append( ExtractFile(idir+'/higgsCombineGluinoResToBoostFullComboT5HH2200.AsymptoticLimits.mH2200.root','2200') );


	#Change this
	# xsecs=[	0.277E+01, 0.385E+00,0.191E+00, 0.985E-01, 0.522E-01, 0.284E-01, 0.157E-01, 0.887E-02, 0.507E-02, 0.293E-02, 0.171E-02, 0.101E-02, 0.598E-03, 0.356E-03]
	xsecs=[2.26585,0.325388, 0.163491, 0.0856418, 0.0460525, 0.0252977, 0.0141903, 0.00810078, 0.00470323, 0.00276133, 0.00163547, 0.000981077, 0.000591918,0.000359318]

	names   = []; names2   = []; names3   = []; names4   = [];
	l_obs   = []; l_obs2   = []; l_obs3   = []; l_obs4   = [];
	l_m2sig = []; l_m2sig2 = []; l_m2sig3 = []; l_m2sig4 = [];
	l_m1sig = []; l_m1sig2 = []; l_m1sig3 = []; l_m1sig4 = [];
	l_exp   = []; l_exp2   = []; l_exp3   = []; l_exp4   = [];
	l_p1sig = []; l_p1sig2 = []; l_p1sig3 = []; l_p1sig4 = [];
	l_p2sig = []; l_p2sig2 = []; l_p2sig3 = []; l_p2sig4 = [];
	count=0; count2=0; count3=0; count4=0;
	BR=1.0 #not needed for gluino since it is fully simulated
	for r in results:
		names.append(r[0]);
		l_m2sig.append(r[1]*xsecs[count]*BR);
		l_m1sig.append(r[2]*xsecs[count]*BR);
		l_exp.append(r[3]*xsecs[count]*BR);
		l_p1sig.append(r[4]*xsecs[count]*BR);
		l_p2sig.append(r[5]*xsecs[count]*BR);
		l_obs.append(r[6]*xsecs[count]*BR);
		count=count+1
	# print "l_exp = ", l_exp
	# print "l_obs = ", l_obs

	for r in results2:
		names2.append(r[0]);
		l_m2sig2.append(r[1]*xsecs[count2]*BR);
		l_m1sig2.append(r[2]*xsecs[count2]*BR);
		l_exp2.append(r[3]*xsecs[count2]*BR);
		l_p1sig2.append(r[4]*xsecs[count2]*BR);
		l_p2sig2.append(r[5]*xsecs[count2]*BR);
		l_obs2.append(r[6]*xsecs[count2]*BR);
		count2=count2+1

	for r in results3:
		names3.append(r[0]);
		l_m2sig3.append(r[1]*xsecs[count3]*BR);
		l_m1sig3.append(r[2]*xsecs[count3]*BR);
		l_exp3.append(r[3]*xsecs[count3]*BR);
		l_p1sig3.append(r[4]*xsecs[count3]*BR);
		l_p2sig3.append(r[5]*xsecs[count3]*BR);
		l_obs3.append(r[6]*xsecs[count3]*BR);
		count3=count3+1

	for r in results4:
		names4.append(r[0]);
		l_m2sig4.append(r[1]*xsecs[count4]*BR);
		l_m1sig4.append(r[2]*xsecs[count4]*BR);
		l_exp4.append(r[3]*xsecs[count4]*BR);
		l_p1sig4.append(r[4]*xsecs[count4]*BR);
		l_p2sig4.append(r[5]*xsecs[count4]*BR);
		l_obs4.append(r[6]*xsecs[count4]*BR);
		count4=count4+1


	a_xax = array('d', []); a_xax2 = array('d', []); a_xax3 = array('d', []); a_xax4 = array('d', []);
	a2_xax = array('d', []);
	a_exp = array('d', []);a_exp2 = array('d', []); a_exp3 = array('d', []); a_exp4 = array('d', []);
	a_obs = array('d', []);
	a_1sig = array('d', []);
	a_2sig = array('d', []);
	#Need to do this a bit more clever
	for i in range(len(names)): a_xax.append( float(names[i]) );
	for i in range(len(names)): a2_xax.append( float(names[i]) );
	for i in range(len(names)-1,-1,-1): a2_xax.append( float(names[i]));
	for i in range(len(l_obs)): a_obs.append( float(l_obs[i]) );
	for i in range(len(l_exp)): a_exp.append( float(l_exp[i]) );
	for i in range(len(l_m2sig)): a_2sig.append( float(l_m2sig[i]) );
	for i in range(len(l_p2sig)-1,-1,-1): a_2sig.append( float(l_p2sig[i]) );
	for i in range(len(l_m1sig)): a_1sig.append( float(l_m1sig[i]) );
	for i in range(len(l_p1sig)-1,-1,-1): a_1sig.append( float(l_p1sig[i]) );


	for i in range(len(names2)): a_xax2.append( float(names2[i]) );
	for i in range(len(l_exp2)): a_exp2.append( float(l_exp2[i]) );

	for i in range(len(names3)): a_xax3.append( float(names3[i]) );
	for i in range(len(l_exp3)): a_exp3.append( float(l_exp3[i]) );

	for i in range(len(names4)): a_xax4.append( float(names4[i]) );
	for i in range(len(l_exp4)): a_exp4.append( float(l_exp4[i]) );

	a_2sig.append(results[0][6])
	a2_xax.append(0.5)

	g_exp = ROOT.TGraph(len(a_xax), a_xax, a_exp)
	g_obs = ROOT.TGraph(len(a_xax), a_xax, a_obs)
	g_1sig = ROOT.TGraph(len(2*a_xax), a2_xax, a_1sig)
	g_2sig = ROOT.TGraph(len(2*a_xax), a2_xax, a_2sig)


	g_exp2 = ROOT.TGraph(len(a_xax2), a_xax2, a_exp2)
	g_exp3 = ROOT.TGraph(len(a_xax3), a_xax3, a_exp3)
	g_exp4 = ROOT.TGraph(len(a_xax4), a_xax4, a_exp4)

	can = ROOT.TCanvas("can","can",1800,1200);
	hrl = ROOT.TH1F("hrl","hrl",36,700,2500);

	# hrl = can.DrawFrame(0,0,6,15);
	hrl.GetXaxis().SetTitle("Gluino mass m_{#tilde{g}} [GeV]");
	hrl.GetXaxis().SetRangeUser(750,2200);
	hrl.GetXaxis().SetLabelSize(0.035)
	hrl.GetXaxis().SetTitleSize(0.04)
	hrl.GetXaxis().SetTitleOffset(1.2)
	hrl.GetYaxis().SetTitle("#sigma_{95% CL} [pb] ");
	#hrl.GetYaxis().SetTitleOffset(0.2)
	hrl.GetYaxis().SetTitleSize(0.04)

	#for i in range(0,15):
	#	if i%3==0:
	#		hrl.GetXaxis().SetBinLabel(i+1,names[i])
	#	if i==14:hrl.GetXaxis().SetBinLabel(i+1, names[i])
	#hrl.GetXaxis().SetBinLabel(2,names[1])
	#hrl.GetXaxis().SetBinLabel(3,names[2])
	#hrl.GetXaxis().SetBinLabel(4,names[3])
	#hrl.GetXaxis().SetBinLabel(5,names[4])
	#hrl.GetXaxis().SetBinLabel(6,names[5])
	#hrl.SetMaximum(1000.);
	hrl.Draw();
	hrl.GetYaxis().SetRangeUser(0.00001,10.);

	#hrl.GetYaxis().SetTitle("UL on #sigma_{95\% CL} ")#x BR(hh#rightarrow bbbb) [fb] ");
	can.SetGrid();
	can.SetLogy();
	txta = ROOT.TLatex(0.2,0.95,"CMS");
	txta.SetNDC();
	txtb = ROOT.TLatex(0.27,0.95,"Preliminary");
	txtb.SetNDC(); txtb.SetTextFont(52); txtb.SetTextSize(0.042);
	txtc = ROOT.TLatex(0.75,0.96,"137 fb^{-1} (13 TeV)");
	# txtc = ROOT.TLatex(0.75,0.96,"35.9 fb^{-1} (13 TeV)");
	#txtc = ROOT.TLatex(0.75,0.96,"150 fb^{-1} (13 TeV)");
	txtc.SetNDC(); txtc.SetTextFont(42); txtc.SetTextSize(0.04);
	txtd = ROOT.TLatex(0.6,0.88,"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}^{2}_{0},  #tilde{#chi}^{2}_{0} #rightarrow H #tilde{#chi}^{1}_{0}");
	txtd.SetNDC(); txtd.SetTextFont(42); txtd.SetTextSize(0.05);
	#f=open("LatestXsecGluGlu.txt", 'r')
	a_stop = array('d', []);
	a_xsec = array('d', []);
	count=0;
	for x in xsecs:
		a_stop.append(float(names[count]));
		# a_xsec.append(x*0.57*0.57)
		a_xsec.append(x)
		count=count+1
	g_xsec=ROOT.TGraph(len(a_stop), a_stop, a_xsec)
	leg = ROOT.TLegend(0.65,0.60,0.9,0.8);
	# leg = ROOT.TLegend(0.15,0.20,0.4,0.4);
	leg.SetFillStyle(1001);
	leg.SetFillColor(0);
	leg.SetBorderSize(1);
	# leg.SetNColumns(2);

	#For Full combo limits
	# leg.AddEntry(g_exp,"expected","l") #E commented out
	# leg.AddEntry(g_obs,"observed","l")#E commented out
	# leg.AddEntry(g_2sig,"expected 2#sigma","f") #E commented out
	# leg.AddEntry(g_1sig,"expected 1#sigma","f") #E commented out
	#

	#For overlaid
	leg.AddEntry(g_exp, "Resolved Only","l")
	leg.AddEntry(g_exp2,"Boosted Only","l")
	leg.AddEntry(g_exp3,"Resolved+Boost","l")
	# leg.AddEntry(g_exp4,"Full Combo","l")



	leg.AddEntry(g_xsec, "Theory cross-section", "l")


	g_1sig.SetFillColor(ROOT.kGreen);
	g_1sig.SetFillStyle(3244);
	g_2sig.SetFillColor(ROOT.kYellow);
	g_2sig.SetFillStyle(3244);
	g_exp.SetLineStyle(2);
	g_exp.SetLineColor(ROOT.kBlack);
	g_exp.SetLineWidth(2);
	g_obs.SetLineWidth(2);

	g_exp.SetLineColor(ROOT.kCyan+2); g_exp.SetLineWidth(2);
	g_exp2.SetLineColor(ROOT.kRed); g_exp2.SetLineWidth(2);
	g_exp3.SetLineColor(ROOT.kBlue); g_exp3.SetLineWidth(2);
	# g_exp4.SetLineColor(ROOT.kBlack); g_exp4.SetLineWidth(2);


	# g_2sig.Draw('f'); #E commented out
	# g_1sig.Draw('fsames'); #E commented out
	#g_1sig.Draw('f');
	# g_obs.Draw('lsames');
	# g_exp.Draw('lsames'); #E commented out


	g_exp.Draw('l');
	g_exp2.Draw('lsames');
	g_exp3.Draw('lsames');
	# g_exp4.Draw('lsames');
	# for i in range(0,100):
	# 	if(g_xsec.Eval(900+i)-g_exp2.Eval(900+i)>0):
	# 		print "Mass %d  Exp Excl %g " %(900+i,g_exp2.Eval(900+i))
	# 		print "Theory Xsec %g " %g_xsec.Eval(900+i)
	#oneLine.Draw("LSAMES");
	txta.Draw();
	txtb.Draw();
	txtc.Draw();
	txtd.Draw();
	leg.Draw();
	g_xsec.SetLineStyle(2);
	g_xsec.SetLineWidth(2);
	# g_xsec.SetLineColor(ROOT.kBlue);
	g_xsec.SetLineColor(ROOT.kBlack);
	g_xsec.Draw("lsame")
	# can.SaveAs('T5HHResults_BoostToResFull.pdf');
	# can.SaveAs('T5HHResults_FullCombo_overlaid.pdf');
	# can.SaveAs('T5HHResults_BoostOnly.pdf');
	can.SaveAs('T5HHResults_ResToBoost.pdf');
