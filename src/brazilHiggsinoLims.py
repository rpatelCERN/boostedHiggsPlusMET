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
	results = []; results2 = []; results3 = []; results4 = []; results5 = []; results_hino = [];

	#BoostedOnly
	results.append( ExtractFile(idir+'/higgsCombineTChiHH200_LSP1_2BoostedH.AsymptoticLimits.mH120.root','200') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH225_LSP1_2BoostedH.AsymptoticLimits.mH120.root','225') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH250_LSP1_2BoostedH.AsymptoticLimits.mH120.root','250') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH275_LSP1_2BoostedH.AsymptoticLimits.mH120.root','275') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH300_LSP1_2BoostedH.AsymptoticLimits.mH120.root','300') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH325_LSP1_2BoostedH.AsymptoticLimits.mH120.root','325') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH350_LSP1_2BoostedH.AsymptoticLimits.mH120.root','350') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH375_LSP1_2BoostedH.AsymptoticLimits.mH120.root','375') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH400_LSP1_2BoostedH.AsymptoticLimits.mH120.root','400') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH425_LSP1_2BoostedH.AsymptoticLimits.mH120.root','425') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH450_LSP1_2BoostedH.AsymptoticLimits.mH120.root','450') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH475_LSP1_2BoostedH.AsymptoticLimits.mH120.root','475') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH500_LSP1_2BoostedH.AsymptoticLimits.mH120.root','500') );
	# results.append( ExtractFile(idir+'/higgsCombineTChiHH525_LSP1_2BoostedH.AsymptoticLimits.mH120.root','525') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH550_LSP1_2BoostedH.AsymptoticLimits.mH120.root','550') );
	# results.append( ExtractFile(idir+'/higgsCombineTChiHH575_LSP1_2BoostedH.AsymptoticLimits.mH120.root','575') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH600_LSP1_2BoostedH.AsymptoticLimits.mH120.root','600') );
	# results.append( ExtractFile(idir+'/higgsCombineTChiHH625_LSP1_2BoostedH.AsymptoticLimits.mH120.root','625') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH650_LSP1_2BoostedH.AsymptoticLimits.mH120.root','650') );
	# results.append( ExtractFile(idir+'/higgsCombineTChiHH675_LSP1_2BoostedH.AsymptoticLimits.mH120.root','675') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH700_LSP1_2BoostedH.AsymptoticLimits.mH120.root','700') );
	# results.append( ExtractFile(idir+'/higgsCombineTChiHH725_LSP1_2BoostedH.AsymptoticLimits.mH120.root','725') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH750_LSP1_2BoostedH.AsymptoticLimits.mH120.root','750') );
	# results.append( ExtractFile(idir+'/higgsCombineTChiHH775_LSP1_2BoostedH.AsymptoticLimits.mH120.root','775') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH800_LSP1_2BoostedH.AsymptoticLimits.mH120.root','800') );
	# results.append( ExtractFile(idir+'/higgsCombineTChiHH825_LSP1_2BoostedH.AsymptoticLimits.mH120.root','825') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH850_LSP1_2BoostedH.AsymptoticLimits.mH120.root','850') );
	# results.append( ExtractFile(idir+'/higgsCombineTChiHH875_LSP1_2BoostedH.AsymptoticLimits.mH120.root','875') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH900_LSP1_2BoostedH.AsymptoticLimits.mH120.root','900') );
	# results.append( ExtractFile(idir+'/higgsCombineTChiHH925_LSP1_2BoostedH.AsymptoticLimits.mH120.root','925') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH950_LSP1_2BoostedH.AsymptoticLimits.mH120.root','950') );
	# results.append( ExtractFile(idir+'/higgsCombineTChiHH975_LSP1_2BoostedH.AsymptoticLimits.mH120.root','975') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH1000_LSP1_2BoostedH.AsymptoticLimits.mH120.root','1000') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH1100_LSP1_2BoostedH.AsymptoticLimits.mH120.root','1100') );
	results.append( ExtractFile(idir+'/higgsCombineTChiHH1200_LSP1_2BoostedH.AsymptoticLimits.mH120.root','1200') );


	# #Combo (ResolvedOnly)
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH200ResOnly.AsymptoticLimits.mH120.root','200') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH225ResOnly.AsymptoticLimits.mH120.root','225') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH250ResOnly.AsymptoticLimits.mH120.root','250') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH275ResOnly.AsymptoticLimits.mH120.root','275') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH300ResOnly.AsymptoticLimits.mH120.root','300') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH325ResOnly.AsymptoticLimits.mH120.root','325') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH350ResOnly.AsymptoticLimits.mH120.root','350') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH375ResOnly.AsymptoticLimits.mH120.root','375') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH400ResOnly.AsymptoticLimits.mH120.root','400') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH425ResOnly.AsymptoticLimits.mH120.root','425') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH450ResOnly.AsymptoticLimits.mH120.root','450') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH475ResOnly.AsymptoticLimits.mH120.root','475') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH500ResOnly.AsymptoticLimits.mH120.root','500') );
	# results2.append( ExtractFile(idir+'/higgsCombineTChiHH525ResOnly.AsymptoticLimits.mH120.root','525') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH550ResOnly.AsymptoticLimits.mH120.root','550') );
	# results2.append( ExtractFile(idir+'/higgsCombineTChiHH575ResOnly.AsymptoticLimits.mH120.root','575') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH600ResOnly.AsymptoticLimits.mH120.root','600') );
	# results2.append( ExtractFile(idir+'/higgsCombineTChiHH625ResOnly.AsymptoticLimits.mH120.root','625') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH650ResOnly.AsymptoticLimits.mH120.root','650') );
	# results2.append( ExtractFile(idir+'/higgsCombineTChiHH675ResOnly.AsymptoticLimits.mH120.root','675') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH700ResOnly.AsymptoticLimits.mH120.root','700') );
	# results2.append( ExtractFile(idir+'/higgsCombineTChiHH725ResOnly.AsymptoticLimits.mH120.root','725') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH750ResOnly.AsymptoticLimits.mH120.root','750') );
	# results2.append( ExtractFile(idir+'/higgsCombineTChiHH775ResOnly.AsymptoticLimits.mH120.root','775') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH800ResOnly.AsymptoticLimits.mH120.root','800') );
	# results2.append( ExtractFile(idir+'/higgsCombineTChiHH825ResOnly.AsymptoticLimits.mH120.root','825') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH850ResOnly.AsymptoticLimits.mH120.root','850') );
	# results2.append( ExtractFile(idir+'/higgsCombineTChiHH875ResOnly.AsymptoticLimits.mH120.root','875') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH900ResOnly.AsymptoticLimits.mH120.root','900') );
	# results2.append( ExtractFile(idir+'/higgsCombineTChiHH925ResOnly.AsymptoticLimits.mH120.root','925') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH950ResOnly.AsymptoticLimits.mH120.root','950') );
	# results2.append( ExtractFile(idir+'/higgsCombineTChiHH975ResOnly.AsymptoticLimits.mH120.root','975') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH1000ResOnly.AsymptoticLimits.mH120.root','1000') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH1100ResOnly.AsymptoticLimits.mH120.root','1100') );
	results2.append( ExtractFile(idir+'/higgsCombineTChiHH1200ResOnly.AsymptoticLimits.mH120.root','1200') );


	# #Combo (ResolvedOnly+BoostedVeto)
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH200Combo.AsymptoticLimits.mH120.root','200') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH225Combo.AsymptoticLimits.mH120.root','225') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH250Combo.AsymptoticLimits.mH120.root','250') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH275Combo.AsymptoticLimits.mH120.root','275') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH300Combo.AsymptoticLimits.mH120.root','300') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH325Combo.AsymptoticLimits.mH120.root','325') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH350Combo.AsymptoticLimits.mH120.root','350') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH375Combo.AsymptoticLimits.mH120.root','375') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH400Combo.AsymptoticLimits.mH120.root','400') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH425Combo.AsymptoticLimits.mH120.root','425') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH450Combo.AsymptoticLimits.mH120.root','450') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH475Combo.AsymptoticLimits.mH120.root','475') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH500Combo.AsymptoticLimits.mH120.root','500') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH550Combo.AsymptoticLimits.mH120.root','550') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH600Combo.AsymptoticLimits.mH120.root','600') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH650Combo.AsymptoticLimits.mH120.root','650') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH700Combo.AsymptoticLimits.mH120.root','700') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH750Combo.AsymptoticLimits.mH120.root','750') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH800Combo.AsymptoticLimits.mH120.root','800') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH850Combo.AsymptoticLimits.mH120.root','850') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH900Combo.AsymptoticLimits.mH120.root','900') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH950Combo.AsymptoticLimits.mH120.root','950') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH1000Combo.AsymptoticLimits.mH120.root','1000') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH1100Combo.AsymptoticLimits.mH120.root','1100') );
	results3.append( ExtractFile(idir+'/higgsCombineTChiHH1200Combo.AsymptoticLimits.mH120.root','1200') );

	# results_hino=results3
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH200Combo.AsymptoticLimits.mH120.root','200') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH250Combo.AsymptoticLimits.mH120.root','250') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH300Combo.AsymptoticLimits.mH120.root','300') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH350Combo.AsymptoticLimits.mH120.root','350') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH400Combo.AsymptoticLimits.mH120.root','400') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH450Combo.AsymptoticLimits.mH120.root','450') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH500Combo.AsymptoticLimits.mH120.root','500') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH550Combo.AsymptoticLimits.mH120.root','550') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH600Combo.AsymptoticLimits.mH120.root','600') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH650Combo.AsymptoticLimits.mH120.root','650') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH700Combo.AsymptoticLimits.mH120.root','700') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH750Combo.AsymptoticLimits.mH120.root','750') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH800Combo.AsymptoticLimits.mH120.root','800') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH850Combo.AsymptoticLimits.mH120.root','850') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH900Combo.AsymptoticLimits.mH120.root','900') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH950Combo.AsymptoticLimits.mH120.root','950') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH1000Combo.AsymptoticLimits.mH120.root','1000') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH1100Combo.AsymptoticLimits.mH120.root','1100') );
	results_hino.append( ExtractFile(idir+'/higgsCombineTChiHH1200Combo.AsymptoticLimits.mH120.root','1200') );

	#200-1450, 1000-1450 is just every 50
	# xsecs=[1335.62,860.597,577.314,400.107,284.855,207.36,153.841,116.006,88.7325,68.6963,53.7702,42.4699,33.8387,27.1867,21.9868,17.9062,14.6677,12.062,9.96406,8.28246,6.89981,5.78355,4.8731,4.09781,3.46143,2.9337,2.4923,2.13679,1.80616,1.55453,1.32692,1.12975,0.968853,0.538005,0.299347]
	xsecs=[1335.62,860.597,577.314,400.107,284.855,207.36,153.841,116.006,88.7325,68.6963,53.7702,42.4699,33.8387,21.9868,14.6677,9.96406,6.89981,4.8731,3.46143,2.4923,1.80616,1.32692,0.968853,0.538005,0.299347]

	#N2N1, 200-1200 (every 50 until 1000, then 1100 and 1200)
	xsec_hino=[244.213,104.252, 50.9994, 27.3286, 15.6691, 9.44017, 5.90757, 3.8167, 2.53015, 1.71418, 1.18113, 0.826366, 0.586211, 0.420556, 0.305935, 0.22285, 0.16428,0.0912469,0.0516263] #this doesn't include the 25 GeV and 75 GeV mass points

	names=[]; l_obs=[]; l_m2sig=[]; l_m1sig=[]; l_exp=[]; l_p1sig=[]; l_p2sig=[]; count=0;
	names2=[]; l_obs2=[]; l_m2sig2=[]; l_m1sig2=[]; l_exp2=[]; l_p1sig2=[]; l_p2sig2=[]; count2=0;
	names3=[]; l_obs3=[]; l_m2sig3=[]; l_m1sig3=[]; l_exp3=[]; l_p1sig3=[]; l_p2sig3=[]; count3=0;
	# names4=[]; l_obs4=[]; l_m2sig4=[]; l_m1sig4=[]; l_exp4=[]; l_p1sig4=[]; l_p2sig4=[]; count4=0;
	# names5=[]; l_obs5=[]; l_m2sig5=[]; l_m1sig5=[]; l_exp5=[]; l_p1sig5=[]; l_p2sig5=[]; count5=0;

	names_hino=[]; count_hino=0;
	# BR=1.0 #BR set in ALPHABET (might be needed for gluino)
	BR=0.5823329*0.5823329

	for r in results:
		names.append(r[0]);
		l_m2sig.append(r[1]*xsecs[count]*BR);
		l_m1sig.append(r[2]*xsecs[count]*BR);
		l_exp.append(r[3]*xsecs[count]*BR);
		l_p1sig.append(r[4]*xsecs[count]*BR);
		l_p2sig.append(r[5]*xsecs[count]*BR);
		l_obs.append(r[6]*xsecs[count]*BR);
		count=count+1

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

	# for r in results4:
	# 	names4.append(r[0]);
	# 	l_m2sig4.append(r[1]*xsecs[count4]*BR);
	# 	l_m1sig4.append(r[2]*xsecs[count4]*BR);
	# 	l_exp4.append(r[3]*xsecs[count4]*BR);
	# 	l_p1sig4.append(r[4]*xsecs[count4]*BR);
	# 	l_p2sig4.append(r[5]*xsecs[count4]*BR);
	# 	l_obs4.append(r[6]*xsecs[count4]*BR);
	# 	count4=count4+1
	#
	# for r in results5:
	# 	names5.append(r[0]);
	# 	l_m2sig5.append(r[1]*xsecs[count5]*BR);
	# 	l_m1sig5.append(r[2]*xsecs[count5]*BR);
	# 	l_exp5.append(r[3]*xsecs[count5]*BR);
	# 	l_p1sig5.append(r[4]*xsecs[count5]*BR);
	# 	l_p2sig5.append(r[5]*xsecs[count5]*BR);
	# 	l_obs5.append(r[6]*xsecs[count5]*BR);
	# 	count5=count5+1
	#

	for r in results_hino:
		names_hino.append(r[0]);
		count_hino=count_hino+1

	a_xax = array('d', []); a2_xax = array('d', []);
	a_exp = array('d', []); a_obs = array('d', []);
	a_1sig = array('d', []); a_2sig = array('d', []);

	a_xax2 = array('d', []); a_exp2 = array('d', []);
	a_xax3 = array('d', []); a_exp3 = array('d', []);
	# a_xax4 = array('d', []); a_exp4 = array('d', []);
	# a_xax5 = array('d', []); a_exp5 = array('d', []);

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
	# for i in range(len(names4)): a_xax4.append( float(names4[i]) );
	# for i in range(len(l_exp4)): a_exp4.append( float(l_exp4[i]) );
	# for i in range(len(names5)): a_xax5.append( float(names5[i]) );
	# for i in range(len(l_exp5)): a_exp5.append( float(l_exp5[i]) );

	a_2sig.append(results[0][6])
	a2_xax.append(0.5)

	g_exp = ROOT.TGraph(len(a_xax), a_xax, a_exp)
	g_obs = ROOT.TGraph(len(a_xax), a_xax, a_obs)
	g_1sig = ROOT.TGraph(len(2*a_xax), a2_xax, a_1sig)
	g_2sig = ROOT.TGraph(len(2*a_xax), a2_xax, a_2sig)

	g_exp2 = ROOT.TGraph(len(a_xax2), a_xax2, a_exp2)
	g_exp3 = ROOT.TGraph(len(a_xax3), a_xax3, a_exp3)
	# g_exp4 = ROOT.TGraph(len(a_xax4), a_xax4, a_exp4)
	# g_exp5 = ROOT.TGraph(len(a_xax5), a_xax5, a_exp5)

	can = ROOT.TCanvas("can","can",1800,1200);
	hrl = ROOT.TH1F("hrl","hrl",20,100,1200);

	hrl.GetXaxis().SetTitle("Higgsino mass m_{#tilde{#chi}^{0}_{1}} [GeV]");
	hrl.GetXaxis().SetRangeUser(200,1200);
	hrl.GetXaxis().SetLabelSize(0.035)
	hrl.GetXaxis().SetTitleSize(0.04)
	hrl.GetXaxis().SetTitleOffset(1.2)
	hrl.GetYaxis().SetTitle("#sigma [fb] #times BR(HH#rightarrow 4b) ");
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
	hrl.GetYaxis().SetRangeUser(0.1,10000);

	can.SetLeftMargin( 0.11 ); can.SetRightMargin( 0.04 );
	can.SetTopMargin( 0.06 ); can.SetBottomMargin( 0.12 );
	can.SetGrid();
	can.SetLogy();
	txta = ROOT.TLatex(0.2,0.95,"CMS");
	txta.SetNDC();
	txtb = ROOT.TLatex(0.27,0.95,"Preliminary");
	txtb.SetNDC(); txtb.SetTextFont(52); txtb.SetTextSize(0.042);
	txtc = ROOT.TLatex(0.75,0.96,"137 fb^{-1} (13 TeV)");
	txtc.SetNDC(); txtc.SetTextFont(42); txtc.SetTextSize(0.04);
	txtd = ROOT.TLatex(0.64,0.88,"pp#rightarrow #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} #rightarrow HH #tilde{G} #tilde{G}");
	txtd.SetNDC(); txtd.SetTextFont(42); txtd.SetTextSize(0.05);
	txte = ROOT.TLatex(0.78,0.82,"m_{#tilde{G}} = 1 GeV");
	txte.SetNDC(); txte.SetTextFont(42); txte.SetTextSize(0.05);
	#f=open("LatestXsecGluGlu.txt", 'r')
	a_stop = array('d', []);
	a_xsec = array('d', []);
	count=0;
	for x in xsecs:
		a_stop.append(float(names[count]));
		a_xsec.append(x*0.5823329*0.5823329) #I believe I need to multiply the theory xsec, maybe
		count=count+1
	g_xsec=ROOT.TGraph(len(a_stop), a_stop, a_xsec)

	a_stop_hino = array('d', []);
	a_xsec_hino = array('d', []);
	count_hino=0;
	for x in xsec_hino:
		a_stop_hino.append(float(names_hino[count_hino]));
		a_xsec_hino.append(x*0.5823329*0.5823329)
		count_hino=count_hino+1
	g_xsec_hino=ROOT.TGraph(len(a_stop_hino), a_stop_hino, a_xsec_hino)

	leg = ROOT.TLegend(0.48,0.48,0.93,0.78);
	leg.SetFillStyle(1001);
	leg.SetFillColor(0);
	leg.SetBorderSize(1);
	# leg.SetNColumns(2);

	# leg.AddEntry(g_exp,"Expected","l") #E commented out
	leg.AddEntry(g_exp, "Boosted Only","l")
	leg.AddEntry(g_exp2,"Resolved Only","l")
	# leg.AddEntry(g_exp3,"Boosted+Resolved","l")
	# leg.AddEntry(g_exp4,"Combo, boosted MET>300","l")
	# leg.AddEntry(g_exp5,"Full Combo","l")

	# leg.AddEntry(g_obs,"Observed","l")#E commented out
	# leg.AddEntry(g_2sig,"Expected 2#sigma","f") #E commented out
	# leg.AddEntry(g_1sig,"Expected 1#sigma","f") #E commented out
  	leg.AddEntry(g_xsec, "Theory cross-section", "l")
  	leg.AddEntry(g_xsec_hino, "Theory cross-section, N2N1 only", "l")

	g_1sig.SetFillColor(ROOT.kGreen); g_1sig.SetFillStyle(3244);
	g_2sig.SetFillColor(ROOT.kYellow); g_2sig.SetFillStyle(3244);

	# g_exp.SetLineStyle(2);
	# g_exp.SetLineColor(ROOT.kBlack);
	# g_exp.SetLineWidth(3);
	# g_obs.SetLineWidth(2);

	g_exp.SetLineColor(ROOT.kRed); g_exp.SetLineWidth(2);
	g_exp2.SetLineColor(ROOT.kBlue); g_exp2.SetLineWidth(2);
	g_exp3.SetLineColor(ROOT.kBlack); g_exp3.SetLineWidth(2); #g_exp3.SetLineStyle(2);
	# g_exp4.SetLineColor(ROOT.kRed); g_exp4.SetLineWidth(2);  g_exp4.SetLineStyle(2);
	# g_exp5.SetLineColor(ROOT.kBlack); g_exp5.SetLineWidth(2);

	# For one line with 1- and 2-sigma bands
	# g_2sig.Draw('f');
	# g_1sig.Draw('fsames');
	# g_obs.Draw('lsames');
	# g_exp.Draw('lsames');

	# For multiple lines
	g_exp.Draw('l');
	g_exp2.Draw('lsames');
	# g_exp3.Draw('lsames');
	# g_exp4.Draw('lsames');
	# g_exp5.Draw('lsames');

	# for i in range(0,100):
	# 	if (g_xsec.Eval(900+i)-g_exp2.Eval(900+i)>0):
	# 		print "Mass %d  Exp Excl %g " %(900+i,g_exp2.Eval(900+i))
	# 		print "Theory Xsec %g " %g_xsec.Eval(900+i)
	#oneLine.Draw("LSAMES");

	txta.Draw();
	txtb.Draw();
	txtc.Draw();
	txtd.Draw();
	txte.Draw();
	leg.Draw();
	g_xsec.SetLineStyle(2);
	g_xsec.SetLineWidth(3);
	# g_xsec.SetLineColor(ROOT.kRed);
	g_xsec.SetLineColor(ROOT.kBlack);
	g_xsec.Draw("lsame")

	#Second xsec (https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVn1n2hino)
	g_xsec_hino.SetLineStyle(3);
	g_xsec_hino.SetLineWidth(2);
	g_xsec_hino.SetLineColor(ROOT.kMagenta);
	g_xsec_hino.Draw("lsame")

	can.SaveAs('TChiHHResults_Combo.pdf');
