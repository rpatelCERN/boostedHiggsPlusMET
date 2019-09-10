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
# ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
# ROOT.setTDRStyle();
# ROOT.gStyle.SetPadLeftMargin(0.16);
# ROOT.gStyle.SetPadRightMargin(0.10);
# ROOT.gStyle.SetPadTopMargin(0.10);
# ROOT.gStyle.SetPalette(1);

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

	#idir = "/eos/uscms/store/user/ntran/SUSY/statInterp/scanOutput/Dec6";
	#idir = "/uscms_data/d2/rgp230/BoostedHPush/NewCommit/CMSSW_7_4_2/src/boostedHiggsPlusMET/datacardsRateParamTest/";
	idir = "./";
	results = []; results2 = []; results3 = []; results4 = []; results5 = [];

	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH200.AsymptoticLimits.mH200.root','200') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH225.AsymptoticLimits.mH225.root','225') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH250.AsymptoticLimits.mH250.root','250') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH275.AsymptoticLimits.mH275.root','275') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH300.AsymptoticLimits.mH300.root','300') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH325.AsymptoticLimits.mH325.root','325') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH350.AsymptoticLimits.mH350.root','350') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH375.AsymptoticLimits.mH375.root','375') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH400.AsymptoticLimits.mH400.root','400') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH425.AsymptoticLimits.mH425.root','425') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH450.AsymptoticLimits.mH450.root','450') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH475.AsymptoticLimits.mH475.root','475') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH500.AsymptoticLimits.mH500.root','500') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH525.AsymptoticLimits.mH525.root','525') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH550.AsymptoticLimits.mH550.root','550') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH575.AsymptoticLimits.mH575.root','575') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH600.AsymptoticLimits.mH600.root','600') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH625.AsymptoticLimits.mH625.root','625') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH650.AsymptoticLimits.mH650.root','650') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH675.AsymptoticLimits.mH675.root','675') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH700.AsymptoticLimits.mH700.root','700') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH725.AsymptoticLimits.mH725.root','725') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH750.AsymptoticLimits.mH750.root','750') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH775.AsymptoticLimits.mH775.root','775') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH800.AsymptoticLimits.mH800.root','800') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH825.AsymptoticLimits.mH825.root','825') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH850.AsymptoticLimits.mH850.root','850') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH875.AsymptoticLimits.mH875.root','875') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH900.AsymptoticLimits.mH900.root','900') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH925.AsymptoticLimits.mH925.root','925') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH950.AsymptoticLimits.mH950.root','950') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH975.AsymptoticLimits.mH975.root','975') );
	results.append( ExtractFile(idir+'/higgsCombineHiggsinoResolvedOnlyTChiHH1000.AsymptoticLimits.mH1000.root','1000') );

	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH200.AsymptoticLimits.mH200.root','200') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH225.AsymptoticLimits.mH225.root','225') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH250.AsymptoticLimits.mH250.root','250') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH275.AsymptoticLimits.mH275.root','275') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH300.AsymptoticLimits.mH300.root','300') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH325.AsymptoticLimits.mH325.root','325') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH350.AsymptoticLimits.mH350.root','350') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH375.AsymptoticLimits.mH375.root','375') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH400.AsymptoticLimits.mH400.root','400') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH425.AsymptoticLimits.mH425.root','425') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH450.AsymptoticLimits.mH450.root','450') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH475.AsymptoticLimits.mH475.root','475') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH500.AsymptoticLimits.mH500.root','500') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH525.AsymptoticLimits.mH525.root','525') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH550.AsymptoticLimits.mH550.root','550') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH575.AsymptoticLimits.mH575.root','575') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH600.AsymptoticLimits.mH600.root','600') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH625.AsymptoticLimits.mH625.root','625') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH650.AsymptoticLimits.mH650.root','650') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH675.AsymptoticLimits.mH675.root','675') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH700.AsymptoticLimits.mH700.root','700') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH725.AsymptoticLimits.mH725.root','725') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH750.AsymptoticLimits.mH750.root','750') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH775.AsymptoticLimits.mH775.root','775') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH800.AsymptoticLimits.mH800.root','800') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH825.AsymptoticLimits.mH825.root','825') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH850.AsymptoticLimits.mH850.root','850') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH875.AsymptoticLimits.mH875.root','875') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH900.AsymptoticLimits.mH900.root','900') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH925.AsymptoticLimits.mH925.root','925') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH950.AsymptoticLimits.mH950.root','950') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH975.AsymptoticLimits.mH975.root','975') );
	results2.append( ExtractFile(idir+'/higgsCombineHiggsinoBoostedOnlyTChiHH1000.AsymptoticLimits.mH1000.root','1000') );

	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH200.AsymptoticLimits.mH200.root','200') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH225.AsymptoticLimits.mH225.root','225') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH250.AsymptoticLimits.mH250.root','250') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH275.AsymptoticLimits.mH275.root','275') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH300.AsymptoticLimits.mH300.root','300') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH325.AsymptoticLimits.mH325.root','325') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH350.AsymptoticLimits.mH350.root','350') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH375.AsymptoticLimits.mH375.root','375') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH400.AsymptoticLimits.mH400.root','400') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH425.AsymptoticLimits.mH425.root','425') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH450.AsymptoticLimits.mH450.root','450') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH475.AsymptoticLimits.mH475.root','475') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH500.AsymptoticLimits.mH500.root','500') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH525.AsymptoticLimits.mH525.root','525') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH550.AsymptoticLimits.mH550.root','550') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH575.AsymptoticLimits.mH575.root','575') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH600.AsymptoticLimits.mH600.root','600') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH625.AsymptoticLimits.mH625.root','625') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH650.AsymptoticLimits.mH650.root','650') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH675.AsymptoticLimits.mH675.root','675') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH700.AsymptoticLimits.mH700.root','700') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH725.AsymptoticLimits.mH725.root','725') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH750.AsymptoticLimits.mH750.root','750') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH775.AsymptoticLimits.mH775.root','775') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH800.AsymptoticLimits.mH800.root','800') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH825.AsymptoticLimits.mH825.root','825') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH850.AsymptoticLimits.mH850.root','850') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH875.AsymptoticLimits.mH875.root','875') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH900.AsymptoticLimits.mH900.root','900') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH925.AsymptoticLimits.mH925.root','925') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH950.AsymptoticLimits.mH950.root','950') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH975.AsymptoticLimits.mH975.root','975') );
	results3.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostOnlyTChiHH1000.AsymptoticLimits.mH1000.root','1000') );


	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH200.AsymptoticLimits.mH200.root','200') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH225.AsymptoticLimits.mH225.root','225') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH250.AsymptoticLimits.mH250.root','250') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH275.AsymptoticLimits.mH275.root','275') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH300.AsymptoticLimits.mH300.root','300') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH325.AsymptoticLimits.mH325.root','325') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH350.AsymptoticLimits.mH350.root','350') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH375.AsymptoticLimits.mH375.root','375') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH400.AsymptoticLimits.mH400.root','400') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH425.AsymptoticLimits.mH425.root','425') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH450.AsymptoticLimits.mH450.root','450') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH475.AsymptoticLimits.mH475.root','475') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH500.AsymptoticLimits.mH500.root','500') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH525.AsymptoticLimits.mH525.root','525') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH550.AsymptoticLimits.mH550.root','550') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH575.AsymptoticLimits.mH575.root','575') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH600.AsymptoticLimits.mH600.root','600') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH625.AsymptoticLimits.mH625.root','625') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH650.AsymptoticLimits.mH650.root','650') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH675.AsymptoticLimits.mH675.root','675') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH700.AsymptoticLimits.mH700.root','700') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH725.AsymptoticLimits.mH725.root','725') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH750.AsymptoticLimits.mH750.root','750') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH775.AsymptoticLimits.mH775.root','775') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH800.AsymptoticLimits.mH800.root','800') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH825.AsymptoticLimits.mH825.root','825') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH850.AsymptoticLimits.mH850.root','850') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH875.AsymptoticLimits.mH875.root','875') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH900.AsymptoticLimits.mH900.root','900') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH925.AsymptoticLimits.mH925.root','925') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH950.AsymptoticLimits.mH950.root','950') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH975.AsymptoticLimits.mH975.root','975') );
	results4.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullComboTChiHH1000.AsymptoticLimits.mH1000.root','1000') );


	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH200.AsymptoticLimits.mH200.root','200') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH225.AsymptoticLimits.mH225.root','225') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH250.AsymptoticLimits.mH250.root','250') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH275.AsymptoticLimits.mH275.root','275') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH300.AsymptoticLimits.mH300.root','300') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH325.AsymptoticLimits.mH325.root','325') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH350.AsymptoticLimits.mH350.root','350') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH375.AsymptoticLimits.mH375.root','375') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH400.AsymptoticLimits.mH400.root','400') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH425.AsymptoticLimits.mH425.root','425') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH450.AsymptoticLimits.mH450.root','450') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH475.AsymptoticLimits.mH475.root','475') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH500.AsymptoticLimits.mH500.root','500') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH525.AsymptoticLimits.mH525.root','525') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH550.AsymptoticLimits.mH550.root','550') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH575.AsymptoticLimits.mH575.root','575') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH600.AsymptoticLimits.mH600.root','600') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH625.AsymptoticLimits.mH625.root','625') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH650.AsymptoticLimits.mH650.root','650') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH675.AsymptoticLimits.mH675.root','675') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH700.AsymptoticLimits.mH700.root','700') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH725.AsymptoticLimits.mH725.root','725') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH750.AsymptoticLimits.mH750.root','750') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH775.AsymptoticLimits.mH775.root','775') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH800.AsymptoticLimits.mH800.root','800') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH825.AsymptoticLimits.mH825.root','825') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH850.AsymptoticLimits.mH850.root','850') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH875.AsymptoticLimits.mH875.root','875') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH900.AsymptoticLimits.mH900.root','900') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH925.AsymptoticLimits.mH925.root','925') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH950.AsymptoticLimits.mH950.root','950') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH975.AsymptoticLimits.mH975.root','975') );
	results5.append( ExtractFile(idir+'/higgsCombineHiggsinoResToBoostFullTest1TChiHH1000.AsymptoticLimits.mH1000.root','1000') );

	xsecs=[1335.62,860.597,577.314,400.107,284.855,207.36,153.841,116.006,88.7325,68.6963,53.7702,42.4699,33.8387,27.1867,21.9868,17.9062,14.6677,12.062,9.96406,8.28246,6.89981,5.78355,4.8731,4.09781,3.46143,2.9337,2.4923,2.13679,1.80616,1.55453,1.32692,1.12975,0.968853]

	names   = []; names2   = []; names3   = []; names4   = []; names5   = [];
	l_obs   = []; l_obs2   = []; l_obs3   = []; l_obs4   = []; l_obs5   = [];
	l_m2sig = []; l_m2sig2 = []; l_m2sig3 = []; l_m2sig4 = []; l_m2sig5 = [];
	l_m1sig = []; l_m1sig2 = []; l_m1sig3 = []; l_m1sig4 = []; l_m1sig5 = [];
	l_exp   = []; l_exp2   = []; l_exp3   = []; l_exp4   = []; l_exp5   = [];
	l_p1sig = []; l_p1sig2 = []; l_p1sig3 = []; l_p1sig4 = []; l_p1sig5 = [];
	l_p2sig = []; l_p2sig2 = []; l_p2sig3 = []; l_p2sig4 = []; l_p2sig5 = [];
	count=0; count2=0; count3=0; count4=0; count5=0;
	BR=0.57*0.57 #BR not needed for gluino - fully simulated
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

	for r in results5:
		names5.append(r[0]);
		l_m2sig5.append(r[1]*xsecs[count5]*BR);
		l_m1sig5.append(r[2]*xsecs[count5]*BR);
		l_exp5.append(r[3]*xsecs[count5]*BR);
		l_p1sig5.append(r[4]*xsecs[count5]*BR);
		l_p2sig5.append(r[5]*xsecs[count5]*BR);
		l_obs5.append(r[6]*xsecs[count5]*BR);
		count5=count5+1

	a_xax = array('d', []); a_xax2 = array('d', []); a_xax3 = array('d', []); a_xax4 = array('d', []); a_xax5 = array('d', []);
	a2_xax = array('d', []);
	a_exp = array('d', []);a_exp2 = array('d', []); a_exp3 = array('d', []); a_exp4 = array('d', []);  a_exp5 = array('d', []);
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

	for i in range(len(names5)): a_xax5.append( float(names5[i]) );
	for i in range(len(l_exp5)): a_exp5.append( float(l_exp5[i]) );

	a_2sig.append(results[0][6])
	a2_xax.append(0.5)

	g_exp = ROOT.TGraph(len(a_xax), a_xax, a_exp)
	g_obs = ROOT.TGraph(len(a_xax), a_xax, a_obs)
	g_1sig = ROOT.TGraph(len(2*a_xax), a2_xax, a_1sig)
	g_2sig = ROOT.TGraph(len(2*a_xax), a2_xax, a_2sig)

	g_exp2 = ROOT.TGraph(len(a_xax2), a_xax2, a_exp2)
	g_exp3 = ROOT.TGraph(len(a_xax3), a_xax3, a_exp3)
	g_exp4 = ROOT.TGraph(len(a_xax4), a_xax4, a_exp4)
	g_exp5 = ROOT.TGraph(len(a_xax5), a_xax5, a_exp5)

	can = ROOT.TCanvas("can","can",1800,1200);
	hrl = ROOT.TH1F("hrl","hrl",36,100,1000);

	# hrl = can.DrawFrame(0,0,6,15);
	hrl.GetXaxis().SetTitle("Higgsino mass m_{#tilde{#chi}^{0}_{1}} [GeV]");
	hrl.GetXaxis().SetRangeUser(200,1000);
	hrl.GetXaxis().SetLabelSize(0.035)
	hrl.GetXaxis().SetTitleSize(0.04)
	hrl.GetXaxis().SetTitleOffset(1.2)
	hrl.GetYaxis().SetTitle("#sigma_{95% CL} x BR(hh#rightarrow bbbb) [fb] ");
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
	hrl.GetYaxis().SetRangeUser(0.1,10000);

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
	txtd = ROOT.TLatex(0.65,0.88,"pp#rightarrow #tilde{#chi}^{1}_{0} #tilde{#chi}^{1}_{0} #rightarrow hh #tilde{G} #tilde{G}");
	txtd.SetNDC(); txtd.SetTextFont(42); txtd.SetTextSize(0.05);
	#f=open("LatestXsecGluGlu.txt", 'r')
	a_stop = array('d', []);
	a_xsec = array('d', []);
	count=0;
	for x in xsecs:
		a_stop.append(float(names[count]));
		a_xsec.append(x*0.57*0.57)
		count=count+1
	g_xsec=ROOT.TGraph(len(a_stop), a_stop, a_xsec)
	leg = ROOT.TLegend(0.65,0.60,0.9,0.8);
	leg.SetFillStyle(1001);
	leg.SetFillColor(0);
	leg.SetBorderSize(1);
	# leg.SetNColumns(2);

	# leg.AddEntry(g_exp,"expected","l") #E commented out
	# leg.AddEntry(g_exp, "Resolved Only","l")
	# leg.AddEntry(g_exp2,"Boosted Only","l")
	# leg.AddEntry(g_exp3,"Resolved+Boost","l")
	leg.AddEntry(g_exp4,"Full Combo","l")
	leg.AddEntry(g_exp5,"Full Combo, test1","l")

	# leg.AddEntry(g_obs,"observed","l")#E commented out
	# leg.AddEntry(g_2sig,"expected 2#sigma","f") #E commented out
	# leg.AddEntry(g_1sig,"expected 1#sigma","f") #E commented out
  	leg.AddEntry(g_xsec, "Theory cross-section", "l")
	#oneLine = ROOT.TF1("oneLine","1",175,550);
	#oneLine.SetLineColor(ROOT.kRed+2);
	#oneLine.SetLineWidth(2);
	#oneLine.SetLineStyle(1);

	# g_1sig.SetFillColor(ROOT.kGreen);
	# g_1sig.SetFillStyle(3244);
	# g_2sig.SetFillColor(ROOT.kYellow);
	# g_2sig.SetFillStyle(3244);
	# g_exp.SetLineStyle(2);
	# g_exp.SetLineColor(ROOT.kBlack);
	# g_exp.SetLineWidth(2);
	# g_obs.SetLineWidth(2);

	g_exp.SetLineColor(ROOT.kCyan+2); g_exp.SetLineWidth(2);
	g_exp2.SetLineColor(ROOT.kRed); g_exp2.SetLineWidth(2);
	g_exp3.SetLineColor(ROOT.kBlue); g_exp3.SetLineWidth(2);
	g_exp4.SetLineColor(ROOT.kBlack); g_exp4.SetLineWidth(2);
	g_exp5.SetLineColor(ROOT.kGreen+2); g_exp5.SetLineWidth(2);


	# g_2sig.Draw('f'); #E commented out
	# g_1sig.Draw('fsames'); #E commented out
	#g_1sig.Draw('f');
	#g_obs.Draw('lsames');
	# g_exp.Draw('lsames'); #E commented out



	# g_exp.Draw('l');
	# g_exp2.Draw('lsames');
	# g_exp3.Draw('lsames');
	# g_exp4.Draw('lsames');

	g_exp4.Draw('l');
	g_exp5.Draw('lsames');
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
	# can.SaveAs('TChiHHResults_FullCombo_RtoB.pdf');
	can.SaveAs('TChiHHResults_FullCombo_Test1.pdf');
	# can.SaveAs('TChiHHResults_ResolvedAndBoosted.pdf');
