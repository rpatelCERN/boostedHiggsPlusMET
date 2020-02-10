import os
from ROOT import *
hino=[]

for h in range(0, 35):
	hino.append(150+h*25)
for h in hino:
	fin=TFile("Skims/tree_TChiHH_HToBB_HToBB_%d_1_MC2016_fast.root" %h,"READ")
	TotalEvents=fin.Get("nEventProc")
	os.system("python QuickDataCardsABCD.py %d %d " %(h, TotalEvents.GetBinContent(1)))
