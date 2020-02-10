from ROOT import *
import os
Bkgs=["QCD","ZJets", "WJets", "TT"]

#Files
# directory = "lowMETbin"
# directory = "BoostedOnly_medB/"
# directory = "BoostToRes/"
directory = "ResToBoost/"
# directory = "ResToBoost_Test3/"
# directory = "ResolvedOnly_medB/"
# directory = "BoostedResolved_medB/"


fboost=TFile(directory+"Boost/Bkg.root", "READ");
fDiAk8=TFile(directory+"DiAK8/Bkg.root", "READ");
fSemiRes=TFile(directory+"SemiResolved/Bkg.root", "READ");
fSingleAk8=TFile(directory+"SingleAK8/Bkg.root", "READ");
fRes4b_lowR=TFile(directory+"/Res4b_lowR/Bkg.root", "READ");
fRes4b_highR=TFile(directory+"/Res4b_highR/Bkg.root", "READ");
fRes3b_lowR=TFile(directory+"/Res3b_lowR/Bkg.root", "READ");
fRes3b_highR=TFile(directory+"/Res3b_highR/Bkg.root", "READ");

BoostedBins=fboost.Get("%s_MET" %Bkgs[0])
DiAk8=fDiAk8.Get("%s_MET" %Bkgs[0])
SemiRes=fSemiRes.Get("%s_MET" %Bkgs[0])
SingleAK8=fSingleAk8.Get("%s_MET" %Bkgs[0])
Res3b_lowR=fRes3b_lowR.Get("%s_MET" %Bkgs[0]) #initialize to get the number of MET bins
Res3b_highR=fRes3b_highR.Get("%s_MET" %Bkgs[0]) #initialize to get the number of MET bins
Res4b_lowR=fRes4b_lowR.Get("%s_MET" %Bkgs[0]) #initialize to get the number of MET bins
Res4b_highR=fRes4b_highR.Get("%s_MET" %Bkgs[0]) #initialize to get the number of MET bins


FullBkgYieldsBoostH=[[]]
FullBkgYieldsDiAk8=[[]]
FullBkgYieldsSemi=[[]]
FullBkgYieldsSingleAk8=[[]]
FullBkgYieldsRes3b_lowR=[[]]
FullBkgYieldsRes3b_highR=[[]]
FullBkgYieldsRes4b_lowR=[[]]
FullBkgYieldsRes4b_highR=[[]]


for b in Bkgs:

	METYields=[]
	BoostedBins=fboost.Get("%s_MET" %b)
	BoostedBins.SetDirectory(0)
	for m in range(1, BoostedBins.GetNbinsX()+1):
		METYields.append(BoostedBins.GetBinContent(m))
	FullBkgYieldsBoostH.append(METYields)

	METYields=[]
	DiAk8=fDiAk8.Get("%s_MET" %b)
	DiAk8.SetDirectory(0)
	for m in range(1, DiAk8.GetNbinsX()+1):
		METYields.append(DiAk8.GetBinContent(m))
	FullBkgYieldsDiAk8.append(METYields)

	METYields=[]
	SemiRes=fSemiRes.Get("%s_MET" %b)
	SemiRes.SetDirectory(0)
	for m in range(1, SemiRes.GetNbinsX()+1):
		METYields.append(SemiRes.GetBinContent(m))
	FullBkgYieldsSemi.append(METYields)

	METYields=[]
	SingleAK8=fSingleAk8.Get("%s_MET" %b)
	SemiRes.SetDirectory(0)
	for m in range(1, SingleAK8.GetNbinsX()+1):
		METYields.append(SingleAK8.GetBinContent(m))
	FullBkgYieldsSingleAk8.append(METYields)

	Res4b_lowR=fRes4b_lowR.Get("%s_MET" %b)
	Res4b_lowR.SetDirectory(0)
	METYields=[]
	for m in range(1,Res4b_lowR.GetNbinsX()+1):
		METYields.append(Res4b_lowR.GetBinContent(m));
	FullBkgYieldsRes4b_lowR.append(METYields)

	Res4b_highR=fRes4b_highR.Get("%s_MET" %b)
	Res4b_highR.SetDirectory(0)
	METYields=[]
	for m in range(1,Res4b_highR.GetNbinsX()+1):
		METYields.append(Res4b_highR.GetBinContent(m));
	FullBkgYieldsRes4b_highR.append(METYields)

	Res3b_lowR=fRes3b_lowR.Get("%s_MET" %b)
	Res3b_lowR.SetDirectory(0)
	METYields=[]
	for m in range(1,Res3b_lowR.GetNbinsX()+1):
		METYields.append(Res3b_lowR.GetBinContent(m));
	FullBkgYieldsRes3b_lowR.append(METYields)

	Res3b_highR=fRes3b_highR.Get("%s_MET" %b)
	Res3b_highR.SetDirectory(0)
	METYields=[]
	for m in range(1,Res3b_highR.GetNbinsX()+1):
		METYields.append(Res3b_highR.GetBinContent(m));
	FullBkgYieldsRes3b_highR.append(METYields)


print FullBkgYieldsRes4b_lowR,FullBkgYieldsRes3b_highR

fboostSig=TFile(directory+"Boost/TChiHH.root", "READ");
fDiAk8Sig=TFile(directory+"DiAK8/TChiHH.root", "READ");
fSemiResSig=TFile(directory+"SemiResolved/TChiHH.root", "READ");
fSingleAk8Sig=TFile(directory+"SingleAK8/TChiHH.root", "READ");

fResSig4b_lowR=TFile(directory+"Res4b_lowR/TChiHH.root", "READ");
fResSig4b_highR=TFile(directory+"Res4b_highR/TChiHH.root", "READ");
fResSig3b_lowR=TFile(directory+"Res3b_lowR/TChiHH.root", "READ");
fResSig3b_highR=TFile(directory+"Res3b_highR/TChiHH.root", "READ");
mSig=[127, 150, 175,200,225, 250, 275, 300, 325,350, 375, 400, 425,450,475,500, 525, 550,575, 600, 625, 650, 675, 700,725, 750,775, 800,825, 850, 875, 900, 925, 950, 975, 1000]


for s in mSig:

	SigHistoBoost=fboostSig.Get("TChiHH%d_MET" %s)
	SigHisto2Ak8=fDiAk8Sig.Get("TChiHH%d_MET" %s)
	SigHistoSemiRes=fSemiResSig.Get("TChiHH%d_MET" %s)
	SigHisto1Ak8=fSingleAk8Sig.Get("TChiHH%d_MET" %s)
	SigHisto3b_lowR=fResSig3b_lowR.Get("TChiHH%d_MET" %s)
	SigHisto3b_highR=fResSig3b_highR.Get("TChiHH%d_MET" %s)
	SigHisto4b_lowR=fResSig4b_lowR.Get("TChiHH%d_MET" %s)
	SigHisto4b_highR=fResSig4b_highR.Get("TChiHH%d_MET" %s)

	CardsBoost=""
	DiJetCardsBoost=""
	SemiResolvedCards=""
	SingleJetCardsBoost=""
	ResolvedCards3b_lowR=""
	ResolvedCards3b_highR=""
	ResolvedCards4b_lowR=""
	ResolvedCards4b_highR=""


	for m in range(0,SigHistoBoost.GetNbinsX()):
		fout=open("HiggsinoBoosted_MET%d.txt" %m,'w')
		fout.seek(0)
		fout.write("# Simple counting experiment, with higgsino signal and a few background processes\n")
		fout.write("imax 1  number of channels\n") #one .txt card for each MET bin
		fout.write("jmax %d  number of backgrounds\n" %len(Bkgs)) #total bkg processes
		fout.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
		fout.write("bin bin%d\n" %m)
		Obs=0;
		for i in range(1, len(Bkgs)+1):Obs=Obs+FullBkgYieldsBoostH[i][m]
		fout.write("observation %g\n" %Obs)
		fout.write("bin ");
		for i in range(0, len(Bkgs)+1):fout.write("bin%d " %m)
		fout.write("\nprocess sig");
		for b in Bkgs: fout.write(" %s" %b)
		fout.write("\nprocess")
		for i in range(0, len(Bkgs)+1):fout.write(" %d" %i)
		Sig=SigHistoBoost.GetBinContent(m+1)
	 	if(m==SigHistoBoost.GetNbinsX()-1):Sig=Sig+SigHistoBoost.GetBinContent(m+2)#ADD OVERFLOW Bin
		fout.write("\nrate %g " %(Sig));
		for i in range(1, len(Bkgs)+1):
			if FullBkgYieldsBoostH[i][m]>0:  fout.write(" %g " %FullBkgYieldsBoostH[i][m])
			else:fout.write(" 0.00001 ")
		fout.close()
		CardsBoost=CardsBoost+" HiggsinoBoosted_MET%d.txt" %m

	for m in range(0,SigHisto2Ak8.GetNbinsX()):
		fout=open("HiggsinoDiAk8_MET%d.txt" %m,'w')
		fout.seek(0)
		fout.write("# Simple counting experiment, with higgsino signal and a few background processes\n")
		fout.write("imax 1  number of channels\n") #one .txt card for each MET bin
		fout.write("jmax %d  number of backgrounds\n" %len(Bkgs)) #total bkg processes
		fout.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
		fout.write("bin bin%d\n" %m)
		Obs=0;
		for i in range(1, len(Bkgs)+1):Obs=Obs+FullBkgYieldsDiAk8[i][m]
		fout.write("observation %g\n" %Obs)
		fout.write("bin ");
		for i in range(0, len(Bkgs)+1):fout.write("bin%d " %m)
		fout.write("\nprocess sig");
		for b in Bkgs: fout.write(" %s" %b)
		fout.write("\nprocess")
		for i in range(0, len(Bkgs)+1):fout.write(" %d" %i)
		Sig=SigHisto2Ak8.GetBinContent(m+1)
	 	if(m==SigHisto2Ak8.GetNbinsX()-1):Sig=Sig+SigHisto2Ak8.GetBinContent(m+2)#ADD OVERFLOW Bin
		fout.write("\nrate %g " %(Sig));
		for i in range(1, len(Bkgs)+1):
			if FullBkgYieldsDiAk8[i][m]>0:  fout.write(" %g " %FullBkgYieldsDiAk8[i][m])
			else:fout.write(" 0.00001 ")
		fout.close()
		DiJetCardsBoost=DiJetCardsBoost+" HiggsinoDiAk8_MET%d.txt" %m

	for m in range(0,SigHistoSemiRes.GetNbinsX()):
		fout=open("HiggsinoSemiResolved_MET%d.txt" %m,'w')
		fout.seek(0)
		fout.write("# Simple counting experiment, with higgsino signal and a few background processes\n")
		fout.write("imax 1  number of channels\n") #one .txt card for each MET bin
		fout.write("jmax %d  number of backgrounds\n" %len(Bkgs)) #total bkg processes
		fout.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
		fout.write("bin bin%d\n" %m)
		Obs=0;
		for i in range(1, len(Bkgs)+1):Obs=Obs+FullBkgYieldsSemi[i][m]
		fout.write("observation %g\n" %Obs)
		fout.write("bin ");
		for i in range(0, len(Bkgs)+1):fout.write("bin%d " %m)
		fout.write("\nprocess sig");
		for b in Bkgs: fout.write(" %s" %b)
		fout.write("\nprocess")
		for i in range(0, len(Bkgs)+1):fout.write(" %d" %i)
		Sig=SigHistoSemiRes.GetBinContent(m+1)
	 	if(m==SigHistoSemiRes.GetNbinsX()-1):Sig=Sig+SigHistoSemiRes.GetBinContent(m+2)#ADD OVERFLOW Bin
		fout.write("\nrate %g " %(Sig));
		for i in range(1, len(Bkgs)+1):
			if FullBkgYieldsSemi[i][m]>0:  fout.write(" %g " %FullBkgYieldsSemi[i][m])
			else:fout.write(" 0.00001 ")
		fout.close()
		SemiResolvedCards=SemiResolvedCards+" HiggsinoSemiResolved_MET%d.txt" %m

	for m in range(0,SigHisto1Ak8.GetNbinsX()):
		fout=open("HiggsinoSingleAk8_MET%d.txt" %m,'w')
		fout.seek(0)
		fout.write("# Simple counting experiment, with higgsino signal and a few background processes\n")
		fout.write("imax 1  number of channels\n") #one .txt card for each MET bin
		fout.write("jmax %d  number of backgrounds\n" %len(Bkgs)) #total bkg processes
		fout.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
		fout.write("bin bin%d\n" %m)
		Obs=0;
		for i in range(1, len(Bkgs)+1):Obs=Obs+FullBkgYieldsSingleAk8[i][m]
		fout.write("observation %g\n" %Obs)
		fout.write("bin ");
		for i in range(0, len(Bkgs)+1):fout.write("bin%d " %m)
		fout.write("\nprocess sig");
		for b in Bkgs: fout.write(" %s" %b)
		fout.write("\nprocess")
		for i in range(0, len(Bkgs)+1):fout.write(" %d" %i)
		Sig=SigHisto1Ak8.GetBinContent(m+1)
	 	if(m==SigHisto1Ak8.GetNbinsX()-1):Sig=Sig+SigHisto1Ak8.GetBinContent(m+2)#ADD OVERFLOW Bin
		fout.write("\nrate %g " %(Sig));
		for i in range(1, len(Bkgs)+1):
			if FullBkgYieldsSingleAk8[i][m]>0:  fout.write(" %g " %FullBkgYieldsSingleAk8[i][m])
			else:fout.write(" 0.00001 ")
		fout.close()
		SingleJetCardsBoost=SingleJetCardsBoost+" HiggsinoSingleAk8_MET%d.txt" %m


	for m in range(0,SigHisto3b_lowR.GetNbinsX()):
		fout=open("HiggsinoRes3b_lowR_MET%d.txt" %m,'w')
		fout.seek(0)
		fout.write("# Simple counting experiment, with higgsino signal and a few background processes\n")
		fout.write("imax 1  number of channels\n") #one .txt card for each MET bin
		fout.write("jmax %d  number of backgrounds\n" %len(Bkgs)) #total bkg processes
		fout.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
		fout.write("bin bin%d\n" %m)
		Obs=0;
		for i in range(1, len(Bkgs)+1):Obs=Obs+FullBkgYieldsRes3b_lowR[i][m]
		fout.write("observation %g\n" %Obs)
		fout.write("bin ");
		for i in range(0, len(Bkgs)+1):fout.write("bin%d " %m)
		fout.write("\nprocess sig");
		for b in Bkgs: fout.write(" %s" %b)
		fout.write("\nprocess")
		for i in range(0, len(Bkgs)+1):fout.write(" %d" %i)
		Sig=SigHisto3b_lowR.GetBinContent(m+1)
	 	if(m==SigHisto3b_lowR.GetNbinsX()-1):Sig=Sig+SigHisto3b_lowR.GetBinContent(m+2)#ADD OVERFLOW Bin
		fout.write("\nrate %g " %(Sig));
		for i in range(1, len(Bkgs)+1):
			#fout.write(" %g " %FullBkgYieldsRes3b[i][m])
			if FullBkgYieldsRes3b_lowR[i][m]>0:  fout.write(" %g " %FullBkgYieldsRes3b_lowR[i][m])
			else:fout.write(" 0.00001 ")
		fout.close()
		ResolvedCards3b_lowR=ResolvedCards3b_lowR+" HiggsinoRes3b_lowR_MET%d.txt" %m

	for m in range(0,SigHisto3b_highR.GetNbinsX()):
		fout=open("HiggsinoRes3b_highR_MET%d.txt" %m,'w')
		fout.seek(0)
		fout.write("# Simple counting experiment, with higgsino signal and a few background processes\n")
		fout.write("imax 1  number of channels\n") #one .txt card for each MET bin
		fout.write("jmax %d  number of backgrounds\n" %len(Bkgs)) #total bkg processes
		fout.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
		fout.write("bin bin%d\n" %m)
		Obs=0;
		for i in range(1, len(Bkgs)+1):Obs=Obs+FullBkgYieldsRes3b_highR[i][m]
		fout.write("observation %g\n" %Obs)
		fout.write("bin ");
		for i in range(0, len(Bkgs)+1):fout.write("bin%d " %m)
		fout.write("\nprocess sig");
		for b in Bkgs: fout.write(" %s" %b)
		fout.write("\nprocess")
		for i in range(0, len(Bkgs)+1):fout.write(" %d" %i)
		Sig=SigHisto3b_highR.GetBinContent(m+1)
	 	if(m==SigHisto3b_highR.GetNbinsX()-1):Sig=Sig+SigHisto3b_highR.GetBinContent(m+2)#ADD OVERFLOW Bin
		fout.write("\nrate %g " %(Sig));
		for i in range(1, len(Bkgs)+1):
			#fout.write(" %g " %FullBkgYieldsRes3b[i][m])
			if FullBkgYieldsRes3b_highR[i][m]>0:  fout.write(" %g " %FullBkgYieldsRes3b_highR[i][m])
			else:fout.write(" 0.00001 ")
		fout.close()
		ResolvedCards3b_highR=ResolvedCards3b_highR+" HiggsinoRes3b_highR_MET%d.txt" %m

	for m in range(0,SigHisto4b_lowR.GetNbinsX()):
		fout=open("HiggsinoRes4b_lowR_MET%d.txt" %m,'w')
		fout.seek(0)
		fout.write("# Simple counting experiment, with higgsino signal and a few background processes\n")
		fout.write("imax 1  number of channels\n") #one .txt card for each MET bin
		fout.write("jmax %d  number of backgrounds\n" %len(Bkgs)) #total bkg processes
		fout.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
		fout.write("bin bin%d\n" %m)
		Obs=0;
		for i in range(1, len(Bkgs)+1):Obs=Obs+FullBkgYieldsRes4b_lowR[i][m]
		fout.write("observation %g\n" %Obs)
		fout.write("bin ");
		for i in range(0, len(Bkgs)+1):fout.write("bin%d " %m)
		fout.write("\nprocess sig");
		for b in Bkgs: fout.write(" %s" %b)
		fout.write("\nprocess")
		for i in range(0, len(Bkgs)+1):fout.write(" %d" %i)
		Sig=SigHisto4b_lowR.GetBinContent(m+1)
	 	if(m==SigHisto4b_lowR.GetNbinsX()-1):Sig=Sig+SigHisto4b_lowR.GetBinContent(m+2)#ADD OVERFLOW Bin
		fout.write("\nrate %g " %(Sig));
		for i in range(1, len(Bkgs)+1):
			#fout.write(" %g " %FullBkgYieldsRes4b[i][m])
			if FullBkgYieldsRes4b_lowR[i][m]>0:  fout.write(" %g " %FullBkgYieldsRes4b_lowR[i][m])
			else:fout.write(" 0.00001 ")
		fout.close()
		ResolvedCards4b_lowR=ResolvedCards4b_lowR+" HiggsinoRes4b_lowR_MET%d.txt" %m

	for m in range(0,SigHisto4b_highR.GetNbinsX()):
		fout=open("HiggsinoRes4b_highR_MET%d.txt" %m,'w')
		fout.seek(0)
		fout.write("# Simple counting experiment, with higgsino signal and a few background processes\n")
		fout.write("imax 1  number of channels\n") #one .txt card for each MET bin
		fout.write("jmax %d  number of backgrounds\n" %len(Bkgs)) #total bkg processes
		fout.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
		fout.write("bin bin%d\n" %m)
		Obs=0;
		for i in range(1, len(Bkgs)+1):Obs=Obs+FullBkgYieldsRes4b_highR[i][m]
		fout.write("observation %g\n" %Obs)
		fout.write("bin ");
		for i in range(0, len(Bkgs)+1):fout.write("bin%d " %m)
		fout.write("\nprocess sig");
		for b in Bkgs: fout.write(" %s" %b)
		fout.write("\nprocess")
		for i in range(0, len(Bkgs)+1):fout.write(" %d" %i)
		Sig=SigHisto4b_highR.GetBinContent(m+1)
	 	if(m==SigHisto4b_highR.GetNbinsX()-1):Sig=Sig+SigHisto4b_highR.GetBinContent(m+2)#ADD OVERFLOW Bin
		fout.write("\nrate %g " %(Sig));
		for i in range(1, len(Bkgs)+1):
			#fout.write(" %g " %FullBkgYieldsRes4b[i][m])
			if FullBkgYieldsRes4b_highR[i][m]>0:  fout.write(" %g " %FullBkgYieldsRes4b_highR[i][m])
			else:fout.write(" 0.00001 ")
		fout.close()
		ResolvedCards4b_highR=ResolvedCards4b_highR+" HiggsinoRes4b_highR_MET%d.txt" %m

	print "signal: TChiHH%s " %s
	# os.system("combineCards.py %s %s > HiggsinoBoostedOnlyTChiHH%s.txt" %(DiJetCardsBoost,CardsBoost,s))
	# os.system("combineCards.py %s %s %s %s > HiggsinoResolvedOnlyTChiHH%s.txt" %(ResolvedCards4b_lowR,ResolvedCards4b_highR,ResolvedCards3b_lowR,ResolvedCards3b_highR,s))


	# os.system("combineCards.py HiggsinoBoostedOnlyTChiHH%s.txt HiggsinoResolvedOnlyTChiHH%s.txt > HiggsinoResToBoostOnlyTChiHH%s.txt" %(s,s,s))
	# os.system("combineCards.py HiggsinoBoostedOnlyTChiHH%s.txt HiggsinoResolvedOnlyTChiHH%s.txt %s %s > HiggsinoResToBoostFullComboTChiHH%s.txt " %(s,s,SemiResolvedCards,SingleJetCardsBoost,s))
	# os.system("combineCards.py HiggsinoBoostedOnlyTChiHH%s.txt HiggsinoResolvedOnlyTChiHH%s.txt %s %s > HiggsinoResToBoostFullTest2TChiHH%s.txt " %(s,s,SemiResolvedCards,SingleJetCardsBoost,s))
	os.system("combineCards.py %s %s %s %s %s %s %s %s > HiggsinoResToBoostFullTest1TChiHH%s.txt " %(DiJetCardsBoost,CardsBoost, ResolvedCards4b_lowR,ResolvedCards4b_highR,ResolvedCards3b_lowR,ResolvedCards3b_highR, SemiResolvedCards,SingleJetCardsBoost,s))


	#Options for combine
	# os.system("combine -M AsymptoticLimits HiggsinoResolvedOnlyTChiHH%s.txt -n HiggsinoResolvedOnlyTChiHH%s -m %s" %(s,s,s)) #Res only
	# os.system("combine -M AsymptoticLimits HiggsinoBoostedOnlyTChiHH%s.txt -n HiggsinoBoostedOnlyTChiHH%s -m %s" %(s,s,s)) #Boosted only
	# os.system("combine -M AsymptoticLimits HiggsinoResToBoostOnlyTChiHH%s.txt -n HiggsinoResToBoostOnlyTChiHH%s -m %s" %(s,s,s)) #Res to Boost, but no combo cases
	# os.system("combine -M AsymptoticLimits HiggsinoResToBoostFullComboTChiHH%s.txt -n HiggsinoResToBoostFullComboTChiHH%s -m %s" %(s,s,s)) #Res to Boost with combo cases
	os.system("combine -M AsymptoticLimits HiggsinoResToBoostFullTest1TChiHH%s.txt -n HiggsinoResToBoostFullTest1TChiHH%s -m %s" %(s,s,s)) #Res to Boost with combo cases, test cuts
